/*
  The file loading is using code from the LV2 Sampler Example Plugin
  Copyright 2011-2012 David Robillard <d@drobilla.net>
  Copyright 2011 Gabriel M. Beddingfield <gabriel@teuton.org>
  Copyright 2011 James Morris <jwm.art.net@gmail.com>

  Resampler function is taken from https://github.com/cpuimage/resampler/
  Copyright (c) 2019 Zhihan Gao

  Permission to use, copy, modify, and/or distribute this software for any
  purpose with or without fee is hereby granted, provided that the above
  copyright notice and this permission notice appear in all copies.

  THIS SOFTWARE IS PROVIDED "AS IS" AND THE AUTHOR DISCLAIMS ALL WARRANTIES
  WITH REGARD TO THIS SOFTWARE INCLUDING ALL IMPLIED WARRANTIES OF
  MERCHANTABILITY AND FITNESS. IN NO EVENT SHALL THE AUTHOR BE LIABLE FOR
  ANY SPECIAL, DIRECT, INDIRECT, OR CONSEQUENTIAL DAMAGES OR ANY DAMAGES
  WHATSOEVER RESULTING FROM LOSS OF USE, DATA OR PROFITS, WHETHER IN AN
  ACTION OF CONTRACT, NEGLIGENCE OR OTHER TORTIOUS ACTION, ARISING OUT OF
  OR IN CONNECTION WITH THE USE OR PERFORMANCE OF THIS SOFTWARE.
*/

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <stdbool.h>
#include <math.h>
#include "fftw3.h"
#include <math.h>
#include <stdlib.h>
#include <string.h>

#ifndef __cplusplus
#    include <stdbool.h>
#endif

#include <sndfile.h>

#include "lv2/lv2plug.in/ns/ext/atom/forge.h"
#include "lv2/lv2plug.in/ns/ext/atom/util.h"
#include "lv2/lv2plug.in/ns/ext/log/log.h"
#include "lv2/lv2plug.in/ns/ext/log/logger.h"
#include "lv2/lv2plug.in/ns/ext/midi/midi.h"
#include "lv2/lv2plug.in/ns/ext/patch/patch.h"
#include "lv2/lv2plug.in/ns/ext/state/state.h"
#include "lv2/lv2plug.in/ns/ext/urid/urid.h"
#include "lv2/lv2plug.in/ns/ext/worker/worker.h"
#include "lv2/lv2plug.in/ns/lv2core/lv2.h"

#include "./uris.h"
#include "./circular_buffer.h"

#define REAL 0
#define IMAG 1

#define MAX_FFT_SIZE 2048

//macro for Volume in DB to a coefficient
#define DB_CO(g) ((g) > -90.0f ? powf(10.0f, (g) * 0.05f) : 0.0f)

enum {
    CABSIM_CONTROL = 0,
    CABSIM_NOTIFY  = 1,
    CABSIM_IN      = 2,
    CABSIM_OUT     = 3,
    ATTENUATE      = 4
};

//static const char* default_sample_file = "Orange_PPC412_V30_412_C_Hi-Gn_121+57_Celestion.wav";

typedef struct {
    SF_INFO  info;      // Info about sample from sndfile
    float*   data;      // ImpulseResponse data in float
    char*    path;      // Path of file
    uint32_t path_len;  // Length of path
} ImpulseResponse;

typedef struct {
    // Features
    LV2_URID_Map*        map;
    LV2_Worker_Schedule* schedule;
    LV2_Log_Log*         log;

    // Forge for creating atoms
    LV2_Atom_Forge forge;

    // Logger convenience API
    LV2_Log_Logger logger;

    // ImpulseResponse
    ImpulseResponse* ir;

    // Ports
    const LV2_Atom_Sequence* control_port;
    LV2_Atom_Sequence*       notify_port;
    float*                   output_port;
    float*                   input_port;

    // Forge frame for notify port (for writing worker replies)
    LV2_Atom_Forge_Frame notify_frame;

    // URIs
    CabsimURIs uris;

    // Current position in run()
    uint32_t frame_offset;

    // Playback state
    sf_count_t frame;
    bool       play;
    bool       new_ir;
    bool       ir_loaded;

    double samplerate;

    //CABSIM ========================================

    float *outbuf;
    float *inbuf;
    float *IR;

    uint32_t prev_buffer_size;
    const float *attenuation;
    uint32_t valid_IR_data_length;

    fftwf_complex *outComplex;
    fftwf_complex *IRout;
    fftwf_complex *convolved;

    fftwf_plan fft;
    fftwf_plan ifft;
    fftwf_plan IRfft;

    int overlap_add_buffers;
    ringbuffer_t overlap_buffer;
} Cabsim;

typedef struct {
    LV2_Atom atom;
    ImpulseResponse*  ir;
} ImpulseResponseMessage;

static uint64_t
Resample_f32(const float *input, float *output, int inSampleRate,
        int outSampleRate, uint64_t inputSize, uint32_t channels)
{
    if (input == NULL)
        return 0;
    uint64_t outputSize = inputSize * outSampleRate / inSampleRate;
    if (output == NULL)
        return outputSize;
    double stepDist = ((double) inSampleRate / (double) outSampleRate);
    const uint64_t fixedFraction = (1LL << 32);
    const double normFixed = (1.0 / (1LL << 32));
    uint64_t step = ((uint64_t) (stepDist * fixedFraction + 0.5));
    uint64_t curOffset = 0;
    for (uint32_t i = 0; i < outputSize; i += 1) {
        for (uint32_t c = 0; c < channels; c += 1) {
            *output++ = (float) (input[c] + (input[c + channels] - input[c]) * (
                        (double) (curOffset >> 32) + ((curOffset & (fixedFraction - 1)) * normFixed)));
        }
        curOffset += step;
        input += (curOffset >> 32) * channels;
        curOffset &= (fixedFraction - 1);
    }
    return outputSize;
}

static sf_count_t
convert_to_mono(float *data, sf_count_t num_input_frames, uint32_t channels)
{
    sf_count_t mono_index = 0;
    for (sf_count_t i = 0; i < num_input_frames * channels; i+=channels) {
        data[mono_index++] = data[i];
    }

    sf_count_t num_output_frames = mono_index;

    return num_output_frames;
}

/**
   Load a new ir and return it.

   Since this is of course not a real-time safe action, this is called in the
   worker thread only.  The ir is loaded and returned only, plugin state is
   not modified.
*/
static ImpulseResponse*
load_ir(Cabsim* self, const char* path, uint32_t path_len)
{
    char* irpath = (char*)malloc(path_len + 1);
    memcpy(irpath, path, path_len);
    irpath[path_len] = 0;

    lv2_log_trace(&self->logger, "Loading ir %s\n", irpath);

    ImpulseResponse* const  ir  = (ImpulseResponse*)calloc(1, sizeof(ImpulseResponse));
    SF_INFO* const info    = &ir->info;
    SNDFILE* const sndfile = sf_open(irpath, SFM_READ, info);

    if (!sndfile || !info->frames) {
        lv2_log_error(&self->logger, "Failed to open ir '%s'\n", irpath);
        free(ir);
        return NULL;
    }

    // Read data
    float* const data = malloc(sizeof(float) * (info->frames * info->channels));
    if (!data) {
        lv2_log_error(&self->logger, "Failed to allocate memory for ir\n");
        return NULL;
    }
    sf_seek(sndfile, 0ul, SEEK_SET);
    sf_read_float(sndfile, data, info->frames * info->channels);
    sf_close(sndfile);

    //When IR has multiple channels, only use first channel
    if (info->channels != 1) {
        info->frames = convert_to_mono(data, info->frames, info->channels);
        info->channels = 1;
    }

    //apply samplerate conversion if needed
    if (info->samplerate == (int)self->samplerate) {
        ir->data = data;
    } else {
        uint64_t targetSampleCount = Resample_f32(data, 0, info->samplerate, (int)self->samplerate, (uint64_t)info->frames, 1);
        float* const resampled_data = malloc(targetSampleCount * sizeof(float));
        info->frames = Resample_f32(data, resampled_data, info->samplerate, (int)self->samplerate, (uint64_t)info->frames, 1);
        free(data);
        ir->data = resampled_data;
    }

    // Fill ir struct and return it
    ir->path     = irpath;
    ir->path_len = path_len;
    return ir;
}

static void
free_ir(Cabsim* self, ImpulseResponse* ir)
{
    if (ir) {
        lv2_log_trace(&self->logger, "Freeing %s\n", ir->path);
        free(ir->path);
        free(ir->data);
        free(ir);
    }
}

/**
   Do work in a non-realtime thread.

   This is called for every piece of work scheduled in the audio thread using
   self->schedule->schedule_work().  A reply can be sent back to the audio
   thread using the provided respond function.
*/
static LV2_Worker_Status
work(LV2_Handle                  instance,
     LV2_Worker_Respond_Function respond,
     LV2_Worker_Respond_Handle   handle,
     uint32_t                    size,
     const void*                 data)
{
    Cabsim*        self = (Cabsim*)instance;
    const LV2_Atom* atom = (const LV2_Atom*)data;
    if (atom->type == self->uris.cab_freeImpulseResponse) {
        // Free old ir
        const ImpulseResponseMessage* msg = (const ImpulseResponseMessage*)data;
        free_ir(self, msg->ir);
    } else {
        // Handle set message (load ir).
        const LV2_Atom_Object* obj = (const LV2_Atom_Object*)data;

        // Get file path from message
        const LV2_Atom* file_path = read_set_file(&self->uris, obj);
        if (!file_path) {
            return LV2_WORKER_ERR_UNKNOWN;
        }

        // Load ir.
        ImpulseResponse* ir = load_ir(self, LV2_ATOM_BODY_CONST(file_path), file_path->size);
        if (ir) {
            // Loaded ir, send it to run() to be applied.
            respond(handle, sizeof(ir), &ir);
        }
    }

    return LV2_WORKER_SUCCESS;
}

/**
   Handle a response from work() in the audio thread.

   When running normally, this will be called by the host after run().  When
   freewheeling, this will be called immediately at the point the work was
   scheduled.
*/
static LV2_Worker_Status
work_response(LV2_Handle  instance,
              uint32_t    size,
              const void* data)
{
    Cabsim* self = (Cabsim*)instance;

    ImpulseResponseMessage msg = { { sizeof(ImpulseResponse*), self->uris.cab_freeImpulseResponse },
        self->ir };

    // Send a message to the worker to free the current ir
    self->schedule->schedule_work(self->schedule->handle, sizeof(msg), &msg);

    // Install the new ir
    self->ir = *(ImpulseResponse*const*)data;

    self->new_ir = true;
    self->ir_loaded = false;

    return LV2_WORKER_SUCCESS;
}

static void
connect_port(LV2_Handle instance,
             uint32_t   port,
             void*      data)
{
    Cabsim* self = (Cabsim*)instance;
    switch (port) {
        case CABSIM_CONTROL:
            self->control_port = (const LV2_Atom_Sequence*)data;
            break;
        case CABSIM_NOTIFY:
            self->notify_port = (LV2_Atom_Sequence*)data;
            break;
        case CABSIM_IN:
            self->input_port  = (float*)data;
            break;
        case CABSIM_OUT:
            self->output_port = (float*)data;
            break;
        case ATTENUATE:
            self->attenuation = (const float*) data;
            break;
        default:
            break;
    }
}

static LV2_Handle
instantiate(const LV2_Descriptor*     descriptor,
            double                    rate,
            const char*               path,
            const LV2_Feature* const* features)
{
    // Allocate and initialise instance structure.
    Cabsim* self = (Cabsim*)calloc(1, sizeof(Cabsim));
    if (!self) {
        return NULL;
    }

    self->samplerate = rate;
    // Get host features
    for (int i = 0; features[i]; ++i) {
        if (!strcmp(features[i]->URI, LV2_URID__map)) {
            self->map = (LV2_URID_Map*)features[i]->data;
        } else if (!strcmp(features[i]->URI, LV2_WORKER__schedule)) {
            self->schedule = (LV2_Worker_Schedule*)features[i]->data;
        } else if (!strcmp(features[i]->URI, LV2_LOG__log)) {
            self->log = (LV2_Log_Log*)features[i]->data;
        }
    }
    if (!self->map) {
        lv2_log_error(&self->logger, "Missing feature urid:map\n");
        goto fail;
    } else if (!self->schedule) {
        lv2_log_error(&self->logger, "Missing feature work:schedule\n");
        goto fail;
    }

    // Map URIs and initialise forge/logger
    map_cabsim_uris(self->map, &self->uris);
    lv2_atom_forge_init(&self->forge, self->map);
    lv2_log_logger_init(&self->logger, self->map, self->log);

    self->outComplex = (fftwf_complex*) fftwf_malloc(sizeof(fftwf_complex)*(MAX_FFT_SIZE));
    self->IRout =  (fftwf_complex*) fftwf_malloc(sizeof(fftwf_complex)*(MAX_FFT_SIZE));
    self->convolved =  (fftwf_complex*) fftwf_malloc(sizeof(fftwf_complex)*(MAX_FFT_SIZE));
    self->outbuf = (float *) calloc((MAX_FFT_SIZE),sizeof(float));
    self->inbuf = (float *) calloc((MAX_FFT_SIZE),sizeof(float));
    self->IR = (float *) calloc((MAX_FFT_SIZE),sizeof(float));

    if (fftwf_import_system_wisdom() != 0)
    {
        self->fft = fftwf_plan_dft_r2c_1d(MAX_FFT_SIZE, self->inbuf, self->outComplex, FFTW_WISDOM_ONLY|FFTW_ESTIMATE);
        self->IRfft = fftwf_plan_dft_r2c_1d(MAX_FFT_SIZE, self->IR, self->IRout, FFTW_WISDOM_ONLY|FFTW_ESTIMATE);
        self->ifft = fftwf_plan_dft_c2r_1d(MAX_FFT_SIZE, self->convolved, self->outbuf, FFTW_WISDOM_ONLY|FFTW_ESTIMATE);
        lv2_log_note(&self->logger, "wisdom file loaded from system\n");
    } else {
        self->fft = fftwf_plan_dft_r2c_1d(MAX_FFT_SIZE, self->inbuf, self->outComplex, FFTW_ESTIMATE);
        self->IRfft = fftwf_plan_dft_r2c_1d(MAX_FFT_SIZE, self->IR, self->IRout, FFTW_ESTIMATE);
        self->ifft = fftwf_plan_dft_c2r_1d(MAX_FFT_SIZE, self->convolved, self->outbuf, FFTW_ESTIMATE);
        lv2_log_warning(&self->logger, "failed to import system wisdom file\n");
    }

    self->new_ir = false;
    self->ir_loaded = false;
    self->prev_buffer_size = UINT32_MAX;

    for (int i = 0; i < MAX_FFT_SIZE; i++)
        self->inbuf[i] = 0.0f;

    return (LV2_Handle)self;

fail:
    free(self);
    return 0;
}

static void
cleanup(LV2_Handle instance)
{
    Cabsim* self = (Cabsim*)instance;

    fftwf_destroy_plan(self->fft);
    fftwf_destroy_plan(self->IRfft);
    fftwf_destroy_plan(self->ifft);
    fftwf_free(self->outComplex);
    fftwf_free(self->IRout);
    fftwf_free(self->convolved);
    free(self->outbuf);
    free(self->inbuf);
    free(self->IR);
    free_ir(self, self->ir);
    free(self);
}

/** Define a macro for converting a gain in dB to a coefficient. */
#define DB_CO(g) ((g) > -90.0f ? powf(10.0f, (g) * 0.05f) : 0.0f)

static void
run(LV2_Handle instance,
    uint32_t   n_frames)
{
    Cabsim*     self   = (Cabsim*)instance;
    CabsimURIs* uris   = &self->uris;
    float*      input  = self->input_port;
    float*      output = self->output_port;

    float *outbuf  = self->outbuf;
    float *inbuf   = self->inbuf;
    float *IR      = self->IR;

    if (n_frames > MAX_FFT_SIZE) {
        // unsupported
        memset(output, 0, sizeof(float)*n_frames);
        return;
    }

    // Set up forge to write directly to notify output port.
    const uint32_t notify_capacity = self->notify_port->atom.size;
    lv2_atom_forge_set_buffer(&self->forge,
            (uint8_t*)self->notify_port,
            notify_capacity);

    // Start a sequence in the notify output port.
    lv2_atom_forge_sequence_head(&self->forge, &self->notify_frame, 0);

    // Read incoming events
    LV2_ATOM_SEQUENCE_FOREACH(self->control_port, ev) {
        self->frame_offset = ev->time.frames;
        if (lv2_atom_forge_is_object_type(&self->forge, ev->body.type)) {
            const LV2_Atom_Object* obj = (const LV2_Atom_Object*)&ev->body;
            if (obj->body.otype == uris->patch_Set) {
                // Get the property and value of the set message
                const LV2_Atom* property = NULL;
                const LV2_Atom* value    = NULL;
                lv2_atom_object_get(obj,
                        uris->patch_property, &property,
                        uris->patch_value,    &value,
                        0);
                if (!property) {
                    lv2_log_error(&self->logger,
                            "patch:Set message with no property\n");
                    continue;
                } else if (property->type != uris->atom_URID) {
                    lv2_log_error(&self->logger,
                            "patch:Set property is not a URID\n");
                    continue;
                }

                const uint32_t key = ((const LV2_Atom_URID*)property)->body;
                if (key == uris->cab_ir) {
                    // ImpulseResponse change, send it to the worker.
                    lv2_log_trace(&self->logger, "Queueing set message\n");
                    self->schedule->schedule_work(self->schedule->handle,
                            lv2_atom_total_size(&ev->body),
                            &ev->body);
                }
            } else {
                lv2_log_trace(&self->logger,
                        "Unknown object type %d\n", obj->body.otype);
            }
        } else {
            lv2_log_trace(&self->logger,
                    "Unknown event type %d\n", ev->body.type);
        }
    }

    // CABSIM =================================================================

    const float attenuation = *self->attenuation;

    const float coef = DB_CO(attenuation);

    uint32_t i, j, m;

    if (n_frames != self->prev_buffer_size) {
        switch (n_frames) {
            case 64:
                self->overlap_add_buffers = 32;
            break;
            case 128:
                self->overlap_add_buffers = 16;
            break;
            case 256:
                self->overlap_add_buffers = 8;
            break;
            case 512:
                self->overlap_add_buffers = 4;
            break;
            default:
                lv2_log_warning(&self->logger, "Non standard buffer size: '%i'. alliasing may happen\n", n_frames);

                self->overlap_add_buffers = 16;
            break;
        }

        self->prev_buffer_size = n_frames;

        //reset ringbuffer
        //max ringbuffer size is defined in ringbuffer.h, it should always be > max_overlap_add_buffers * MAX_FFT_SZIE
        ringbuffer_clear(&self->overlap_buffer, self->overlap_add_buffers * MAX_FFT_SIZE);

        memset(self->outbuf, 0, MAX_FFT_SIZE*sizeof(float));
        memset(self->inbuf, 0, MAX_FFT_SIZE*sizeof(float));
    }

    //copy inputbuffer and IR buffer with zero padding.
    if (self->new_ir)
    {
        self->valid_IR_data_length = (MAX_FFT_SIZE - n_frames) < self->ir->info.frames ?
                                  (MAX_FFT_SIZE - n_frames) : self->ir->info.frames;
        memcpy(IR, self->ir->data, self->valid_IR_data_length*sizeof(float));
        memset(IR+self->valid_IR_data_length, 0, (MAX_FFT_SIZE - self->valid_IR_data_length)*sizeof(float));

        fftwf_execute(self->IRfft);

        lv2_log_trace(&self->logger, "Responding to get request\n");
        lv2_atom_forge_frame_time(&self->forge, self->frame_offset);
        write_set_file(&self->forge, &self->uris,
                self->ir->path,
                self->ir->path_len);

        self->new_ir = false;
        self->ir_loaded = true;

        ringbuffer_clear(&self->overlap_buffer, self->overlap_add_buffers * MAX_FFT_SIZE);

        memset(self->outbuf, 0, MAX_FFT_SIZE*sizeof(float));
        memset(self->inbuf, 0, MAX_FFT_SIZE*sizeof(float));
    }

    for (i = 0; i < n_frames; i++)
        inbuf[i] = input[i] * coef * 0.2f;

    fftwf_execute(self->fft);

    if (self->ir_loaded) {

        //complex multiplication
        for (m = 0; m < self->valid_IR_data_length; m++) {
            self->convolved[m][REAL] = self->outComplex[m][REAL] * self->IRout[m][REAL] - self->outComplex[m][IMAG] * self->IRout[m][IMAG];
            self->convolved[m][IMAG] = self->outComplex[m][REAL] * self->IRout[m][IMAG] + self->outComplex[m][IMAG] * self->IRout[m][REAL];
        }
        memset(self->convolved+m, 0, (MAX_FFT_SIZE-m)*sizeof(fftwf_complex));

        fftwf_execute(self->ifft);

        for (j = 0; j < MAX_FFT_SIZE; j++)
            ringbuffer_push_sample(&self->overlap_buffer , outbuf[j] / MAX_FFT_SIZE);

        //normalize output with overlap add.
        for (j = 0; j < n_frames; j++) {
            float overlap_value = 0.0f;

            for (uint8_t Oa = 0; Oa < self->overlap_add_buffers-1; Oa++)
                overlap_value += ringbuffer_get_relative_val(&self->overlap_buffer, (Oa * MAX_FFT_SIZE) + ((self->overlap_add_buffers - Oa - 1) * n_frames) + j + 1);

            output[j] = ((outbuf[j] / MAX_FFT_SIZE) + overlap_value);
        }
    } else {
        memset(output, 0, sizeof(float)*n_frames);
    }
}

static LV2_State_Status
save(LV2_Handle                instance,
     LV2_State_Store_Function  store,
     LV2_State_Handle          handle,
     uint32_t                  flags,
     const LV2_Feature* const* features)
{
    Cabsim* self = (Cabsim*)instance;

    if (!self->ir) {
        return LV2_STATE_SUCCESS;
    }

    LV2_State_Map_Path* map_path = NULL;
    for (int i = 0; features[i]; ++i) {
        if (!strcmp(features[i]->URI, LV2_STATE__mapPath)) {
            map_path = (LV2_State_Map_Path*)features[i]->data;
        }
    }

    if (map_path) {
        char* apath = map_path->abstract_path(map_path->handle, self->ir->path);
        store(handle,
                self->uris.cab_ir,
                apath,
                strlen(self->ir->path) + 1,
                self->uris.atom_Path,
                LV2_STATE_IS_POD | LV2_STATE_IS_PORTABLE);
        free(apath);
        return LV2_STATE_SUCCESS;
    } else {
        return LV2_STATE_ERR_NO_FEATURE;
    }
}

static LV2_State_Status
restore(LV2_Handle                  instance,
        LV2_State_Retrieve_Function retrieve,
        LV2_State_Handle            handle,
        uint32_t                    flags,
        const LV2_Feature* const*   features)
{
    Cabsim* self = (Cabsim*)instance;

    size_t   size;
    uint32_t type;
    uint32_t valflags;

    const void* value = retrieve(
            handle,
            self->uris.cab_ir,
            &size, &type, &valflags);

    if (value) {
        const char* path = (const char*)value;
        ImpulseResponse *ir = load_ir(self, path, size);
        if (ir) {
            lv2_log_trace(&self->logger, "Restoring file %s\n", path);
            free_ir(self, self->ir);
            self->ir = ir;
            self->new_ir = true;
        } else {
            lv2_log_error(&self->logger, "File %s couldn't be loaded\n", path);
            return LV2_STATE_ERR_UNKNOWN;
        }
    }

    return LV2_STATE_SUCCESS;
}

static const void*
extension_data(const char* uri)
{
    static const LV2_State_Interface  state  = { save, restore };
    static const LV2_Worker_Interface worker = { work, work_response, NULL };
    if (!strcmp(uri, LV2_STATE__interface)) {
        return &state;
    } else if (!strcmp(uri, LV2_WORKER__interface)) {
        return &worker;
    }
    return NULL;
}

static const LV2_Descriptor descriptor = {
    CABSIM_URI,
    instantiate,
    connect_port,
    NULL,  // activate,
    run,
    NULL,  // deactivate,
    cleanup,
    extension_data
};

LV2_SYMBOL_EXPORT
const LV2_Descriptor* lv2_descriptor(uint32_t index)
{
    switch (index) {
        case 0:
            return &descriptor;
        default:
            return NULL;
    }
}
