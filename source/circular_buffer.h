#ifndef CIRCULAR_BUFFER_H
#define CIRCULAR_BUFFER_H

#include <stdint.h>

#define MAX_BUFFER_SIZE 65536

typedef struct RINGBUFFER_T {
	float m_buffer[MAX_BUFFER_SIZE];
	uint32_t S;
	uint32_t m_size;
	uint32_t m_front;
	uint32_t m_back;
	float power;
} ringbuffer_t;

void ringbuffer_clear(ringbuffer_t *buffer, uint32_t size);
void ringbuffer_push(ringbuffer_t *buffer);
void ringbuffer_push_sample(ringbuffer_t *buffer, const float x);
void ringbuffer_pop(ringbuffer_t *buffer);
void ringbuffer_back_erase(ringbuffer_t *buffer, const uint32_t n);
void ringbuffer_front_erase(ringbuffer_t *buffer, const uint32_t n);
int ringbuffer_peek_index(ringbuffer_t *buffer);
float ringbuffer_push_and_calculate_power(ringbuffer_t *buffer, const float input);
float ringbuffer_front(ringbuffer_t *buffer);
float ringbuffer_back(ringbuffer_t *buffer);
float ringbuffer_get_relative_val(ringbuffer_t *buffer, uint32_t index);
float ringbuffer_get_index_val(ringbuffer_t *buffer, uint32_t index);
int ringbuffer_empty(ringbuffer_t *buffer);
int ringbuffer_full(ringbuffer_t *buffer);
float * ringbuffer_get_first_pointer(ringbuffer_t *buffer);

#endif // __RINGBUFFER_H__
