# mod-cabsim-IR-loader

An LV2 cabinet simulator plugin that loads impulse response (IR) files.

This plugin is specifically created for handling speaker cabinet IRs,
this plugin is not optimized for handling larger files like reverb IRs.

Currently it only uses the first 42.7 ms (2048 samples at 48 kHz sampling rate) of the loaded IR file.
IR files at different sample rates are resampled to 48 kHz by the plugin.
It is recommended to trim any silence at the start of the IR file for optimal results.

Default IR file provided by forward audio.
