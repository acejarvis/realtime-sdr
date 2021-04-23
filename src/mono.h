/*
Comp Eng 3DY4 (Computer Systems Integration Project)

Department of Electrical and Computer Engineering
McMaster University
Ontario, Canada
*/

#ifndef DY4_MONO_H
#define DY4_MONO_H

// add headers as needed
#include <vector>
#include <math.h>
#include <chrono>
#include <queue>
#include <mutex>
#include <condition_variable>

// declaration of a function prototypes
void monoPath(std::queue<std::vector<float>> &path_queue, std::condition_variable &path_cvar,
              std::vector<float> &fm_dem, std::vector<float> &audio_coeff, std::vector<float> &temp_audio, int audio_decim, int mode);

#endif // DY4_MONO_H