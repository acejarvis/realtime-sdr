/*
Comp Eng 3DY4 (Computer Systems Integration Project)

Department of Electrical and Computer Engineering
McMaster University
Ontario, Canada
*/

#ifndef DY4_STEREO_H
#define DY4_STEREO_H

// add headers as needed
#include <vector>
#include <math.h>
#include <chrono>
#include <queue>
#include <thread>
#include <mutex>
#include <condition_variable>

// declaration of a function prototypes
void stereoPath(std::queue<std::vector<float>> &path_queue, std::mutex &path_mutex, std::condition_variable &path_cvar,
                std::vector<float> &fm_dem, std::vector<float> &temp_path1, std::vector<float> &temp_path2, std::vector<float> &temp_path_mono,
                int mode, float rf_Fs, const std::vector<float> &audio_coeff, std::vector<float> &path1_coeff, std::vector<float> &path2_coeff,
                float &integrator, float &phaseEst, float &feedbackI, float &feedbackQ, float &nco_buf, int &trigOffset);

#endif // DY4_STEREO_H