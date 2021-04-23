/*
Comp Eng 3DY4 (Computer Systems Integration Project)

Department of Electrical and Computer Engineering
McMaster University
Ontario, Canada
*/

#include "dy4.h"
#include "filter.h"
#include "stereo.h"
#include <math.h>
#include <chrono>
#include <queue>
#include <mutex>
#include <condition_variable>

void monoPath(std::queue<std::vector<float>> &path_queue, std::condition_variable &path_cvar, std::vector<float> &fm_dem, std::vector<float> &audio_coeff, std::vector<float> &temp_audio, int audio_decim, int mode)
{
    // mono
    std::vector<float> audio_filt;
    if (mode == 0)
    {
        downConvolveFIR(audio_filt, fm_dem, audio_coeff, temp_audio, audio_decim);
    }
    else if (mode == 1)
    {
        int demod_count = 0;
        std::vector<float> upsampled_fm_dem(24 * fm_dem.size(), 0.0);
        for (auto i = 0; i < fm_dem.size(); i++)
        {
            upsampled_fm_dem[i * 24] = 24 * fm_dem[i];
        }
        upDownConvolveFIR(audio_filt, upsampled_fm_dem, audio_coeff, temp_audio, 24, audio_decim * 25);
    }

    path_queue.push(audio_filt);
    path_cvar.notify_one();
}