/*
Comp Eng 3DY4 (Computer Systems Integration Project)

Department of Electrical and Computer Engineering
McMaster University
Ontario, Canada
*/

#include "dy4.h"
#include "filter.h"
#include <math.h>
#include <queue>
#include <thread>
#include <mutex>
#include <condition_variable>

void stereoPath(std::queue<std::vector<float>> &path_queue, std::mutex &path_mutex, std::condition_variable &path_cvar,
                std::vector<float> &fm_dem, std::vector<float> &temp_path1, std::vector<float> &temp_path2, std::vector<float> &temp_path_mono,
                int mode, float rf_Fs, const std::vector<float> &audio_coeff, std::vector<float> &path1_coeff, std::vector<float> &path2_coeff,
                float &integrator, float &phaseEst, float &feedbackI, float &feedbackQ, float &nco_buf, int &trigOffset)
{
    std::vector<float> path1_filt, path1_out, path2_out;

    // --- Stereo Carrier Recovery ---
    // Band-pass filter 18.5KHz to 19.5KHz
    downConvolveFIR(path1_filt, fm_dem, path1_coeff, temp_path1, 1);
    // PLL and NCO
    float PLL_freq = 19e3;
    fmPLL(path1_filt, path1_out, PLL_freq, rf_Fs, integrator, phaseEst, feedbackI, feedbackQ, nco_buf, trigOffset, 2.0, 0.0, 0.005);

    // --- Stereo Channel Extraction ---
    // Band-pass filter 22KHz to 54KHz
    downConvolveFIR(path2_out, fm_dem, path2_coeff, temp_path2, 1);

    // --- Stereo Processing ---
    // Mixer
    std::vector<float> mixed_data(path1_filt.size(), 0.0);
    for (int i = 0; i < path1_filt.size(); i++)
    {
        mixed_data[i] = path1_out[i] * path2_out[i];
    }

    // Digital Filtering Sample Rate Conversion
    std::vector<float> audio_filt;
    float audio_decim = 5;
    if (mode == 0)
    {
        downConvolveFIR(audio_filt, mixed_data, audio_coeff, temp_path_mono, audio_decim);
    }
    else if (mode == 1)
    {
        int demod_count = 0;
        std::vector<float> upsampled_fm_dem(24 * fm_dem.size(), 0.0);
        for (auto i = 0; i < fm_dem.size(); i++)
        {
            upsampled_fm_dem[i * 24] = fm_dem[i];
        }
        upDownConvolveFIR(audio_filt, upsampled_fm_dem, audio_coeff, temp_path_mono, 24, audio_decim * 25);
    }

    // Stereo Combiner
    std::vector<float> stereo_L(audio_filt.size(), 0.0);
    std::vector<float> stereo_R(audio_filt.size(), 0.0);
    // wait for mono to finish
    std::unique_lock<std::mutex> path_lock(path_mutex);
    if (path_queue.empty())
        path_cvar.wait(path_lock);

    std::vector<float> mono_data = path_queue.front();
    path_queue.pop();

    for (int j = 0; j < audio_filt.size(); j++)
    {
        stereo_L[j] = (audio_filt[j] + mono_data[j]) / 2;
        stereo_R[j] = (audio_filt[j] - mono_data[j]) / 2;
    }
    std::vector<short int> audio_out(stereo_L.size() * 2, 0.0);

    // audio output at 48K samples/sec
    int stereo_count = 0;
    for (unsigned int k = 0; k < stereo_L.size(); k++)
    {
        if (std::isnan(stereo_L[k]))
            audio_out[stereo_count] = 0;
        else
            audio_out[stereo_count] = static_cast<short int>(stereo_L[k] * 16384);

        stereo_count++;

        if (std::isnan(stereo_R[k]))
            audio_out[stereo_count] = 0;
        else
            audio_out[stereo_count] = static_cast<short int>(stereo_R[k] * 16384);

        stereo_count++;
    }
    fwrite(&audio_out[0], sizeof(short int), audio_out.size(), stdout);

    path_lock.unlock();
    path_cvar.notify_one();
}