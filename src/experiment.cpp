/*
Comp Eng 3DY4 (Computer Systems Integration Project)

Copyright by Nicola Nicolici
Department of Electrical and Computer Engineering
McMaster University
Ontario, Canada
*/

#include "dy4.h"
#include "filter.h"
#include "fourier.h"
#include "genfunc.h"
#include "iofunc.h"
#include "logfunc.h"
#include "rfFrontend.h"
#include "mono.h"
#include "stereo.h"
#include "rds.h"
#include <chrono>
#include <queue>
#include <thread>
#include <mutex>
#include <condition_variable>

// Filter parameters
float rf_Fs = 2.4e6;
float rf_Fc = 100e3;
unsigned short int rf_taps = 151;
int rf_decim = 10;

float audio_Fs = 48e3;
float audio_Fc = 16e3;
unsigned short int audio_taps = 151;
int audio_decim = 5;

// Mode select upon user's input
int modeSelect(int argc, char *argv[])
{
    int mode = 0;
    if (argc < 2)
    {
        std::cerr << "Operating in default mode 0" << std::endl;
    }
    else if (argc == 2)
    {
        mode = atoi(argv[1]);
        if (mode != 1)
        {
            std::cerr << "Wrong mode" << mode << std::endl;
            exit(1);
        }
    }
    else
    {
        std::cerr << "Usage: " << argv[0] << std::endl;
        std::cerr << "Or " << std::endl;
        std::cerr << "Usage: " << argv[0] << " 1" << std::endl;
    }
    return mode;
}

// RF frontend thread
void frontendThread(std::queue<std::vector<float>> &rf_queue, std::mutex &rf_mutex, std::condition_variable &rf_cvar, int mode)
{
    int BLOCK_SIZE = (mode == 0) ? 1024 * rf_decim * audio_decim : 1000 * rf_decim * audio_decim;

    if (mode == 1)
        rf_Fs = 2.5e6;

    // Filter coefficients
    std::vector<float> rf_coeff;
    impulseResponseLPF(rf_Fs, rf_Fc, rf_taps, rf_coeff);
    // Filter buffers for rf frontend block-processing
    std::vector<float> temp_I(rf_taps - 1, 0.0);
    std::vector<float> temp_Q(rf_taps - 1, 0.0);
    std::vector<float> prev_IQ(2, 0.0);

    int block_id = 0;
    for (unsigned int block_id = 0; block_id < MAX_BLOCKS; block_id++)
    {
        // rf frontend process
        std::vector<float> fm_dem;
        rfFrontend(fm_dem, prev_IQ, rf_coeff, temp_I, temp_Q, rf_decim, BLOCK_SIZE);

        // lock thread if limit of queue blocks reached
        std::unique_lock<std::mutex> rf_lock(rf_mutex);
        if (rf_queue.size() == QUEUE_BLOCKS)
            rf_cvar.wait(rf_lock);
        rf_queue.push(fm_dem);
        rf_lock.unlock();
        rf_cvar.notify_one();
    }
    std::cerr << "Maximum block limit of " << MAX_BLOCKS << " reached" << std::endl;
    exit(1);
}

// Audio and RDS path thread
void pathThread(std::queue<std::vector<float>> &rf_queue, std::mutex &rf_mutex, std::condition_variable &rf_cvar, int audio_Fc, int mode)
{
    // mono path Block status
    std::vector<float> temp_audio(audio_taps - 1, 0.0);

    // stereo path Block status
    std::vector<float> path1_coeff, path2_coeff;
    impulseResponseBPF(rf_Fs, 19.5e3, 18.5e3, rf_taps, path1_coeff);
    impulseResponseBPF(rf_Fs, 54e3, 22e3, rf_taps, path2_coeff);

    // stereo PLL buffers for block-processing
    float integrator = 0.0;
    float phaseEst = 0.0;
    float feedbackI = 1.0;
    float feedbackQ = 0.0;
    float nco_buf = 1.0;
    int trigOffset = 0;

    // Filter buffers for stereo block-processing
    std::vector<float> temp_path1(rf_taps - 1, 0.0);
    std::vector<float> temp_path2(rf_taps - 1, 0.0);
    std::vector<float> temp_path_mono(rf_taps - 1, 0.0);

    // Audio low-pass filter coefficients
    std::vector<float> audio_coeff;
    if (mode == 0)
        impulseResponseLPF(24e4, audio_Fc, audio_taps, audio_coeff);
    else if (mode == 1)
        impulseResponseLPF(6e6, audio_Fc, audio_taps * 24, audio_coeff);

    // RDS PLL buffers for block-processing
    float rds_integrator = 0.0;
    float rds_phaseEst = 0.0;
    float rds_feedbackI = 1.0;
    float rds_feedbackQ = 0.0;
    float rds_nco_buf_i = 1.0;
    float rds_nco_buf_q = 1.0;
    int rds_trigOffset = 0;
    int up_decim = 19;
    int down_decim = 80;

    // Filter buffers for RDS block-processing
    std::vector<float> temp_after_3k_i(rf_taps - 1, 0.0);
    std::vector<float> temp_after_3k_q(rf_taps - 1, 0.0);
    std::vector<float> temp_channel_BPF(rf_taps - 1, 0.0);
    std::vector<float> temp_CR_BPF(rf_taps - 1, 0.0);
    std::vector<float> temp_Demod_LPF_i(rf_taps * up_decim - 1, 0.0);
    std::vector<float> temp_Demod_LPF_q(rf_taps * up_decim - 1, 0.0);
    std::vector<float> temp_rrc_i(rf_taps - 1, 0.0);
    std::vector<float> temp_rrc_q(rf_taps - 1, 0.0);

    // Filter coefficients
    std::vector<float> rds1_coeff, rds2_coeff, rrc_coeff, LP3k_coeff, LP16k_coeff; 
    impulseResponseBPF(rf_Fs / 10, 60e3, 54e3, rf_taps, rds1_coeff); // Band-pass Filter 54KHz to 60KHz
    impulseResponseBPF(rf_Fs / 10, 114.5e3, 113.5e3, rf_taps, rds2_coeff); // Band-pass Filter 113.5KHz to 114.5KHz
    impulseResponseLPF(rf_Fs / 10, 3e3, rf_taps, LP3k_coeff); // Low-pass Filter Fc=3KHz
    impulseResponseLPF(rf_Fs / 10 * up_decim, 16e3, rf_taps * up_decim, LP16k_coeff); // Low-pass Filter Fc=16KHz
    impulseResponseRRC(57e3, rf_taps, rrc_coeff); // Root-raised-cosine Filter

    // RDS data processing parameters
    int count_RRC = 0; // 3-block data counter
    int phase_index_i = 0; // phase index for data recovery
    int phase_index_q = 0; 
    int deco_buf_i = 0; // decoding buffer
    int multi_resize_count = 0; // RDS block counter
    std::vector<int> multi_prepared; // RDS bitstream container
    std::vector<float> RRC_3block_i(608 * 3, 0.0); // 3-block data buffer
    std::vector<float> RRC_3block_q(608 * 3, 0.0);
    int syn_base = 0; // bitstream offset
    int syn_base_buf = 0;

    int block_id = 0;
    for (unsigned int block_id = 0;; block_id++)
    {
        // std::cerr << "********* Processing block "<< block_id << std::endl;
        // multi-threading with RF Frontend
        std::unique_lock<std::mutex> rf_lock(rf_mutex);
        if (rf_queue.empty())
            rf_cvar.wait(rf_lock);
        std::vector<float> fm_dem = rf_queue.front();
        rf_queue.pop();

        // Audio & RDS multi-threading params
        std::queue<std::vector<float>> path_queue;
        std::mutex path_mutex;
        std::condition_variable path_cvar;

        std::thread t_mono(monoPath, std::ref(path_queue), std::ref(path_cvar),
                           std::ref(fm_dem), std::ref(audio_coeff), std::ref(temp_audio), std::ref(audio_decim), std::ref(mode));

        std::thread t_stereo(stereoPath, std::ref(path_queue), std::ref(path_mutex), std::ref(path_cvar),
                             std::ref(fm_dem), std::ref(temp_path1), std::ref(temp_path2), std::ref(temp_path_mono), std::ref(mode), std::ref(rf_Fs), std::ref(audio_coeff), std::ref(path1_coeff), std::ref(path2_coeff),
                             std::ref(integrator), std::ref(phaseEst), std::ref(feedbackI), std::ref(feedbackQ), std::ref(nco_buf), std::ref(trigOffset = 0));

        // RDS for mode 0 only
        if (mode == 0)
        {
            std::thread t_rds(rdsPath, std::ref(path_cvar), std::ref(fm_dem), std::ref(rf_Fs), std::ref(rds_integrator),
                              std::ref(rds_phaseEst), std::ref(rds_feedbackI), std::ref(rds_feedbackQ), std::ref(rds_nco_buf_i), std::ref(rds_nco_buf_q), std::ref(rds_trigOffset = 0),
                              std::ref(rds1_coeff), std::ref(rds2_coeff), std::ref(rrc_coeff), std::ref(LP3k_coeff), std::ref(LP16k_coeff),
                              std::ref(temp_after_3k_i), std::ref(temp_after_3k_q), std::ref(temp_channel_BPF), std::ref(temp_CR_BPF), std::ref(temp_Demod_LPF_i), std::ref(temp_Demod_LPF_q), std::ref(temp_rrc_i), std::ref(temp_rrc_q), std::ref(RRC_3block_i), std::ref(RRC_3block_q),
                              std::ref(block_id), std::ref(count_RRC), std::ref(phase_index_i), std::ref(phase_index_q), std::ref(deco_buf_i), std::ref(multi_resize_count), std::ref(multi_prepared), std::ref(syn_base), std::ref(syn_base_buf));
            t_rds.join();
        }

        t_mono.join();
        t_stereo.join();

        rf_lock.unlock();
        rf_cvar.notify_one();
    }
}

int main(int argc, char *argv[])
{
    int mode = modeSelect(argc, argv);

    // multi-threading params
    std::queue<std::vector<float>> rf_queue;
    std::mutex rf_mutex;
    std::condition_variable rf_cvar;

    std::thread t_frontend(frontendThread, std::ref(rf_queue), std::ref(rf_mutex), std::ref(rf_cvar), std::ref(mode));
    std::thread t_path(pathThread, std::ref(rf_queue), std::ref(rf_mutex), std::ref(rf_cvar), std::ref(audio_Fc), std::ref(mode));

    t_frontend.join();
    t_path.join();

    return 0;
}