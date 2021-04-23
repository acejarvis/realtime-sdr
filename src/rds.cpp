/*
Comp Eng 3DY4 (Computer Systems Integration Project)

Department of Electrical and Computer Engineering
McMaster University
Ontario, Canada
*/

#include "dy4.h"
#include "filter.h"
#include "logfunc.h"
#include <math.h>
#include <mutex>
#include <condition_variable>
#include <algorithm>
#include <functional>

void rdsPath(std::condition_variable &path_cvar,
             std::vector<float> &fm_dem,
             float rf_Fs, float &integrator, float &phaseEst, float &feedbackI, float &feedbackQ, float &nco_buf_i, float &nco_buf_q, int &trigOffset,
             std::vector<float> &rds1_coeff, std::vector<float> &rds2_coeff, std::vector<float> &rrc_coeff, std::vector<float> &LP3k_coeff, std::vector<float> &LP16k_coeff,
             std::vector<float> &temp_after_3k_i, std::vector<float> &temp_after_3k_q, std::vector<float> &temp_channel_BPF, std::vector<float> &temp_CR_BPF, std::vector<float> &temp_Demod_LPF_i, std::vector<float> &temp_Demod_LPF_q,
             std::vector<float> &temp_rrc_i, std::vector<float> &temp_rrc_q, std::vector<float> &RRC_3block_i, std::vector<float> &RRC_3block_q,
             int block_id, int &count_RRC, int &phase_index_i, int &phase_index_q, int &deco_buf_i, int &multi_resize_count, std::vector<int> &multi_prepared, int &syn_base, int &syn_base_buf)
{
    // --- RDS channel Extraction ---
    // Band-pass Filter 54KHz to 60KHz

    std::vector<float> rds_channel;
    downConvolveFIR(rds_channel, fm_dem, rds1_coeff, temp_channel_BPF, 1);

    // --- RDS Carrier Recovery ---
    // Squaring Nonlinearity
    std::vector<float> post_channel(rds_channel.size(), 0.0);
    transform(rds_channel.begin(), rds_channel.end(), rds_channel.begin(), post_channel.begin(), std::multiplies<float>());

    // Band-pass Filter 113.5KHz to 114.5KHz

    std::vector<float> filter_out;
    downConvolveFIR(filter_out, post_channel, rds2_coeff, temp_CR_BPF, 1);

    // PLL & NCO
    std::vector<float> carrier_out_i, carrier_out_q;
    float PLL_freq = 114e3;
    fmPLL(filter_out, carrier_out_i, carrier_out_q, PLL_freq, rf_Fs / 10, integrator, phaseEst, feedbackI, feedbackQ, nco_buf_i, nco_buf_q, trigOffset, 0.5, 0.0, 0.01);

    // --- RDS Demodulation ---
    // Mixer
    std::vector<float> mixed_data_i(rds_channel.size(), 0.0);
    std::vector<float> mixed_data_q(rds_channel.size(), 0.0);
    transform(carrier_out_i.begin(), carrier_out_i.end(), rds_channel.begin(), mixed_data_i.begin(), std::multiplies<float>());
    transform(carrier_out_q.begin(), carrier_out_q.end(), rds_channel.begin(), mixed_data_q.begin(), std::multiplies<float>());

    // Low-pass Filter Fc=3KHz
    std::vector<float> after_3k_i, after_3k_q;
    downConvolveFIR(after_3k_i, mixed_data_i, LP3k_coeff, temp_after_3k_i, 1);
    downConvolveFIR(after_3k_q, mixed_data_q, LP3k_coeff, temp_after_3k_q, 1);

    // Rational resampler
    int up_decim = 19;
    int down_decim = 80;
    std::vector<float> upsampled_data_i(up_decim * after_3k_i.size(), 0.0);
    std::vector<float> upsampled_data_q(up_decim * after_3k_q.size(), 0.0);
    // upsample by 19
    for (auto i = 0; i < after_3k_i.size(); i++)
    {
        upsampled_data_i[i * up_decim] = after_3k_i[i];
        upsampled_data_q[i * up_decim] = after_3k_q[i];
    }
    // Low-pass Filter Fc=16KHz downsampled by 80 
    std::vector<float> resampled_data_i, resampled_data_q;
    upDownConvolveFIR(resampled_data_i, upsampled_data_i, LP16k_coeff, temp_Demod_LPF_i, up_decim, down_decim);
    upDownConvolveFIR(resampled_data_q, upsampled_data_q, LP16k_coeff, temp_Demod_LPF_q, up_decim, down_decim);

    // Root-raised-cosine Filter
    std::vector<float> rrc_out_i, rrc_out_q;
    downConvolveFIR(rrc_out_i, resampled_data_i, rrc_coeff, temp_rrc_i, 1);
    downConvolveFIR(rrc_out_q, resampled_data_q, rrc_coeff, temp_rrc_q, 1);

    // uncomment the following code to plot rrc
    // if (block_id == 18)
    // {
    //     std::vector<float> sample;
    //     genIndexVector(sample, rrc_out_i.size());
    //     logVector("sample", sample, rrc_out_i);
    // }

    // Clock & data recovery output
    std::vector<float> i_recover(RRC_3block_i.size() / 24, 0.0);
    std::vector<float> q_recover(RRC_3block_q.size() / 24, 0.0);

    // Manchaster boolean output
    std::vector<int> after_man_bool_i(38, 0);

    // Decoding output
    std::vector<int> after_deco_i(38, 0);

    // Parity check matrix
    int Parity_check[26][10] = {{1, 0, 0, 0, 0, 0, 0, 0, 0, 0},
                                {0, 1, 0, 0, 0, 0, 0, 0, 0, 0},
                                {0, 0, 1, 0, 0, 0, 0, 0, 0, 0},
                                {0, 0, 0, 1, 0, 0, 0, 0, 0, 0},
                                {0, 0, 0, 0, 1, 0, 0, 0, 0, 0},
                                {0, 0, 0, 0, 0, 1, 0, 0, 0, 0},
                                {0, 0, 0, 0, 0, 0, 1, 0, 0, 0},
                                {0, 0, 0, 0, 0, 0, 0, 1, 0, 0},
                                {0, 0, 0, 0, 0, 0, 0, 0, 1, 0},
                                {0, 0, 0, 0, 0, 0, 0, 0, 0, 1},
                                {1, 0, 1, 1, 0, 1, 1, 1, 0, 0},
                                {0, 1, 0, 1, 1, 0, 1, 1, 1, 0},
                                {0, 0, 1, 0, 1, 1, 0, 1, 1, 1},
                                {1, 0, 1, 0, 0, 0, 0, 1, 1, 1},
                                {1, 1, 1, 0, 0, 1, 1, 1, 1, 1},
                                {1, 1, 0, 0, 0, 1, 0, 0, 1, 1},
                                {1, 1, 0, 1, 0, 1, 0, 1, 0, 1},
                                {1, 1, 0, 1, 1, 1, 0, 1, 1, 0},
                                {0, 1, 1, 0, 1, 1, 1, 0, 1, 1},
                                {1, 0, 0, 0, 0, 0, 0, 0, 0, 1},
                                {1, 1, 1, 1, 0, 1, 1, 1, 0, 0},
                                {0, 1, 1, 1, 1, 0, 1, 1, 1, 0},
                                {0, 0, 1, 1, 1, 1, 0, 1, 1, 1},
                                {1, 0, 1, 0, 1, 0, 0, 1, 1, 1},
                                {1, 1, 1, 0, 0, 0, 1, 1, 1, 1},
                                {1, 1, 0, 0, 0, 1, 1, 0, 1, 1}};
    // Offset word
    std::vector<int> multi_result(10, 0);
    //  Syndrome word
    std::vector<int> syn_A({1, 1, 1, 1, 0, 1, 1, 0, 0, 0});
    std::vector<int> syn_B({1, 1, 1, 1, 0, 1, 0, 1, 0, 0});
    std::vector<int> syn_C({1, 0, 0, 1, 0, 1, 1, 1, 0, 0});
    std::vector<int> syn_Cprime({1, 1, 1, 1, 0, 0, 1, 1, 0, 0});
    std::vector<int> syn_D({1, 0, 0, 1, 0, 1, 1, 0, 0, 0});

    if (block_id != 0)
    {
        // collect 3 blocks (38 bits)
        if (count_RRC < 3)
        {
            for (auto i = 608 * count_RRC; i < (608 * (count_RRC + 1)); i++)
            {
                RRC_3block_i[i] = rrc_out_i[i % 608];
                RRC_3block_q[i] = rrc_out_q[i % 608];
            }
            // Clock recovery, find the best phase index that having least HHLLs
            if (block_id == 1)
            {
                int error_count_min = 0;
                for (auto i = 0; i < 24; i++) // recovery for every 24 signal points
                {
                    int error_count = 0;
                    for (auto j = 0; j < 25; j += 2) // 25 data points in one block
                    {
                        if ((RRC_3block_i[j * 24 + i] * RRC_3block_i[(j + 1) * 24 + i]) >= 0) // HHLL
                            error_count++;
                    }
                    if (i == 0)
                        error_count_min = error_count;
                    else
                    {
                        if (error_count < error_count_min)
                        {
                            error_count_min = error_count;
                            phase_index_i = i;
                            phase_index_q = i;
                        }
                    }
                }
            }
        }
        count_RRC++;

        // process every 3 blocks
        if (count_RRC == 3)
        {
            count_RRC = 0;
            // Data recovery
            for (int i = phase_index_i; i < RRC_3block_i.size(); i += 24)
                i_recover[(i - phase_index_i) / 24] = RRC_3block_i[i];

            for (int i = phase_index_q; i < RRC_3block_q.size(); i += 24)
                q_recover[(i - phase_index_i) / 24] = RRC_3block_q[i];

            // uncomment the following code to plot a constellation diagram
            // if (block_id == 18)
            // {
            //     logVector("rrc_output", i_recover, q_recover);
            // }

            // --- RDS Data Processing ---
            //Manchester
            for (int i = 0; i < i_recover.size(); i += 2)
            {
                if (i_recover[i] < i_recover[i + 1])
                    after_man_bool_i[i / 2] = 0;
                else
                    after_man_bool_i[i / 2] = 1;
            }

            //Differential decoding
            // first bit buffered
            if (deco_buf_i == after_man_bool_i[0])
                after_deco_i[0] = 0;
            else
                after_deco_i[0] = 1;
            // common case bits
            for (int i = 0; i < after_man_bool_i.size() - 1; i++)
            {
                if (after_man_bool_i[i] == after_man_bool_i[i + 1])
                    after_deco_i[i + 1] = 0;
                else
                    after_deco_i[i + 1] = 1;
            }
            deco_buf_i = after_man_bool_i[-1];

            // extend the RDS bitstream container after each RDS block
            multi_resize_count++;
            multi_prepared.resize(after_deco_i.size() * multi_resize_count);
            for (int i = after_deco_i.size() * (multi_resize_count - 1); i < multi_prepared.size(); i++)
                multi_prepared[i] = after_deco_i[i % after_deco_i.size()];

            // Frame Synchronization
            while ((multi_prepared.size() - syn_base) >= 26)
            {
                // matrix multiplication
                for (int i = 0; i < 10; i++)
                {
                    for (int j = 0; j < 26; j++)
                    {
                        multi_result[i] += multi_prepared[syn_base + j] * Parity_check[j][i];
                    }
                    multi_result[i] = multi_result[i] % 2;
                }

                // error detection and syndromes computation
                if (multi_result == syn_A)
                {
                    if ((syn_base - syn_base_buf) / 26 == 0)
                    {
                        std::cerr << "Bit pos " << syn_base << " Block type "
                                  << "A" << std::endl;
                        syn_base_buf = syn_base;
                    }
                    else
                        std::cerr << " Probable false positive with A syndrome starting at position " << syn_base << std::endl;
                }
                else if (multi_result == syn_B)
                {
                    if ((syn_base - syn_base_buf) / 26 == 0)
                    {
                        std::cerr << "Bit pos " << syn_base << " Block type "
                                  << "B" << std::endl;
                        syn_base_buf = syn_base;
                    }
                    else
                        std::cerr << " Probable false positive with B syndrome starting at position " << syn_base << std::endl;
                }
                else if (multi_result == syn_C)
                {
                    if ((syn_base - syn_base_buf) / 26 == 0)
                    {
                        std::cerr << "Bit pos " << syn_base << " Block type "
                                  << "C" << std::endl;
                        syn_base_buf = syn_base;
                    }
                    else
                        std::cerr << " Probable false positive with C syndrome starting at position " << syn_base << std::endl;
                }
                else if (multi_result == syn_Cprime)
                {
                    if ((syn_base - syn_base_buf) / 26 == 0)
                    {
                        std::cerr << "Bit pos " << syn_base << " Block type "
                                  << "C'" << std::endl;
                        syn_base_buf = syn_base;
                    }
                    else
                        std::cerr << " Probable false positive with C' syndrome starting at position " << syn_base << std::endl;
                }
                else if (multi_result == syn_D)
                {
                    if ((syn_base - syn_base_buf) / 26 == 0)
                    {
                        std::cerr << "Bit pos " << syn_base << " Block type "
                                  << "D" << std::endl;
                        syn_base_buf = syn_base;
                    }
                    else
                        std::cerr << " Probable false positive with D syndrome starting at position " << syn_base << std::endl;
                }
                syn_base++;
            }
        }
    }
    path_cvar.notify_one();
}