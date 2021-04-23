/*
Comp Eng 3DY4 (Computer Systems Integration Project)

Department of Electrical and Computer Engineering
McMaster University
Ontario, Canada
*/

#ifndef DY4_RDS_H
#define DY4_RDS_H

// add headers as needed
#include <vector>
#include <math.h>
#include <chrono>
#include <queue>
#include <mutex>
#include <condition_variable>

// declaration of a function prototypes
void rdsPath(std::condition_variable &path_cvar,
             std::vector<float> &fm_dem,
             float rf_Fs, float &integrator, float &phaseEst, float &feedbackI, float &feedbackQ, float &nco_buf_i, float &nco_buf_q, int &trigOffset,
             std::vector<float> &rds1_coeff, std::vector<float> &rds2_coeff, std::vector<float> &rrc_coeff, std::vector<float> &LP3k_coeff, std::vector<float> &LP16k_coeff,
             std::vector<float> &temp_after_3k_i, std::vector<float> &temp_after_3k_q, std::vector<float> &temp_channel_BPF, std::vector<float> &temp_CR_BPF, std::vector<float> &temp_Demod_LPF_i, std::vector<float> &temp_Demod_LPF_q, 
             std::vector<float> &temp_rrc_i, std::vector<float> &temp_rrc_q, std::vector<float> &RRC_3block_i, std::vector<float> &RRC_3block_q,
             int block_id, int &count_RRC, int &phase_index_i, int &phase_index_q, int &deco_buf_i, int &multi_resize_count, std::vector<int> &multi_prepared, int &syn_base, int &syn_base_buf);

#endif // DY4_RDS_H