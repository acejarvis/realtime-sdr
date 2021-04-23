/*
Comp Eng 3DY4 (Computer Systems Integration Project)

Department of Electrical and Computer Engineering
McMaster University
Ontario, Canada
*/

#ifndef DY4_RFFRONTEND_H
#define DY4_RFFRONTEND_H

// add headers as needed
#include <vector>

// declaration of a function prototypes
void readStdinBlockData(unsigned int num_samples, std::vector<float> &block_data);
void rfFrontend(std::vector<float> &fm_dem, std::vector<float> &prev_IQ, std::vector<float> &rf_coeff, std::vector<float> &temp_I, std::vector<float> &temp_Q, int rf_decim, int BLOCK_SIZE);

#endif // DY4_RFFRONTEND_H