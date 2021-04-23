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

// read block data from stdin
void readStdinBlockData(unsigned int num_samples, std::vector<float> &block_data)
{
    std::vector<char> raw_data(num_samples);
    std::cin.read(reinterpret_cast<char *>(&raw_data[0]), num_samples * sizeof(char));
    for (unsigned int k = 0; k < num_samples; k++)
    {
        //automatically normalizes the data to the range -1 tp +1
        block_data[k] = float(((unsigned char)raw_data[k] - 128) / 128.0);
    }
}

// process RF Frontend and output to stdout
void rfFrontend(std::vector<float> &fm_dem, std::vector<float> &prev_IQ, std::vector<float> &rf_coeff, std::vector<float> &temp_I, std::vector<float> &temp_Q, int rf_decim, int BLOCK_SIZE)
{
    std::vector<float> block_data(BLOCK_SIZE);
    readStdinBlockData(BLOCK_SIZE, block_data);
    if ((std::cin.rdstate()) != 0)
    {
        std::cerr << "End of input stream reached" << std::endl;
        exit(1);
    }

    std::vector<float> I;
    I.resize(block_data.size() / 2, 0.0);
    std::vector<float> Q;
    Q.resize(block_data.size() / 2, 0.0);
    for (auto k = 0; k < I.size(); k++)
    {
        I[k] = block_data[k * 2];
        Q[k] = block_data[k * 2 + 1];
    }

    std::vector<float> I_filt;
    std::vector<float> Q_filt;

    downConvolveFIR(I_filt, I, rf_coeff, temp_I, rf_decim);
    downConvolveFIR(Q_filt, Q, rf_coeff, temp_Q, rf_decim);

    //FM demodulator
    fmDemodNoArctan(I_filt, Q_filt, prev_IQ, fm_dem);
}