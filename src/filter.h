/*
Comp Eng 3DY4 (Computer Systems Integration Project)

Department of Electrical and Computer Engineering
McMaster University
Ontario, Canada
*/

#ifndef DY4_FILTER_H
#define DY4_FILTER_H

// add headers as needed
#include <iostream>
#include <vector>

// declaration of a function prototypes
// function to compute the impulse response "h" for a Low-pass filter
void impulseResponseLPF(float Fs, float Fc, unsigned short int num_taps, std::vector<float> &h);
// function to compute the impulse response "h" for a Band-pass filter
void impulseResponseBPF(float Fs, float Fe, float Fb, unsigned short int num_taps, std::vector<float> &h);
// function to compute the impulse response "h" for a Root-raised Cosine filter
void impulseResponseRRC(float Fs, unsigned short int num_taps, std::vector<float> &h);
// function to process FM Demodulation without arctan function (computational friendly)
void fmDemodNoArctan(const std::vector<float> &I, const std::vector<float> &Q, std::vector<float> &prev_IQ, std::vector<float> &fm_demod);
// function to filter signal with upsampling and downsampling
void upDownConvolveFIR(std::vector<float> &y, const std::vector<float> &x, const std::vector<float> &h, std::vector<float> &temp, int up_decim, int down_decim);
// function to filter signal with downsampling only
void downConvolveFIR(std::vector<float> &y, const std::vector<float> &x, const std::vector<float> &h, std::vector<float> &temp, int decim);
// function for phase-locked loop with mixed output
void fmPLL(std::vector<float> &pllIn, std::vector<float> &ncoOut, float freq, float Fs,
           float &integrator, float &phaseEst, float &feedbackI, float &feedbackQ, float &nco_buf, int &trigOffset,
           float ncoScale, float phaseAdjust, float normBandwidth);
// function for phase-locked loop with I/Q output
void fmPLL(std::vector<float> &pllIn, std::vector<float> &ncoOut_i, std::vector<float> &ncoOut_q, float freq, float Fs,
           float &integrator, float &phaseEst, float &feedbackI, float &feedbackQ, float &nco_buf_i, float &nco_buf_q, int &trigOffset,
           float ncoScale, float phaseAdjust, float normBandwidth);

#endif // DY4_FILTER_H
