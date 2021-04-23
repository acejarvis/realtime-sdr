/*
Comp Eng 3DY4 (Computer Systems Integration Project)

Department of Electrical and Computer Engineering
McMaster University
Ontario, Canada
*/

#include "dy4.h"
#include "filter.h"
#include <iostream>
#include <fstream>
#include <vector>
#include <iomanip>
#include <complex>
#include <cmath>
#include <chrono>

// function to compute the impulse response "h" for a Low-pass filter
void impulseResponseLPF(float Fs, float Fc, unsigned short int num_taps, std::vector<float> &h)
{
    // bring your own functionality
    // allocate memory for the impulse response
    h.resize(num_taps, 0.0);

    float norm_cutoff = Fc / Fs * 2;

    for (auto i = 0; i < num_taps; i++)
    {
        if (i == (num_taps - 1) / 2)
        {
            h[i] = norm_cutoff;
        }
        else
        {
            h[i] = norm_cutoff * std::sin(PI * norm_cutoff * (i - (num_taps - 1) / 2)) / (PI * norm_cutoff * (i - (num_taps - 1) / 2));
        }

        h[i] *= pow(std::sin(i * PI / num_taps), 2.0);
    }
}

// function to compute the impulse response "h" for a Band-pass filter
void impulseResponseBPF(float Fs, float Fe, float Fb, unsigned short int num_taps, std::vector<float> &h)
{
    // bring your own functionality
    // allocate memory for the impulse response
    h.resize(num_taps, 0.0);

    float norm_center = ((Fe + Fb) / 2) / (Fs / 2);
    float norm_pass = (Fe - Fb) / (Fs / 2);

    for (auto i = 0; i < num_taps; i++)
    {
        if (i == (num_taps - 1) / 2)
        {
            h[i] = norm_pass;
        }
        else
        {
            h[i] = norm_pass * std::sin(PI * (norm_pass / 2) * (i - (num_taps - 1) / 2)) / (PI * (norm_pass / 2) * (i - (num_taps - 1) / 2));
        }

        h[i] *= std::cos(i * PI * norm_center);
        h[i] *= pow(std::sin(i * PI / num_taps), 2.0);
    }
}

// function to compute the impulse response "h" for a Root-raised Cosine filter
void impulseResponseRRC(float Fs, unsigned short int num_taps, std::vector<float> &h)
{
    float T_symbol = 1 / 2375.0;
    float beta = 0.90;
    h.resize(num_taps, 0.0);
    for (auto i = 0; i < num_taps; i++)
    {
        float t = (i - num_taps / 2) / Fs;
        if (t == 0.0)
        {
            h[i] = 1.0 + beta * (4 / PI - 1);
        }
        else if (t == -T_symbol / (4 * beta) || t == T_symbol / (4 * beta))
        {
            h[i] = (beta / sqrt(2)) * (((1 + 2 / PI) *
                                        (std::sin(PI / (4 * beta)))) +
                                       ((1 - 2 / PI) * (std::cos(PI / (4 * beta)))));
        }
        else
        {
            h[i] = (std::sin(PI * t * (1 - beta) / T_symbol) +
                    4 * beta * (t / T_symbol) * std::cos(PI * t * (1 + beta) / T_symbol)) /
                   (PI * t * (1 - (4 * beta * t / T_symbol) * (4 * beta * t / T_symbol)) / T_symbol);
        }
    }
}

// function to process FM Demodulation without arctan function (computational friendly)
void fmDemodNoArctan(const std::vector<float> &I, const std::vector<float> &Q, std::vector<float> &prev_IQ, std::vector<float> &fm_demod)
{
    fm_demod.resize(I.size(), 0.0);
    std::vector<float> sI, sQ;
    sI.resize(I.size(), 0.0);
    sQ.resize(I.size(), 0.0);

    for (auto i = 0; i < sI.size(); i++)
    {
        sI[i] = I[i] * I[i];
        sQ[i] = Q[i] * Q[i];
    }

    std::vector<float> Inew, Qnew;
    Inew.resize(I.size() + 1, 0.0);
    Qnew.resize(Q.size() + 1, 0.0);
    for (auto j = 0; j < Inew.size(); j++)
    {
        if (j == 0)
        {
            Inew[j] = prev_IQ[0];
            Qnew[j] = prev_IQ[1];
        }
        else
        {
            Inew[j] = I[j - 1];
            Qnew[j] = Q[j - 1];
        }
    }

    float current_phase;
    for (auto k = 0; k < Inew.size() - 1; k++)
    {
        if (sI[k] + sQ[k] == 0)
        {
            current_phase = 0;
        }
        else
        {
            current_phase = (Inew[k + 1] * (Qnew[k + 1] - Qnew[k]) - Qnew[k + 1] * (Inew[k + 1] - Inew[k])) / (sI[k] + sQ[k]);
        }
        fm_demod[k] = current_phase;
        //store previous I,Q for block processing
        prev_IQ[0] = Inew[k + 1];
        prev_IQ[1] = Qnew[k + 1];
    }
}

// function to filter signal with upsampling and downsampling
void upDownConvolveFIR(std::vector<float> &y, const std::vector<float> &x, const std::vector<float> &h, std::vector<float> &temp, int up_decim, int down_decim)
{
    // bring your own functionality
    // allocate memory for the output (filtered) data
    y.resize(x.size() / down_decim, 0.0);
    int N_taps = h.size();
    std::vector<int> phase;
    if (up_decim == 24)
        phase = {23, 18, 13, 8, 3, 22, 17, 12, 7, 2, 21, 16, 11, 6, 1, 20, 15, 10, 5, 0, 19, 14, 9, 4};
    else
        phase = {18, 14, 10, 6, 2, 17, 13, 9, 5, 1, 16, 12, 8, 4, 0, 15, 11, 7, 3};

    while (true)
    {
        //k is used to write elements into the array
        int k = 0;
        //to locate the beginning of the array for each data
        int l_base = 0;
        int m = 0;
        int count = 0;
        // phase array for the first j to be calculated

        for (auto i = 0; i < x.size(); i += down_decim)
        {
            int l = (up_decim == 24) ? 6 : 12; // temp offset
            for (auto j = phase[(i / down_decim) % up_decim]; j < N_taps; j += 24)
            {
                if (l_base + l < N_taps - 1)
                {
                    y[i / down_decim] = y[i / down_decim] + temp[l_base + l] * h[N_taps - j - 1];
                    l += 24;
                }
                else
                {
                    if ((i + j - N_taps + 1) % up_decim == 0 && (i + j - N_taps + 1) != 0)
                    {
                        y[i / down_decim] = y[i / down_decim] + x[i + j + 1 - N_taps] * h[N_taps - j - 1];
                    }
                }
            }
            y[i / down_decim] *= up_decim;

            //the beginning of the array needs to increase after finish a y data
            l_base = l_base + down_decim;

            if (x.size() - count < N_taps)
            {
                temp[k] = x[count];
                k++;
            }
            m = i;
            count++;
        }
        if (m == x.size() - down_decim)
        {
            break;
        }
    }
}

// function to filter signal with downsampling only
void downConvolveFIR(std::vector<float> &y, const std::vector<float> &x, const std::vector<float> &h, std::vector<float> &temp, int decim)
{
    // bring your own functionality
    // allocate memory for the output (filtered) data
    y.resize(x.size() / decim, 0.0);
    int N_taps = h.size();

    while (true)
    {
        //k is used to write elements into the array
        int k = 0;
        //to locate the beginning of the array for each data
        int l_base = 0;
        int m = 0;
        int count = 0;

        for (auto i = 0; i < x.size(); i += decim)
        {

            int l = 0;

            for (auto j = 0; j < N_taps; j++)
            {
                if (l_base + l < N_taps - 1)
                {
                    y[i / decim] = y[i / decim] + temp[l_base + l] * h[N_taps - j - 1];
                    l++;
                }
                else
                {
                    y[i / decim] = y[i / decim] + x[i + j + 1 - N_taps] * h[N_taps - j - 1];
                }
            }
            //the beginning of the array needs to increase after finish a y data
            l_base = l_base + decim;

            if ((x.size() - count) < N_taps)
            {
                temp[k] = x[count];
                k++;
            }
            m = i;
            count++;
        }
        if (m == x.size() - decim)
        {
            break;
        }
    }
}

// function for phase-locked loop with mixed output
void fmPLL(std::vector<float> &pllIn, std::vector<float> &ncoOut, float freq, float Fs, float &integrator, float &phaseEst, float &feedbackI, float &feedbackQ, float &nco_buf, int &trigOffset, float ncoScale = 1.0, float phaseAdjust = 0.0, float normBandwidth = 0.01)
{
    //scale factors for proportional/integrator terms
    //these scale factors were derived assuming the following:
    //damping factor of 0.707 (1 over square root of 2)
    //there is no oscillator gain and no phase detector gain
    float Cp = 2.666;
    float Ci = 3.555;

    //gain for the proportional term
    float Kp = (normBandwidth)*Cp;
    //gain for the integrator term
    float Ki = (normBandwidth * normBandwidth) * Ci;
    ncoOut.resize(pllIn.size() + 1, 0.0);

    ncoOut[0] = nco_buf;

    float errorI, errorQ, errorD, triArg;

    for (auto i = 0; i < pllIn.size(); i++)
    {
        //phase detector
        errorI = pllIn[i] * feedbackI;        // complex conjugate of the
        errorQ = pllIn[i] * (-1) * feedbackQ; //feedback complex exponential
        errorD = std::atan2(errorQ, errorI);

        //loop filter
        integrator += Ki * errorD;
        //update phase estimate
        phaseEst += (Kp * errorD) + integrator;
        //internal oscillator
        triArg = 2 * PI * (freq / Fs) * (trigOffset + i + 1) + phaseEst;
        feedbackI = std::cos(triArg);
        feedbackQ = std::sin(triArg);
        ncoOut[i + 1] = std::cos(triArg * ncoScale + phaseAdjust);
    }
    trigOffset += pllIn.size();
    nco_buf = ncoOut[-1];
    ncoOut.pop_back();
}

// function for phase-locked loop with I/Q output
void fmPLL(std::vector<float> &pllIn, std::vector<float> &ncoOut_i, std::vector<float> &ncoOut_q, float freq, float Fs, float &integrator, float &phaseEst, float &feedbackI, float &feedbackQ, float &nco_buf_i, float &nco_buf_q, int &trigOffset, float ncoScale = 1.0, float phaseAdjust = 0.0, float normBandwidth = 0.01)
{
    float Cp = 2.666;
    float Ci = 3.555;

    //gain for the proportional term
    float Kp = (normBandwidth)*Cp;
    //gain for the integrator term
    float Ki = (normBandwidth * normBandwidth) * Ci;
    ncoOut_i.resize(pllIn.size() + 1, 0.0);
    ncoOut_q.resize(pllIn.size() + 1, 0.0);

    ncoOut_i[0] = nco_buf_i;
    ncoOut_q[0] = nco_buf_q;

    float errorI, errorQ, errorD, triArg;

    for (auto i = 0; i < pllIn.size(); i++)
    {
        //phase detector
        errorI = pllIn[i] * feedbackI;        // complex conjugate of the
        errorQ = pllIn[i] * (-1) * feedbackQ; //feedback complex exponential
        errorD = std::atan2(errorQ, errorI);
        //std::cerr << errorD << std::endl;

        //loop filter
        integrator += Ki * errorD;
        //update phase estimate
        phaseEst += (Kp * errorD) + integrator;
        //internal oscillator
        triArg = 2 * PI * (freq / Fs) * (trigOffset + i + 1) + phaseEst;
        feedbackI = std::cos(triArg);
        feedbackQ = std::sin(triArg);
        ncoOut_i[i + 1] = std::cos(triArg * ncoScale + phaseAdjust);
        ncoOut_q[i + 1] = std::sin(triArg * ncoScale + phaseAdjust);
    }
    trigOffset += pllIn.size();
    nco_buf_i = ncoOut_i[ncoOut_i.size() - 1];
    nco_buf_q = ncoOut_q[ncoOut_q.size() - 1];
    ncoOut_i.pop_back();
    ncoOut_q.pop_back();
}