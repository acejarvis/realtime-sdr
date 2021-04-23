/*
Comp Eng 3DY4 (Computer Systems Integration Project)

Department of Electrical and Computer Engineering
McMaster University
Ontario, Canada
*/

// source code for Fourier-family of functions
#include "dy4.h"
#include "fourier.h"

// just DFT function (no FFT yet)
void DFT(const std::vector<float> &x, std::vector<std::complex<float>> &Xf)
{
    Xf.resize(x.size(), static_cast<std::complex<float>>(0, 0));
    for (auto m = 0; m < Xf.size(); m++)
    {
        for (auto k = 0; k < x.size(); k++)
        {
            std::complex<float> expval(0, -2 * PI * (k * m) / x.size());
            Xf[m] += x[k] * std::exp(expval);
        }
    }
}

// function to compute the magnitude values in a complex vector
void computeVectorMagnitude(const std::vector<std::complex<float>> &Xf, std::vector<float> &Xmag)
{
    // only the positive frequencies
    Xmag.resize(Xf.size(), static_cast<float>(0));
    for (auto i = 0; i < Xf.size(); i++)
    {
        Xmag[i] = std::abs(Xf[i]) / Xf.size();
    }
}

// add your own code to estimate the PSD
void estimatePSD(const std::vector<float> &samples, std::vector<float> &psd_est, float Fs, std::vector<float> &freq)
{
    int freq_bins = NFFT;
    float df = Fs / freq_bins;
    freq.resize(freq_bins / 2, static_cast<float>(0));
    // create x array for plot
    for (auto i = 0; i < freq_bins / 2; i++)
    {
        freq[i] = i * df;
    }

    std::vector<float> hann;
    hann.resize(freq_bins, static_cast<float>(0));
    for (auto i = 0; i < hann.size(); i++)
    {
        hann[i] = pow(sin(i * PI / freq_bins), 2);
    }
    int no_segments = floor(samples.size() / freq_bins);

    std::vector<float> psd_list;
    psd_list.resize(freq_bins * no_segments, static_cast<float>(0));

    for (auto k = 0; k < no_segments; k++)
    {
        std::vector<float> windowed_samples;
        windowed_samples.resize(freq_bins, static_cast<float>(0));
        std::vector<float> psd_seg;
        psd_seg.resize(freq_bins / 2, static_cast<float>(0));
        std::vector<std::complex<float>> Xf;
        Xf.resize(freq_bins, static_cast<std::complex<float>>(0, 0));
        int j = 0;
        for (auto m = k * freq_bins; m < (k + 1) * freq_bins; m++)
        {
            windowed_samples[j] = samples[m] * hann[j];
            j++;
        }
        DFT(windowed_samples, Xf);

        for (auto n = 0; n < freq_bins / 2; n++)
        {
            psd_seg[n] = 2 * (1 / (Fs * freq_bins / 2)) * pow(abs(Xf[n]), 2);
            psd_seg[n] = 10 * log10(psd_seg[n]);
            psd_list[n + k * freq_bins / 2] = psd_seg[n];
        }
    }

    psd_est.resize(freq_bins / 2, static_cast<float>(0));

    for (auto j = 0; j < freq_bins / 2; j++)
    {
        for (auto l = 0; l < no_segments; l++)
        {
            psd_est[j] += psd_list[j + l * freq_bins / 2];
        }
        psd_est[j] /= no_segments;
    }
}