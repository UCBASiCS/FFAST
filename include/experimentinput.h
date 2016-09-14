#ifndef EXPERIMENTINPUT_H
#define EXPERIMENTINPUT_H

#include "input.h"
#include "utils.h"

#include <fftw3.h>
#include <set>
#include <unordered_map>
#include <vector>

class ExperimentInput: public Input
{
private:
    std::unordered_map<int,ffast_complex> timeSignal;
    std::set<int> neededSamples;

    // for the experiment mode
    ffast_complex* frequencySignal;
    std::unordered_map<int,ffast_complex> nonZeroFrequencies;
    ffast_real signalMagnitude;
    ffast_real realSNR;
    ffast_real noiseStdDeviation;
    fftw_plan plan;

public:
    ExperimentInput(Chrono* newChrono, const Config* newConfig);
    ~ExperimentInput();
    void process();
    void process(std::vector<int> delays);
    const std::unordered_map<int,ffast_complex>& getTimeSignal() const;

    // for the experiment mode
    ffast_complex* getFrequencySignal() const;
    const std::unordered_map<int,ffast_complex>& getNonZeroFrequencies() const;
    ffast_real getSignalMagnitude() const;
    ffast_real getRealSNR() const;
    ffast_real getRealSNRdB() const;
    ffast_complex getSignalAtIndex(int k) const;

private:
    void findNeededSamples(std::vector<int> delays);

    // for the experiment mode
    ffast_real distribution(ffast_real _urand, std::vector<double> F);
    void generateNonZeroFrequencies();
    ffast_real getRandomPhase() const;
    void computeFrequencySignal();
    void addNoise();
    void frequencyToTime();
    void frequencyToTimeUsingFFT(std::vector<int> delays);
    void applyQuantization(int bitsNb);
};

#endif // EXPERIMENTINPUT_H
