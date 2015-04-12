#ifndef FRONTEND_H
#define FRONTEND_H

#include "input.h"
#include "step.h"
#include "utils.h"

#include <fftw3.h>
#include <set>
#include <vector>

class FrontEnd: public Step
{
private:
    std::unordered_map<int,ffast_complex> signal;

    ffast_complex** observationMatrix;
    std::set<int> usedSamples; // A set to keep track of samples used by FFAST
    int signalLength;
    const Input* input;

    bool countSamplesDone;  // Flag to indicate if the number of samples used by FFAST are counted.
    ffast_real* samplingPeriods;
    ffast_complex* subSampledSignal;
    ffast_complex* DFTresults;
    fftw_plan* plans;
    std::vector<int> delays;
    ffast_real window(int i);

public:
    FrontEnd(Chrono* newChrono, const Config* newConfig, const Input* newInput);
    ~FrontEnd();
    void process();
    ffast_complex** getObservationMatrix() const;
    int getUsedSamplesNb() const;
    const std::vector<int> getDelays() const;
    const std::unordered_map<int,ffast_complex>& getTimeSignal() const;

private:
    void computeDelays();
};

#endif // FRONTEND_H
