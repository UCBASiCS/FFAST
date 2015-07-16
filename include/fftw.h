#ifndef FFTW_H
#define FFTW_H

#include "input.h"
#include "step.h"
#include "utils.h"

#include <fftw3.h>

class FFTW: public Step
{
private:
    ffast_complex* signal;
    
    ffast_complex* outputSignal;

    ffast_real FFTWscalingFactor;
    fftw_plan plan;

public:
    FFTW(Chrono* newChrono, const Config* newConfig, ffast_complex *newSignal);
    ~FFTW();
    void process();
};

#endif // FFTW_H
