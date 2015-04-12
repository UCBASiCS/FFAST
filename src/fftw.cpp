#include "fftw.h"

#include <iostream>

FFTW::FFTW(Chrono* newChrono, const Config* newConfig, ffast_complex* newSignal): Step(newChrono, newConfig), signal(newSignal)
{
    outputSignal = (ffast_complex*) fftw_malloc(config->getSignalLength() * sizeof(ffast_complex));

    FFTWscalingFactor = sqrt((ffast_real) config->getSignalLength());

    plan = fftw_plan_dft_1d(config->getSignalLength(), reinterpret_cast<fftw_complex*>(signal),
            reinterpret_cast<fftw_complex*>(outputSignal), FFTW_FORWARD, config->getFFTWstrategy());
}

FFTW::~FFTW()
{
    fftw_destroy_plan(plan);

    free(outputSignal);
}

void FFTW::process()
{
    chrono->start("FFTW");

    fftw_execute(plan);

    chrono->stop("FFTW");
}
