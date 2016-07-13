#include "frontend.h"

#include <iostream>
#include <math.h>
#include <stdlib.h>

FrontEnd::FrontEnd(Chrono* newChrono, const Config* newConfig, const Input* newInput): Step(newChrono, newConfig),
    input(newInput), countSamplesDone(false)
{
    signalLength = config->getSignalLength();
    observationMatrix = (ffast_complex**) malloc(config->getBinsSum() * sizeof(ffast_complex*));

    for (int binAbsoluteIndex=0; binAbsoluteIndex<config->getBinsSum(); binAbsoluteIndex++)
    {
        observationMatrix[binAbsoluteIndex] = (ffast_complex*) malloc(config->getDelaysNb() * sizeof(ffast_complex));
    }

    // We compute short DFTs in place. Hence only one common malloc is needed for all the stages.
    subSampledSignal = (ffast_complex*) fftw_malloc(config->getBiggestBin() * sizeof(ffast_complex));
    DFTresults = (ffast_complex*) fftw_malloc(config->getBiggestBin() * sizeof(ffast_complex));

    samplingPeriods = (ffast_real*) malloc(config->getBinsNb() * sizeof(ffast_real));
    plans = (fftw_plan*) malloc(config->getBinsNb() * sizeof(fftw_plan));

    for (int stage=0; stage<config->getBinsNb(); stage++)
    {
        samplingPeriods[stage] = signalLength/config->getBinSize(stage);
        plans[stage] = fftw_plan_dft_1d(config->getBinSize(stage), reinterpret_cast<fftw_complex*>(subSampledSignal),
                    reinterpret_cast<fftw_complex*>(DFTresults), FFTW_FORWARD, config->getFFTWstrategy());
    }

    computeDelays();
}

FrontEnd::~FrontEnd()
{
    for (int binAbsoluteIndex=0; binAbsoluteIndex<config->getBinsSum(); binAbsoluteIndex++)
    {
        free(observationMatrix[binAbsoluteIndex]);
    }

    for (int stage=0; stage<config->getBinsNb(); stage++)
    {
        fftw_destroy_plan(plans[stage]);
    }

    free(plans);
    free(observationMatrix);
    free(subSampledSignal);
    free(DFTresults);
    free(samplingPeriods);
}

void FrontEnd::process()
{
    chrono->start("FrontEnd");

    int sampleIndex;
    int binAbsoluteIndex;
    signal = input->getTimeSignal();
    ffast_real sqrtStageSamplingPeriod;

    for (int stage = 0; stage < config->getBinsNb(); stage++)
    {
        int delayIndex = 0;
        sqrtStageSamplingPeriod = sqrt(samplingPeriods[stage]);

        for(auto delayIterator=delays.begin(); delayIterator != delays.end(); ++delayIterator)
        {
            // subsample the signal
            for (int binRelativeIndex=0; binRelativeIndex<config->getBinSize(stage); binRelativeIndex++)
            {
                // samples at the subsampled signal
                sampleIndex = ((int) (*delayIterator + binRelativeIndex*samplingPeriods[stage])) % signalLength;

                // Just counting the number of samples once if the -c flag is given
                if (!countSamplesDone && config->needToCountSamples())
                {
                    usedSamples.insert(sampleIndex);
                }
                // subSampledSignal[binRelativeIndex] = sqrt(samplingPeriods[stage])*input->getSignalAtIndex(sampleIndex)*window(sampleIndex);
                subSampledSignal[binRelativeIndex] = sqrtStageSamplingPeriod*input->getSignalAtIndex(sampleIndex)*window(sampleIndex);
            }

            fftw_execute(plans[stage]);

            binAbsoluteIndex = config->getBinOffset(stage);

            // occupy the observation matrix
            for (int binRelativeIndex=0; binRelativeIndex<config->getBinSize(stage); binRelativeIndex++)
            {
                observationMatrix[binAbsoluteIndex][delayIndex] = DFTresults[binRelativeIndex]/sqrt(config->getBinSize(stage));

                binAbsoluteIndex++;
            }

            delayIndex++;
        }
    }

    countSamplesDone = true;

    chrono->stop("FrontEnd");
}

ffast_complex** FrontEnd::getObservationMatrix() const
{
    return observationMatrix;
}

int FrontEnd::getUsedSamplesNb() const
{
    return usedSamples.size();
}

const std::vector<int> FrontEnd::getDelays() const
{
    return delays;
}

void FrontEnd::computeDelays()
{
    if (config->isNoisy() || config->applyWindow())
    {
        std::set<int> tempDelaysSet;
        if (config->needToUseMaximumLikelihoodDetection())
        {
            // We use ML detection only when the delay samples are randomly placed.
            while((int) tempDelaysSet.size() < config->getDelaysNb())
            {
                int tempDelay = ((int) floor(signalLength*((ffast_real) drand48()))) % signalLength;

                tempDelaysSet.insert(tempDelay);
            }
            delays.resize(config->getDelaysNb());
            std::copy(tempDelaysSet.begin(), tempDelaysSet.end(), delays.begin());
        }
        else {
            // We use FFAST_Search
	    int r = 0;

	    unsigned long long int skip = 1;
            for (int i = 0; i < config->getChainsNb(); ++i)
            {
		r = ((int) floor(signalLength*((ffast_real) drand48()))) % signalLength;
                for (int j = 0; j < config->getDelaysPerBunchNb(); ++j)
                {
                    delays.push_back(positiveMod(r+j*skip, signalLength));
                }
                skip *= 2; // 2^i
            }
        }
    }
    else
    {
        for (int delay=0; delay<config->getDelaysNb(); delay++)
        {
            delays.push_back(delay);
        }
    }
}

ffast_real FrontEnd::window(int i) {
    if(config->applyWindow())
    {
        // Blackman-Nuttall window
        double a0 = 0.3635819;
        double a1 = 0.4891775;
        double a2 = 0.1365995;
        double a3 = 0.0106411;
        return (ffast_real) 2 * ( a0
			          - a1 * cos((2*M_PI*i) / (signalLength-1))
		                  + a2 * cos((4*M_PI*i) / (signalLength-1))
			          - a3 * cos((6*M_PI*i) / (signalLength-1)) );
    }
    else
    {
        return 1;
    }
}

const std::unordered_map<int,ffast_complex>& FrontEnd::getTimeSignal() const
{
    return input->getTimeSignal();
}
