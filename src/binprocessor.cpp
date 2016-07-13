#include "binprocessor.h"

#include <algorithm>
#include <cmath>
#include <iostream>
#include <limits>
#include <math.h>
#include <stdlib.h>

BinProcessor::BinProcessor(Chrono *newChrono, const Config* newConfig, ffast_complex** newObservationMatrix, const std::vector<int>& newDelays):
    Step(newChrono, newConfig), observationMatrix(newObservationMatrix), delays(newDelays)
{
    signalLength     = config->getSignalLength();
    delaysNb         = config->getDelaysNb();
    chainsNb         = config->getChainsNb();
    delaysPerBunchNb = config->getDelaysPerBunchNb();

    MLdetection         = config->needToUseMaximumLikelihoodDetection();
    twoPiBySignalLength = (2*M_PI)/signalLength;
    signalLengthByTwoPi = 1/twoPiBySignalLength;

    signalVector    = (ffast_complex*) malloc(delaysNb * sizeof(ffast_complex));
    thresholds      = (ffast_real*) malloc(config->getBinsNb() * sizeof(ffast_real));
    directionVector = (ffast_complex*) malloc(delaysNb * sizeof(ffast_complex));

    computeThresholds();

    if (!MLdetection)
    {
        angles  = (ffast_real*) malloc((delaysPerBunchNb-1) * sizeof(ffast_real));
        weights = (ffast_real*) malloc((delaysPerBunchNb-1) * sizeof(ffast_real));
        computeWeights();
    }
}

BinProcessor::~BinProcessor()
{
    free(signalVector);
    free(thresholds);
    free(directionVector);

    if (!MLdetection)
    {
        free(angles);
        free(weights);
    }
}

void BinProcessor::process()
{
    if (MLdetection)
    {
        // MLprocess tries out all possible locations and estimates bin signal
        MLprocess();
    }
    else
    {
        computeLocation();
        estimateBinSignal();
    }
}

void BinProcessor::adjustTo(int newBinAbsoluteIndex, int newBinRelativeIndex, int newStage)
{
    binAbsoluteIndex = newBinAbsoluteIndex;
    binRelativeIndex = newBinRelativeIndex;
    stage   = newStage;
    binSize = config->getBinSize(stage);
}

bool BinProcessor::isSingleton()
{
    // a singleton is not zero-ton
    if (isZeroTon())
    {
        // The bin is a 'zero-ton'
        return false;
    }

    process();
    
    // if the bin is a singleton, when we peel the frequency from the bin
    // the remaining power should be coming from noise only
    // also we need the singleton to have energy larger than a threshold (minimumEnergy)
    if (noise <= thresholds[stage] && std::norm(amplitude) > minimumEnergy)
    {
        /* for checking the thresholds
        std::cout << "---------" << std::endl;
        std::cout << "remaining noise " << noise << " threshold: " << thresholds[stage] << std::endl;
        std::cout << "signal energy " << std::norm(amplitude) << " minimum energy: " << minimumEnergy << std::endl;
        */
        // The bin is a singleton
        return true;
    }
    /* for checking the thresholds
    else
    {
        std::cout << "---------" << std::endl;
        std::cout << "remaining noise " << noise << " threshold: " << thresholds[stage] << std::endl;
        std::cout << "signal energy " << std::norm(amplitude) << " minimum energy: " << minimumEnergy << std::endl;
    }
    */

    // The bin is a 'multi-ton'
    return false;
}

int BinProcessor::getLocation() const
{
    return location;
}

ffast_complex BinProcessor::getAmplitude() const
{
    return amplitude;
}

ffast_complex BinProcessor::getSignal(int delayIndex) const
{
    return signalVector[delayIndex];
}

bool BinProcessor::isZeroTon() const
{
    ffast_real energy = 0;

    for (int delayIndex=0; delayIndex<delaysNb; delayIndex++)
    {
        energy += std::norm(observationMatrix[binAbsoluteIndex][delayIndex]);
    }

    if (energy <= thresholds[stage])
    {
        return true;
    }

    return false;
}

void BinProcessor::MLprocess()
{
    ffast_real MLnoise = std::numeric_limits<ffast_real>::infinity();
    int MLlocation = 0;

    /*
        The balls (spectrum of the original signal) that can be in this bin are
        at indices binRelativeIndex + n*numberOfBinsAtStage.
        Here, we find the best possible ball that matches the signal at this bin.
    */
    location = binRelativeIndex;

    while (location < signalLength)
    {
        estimateBinSignal();

        if (noise < MLnoise)
        {
            MLnoise = noise;
            MLlocation = location;
        }

        // binSize is equal to the number of bins at the stage
        location += binSize;
    }

    location = MLlocation;
    estimateBinSignal();
}

void BinProcessor::estimateBinSignal()
{
    int delayIndex = 0;
    amplitude      = 0;

    for(auto delayIterator=delays.begin(); delayIterator != delays.end(); ++delayIterator)
    {
        // phase picked up by a tone at 'location' due to the delay
        directionVector[delayIndex] = std::polar((ffast_real) 1, (ffast_real) twoPiBySignalLength*location*(*delayIterator));
        // remove the picked up phase
        amplitude += std::conj(directionVector[delayIndex]) * observationMatrix[binAbsoluteIndex][delayIndex];

        delayIndex++;
    }

    // average the found amplitudes
    amplitude /= delaysNb;

    noise = 0;

    for (int delayIndex=0; delayIndex<delaysNb; delayIndex++)
    {
        // compute the bin signal that would result from a tone with such amplitude and location
        signalVector[delayIndex] = amplitude*directionVector[delayIndex];
        // the residual between the actual bin measurement and the computed one is the noise
        noise += std::norm(observationMatrix[binAbsoluteIndex][delayIndex]-signalVector[delayIndex]);
    }
}

void BinProcessor::computeLocation()
{
    ffast_real location_bis = 0;
    ffast_real nwrap = 1;
    ffast_real tempLocation;
    ffast_real loc_update;
    ffast_real r;

    for (int i = 0; i < chainsNb; ++i)
    {
        // Kay's method
        tempLocation = fmod(signalLengthByTwoPi*getOmega(i), signalLength);
        tempLocation = tempLocation > 0 ? tempLocation : tempLocation+signalLength;

	// Update Location
        loc_update = tempLocation/nwrap - location_bis;
        r =  loc_update - round((loc_update * nwrap) / ((ffast_real) signalLength)) * (((ffast_real) signalLength)/nwrap);
        location_bis += r;
        nwrap *= 2; // equal to 2^i
    }

    location = positiveMod(
                (int) round((location_bis-binRelativeIndex)/binSize)*binSize + binRelativeIndex,
                signalLength
            );
}

// i corresponds to the delay-bunch index
ffast_real BinProcessor::getOmega(int i) const
{
    ffast_real omega = 0;
    bool needToShift = false;

    // needToShift flag is used to distinguish between the case when omega is close to 0 or close to Pi
    // These two points require different averaging process.
    for (int delayIndex=0; delayIndex<delaysPerBunchNb-1; delayIndex++)
    {
        angles[delayIndex] = std::arg(std::conj(observationMatrix[binAbsoluteIndex][i*delaysPerBunchNb + delayIndex])
                                  * observationMatrix[binAbsoluteIndex][i*delaysPerBunchNb + delayIndex+1]);

        if (!needToShift && angles[delayIndex] < -M_PI/2)
        {
           needToShift = true;
        }

        omega += weights[delayIndex] * angles[delayIndex];
    }

    if (needToShift)
    {
        for (int delayIndex=0; delayIndex<delaysPerBunchNb-1; delayIndex++)
        {
            if (angles[delayIndex] < 0)
            {
                omega += weights[delayIndex]*2*M_PI;
            }
        }
    }

    if (omega < 0)
    {
        omega += 2*M_PI;
    }

    return omega;
}

// for the fast search linear estimator
void BinProcessor::computeWeights()
{
    ffast_real baseWeight = 6.0 / ((ffast_real) delaysPerBunchNb*(delaysPerBunchNb*delaysPerBunchNb - 1.0));

    for (int delayIndex=0; delayIndex<delaysPerBunchNb-1; delayIndex++)
    {
        weights[delayIndex] = baseWeight*(delayIndex+1)*(delaysPerBunchNb-(delayIndex+1));
    }
}

void BinProcessor::computeThresholds()
{
    // for computing the histogram for zeroton-singleton-multiton test
    std::vector<double> energyBins;
    double noiseEstimation = 0;
    double tempEnergy = 0;
    int    energyHistogramBinsCounted = 0;

    if ( config->isNoisy() || config->applyWindow() )
    {
        // go over the stages
        for (int stage = 0; stage < config->getBinsNb(); ++stage)
        {
            // go over the bins in the stage
            for (int i = 0; i < config->getBinSize(stage); ++i)
            {
                // this is the energy of a bin measurement
                tempEnergy = 0;
                // go over the delays (go over bin measurements)
                for (int j = 0; j < config->getDelaysNb(); ++j)
                {
                    // std::norm calculates norm-squared
                    tempEnergy += std::norm(observationMatrix[config->getBinOffset(stage)+i][j]);
                }
                energyBins.push_back(tempEnergy);
            }
        }

        // make histogram below
        // sort the energy bins
        std::sort (energyBins.begin(), energyBins.end());

        double oneOverEta = (double)config->getSignalSparsityPeeling()/((double)config->getBinsSum()/(double)config->getBinsNb());
        double maxValueWanted = round(energyBins.size()/M_E);

        for (int i = 1; i <= 2; ++i)
        {
            maxValueWanted += energyBins.size()*exp(-oneOverEta)*pow(oneOverEta,i)/tgamma(i+1);
        }

        /*
        for (int i = 0; i < round(energyBins.size()*exp(-oneOverEta)); ++i)
        {
            std::cout << "energy bins" << i << ": " << energyBins[i] << std::endl;
            noiseEstimation += energyBins[i];
            energyHistogramBinsCounted++;
        }
        */

        bool noiseLevelCrossed = false;
        while (!noiseLevelCrossed)
        {
            // std::cout << "energy bins" << energyHistogramBinsCounted << ": " << energyBins[energyHistogramBinsCounted] << std::endl;
            noiseEstimation += energyBins[energyHistogramBinsCounted];
            energyHistogramBinsCounted++;
            // if the bin energy is 10 times the previous one, declare it has signal in it
            if ( energyBins[energyHistogramBinsCounted] / energyBins[energyHistogramBinsCounted-1] >= 10 )
            {
                noiseLevelCrossed = true;
            }
        }


        // CHOICE 1: this corresponds to the average energy of the zero-ton bins
        noiseEstimation /= energyHistogramBinsCounted;

        // CHOICE 2: this corresponds to the maximum energy of the zero-ton bins
        // noiseEstimation = energyBins[ round( energyBins.size() / (M_E) ) ];
    }

    // Minimum energy for the signal to be accepted as non-zero
    // 10 percent of the minimum SNR

    if( !(config->isNoisy() || config->applyWindow()) ) // noiseless ongrid
    {
        minimumEnergy = pow(10,-8);
    }
    else if ( config->isNoisy() && !config->applyWindow() ) // noisy ongrid
    {
        /*
            We want a threshold on the minimum signal energy to eliminate false
            detections. A simple one would be of the sort:
            minimumEnergy = 0.1 * noiseEstimation * pow(10,config->getSNRdB()/10);
            However, for large signal to noise ratio, this is not working since there
            is always a noise floor around 1e-18 due to numerical issues, hence when
            multiplied with high SNR it blows up. Hence, for small SNR we do that, but
            for large SNR we clip it.
            We chose it to clip to 1000 times the noise floor.
        */
        minimumEnergy = std::min( 0.1 * noiseEstimation * pow(10,config->getSNRdB()/10) , 1000 * noiseEstimation );
    }

    if ( config->applyWindow() ) // offgrid
    {
        minimumEnergy = 0.1 * config->getMinFourierMagnitude();
    }

    // std::cout << "estimated noise: " << noiseEstimation << std::endl;

    /*
        Base threshold for the remaining bin signal to be considered as
        zero-ton after the signal is peeled from it. This value, in theory, is zero
        for noiseless simulations. However, due to machine precision, it needs to be positive
        value. For noiseless simulations we chose it to be 1e-13. If this value is high,
        we start to get false detections. One can decrease the thrsholds or increase the
        number of delays to reduce the false detections.
    */
    ffast_real baseThreshold = pow(10, -13);

    ffast_real factor;

    // Thresholds are obtained using tailbounds of Gaussian random variable
    if (delaysNb < 10)
    {
        factor = 4;
    }
    else if (delaysNb < 20)
    {
        factor = 3;
    }
    else if (delaysNb < 50)
    {
        factor = 2;
    }
    else
    {
        factor = 1.5;
    }

    for (int stage=0; stage<config->getBinsNb(); stage++)
    {
        thresholds[stage] = baseThreshold;

        if ( config->applyWindow() || config->isNoisy() )
        {
            thresholds[stage] = pow(10,-10) + (ffast_real) factor * noiseEstimation;
        }
    }

}
