#include "backend.h"

#include <algorithm>
#include <iostream>
#include <cmath>
#include <stdlib.h>
#include <stdio.h>
#include <fstream>

BackEnd::BackEnd(Chrono* newChrono, const Config* newConfig, const FrontEnd* newFrontEnd):
    Step(newChrono, newConfig), frontEnd(newFrontEnd)
{
    observationMatrix = frontEnd->getObservationMatrix();
    changed = (bool*) malloc(config->getBinsSum() * sizeof(bool));
}

BackEnd::~BackEnd()
{
    free(changed);
}

void BackEnd::process()
{
    chrono->start("BackEnd");
    binProcessor = new BinProcessor(chrono, config, observationMatrix, frontEnd->getDelays());

    int binAbsoluteIndex;
    bool singletonFound = true;
    initialize();

    while(singletonFound)
    {       
        singletonFound = false;

        for (int stage = 0; stage < config->getBinsNb(); stage++)
        {
            binAbsoluteIndex = config->getBinOffset(stage);

            for (int binRelativeIndex = 0; binRelativeIndex < config->getBinSize(stage); binRelativeIndex++)
            {
                // Copies the values of bin index and related variables in binProcessor object.
                binProcessor->adjustTo(binAbsoluteIndex, binRelativeIndex, stage);

                if (   
			changed[binAbsoluteIndex]
                        && binProcessor->isSingleton()
                        && decodedFrequencies.count(binProcessor->getLocation()) == 0
                   )
                {
                    singletonFound = true;
                    decodedFrequencies[binProcessor->getLocation()] = binProcessor->getAmplitude();
                    peelFrom(binProcessor->getLocation());
                }

                changed[binAbsoluteIndex] = false;

                if ((int) decodedFrequencies.size() == config->getSignalSparsityPeeling())
                {
                    goto stopPeeling;
                }

                binAbsoluteIndex++;
            }
        }
    }

    stopPeeling:

    if ( config->applyWindow() )
    {
        getClusteredFrequencies();
    }
    else
    {
        for (auto it = decodedFrequencies.cbegin(); it != decodedFrequencies.cend(); ++it)
        {
            realFrequenciesIndices.push_back(it->first);
        }
    }

    delete binProcessor;

    chrono->stop("BackEnd");
}

const std::vector<double>& BackEnd::getRealFrequenciesIndices() const
{
    return realFrequenciesIndices;
}

const std::unordered_map<int,ffast_complex>& BackEnd::getDecodedFrequencies() const
{
    return decodedFrequencies;
}

void BackEnd::swapDecodedFrequencies(std::unordered_map<int,ffast_complex> newDecodedFrequencies)
{
    decodedFrequencies.swap(newDecodedFrequencies);
}

void BackEnd::initialize()
{
    decodedFrequencies.clear();
    realFrequenciesIndices.clear();

    for (int binAbsoluteIndex = 0; binAbsoluteIndex < config->getBinsSum(); binAbsoluteIndex++)
    {
        changed[binAbsoluteIndex] = true;
    }
}

void BackEnd::peelFrom(int location)
{
    int hash;
    for (int stage = 0; stage < config->getBinsNb(); stage++)
    {
        hash = ( location % config->getBinSize(stage) ) + config->getBinOffset(stage);
        for (int delayIndex = 0; delayIndex < config->getDelaysNb(); delayIndex++)
        {
            observationMatrix[hash][delayIndex] -= binProcessor->getSignal(delayIndex);
        }

        changed[hash] = true;
    }
}

void BackEnd::getClusteredFrequencies()
{
    // DecodedFrequencies is an unordered_map hence the elements are
    // not ordered with respect to its key
    // We order them here
    std::map<int,ffast_complex> decodedFrequenciesOrdered;
    
    double energy;
    double totalEnergy = 0;
    double peak = 0;
    ffast_complex amplitude=0;

    // Order the decoded frequencies with respect to their locations
    for (auto it = decodedFrequencies.cbegin(); it != decodedFrequencies.cend(); ++it)
    {
        decodedFrequenciesOrdered.emplace(it->first,it->second);
    }

    // Clear the decoded frequencies map
    decodedFrequencies.clear();

    // Ratio of the original signal length to the truncated one
    double ratio = ((double) config->getSignalLengthOriginal())/((double) config->getSignalLength());
    
    int nb = config->getSignalSparsity() * 100;

    for (auto it = decodedFrequenciesOrdered.cbegin(); it != decodedFrequenciesOrdered.cend(); ++it)
    {
        energy = std::norm(it->second);
        peak  += energy * (it->first);
        totalEnergy += energy;
        amplitude   += (it->second) * energy;
	// if the cluster has a band of two zeros afterwards
        if( decodedFrequenciesOrdered.count(it->first+1) + decodedFrequenciesOrdered.count(it->first+2) ==  0 )
        {
	    // weighted average location
            realFrequenciesIndices.push_back( ratio * peak / totalEnergy );
	    
	    // round to the nearest integer location
            decodedFrequencies[(int) round(ratio * peak / totalEnergy)] = std::polar(sqrt(totalEnergy),0.0);

            // Begin Average Method
            decodedFrequencies[(int) round(ratio * peak / totalEnergy)] = 0;
	    
            for (int i = 0; i < nb; ++i)
            {
                decodedFrequencies[(int) round(ratio*peak/totalEnergy)] += 
		    (frontEnd->getTimeSignal().at(i)) * std::polar(1.0,-2*M_PI*(i)*(peak/totalEnergy)/config->getSignalLength());
            }
            decodedFrequencies[(int) round(ratio*peak/totalEnergy)] = 
		std::polar(sqrt(totalEnergy),std::arg(decodedFrequencies[(int) round(ratio*peak/totalEnergy)]));
            // End Average Method

            amplitude = 0;
            totalEnergy = 0;
            peak = 0;
        }
    }
}
