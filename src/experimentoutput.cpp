#include "experimentoutput.h"

#include <algorithm>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <stdlib.h>
#include <set>
#include <unordered_map>

ExperimentOutput::ExperimentOutput(Chrono* newChrono, const Config* newConfig, const ExperimentInput* newInput):
Output(newChrono, newConfig), input(newInput)
{
    binningFailuresNb = 0;
    fullRecoveriesNb = 0;
}

void ExperimentOutput::process()
{
    chrono->start("Output");

    std::unordered_map<int,ffast_complex> decodedFrequencies;

    // If the signal was off-grid adjust the output to the original signal
    if ( config->getSignalLengthOriginal() != config->getSignalLength() )
    {
        for (auto it = backEnd->getDecodedFrequencies().cbegin(); it != backEnd->getDecodedFrequencies().cend(); ++it)
        {
            if (input->getNonZeroFrequencies().count(it->first+1) > 0)
            {
                decodedFrequencies[it->first+1] = it->second;
            }
            else if (input->getNonZeroFrequencies().count(it->first-1) > 0)
            {
                decodedFrequencies[it->first-1] = it->second;
        }
        else if (input->getNonZeroFrequencies().count(it->first+2) > 0)
        {
            decodedFrequencies[it->first+2] = it->second;
        }
        else if (input->getNonZeroFrequencies().count(it->first-2) > 0)
        {
            decodedFrequencies[it->first-2] = it->second;
        }
        else
        {
            decodedFrequencies[it->first] = it->second;
        }
    }
    backEnd->swapDecodedFrequencies(decodedFrequencies);
    decodedFrequencies.clear();
}

    /*
        we first check if there are false detections. If there are 
        not any false detections we check if the peeling got stuck 
        due to the graph. If so we declare a binning failure and 
        do not count missed locations. If it was not a binning failure,
        then we count the missed locations.
    */
        checkFalseDetections();
        checkMissedLocationsOrBinningFailure();
        checkFullRecovery();
        computeMSEErrorAmplitude();
        computeMSEerror();

        chrono->stop("Output");

        if ( config->isVerbose() && config->getIterations() >= 2 )
        {
            displayIterationResults();
        }

    }

    void ExperimentOutput::displayIterationResults() const
    {
        std::cout << std::endl << "### Iteration results ###" << std::endl << MSEerrorAmplitude << " -> MSE error Amplitude" << std::endl;

        if (falseDetectionsIterationNb > 0)
        {
            std::cout << falseDetectionsIterationNb << " -> false detections" << std::endl;
        }
        else if (missedLocationsIterationNb > 0)
        {
            std::cout << missedLocationsIterationNb << " -> missed locations" << std::endl;
        }
        else if (binningFailure)
        {
         std::cout << "Binning failure" << std::endl;
     }
     else
     {
        std::cout << "Full recovery" << std::endl;
    }
}

void ExperimentOutput::displayGlobalResults()
{
    setGlobalResultsFigures();

    std::cout 	<< std::endl << "<===== FFAST RESULTS =====>" << std::endl;
    std::cout 	<< std::setprecision(4)  << averageMSEerrorAmplitude << " -> MSE error amplitude" << std::endl;
    // std::cout   << maxMSEerrorAmplitude << " -> maximum MSE amplitude"   << std::endl;
    std::cout 	<< "average freq error is " << averageMSEerror*100 << " percent of 2PI/N = " << 2*M_PI/config->getSignalLengthOriginal() << std::endl;
    // std::cout   << maxMSEerror          << " -> max MSE error frequencies"    << std::endl;
    std::cout   << averageFalseDetectionsNb << " -> average false detections" << std::endl;
    std::cout   << maxFalseDetectionsNb     << " -> maximum false detections" << std::endl;
    std::cout   << averageMissedLocationsNb << " -> average missed locations" << std::endl;
    std::cout   << maxMissedLocationsNb     << " -> maximum missed locations" << std::endl;
              // << probabilityOfFalseDetection << " -> probability of one or more false detection" << std::endl
              // << probabilityOfMissedLocation << " -> probability of one or more missed location" << std::endl
    std::cout << probabilityOfBinningFailure << " -> probability of binning failure" << std::endl;
    // if there were binning failures, list the sizes and counts of trapping sets
    if ( probabilityOfBinningFailure != 0 )
    {
        std::cout << "trapping sets (size, count): ";
        for ( auto it = trappingSetSizeCounter.cbegin(); it != trappingSetSizeCounter.cend(); ++it )
        {
            std::cout << "( " << it->first << " , " << it->second << " ), ";
        }
        std::cout << std::endl;
    }
    std::cout  << probabilityOfFullRecovery   << " -> probability of full recovery"   << std::endl;
}

int ExperimentOutput::getFullRecoveriesNb() const
{
    return fullRecoveriesNb;
}

void ExperimentOutput::checkFalseDetections()
{
    // number of false detections in the current iteration of the algorithm
    falseDetectionsIterationNb = 0;

    // go over the decoded frequencies by the backEnd
    for (auto it = backEnd->getDecodedFrequencies().cbegin(); it != backEnd->getDecodedFrequencies().cend(); ++it)
    {
        // if decoded frequency is a false detection
        if (input->getNonZeroFrequencies().count(it->first) == 0)
        {
            falseDetectionsIterationNb++;
        }
    }

    // addd the results of current iteration to the global couter
    falseDetectionsNbs.push_back(falseDetectionsIterationNb);
}

void ExperimentOutput::checkMissedLocationsOrBinningFailure()
{
    // indices of missed locations
    std::set<int> missedLocations;
    binningFailure = false;


    // if there are no false detections in the current iteration
    if (falseDetectionsIterationNb == 0)
    {
        // go over the true frequencies
        for (auto it = input->getNonZeroFrequencies().cbegin(); it != input->getNonZeroFrequencies().cend(); ++it)
        {
            // if the frequency is missing in the decoded frequencies
            if (backEnd->getDecodedFrequencies().count(it->first) == 0)
            {
                // add the missing location to the missed locations
                missedLocations.insert(it->first);
            }
        }

        // if there is a missing location check if it is due to binning failure
        if (!missedLocations.empty())
        {
            binningFailure = isBinningFailure(missedLocations);
        }
    }

    if (binningFailure)
    {
        binningFailuresNb++;
        missedLocationsIterationNb = 0;
        
        // if it is a binning failure, increment the counter for trapping set size
        // that is held with the unordered_map trappingSetSizeCounter
        auto it = trappingSetSizeCounter.find( missedLocations.size() );
        if( it != trappingSetSizeCounter.end() )  
        {
            it->second =  it->second + 1;
        }
        else
        {
            trappingSetSizeCounter[ missedLocations.size() ] = 1;
        }
    }
    else
    {
        missedLocationsIterationNb = missedLocations.size();
    }

    missedLocationsNbs.push_back(missedLocationsIterationNb);
}

void ExperimentOutput::checkFullRecovery()
{
    if (falseDetectionsIterationNb == 0 && missedLocationsIterationNb == 0 && binningFailure == false)
    {
       fullRecovery = true;
       fullRecoveriesNb++;
   }
   else
   {
    fullRecovery = false;
}
}

// MSE is calculated if we are successful
void ExperimentOutput::computeMSEErrorAmplitude()
{
    if (falseDetectionsIterationNb == 0 && missedLocationsIterationNb == 0 && binningFailure == false)
    {
        ffast_real denominator = 0, numerator = 0;
        for (auto i = input->getNonZeroFrequencies().cbegin(); i != input->getNonZeroFrequencies().cend(); ++i)
        {
            numerator+=std::norm(i->second-backEnd->getDecodedFrequencies().at(i->first));
            denominator+=std::norm(i->second);
        }
        MSEerrorsAmplitude.push_back(numerator/denominator);
    }
}

void ExperimentOutput::computeMSEerror()
{
    if (falseDetectionsIterationNb == 0 && missedLocationsIterationNb == 0 && binningFailure == false)
    {
        ffast_real MSEerror = 0;
        for (auto i = backEnd->getRealFrequenciesIndices().cbegin(); i != backEnd->getRealFrequenciesIndices().cend(); ++i)
        {
            bool found = false;
            auto it = input->getNonZeroFrequencies().cbegin();
            while (!found && it != input->getNonZeroFrequencies().cend())
            {
                if(std::abs((double)it->first-(*i))<=3)
                {
                    MSEerror += std::abs( it->first-(*i) );
                    found = true;
                }
                ++it;
            }
        }
        MSEerrors.push_back(MSEerror/config->getSignalSparsity());
    }
}

bool ExperimentOutput::isBinningFailure(std::set<int> missedLocations) const
{
    /*
        if at least one of the missed locations appear in a singleton bin, that is,
        it could have been peeled, declare NOT a bining failure. If all of the missed
        locations appear in multiton bins, than it is a binning failure.
    */
    for (int stage=0; stage<config->getBinsNb(); stage++)
    {
        std::vector<int> binStatus(config->getBinSize(stage), 0);

        for (auto it = missedLocations.cbegin(); it != missedLocations.cend(); ++it)
        {
            binStatus[*it % config->getBinSize(stage)]++;
        }

        for (unsigned int index=0; index<binStatus.size(); index++)
        {
            if (binStatus[index] == 1)
            {
                return false;
            }
        }
    }
    return true;
}

void ExperimentOutput::setGlobalResultsFigures()
{
    averageMSEerrorAmplitude    = 0;
    maxMSEerrorAmplitude        = 0;
    averageMissedLocationsNb    = 0;
    maxMissedLocationsNb        = 0;
    averageFalseDetectionsNb    = 0;
    maxFalseDetectionsNb        = 0;
    averageMSEerror             = 0;
    maxMSEerror                 = 0;
    int falseDetectionsNb       = 0;
    int missedLocationsNb       = 0;

    setAverageAndMax(MSEerrorsAmplitude, averageMSEerrorAmplitude, maxMSEerrorAmplitude);
    setAverageAndMax(missedLocationsNbs, averageMissedLocationsNb, maxMissedLocationsNb, &missedLocationsNb);
    setAverageAndMax(falseDetectionsNbs, averageFalseDetectionsNb, maxFalseDetectionsNb, &falseDetectionsNb);
    setAverageAndMax(MSEerrors, averageMSEerror, maxMSEerror);

    probabilityOfFalseDetection = ((ffast_real) falseDetectionsNb)/ config->getIterations();
    probabilityOfMissedLocation = ((ffast_real) missedLocationsNb)/ config->getIterations();
    probabilityOfBinningFailure = ((ffast_real) binningFailuresNb)/ config->getIterations();
    probabilityOfFullRecovery = ((ffast_real) fullRecoveriesNb)/ config->getIterations();
}

template <class tempVector>
    void ExperimentOutput::setAverageAndMax(tempVector& array, ffast_real& average, ffast_real& max, int* counter)
    {
        for (unsigned int i=0; i<array.size(); i++)
        {
            average += array[i];

            if (max < array[i])
            {
                max = array[i];
            }

            if (counter != NULL && array[i] > 0)
            {
                (*counter)++;
            }
        }
        if (array.size() > 0)
        {
            average /= array.size();
        }
    }
