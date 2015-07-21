#include "config.h"
#include "lcm.h"
#include "split.h"

#include <algorithm>
#include <assert.h>
#include <fftw3.h>
#include <fstream>
#include <getopt.h>
#include <iostream>
#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <limits>



Config::Config(int newArgc, char** newArgv):argc(newArgc), argv(newArgv), helpDisplayed(false)
{
    setDefaultOptions();
    setOptionsFromCommandLine();
    if(bins.size() == 0)
    {
        proposedBins();
    }

    if (applyWindowVar)
    {
      signalSparsityPeeling = 3*signalSparsity;
    }

    maxSNRdB = (maxSNRdB < SNRdB) ? SNRdB : maxSNRdB;

    setBinOffsetsAndBinsSum();
    
    if(distribution.size() == 0)
    {
	preprocessDistribution((char*) "1");
    }

    // if the delays are not given as input by the user
    // and signal is noisy and/or we apply the window
    if (defaultDelays && (noisy || applyWindowVar))
    {
	double effectiveSNR = pow(10,getSNRdB()/10);

	if ( applyWindowVar )
	{
	    double mainLobePower = 2 * (  pow(sin(M_PI/2)/sin(M_PI/2/signalLength),2) 
		                        + pow(sin(3*M_PI/2)/sin(3*M_PI/2/signalLength),2) 
				        + pow(sin(5*M_PI/2)/sin(5*M_PI/2/signalLength),2)  );
	    double signalToNoiseContributionFromOffgrid = mainLobePower / ( pow(signalLength,2) - mainLobePower );
	    
	    offgridSNRdB = 10 * log10(signalToNoiseContributionFromOffgrid);

	    effectiveSNR = 1/( 1/effectiveSNR + 1/signalToNoiseContributionFromOffgrid );

	}
	
	
	// double delayScaling = 3 * pow(10,getSNRdB()/10)/(1 + 4*pow(10,getSNRdB()/20));
	
	double delayScaling = 3 * effectiveSNR / ( 1 + 4 * sqrt(effectiveSNR) );

	chainsNb = ceil( log(signalLength) / sqrt(delayScaling) );
	
	if (!maximumLikelihood) // Fast search
	{
	    delaysPerBunchNb = 2 * ceil( pow(log(signalLength),1.0/3.0) / sqrt(delayScaling) );
	}

	// the number of chains should be at least 1
	chainsNb = (chainsNb >= 1) ? chainsNb : 1;
	// the number of delays per chain should be at least equal to 2
	delaysPerBunchNb = (delaysPerBunchNb >= 2) ? delaysPerBunchNb : 2;

	if ( applyWindowVar )
	{

	    chainsNb *= (1 + log(signalSparsity)/log(10) );
	}
    }

    delaysNb = chainsNb * delaysPerBunchNb;
    
    if (!noisy && !applyWindowVar)
    {
	chainsNb = 1;
	delaysPerBunchNb = delaysNb;
    }

    assert( signalLength >= signalSparsity );
    assert( delaysNb >= 2 );
    assert( delaysNb <= signalLength/getBiggestBin() );

    // in fast search using Kay's method, we need at least
    // 2 delays within each bunch to calculate the locations
    if (!maximumLikelihood)
    {
	assert( delaysPerBunchNb >= 2 );
    }

}

void Config::setDefaultOptions()
{
    outputFile = (char*) "ffastOutput.txt";
    
    signalLength = 124950;
    signalLengthOriginal = signalLength;
    signalSparsityPeeling = 40;
    signalSparsity = 40;
    lengthFactor = 1; // n = LCM(Bins)*lengthfactor;
    maximumLikelihood = false;
    countSamples = true;
    delaysPerBunchNb = 2;
    chainsNb = 1;
    offgridSNRdB = 100;
    minFourierMagnitude = 1;

    /* for experiment mode */
    iterations = 1;
    experimentMode = true;
    // Phase = 0 implies phase of non-zero 
    // coefficients is uniformly random in [0,2*pi]. 
    phasesNb = 0; 
    FFTWstrategy = FFTW_ESTIMATE;
    compareWithFFTW = false;
    displayIterationTime = false;
    noisy = false;
    SNRdB = 50;
    maxSNRdB = -std::numeric_limits<float>::infinity();
    verbose = false;
    reconstructSignalInBackEnd = false;
    offGrid = false;
    defaultDelays = true;
}

bool Config::isOffTheGrid() const
{
    return offGrid;
}

bool Config::isExperimentMode() const
{
    return experimentMode;
}

float Config::getMinFourierMagnitude() const
{
    return minFourierMagnitude;
}

int Config::getIterations() const
{
    return iterations;
}

int Config::getSignalLength() const
{
    return signalLength;
}

bool Config::needToUseMaximumLikelihoodDetection() const
{
    return maximumLikelihood;
}

int Config::getSignalSparsity() const
{
    return signalSparsity;
}

int Config::getSignalSparsityPeeling() const
{
    return signalSparsityPeeling;
}

int Config::getBinSize(int stage) const
{
    return bins[stage];
}

int Config::getBinOffset(int stage) const
{
    return binOffsets[stage];
}

int Config::getBinsSum() const
{
    return binsSum;
}

int Config::getBiggestBin() const
{
    return bins[bins.size()-1];
}

int Config::getBinsNb() const
{
    return bins.size();
}

int Config::getDelaysNb() const
{
    return delaysNb;
}

int Config::getChainsNb() const
{
    return chainsNb;
}

int Config::getDelaysPerBunchNb() const
{
    return delaysPerBunchNb;
}

bool Config::isNoisy() const
{
    return noisy;
}

float Config::getSNRdB() const
{
    return SNRdB;
}

float Config::getOffgridSNRdB() const
{
    return offgridSNRdB;
}

float Config::getMaxSNRdB() const
{
    return maxSNRdB;
}

int Config::getPhasesNb() const
{
    return phasesNb;
}


bool Config::isVerbose() const
{
    return verbose;
}

bool Config::needToCompareWithFFTW() const
{
    return compareWithFFTW;
}

bool Config::needToDisplayIterationTime() const
{
    return displayIterationTime;
}

bool Config::needToReconstructSignalInBackEnd() const
{
    return reconstructSignalInBackEnd;
}

bool Config::needToCountSamples() const
{
    return countSamples;
}

int Config::getFFTWstrategy() const
{
    return FFTWstrategy;
}

void Config::display() const
{
    std::cout << std::endl << "<=======================>" << std::endl;
    
    if( isExperimentMode() )
    {
	
	std::cout << "Running experiment mode" << std::endl;
    }

    std::cout << std::endl << "<===== INPUT PROFILE =====>" << std::endl;

    std::cout << signalLength << " -> signal length" << std::endl;
    std::cout << signalSparsity << " -> signal sparsity" << std::endl;
    if (noisy)
    {
        std::cout << SNRdB << " dB" << " -> signal-to-noise ratio (SNR)" << std::endl;
    }
    else
    {
        std::cout << "Noiseless signal" << std::endl;
    }
    if (phasesNb == 0)
    {
        std::cout << "Random phase"  << std::endl;
    }
    else
    {
        std::cout << phasesNb << " -> possible phase";

        if (phasesNb > 1)
        {
            std::cout << "s";
        }

        std::cout << " for the non-zero frequencies in the input signal" << std::endl;
    }

    std::cout << std::endl << "<===== FFAST CONFIG =====>" << std::endl;
    
    std::cout << delaysNb << " -> delays" << std::endl;

    if( signalLengthOriginal != signalLength )
    {
	std::cout << "The signal length is not Chinese Remainder Theorem friendly." << std::endl; 
	std::cout << "The last " << (signalLengthOriginal - signalLength) << " samples will be deleted" << std::endl;
    }

    for (unsigned int stage=0; stage<bins.size(); stage++)
    {
        std::cout << bins[stage] << " ";
    }
    std::cout << "-> bins" << std::endl;

    if (lengthFactor > 1)
    {
        std::cout << lengthFactor << " -> length factor"  << std::endl;
    }

    if (iterations > 1)
    {
        std::cout << iterations << " -> iterations"  << std::endl;
    }

    

    

    if (maximumLikelihood)
    {
        std::cout << "Maximum likelihood detection enabled"  << std::endl;
    }

    if (reconstructSignalInBackEnd)
    {
        std::cout << "The back-end will reconstruct the full length frequency signal" << std::endl;
    }
    
    if (FFTWstrategy == FFTW_MEASURE)
    {
        std::cout << "Optimized FFTW plan generation enabled" << std::endl;
    }

    if (verbose)
    {
        std::cout << "Verbose mode enabled" << std::endl;

        if (displayIterationTime)
        {
            std::cout << "Prompting execution time for each iteration" << std::endl;
        }
    }

    if (compareWithFFTW)
    {
        std::cout << "Comparing with FFTW" << std::endl;
    }

    if(distribution.size() > 2) 
    {
        std::cout << "Non uniform distribution -> ";

        for (unsigned int i = 0; i < distribution.size(); ++i)
        {
            std::cout << " " << distribution[i];
        }
        
        std::cout << std::endl;
    }
}

void Config::help()
{
    helpDisplayed = true;

    std::cout   << std::endl << "<===== HELP =====>"
                << std::endl << std::endl
                << " [-a or --experiment]" << std::endl
                << "     Run experiment mode."  
                << std::endl << std::endl
                << " [-b BINS or --bins BINS]" << std::endl
                << "     Set the bins to use." << std::endl
                << "     Example: -b \"125 128 243\"."
                << std::endl << std::endl
                << " [-c or --samples]" << std::endl
                << "     Do not count the number of signal samples used in the front-end and display it in the results section."
                << std::endl << std::endl
                << " [-d NUM or --delays NUM]" << std::endl
                << "     Set the number of delays to use in the back-end." << std::endl
                << "     Note: the following assertions should be true (d >= 2) and (d < n/max(Bins))."
                << std::endl << std::endl
                << " [-e NUM or --chains NUM]" << std::endl
                << "     Set the number of chains to use in the back-end for fast search."
                << std::endl << std::endl
                << " [-f FNAME or --file FNAME]" << std::endl
                << "     Input file to use" << std::endl
                << "     Example: -f \"timeSignal.txt\"."
                << std::endl << std::endl
                << " [-g NUM or --minmagnitude NUM]" << std::endl
                << "     Minimum Fourier magnitude that will be recovered. Strongly advised to be entered in offgrid setting."
                << std::endl << std::endl
                << " [-h or --help]" << std::endl
                << "     Displays help."
                << std::endl << std::endl
                << " [-i NUM or --iterations NUM]" << std::endl
                << "     Set the number of iterations." << std::endl
                << "     When the algorithm is be executed several times, the results will be averaged."
                << std::endl << std::endl
                << " [-k NUM or --sparsity NUM]" << std::endl
                << "     Set the signal sparsity. The input signal will have k uniformly distributed non zero frequency." << std::endl
                << "     Note: the following assertions should be true (k <= n) and (k <= 2*min(Bins))."
                << std::endl << std::endl
                << " [-l or --ml]" << std::endl
                << "     Enable the maximum likelihood detection." << std::endl
                << "     Note: This may be required when the snr is really low but it will dramatically slow down the execution."
                << std::endl << std::endl
                << " [-m NUM or --factor NUM]" << std::endl
                << "     Set the length factor. It is useful for having small bins and large signal length." << std::endl
                << "     Note: The length of the signal is equal to lcm(Bins)*(length factor)."
                << std::endl << std::endl
                << " [-n NUM or --length NUM]" << std::endl
                << "     Set the signal length."
                << std::endl << std::endl
                << " [-o or --optimize]" << std::endl
                << "     Create optimized FFTW plans." << std::endl
                << "     Note: This option should be used only with small n (less than 5000) since creating plans for big signals is slow."
                << std::endl << std::endl
                << " [-r or --reconstruct]" << std::endl
                << "     Force the back-end to reconstruct the full n-length frequency signal." << std::endl
                << "     Note: this can increase the execution time of the FFAST algorithm."
                << std::endl << std::endl
                << " [-s NUM or --snr NUM]" << std::endl
                << "     Add noise to the input signal in frequency domain and specify the SNR in dB." << std::endl
                << "     Example: -s 15 will generate a noisy signal with a SNR equal to 15 dB."
                << std::endl << std::endl
                << " [-v or --verbose]" << std::endl
                << "     Enable the verbose mode: the details of each iteration will be displayed."
                << std::endl << std::endl
                << " [-u DIST or --distribution DIST]" << std::endl
                << "     Enable the non-uniform distribution mode: input the weights of your distribution." << std::endl
                << "     Example: -u \"3 2 1\" will divide your frequency domain in 3, the probability to be in the first third will be 1/2, the next 1/3 and the last 1/6"
                << std::endl << std::endl
                << " [-w or --fftw]" << std::endl
                << "     Compare the execution time of FFAST with FFTW." << std::endl
                << "     Note: FFTW will be exectuted at each iterarion and this will lead to slower execution." 
                << std::endl << std::endl
                << " [-x BINF or --bins_with_factor BINF]" << std::endl
                << "     Set the bins and the factor length at the same time." << std::endl
                << "     The syntax is the same than for setting bin but the last value is used for the factor." << std::endl
                << "     Example: -x \"49 50 51 10\" will set bins to \"49 50 51\" and length factor to 10."
                << std::endl << std::endl
                << " [-z FNAME or --write FNAME]" << std::endl
                << "     Write the output on the choosen file"
                << std::endl;
}

bool Config::isHelpDisplayed() const
{
    return helpDisplayed;
}

void Config::setOptionsFromCommandLine()
{
    const struct option longOptions[] =
    {	
	{"experiment",	 no_argument,	    NULL, 'a'},
	{"bins",	 required_argument, NULL, 'b'},
        {"samples",      no_argument,       NULL, 'c'},
        {"delays",       required_argument, NULL, 'd'},
        {"chains",       required_argument, NULL, 'e'},
	{"file",         required_argument, NULL, 'f'},
	{"minmagnitude", required_argument, NULL, 'g'},
	{"write",        required_argument, NULL, 'w'},
        {"help",         no_argument,       NULL, 'h'},
        {"iterations",   required_argument, NULL, 'i'},
        {"sparsity",     required_argument, NULL, 'k'},
        {"ml",           no_argument,       NULL, 'l'},
        {"factor",       required_argument, NULL, 'm'},
        {"length",       required_argument, NULL, 'n'},
        {"optimize",     no_argument,       NULL, 'o'},
        {"reconstruct",  no_argument,       NULL, 'r'},
        {"snr",          required_argument, NULL, 's'},
        {"distribution", required_argument, NULL, 'u'},
        {"verbose",      no_argument,       NULL, 'v'},
        {"fftw",         no_argument,       NULL, 'w'},
        {"maxsnr",       required_argument, NULL, 'x'},
        {0, 0, 0, 0}
    };

    int option = 0;

    while(option != -1)
    {
        option = getopt_long(argc, argv, "acf:d:e:hg:i:k:b:lm:n:op:rs:tu:vwx:", longOptions, NULL); //abfgjuxyz

        switch (option)
        {
	case 'a':
	    experimentMode = true;
	    break;
	    
	case 'b':
            setBins(optarg);
            break;

        case 'n':
            signalLength = atoi(optarg);
            signalLengthOriginal = atoi(optarg);
            break;

        case 'c':
            countSamples = false;
            break;

        case 'd':
            delaysPerBunchNb = atoi(optarg);
            defaultDelays = false;
            break;

        case 'e':
            chainsNb = atoi(optarg);
            defaultDelays = false;
            break;

	case 'f':
            inputFile = optarg;
            break;

        case 'h':
            help();
            break;

	case 'g':
	    minFourierMagnitude = strtof(optarg,0);
	    break;

        case 'i':
            iterations = atoi(optarg);
            break;

        case 'k':
            signalSparsity = atoi(optarg);
            signalSparsityPeeling = atoi(optarg);
            break;

        case 'l':
            maximumLikelihood = true;
            break;

        case 'm':
            lengthFactor = atoi(optarg);
            break;

        case 'o':
            FFTWstrategy = FFTW_MEASURE;
            break;

        case 'r':
            reconstructSignalInBackEnd = true;
            break;

        case 's':
            noisy = true;
            SNRdB = strtof(optarg, 0);
            break;

        case 'x':
            maxSNRdB = strtof(optarg, 0);
            break;

        case 'u':
            preprocessDistribution(optarg);
            break;

	case 'z':
            outputFile = optarg;
            break;

        case 'v':
            verbose = true;
            break;

        case 'w':
            compareWithFFTW = true;
            break;
        }
    }
}

void Config::setBins(const char *newBins)
{
    // split: divides a string into multiple strings that are separated by space ' '.
    std::vector<std::string> tempBins = split(newBins, ' ');
    bins.clear();
    for(unsigned int stage=0; stage<tempBins.size(); stage++)
    {
        int value = atoi(tempBins[stage].c_str());
        bins.push_back(value);
    }
}

void Config::setBinOffsetsAndBinsSum()
{
    std::sort(bins.begin(), bins.end());

    binsSum = 0;

    for (unsigned int stage=0; stage<bins.size(); stage++)
    {
        binOffsets.push_back(binsSum);
        binsSum += bins[stage];
    }
}

std::vector<double> Config::getDistribution() const
{
  return distribution;
}

void Config::preprocessDistribution(char* newDistribution)
{
    std::vector<std::string> tempDistribution = split(newDistribution, ' ');

    distribution.clear();
    distribution.push_back(0);

    int l = tempDistribution.size()+1;

    for (int i = 1; i < l; ++i)
    {
	distribution.push_back(distribution[i-1]+atof(tempDistribution[i-1].c_str()));
    }

    for (int i = 1; i < l; ++i) 
    {
	distribution[i] = distribution[i]/distribution[l-1];
    }
}

void Config::proposedBins(int newRangeSemiLength) {
    int F = pow(signalLength,1.0/3.0);
    int rangeSemiLength = F > newRangeSemiLength ? newRangeSemiLength : F-2;
    int bestLength = pow(F-rangeSemiLength,3);
    int p,r,f3;
    std::vector<int> bestBins,testedBins;

    for (int f1 = F-rangeSemiLength; f1 <= F+rangeSemiLength; ++f1) {
        testedBins.push_back(f1);
        p = f1*(F-rangeSemiLength); // p = f1.f2

        for (int f2 = F-rangeSemiLength; f2 <= F+rangeSemiLength; ++f2) {
            testedBins.push_back(f2);
            r = floor(((double) signalLength)/((double) p));

            for(int n = p*r; n > bestLength; n-=p) {
                f3 = n/p;
                if(f3 <= F+rangeSemiLength && f3 >= F-rangeSemiLength) {
                    testedBins.push_back(f3);

                    if(areCoprime(testedBins)) {
                        bestLength = n;
                        bestBins = testedBins;
                    }

                    testedBins.pop_back();
                }
            }
            testedBins.pop_back();
            p += f1;
        }
        testedBins.pop_back();
    }

    bins = bestBins;
    applyWindowVar = (signalLength != bestLength) || applyWindowVar;

    signalLength = bestLength;
}

bool Config::applyWindow() const
{
    return applyWindowVar;
}

int Config::getSignalLengthOriginal() const
{
    return signalLengthOriginal;
}

char* Config::getInputFile() const
{
    return inputFile;
}

char* Config::getOutputFile() const
{
    return outputFile;
}
