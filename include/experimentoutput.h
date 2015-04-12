#ifndef EXPERIMENTOUTPUT_H
#define EXPERIMENTOUTPUT_H

#include "frontend.h"
#include "backend.h"
#include "step.h"
#include "utils.h"
#include "output.h"
#include "experimentinput.h"

#include <fftw3.h>
#include <set>
#include <vector>
#include <map>
#include <unordered_map>

class ExperimentOutput: public Output
{
private:
    const ExperimentInput* input;

    ffast_complex* inputFrequencySignal;
    ffast_complex* outputFrequencySignal;

    int falseDetectionsIterationNb;
    int missedLocationsIterationNb;
    bool binningFailure;
    ffast_real MSEerrorAmplitude;
    bool fullRecovery;

    std::vector<unsigned int> falseDetectionsNbs;
    std::vector<unsigned int> missedLocationsNbs;
    std::vector<double> MSEerrors;
    int binningFailuresNb;
    std::vector<ffast_real> MSEerrorsAmplitude;
    int fullRecoveriesNb;

    ffast_real averageFalseDetectionsNb;
    ffast_real maxFalseDetectionsNb;
    ffast_real averageMissedLocationsNb;
    ffast_real maxMissedLocationsNb;
    ffast_real averageMSEerrorAmplitude;
    ffast_real maxMSEerrorAmplitude;
    ffast_real averageMSEerror;
    ffast_real maxMSEerror;

    ffast_real probabilityOfFalseDetection;
    ffast_real probabilityOfMissedLocation;
    ffast_real probabilityOfBinningFailure;
    ffast_real probabilityOfFullRecovery;

public:
    ExperimentOutput(Chrono* newChrono, const Config* newConfig, const ExperimentInput* newInput);
    ~ExperimentOutput() {}
    void process();
    void displayGlobalResults();
    void displayIterationResults() const;
    int getFullRecoveriesNb() const;
    void writeOutputToFile();

private:
    void checkFalseDetections();
    void checkMissedLocationsOrBinningFailure();
    void checkFullRecovery();
    void computeMSEErrorAmplitude();
    void computeMSEerror();
    bool isBinningFailure(std::set<int> missedLocations) const;
    void setGlobalResultsFigures();
    template <class tempVector>
    void setAverageAndMax(tempVector& array, ffast_real& average, ffast_real& max, int* counter = NULL);
    ffast_complex* getClusteredFrequencies();
};

#endif // EXPERIMENTOUTPUT_H
