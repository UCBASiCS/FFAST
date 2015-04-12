#ifndef BINPROCESSOR_H
#define BINPROCESSOR_H

#include "step.h"
#include "utils.h"

#include <vector>

class BinProcessor: public Step
{
private:
    ffast_complex** observationMatrix;
    int signalLength;
    int delaysNb;           // in draft paper : D=CN
    int chainsNb;           // in draft paper : C=noisy ? log(log(n)) : 1
    int delaysPerBunchNb;   // in draft paper : N=input[def:2]

    ffast_real minimumEnergy;

    int location;
    ffast_complex amplitude;
    ffast_complex* signalVector;

    int binAbsoluteIndex;
    int binRelativeIndex;
    int stage;
    int binSize;
    bool MLdetection;
    ffast_real twoPiBySignalLength;
    ffast_real signalLengthByTwoPi;
    ffast_real noise;
    ffast_real* thresholds;
    ffast_real* weights;
    ffast_real* angles;
    ffast_complex* directionVector;
    std::vector<int> delays;

public:
    BinProcessor(Chrono* newChrono, const Config* newConfig, ffast_complex** newObservationMatrix, const std::vector<int>& newDelays);
    ~BinProcessor();
    void process();
    void adjustTo(int newBinAbsoluteIndex, int newBinRelativeIndex, int newStage);
    bool isSingleton();
    int getLocation() const;
    ffast_complex getAmplitude() const;
    ffast_complex getSignal(int delayIndex) const;

private:
    bool isZeroTon() const;
    void MLprocess();
    void estimateBinSignal();
    void computeLocation();
    ffast_real getOmega(int i) const;
    void computeWeights();
    void computeThresholds();
};

#endif // BINPROCESSOR_H
