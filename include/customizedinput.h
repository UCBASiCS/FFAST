#ifndef CUSTOMIZEDINPUT_H
#define CUSTOMIZEDINPUT_H

#include "input.h"
#include "step.h"
#include "utils.h"
#include <set>
#include <unordered_map>

class CustomizedInput: public Input
{
private:
    std::unordered_map<int,ffast_complex> timeSignal;
    std::set<int> neededSamples;

public:
    CustomizedInput(Chrono* newChrono, const Config* newConfig);
    ~CustomizedInput();
    void process();
    void process(std::vector<int> delays);
    const std::unordered_map<int,ffast_complex>& getTimeSignal() const;
    ffast_complex getSignalAtIndex(int k) const;

private:
    void convertFileToTimeSignal();
    void findNeededSamples(std::vector<int> delays);
};

#endif // INPUT_H

