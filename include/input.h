#ifndef INPUT_H
#define INPUT_H

#include "step.h"
#include "utils.h"

#include <fftw3.h>
#include <set>
#include <unordered_map>
#include <vector>

class Input: public Step
{
public:
    Input(Chrono* newChrono, const Config* newConfig);
    virtual ~Input() {}
    virtual const std::unordered_map<int,ffast_complex>& getTimeSignal() const = 0;
    virtual void process() {}
    virtual void process(std::vector<int> delays) {}
    virtual ffast_complex getSignalAtIndex(int k) const = 0;
};

#endif // INPUT_H
