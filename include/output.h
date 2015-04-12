#ifndef OUTPUT_H
#define OUTPUT_H

#include "frontend.h"
#include "backend.h"
#include "step.h"
#include "utils.h"

#include <fftw3.h>
#include <set>
#include <vector>
#include <map>
#include <unordered_map>

class Output: public Step
{
protected:
    BackEnd* backEnd;

public:
    Output(Chrono* newChrono, const Config* newConfig);
    virtual ~Output() {}
    virtual void process() {}
    virtual void displayGlobalResults() {}
    void setBackEnd(BackEnd* newBackend);
};

#endif // OUTPUT_H
