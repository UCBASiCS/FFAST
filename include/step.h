#ifndef STEP_H
#define STEP_H

#include "chrono.h"
#include "config.h"

class Step
{
protected:
    Chrono* chrono;
    const Config* config;

public:
    Step(Chrono* newChrono, const Config* newConfig): chrono(newChrono), config(newConfig) {}
    virtual ~Step() {}
    virtual void process() {}
};

#endif // STEP_H
