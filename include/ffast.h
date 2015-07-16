#ifndef FFAST_H
#define FFAST_H

#include "backend.h"
#include "chrono.h"
#include "config.h"
#include "fftw.h"
#include "frontend.h"
#include "input.h"
#include "output.h"

class FFAST
{
public:
    const Config* config;

private:
    Chrono* chrono;
    Input* input;
    FrontEnd* frontEnd;
    BackEnd* backEnd;
    Output* output;

    int iteration;

public:
    FFAST(const Config* newConfig, Input* newInput, Output* newOutput);
    ~FFAST();
    void process();
    void displayResults() const;
    const std::vector<int> getDelays() const;

private:
    void displayExecutionTime(double frontEndTime, double backEndTime, double globalTime) const;
};

#endif // FFAST_H
