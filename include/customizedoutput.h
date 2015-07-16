#ifndef CUSTOMIZEDOUTPUT_H
#define CUSTOMIZEDOUTPUT_H

#include "frontend.h"
#include "backend.h"
#include "step.h"
#include "utils.h"
#include "output.h"

#include <fftw3.h>
#include <set>
#include <vector>
#include <map>
#include <unordered_map>

class CustomizedOutput: public Output
{
public:
    CustomizedOutput(Chrono* newChrono, const Config* newConfig);
    ~CustomizedOutput() {}
    void process();
    void displayGlobalResults();
};

#endif // CUSTOMIZEDOUTPUT_H
