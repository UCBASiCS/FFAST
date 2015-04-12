#include "output.h"

#include <algorithm>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <stdlib.h>
#include <set>
#include <unordered_map>

Output::Output(Chrono* newChrono, const Config* newConfig): Step(newChrono, newConfig)
{
}

void Output::setBackEnd(BackEnd* newBackEnd)
{
    backEnd = newBackEnd;
}

