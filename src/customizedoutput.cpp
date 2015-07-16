#include "customizedoutput.h"

#include <algorithm>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <stdlib.h>
#include <set>
#include <unordered_map>

CustomizedOutput::CustomizedOutput(Chrono* newChrono, const Config* newConfig):
    Output(newChrono, newConfig)
{
}

void CustomizedOutput::process()
{
    chrono->start("Output");

    std::unordered_map<int,ffast_complex> decodedFrequencies;

    decodedFrequencies = backEnd->getDecodedFrequencies();

    std::ofstream myfile (config->getOutputFile());
    if (myfile.is_open())
    {
        myfile.precision(15);
        for (auto i = decodedFrequencies.cbegin(); i != decodedFrequencies.cend(); ++i)
            myfile << i->first << ":" << i->second << std::endl;
        myfile.close();
    }

    std::cout << "CUSTOMIZED OUTPUT LINE 33" << std::endl;
    
    chrono->stop("Output");
}

void CustomizedOutput::displayGlobalResults()
{
    std::cout << "written file" << std::endl;
}
