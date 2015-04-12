#include "ffast.h"
#include "customizedoutput.h"

#include <iomanip>
#include <iostream>


FFAST::FFAST(const Config* newConfig, Input* newInput, Output* newOutput): config(newConfig), input(newInput), output(newOutput)
{
    chrono    = new Chrono();
    frontEnd  = new FrontEnd(chrono, config, input);
    backEnd   = new BackEnd(chrono, config, frontEnd);
    output->setBackEnd(backEnd);
    iteration = 0;
}

FFAST::~FFAST()
{
    delete output;
    delete backEnd;
    delete frontEnd;
    delete input;
    delete config;
    delete chrono;
}

void FFAST::process()
{
    chrono->start("Global");

    iteration++;

    frontEnd->process();
    backEnd->process();

    chrono->stop("Global");
}

void FFAST::displayResults() const
{

    if ( config->needToCountSamples() )
    {
        std::cout << std::endl << "<===== NUMBER OF SAMPLES =====>" << std::endl
                  << frontEnd->getUsedSamplesNb() << " -> used samples" << std::endl;
	std::cout << std::setprecision(3) << frontEnd->getUsedSamplesNb()/((double) config->getSignalLengthOriginal())*100 << "\% samples touched" << std::endl;
    }
    
    std::cout << std::endl << "<===== AVERAGE TIME (in seconds) =====>" << std::endl;
    displayExecutionTime(chrono->average("FrontEnd"), chrono->average("BackEnd"), chrono->average("Global"));
    output->displayGlobalResults();
}

void FFAST::displayExecutionTime(double frontEndTime, double backEndTime, double globalTime) const
{
    double FFASTtime = frontEndTime+backEndTime;

    std::cout << std::setprecision(3) << std::scientific
              << frontEndTime << " -> execution of the front-end"  << std::endl
              << backEndTime << " -> execution of the back-end"  << std::endl
              << FFASTtime << " -> total execution of FFAST"  << std::endl;

    std::cout << globalTime << " -> global execution" << std::endl;

    std::cout.unsetf(std::ios_base::floatfield);
}

const std::vector<int> FFAST::getDelays() const
{
    return frontEnd->getDelays();
}
