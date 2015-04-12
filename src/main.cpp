#include "ffast.h"
#include "customizedinput.h"
#include "experimentinput.h"
#include "customizedoutput.h"
#include "experimentoutput.h"
#include "fftw.h"

#include <stdlib.h>
#include <time.h>
// some additional includes
#include <iomanip>
#include <iostream>
#include <complex.h>
#include <fftw3.h>

int main(int argc, char** argv)
{
    // Sets a random seed for each run of FFAST.
    srand48((long int) time (NULL));
    
    // Configurations are set through the command line arguments
    Config* configurations = new Config(argc, argv);
    // This is here measuring how long it takes to generate the input and output
    // and FFTW if it is asked
    Chrono* mChrono = new Chrono();

    // Input object
    Input* input;
    // Output object
    Output* output;
     
    // In the experiment mode, program uses randomly generated input and
    // compare the FFAST output to the input.
    // In the customized mode, program takes input from a file and outputs
    // the recovered signal to a file. 
    // Note that input and output to FFAST is made modular and for using FFAST
    // in other applications changing these objects would be enough.
    if ( configurations->isExperimentMode() )
    {
	input  = new ExperimentInput(mChrono, configurations);
	// output = new CustomizedOutput(mChrono, configurations);
	output = new ExperimentOutput(mChrono, configurations, dynamic_cast<ExperimentInput*>(input));
    }
    else
    {
	input  = new CustomizedInput(mChrono, configurations);    
	output = new CustomizedOutput(mChrono, configurations);
    }

    FFAST* ffast = new FFAST(configurations, input, output);
    
    if (!configurations->isHelpDisplayed())
    {
        int iterations = configurations->getIterations();

        configurations->display();

        for (int iteration=0; iteration<iterations; iteration++)
        {
	    input->process(ffast->getDelays());
            ffast->process();
	    output->process();
        }
        
        ffast->displayResults();
    }

    	
    if (configurations->needToCompareWithFFTW())
    {
	input->process();
	const std::unordered_map<int,ffast_complex> mymap = input->getTimeSignal();

	ffast_complex* inputSignal = (ffast_complex*) fftw_malloc(configurations->getSignalLength() * sizeof(ffast_complex));
	
	for ( auto it = mymap.begin(); it!= mymap.end(); ++it )
	{
	    inputSignal[it->first] = it->second;
	}
      
	FFTW* fftw = new FFTW(mChrono, configurations, inputSignal);
	fftw->process();
	
	std::cout << std::endl;
	std::cout << "<===== FFTW =====>" << std::endl;
	std::cout << std::setprecision(3) << std::scientific
                  << mChrono->average("FFTW") << " -> execution of FFTW"  << std::endl;

	delete fftw;
    }

    std::cout << std::endl;
	
    delete ffast;

    exit(EXIT_SUCCESS);
}

