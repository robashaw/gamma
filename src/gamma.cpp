/*
*
*   PURPOSE: The main program loop for the molecular suite of ab initio quantum
*            chemistry routines.
* 
*   DATE              AUTHOR                CHANGES
*   ===========================================================================
*   23/09/15          Robert Shaw           Original code.
*
*/

#include <iostream>
#include <fstream>
#include <string>
#include "mpi.h"
#include "ioutil.hpp"
#include "ProgramController.hpp"

int main (int argc, char* argv[])
{
	MPI_Init(&argc, &argv);
	int np; 
	MPI_Comm_size(MPI_COMM_WORLD, &np);
	
	if (argc == 1) { 
		std::cerr << "You must provide an input file.\nUSAGE: gamma inputfile.mol" << std::endl;
	} else {
		// Get the input file name, and create the error and output file names
		std::string ifname = argv[1];
		std::string name = ifname;
		std::size_t pos = name.find('.');
		if (pos != std::string::npos) { 
			name.erase(pos, name.length());
		}    
		std::string efname = name + ".log";
		std::string ofname = name + ".out";

		// Open the file streams
		std::ifstream input(ifname);
		if (!input.is_open()) {
			std::cerr << "Could not open input file." << std::endl;
		} else {
			std::ofstream output, err; 
			bool ochk = open_outfile(ofname, output);
			bool echk = open_outfile(efname, err); 
			
			if (ochk && echk) {
      
				// Create the program controller
				std::shared_ptr<ProgramController> control = std::make_shared<ProgramController>(input, output, err, name);
				control->run(); 
		   
				// Close file streams
				output.close();
				err.close();
				
			} else {
				std::cerr << "Could not open output files." << std::endl; 
			}
			input.close();
		}
	}
	MPI_Finalize();
	return 0;
}

