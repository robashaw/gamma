/* 
 * 	Copyright (c) 2017 Robert Shaw
 *
 * 	Permission is hereby granted, free of charge, to any person obtaining
 *	a copy of this software and associated documentation files (the
 * 	"Software"), to deal in the Software without restriction, including
 * 	without limitation the rights to use, copy, modify, merge, publish,
 * 	distribute, sublicense, and/or sell copies of the Software, and to
 * 	permit persons to whom the Software is furnished to do so, subject to
 *	the following conditions:
 *
 *	The above copyright notice and this permission notice shall be
 * 	included in all copies or substantial portions of the Software.
 *
 *	THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
 *	EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
 *	MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
 *	NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE
 *	LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION
 *	OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION
 *	WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
 */

#ifndef BASISREADERHEADERDEF
#define BASISREADERHEADERDEF

// Includes
#include <string>
#include <fstream>
#include <map>
#include "eigen_wrapper.hpp"
#include <libecpint/gshell.hpp>
#include <libecpint/ecp.hpp>
#include "pugixml.hpp"

// Declare forward dependencies
class Basis;

/*! \class BasisReader
 	\brief Reads in basis sets from the basis set library. 
	
	Given a list of atom types and their basis names, BasisReader reads the basis sets in 
	from the basis set library and adds these shells to a Basis object. This covers
	both orbital and auxiliary sets, and ECPs. 
 */
class BasisReader
{
private:
  /*! A map of atom types (by atomic number) and their basis names. The zero element
   	  is taken to be the default basis, which is STO-3G if not otherwise specified.
   */
  std::map<int, std::string> names; 
  
  std::string libpath; //< The path to the basis library 
  
  /*! Maps of all possible basis set names to the files in the basis set library.
   	  These lists are read in from list files in the basis set library when the 
      BasisReader is created, and include synonyms. 
   */
  std::map<std::string, int> obs_list, jk_list, mp2_list, ecp_list; 
  
  /*! Reads in the basis set list files of type defined by name. These are the .list
      files in the basis set library share. It converts these into a map from name
  	  to file index.
  	
	  @param name - the name of the basis type; currently "basis", "jkfit", "mp2fit", or "ecp".
  	  @return A map of names to file indices, to be stored and read by read_basis. 
   */
  std::map<std::string, int> read_basis_list(std::string name); 
  
  /*! Reads the basis set for a given atom name from the corresponding file given a map
  	  of basis names to file indices. 
  	  @param atom - the atomic symbol (lower case).
  	  @param pos - the xyz coordinates of the atom in question.
  	  @param basis - the name of the basis set.
  	  @param name - the name of the basis type ("basis", "jkfit", "mp2fit", or "ecp")
  	  @param basis_list - reference to the map of basis names to file indices.
      @return A vector of GaussianShells corresponding to the basis set for the given atom.
   */
  std::vector<libecpint::GaussianShell> read_basis(std::string atom, double* pos,
  	 									std::string basis, std::string name, 
  										std::map<std::string, int>& basis_list); 

public: 
	
  /*! Creates a BasisReader object. 
	  @param ns - a map of atom types (as atomic numbers) to basis set names.
   */ 
  BasisReader(std::map<int, std::string> ns); // Constructor

  /*! Reads an ECP basis set from the basis set library. 
      @param q - the atomic number of the atom with the ECP. 
      @param ecpset - the ECPBasis to which the shells should be added when read in.
      @param pos - the xyz coordinates of the atom.
      @return An ECP object containing the ECP shells for the specified atom type and position.
   */
  libecpint::ECP readECP(int q, libecpint::ECPBasis& ecpset, double *pos); 
  
  /*! Reads an orbital or auxiliary basis for a given atom type and adds the shells to a Basis object.
      @param b - the Basis to which the shells should be added when read in. 
      @param q - the atomic number of the atom.
      @param pos - the xyz coordinates of the atom.
      @param atom - the index of the atom.
      @param type - the basis set type (0 = orbital, 1 = jkfit, 2 = mp2fit)
   */ 
  void readShellBasis(Basis& b, int q, double *pos, int atom, int type = 0);
};

#endif
