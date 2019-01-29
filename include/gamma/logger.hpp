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

#ifndef LOGGERHEADERDEF
#define LOGGERHEADERDEF

// Exceptional use of preprocessor here for fundamental constants
// purely because they aren't always defined by every compiler
#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

// Includes
#include <string>
#include <fstream>
#include <iostream>
#include "basis.hpp"
#include "atom.hpp"
#include "molecule.hpp"
#include "eigen_wrapper.hpp"
#include <libecpint/ecp.hpp>
#include <vector>
#include <chrono>
#include <libint2.hpp>

// Declare forward dependencies
class Error;
class Fock; 

// Begin class declaration
class Logger
{
private:
  ProgramController& control; 
  std::ofstream& outfile;
  std::ofstream intfile, optfile;
  std::ostream& errstream;
  Error* errs;
  int nerr;
  std::chrono::steady_clock::time_point last_time, first_time;
  
public:
  // Conversion factors
  static const double RTOCM;
  static const double RTOGHZ;
  static const double TOKCAL;
  static const double TOKJ;
  static const double TOEV;
  static const double TOBOHR;
  static const double TOANG;
  static const double TOWAVENUMBERS; 
  
  // Constructor/destructor
  Logger(ProgramController &control, std::ofstream& out, std::ostream& e);
  ~Logger(); // Delete the various arrays
  // Accessors
  Logger& operator=(const Logger& other);
  
  void init_intfile();
  void init_optfile(); 
  
  std::ofstream& getIntFile() { return intfile; }
  std::ofstream& getOptFile() { return optfile; }
  
  // Overloaded print functions
  void print(const std::string& msg) const; // Print a string message
  // Print out a vector with precision digits, either horizontally or vertically
  void print(const Vector& v, int digits = 6, bool vertical = false) const; 
  void print(const Matrix& m, int digits = 6) const; // Matrix with precision digits
  void print(Basis& b, SharedMolecule m, bool full = false) const; // Basis set - spec., no. of bfs, etc.
  void print(libecpint::ECPBasis& b, bool full = false) const; // ECP Basis set - spec., no. of bfs, etc.
  void print(libecpint::ECP& ecp, bool full = false) const; 
  void print(std::vector<libint2::Shell>& basis, std::vector<int>& shellAtomList, 
					SharedMolecule m, bool full = false) const; 
  void print(const Atom& a) const; // Atom - i.e coords, etc.
  // Print out the details of the molecule, including inertial data if wanted.
  void print(Molecule& mol, bool inertia = false) const; 

  // Print out an iteration
  void iteration(int iter, double energy, double delta, double dd);
  void iterationCC(int iter, double energy, double delta_e, double delta_s, double delta_d, double t_interm, double t_amps, double t_iter);
  void initIteration();
  void initIterationCC();
  void initALMOTable(); 
  void initCTTable(); 
  void ALMORow(int f1, int f2, double sep, double edisp, double edispexch, double eionic, double ebsse, double eintra); 
  void CTRow(int f1, int f2, double en); 
  void orbitals(const Vector& eps, int nel, bool one = false);
  void coefficient_matrix(const Vector& eps, int nel, const Matrix& coeffs, bool one = false); 
  void frequencies(const Vector& freqs, const Matrix& modes, bool printmodes); 
  void printDomains(Fock& f); 
  
  void initIterationOpt(); 
  void optg_dump(int iter, Vector& grad, Vector& step, SharedMolecule m, Matrix& hessian,
  				double trust, double delta_e, double grad_norm, double step_norm, 
				double energy, double expect);
  void optx_dump(int iter, Vector& grad, Vector& step, SharedMolecule m, Matrix& hessian,
			    double trust, double delta_e, double grad_norm, double step_norm, 
			  	double energy, double expect, std::vector<int>& activex);
				
  void mo_map(Vector& coeffs, SharedMolecule m, int fineness, std::string& filename);  
  
  // Specific logging formats
  void title(const std::string& msg) const;
  void result(const std::string& msg) const;
  void result(const std::string& name, const double value, const std::string& units) const;
  void error(Error& e);
  void localTime();
  void globalTime();
  void errTime();
  double getLocalTime();
  double getGlobalTime();
  void flush();
  void init();
  void finalise();
};
 
#endif

