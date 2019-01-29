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

#ifndef SCFHEADERDEF
#define SCFHEADERDEF

#include "fock.hpp"
#include "eigen_wrapper.hpp"
#include "molecule.hpp"
#include "diis.hpp"


// Declare forward dependencies
class IntegralEngine;
class Command; 

// Begin class
class SCF
{
private:
	Command& cmd; 
	SharedMolecule molecule;
	Fock& focker;
	DIISEngine diis;
	double energy, last_energy, one_E, two_E, error, last_err;
public:
	// Constructor
	SCF(Command& c, SharedMolecule m, Fock& f);
	// Routines
	void calcE();
	double getEnergy() const { return energy; } 
	double calcE(const Matrix& hcore, const Matrix& dens, const Matrix& fock); 
	Vector calcErr(const Matrix& F, const Matrix& D, const Matrix& S, const Matrix& orthog);
	Vector calcErr();
	bool testConvergence(double val);
	void rhf(bool print = true);
	void uhf(bool print = true);
	void uhf_internal(bool print, UnrestrictedFock& ufocker);
};

#endif
