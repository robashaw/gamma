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

#ifndef OPTIMISER_HEADER
#define OPTIMISER_HEADER	


#include "molecule.hpp"
#include "fock.hpp"
#include "eigen_wrapper.hpp"
#include "ProgramController.hpp"
#include <vector>

double rfo(Vector& dx, Vector& grad, Matrix& hessian, double alpha, double stepsize); 
double trust_newton(Vector& dx, Vector& grad, Matrix& hessian, double stepsize); 
std::pair<double, double> rfo_prime(Vector& dy, Vector &grad, Matrix& hessian, double alpha); 

struct RHFCalculator { 
	Command& cmd; 
	SharedMolecule mol; 
	double old_e, delta_e, grad_norm, energy; 
	int iter; 
	
	RHFCalculator(Command& _cmd, SharedMolecule _m) : cmd(_cmd), mol(_m), old_e(0.0), iter(0) { }
	
	double operator()(const Vector& dx, Vector& grad, Matrix& hessian);
};

struct RHFOptimiser {

	Command& cmd; 
	SharedMolecule mol;
	RHFCalculator calc; 
	
	RHFOptimiser(Command& _cmd, SharedMolecule _m) : cmd(_cmd), mol(_m), calc(_cmd, _m) { }
	
	void optimise(); 
	void frequencies(Matrix& hessian); 
	void exponents(); 
	
}; 
#endif
