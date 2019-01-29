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

#include "eigen_wrapper.hpp"
#include <vector> 

/*! Iteratively solves the CPHF equations via
	    U(n+1) = Q + V * U(n)
	until ||U(n+1) - U(n)|| is less than a tolerance or the maximum number
	of iterations has been reached.

	@param Q - the source term matrix
	@param V - the gradient term matrix
	@param tol - the desired convergence tolerance of the change in U.
	@param maxiter - the maximum number of iterations allowed to reach convergence.
	@return The solution vector, U
 */
Vector cphf_single_solver(Vector& Q, Matrix& V, double tol = 1e-6, int maxiter = 10); 

/*! Solves the CPHF equations for a group of systems, calling
 	cphf_single_solver for each system. 
	
	@param Q - a vector of source matrices, in same order as V.
	@param V - a vector of gradient matrices, in same order as Q.
	@param tol - the desired convergence tolerance for the change in the solution vector.
	@param maxiter - the maximum number of iterations allowed to reach convergence.
	@return A vector of solution vectors, U, for each system, in the same order as Q and V. 
 */
std::vector<Vector> cphf_group_solver(std::vector<Vector>& Q, Matrix& V, double tol = 1e-6, int maxiter = 10); 

