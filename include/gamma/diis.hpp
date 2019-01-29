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

#ifndef DIISHEADERDEF
#define DIISHEADERDEF

#include <vector> 
#include "eigen_wrapper.hpp"

/*! \class DIISEngine
	\brief Performs generic DIIS extrapolations.

	DIISEngine stores error vectors from any DIIS extrapolation procedure
	(e.g. SCF or CC), and computes the extrapolation weights from them. 
	To do this it uses a damped solver that ensures the problem does not 
	become ill-conditioned; note that this means not all previous errors
	are necessarily used in each extrapolation. 
 */
class DIISEngine
{
private:
	int maxDiis; ///< Maximum number of error vectors to store.
	
	/*! True if DIIS extrapolation is active. Error vectors will be stored
	 	but weights not calculated if this is set to false. 
	 */
	bool useDiis; 
	
	std::vector<Vector> errs; ///< List of previous error vectors, in order.
	double damping_factor; ///< Factor by which eigenvalues of B matrix are damped.
	
	Matrix lastB; ///< Previous B matrix, for damping procedure. 
	
public:
	
	/// Creates an empty DIISEngine. 
	DIISEngine();
	
	/*! Initialises the DIISEngine.
	 	@param _maxDiis - the maximum number of error vectors to use in the DIIS extrapolation.
		@param _useDiis - whether the DIIS extrapolation is active straight away.
		@param _damp - the damping factor to be used in solving the matrix equations.
	 */
	void init(int _maxDiis, bool _useDiis, double _damp = 0.02);
	
	/*! Switches the DIIS extrapolation on or off.
		@param on - True will set DIIS to active, False to inactive.
	 */
	void use(bool on) { useDiis = on; }
	
	/*! Computes the DIIS weights, given a new set of errors. 
		@param errors - a collection of one or more error vectors, which will be concatenated and added to errs.
		@return A vector of the DIIS weights. 
	 */
	Vector compute(std::vector<Vector> &errors);
	
	/*! Solves the DIIS equations using a damping procedure. 
		@param B - the B matrix.
		@param start - the index of the eigenvalue to start with (eigenvalues in ascending order).
		@return A vector of the weights, with the last element the rank of the problem. 
	 */
	Vector solve(Matrix &B, int start);
};

#endif
