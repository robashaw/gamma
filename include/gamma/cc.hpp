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

#ifndef CCHEADERDEF
#define CCHEADERDEF

#include "tensor4.hpp"
#include <ctf.hpp>
#include "eigen_wrapper.hpp"
#include "mp2.hpp"
#include "diis.hpp"
#include <vector>
#include "ProgramController.hpp"

class IntegralEngine;

/*! \class CCSD
 	\brief Coupled-cluster singles and doubles container class.

	Calculates restricted CCSD and (T) energies for a Molecule using the 
	Cyclops Tensor Framework (and DIIS extrapolation - under construction). 
 */
class CCSD
{
private:
	int N; ///< The total number of orbitals
	int nocc; ///< The number of occupied orbitals.
	int iter; ///< Number of the current iteration.
	Command& cmd; ///< Reference to the Command with Options for calculation.
	MP2& mp2; ///< Reference to the MP2 object for transformed integrals.
	
	DIISEngine diis; ///< Container for DIIS extrapolation of the amplitudes.
	bool doDiis; ///< True if DIIS is to be used, False otherwise.
	int maxDiis; ///< Maximum number of error vectors to be used in DIIS extrapolation (default 4).
	std::vector<Matrix> singles_cache; ///< Cache of (up to maxDiis) previous singles amplitudes.
	std::vector<S4OddTensor4> doubles_cache; ///< Cache of (up to maxDiis) previous doubles amplitudes.

	double energy; ///< Current iteration CCSD energy in Hartree. 
	double delta_e; ///< The change in CCSD energy from the last iteration.
	double triples_energy; ///< The (T) energy in Hartree.
	double delta_singles; ///< Change in Frobenius norm of singles amplitudes from the last iteration.
	double delta_doubles; ///< Change in Frobenius norm of doubles amplitudes from the last iteration.
	bool withTriples; ///< True if (T) energy is to be calculated, false otherwise.
	
public:
	
	/*! Creates a CCSD object. 
	    @param c - the Command with the Options for the calculation.
	 	@param _mp2 - the MP2 object with the integrals and Fock matrix for the calculation.
	 */
	CCSD(Command& c, MP2& _mp2);
	
	/*! Calculates the (T) energy, given integrals and amplitudes.
		@param V - the container with the MO-basis integrals as tensors.
		@param T - the container with the singles and doubles CCSD amplitudes as tensors.
	 */
	void calculateTriples(Integrals& V, Amplitudes& T);
	
	/*! Calculates the DIIS error vector of the amplitudes. UNDER CONSTRUCTION.
	 	@param newSingles - the new singles amplitudes from this iteration.
		@param newDoubles - the new doubles amplitudes from this iteration.
	 */
	void calculateError(Matrix& newSingles, S4OddTensor4& newDoubles);
	
	/// Performs a single CCSD iteration, calculating the new amplitudes - SLOW VERSION.  
	void compute();
	
  	double getEnergy() const { return energy; } ///< @return The CCSD energy in Hartree. 
  	double getETriples() const { return triples_energy; } ///< @return The (T) energy in Hartree.
	
	/*! Performs a single CCSD iteration, updating the singles and doubles amplitudes.
		This is the fast version using the CTF, but it is very memory hungry. 
	
		@param V - the container with the transformed integrals.
		@param aT - the container with the amplitudes to be updated.
		@param sched_nparts - MPI variable, not currently used.
	 */
	void ccsd_iteration(Integrals& V, Amplitudes& aT, int sched_nparts = 0);
};

#endif 
