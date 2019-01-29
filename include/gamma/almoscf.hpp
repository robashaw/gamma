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

#ifndef ALMOSCFHEADERDEF
#define ALMOSCFHEADERDEF

#include "fock.hpp"
#include "eigen_wrapper.hpp"
#include "molecule.hpp"
#include "diis.hpp"
#include <vector>
#include "integrals.hpp"
#include "ProgramController.hpp"
#include <ctf.hpp>

/*! \class ALMOSCF
	\brief Container class for performing calculations using absolutely localised molecular orbitals.
	
	Given a molecule and a set of fragments, this allows for the calculation of SCF (+RPA) calculations
	on restricted or unrestricted systems.
 */
class ALMOSCF
{
private:
	SharedMolecule molecule; ///< Shared pointer to the Molecule, with information about fragments.
	Command& cmd; ///< Command with Options for this calculation.
	Fock& focker; ///< Fock object for the system-wide integrals, core Hamiltonian, and Fock matrix.
	std::vector<FockFragment> fragments; ///< Subsystem Fock objects for a restricted calculation.
	std::vector<UnrestrictedFockFragment> ufragments; ///< Subsystem Fock objects for an unrestricted calculation.
	std::vector<IntegralEngine> ints; ///< IntegralEngine for each fragment.
	DIISEngine diis; ///< Container for DIIS extrapolation of the ALMO SCF Fock matrices.
	
	double dimer_energy; ///< Total energy of the "dimer" (macrosystem) in Hartree.
	double e_frz; ///< Frozen orbital contribution to the energy. 
	double e_pol; ///< Polarization contribution to the energy after ALMO relaxation.
	double e_ct; ///< Charge transfer contribution to the energy.
	double e_int; ///< Intramolecular RPA energy
	double e_disp; ///< Dispersion and exchange-dispersion contribution to the energy.
	double e_pert_2; ///< 2nd-order perturbative correction to the ALMO SCF energy.
	double e_pert_4; ///< 4th-order perturbative correction to the ALMO SCF energy.
	double e_mon_rpa; ///< Cumulative RPA energy for each monomer.
	double delta_e; ///< Change in energy in the most recent SCF iteration.
	double delta_d; ///< Change in density matrix Frobenius norm in the most recent SCF iteration.
	
	std::vector<double> monomer_energies; ///< Ordered list of monomer total energies in Hartree.
	
	/// Basic information about the fragment subsystems, used in the local density-fitting procedure.
	std::vector<FragmentInfo> finfo, finfo_alpha, finfo_beta;
	 
	int nfrags; ///< The number of fragments in the system.
	int MAX; ///< The maximum number of DIIS error vectors to store.
	
	/// The total ALMO density matrix for the system.
	Matrix P, P_alpha, P_beta;
	/// The inverse overlap metric for the occupied subsystem.
	Matrix sigma, sigma_alpha, sigma_beta; 
public:
	/*! Creates an ALMOSCF container for the given Molecule, 
	 	using the options specified in the Command, and the given Fock object.
		@param c - the Command with the Options for this calculation.
		@param m - the SharedMolecule to perform the calculation on.
		@param f - the Fock object to be used.
	 */
	ALMOSCF(Command& c, SharedMolecule m, Fock& f);
	
	/*! Perturbatively correct the restricted ALMO SCF energy for charge transfer to 2nd or 4th order.
		@param order4 - performs 4th order perturbative correction if True, only 2nd order otherwise.
	 */
	void rperturb(bool order4 = false);
	
	/*! Perturbatively correct the restricted ALMO SCF energy for charge transfer to infinite order.
		Gives a directional, fragment pairwise decomposition of the charge transfer energy. 
	 */
	void rinf(); 
	
	/*! Perturbatively correct the unrestricted ALMO SCF energy for charge transfer to 2nd or 4th order.
		@param order4 - performs 4th order perturbative correction if True, only 2nd order otherwise.
	 */
	void uperturb(bool order4 = false);
	
	/*! Perturbatively correct the unrestricted ALMO SCF energy for charge transfer to infinite order.
		Gives a directional, fragment pairwise decomposition of the charge transfer energy. 
		UNDER CONSTRUCTION.
 	 */
	void uinf(); 
	
	/*! Performs the monomer SCF calculations.
		Makes the fragment Fock objects, calculates the monomer energies, and stores information about
		the fragments in the FragmentInfo objects.
		@param unrestricted - set to True if unrestricted calculations necessary; defaults to False.
	 */
	void setFragments(bool unrestricted = false);
	
	/// Returns the total energy of the system in Hartree.
	double getDimerEnergy() const { return dimer_energy; } 
	
	/// Returns a reference to the ordered list of fragment energies. 
	std::vector<double>& getMonomerEnergies() { return monomer_energies; }
	
	/// Preforms a restricted ALMO SCF calculation.
	void rscf();
	
	/// Performs an unrestricted ALMO SCF calculation.
	void uscf();
	
	/// Performs a single, restricted ALMO SCF iteration.
	void rcompute();
	
	/// Performs a single, unrestricted ALMO SCF iteration.
	void ucompute();
	
	/*! Computes the full density-fitted Fock matrix, correcting for the local exchange approximation.
	    For the restricted case only. It is usually a very small correction, and can be very expensive.
		@return The correction to the total energy from using the full Fock matrix, in Hartree.
	 */ 
	double r_energy_df(); 
	
	/*! Builds the unrestricted ALMO density matrix for either alpha or beta electrons.
	 	@param alpha - calculates for the alpha electrons if True, beta if False, and stores results accordingly.
		@return The change in Frobenius norm of the alpha/beta density matrix. 
	 */
	double makeDens(bool alpha);
	
	void pairwise_rpa(std::vector<std::pair<int, int>>& pairs, std::vector<std::pair<int, int>>& offsets, int L, int R, std::atomic<double>& result, bool withcore, CTF::World& dw); 
};

#endif
