/*
*
*   PURPOSE: To implement class IntegralEngine, which calculates, stores,
*            and processes the molecular integrals needed in ab initio
*            quantum chemistry calculations.
*
*   DATE         AUTHOR           CHANGES
*   ======================================================================
*   02/09/15     Robert Shaw      Original code.
*   03/09/15     Robert Shaw      Merged overlap and kinetic integrals.
*   04/09/15     Robert Shaw      Now supports general contracted.
*   06/09/15     Robert Shaw      Nuclear attraction ints for prims.
*   07/09/15     Robert Shaw      formNucAttract() now works, as does
*                                 makeSpherical(ints, lnums)
*   08/09/15     Robert Shaw      Auxiliary 2e- integrals, twoe
*   09/09/15     Robert Shaw      Shell 2e- integrals, up to (m0|pq)
*                                 - need to sphericalise, increment, 
*                                 then sphericalise again.
*   10/09/15     Robert Shaw      Finished shell twoe. Got rid of 
*   							   makeContracted and makeSpherical
*                                 for 2e ints, as absorbed into twoe.
*/

#include "error.hpp"
#include "integrals.hpp"
#include "mathutil.hpp"
#include "basis.hpp"
#include "ioutil.hpp"
#include "logger.hpp"
#include "tensor4.hpp"
#include <libecpint/ecp.hpp>
#include <libecpint/ecpint.hpp>
#include <libecpint/multiarr.hpp>
#include <libecpint/gshell.hpp>
#include "ProgramController.hpp"

#include <cmath>
#include <iomanip>
#include <string>
#include <thread>   

// Constructor
IntegralEngine::IntegralEngine(SharedMolecule m, bool print) : molecule(m)
{
	// Calculate sizes
	int natoms = molecule->getNAtoms();
	int M = nbasis(molecule->getBasis().getIntShells()); 
	int N = ncart(molecule->getBasis().getIntShells());
	// Cartesian is easy - there are (N^2+N)/2
	// unique 1e integrals and ([(N^2+N)/2]^2 + (N^2+N)/2)/2
	// unique 2e integrals
	int ones = (N*(N+1));
	sizes.resize(4);
	sizes[0] = ones;
	sizes[1] = (ones*(ones+1))/4;
  
	ones = (M*(M+1));
	sizes[2] = ones;
	sizes[3] = (ones*(ones+1))/4;
	

	if (print) {
		molecule->control->log.title("INTEGRAL GENERATION");
		molecule->control->log.print("Forming the one electron integrals\n");
	}
	
	auto &shells = molecule->getBasis().getIntShells();
	std::vector<Atom> atoms;
	for (int i = 0; i < molecule->getNAtoms(); i++) atoms.push_back(molecule->getAtom(i));
	
	sints = compute_1body_ints(shells, libint2::Operator::overlap);
	tints = compute_1body_ints(shells, libint2::Operator::kinetic);
	naints = compute_1body_ints(shells, libint2::Operator::nuclear, atoms);

	//std::cout << naints + tints << std::endl << std::endl;
	if(molecule->getBasis().hasECPS()) {
		buildTransMat();
		naints = naints + compute_ecp_ints(shells);
	}
	//std::cout << std::endl << naints + tints << std::endl;
	
	if (print) {
		molecule->control->log.print("One electron integrals complete\n");
		molecule->control->log.localTime();
	}
	
	Vector ests = getEstimates();
	if ( molecule->control->get_option<double>("memory") < ests(3) && !molecule->control->get_option<bool>("direct")) {
		Error e("MEMERR", "Not enough memory for ERIs.");
		if (print) {
			molecule->control->log.error(e);
			molecule->control->set_option<bool>("direct", true);
		}
	}
	
	prescreen = compute_schwarz_ints<>(shells);
	if (prescreen.rows() < 10 && print) { 
		molecule->control->log.print("Forming the two electron repulsion integrals.\n");
		molecule->control->log.print("PRESCREENING MATRIX:\n");
		molecule->control->log.print(prescreen);
		molecule->control->log.print("\n\n");
	}
		
	if ( molecule->control->get_option<bool>("direct") ){
		if(print) molecule->control->log.print("Two electron integrals to be calculated on the fly.\n");
	} else { // Check memory requirements
		twoints = compute_eris(shells);
		
		if(print) molecule->control->log.print("Two electron integrals completed.\n");
		std::string mem = "Approximate memory usage = ";
		mem += std::to_string(ests(3));
		mem += " MB\n";
		if(print) {
			molecule->control->log.print(mem);
			molecule->control->log.localTime();
		}
		
		if (molecule->control->get_option<bool>("printeris")) {
			molecule->control->log.print("Writing ERIs to file.\n");
			printERI(molecule->control->log.getIntFile(), M);
		}	
	} 
	
	molecule->control->log.flush();
}

IntegralEngine::IntegralEngine(SharedMolecule m, const IntegralEngine& ints, int start, int finish) : molecule(m)
{
	int nbfs = finish - start; 
	
	int ones = (nbfs*(nbfs+1));
	sizes.resize(4);
	sizes[0] = ones;
	sizes[1] = (ones*(ones+1))/4;
	sizes[2] = ones;
	sizes[3] = (ones*(ones+1))/4;

	Matrix S = ints.getOverlap();
	sints = S.block(start, start, nbfs, nbfs); 
	Matrix T = ints.getKinetic();
	tints = T.block(start, start, nbfs, nbfs);
		
	auto &shells = molecule->getBasis().getIntShells();
	std::vector<Atom> atoms;
	for (int i = 0; i < molecule->getNAtoms(); i++) atoms.push_back(molecule->getAtom(i));
	
	naints = compute_1body_ints(shells, libint2::Operator::nuclear, atoms);
	
	if(molecule->getBasis().hasECPS()) {
		naints = naints + compute_ecp_ints(shells);
	}
	
	prescreen = compute_schwarz_ints<>(shells);
	
	if ( !molecule->control->get_option<bool>("direct") ) {
		twoints.assign(nbfs, 0.0);
		for (int i = 0; i < nbfs; i++)
			for (int j = 0; j <= i; j++)
				for (int k = 0; k <= i; k++)
					for (int l = 0; l <= k; l++)
						twoints(i, j, k, l) = ints.getERI(i+start, j+start, k+start, l+start);
	}
}

IntegralEngine::IntegralEngine(const IntegralEngine& other) : molecule(other.molecule) {
	sints = other.sints;
	tints = other.tints;
	transmat = other.transmat; 
	naints = other.naints;
	prescreen = other.prescreen;
	sizes = other.sizes;
	twoints = other.twoints;
}

IntegralEngine::~IntegralEngine() {

}

// Accessors

// Return estimates of the memory that will be needed by the 
// one and two electron integrals. Returns as:
// [1e cart, 2e cart, 1e spher, 2e spher]
Vector IntegralEngine::getEstimates() const
{
	Vector estimates(4);
	// The amount of memory is roughly the number of integrals times the size
	// of a double in memory
	double TOMB = 1.0/(1024.0 * 1024.0);
	estimates[0] = TOMB*sizeof(double)*sizes(0);
	estimates[1] = TOMB*sizeof(double)*sizes(1);
	estimates[2] = TOMB*sizeof(double)*sizes(2);
	estimates[3] = TOMB*sizeof(double)*sizes(3);
	return estimates;
}

// Return a particular integral from twoints (taking into account symmetries)
double IntegralEngine::getERI(int i, int j, int k, int l) const
{
	return twoints(i, j, k, l);
} 

// Print a sorted list of ERIs to ostream output
void IntegralEngine::printERI(std::ostream& output, int NSpher) const
{
	output << " TWO ELECTRON INTEGRALS: \n";
	// print them out
	int icount = 0;
	int scount = 0;
	for (int c1 = 0; c1 < NSpher; c1++){
		for (int c2 = 0; c2 < c1+1; c2++){
			for (int c3 = 0; c3 < c1+1; c3++){
				for(int c4 = 0; c4 < c3+1; c4++){
					icount++;
					double multiplier = 0.125;
					if (c1!=c2) { multiplier *= 2.0; }
					if (c3!=c4) { multiplier *= 2.0; }
					if ((c1+c2) != (c3+c4)) { multiplier*=2.0; }
					else if ((c1!=c3) && (c1 != c4)) { multiplier *=2.0; }
					else if ((c2!=c3) && (c2 != c4)) { multiplier *=2.0; }
					if (fabs(getERI(c4, c3, c2, c1)) < molecule->control->get_option<double>("thrint")) { scount++; multiplier = 0; }
					output << std::setw(6) << c1+1;
					output << std::setw(6) << c2+1;
					output << std::setw(6) << c3+1;
					output << std::setw(6) << c4+1;
					output << std::setw(20) << multiplier*getERI(c4, c3, c2, c1);
					output << "\n";
				}
			}
		}
	}			
	output << "N 2e Ints: " << icount << "\n";
	output << "N insig. Ints: " << scount << "\n";
}

// Contract a set of 1e- integrals
// Assumes that integrals are ordered as: 00, 01, 02, ..., 10, 11, 12, ...,
// and so on, where the first number refers to the index of c1, and the second, 
// that of c2.
double IntegralEngine::makeContracted(Vector& c1, Vector& c2, Vector& ints) const
{
	double integral = 0.0;
	int N1 = c1.size();
	int N2 = c2.size();
	// Loop over contraction coefficients
	for (int i = 0; i < N1; i++){
		for (int j = 0; j < N2; j++){
			integral += c1(i)*c2(j)*ints(i*N2+j);
		}
	}
	return integral;
}

// Sphericalise a matrix of 1e- integrals (ints)
// where the cols have angular momenta lnums.
// Returns matrix of integrals in canonical order
Matrix IntegralEngine::makeSpherical(const Matrix& ints) const
{
	// Calculate the size of matrix needed
	return transmat * ints * transmat.transpose(); 
}


void IntegralEngine::buildTransMat() 
{
	auto &shells = molecule->getBasis().getIntShells();
	
	int natoms = molecule->getNAtoms();
	int ncar = ncart(shells); // No. of cartesian basis functions
	int nspher = nbasis(shells); // No. of spherical basis functions
	
	transmat = Matrix::Zero(nspher, ncar);
	int row = 0; int col_offset = 0; 
	for (auto s : shells) {
		int lam = s.contr[0].l; 
		switch(lam) {
			case 0: { // s-type
				transmat(row++, col_offset) = 1.0;
				break; 
			}
			case 1: { // l-type 
				transmat(row++, col_offset+1) = 1.0;
				transmat(row++, col_offset+2) = 1.0;
				transmat(row++, col_offset) = 1.0;
				break; 
			}
			case 2: { // d-type
				transmat(row++, col_offset + 1) = std::sqrt(3.0); 
				transmat(row++, col_offset + 4) = std::sqrt(3.0);
				transmat(row, col_offset) = -0.5;
				transmat(row, col_offset + 3) = -0.5;
				transmat(row++, col_offset + 5) = 1.0;
				transmat(row++, col_offset + 2) = std::sqrt(3.0);
				transmat(row, col_offset) = 0.5*std::sqrt(3.0);
				transmat(row++, col_offset + 3) = -0.5*std::sqrt(3.0); 
				break;
			}
			case 3: { // f-type
				transmat(row, col_offset+1) = 0.75*std::sqrt(10.0);
				transmat(row++, col_offset+6) = -0.25*std::sqrt(10.0);
				transmat(row++, col_offset+4) = std::sqrt(15.0);
				transmat(row, col_offset+1) = -0.25*sqrt(6.0); 
				transmat(row, col_offset + 6) = -0.25*sqrt(6.0);
				transmat(row++, col_offset + 8) = std::sqrt(6.0);
				transmat(row, col_offset+2) = -1.5;
				transmat(row, col_offset+7) = -1.5;
				transmat(row++, col_offset+9) = 1.0;
				transmat(row, col_offset) = -0.25*std::sqrt(6.0);
				transmat(row, col_offset+3) = -0.25*std::sqrt(6.0);
				transmat(row++, col_offset+5) = std::sqrt(6.0); 
				transmat(row, col_offset+2) = 0.5*std::sqrt(15.0);
				transmat(row++, col_offset+7) = -0.5*std::sqrt(15.0);
				transmat(row, col_offset) = 0.25 * std::sqrt(10.0);
				transmat(row++, col_offset+3) = -0.75*std::sqrt(10.0); 
				break;
			}
			case 4: { // g-type
				break;
			}
			case 5: { // h-type
				break;
			}
			
		}
		col_offset += (lam+1)*(lam+2)/2;
	}
	
}

size_t IntegralEngine::nbasis(const std::vector<libint2::Shell>& shells) {
	size_t n = 0;
	for (const auto& shell: shells)
		n += shell.size();
	return n;
}

size_t IntegralEngine::ncart(const std::vector<libint2::Shell>& shells) {
	size_t n = 0;
	for (const auto& shell: shells)
		n += shell.cartesian_size();
	return n;
}

size_t IntegralEngine::max_nprim(const std::vector<libint2::Shell>& shells) {
	size_t n = 0;
	for (auto shell: shells)
		n = std::max(shell.nprim(), n);
	return n;
}

int IntegralEngine::max_l(const std::vector<libint2::Shell>& shells) {
	int l = 0;
	for (auto shell: shells)
		for (auto c: shell.contr)
			l = std::max(c.l, l);
	return l;
}

std::vector<size_t> IntegralEngine::map_shell_to_basis_function(const std::vector<libint2::Shell>& shells) {
	std::vector<size_t> result;
	result.reserve(shells.size());

	size_t n = 0;
	for (auto shell: shells) {
		result.push_back(n);
		n += shell.size();
	}

	return result;
}

std::vector<size_t> IntegralEngine::map_shell_to_cart_basis_function(const std::vector<libint2::Shell>& shells) {
	std::vector<size_t> result;
	result.reserve(shells.size());

	size_t n = 0;
	for (auto shell: shells) {
		result.push_back(n);
		n += shell.cartesian_size();
	}

	return result;
}

std::vector<long> IntegralEngine::map_shell_to_atom(const std::vector<Atom>& atoms, const std::vector<libint2::Shell>& shells) {
	std::vector<long> result;
	result.reserve(shells.size());
	for(const auto& s: shells) {
		auto a = std::find_if(atoms.begin(), atoms.end(), [&s](const Atom& a){ return s.O[0] == a.getX() && s.O[1] == a.getY() && s.O[2] == a.getZ(); } );
		result.push_back( a != atoms.end() ? a - atoms.begin() : -1);
	}
	return result;
}

Matrix IntegralEngine::compute_1body_ints(const std::vector<libint2::Shell>& shells,
libint2::Operator obtype,
const std::vector<Atom>& atoms)
{
	using libint2::Shell;
	using libint2::Engine;
	using libint2::Operator;
	
	int nthreads = molecule->control->get_option<int>("nthreads"); 

	const auto n = nbasis(shells);
	EMatrix result(n, n); 
	int nshells = shells.size(); 

	std::vector<libint2::Engine> engines(nthreads);
	engines[0] = libint2::Engine(obtype, max_nprim(shells), max_l(shells), 0);
	// pass operator params to the engine, e.g.
	// nuclear attraction ints engine needs to know where the charges sit ...
	// the nuclei are charges in this case; in QM/MM there will also be classical
	// charges
	if (obtype == Operator::nuclear) {
		std::vector<std::pair<double,std::array<double,3>>> q;
		for(const auto& atom : atoms) {
			q.push_back( {static_cast<double>(atom.getEffectiveCharge()), {{atom.getX(), atom.getY(), atom.getZ()}}} );
		}
			
		engines[0].set_params(q);
	}
	
	for (size_t i = 1; i != nthreads; ++i) {
		engines[i] = engines[0];
	}                                                                                              

	auto shell2bf = map_shell_to_basis_function(shells);

	// buf[0] points to the target shell set after every call  to engine.compute()                                                                                                    
	
	auto compute = [&](int thread_id) {

		auto& engine = engines[thread_id]; 
		const auto& buf = engine.results();

		// loop over unique shell pairs, {s1,s2} such that s1 >= s2
		// this is due to the permutational symmetry of the real integrals over
		// Hermitian operators: (1|2) = (2|1)
		for(auto s1=0, s12=0; s1!=shells.size(); ++s1) {

			auto bf1 = shell2bf[s1]; // first basis function in this shell                                                                                                                           
			auto n1 = shells[s1].size();

			for(auto s2=0; s2<=s1; ++s2, ++s12) {

				if (s12 % nthreads != thread_id) continue; 
				
				auto bf2 = shell2bf[s2];
				auto n2 = shells[s2].size();

				// compute shell pair; return is the pointer to the buffer                                                                                                                             
				engine.compute(shells[s1], shells[s2]);

				// "map" buffer to a const Eigen Matrix, and copy it to the corresponding blocks of the result                                                                                         
				Eigen::Map<const EMatrix> buf_mat(buf[0], n1, n2);
				result.block(bf1, bf2, n1, n2) = buf_mat;
				if (s1 != s2) // if s1 >= s2, copy {s1,s2} to the corresponding {s2,s1} block, note the transpose!                                                                                     
					result.block(bf2, bf1, n2, n1) = buf_mat.transpose();

			}
		}
	};  // compute lambda

	parallel_do(compute, nthreads);

	return result;                                                                     
}

EMatrix IntegralEngine::compute_shellblock_norm(const std::vector<libint2::Shell> &shells, const Matrix& A) {
	const auto nsh = shells.size();
	EMatrix Ash(nsh, nsh);
	EMatrix EA = A;
	
	auto shell2bf = map_shell_to_basis_function(shells);
	for (size_t s1 = 0; s1 != nsh; ++s1) {
		const auto& s1_first = shell2bf[s1];
		const auto& s1_size = shells[s1].size();
		for (size_t s2 = 0; s2 != nsh; ++s2) {
			const auto& s2_first = shell2bf[s2];
			const auto& s2_size = shells[s2].size();

			Ash(s1, s2) = EA.block(s1_first, s2_first, s1_size, s2_size).lpNorm<Eigen::Infinity>();
		}
	}

	return Ash;
}

template <libint2::Operator Kernel>
Matrix IntegralEngine::compute_schwarz_ints( const std::vector<libint2::Shell> &bs1, const std::vector<libint2::Shell> _bs2, 
bool use_2norm, typename libint2::operator_traits<Kernel>::oper_params_type params) {
	const std::vector<libint2::Shell>& bs2 = (_bs2.empty() ? bs1 : _bs2);
	const auto nsh1 = bs1.size();
	const auto nsh2 = bs2.size();
	const auto bs1_equiv_bs2 = (&bs1 == &bs2);

	Matrix K = Matrix::Zero(nsh1, nsh2);

	// construct the 2-electron repulsion integrals engine
	using libint2::Engine;
	int nthreads = molecule->control->get_option<int>("nthreads"); 

	std::vector<Engine> engines(nthreads);
	// !!! very important: cannot screen primitives in Schwarz computation !!!
	auto epsilon = 0.;
	engines[0] = Engine(Kernel, std::max(max_nprim(bs1), max_nprim(bs2)),
	std::max(max_l(bs1), max_l(bs2)), 0, epsilon, params);
	for (size_t i = 1; i != nthreads; ++i) {
		engines[i] = engines[0];
	}
	
	auto compute = [&](int thread_id) {

		const auto& buf = engines[thread_id].results();

		// loop over permutationally-unique set of shells
		for (auto s1 = 0l, s12 = 0l; s1 != nsh1; ++s1) {
			auto n1 = bs1[s1].size();  // number of basis functions in this shell

			auto s2_max = bs1_equiv_bs2 ? s1 : nsh2 - 1;
			for (auto s2 = 0; s2 <= s2_max; ++s2, ++s12) {
				if (s12 % nthreads != thread_id) continue;

				auto n2 = bs2[s2].size();
				auto n12 = n1 * n2;

				engines[thread_id].compute2<Kernel, libint2::BraKet::xx_xx, 0>(bs1[s1], bs2[s2],
				bs1[s1], bs2[s2]);
				assert(buf[0] != nullptr &&
					"to compute Schwarz ints turn off primitive screening");

				// to apply Schwarz inequality to individual integrals must use the diagonal elements
				// to apply it to sets of functions (e.g. shells) use the whole shell-set of ints here
				Eigen::Map<const Matrix> buf_mat(buf[0], n12, n12);
				auto norm2 = use_2norm ? buf_mat.norm()
					: buf_mat.lpNorm<Eigen::Infinity>();
				K(s1, s2) = std::sqrt(norm2);
				if (bs1_equiv_bs2) K(s2, s1) = K(s1, s2);
			}
		}
	};  // thread lambda

	parallel_do(compute, nthreads);
	   
	return K;
}
 
shellpair_list_t IntegralEngine::compute_shellpair_list(const std::vector<libint2::Shell>& bs1, 
const std::vector<libint2::Shell> _bs2, double threshold) {
 		
	const std::vector<libint2::Shell> &bs2 = (_bs2.empty() ? bs1 : _bs2);
	const auto nsh1 = bs1.size();
	const auto nsh2 = bs2.size();
	const auto bs1_equiv_bs2 = (&bs1 == &bs2);

	// construct the 2-electron repulsion integrals engine
	using libint2::Engine;
	using libint2::Operator;
	Engine engine(Operator::overlap, std::max(max_nprim(bs1), max_nprim(bs2)),
	std::max(max_l(bs1), max_l(bs2)), 0);

	shellpair_list_t result;
	const auto& buf = engine.results();

	// loop over permutationally-unique set of shells
	for (auto s1 = 0l, s12 = 0l; s1 != nsh1; ++s1) {
		if (result.find(s1) == result.end())
			result.insert(std::make_pair(s1, std::vector<size_t>()));

		auto n1 = bs1[s1].size();  // number of basis functions in this shell

		auto s2_max = bs1_equiv_bs2 ? s1 : nsh2 - 1;
		for (auto s2 = 0; s2 <= s2_max; ++s2, ++s12) {
	
			auto on_same_center = (bs1[s1].O == bs2[s2].O);
			bool significant = on_same_center;
			if (not on_same_center) {
				auto n2 = bs2[s2].size();
				engine.compute(bs1[s1], bs2[s2]);
				Eigen::Map<const EMatrix> buf_mat(buf[0], n1, n2);
				auto norm = buf_mat.norm();
				significant = (norm >= threshold);
			}

			if (significant) {
				result[s1].emplace_back(s2);
			}
		}
	}

	// resort shell list in increasing order, i.e. result[s][s1] < result[s][s2]
	// if s1 < s2
	for (auto s1 = 0l; s1 != nsh1; ++s1) {
		auto& list = result[s1];
		std::sort(list.begin(), list.end());
	}

	return result;
}
 
Matrix IntegralEngine::compute_ecp_ints(const std::vector<libint2::Shell>& shells, int deriv_order) {
	const auto n = ncart(shells);
	Matrix ecps = Matrix::Zero(n, n);
	
	// Initialise ecp integral engine
	molecule->control->log.print("\nIntialising ECP integral calculations...\n");
	libecpint::ECPIntegral ecpint(molecule->getBasis().getMaxL(), molecule->getECPBasis().getMaxL(), deriv_order);
	molecule->control->log.localTime();
	
	libecpint::ECPBasis& ecpset = molecule->getECPBasis();
	auto shell2bf = map_shell_to_cart_basis_function(shells);
	
	// loop over shells
	for(auto s1=0; s1!=shells.size(); ++s1) {

		auto bf1 = shell2bf[s1];
		auto n1 = shells[s1].size();
		
		double A[3] = { shells[s1].O[0], shells[s1].O[1], shells[s1].O[2] };
		libecpint::GaussianShell shellA(A, shells[s1].contr[0].l);
		for (auto c : shells[s1].contr)
			for (int i = 0; i < c.coeff.size(); i++)
				shellA.addPrim(shells[s1].alpha[i], c.coeff[i]);

		for(auto s2=0; s2<=s1; ++s2) {

			auto bf2 = shell2bf[s2];
			auto n2 = shells[s2].size();
			
			double B[3] = { shells[s2].O[0], shells[s2].O[1], shells[s2].O[2] };
			libecpint::GaussianShell shellB(B, shells[s2].contr[0].l); 
			for (auto c : shells[s2].contr)
				for (int i = 0; i < c.coeff.size(); i++)
					shellB.addPrim(shells[s2].alpha[i], c.coeff[i]);
			
			libecpint::TwoIndex<double> tempValues;
			Matrix shellPairInts = Matrix::Zero(shellA.ncartesian(), shellB.ncartesian());
			for (int i = 0; i < ecpset.getN(); i++) {
				ecpint.compute_shell_pair(ecpset.getECP(i), shellA, shellB, tempValues);
				for (int na = 0; na < shellA.ncartesian(); na++) {
					for (int nb = 0; nb < shellB.ncartesian(); nb++) shellPairInts(na, nb) += tempValues(na, nb);
				} 
			}
			
			//shellPairInts.print(); std::cout << "\n\n";
			for (int i = bf1; i < bf1 + shellPairInts.rows(); i++) {
				for (int j = bf2; j < bf2 + shellPairInts.cols(); j++) {
					ecps(i, j) = shellPairInts(i-bf1, j-bf2);
					ecps(j, i) = ecps(i, j);
				}
			}
			
			bf2 += n2; 
		}
		
		bf1 += n1; 
	}
	
	//std::cout << makeSpherical(ecps) << std::endl; 
	
	return makeSpherical(ecps);
}
 
