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

#ifndef INTEGRALSHEADERDEF
#define INTEGRALSHEADERDEF

// Includes
#include "eigen_wrapper.hpp"
#include "molecule.hpp"
#include <iostream>
#include "tensor4.hpp"
#include <libint2.hpp>
#include <vector>
#include <array>
#include <unordered_map>
#include <iostream>

using shellpair_list_t = std::unordered_map<size_t, std::vector<size_t>>;

// Declare forward dependencies
class Atom;
class Tensor6;

#define PAIR_INDEX(x, y) (x < y ? (y*(y+1))/2 + x : (x*(x+1))/2 + y)
#define TRIPLE_INDEX(x, y, k, n) (x < y ? (k*n + (y*(y+1))/2 + x) : (k*n + (x*(x+1))/2 + y))

class DFBlocks {
private:
	std::vector<Matrix> eris;
	std::vector<bool> zeroList; 
	
	int nfrags, npairs, ntriples; 

public:
	DFBlocks() : nfrags(-1), npairs(-1), ntriples(-1) { } 
	DFBlocks(int nfrags); 
	 
	Matrix& operator()(int f1, int f2, int fk);
	bool isZero(int f1, int f2, int fk) const;
	void setZero(int f1, int f2, int fk); 
	void setNonZero(int f1, int f2, int fk);  
	
	int getNFrags() const { return nfrags; }
	
	void resize(int nfrags); 
}; 

//Begin class declaration
class IntegralEngine
{
private:
	SharedMolecule molecule;
	Matrix sints;
	Matrix tints;
	Matrix naints;
	Matrix prescreen;
	Matrix transmat; 
	iVector sizes;
	S8EvenTensor4 twoints;
public:
	IntegralEngine(SharedMolecule m, bool print = true); //Constructor
	IntegralEngine(SharedMolecule m, const IntegralEngine& i, int start, int finish); 
	IntegralEngine(const IntegralEngine& other);
	~IntegralEngine();

	// Accessors
	Vector getEstimates() const;
	double getOverlap(int i, int j) const { return sints(i, j); }
	Matrix& getOverlap() { return sints; }
	Matrix getOverlap() const { return sints; }
	double getKinetic(int i, int j) const { return tints(i, j); }
	Matrix& getKinetic() { return tints; }
	Matrix getKinetic() const { return tints; }
	double getNucAttract(int i, int j) const { return naints(i, j); }
	Matrix& getNucAttract() { return naints; }
	Matrix getNucAttract() const { return naints; }
	Matrix& getTransMat() { return transmat; }
	double getERI(int i, int j, int k, int l) const;
	S8EvenTensor4 getERI() const { return twoints; }
	Matrix& getPrescreen() { return prescreen; }
	void clearTwoInts() { twoints.resize(0); }

	// Intrinsic routines
	void printERI(std::ostream& output, int NSpher) const;
	double makeContracted(Vector& c1, Vector& c2, Vector& ints) const;
	Matrix makeSpherical(const Matrix& ints) const;
	void buildTransMat(); 

	size_t nbasis(const std::vector<libint2::Shell>& shells);
	size_t ncart(const std::vector<libint2::Shell>& shells);
	size_t max_nprim(const std::vector<libint2::Shell>& shells);
	int max_l(const std::vector<libint2::Shell>& shells);
	std::vector<size_t> map_shell_to_basis_function(const std::vector<libint2::Shell>& shells);
	std::vector<size_t> map_shell_to_cart_basis_function(const std::vector<libint2::Shell>& shells);
	std::vector<long> map_shell_to_atom(const std::vector<Atom>& atoms, const std::vector<libint2::Shell>& shells);
  
	shellpair_list_t compute_shellpair_list(const std::vector<libint2::Shell>& bs1, 
	const std::vector<libint2::Shell> _bs2 = std::vector<libint2::Shell>(), double threshold = 1e-12);
  
	EMatrix compute_shellblock_norm(const std::vector<libint2::Shell> &shells, const Matrix& A);
  
	template <libint2::Operator Kernel = libint2::Operator::coulomb>
	Matrix compute_schwarz_ints( const std::vector<libint2::Shell> &_bs1, 
	const std::vector<libint2::Shell> _bs2 = std::vector<libint2::Shell>(), bool use_2norm = false,
	typename libint2::operator_traits<Kernel>::oper_params_type params = libint2::operator_traits<Kernel>::default_params());

	Matrix compute_1body_ints(const std::vector<libint2::Shell>& shells,
	libint2::Operator t,
	const std::vector<Atom>& atoms = std::vector<Atom>());
	
	Matrix compute_ecp_ints(const std::vector<libint2::Shell>& shells, int deriv_order = 0);
	

	S8EvenTensor4 compute_eris(const std::vector<libint2::Shell>& shells);
	void compute_eris_3index(const std::vector<libint2::Shell>& obs, const std::vector<libint2::Shell>& auxbs, Matrix& eris);
	void compute_eris_3index(const std::vector<libint2::Shell>& obs, const std::vector<libint2::Shell>& auxbs, SparseMatrix& eris);
	void compute_eris_3index(const std::vector<libint2::Shell>& obs1, const std::vector<libint2::Shell>& obs2,
		 const std::vector<libint2::Shell>& auxbs, Matrix& eris);
	Matrix compute_eris_2index(const std::vector<libint2::Shell>& auxbs);
  
	template<libint2::Operator obtype>
	std::vector<EMatrix> compute_1body_ints_deriv(unsigned deriv_order, const std::vector<libint2::Shell>& obs, const std::vector<Atom> &atoms)
	{
		const auto n = nbasis(obs);
		const auto nshells = obs.size();
		constexpr auto nopers = libint2::operator_traits<obtype>::nopers;
		const auto nresults = nopers * libint2::num_geometrical_derivatives(atoms.size(), deriv_order);
		typedef std::vector<EMatrix> result_type;
		result_type result(nresults);
		for (auto& r : result) r = EMatrix::Zero(n, n);

		using libint2::Engine;
		using libint2::Operator; 
	
		// construct the 1-body integrals engine
		Engine engine(obtype, max_nprim(obs), max_l(obs), deriv_order);

		// nuclear attraction ints engine needs to know where the charges sit ...
		// the nuclei are charges in this case; in QM/MM there will also be classical
		// charges
		if (obtype == Operator::nuclear) {
			std::vector<std::pair<double,std::array<double,3>>> q;
			for(const auto& atom : atoms) {
				q.push_back( {static_cast<double>(atom.getEffectiveCharge()), {{atom.getX(), atom.getY(), atom.getZ()}}} );
			}
			engine.set_params(q);
		}
	 
		auto shell2bf = map_shell_to_basis_function(obs);
		auto shell2atom = map_shell_to_atom(atoms, obs);
		
		const auto natoms = atoms.size();
		const auto two_times_ncoords = 6*natoms;
		const auto nderivcenters_shset = 2 + ((obtype == Operator::nuclear) ? natoms : 0);

		const auto& buf = engine.results();
		shellpair_list_t obs_shellpair_list = compute_shellpair_list(obs);
		
		// loop over unique shell pairs, {s1,s2} such that s1 >= s2
		// this is due to the permutational symmetry of the real integrals over
		// Hermitian operators: (1|2) = (2|1)
		for (auto s1 = 0l, s12 = 0l; s1 != nshells; ++s1) {
			auto bf1 = shell2bf[s1];  // first basis function in this shell
			auto n1 = obs[s1].size();
			auto atom1 = shell2atom[s1];
			assert(atom1 != -1);

			auto s1_offset = s1 * (s1+1) / 2;
			for (auto s2: obs_shellpair_list[s1]) {
				auto s12 = s1_offset + s2;

				auto bf2 = shell2bf[s2];
				auto n2 = obs[s2].size();
				auto atom2 = shell2atom[s2];

				auto n12 = n1 * n2;

				// compute shell pair; return is the pointer to the buffer
				engine.compute(obs[s1], obs[s2]);

				// "copy" lambda copies shell set \c idx to the operator matrix with
				// index \c op
				auto add_shellset_to_dest = [&](std::size_t op, std::size_t idx,
				double scale = 1.0) {
					// "map" buffer to a const Eigen Matrix, and copy it to the
					// corresponding blocks of the result
					Eigen::Map<const EMatrix> buf_mat(buf[idx], n1, n2);
					if (scale == 1.0)
						result[op].block(bf1, bf2, n1, n2) += buf_mat;
					else
						result[op].block(bf1, bf2, n1, n2) += scale * buf_mat;
					if (s1 != s2) {  // if s1 >= s2, copy {s1,s2} to the corresponding
					// {s2,s1} block, note the transpose!
					if (scale == 1.0)
						result[op].block(bf2, bf1, n2, n1) += buf_mat.transpose();
					else
						result[op].block(bf2, bf1, n2, n1) += scale * buf_mat.transpose();
				}
			};

			switch (deriv_order) {
				case 0: {
					for (std::size_t op = 0; op != nopers; ++op) {
						add_shellset_to_dest(op, op);
					}
				} break;

				// map deriv quanta for this shell pair to the overall deriv quanta
				//
				// easiest to explain with example:
				// in sto-3g water shells 0 1 2 sit on atom 0, shells 3 and 4 on atoms
				// 1 and 2 respectively
				// each call to engine::compute for nuclear ints will return
				// derivatives
				// with respect to 15 coordinates, obtained as 3 (x,y,z) times 2 + 3 =
				// 5 centers
				// (2 centers on which shells sit + 3 nuclear charges)
				// (for overlap, kinetic, and emultipole ints we there are only 6
				// coordinates
				//  since the operator is coordinate-independent, or derivatives with
				//  respect to
				//  the operator coordinates are not computed)
				//

				case 1: {
					std::size_t shellset_idx = 0;
					for (auto c = 0; c != nderivcenters_shset; ++c) {
						auto atom = (c == 0) ? atom1 : ((c == 1) ? atom2 : c - 2);
						auto op_start = 3 * atom * nopers;
						auto op_fence = op_start + nopers;
						for (auto xyz = 0; xyz != 3; ++xyz, op_start += nopers, op_fence += nopers) {
							for (unsigned int op = op_start; op != op_fence; ++op, ++shellset_idx) {
								add_shellset_to_dest(op, shellset_idx);
							}
						}
					}
				} break;

				case 2: {
					//
					// must pay attention to symmetry when computing 2nd and higher-order derivs
					// e.g. d2 (s1|s2) / dX dY involves several cases:
					// 1. only s1 (or only s2) depends on X AND Y (i.e. X and Y refer to same atom) =>
					//    d2 (s1|s2) / dX dY = (d2 s1 / dX dY | s2)
					// 2. s1 depends on X only, s2 depends on Y only (or vice versa) =>
					//    d2 (s1|s2) / dX dY = (d s1 / dX | d s2 / dY)
					// 3. s1 AND s2 depend on X AND Y (i.e. X and Y refer to same atom) =>
					//    case A: X != Y
					//    d2 (s1|s2) / dX dY = (d2 s1 / dX dY | s2) + (d s1 / dX | d s2 / dY)
					//      + (d s1 / dY | d s2 / dX) + (s1| d2 s2 / dX dY )
					//    case B: X == Y
					//    d2 (s1|s2) / dX2 = (d2 s1 / dX2 | s2) + 2 (d s1 / dX | d s2 / dX)
					//      + (s1| d2 s2 / dX2 )

					// computes upper triangle index
					// n2 = matrix size times 2
					// i,j = (unordered) indices
#define upper_triangle_index(n2, i, j) (std::min((i), (j))) * ((n2) - (std::min((i), (j))) - 1) / 2 + (std::max((i), (j)))

					// look over shellsets in the order in which they appear
					std::size_t shellset_idx = 0;
					for (auto c1 = 0; c1 != nderivcenters_shset; ++c1) {
						auto a1 = (c1 == 0) ? atom1 : ((c1 == 1) ? atom2 : c1 - 2);
						auto coord1 = 3 * a1;
						for (auto xyz1 = 0; xyz1 != 3; ++xyz1, ++coord1) {

							for (auto c2 = c1; c2 != nderivcenters_shset; ++c2) {
								auto a2 = (c2 == 0) ? atom1 : ((c2 == 1) ? atom2 : c2 - 2);
								auto xyz2_start = (c1 == c2) ? xyz1 : 0;
								auto coord2 = 3 * a2 + xyz2_start;
								for (auto xyz2 = xyz2_start; xyz2 != 3; ++xyz2, ++coord2) {

									double scale = (coord1 == coord2 && c1 != c2) ? 2.0 : 1.0;

									const auto coord12 =
										upper_triangle_index(two_times_ncoords, coord1, coord2);
									auto op_start = coord12 * nopers;
									auto op_fence = op_start + nopers;
									for (auto op = op_start; op != op_fence; ++op, ++shellset_idx) {
										add_shellset_to_dest(op, shellset_idx, scale);
									}
								}
							}
						}
					}
				} break;
#undef upper_triangle_index

				default: {
					assert(false && "not yet implemented");

					using ShellSetDerivIterator = libint2::FixedOrderedIntegerPartitionIterator<std::vector<unsigned int>>;
					ShellSetDerivIterator shellset_diter(deriv_order, nderivcenters_shset);
					while (shellset_diter) {
						const auto& deriv = *shellset_diter;
					}
				}
			}  // copy shell block switch
		}
	}  // s2 <= s1
	return result;
}
	  

};

#endif  
    
