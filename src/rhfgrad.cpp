#include "fock.hpp"
#include "error.hpp"
#include "atom.hpp"
#include "basis.hpp"
#include "tensor4.hpp"
#include <iostream>
#include <iomanip>
#include <cmath>
#include <libint2.hpp>
#include "logger.hpp"
#include "ProgramController.hpp"
#include "mp2.hpp"
#include "cphf.hpp"

void Fock::compute_forces(const std::vector<Atom> &atoms, int nocc) {
	using libint2::Shell;
	using libint2::Operator;

	std::vector<Shell>& shells = molecule->getBasis().getIntShells();
	
	Matrix F1 = Matrix::Zero(atoms.size(), 3);
	Matrix F_Pulay = Matrix::Zero(atoms.size(), 3);
	
	// One-body contributions to forces
	auto T1 = integrals.compute_1body_ints_deriv<Operator::kinetic>(1, shells, atoms);
	auto V1 = integrals.compute_1body_ints_deriv<Operator::nuclear>(1, shells, atoms);
	for (auto atom = 0, i = 0; atom != atoms.size(); ++atom) {
		for (auto xyz = 0; xyz != 3; ++xyz, ++i) {
			auto force =  (T1[i] + V1[i]).cwiseProduct(dens).sum();
			F1(atom, xyz) += force;
		}
	}
	
	// Pulay force
	EMatrix evals_occ = EMatrix::Zero(nocc, nocc);
	for (int i = 0; i < nocc; ++i) evals_occ(i, i) = eps[i];
	EMatrix C_occ = CP.block(0, 0, dens.rows(), nocc); 
	EMatrix W = C_occ * evals_occ * C_occ.transpose();

	auto S1 = integrals.compute_1body_ints_deriv<Operator::overlap>(1, shells, atoms);
	for (auto atom = 0, i = 0; atom != atoms.size(); ++atom) {
		for (auto xyz = 0; xyz != 3; ++xyz, ++i) {
			auto force = 2 * S1[i].cwiseProduct(W).sum();
			F_Pulay(atom, xyz) -= force;
		}
	}
	
	// Two-body contrtibutions to forces
	Matrix F2 = Matrix::Zero(atoms.size(), 3);
	auto G1 = compute_2body_fock_deriv<1>(atoms, dens);
	for (auto atom = 0, i = 0; atom != atoms.size(); ++atom) {
		for (auto xyz = 0; xyz != 3; ++xyz, ++i) {
			// identity prefactor since E(HF) = trace(H + F, D) = trace(2H + G, D)
			auto force = 0.25 * G1[i].cwiseProduct(dens).sum();
			F2(atom, xyz) += force;
		}
	}
	
	// Compute nuclear repulsion forces
	Matrix FN = Matrix::Zero(atoms.size(), 3);
	for (auto a1 = 1; a1 != atoms.size(); ++a1) {
		const auto& atom1 = atoms[a1];
		for (auto a2 = 0; a2 < a1; ++a2) {
			const auto& atom2 = atoms[a2];

			auto x12 = atom1.getX() - atom2.getX();
			auto y12 = atom1.getY() - atom2.getY();
			auto z12 = atom1.getZ() - atom2.getZ();
			auto r12_2 = x12 * x12 + y12 * y12 + z12 * z12;
			auto r12 = sqrt(r12_2);
			auto r12_3 = r12 * r12_2;

			auto z1z2_over_r12_3 =
				atom1.getEffectiveCharge() * atom2.getEffectiveCharge() / r12_3;

			auto fx = -x12 * z1z2_over_r12_3;
			auto fy = -y12 * z1z2_over_r12_3;
			auto fz = -z12 * z1z2_over_r12_3;
			FN(a1, 0) += fx;
			FN(a1, 1) += fy;
			FN(a1, 2) += fz;
			FN(a2, 0) -= fx;
			FN(a2, 1) -= fy;
			FN(a2, 2) -= fz;
		}
	}
	
	forces = F1 + F_Pulay + F2 + FN;
}

void Fock::compute_hessian(const std::vector<Atom> &atoms, int nocc) {
	
	dens /= 2.0; 
	
	const auto N = atoms.size();
	const auto ncoords = 3 * N;
	const auto nelem = ncoords * (ncoords +1) / 2;
	
	using libint2::Shell;
	using libint2::Operator;

	std::vector<Shell>& shells = molecule->getBasis().getIntShells();
	
	// Calculate W
	
	EMatrix evals_occ = EMatrix::Zero(nocc, nocc);
	for (int i = 0; i < nocc; ++i) evals_occ(i, i) = eps[i];
	EMatrix C_occ = CP.block(0, 0, dens.rows(), nocc); 
	EMatrix W = C_occ * evals_occ * C_occ.transpose();
	
	// Calculate first derivative matrices
	
	auto T1 = integrals.compute_1body_ints_deriv<Operator::kinetic>(1, shells, atoms);
	auto V1 = integrals.compute_1body_ints_deriv<Operator::nuclear>(1, shells, atoms);
	auto S1 = integrals.compute_1body_ints_deriv<Operator::overlap>(1, shells, atoms);
	auto G1 = compute_2body_fock_deriv<1>(atoms, dens);
	
	// Determine contributions to forces
	forces = Matrix::Zero(N, 3);
	for (auto atom = 0, i = 0; atom != N; ++atom) {
		for (auto xyz = 0; xyz != 3; ++xyz, ++i) {
			auto h_force = 2.0*(T1[i] + V1[i]).cwiseProduct(dens).sum();
			auto w_force = 2.0*S1[i].cwiseProduct(W).sum();
			auto g_force = G1[i].cwiseProduct(dens).sum();
			forces(atom, xyz) += h_force + g_force - w_force; 
		}
	}
	
	// Nuclear repulsion derivatives
	Matrix FN = Matrix::Zero(N, 3);
	for (auto a1 = 1; a1 != N; ++a1) {
		const auto& atom1 = atoms[a1];
		for (auto a2 = 0; a2 < a1; ++a2) {
			const auto& atom2 = atoms[a2];

			auto x12 = atom1.getX() - atom2.getX();
			auto y12 = atom1.getY() - atom2.getY();
			auto z12 = atom1.getZ() - atom2.getZ();
			auto r12_2 = x12 * x12 + y12 * y12 + z12 * z12;
			auto r12 = sqrt(r12_2);
			auto r12_3 = r12 * r12_2;

			auto z1z2_over_r12_3 =
				atom1.getEffectiveCharge() * atom2.getEffectiveCharge() / r12_3;

			auto fx = -x12 * z1z2_over_r12_3;
			auto fy = -y12 * z1z2_over_r12_3;
			auto fz = -z12 * z1z2_over_r12_3;
			FN(a1, 0) += fx;
			FN(a1, 1) += fy;
			FN(a1, 2) += fz;
			FN(a2, 0) -= fx;
			FN(a2, 1) -= fy;
			FN(a2, 2) -= fz;
		}
	}
	
	forces += FN; 
	
	// Calculate second derivative matrices
	auto T2 = integrals.compute_1body_ints_deriv<Operator::kinetic>(2, shells, atoms);
	auto V2 = integrals.compute_1body_ints_deriv<Operator::nuclear>(2, shells, atoms);
	auto S2 = integrals.compute_1body_ints_deriv<Operator::overlap>(2, shells, atoms);
	auto G2 = compute_2body_fock_deriv<2>(atoms, dens);
	
	// Determine contribution to non-response part of the Hessian
	
	Matrix H_nr = Matrix::Zero(ncoords, ncoords); 
	for (auto row = 0, i = 0; row != ncoords; ++row) {
		for (auto col = row; col != ncoords; ++col, ++i) {
			auto h_hess = 2.0 * (T2[i] + V2[i]).cwiseProduct(dens).sum();
			auto w_hess = 2.0 * S2[i].cwiseProduct(W).sum();
			auto g_hess = G2[i].cwiseProduct(dens).sum();
			H_nr(row, col) += h_hess + g_hess - w_hess; 
		}
	}
	
	// Nuclear repulsion hessian
	Matrix HN = Matrix::Zero(ncoords, ncoords);
	for (auto a1 = 1; a1 != atoms.size(); ++a1) {
		const auto& atom1 = atoms[a1];
		for (auto a2 = 0; a2 < a1; ++a2) {
			const auto& atom2 = atoms[a2];

			auto x12 = atom1.getX() - atom2.getX();
			auto y12 = atom1.getY() - atom2.getY();
			auto z12 = atom1.getZ() - atom2.getZ();
			auto x12_2 = x12 * x12;
			auto y12_2 = y12 * y12;
			auto z12_2 = z12 * z12;
			auto r12_2 = x12 * x12 + y12 * y12 + z12 * z12;
			auto r12 = sqrt(r12_2);
			auto r12_5 = r12 * r12_2 * r12_2;

			auto z1z2_over_r12_5 = atom1.getEffectiveCharge() * atom2.getEffectiveCharge() / r12_5;

			HN(3*a1 + 0, 3*a1 + 0) += z1z2_over_r12_5 * (3*x12_2 - r12_2);
			HN(3*a1 + 1, 3*a1 + 1) += z1z2_over_r12_5 * (3*y12_2 - r12_2);
			HN(3*a1 + 2, 3*a1 + 2) += z1z2_over_r12_5 * (3*z12_2 - r12_2);
			HN(3*a1 + 0, 3*a1 + 1) += z1z2_over_r12_5 * (3*x12*y12);
			HN(3*a1 + 0, 3*a1 + 2) += z1z2_over_r12_5 * (3*x12*z12);
			HN(3*a1 + 1, 3*a1 + 2) += z1z2_over_r12_5 * (3*y12*z12);

			HN(3*a2 + 0, 3*a2 + 0) += z1z2_over_r12_5 * (3*x12_2 - r12_2);
			HN(3*a2 + 1, 3*a2 + 1) += z1z2_over_r12_5 * (3*y12_2 - r12_2);
			HN(3*a2 + 2, 3*a2 + 2) += z1z2_over_r12_5 * (3*z12_2 - r12_2);
			HN(3*a2 + 0, 3*a2 + 1) += z1z2_over_r12_5 * (3*x12*y12);
			HN(3*a2 + 0, 3*a2 + 2) += z1z2_over_r12_5 * (3*x12*z12);
			HN(3*a2 + 1, 3*a2 + 2) += z1z2_over_r12_5 * (3*y12*z12);

			HN(3*a2 + 0, 3*a1 + 0) -= z1z2_over_r12_5 * (3*x12_2 - r12_2);
			HN(3*a2 + 1, 3*a1 + 1) -= z1z2_over_r12_5 * (3*y12_2 - r12_2);
			HN(3*a2 + 2, 3*a1 + 2) -= z1z2_over_r12_5 * (3*z12_2 - r12_2);
			HN(3*a2 + 1, 3*a1 + 0) -= z1z2_over_r12_5 * (3*y12*x12);
			HN(3*a2 + 2, 3*a1 + 0) -= z1z2_over_r12_5 * (3*z12*x12);
			HN(3*a2 + 2, 3*a1 + 1) -= z1z2_over_r12_5 * (3*z12*y12);
			HN(3*a2 + 0, 3*a1 + 1) -= z1z2_over_r12_5 * (3*x12*y12);
			HN(3*a2 + 0, 3*a1 + 2) -= z1z2_over_r12_5 * (3*x12*z12);
			HN(3*a2 + 1, 3*a1 + 2) -= z1z2_over_r12_5 * (3*y12*z12);
		}
	}
	
	H_nr += HN; 
	
	// Transform integrals and derivative matrices to MO basis 
	
	MP2 mp2(*this); 
	mp2.transformIntegrals(false); 
	S8EvenTensor4& moInts = mp2.getMOInts();
	
	std::vector<EMatrix> S1_mo;
	std::vector<EMatrix> H1_mo;
	std::vector<EMatrix> G1_mo; 
	for (int i = 0; i < ncoords; i++) {
		auto stilde = CP.transpose() * S1[i] * CP; 
		auto htilde = CP.transpose() * (T1[i] + V1[i]) * CP; 
		auto gtilde = CP.transpose() * G1[i] * CP;
		S1_mo.push_back(stilde);
		H1_mo.push_back(htilde); 
		G1_mo.push_back(gtilde); 
	}
	
	// Calculate Q(1) and V for CPHF 
	int nvirt = CP.rows() - nocc; 
	int nix =  nvirt * nocc; 
	
	Matrix Vaibj = Matrix::Zero(nix, nix); 
	int row = 0; int col = 0;
	for (int a = 0; a < nvirt; a++) {
		for (int i = 0; i < nocc; i++) {
			col = 0; 
			
			auto dia = eps[i] - eps[a+nocc]; 
			
			for (int b = 0; b < nvirt; b++) {
				for (int j = 0; j < nocc; j++) {
					Vaibj(row, col) = -moInts(a+nocc, b+nocc, i, j) - moInts(a+nocc, j, i, b+nocc); 
					Vaibj(row, col) += 4.0 * moInts(a+nocc, i, b+nocc, j); 
					Vaibj(row, col) /= dia; 
					col++; 
				}
			}
			row++; 
		}
	} 
	
	std::vector<Vector> Q1; 
	for (int nc = 0; nc < ncoords; nc++) {
		Vector qai = Vector::Zero(nix); 
		
		row = 0;
		for (int a = 0; a < nvirt; a++) {
			for (int i = 0; i < nocc; i++) {
				qai[row] = H1_mo[nc](a+nocc, i); 
				qai[row] -= S1_mo[nc](a+nocc, i) * eps[i]; 
				
				for (int k = 0; k < nocc; k++)
					for (int l = 0; l < nocc; l++)
						qai[row] -= S1_mo[nc](k, l) * (2.0 * moInts(a+nocc, i, l, k) - moInts(a+nocc, k, l, i)); 
				
				qai[row] += G1_mo[nc](a+nocc, i); 
				qai[row] /= eps[i] - eps[a+nocc]; 
				
				row++;
			}				
		}
		
		Q1.push_back(qai); 
	}
	
	// Call CPHF solver
	std::vector<Vector> uai = cphf_group_solver(Q1, Vaibj, 1e-4, 30); 
	
	// Calculate perturbed Fock, density, and weighted density matrices, and transform to AO basis
	std::vector<Matrix> P1, W1;
	int nmos = nocc+nvirt; 
	for (int nc = 0; nc < ncoords; nc++) {
		Matrix ftilde = H1_mo[nc] + G1_mo[nc]; 
		
		for (int q = 0; q < nmos; q++) {
			for (int p = 0; p < nmos; p++) {
				
				for (int i = 0; i < nocc; i++) {
					
					for (int j = 0; j < nocc; j++) 
						ftilde(q, p) -= S1_mo[nc](j, i) * ( 2.0 * moInts(q, p, j, i) - moInts(q, i, j, p)); 
					
					for (int a = nocc; a < nmos; a++)
						ftilde(q, p) += uai[nc][(a-nocc)*nocc+i] * ( 4.0 * moInts(q, p, i, a) - moInts(q, i, a, p) - moInts(q, a, p, i) );
					
				}
			}
		}
		
		Matrix ptilde = Matrix::Zero(nmos, nmos);
		Matrix wtilde = Matrix::Zero(nmos, nmos); 
		
		for (int i = 0; i < nocc; i++) {
			for (int j = 0; j < nocc; j++) {
				ptilde(i, j) = -S1_mo[nc](j, i); 
				wtilde(i, j) = ftilde(j, i) - (eps[i] + eps[j])*S1_mo[nc](j, i); 
			}
			
			for (int a = nocc; a < nmos; a++) {
				ptilde(i, a) = 2.0*uai[nc]((a-nocc)*nocc+i); 
				wtilde(i, a) = eps[i] * ptilde(i, a); 
			}
		} 
		
		ptilde = CP * ptilde * CP.transpose(); 
		wtilde = CP * wtilde * CP.transpose();
		P1.push_back(ptilde);
		W1.push_back(wtilde); 
	} 
	
	// Determine response contributions to Hessian
	Matrix H_r = Matrix::Zero(ncoords, ncoords); 
	for (auto i = 0; i < ncoords; i++) {
		for (auto j = i; j < ncoords; j++) {
					
			for (int mu = 0; mu < nmos; mu++) {
				for (int nu = 0; nu < nmos; nu++) {
					
					H_r(i, j) += P1[j](mu, nu) * (T1[i](mu, nu) + V1[i](mu, nu)); 
					H_r(i, j) += P1[j](mu, nu) * G1[i](mu, nu); 
					H_r(i, j) -= W1[j](mu, nu) * S1[i](mu, nu); 
 					
				}
			}
			
		}
	}
	
	hessian = H_nr + 2.0 * H_r; 
	for (auto i = 0; i < ncoords; i++)
		for (auto j = 0; j < i; j++)
			hessian(i, j) = hessian(j, i);
	
}


void Fock::compute_hessian_numeric(const std::vector<Atom> &atoms, int nocc, Command& cmd) {
	
	const auto N = atoms.size();
	const auto ncoords = 3 * N;
	const double h = 0.01; 
	
	hessian = Matrix::Zero(ncoords, ncoords); 
	
	for (int i = 0; i < N; i++) {
		
		for (int xyz = 0; xyz < 3; xyz++) {
			double deltax[3] = {0.0, 0.0, 0.0}; 
			
			deltax[xyz] = h; 
			molecule->getAtom(i).translate(deltax[0], deltax[1], deltax[2]); 
			molecule->updateBasisPositions();
			
			IntegralEngine ints1(molecule, false);
			Fock f1(cmd, ints1, molecule);
			SCF hf1(cmd, molecule, f1); 
			hf1.rhf(false);   
			
			std::vector<Atom> atomlist; 
			for (int i = 0; i < molecule->getNAtoms(); i++) atomlist.push_back(molecule->getAtom(i));
			f1.compute_forces(atomlist, nocc); 
			
			deltax[xyz] = -2.0 * h; 
			molecule->getAtom(i).translate(deltax[0], deltax[1], deltax[2]); 
			molecule->updateBasisPositions();
			
			IntegralEngine ints2(molecule, false);
			Fock f2(cmd, ints2, molecule);
			SCF hf2(cmd, molecule, f2); 
			hf2.rhf(false);   
			
			atomlist.clear();
			for (int i = 0; i < molecule->getNAtoms(); i++) atomlist.push_back(molecule->getAtom(i));
			f2.compute_forces(atomlist, nocc); 

			Matrix hi = (f1.getForces() - f2.getForces()) / (2.0 * h); 
			
			for (int j = 0; j < N; j++) 
				for (int abc = 0; abc < 3; abc++) 
					hessian(3*i + xyz, 3*j + abc) = fabs(hi(j, abc)) > 1e-12 ? hi(j, abc) : 0.0; 
			
			deltax[xyz] = h; 
			molecule->getAtom(i).translate(deltax[0], deltax[1], deltax[2]);
			molecule->updateBasisPositions();
		}
	}
	
	for (int i = 0; i < ncoords; i++) 
		for (int j = 0; j < i; j++)
			hessian(i, j) = hessian(j, i); 
}

template <unsigned deriv_order>
std::vector<EMatrix> Fock::compute_2body_fock_deriv(const std::vector<Atom> &atoms, const EMatrix& D)
{
	using libint2::Shell;	
	using libint2::Engine;	
	using libint2::Operator;		
	
	std::vector<Shell>& shells = molecule->getBasis().getIntShells();
	const Matrix& Schwarz = integrals.getPrescreen();
	const auto n = integrals.nbasis(shells);
	const auto nshells = shells.size();
	const auto nderiv = libint2::num_geometrical_derivatives(atoms.size(), deriv_order);
	const auto nderiv_shellset = libint2::num_geometrical_derivatives(atoms.size(), deriv_order);
	const auto ncoords_times_2 = (atoms.size() * 3) * 2;
	
	std::vector<EMatrix> G(nderiv, EMatrix::Zero(n, n));
	
	const auto do_schwarz_screen = Schwarz.cols() != 0 && Schwarz.rows() != 0;
	EMatrix D_shblk_norm = integrals.compute_shellblock_norm(shells, dens);
	auto fock_precision = precision;
	auto maxnprim =  integrals.max_nprim(shells);
	auto maxnprim4 = maxnprim * maxnprim * maxnprim * maxnprim;
	auto engine_precision = std::min(fock_precision / D_shblk_norm.maxCoeff(), std::numeric_limits<double>::epsilon()) / maxnprim4;
	
	Engine engine(Operator::coulomb, maxnprim, integrals.max_l(shells), deriv_order);
	engine.set_precision(engine_precision);
	
	auto shell2bf = integrals.map_shell_to_basis_function(shells); 
	auto shell2atom = integrals.map_shell_to_atom(atoms, shells);
	
	const auto& buf = engine.results();
	
	size_t shell_atoms[4];
	auto obs_shellpair_list = integrals.compute_shellpair_list(shells);
	
	// loop over permutationally-unique set of shells
	for (auto s1 = 0l, s1234 = 0l; s1 != nshells; ++s1) {
		auto bf1_first = shell2bf[s1];  // first basis function in this shell
		auto n1 = shells[s1].size();       // number of basis functions in this shell
		shell_atoms[0] = shell2atom[s1];

		for (const auto& s2 : obs_shellpair_list[s1]) {
			auto bf2_first = shell2bf[s2];
			auto n2 = shells[s2].size();
			shell_atoms[1] = shell2atom[s2];

			const auto Dnorm12 = do_schwarz_screen ? D_shblk_norm(s1, s2) : 0.0;

			for (auto s3 = 0; s3 <= s1; ++s3) {
				auto bf3_first = shell2bf[s3];
				auto n3 = shells[s3].size();
				shell_atoms[2] = shell2atom[s3];

				const auto Dnorm123 = do_schwarz_screen ? std::max(D_shblk_norm(s1, s3),
				std::max(D_shblk_norm(s2, s3), Dnorm12)) : 0.0;

				const auto s4_max = (s1 == s3) ? s2 : s3;
				for (const auto& s4 : obs_shellpair_list[s3]) {
					if (s4 > s4_max)
						break;

					const auto Dnorm1234 = do_schwarz_screen ? std::max(D_shblk_norm(s1, s4),
					std::max(D_shblk_norm(s2, s4), std::max(D_shblk_norm(s3, s4), Dnorm123))) : 0.0;

					if (do_schwarz_screen && Dnorm1234 * Schwarz(s1, s2) * Schwarz(s3, s4) < fock_precision)
						continue;

					auto bf4_first = shell2bf[s4];
					auto n4 = shells[s4].size();
					shell_atoms[3] = shell2atom[s4];

					const auto n1234 = n1 * n2 * n3 * n4;

					// compute the permutational degeneracy (i.e. # of equivalents) of
					// the given shell set
					auto s12_deg = (s1 == s2) ? 1.0 : 2.0;
					auto s34_deg = (s3 == s4) ? 1.0 : 2.0;
					auto s12_34_deg = (s1 == s3) ? (s2 == s4 ? 1.0 : 2.0) : 2.0;
					auto s1234_deg = s12_deg * s34_deg * s12_34_deg;

					// computes contribution from shell set \c idx to the operator matrix with
					// index \c op
					auto add_shellset_to_dest = [&](std::size_t op, std::size_t idx, int coord1, int coord2, double scale = 1.0) {
						auto& g = G[op];
						auto shset = buf[idx];
						const auto weight = scale * s1234_deg;

						for (auto f1 = 0, f1234 = 0; f1 != n1; ++f1) {
							const auto bf1 = f1 + bf1_first;
							for (auto f2 = 0; f2 != n2; ++f2) {
								const auto bf2 = f2 + bf2_first;
								for (auto f3 = 0; f3 != n3; ++f3) {
									const auto bf3 = f3 + bf3_first;
									for (auto f4 = 0; f4 != n4; ++f4, ++f1234) {
										const auto bf4 = f4 + bf4_first;

										const auto value = shset[f1234];
										const auto wvalue = value * weight;

										g(bf1, bf2) += D(bf3, bf4) * wvalue;
										g(bf3, bf4) += D(bf1, bf2) * wvalue;
										g(bf1, bf3) -= 0.25 * D(bf2, bf4) * wvalue;
										g(bf2, bf4) -= 0.25 * D(bf1, bf3) * wvalue;
										g(bf1, bf4) -= 0.25 * D(bf2, bf3) * wvalue;
										g(bf2, bf3) -= 0.25 * D(bf1, bf4) * wvalue;
									}
								}
							}
						}
					};

					engine.compute2<Operator::coulomb, libint2::BraKet::xx_xx, deriv_order>(shells[s1], shells[s2], shells[s3], shells[s4]);
					if (buf[0] == nullptr)
						continue; // if all integrals screened out, skip to next quartet

					switch (deriv_order) {
						case 0: {
							int coord1 = 0, coord2 = 0;
							add_shellset_to_dest(0, 0, coord1, coord2);
						} break;

						case 1: {
							for (auto d = 0; d != 12; ++d) {
								const int a = d / 3;
								const int xyz = d % 3;

								auto coord = shell_atoms[a] * 3 + xyz;
								auto& g = G[coord];

								int coord1 = 0, coord2 = 0;

								add_shellset_to_dest(coord, d, coord1, coord2);

							}  // d \in [0,12)
						} break;

						case 2: {
							// computes upper triangle index
							// n2 = matrix size times 2
							// i,j = (unordered) indices
#define upper_triangle_index(n2, i, j) (i < j ? i : j) * (n2 - (i < j ? i : j) - 1) / 2 + (i > j ? i : j)
							// look over shellsets in the order in which they appear
							std::size_t shellset_idx = 0;
							for (auto c1 = 0; c1 != 4; ++c1) {
								auto a1 = shell_atoms[c1];
								auto coord1 = 3 * a1;
								for (auto xyz1 = 0; xyz1 != 3; ++xyz1, ++coord1) {
									for (auto c2 = c1; c2 != 4; ++c2) {
										auto a2 = shell_atoms[c2];
										auto xyz2_start = (c1 == c2) ? xyz1 : 0;
										auto coord2 = 3 * a2 + xyz2_start;
										for (auto xyz2 = xyz2_start; xyz2 != 3; ++xyz2, ++coord2) {
											double scale = (coord1 == coord2 && c1 != c2) ? 2.0 : 1.0;

											const auto coord12 = upper_triangle_index(ncoords_times_2, coord1, coord2);
											add_shellset_to_dest(coord12, shellset_idx, coord1, coord2, scale);
											++shellset_idx;
										}
									}
								}
							}
						} break;
#undef upper_triangle_index

						default:
						assert(deriv_order <= 2 &&
							"support for 3rd and higher derivatives of the Fock "
								"matrix not yet implemented");
					}
				}
			}
		}
	}

	std::vector<EMatrix> GG(nderiv);
	for (auto d = 0; d != nderiv; ++d) {
		GG[d] = 0.5 * (G[d] + G[d].transpose());
	}

	return GG;
}

Vector Fock::compute_xgrad(double fx, Matrix& xhessian, std::vector<int>& activex, Command& cmd) {
	
	bool withmp2 = cmd.get_option<bool>("mp2"); 
	
	Basis& b = molecule->getBasis(); 
	const auto nexp = activex.size();
	const double h = 0.01; 
	
	Vector xgrad = Vector::Zero(nexp);
	xhessian = Matrix::Zero(nexp, nexp); 
	
	double currex, currey, fxp, fyp; 
	int xctr = 0;  
	for (int x : activex) {
		currex = b.getExp(x); 
	
		b.setExp(x, currex+h); 
			
		IntegralEngine ints1(molecule, false);
		Fock f1(cmd, ints1, molecule);
		SCF hf1(cmd, molecule, f1); 
		hf1.rhf(false);   
		fxp = hf1.getEnergy(); 
		
		if (withmp2) {
			MP2 mp2(f1);
			mp2.tensormp2(false);
			fxp += mp2.getEnergy(); 
		}
		
		xgrad[xctr] = fxp; 	
		
		int yctr = xctr;
		for (int y : activex) {
			if (y >= x) {
				currey = b.getExp(y); 
			
				b.setExp(y, currey+h); 
			
				IntegralEngine ints2(molecule, false);
				Fock f2(cmd, ints2, molecule);
				SCF hf2(cmd, molecule, f2); 
				hf2.rhf(false); 
				fyp = hf2.getEnergy(); 
				
				if (withmp2) {
					MP2 mp2(f2);
					mp2.tensormp2(false);
					fyp += mp2.getEnergy(); 
				}
			
				xhessian(xctr, yctr++) = fyp; 
			
				b.setExp(y, currey); 
			}
		}
		
		b.setExp(x, currex); 
		xctr++; 
	}

	double h2 = h*h; 
	for (int x = 0; x < nexp; x++) {
		for (int y = x; y < nexp; y++) {
			xhessian(x, y) += fx - xgrad[x] - xgrad[y]; 
			xhessian(x, y) /= h2; 
			xhessian(y, x) = xhessian(x, y); 
		}
		xgrad[x] -= fx;
		xgrad[x] /= h; 
	}

	return xgrad; 
}
