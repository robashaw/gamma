#include "fock.hpp"
#include "error.hpp"
#include "atom.hpp"
#include "tensor4.hpp"
#include <iostream>
#include <iomanip>
#include <cmath>
#include <libint2.hpp>
#include "logger.hpp"
#include "ProgramController.hpp"
#include "eigen_util.hpp"
#include <string>
#include <future>
#include <thread>
#include "ioutil.hpp"


typedef std::vector<libint2::Shell> BasisSet; 
 
void Fock::makeJK(Matrix& P, double multiplier)
{
	
	if (density_fitted){
		formJKdf(P, multiplier); 
	} else if (twoints) {
		formJK(P, multiplier); 
	} else if (direct) {
		formJKdirect(integrals.getPrescreen(), P, multiplier);
	} else {
		try {
			formJKfile();
		} catch (Error e) {
			molecule->control->log.error(e);
		}
	}
}

// Form the 2J-K matrix, given that twoints is stored in memory
void Fock::formJK(Matrix& P, double multiplier)
{
	
	jints = Matrix::Zero(nbfs, nbfs);
	kints = Matrix::Zero(nbfs, nbfs);
	
	for (int u = 0; u < nbfs; ++u){
		for (int v = 0; v < nbfs ; v++){
			for (int s = 0; s < nbfs; s++){
				for (int l = 0; l < nbfs; l++){
					jints(u, v) += multiplier * P(s, l)*integrals.getERI(u, v, s, l);
					kints(u, v) += multiplier * P(s, l)*integrals.getERI(u, l, s, v);
				}
			}
		}
	}
	
	jkints = jints - 0.5*kints;
}

void Fock::formJKlocal(Matrix& Cocc, const Matrix& sigmainv, Matrix& Pt, std::vector<FragmentInfo>& finfo, bool from_file) {
	
	if (from_file)
		compute_2body_fock_df_local_file(Cocc, sigmainv, Pt, finfo, jints, kints, lmo_domains, ao_domains, fit_domains);
	else
		compute_2body_fock_df_local(Cocc, sigmainv, Pt, finfo, jints, kints, lmo_domains, ao_domains, fit_domains); 
	jkints = jints - kints; 

}

// Form JK using integral direct methods
void Fock::formJKdirect(const Matrix& Schwarz, Matrix& D, double multiplier)
{
	using libint2::Shell;
	using libint2::Engine;
	using libint2::Operator;
	
	std::vector<Shell>& shells = molecule->getBasis().getIntShells();
	const auto n = integrals.nbasis(shells);
	
	const auto do_schwarz_screen = Schwarz.cols() != 0 && Schwarz.rows() != 0;
	EMatrix D_shblk_norm = integrals.compute_shellblock_norm(shells, D);
	auto fock_precision = precision;
	auto maxnprim =  integrals.max_nprim(shells);
	auto maxnprim4 = maxnprim * maxnprim * maxnprim * maxnprim;
	auto engine_precision = std::min(fock_precision / D_shblk_norm.maxCoeff(), std::numeric_limits<double>::epsilon()) / maxnprim4;
	
	int nthreads = molecule->control->get_option<int>("nthreads");

	std::vector<Matrix> J(nthreads, Matrix::Zero(n, n));
	std::vector<Matrix> K(nthreads, Matrix::Zero(n, n));

	std::vector<Engine> engines(nthreads);
	engines[0] =
		Engine(Operator::coulomb, integrals.max_nprim(shells), integrals.max_l(shells), 0);
	engines[0].set_precision(engine_precision);  // shellset-dependent precision control
	// will likely break positive
	// definiteness
	// stick with this simple recipe
	for (size_t i = 1; i != nthreads; ++i) {
		engines[i] = engines[0];
	}
	
	auto shell2bf = integrals.map_shell_to_basis_function(shells); 
	
	
	auto compute = [&](int thread_id) {
		auto& engine = engines[thread_id]; 
		const auto& buf = engine.results();
		auto& j = J[thread_id];
		auto& k = K[thread_id]; 
	
		for(auto s1=0l, s1234=0l; s1!=shells.size(); ++s1) {

			auto bf1_first = shell2bf[s1]; // first basis function in this shell
			auto n1 = shells[s1].size();   // number of basis functions in this shell

			for(auto s2=0; s2<=s1; ++s2) {

				auto bf2_first = shell2bf[s2];
				auto n2 = shells[s2].size();
		
				const auto Dnorm12 = do_schwarz_screen ? D_shblk_norm(s1, s2) : 0.0;

				for(auto s3=0; s3<=s1; ++s3) {

					auto bf3_first = shell2bf[s3];
					auto n3 = shells[s3].size();
		  
					const auto Dnorm123 = do_schwarz_screen ? std::max(D_shblk_norm(s1, s3),
					std::max(D_shblk_norm(s2, s3), Dnorm12)) : 0.0;

					const auto s4_max = (s1 == s3) ? s2 : s3;
					for(auto s4=0; s4<=s4_max; ++s4, ++s1234) {
					
						if (s1234 % nthreads != thread_id) continue;
			
						const auto Dnorm1234 = do_schwarz_screen ? std::max(D_shblk_norm(s1, s4),
						std::max(D_shblk_norm(s2, s4), std::max(D_shblk_norm(s3, s4), Dnorm123))) : 0.0;
						
						if(do_schwarz_screen && Dnorm1234 * Schwarz(s1, s2) * Schwarz(s3, s4) < fock_precision) 
							continue;
					
						auto bf4_first = shell2bf[s4];
						auto n4 = shells[s4].size();

						// compute the permutational degeneracy (i.e. # of equivalents) of the given shell set
						auto s12_deg = (s1 == s2) ? 1.0 : 2.0;
						auto s34_deg = (s3 == s4) ? 1.0 : 2.0;
						auto s12_34_deg = (s1 == s3) ? (s2 == s4 ? 1.0 : 2.0) : 2.0;
						auto s1234_deg = s12_deg * s34_deg * s12_34_deg;

						engine.compute(shells[s1], shells[s2], shells[s3], shells[s4]);
						const auto* buf_1234 = buf[0];
						if (buf_1234 == nullptr)
						  	continue; // if all integrals screened out, skip to next quartet
						
						for(auto f1=0, f1234=0; f1!=n1; ++f1) {
							const auto bf1 = f1 + bf1_first;
							for(auto f2=0; f2!=n2; ++f2) {
								const auto bf2 = f2 + bf2_first;
								for(auto f3=0; f3!=n3; ++f3) {
									const auto bf3 = f3 + bf3_first;
									for(auto f4=0; f4!=n4; ++f4, ++f1234) {
										const auto bf4 = f4 + bf4_first;

										const auto value = buf_1234[f1234];

										const auto value_scal_by_deg = value * s1234_deg;

										j(bf1,bf2) += multiplier * D(bf3,bf4) * value_scal_by_deg;
										j(bf3,bf4) += multiplier * D(bf1,bf2) * value_scal_by_deg;
										k(bf1,bf3) += 0.5 * multiplier * D(bf2,bf4) * value_scal_by_deg;
										k(bf2,bf4) += 0.5 * multiplier * D(bf1,bf3) * value_scal_by_deg;
										k(bf1,bf4) += 0.5 * multiplier * D(bf2,bf3) * value_scal_by_deg;
										k(bf2,bf3) += 0.5 * multiplier * D(bf1,bf4) * value_scal_by_deg;
										
									}
								}
							}
						}

					}
				}
			}
		}
	};
	parallel_do(compute, nthreads); 
	
	// symmetrize the result
    for (size_t i = 1; i != nthreads; ++i) {
      J[0] += J[i];
	  K[0] += K[i]; 
    }
	J[0] = 0.25 * (J[0] + J[0].transpose()).eval();
	K[0] = 0.25 * (K[0] + K[0].transpose()).eval();
	jkints = J[0] - 0.5*K[0];
}

// Form the JK matrix from two electron integrals stored on file
void Fock::formJKfile()
{
}
	
void Fock::formJKdf(Matrix& Cocc, double multiplier) {
	compute_2body_fock_df(Cocc, jints, kints);
	jkints = multiplier * (2.0 * jints - kints);  
}	

void Fock::makeFock(Matrix& P, double multiplier)
{
	if (not density_fitted) {
		focka = fockinc; 
		if (not started_incremental &&
		rms_error < incremental_threshold) {
			started_incremental = true;
			reset_incremental = false;
			next_reset = rms_error / 1e1;
			last_reset = iter - 1; 
		}
		if (reset_incremental || not started_incremental) {
			focka = hcore;
			dens_diff = dens;
		}
		if (reset_incremental && started_incremental) {
			reset_incremental = false;
			last_reset = iter; 
			next_reset = rms_error / 1e1;
		}
	} else 
		focka = hcore; 
	
	makeJK(P, multiplier);  
	focka += jkints; 
	if (not density_fitted) fockinc = focka; 
	addDiis(); 
}

void Fock::addDiis() {
	if (diis) { // Archive for averaging
		if (iter >= MAX) {
			focks.erase(focks.begin());
		}
		focks.push_back(focka);
		iter++;
	}		
}

void Fock::makeFock() {
	if (not density_fitted)
		makeFock(dens_diff); 
	else {
		Matrix cocc = CP.block(0, 0, nbfs, nocc); 
		makeFock(cocc); 
	}
}

Vector FockFragment::buildFock(Matrix& qfq, Matrix& qfp, Matrix& pfp, bool alpha) 
{
	int nbfs = end - start; 
		
	focka = qfp.block(start, start, nbfs, nbfs) * Sxx; 
	Matrix e = focka - focka.transpose(); 
	Vector err(Eigen::Map<Vector>(e.data(), e.cols()*e.rows()));
		
	focka = focka + focka.transpose() + qfq.block(start, start, nbfs, nbfs) + Sxx * pfp.block(start, start, nbfs, nbfs) * Sxx; 
		
	if (diis) { // Archive for averaging
		if (iter >= MAX) {
			focks.erase(focks.begin());
		}
		focks.push_back(focka);
	}		
	
	return err; 
}
// Form the 2J-K matrix, given that twoints is stored in memory
void UnrestrictedFock::formJK(Matrix& Pa, Matrix& Pb, double multiplier)
{
	kints_alpha = Matrix::Zero(nbfs, nbfs);
	kints_beta = Matrix::Zero(nbfs, nbfs); 
	jints_alpha = Matrix::Zero(nbfs, nbfs);
	jints_beta = Matrix::Zero(nbfs, nbfs);
	for (int u = 0; u < nbfs; u++){
		for (int v = 0; v < nbfs ; v++){
			for (int s = 0; s < nbfs; s++){
				for (int l = 0; l < nbfs; l++){
					auto ival = multiplier * integrals.getERI(u, v, s, l); 
					jints_alpha(u, v) += Pa(s, l) * ival;
					jints_beta(u, v) += Pb(s, l) * ival;
					ival =  multiplier * integrals.getERI(u, l, s, v);
					kints_alpha(u, v) += Pa(s, l) * ival;
					kints_beta(u, v) += Pb(s, l) * ival;
				}
			}
		}
	}
}

// Form JK using integral direct methods
void UnrestrictedFock::formJKdirect(const Matrix& Schwarz, Matrix& Da, Matrix& Db, double multiplier)
{
	using libint2::Shell;
	using libint2::Engine;
	using libint2::Operator;
		
	std::vector<Shell>& shells = molecule->getBasis().getIntShells();
	const auto n = integrals.nbasis(shells);
	kints_alpha = Matrix::Zero(n, n);
	kints_beta = Matrix::Zero(n, n); 
	jints_alpha = Matrix::Zero(n, n);
	jints_beta = Matrix::Zero(n, n);
		
	const auto do_schwarz_screen = Schwarz.cols() != 0 && Schwarz.rows() != 0;
	Matrix D_tot = Da + Db; 
	EMatrix D_shblk_norm = integrals.compute_shellblock_norm(shells, D_tot);
	auto fock_precision = precision;
	auto maxnprim =  integrals.max_nprim(shells);
	auto maxnprim4 = maxnprim * maxnprim * maxnprim * maxnprim;
	auto engine_precision = std::min(fock_precision / D_shblk_norm.maxCoeff(), std::numeric_limits<double>::epsilon()) / maxnprim4;
		
	Engine engine(Operator::coulomb, integrals.max_nprim(shells), integrals.max_l(shells), 0);
	engine.set_precision(engine_precision);
	auto shell2bf = integrals.map_shell_to_basis_function(shells); 
		
	const auto& buf = engine.results();
		
	for(auto s1=0; s1!=shells.size(); ++s1) {
		auto bf1_first = shell2bf[s1]; // first basis function in this shell
		auto n1 = shells[s1].size();   // number of basis functions in this shell
		for(auto s2=0; s2<=s1; ++s2) {
			auto bf2_first = shell2bf[s2];
			auto n2 = shells[s2].size();
			
			const auto Dnorm12 = do_schwarz_screen ? D_shblk_norm(s1, s2) : 0.0;
			for(auto s3=0; s3<=s1; ++s3) {
				auto bf3_first = shell2bf[s3];
				auto n3 = shells[s3].size();
			  
				const auto Dnorm123 = do_schwarz_screen ? std::max(D_shblk_norm(s1, s3),
				std::max(D_shblk_norm(s2, s3), Dnorm12)) : 0.0;
				const auto s4_max = (s1 == s3) ? s2 : s3;
				for(auto s4=0; s4<=s4_max; ++s4) {
					const auto Dnorm1234 = do_schwarz_screen ? std::max(D_shblk_norm(s1, s4),
					std::max(D_shblk_norm(s2, s4), std::max(D_shblk_norm(s3, s4), Dnorm123))) : 0.0;
							
					if(do_schwarz_screen && Dnorm1234 * Schwarz(s1, s2) * Schwarz(s3, s4) < fock_precision) 
						continue;
						
					auto bf4_first = shell2bf[s4];
					auto n4 = shells[s4].size();
					// compute the permutational degeneracy (i.e. # of equivalents) of the given shell set
					auto s12_deg = (s1 == s2) ? 1.0 : 2.0;
					auto s34_deg = (s3 == s4) ? 1.0 : 2.0;
					auto s12_34_deg = (s1 == s3) ? (s2 == s4 ? 1.0 : 2.0) : 2.0;
					auto s1234_deg = s12_deg * s34_deg * s12_34_deg;
					engine.compute(shells[s1], shells[s2], shells[s3], shells[s4]);
					const auto* buf_1234 = buf[0];
					if (buf_1234 == nullptr)
						continue; // if all integrals screened out, skip to next quartet
					for(auto f1=0, f1234=0; f1!=n1; ++f1) {
						const auto bf1 = f1 + bf1_first;
						for(auto f2=0; f2!=n2; ++f2) {
							const auto bf2 = f2 + bf2_first;
							for(auto f3=0; f3!=n3; ++f3) {
								const auto bf3 = f3 + bf3_first;
								for(auto f4=0; f4!=n4; ++f4, ++f1234) {
									const auto bf4 = f4 + bf4_first;
									const auto value = buf_1234[f1234];
									const auto value_scal_by_deg = multiplier * value * s1234_deg;
									jints_alpha(bf1,bf2) += Da(bf3, bf4) * value_scal_by_deg;
									jints_alpha(bf3,bf4) += Da(bf1, bf2) * value_scal_by_deg;
									kints_alpha(bf1,bf3) += 0.5 * Da(bf2,bf4) * value_scal_by_deg;
									kints_alpha(bf2,bf4) += 0.5 * Da(bf1,bf3) * value_scal_by_deg;
									kints_alpha(bf1,bf4) += 0.5 * Da(bf2,bf3) * value_scal_by_deg;
									kints_alpha(bf2,bf3) += 0.5 * Da(bf1,bf4) * value_scal_by_deg;
										
									jints_beta(bf1,bf2) += Db(bf3, bf4) * value_scal_by_deg;
									jints_beta(bf3,bf4) += Db(bf1, bf2) * value_scal_by_deg;
									kints_beta(bf1,bf3) += 0.5 * Db(bf2,bf4) * value_scal_by_deg;
									kints_beta(bf2,bf4) += 0.5 * Db(bf1,bf3) * value_scal_by_deg;
									kints_beta(bf1,bf4) += 0.5 * Db(bf2,bf3) * value_scal_by_deg;
									kints_beta(bf2,bf3) += 0.5 * Db(bf1,bf4) * value_scal_by_deg;
								}
							}
						}
					}
				}
			}
		}
	}
	
	// symmetrize the result
	kints_alpha = 0.25 * (kints_alpha + kints_alpha.transpose()).eval();
	kints_beta = 0.25 * (kints_beta + kints_beta.transpose()).eval();
	jints_alpha = 0.25 * (jints_alpha + jints_alpha.transpose()).eval();
	jints_beta = 0.25 * (jints_beta + jints_beta.transpose()).eval();
}

void UnrestrictedFock::formJKdf(Matrix& ca_occ, Matrix& cb_occ, double multiplier) { 
	compute_2body_fock_df(ca_occ, jints_alpha, kints_alpha);
	compute_2body_fock_df(cb_occ, jints_beta, kints_beta); 
	jints_alpha *= multiplier;
	jints_beta  *= multiplier;
	kints_alpha *= multiplier;
	kints_beta  *= multiplier; 
}

void UnrestrictedFock::formJKlocal(Matrix& ca, Matrix& cb, const Matrix& sa, const Matrix& sb, Matrix& Pa, Matrix& Pb, std::vector<FragmentInfo>& finfoa, std::vector<FragmentInfo>& finfob) {
	
	if (xyK.rows() == 0) {
		BasisSet& dfbs = molecule->getBasis().getJKShells();
		int ndf = integrals.nbasis(dfbs); 
		
		Matrix V = integrals.compute_eris_2index(dfbs);
		
		EigenSolver esa(sa);
		Matrix cat = ca * esa.operatorSqrt(); 
		build_domains(cat, V, finfoa, lmo_alpha, ao_alpha, fit_alpha); 
		
		EigenSolver esb(sb);
		Matrix cbt = cb * esb.operatorSqrt(); 
		build_domains(cbt, V, finfob, lmo_beta, ao_beta, fit_beta);	
	
		Matrix Pt = Pa + Pb; 
		build_blocked_eris(finfoa, V, Pt); 
			
		auto llt = V.selfadjointView<Eigen::Lower>().llt();
		auto& L = llt.matrixL(); 
		Matrix I(ndf, ndf);
		I.setIdentity(); 
		Linv = L.solve(I).transpose().sparseView(1e-12);
		Linv2 = Linv * Linv.transpose();
		
		xyK.resize(1, 1);	
	}
	
	compute_2body_fock_df_local(ca, sa, Pa, finfoa, jints_alpha, kints_alpha, lmo_alpha, ao_alpha, fit_alpha); 
	compute_2body_fock_df_local(cb, sb, Pb, finfob, jints_beta, kints_beta, lmo_beta, ao_beta, fit_beta); 
	
	jints_alpha *= 0.5;
	jints_beta *= 0.5;
}

// Make the JK matrix, depending on how two electron integrals are stored/needed
void UnrestrictedFock::makeJK()
{
	if (twoints){
		formJK(dens_alpha, dens_beta); 
	} else if (density_fitted) {
		Matrix ca_occ = CP_alpha.block(0, 0, nbfs, nalpha); 
		Matrix cb_occ = CP_beta.block(0, 0, nbfs, nbeta); 
		formJKdf(ca_occ, cb_occ); 
	} else if (direct) {
		formJKdirect(integrals.getPrescreen(), dens_alpha, dens_beta);
	} else {
		try {
			formJKfile();
		} catch (Error e) {
			molecule->control->log.error(e);
		}
	}
}

void UnrestrictedFock::makeFock()
{
	fock_alpha_ao = hcore + jints_alpha + jints_beta;
	fock_beta_ao = fock_alpha_ao - kints_beta; 
	fock_alpha_ao -= kints_alpha; 
	if (diis) { // Archive for averaging
		if (iter >= MAX) {
			alpha_focks.erase(alpha_focks.begin());
			beta_focks.erase(beta_focks.begin());
		}
		alpha_focks.push_back(fock_alpha_ao);
		beta_focks.push_back(fock_beta_ao);
		iter++;
	}		
}

Vector UnrestrictedFockFragment::buildFock(Matrix& qfq, Matrix& qfp, Matrix& pfp, bool alpha) 
{
	int nbfs = end - start; 
	
	focka = qfp.block(start, start, nbfs, nbfs) * Sxx; 
	
	Matrix e = focka - focka.transpose(); 
	Vector err(Eigen::Map<Vector>(e.data(), e.cols()*e.rows()));
	
	focka = focka + focka.transpose() + qfq.block(start, start, nbfs, nbfs) + Sxx * pfp.block(start, start, nbfs, nbfs) * Sxx; 
	
	if (alpha) {
		fock_alpha_ao = focka;
		if (diis) { // Archive for averaging
			if (iter >= MAX)
				alpha_focks.erase(alpha_focks.begin());
			alpha_focks.push_back(fock_alpha_ao);
		}	
	} else { 
		fock_beta_ao = focka;
		if (diis) { // Archive for averaging
			if (iter >= MAX)
				beta_focks.erase(beta_focks.begin());
			beta_focks.push_back(fock_beta_ao);
		}	
	}
	return err; 
}

// an Fock builder that can accept densities expressed a separate basis
Matrix Fock::compute_2body_fock_general(const BasisSet& obs, const Matrix& D, const BasisSet& D_bs,
bool D_is_shelldiagonal, double precision)
{
		
	const auto n = integrals.nbasis(obs);
	const auto nshells = obs.size();
	const auto n_D = integrals.nbasis(D_bs);
	assert(D.cols() == D.rows() && D.cols() == n_D);

	// construct the 2-electron repulsion integrals engine
	using libint2::Engine;
	int nthreads = molecule->control->get_option<int>("nthreads");

	std::vector<Matrix> G(nthreads, Matrix::Zero(n, n));
	
	std::vector<Engine> engines(nthreads);
	engines[0] =
		Engine(libint2::Operator::coulomb,
	std::max(integrals.max_nprim(obs), integrals.max_nprim(D_bs)),
	std::max(integrals.max_l(obs), integrals.max_l(D_bs)), 0);
	engines[0].set_precision(precision);  // shellset-dependent precision control
	// will likely break positive
	// definiteness
	// stick with this simple recipe
	for (size_t i = 1; i != nthreads; ++i) {
		engines[i] = engines[0];
	}

	auto shell2bf = integrals.map_shell_to_basis_function(obs);
	auto shell2bf_D = integrals.map_shell_to_basis_function(D_bs);

	auto compute = [&](int thread_id) {
		auto& engine = engines[thread_id]; 
		const auto& buf = engine.results();

		// loop over permutationally-unique set of shells
		for (auto s1 = 0l, s1234 = 0l; s1 != nshells; ++s1) {
			auto bf1_first = shell2bf[s1];  // first basis function in this shell
			auto n1 = obs[s1].size();       // number of basis functions in this shell

			for (auto s2 = 0; s2 <= s1; ++s2) {
				auto bf2_first = shell2bf[s2];
				auto n2 = obs[s2].size();

				for (auto s3 = 0; s3 < D_bs.size(); ++s3) {
					auto bf3_first = shell2bf_D[s3];
					auto n3 = D_bs[s3].size();

					auto s4_begin = D_is_shelldiagonal ? s3 : 0;
					auto s4_fence = D_is_shelldiagonal ? s3 + 1 : D_bs.size();

					for (auto s4 = s4_begin; s4 != s4_fence; ++s4, ++s1234) {
						
						if (s1234 % nthreads != thread_id) continue; 
						
						auto bf4_first = shell2bf_D[s4];
						auto n4 = D_bs[s4].size();

						// compute the permutational degeneracy (i.e. # of equivalents) of
						// the given shell set
						auto s12_deg = (s1 == s2) ? 1.0 : 2.0;

						if (s3 >= s4) {
							auto s34_deg = (s3 == s4) ? 1.0 : 2.0;
							auto s1234_deg = s12_deg * s34_deg;
							// auto s1234_deg = s12_deg;
							engine.compute2<libint2::Operator::coulomb, libint2::BraKet::xx_xx, 0>(
								obs[s1], obs[s2], D_bs[s3], D_bs[s4]);
							const auto* buf_1234 = buf[0];
							if (buf_1234 != nullptr) {
								for (auto f1 = 0, f1234 = 0; f1 != n1; ++f1) {
									const auto bf1 = f1 + bf1_first;
									for (auto f2 = 0; f2 != n2; ++f2) {
										const auto bf2 = f2 + bf2_first;
										for (auto f3 = 0; f3 != n3; ++f3) {
											const auto bf3 = f3 + bf3_first;
											for (auto f4 = 0; f4 != n4; ++f4, ++f1234) {
												const auto bf4 = f4 + bf4_first;

												const auto value = buf_1234[f1234];
												const auto value_scal_by_deg = value * s1234_deg;
												G[thread_id](bf1, bf2) += 2.0 * D(bf3, bf4) * value_scal_by_deg;
											}
										}
									}
								}
							}
						}

						engine.compute2<libint2::Operator::coulomb, libint2::BraKet::xx_xx, 0>(
							obs[s1], D_bs[s3], obs[s2], D_bs[s4]);
						const auto* buf_1324 = buf[0];
						if (buf_1324 == nullptr)
							continue; // if all integrals screened out, skip to next quartet

						for (auto f1 = 0, f1324 = 0; f1 != n1; ++f1) {
							const auto bf1 = f1 + bf1_first;
							for (auto f3 = 0; f3 != n3; ++f3) {
								const auto bf3 = f3 + bf3_first;
								for (auto f2 = 0; f2 != n2; ++f2) {
									const auto bf2 = f2 + bf2_first;
									for (auto f4 = 0; f4 != n4; ++f4, ++f1324) {
										const auto bf4 = f4 + bf4_first;

										const auto value = buf_1324[f1324];
										const auto value_scal_by_deg = value * s12_deg;
										G[thread_id](bf1, bf2) -= D(bf3, bf4) * value_scal_by_deg;
									}
								}
							}
						}
					}
				}
			}
		}
	};
	
	parallel_do(compute, nthreads);

	for (int i = 1; i != nthreads; ++i)
		G[0] += G[i]; 
	
	// symmetrize the result and return
	return 0.5 * (G[0] + G[0].transpose());		
}

void Fock::compute_2body_fock_df(const Matrix& Cocc, Matrix& j, Matrix& k) {

	BasisSet& obs = molecule->getBasis().getIntShells(); 
	BasisSet& dfbs = molecule->getBasis().getJKShells(); 
  
	const auto n = integrals.nbasis(obs); 
	const auto ndf = integrals.nbasis(dfbs); 
	const int ncocc = Cocc.cols(); 
	
	Logger& log = molecule->control->log; 

	// using first time? compute 3-center ints and transform to inv sqrt
	// representation
	if (xyK.rows() == 0) {
		
		integrals.compute_eris_3index(obs, dfbs, xyK); 

		Matrix V = integrals.compute_eris_2index(dfbs); 
		LLT V_LLt(V);
		Matrix I = Matrix::Identity(ndf, ndf);
		auto L = V_LLt.matrixL();
		Matrix denseLinv = L.solve(I).transpose(); 
		Linv = denseLinv.sparseView(1e-12); 

		Vector row; 
		for (int x = 0; x < n; x++) {
			for (int y = 0; y<= x; y++) {
				int xy = (x*(x+1))/2+y; 
				row = xyK.row(xy); 
			
				for (int K = 0; K < ndf; K++)
					xyK(xy, K) = row.dot(denseLinv.col(K)); 
			}
		} 
	}  // if (xyK.size() == 0)

	int nthreads = molecule->control->get_option<int>("nthreads"); 
	nthreads = std::min<int>(nthreads, n); 
	std::vector<int> limits = bounds(nthreads, n); 
	std::vector<std::thread> workers; 
	
	Matrix xiK = Matrix::Zero(n*ncocc, ndf); 
	
	auto first_trans = [&](int L, int R) {
		for (int x = L; x < R; x++) {
			for (int i = 0; i < ncocc; i++) {
				for (int K = 0; K < ndf; K++) {
					for (int y = 0; y < n; y++) {
						int X = std::max(x, y);
						int Y = std::min(x, y);
						xiK(x*ncocc+i, K) += xyK((X*(X+1))/2+Y, K) * Cocc(y, i); 
					}
				}
			}
		}
	};
	
	for (int thrd = 0; thrd < nthreads - 1; ++thrd)
		workers.push_back(std::thread(first_trans, limits[thrd], limits[thrd+1]));
	
	first_trans(limits[nthreads-1], limits[nthreads]); 
	
	for (auto& t : workers)
		t.join(); 
	
	workers.clear(); 

	k = Matrix::Zero(n, n); 

	auto second_trans = [&](int L, int R) {
		for (int x = L; x < R; x++) {
			for (int y = 0; y <= x; y++) {
				for (int i = 0; i < ncocc; i++)
					for (int K = 0; K < ndf; K++)
						k(x, y) += xiK(x*ncocc+i, K) * xiK(y*ncocc+i, K); 
				k(y, x) = k(x, y); 
			}
		}	
	};
	
	for (int thrd = 0; thrd < nthreads - 1; ++thrd)
		workers.push_back(std::thread(second_trans, limits[thrd], limits[thrd+1]));
	
	second_trans(limits[nthreads-1], limits[nthreads]); 
	
	for (auto& t : workers)
		t.join(); 
	
	workers.clear(); 

	// compute Coulomb
	Vector Jtmp = Vector::Zero(ndf); 
	for (int K = 0; K < ndf; K++) 
		for (int x = 0; x < n; x++)
			for (int i = 0; i < ncocc; i++)
				Jtmp[K] += xiK(x*ncocc+i, K) * Cocc(x, i); 
	xiK.resize(0, 0); 
  
	j = Matrix::Zero(n, n);
	for (int x = 0; x < n; x++) {
		for (int y = 0; y <= x; y++) { 
			for (int K = 0; K < ndf; K++)  
				j(x, y) += xyK((x*(x+1))/2+y, K) * Jtmp[K]; 
			j(y, x) = j(x, y); 
		}
	}
}

void Fock::build_domains(Matrix& Cocc, Matrix& V, std::vector<FragmentInfo>& finfo,
std::vector<Domain>& lmod, std::vector<Domain>& aod, std::vector<Domain>& fitd) {
	int i_offset = 0; 
	for (int B = 0; B < finfo.size(); B++) {
		for (int i = i_offset; i < i_offset+finfo[B].occ; i++) {
			Domain d; 
			
			std::vector<int> notcentres; 
			int mu_offset = 0;
			for (int A = 0; A < finfo.size(); A++) {
				bool add = false; 
				if (A == B) add = true; 
				else {	
					double sum = 0.0; 
					for (int mu = mu_offset; mu < mu_offset + finfo[A].occ; mu++) 
						sum += Cocc(mu, i) * Cocc(mu, i); 
					add = sum > finfo[0].mo_thresh; 
				}
			
				if (add) { 
					d.starts.push_back(finfo[A].start); 
					d.sizes.push_back(finfo[A].nbfs); 
					d.centres.push_back(A); 
				} else 
					notcentres.push_back(A); 
			
				mu_offset += finfo[A].nbfs; 
			}
			
			Domain dao; 
			dao.starts = d.starts;
			dao.sizes = d.sizes; 
			dao.centres = d.centres; 
			mu_offset = 0; 
			for (auto A : d.centres) {
				std::vector<int> newnc; 
				for (int b = 0; b < notcentres.size(); b++){
					int C = notcentres[b]; 
					double sep = (finfo[C].com - finfo[A].com).norm(); 
					if (finfo[C].radius > sep) {
						dao.centres.push_back(C); 
						dao.starts.push_back(finfo[C].start);
						dao.sizes.push_back(finfo[C].nbfs); 
					} else 
						newnc.push_back(C); 
				} 
				notcentres = newnc; 
			}
		
			d.sumsizes();
			dao.sumsizes(); 
			lmod.push_back(d); 
			aod.push_back(dao); 
		}
		i_offset += finfo[B].occ; 
	}
	
	//Kblocks.resize(finfo.size()); 
	
	EigenSolver es2(integrals.getOverlap()); 
	Matrix SC = es2.operatorSqrt() * Cocc; 
	i_offset = 0;
	for (auto& f1 : finfo) { 
		 
		for (int i = i_offset; i < i_offset + f1.occ; i++) {
			Domain d;
			
			int mu_offset = 0; int center = 0; 
			for (auto& f2 : finfo) { 
				
				double sep = (f2.com - f1.com).norm(); 
				
				double INi = 0.0; 
				for (int mu = mu_offset; mu < mu_offset + f2.nbfs; mu++) {
					INi += SC(mu, i) * SC(mu, i);	
				}
				
				if (INi > finfo[0].fit_thresh || sep < finfo[0].r_thresh) {
					d.centres.push_back(center);
					d.sizes.push_back(f2.naux); 
					d.starts.push_back(f2.auxstart); 
					
					//Kblocks[center].push_back(i); 
				}
			
				center++; 
				mu_offset += f2.nbfs; 
			}
			
			d.sumsizes();
			fitd.push_back(d); 
		}
		
		i_offset += f1.occ; 
	}
	
	/*for (auto& klist : Kblocks) {
	std::sort(klist.begin(), klist.end());
	klist.erase(std::unique(klist.begin(), klist.end()), klist.end());
	}*/
	
	for (int i = 0; i < nocc; i++) {
		
		auto& lmo_d = lmod[i]; 
		auto& fit_d = fitd[i];
		int fitsize = fit_d.centres.size();  
		
		Matrix GAB = Matrix::Zero(fit_d.totalsize, fit_d.totalsize); 
		int A = 0; 
		for (int a = 0; a < fitsize; a++) {
			for (int Ax = fit_d.starts[a]; Ax < fit_d.starts[a] + fit_d.sizes[a]; Ax++) {
				
				int B = 0; 
				for (int b = 0; b < fitsize; b++) {
					for (int Bx = fit_d.starts[b]; Bx < fit_d.starts[b] + fit_d.sizes[b]; Bx++) {
						GAB(A, B) = V(Ax, Bx); 
						B++; 
					}
				}
				
				A++; 
			}
		}
		
		Eigen::LLT<Matrix> G_LLt(GAB);
		Matrix I = Matrix::Identity(fit_d.totalsize, fit_d.totalsize);
		auto GL = G_LLt.matrixL();
		lmo_d.G = GL.solve(I).transpose();
		
	}
}

void Fock::compute_2body_fock_df_local(Matrix& Cocc, const Matrix& sigmainv, Matrix& Pt, std::vector<FragmentInfo>& finfo, Matrix& j, Matrix& k,
std::vector<Domain>& lmod, std::vector<Domain>& aod, std::vector<Domain>& fitd) {

	BasisSet& obs = molecule->getBasis().getIntShells(); 
	BasisSet& dfbs = molecule->getBasis().getJKShells(); 
  
	const auto n = integrals.nbasis(obs); 
	const auto ndf = integrals.nbasis(dfbs); 
	
	Logger& log = molecule->control->log; 

	EigenSolver es(sigmainv); 
	Cocc = Cocc * es.operatorSqrt(); 
	
	// using first time? compute 3-center ints and transform to inv sqrt
	// representation
	if (xyK.rows() == 0) {
		 
		Matrix V = integrals.compute_eris_2index(dfbs);
		build_blocked_eris(finfo, V, Pt); 
	
		auto llt = V.selfadjointView<Eigen::Lower>().llt();
		auto& L = llt.matrixL(); 
		Matrix I(ndf, ndf);
		I.setIdentity(); 
		Linv = L.solve(I).transpose().sparseView(1e-12);
		Linv2 = Linv * Linv.transpose();
		
		build_domains(Cocc, V, finfo, lmod, aod, fitd);
		xyK.resize(1, 1); 
		
	}  // if (xyK.size() == 0) 
	
	j = Matrix::Zero(n, n); 
	k = Matrix::Zero(n, n);
	// compute exchange
	int nfrags = finfo.size(); 

	int nthreads = molecule->control->get_option<int>("nthreads"); 
	nthreads = std::min<int>(nthreads, nocc); 
	std::vector<int> limits = bounds(nthreads, nocc); 

	auto exchange = [&](int L, int R) {
		for (int i = L; i < R; i++) {
			auto& lmo_d = lmod[i]; 
			auto& ao_d = aod[i];
			auto& fit_d = fitd[i]; 
		
			int nmo = lmo_d.centres.size();
			int nao = ao_d.centres.size();
			int nfit = fit_d.centres.size(); 
		
			Matrix uA = Matrix::Zero(ao_d.totalsize, fit_d.totalsize); 
		
			int Al = 0; 
			int fitstart, aostart, mostart, fitsize, aosize, mosize, f1, f2, fk; 
			for (int fit = 0; fit < nfit; fit++) {
				fitstart = fit_d.starts[fit]; 
				fitsize = fit_d.sizes[fit]; 
				fk = fit_d.centres[fit]; 
			
				for (int A = fitstart; A < fitstart + fitsize; A++) {
				
					int ul = 0;
					for (int ao = 0; ao < nao; ao++) {
						aostart = ao_d.starts[ao]; 
						aosize = ao_d.sizes[ao]; 
						f1 = ao_d.centres[ao];
					
						for (int u = 0; u < aosize; u++) {
							for (int mo = 0; mo < nmo; mo++) {
								mostart = lmo_d.starts[mo];
								mosize = lmo_d.sizes[mo]; 
								f2 = lmo_d.centres[mo]; 
						
								if (!blocked_xyK.isZero(f1, f2, fk)) {
									Matrix& xyf = blocked_xyK(f1, f2, fk); 
								
									for (int nu = 0; nu < mosize; nu++) {
										int U = f1 < f2 ? nu * aosize : u * mosize;
										int NU = f1 < f2 ? u : nu; 
										uA(ul, Al) += xyf(U + NU, A-fitstart) * Cocc(nu+mostart, i);
									} 
								}
							}
						
							ul++; 
						
						}
					}
		
					Al++; 
				}
			}
		
			uA *= lmo_d.G; 
			uA = uA * uA.transpose(); 
		
			int ao1start, ao2start, ao1size, ao2size; 
			int u = 0; 
			for (int ao1 = 0; ao1 < nao; ao1++) {
				ao1start = ao_d.starts[ao1];
				ao1size = ao_d.sizes[ao1]; 
			
				for (int mu = ao1start; mu < ao1start + ao1size; mu++) {
				
					int v = 0; 
					for (int ao2 = 0; ao2 < nao; ao2++) {
						ao2start = ao_d.starts[ao2]; 
						ao2size = ao_d.sizes[ao2]; 
					 
						for (int nu = ao2start; nu < ao2start + ao2size; nu++) {
							k(mu, nu) += uA(u, v); 
						
							v++; 
						}
					}
					u++; 
				}
			}
	
		}
	
	};
	
	std::vector<std::thread> workers; 
	for (int thrd = 0; thrd < nthreads - 1; ++thrd)
		workers.push_back(std::thread(exchange, limits[thrd], limits[thrd+1]));
	
	exchange(limits[nthreads-1], limits[nthreads]); 
	
	for (auto& t : workers)
		t.join(); 

	workers.clear(); 
	
	// compute Coulomb
	
	nthreads = std::min<int>(nthreads, nfrags); 
	limits = bounds(nthreads, nfrags); 
	
	std::vector<Vector> Pks(nthreads);
	for (Vector& pk : Pks) pk = Vector::Zero(ndf); 
	
	auto coulomb_first = [&](int L, int R, int thread) {
		int nao1, nao2, naux, f1start, f2start, fkstart; 
		for (int f1 = L; f1 < R; ++f1) {
			nao1 = finfo[f1].nbfs; 
			f1start = finfo[f1].start; 
		
			for (int f2 = 0; f2 <= f1; f2++) {
				nao2 = finfo[f2].nbfs; 
				f2start = finfo[f2].start;
			
				for (int fk = 0; fk < nfrags; fk++) { 
					fkstart = finfo[fk].auxstart; 
					naux = finfo[fk].naux; 
				
					if (!blocked_xyK.isZero(f1, f2, fk)) {
						Matrix& xyf = blocked_xyK(f1, f2, fk); 
				
						Matrix Ptf = Pt.block(f1start, f2start, nao1, nao2); 
						Vector Pvf(nao1 * nao2); 
						for (int x = 0; x < nao1; x++) {
							for (int y = 0; y < nao2; y++)
								Pvf[x*nao2+y] = Ptf(x, y); 
						}
			
						if (f1==f2) 
							Pks[thread].segment(fkstart, naux) += 0.5 * xyf.transpose() * Pvf; 
						else 
							Pks[thread].segment(fkstart, naux) += xyf.transpose() * Pvf; 
					}
				}
			}
		}
	};
	
	for (int thrd = 0; thrd < nthreads - 1; ++thrd)
		workers.push_back(std::thread(coulomb_first, limits[thrd], limits[thrd+1], thrd)); 
	
	coulomb_first(limits[nthreads-1], limits[nthreads], nthreads-1); 
	
	for (auto& t : workers)
		t.join();
	
	Vector Pk = Vector::Zero(ndf);  
	for (auto& pk : Pks) Pk += pk; 
	
	workers.clear(); 
	
	Pk = Linv2 * Pk;
	
	auto coulomb_second = [&](int L, int R) {
		int nao1, nao2, naux, f1start, f2start, fkstart; 
		for (int f1 = L; f1 < R; ++f1) {
			nao1 = finfo[f1].nbfs; 
			f1start = finfo[f1].start; 
		
			for (int f2 = 0; f2 <= f1; f2++) {
				nao2 = finfo[f2].nbfs; 
				f2start = finfo[f2].start; 
			
				for (int fk = 0; fk < nfrags; fk++) {
					fkstart = finfo[fk].auxstart;
					naux = finfo[fk].naux; 
				
					if (!blocked_xyK.isZero(f1, f2, fk)) {
						Matrix& xyf = blocked_xyK(f1, f2, fk); 
				
						Vector Pvf = 4.0 * xyf * Pk.segment(fkstart, naux);
						if (f1 == f2) {
							Pvf *= 0.5; 
							for (int x = 0; x < nao1; x++)
								Pvf[x*nao2+x] *= 2.0; 
						}
				
						for (int x = 0; x < nao1; x++) { 
							for (int y = 0; y < nao2; y++) {
						
								j(x+f1start, y+f2start) += Pvf[x * nao2 + y]; 
								j(y+f2start, x+f1start) = j(x+f1start, y+f2start);
							
							}
						}
					}
				}
			}
		}
	};
	
	for (int thrd = 0; thrd < nthreads - 1; ++thrd)
		workers.push_back(std::thread(coulomb_second, limits[thrd], limits[thrd+1])); 
	
	coulomb_second(limits[nthreads-1], limits[nthreads]); 
	
	for (auto& t : workers)
		t.join();
}

void Fock::build_blocked_eris(std::vector<FragmentInfo>& finfo, Matrix& JPQ, Matrix& Pt) {
	
	BasisSet& obs = molecule->getBasis().getIntShells(); 
	BasisSet& dfbs = molecule->getBasis().getJKShells(); 
	Matrix &psc = integrals.getPrescreen(); 
	
	int nfrags = finfo.size(); 
	double THRESHOLD = molecule->control->get_option<double>("thrint");
	blocked_xyK.resize(nfrags);  
	
	int nobs1 = 0;  
	for (int f1 = 0; f1 < finfo.size(); f1++) {
		auto& info1 = finfo[f1];
		
		int nobs2 = 0; 
		for (int f2 = 0; f2 <= f1; f2++) {
			auto& info2 = finfo[f2]; 
			double infnorm = psc.block(nobs1, nobs2, info1.nshells, info2.nshells).lpNorm<Eigen::Infinity>();
			double ptnorm = Pt.block(info1.start, info2.start, info1.nbfs, info2.nbfs).lpNorm<Eigen::Infinity>();  
			
			int nabsk = 0; int ndf = 0;
			for (int fk = 0; fk < finfo.size(); fk++) {
				auto& infok = finfo[fk]; 
				double infjpq = std::sqrt(JPQ.block(ndf, ndf, infok.naux, infok.naux).lpNorm<Eigen::Infinity>()); 
				double rxyz = (info1.com - infok.com).norm() + (info2.com - infok.com).norm();
				rxyz *= 0.5;  
				
				double estimate = rxyz < 1.0 ? 1.0 : 2.0 * infnorm * ptnorm * infjpq / rxyz;
				//estimate = rxyz > RTHRESH ? estimate / rxyz : estimate * infjpq; 
				
				if (estimate > THRESHOLD) {
					blocked_xyK.setNonZero(f1, f2, fk); 
					
					BasisSet obs1, obs2, dfset;
					std::copy(obs.begin()+nobs1, obs.begin()+nobs1+info1.nshells, std::back_inserter(obs1)); 
					std::copy(obs.begin()+nobs2, obs.begin()+nobs2+info2.nshells, std::back_inserter(obs2));
					std::copy(dfbs.begin() + nabsk, dfbs.begin()+nabsk+infok.ndfshells, std::back_inserter(dfset)); 
								
					integrals.compute_eris_3index(obs1, obs2, dfset, blocked_xyK(f1, f2, fk)); 
				} 
				
				nabsk += infok.ndfshells; 
				ndf += infok.naux; 
			}
			
			nobs2 += info2.nshells;
		}
		
		nobs1 += info1.nshells;
	}
	
	/*int nzero = 0; int ntotal = 0; 
	for (int f1 = 0; f1 < finfo.size(); f1++) {
	for (int f2 = 0; f2 <= f1; f2++) {
	for (int fk = 0; fk < finfo.size(); fk++) {
	if (blocked_xyK.isZero(f1, f2, fk)) nzero++; 
	ntotal++; 
	}
	}
	}*/
	//std::cout << nzero << " " << ntotal << std::endl; 
	
}


void read_intfile(int K, Matrix& xyf) {
	std::string filename = "V" + std::to_string(K) + ".ints";
	eigenutil::read_binary(filename.c_str(), xyf); 
}

void Fock::compute_2body_fock_df_local_file(Matrix& Cocc, const Matrix& sigmainv, Matrix& Pt, std::vector<FragmentInfo>& finfo, Matrix& j, Matrix& k,
std::vector<Domain>& lmod, std::vector<Domain>& aod, std::vector<Domain>& fitd) {

	BasisSet& obs = molecule->getBasis().getIntShells(); 
	BasisSet& dfbs = molecule->getBasis().getJKShells(); 
  
	const auto n = integrals.nbasis(obs); 
	const auto ndf = integrals.nbasis(dfbs); 
	
	Logger& log = molecule->control->log; 

	EigenSolver es(sigmainv); 
	Cocc = Cocc * es.operatorSqrt(); 
	
	// using first time? compute 3-center ints and transform to inv sqrt
	// representation
	if (xyK.rows() == 0) {
		
		Matrix V = integrals.compute_eris_2index(dfbs);
	
		auto L = V.selfadjointView<Eigen::Lower>().llt().matrixL(); 
		Matrix I(ndf, ndf);
		I.setIdentity(); 
		Linv = L.solve(I).transpose().sparseView(1e-12);
		Linv2 = Linv * Linv.transpose();   
		
		build_domains(Cocc, V, finfo, lmod, aod, fitd); 
		
		int aux = 0; 
		int fn = 0; 
		for (auto& f : finfo) {
			BasisSet dfset;
			std::copy(dfbs.begin() + aux, dfbs.begin() + aux + f.ndfshells, std::back_inserter(dfset)); 
			
			Matrix xyf; 
			integrals.compute_eris_3index(obs, dfset, xyf);
			std::string filename = "V"+ std::to_string(fn) + ".ints"; 
			eigenutil::write_binary(filename.c_str(), xyf); 
			
			fn++;
			aux += f.ndfshells; 
		}
		
		xyK = Matrix::Zero(1, 1); 
		
	}  // if (xyK.size() == 0) 
	
	j = Matrix::Zero(n, n); 
	k = Matrix::Zero(n, n); 
	
	std::vector<Matrix> xyf(fitd[0].centres.size()); 
	for (int f = 0; f < xyf.size(); f++)
		read_intfile(fitd[0].centres[f], xyf[f]); 
	
	// compute exchange
	int nfrags = finfo.size(); 
	double start = molecule->control->log.getGlobalTime(); 
	for (int i = 0; i < nocc; i++) {
	
		auto& lmo_d = lmod[i]; 
		auto& ao_d = aod[i];
		auto& fit_d = fitd[i]; 
		
		int nmo = lmo_d.centres.size();
		int nao = ao_d.centres.size();
		int nfit = fit_d.centres.size(); 
		
		// Read in integrals
		std::vector<std::future<void>> futures; 
		std::vector<Matrix> xyf_next;
		if (i < nocc-1) { 
			int nextnfit = fitd[i+1].centres.size();
			xyf_next.resize(nextnfit); 
			for (int f = 0; f < nextnfit; f++)
				futures.push_back(std::async(read_intfile, fitd[i+1].centres[f], std::ref(xyf_next[f])));
		} 
		
		Matrix uA = Matrix::Zero(ao_d.totalsize, fit_d.totalsize); 
		
		int Al = 0; 
		int fitstart, aostart, mostart, fitsize, aosize, mosize; 
		for (int fit = 0; fit < nfit; fit++) {
			fitstart = fit_d.starts[fit]; 
			fitsize = fit_d.sizes[fit];
						
			for (int A = fitstart; A < fitstart + fitsize; A++) {
				
				int ul = 0;
				for (int ao = 0; ao < nao; ao++) {
					aostart = ao_d.starts[ao]; 
					aosize = ao_d.sizes[ao]; 
					
					for (int u = aostart; u < aostart + aosize; u++) {
						for (int mo = 0; mo < nmo; mo++) {
							mostart = lmo_d.starts[mo];
							mosize = lmo_d.sizes[mo]; 

							for (int nu = mostart; nu < mostart + mosize; nu++) {
								int U = std::max(u, nu); 
								int NU = std::min(u, nu); 
								uA(ul, Al) += xyf[fit]((U*(U+1))/2+NU, A-fitstart) * Cocc(nu, i);
							} 
							
						}
						
						ul++; 
						
					}
				}
		
				Al++; 
			}
		}
		
		
		uA *= lmo_d.G; 
		uA = uA * uA.transpose(); 
		
		int ao1start, ao2start, ao1size, ao2size; 
		int u = 0; 
		for (int ao1 = 0; ao1 < nao; ao1++) {
			ao1start = ao_d.starts[ao1];
			ao1size = ao_d.sizes[ao1]; 
			
			for (int mu = ao1start; mu < ao1start + ao1size; mu++) {
				
				int v = 0; 
				for (int ao2 = 0; ao2 < nao; ao2++) {
					ao2start = ao_d.starts[ao2]; 
					ao2size = ao_d.sizes[ao2]; 
					 
					for (int nu = ao2start; nu < ao2start + ao2size; nu++) {

						k(mu, nu) += uA(u, v); 
						
						v++; 
					}
				}
				u++; 
			}
		}
		
		if (i < nocc-1) { 
			for (auto& future : futures) 
				future.get(); 
			xyf = xyf_next; 
		}
	}
	double end = molecule->control->log.getGlobalTime(); 
	//std::cout << "Exchange: " << end - start << " seconds,";
	 
	// compute Coulomb
	
	start = molecule->control->log.getGlobalTime(); 
	Vector Pv((n*(n+1))/2);
	for (int x = 0; x < n; x++)
		for (int y = 0; y <= x; y++)
			Pv[(x*(x+1))/2 + y] = x == y ? 0.5 * Pt(x ,y) : Pt(x, y); 
	
	Vector Pk(ndf);
	int aux = 0;  
	Matrix xyfrag; 
	read_intfile(0, xyfrag); 
	for (int f = 0; f < nfrags-1; f++) {
		Matrix xyf_next;
		std::future<void> future(std::async(read_intfile, f+1, std::ref(xyf_next))); 
		Pk.segment(aux, finfo[f].naux) = xyfrag.transpose() * Pv;  
		aux += finfo[f].naux; 
		future.get();
		xyfrag = xyf_next; 
	}
	Pk.segment(aux, finfo[nfrags-1].naux) = xyfrag.transpose() * Pv; 
	
	Pk = Linv2 * Pk;  
	
	aux = 0;
	Pv = Vector::Zero((n*(n+1))/2);
	read_intfile(0, xyfrag); 
	for (int f = 0; f < nfrags-1; f++) {
		Matrix xyf_next; 
		std::future<void> future(std::async(read_intfile, f+1, std::ref(xyf_next))); 
		Pv += 4.0 * xyfrag * Pk.segment(aux, finfo[f].naux);  
		aux += finfo[f].naux; 
		future.get();
		xyfrag = xyf_next; 
	}
	Pv += 4.0 * xyfrag * Pk.segment(aux, finfo[nfrags-1].naux); 
	
	for (int x = 0; x < n; x++)
	for (int y = 0; y <=x; y++) {
		j(x, y) += Pv[(x*(x+1))/2+y];
		j(y, x) = j(x, y); 
	}
	end = molecule->control->log.getGlobalTime(); 	
	//std::cout << "Coulomb: " << end - start << " seconds" << std::endl << std::flush; 
	
	//std::cout << std::endl << std::flush; 
}
