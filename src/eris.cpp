#include "error.hpp"
#include "integrals.hpp"
#include "mathutil.hpp"
#include "basis.hpp"
#include "logger.hpp"
#include "tensor4.hpp"
#include "ioutil.hpp"
#include <libecpint/multiarr.hpp>
#include <libecpint/gshell.hpp>
#include "ProgramController.hpp"

#include <cmath>
#include <iomanip>
#include <string>
#include <thread>  

S8EvenTensor4 IntegralEngine::compute_eris(const std::vector<libint2::Shell>& shells) {
	using libint2::Shell;
	using libint2::Engine;
	using libint2::Operator;

	const auto n = nbasis(shells);
	int nthreads = molecule->control->get_option<int>("nthreads");
	S8EvenTensor4 eris(n, 0.0);

	std::vector<Engine> engines(nthreads);
	engines[0] = Engine(Operator::coulomb, max_nprim(shells), max_l(shells), 0);
	for (size_t i = 1; i != nthreads; ++i) {
		engines[i] = engines[0];
	}

	auto shell2bf = map_shell_to_basis_function(shells);

	auto compute = [&](int thread_id) {
		auto& engine = engines[thread_id]; 
		const auto& buf = engine.results();

		// loop over shell quartets
		for (auto s1=0l, s1234=0l; s1 != shells.size(); ++s1) {
    
			auto bf1_first = shell2bf[s1];
			auto n1 = shells[s1].size();

			for (auto s2=0; s2 <= s1; ++s2) {
      
				auto bf2_first = shell2bf[s2];
				auto n2 = shells[s2].size();
      
				for (auto s3=0; s3 <= s1; ++s3) {
	
					auto bf3_first = shell2bf[s3];
					auto n3 = shells[s3].size();

					const auto s4_max = (s1 == s3) ? s2 : s3;
					for(auto s4=0; s4<=s4_max; ++s4, ++s1234) {
	  
						if (s1234 % nthreads != thread_id) continue; 
					
						auto bf4_first = shell2bf[s4];
						auto n4 = shells[s4].size();

						engine.compute(shells[s1], shells[s2], shells[s3], shells[s4]);

						const auto* buf_1234 = buf[0];
						if (buf_1234 == nullptr)
							continue;

						for(auto f1=0, f1234=0; f1!=n1; ++f1) {
            
							const auto bf1 = f1 + bf1_first;
            
							for(auto f2=0; f2!=n2; ++f2) {
              
								const auto bf2 = f2 + bf2_first;
             
								for(auto f3=0; f3!=n3; ++f3) {
                
									const auto bf3 = f3 + bf3_first;
             
									for(auto f4=0; f4!=n4; ++f4, ++f1234) {
                  
										const auto bf4 = f4 + bf4_first;
										eris(bf1, bf2, bf3, bf4) = buf_1234[f1234];
										eris(bf3, bf4, bf1, bf2) = buf_1234[f1234];
		  
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
	return eris;
}

void IntegralEngine::compute_eris_3index(const std::vector<libint2::Shell>& obs, const std::vector<libint2::Shell>& auxbs, Matrix& eris)
{
	
	using libint2::Shell;
	using libint2::Engine;
	using libint2::Operator;
	using libint2::BraKet; 

	const auto n_obs = nbasis(obs);
	const auto n_abs = nbasis(auxbs);
	int ndim = (n_obs * (n_obs + 1)) / 2;  
	eris = Matrix::Zero(ndim, n_abs);
	
	int nthreads = molecule->control->get_option<int>("nthreads");

	std::vector<Engine> engines(nthreads);
	engines[0] =
		Engine(Operator::coulomb, 
	std::max(max_nprim(obs), max_nprim(auxbs)),
	std::max(max_l(obs), max_l(auxbs)), 0, 
	std::numeric_limits<double>::epsilon(), 
	libint2::operator_traits<Operator::coulomb>::default_params(),
	BraKet::xs_xx);
	for (size_t i = 1; i != nthreads; ++i) {
		engines[i] = engines[0];
	}

	auto shell2bf_obs = map_shell_to_basis_function(obs);
	auto shell2bf_abs = map_shell_to_basis_function(auxbs);


	auto compute = [&](int thread_id) {
		auto& engine = engines[thread_id];
		const auto& buf = engine.results();
	
		// loop over shell quartets
		for (auto s1=0l, s123=0l; s1 != auxbs.size(); ++s1) {
    
			auto bf1_first = shell2bf_abs[s1];
			auto n1 = auxbs[s1].size();

			for (auto s2=0; s2 != obs.size(); ++s2) {
      
				auto bf2_first = shell2bf_obs[s2];
				auto n2 = obs[s2].size();
      
				for (auto s3=0; s3 <= s2; ++s3, ++s123) {
	
					if (s123 % nthreads != thread_id) continue; 
				
					auto bf3_first = shell2bf_obs[s3];
					auto n3 = obs[s3].size();

					engine.compute(auxbs[s1], obs[s2], obs[s3]);

					const auto* buf_123 = buf[0];
					if (buf_123 == nullptr)
						continue;

					for(auto f1=0, f123=0; f1!=n1; ++f1) {
            
						const auto bf1 = f1 + bf1_first;
            
						for(auto f2=0; f2!=n2; ++f2) {
              
							const auto bf2 = f2 + bf2_first;
							int ix = (bf2 * (bf2 + 1)) / 2; 
			 
							for(auto f3=0; f3!=n3; ++f3, ++f123) {
                
								const auto bf3 = f3 + bf3_first; 
								if (bf2 >= bf3)
									eris(ix+bf3, bf1)= buf_123[f123];							
							}
						}
					}
				}
			}
		}		
	};
	
	parallel_do(compute, nthreads); 
							
} 

void IntegralEngine::compute_eris_3index(const std::vector<libint2::Shell>& obs, const std::vector<libint2::Shell>& auxbs, SparseMatrix& eris)
{
	
	using libint2::Shell;
	using libint2::Engine;
	using libint2::Operator;
	using libint2::BraKet; 

	const auto n_obs = nbasis(obs);
	const auto n_abs = nbasis(auxbs);
	int ndim = (n_obs * (n_obs + 1)) / 2;  

	Engine engine(Operator::coulomb, 
	std::max(max_nprim(obs), max_nprim(auxbs)),
	std::max(max_l(obs), max_l(auxbs)), 0, 
	std::numeric_limits<double>::epsilon(), 
	libint2::operator_traits<Operator::coulomb>::default_params(),
	BraKet::xs_xx);

	auto shell2bf_obs = map_shell_to_basis_function(obs);
	auto shell2bf_abs = map_shell_to_basis_function(auxbs);
	
	using T = Eigen::Triplet<double>; 
	std::vector<T> tripletList; 
	eris.resize(ndim, n_abs);

	const auto& buf = engine.results();
	
	// loop over shell quartets
	for (auto s1=0; s1 != auxbs.size(); ++s1) {
    
		auto bf1_first = shell2bf_abs[s1];
		auto n1 = auxbs[s1].size();

		for (auto s2=0; s2 != obs.size(); ++s2) {
      
			auto bf2_first = shell2bf_obs[s2];
			auto n2 = obs[s2].size();
      
			for (auto s3=0; s3 <= s2; ++s3) {
	
				auto bf3_first = shell2bf_obs[s3];
				auto n3 = obs[s3].size();

				engine.compute(auxbs[s1], obs[s2], obs[s3]);

				const auto* buf_123 = buf[0];
				if (buf_123 == nullptr)
					continue;

				for(auto f1=0, f123=0; f1!=n1; ++f1) {
            
					const auto bf1 = f1 + bf1_first;
            
					for(auto f2=0; f2!=n2; ++f2) {
              
						const auto bf2 = f2 + bf2_first;
						int ix = (bf2 * (bf2 + 1)) / 2; 
			 
						for(auto f3=0; f3!=n3; ++f3, ++f123) {
                
							const auto bf3 = f3 + bf3_first; 
							if (bf2 >= bf3)
								tripletList.push_back(T(ix+bf3, bf1, buf_123[f123]));							
						}
					}
				}
			}
		}
	}						
	
	eris.setFromTriplets(tripletList.begin(), tripletList.end()); 
	eris.makeCompressed();
} 

void IntegralEngine::compute_eris_3index(const std::vector<libint2::Shell>& obs1, const std::vector<libint2::Shell>& obs2,
const std::vector<libint2::Shell>& auxbs, Matrix& eris)
{
	
	using libint2::Shell;
	using libint2::Engine;
	using libint2::Operator;
	using libint2::BraKet; 

	const auto n_obs_1 = nbasis(obs1);
	const auto n_obs_2 = nbasis(obs2);
	const auto n_abs = nbasis(auxbs);
	
	int ndim = n_obs_1 * n_obs_2; 
	eris = Matrix::Zero(ndim, n_abs);
	int nthreads = molecule->control->get_option<int>("nthreads");

	std::vector<Engine> engines(nthreads);
	engines[0] =
		Engine(Operator::coulomb, 
	std::max(max_nprim(obs1), std::max(max_nprim(obs2), max_nprim(auxbs))),
	std::max(max_l(obs1), std::max(max_l(obs2), max_l(auxbs))), 0, 
	std::numeric_limits<double>::epsilon(), 
	libint2::operator_traits<Operator::coulomb>::default_params(),
	BraKet::xs_xx);
	for (size_t i = 1; i != nthreads; ++i) {
		engines[i] = engines[0];
	}

	auto shell2bf_obs_1 = map_shell_to_basis_function(obs1);
	auto shell2bf_obs_2 = map_shell_to_basis_function(obs2);
	auto shell2bf_abs = map_shell_to_basis_function(auxbs);


	auto compute = [&](int thread_id) {
		auto& engine = engines[thread_id];
		const auto& buf = engine.results();
	
		// loop over shell quartets
		for (auto s1=0l, s123=0l; s1 != auxbs.size(); ++s1) {
    
			auto bf1_first = shell2bf_abs[s1];
			auto n1 = auxbs[s1].size();

			for (auto s2=0; s2 != obs1.size(); ++s2) {
      
				auto bf2_first = shell2bf_obs_1[s2];
				auto n2 = obs1[s2].size();
      
				for (auto s3=0; s3 != obs2.size(); ++s3, ++s123) {
					
					if (s123 % nthreads != thread_id) continue;
					
					auto bf3_first = shell2bf_obs_2[s3];
					auto n3 = obs2[s3].size();

					engine.compute(auxbs[s1], obs1[s2], obs2[s3]);

					const auto* buf_123 = buf[0];
					if (buf_123 == nullptr)
						continue;

					for(auto f1=0, f123=0; f1!=n1; ++f1) {
            
						const auto bf1 = f1 + bf1_first;
            
						for(auto f2=0; f2!=n2; ++f2) {
              
							const auto bf2 = f2 + bf2_first;
			 
							for(auto f3=0; f3!=n3; ++f3, ++f123) {
                
								const auto bf3 = f3 + bf3_first; 
								eris(bf2*n_obs_2 + bf3, bf1)= buf_123[f123];							
							}
						}
					}
				}
			}
		}	
	};	
	
	parallel_do(compute, nthreads); 				
} 
						
Matrix IntegralEngine::compute_eris_2index(const std::vector<libint2::Shell>& auxbs) {
	using libint2::Shell;
	using libint2::Engine;
	using libint2::Operator;
	using libint2::BraKet; 

	const auto n_abs = nbasis(auxbs); 
	Matrix eris = Matrix::Zero(n_abs, n_abs);
	int nthreads = molecule->control->get_option<int>("nthreads");

	std::vector<Engine> engines(nthreads);
	engines[0] =
		Engine(libint2::Operator::coulomb, max_nprim(auxbs), max_l(auxbs), 0);
	engines[0].set_braket(BraKet::xs_xs);
	for (size_t i = 1; i != nthreads; ++i) {
		engines[i] = engines[0];
	}

	auto shell2bf_abs = map_shell_to_basis_function(auxbs);
	int nshells = auxbs.size(); 
	
	auto compute = [&](int thread_id) {

		auto& engine = engines[thread_id];
		const auto& buf = engine.results();

		// loop over unique shell pairs, {s1,s2} such that s1 >= s2
		// this is due to the permutational symmetry of the real integrals over
		// Hermitian operators: (1|2) = (2|1)
		for (auto s1=0, s12=0; s1 != auxbs.size(); ++s1) {
    
			auto bf1_first = shell2bf_abs[s1];
			auto n1 = auxbs[s1].size();

			for (auto s2=0; s2 <= s1; ++s2, ++s12) {
      		  
			  	if (s12 % nthreads != thread_id) continue; 
				
				auto bf2_first = shell2bf_abs[s2];
				auto n2 = auxbs[s2].size();

				engine.compute(auxbs[s1], auxbs[s2]);

				const auto* buf_12 = buf[0];
				if (buf_12 == nullptr)
					continue;

				for(auto f1=0, f12=0; f1!=n1; ++f1) {
            
					const auto bf1 = f1 + bf1_first;
            
					for(auto f2=0; f2!=n2; ++f2, ++f12) {
              
						const auto bf2 = f2 + bf2_first;
           
             
						eris(bf1, bf2) = buf_12[f12];
						eris(bf2, bf1) = buf_12[f12];
		  
							
					}
				}
			}
		}
	};  // compute lambda

	parallel_do(compute, nthreads);
	
	return eris;
	/*using libint2::Shell;
	using libint2::Engine;
	using libint2::Operator;
	using libint2::BraKet; 

	const auto n_abs = nbasis(auxbs); 
	Matrix eris = Matrix::Zero(n_abs, n_abs);

	Engine engine(Operator::coulomb, max_nprim(auxbs), max_l(auxbs), 0);
	engine.set_braket(BraKet::xs_xs);

	auto shell2bf_abs = map_shell_to_basis_function(auxbs);

	const auto& buf = engine.results();

	// loop over shell quartets
	for (auto s1=0; s1 != auxbs.size(); ++s1) {
    
		auto bf1_first = shell2bf_abs[s1];
		auto n1 = auxbs[s1].size();

		for (auto s2=0; s2 <= s1; ++s2) {
      
			auto bf2_first = shell2bf_abs[s2];
			auto n2 = auxbs[s2].size();

			engine.compute(auxbs[s1], auxbs[s2]);

			const auto* buf_12 = buf[0];
			if (buf_12 == nullptr)
				continue;

			for(auto f1=0, f12=0; f1!=n1; ++f1) {
            
				const auto bf1 = f1 + bf1_first;
            
				for(auto f2=0; f2!=n2; ++f2, ++f12) {
              
					const auto bf2 = f2 + bf2_first;
           
             
					eris(bf1, bf2) = buf_12[f12];
					eris(bf2, bf1) = buf_12[f12];
		  
							
				}
			}
		}
	}
	
	return eris;	*/				
}

DFBlocks::DFBlocks(int _nfrags) {
	resize(_nfrags); 
}

void DFBlocks::resize(int _nfrags) {
	
	nfrags = _nfrags; 
	npairs = (nfrags*(nfrags+1))/2; 
	ntriples = npairs*nfrags; 
	eris.clear(); 
	zeroList.clear();
	
	eris.resize(ntriples);
	for (int n = 0; n < ntriples; n++)
		zeroList.push_back(true); 
	
}

Matrix& DFBlocks::operator()(int f1, int f2, int fk) {
	return eris[TRIPLE_INDEX(f1, f2, fk, npairs)]; 
} 

bool DFBlocks::isZero(int f1, int f2, int fk) const {
	return zeroList[TRIPLE_INDEX(f1, f2, fk, npairs)]; 
}

void DFBlocks::setZero(int f1, int f2, int fk) {
	zeroList[TRIPLE_INDEX(f1, f2, fk, npairs)] = true; 
}

void DFBlocks::setNonZero(int f1, int f2, int fk) {
	zeroList[TRIPLE_INDEX(f1, f2, fk, npairs)] = false; 
}


