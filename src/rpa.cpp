
#include "rpa.hpp"
#include "eigen_wrapper.hpp"
#include <ctf.hpp>
#include <libint2.hpp>
#include "logger.hpp"

// Constructor
RPA::RPA(Command& c, Fock& _focker, int _N, int _nocc, CTF::World& _dw) : cmd(c), focker(_focker), N(_N), nocc(_nocc), nvirt(_N - _nocc), dw(_dw)
{
	energy = 0.0;
}

void RPA::compute(bool print) {
	
	Logger& log = focker.getMolecule()->control->log; 
	bool sosex = cmd.get_option<bool>("sosex"); 
	
	auto &shells = focker.getMolecule()->getBasis().getIntShells();
	std::vector<CTF::Tensor<> > V;
	
	Matrix cp_occ = focker.getCP().block(0, 0, N, nocc); 
	Matrix cp_virt = focker.getCP().block(0, nocc, N, nvirt); 
	
	if (print) log.print("Transforming integrals"); 
	if (cmd.get_option<bool>("longrange")) {
		if(print) log.print("Using long-range potential, with mu = " + std::to_string(cmd.get_option<double>("mu"))); 
		eris(shells, V, cp_occ, cp_virt, true);
	} else { 
		eris(shells, V, cp_occ, cp_virt); 
	}
	if (print) log.localTime();  
	
	if (print) log.print("\nForming excitation matrices"); 
	int64_t sz1, sz2, *i1, *i2;
	double *viajb;
	V[0].read_local(&sz1, &i1, &viajb);
	
	int dim = nvirt*nocc; 
	int nvnono = nocc*dim; 
	Matrix A = Matrix::Zero(dim, dim); 
	Matrix B = Matrix::Zero(dim, dim);
	Matrix K = Matrix::Zero(dim, dim); 
	
	Matrix F = focker.getCP().transpose() * focker.getFockAO() * focker.getCP(); 
	
	int axis = nocc*(nocc+1);
	axis /= 2; 
	int ix1, ix2;  
	for (int i = 0; i < nocc; i++) {
		for (int a = 0; a < nvirt; a++) {
			ix1 = i*nvirt+a;
			
			for (int j = 0; j < nocc; j++) {
				for (int b = 0; b < nvirt; b++) {
					ix2 = j*nvirt + b;
					
					K(ix1, ix2) = 2.0 * viajb[i + a*nocc + j*dim + b*nvnono]; 
					A(ix1, ix2) = K(ix1, ix2);
					B(ix1, ix2) = K(ix1, ix2);
					if (sosex) {
						B(ix1, ix2) -= viajb[j + a*nocc + i*dim + b*nvnono]; 
					}
					if (i==j) A(ix1, ix2) += F(a+nocc, b+nocc); 
					if (a==b) A(ix1, ix2) -= F(i, j); 
				}
			}
		}
	}
	delete i1;
	delete viajb;
	if (print) log.localTime(); 
	
	if (print) log.print("\nSolving Riccatti equations");
	if (cmd.get_option<bool>("longrange") || cmd.get_option<bool>("iterative")) {
		
		Matrix Ap = A; 
		for (int i = 0; i < dim ;i++) Ap(i, i) = 0.0; 
		
		Matrix T = -K; 
		
		double delta = 1.0;
		int iter = 0;
		Matrix newT; 
		while (delta > 1e-4 && iter < 15) {
			newT = -(K + T*Ap);  
			newT -= (T*K + Ap)*T; 
			
			for (int p = 0; p < dim; p++)
				for (int q = 0; q < dim; q++)
					newT(p, q) /= A(p, p) + A(q, q); 
			
			delta = (newT - T).norm(); 
			iter++; 
			T = newT; 
		}
		
		energy = 0.5 * (B*T).trace(); 
		
	} else {
		Matrix M = (A+B)*(A-B); 
		Eigen::EigenSolver<Matrix> solver(M);
		const Eigen::VectorXcd& evals = solver.eigenvalues(); 
		Vector omega = Vector::Zero(evals.size()); 
	
		energy = 0.0; 
		for (int i = 0; i < omega.size(); i++) { 
			std::complex<double> ei = evals[i]; 
			omega[i] = std::sqrt(ei.real()); 
			energy += omega[i] - A(i, i); 
		}	
		if(print) log.print("\nRPA excitation energies: ");
		if(print) log.print(omega); 
	
		energy /= 2.0;
		if (sosex) energy /= 2.0; 
	}
	if (print) log.localTime(); 
}

void RPA::eris(const std::vector<libint2::Shell>& shells, std::vector<CTF::Tensor<> >& moInts, Matrix& T, Matrix& V, bool longrange) {
	using libint2::Shell;
	using libint2::Engine;
	using libint2::Operator;

	IntegralEngine& integrals = focker.getIntegrals(); 

	const auto n = integrals.nbasis(shells);
	int l2 = (n*(n+1))/2;
	
	int aoshape[4] = {N, N, N, N};  
	int iajbshape[4] = {nocc, nvirt, nocc, nvirt};
	int ijabshape[4] = {nocc, nocc, nvirt, nvirt}; 
	int sysy[4] = {SY, NS, SY, NS}; 
	int nsns[4] = {NS, NS, NS, NS};
	int syns[4] = {SY, NS, NS, NS}; 
	int nssy[4] = {NS, NS, SY, NS};
	
	CTF::Tensor<> aoInts(4, aoshape, sysy, dw); 
	
	int64_t sz, *ix; 
	double *values; 
	aoInts.read_local(&sz, &ix, &values); 
	
	double mu = cmd.get_option<double>("mu"); 
	
	Engine* engine; 
	
	if (longrange) {
		engine = new Engine(Operator::erf_coulomb, integrals.max_nprim(shells), integrals.max_l(shells), 0);
		engine->set_params(mu); 
	} else {
		engine = new Engine(Operator::coulomb, integrals.max_nprim(shells), integrals.max_l(shells), 0); 
	}
	
	auto shell2bf = integrals.map_shell_to_basis_function(shells);

	const auto& buf = engine->results();

	// loop over shell quartets
	for (auto s1=0; s1 != shells.size(); ++s1) {
    
		auto bf1_first = shell2bf[s1];
		auto n1 = shells[s1].size();

		for (auto s2=0; s2 <= s1; ++s2) {
      
			auto bf2_first = shell2bf[s2];
			auto n2 = shells[s2].size();
      
			for (auto s3=0; s3 <= s1; ++s3) {
	
				auto bf3_first = shell2bf[s3];
				auto n3 = shells[s3].size();

				const auto s4_max = (s1 == s3) ? s2 : s3;
				for(auto s4=0; s4<=s4_max; ++s4) {
	  
					auto bf4_first = shell2bf[s4];
					auto n4 = shells[s4].size();

					engine->compute(shells[s1], shells[s2], shells[s3], shells[s4]);

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
									int mu = std::min(bf1, bf2);
									int nu = std::max(bf1, bf2);
									int sigma = std::min(bf3, bf4);
									int lambda = std::max(bf3, bf4); 
									int l1 = (nu*(nu+1))/2;
									int l3 = (lambda*(lambda+1))/2; 
									values[mu + l1 + sigma*l2 + l3*l2] = buf_1234[f1234]; 
									values[sigma + l3 + mu*l2 + l1*l2] = buf_1234[f1234]; 
								}
							}
						}
					}
				}
			}
		}
	}
	aoInts.write(sz, ix, values);
	delete ix; 
	delete values; 
	delete engine;
	
	int64_t sz1, sz2, *ix1, *ix2;
	double *v1, *v2; 
	
	CTF::Matrix<> cp_occ(N, nocc, NS, dw); 
	CTF::Matrix<> cp_virt(N, nvirt, NS, dw);  

	int ctr = 0; 
	cp_occ.read_local(&sz1, &ix1, &v1);
	for (int i = 0; i < nocc; i++)
		for (int mu = 0; mu < N; mu++) 
			v1[ctr++] = T(mu, i); 
	cp_occ.write(sz1, ix1, v1); 
	
	delete ix1;
	delete v1; 
	
	ctr = 0; 
	cp_virt.read_local(&sz2, &ix2, &v2);
	for (int a = 0; a < nvirt; a++)
		for (int mu = 0; mu < N; mu++)
			v2[ctr++] = V(mu, a); 
	cp_virt.write(sz2, ix2, v2); 
	
	delete ix2; 
	delete v2;
	
	CTF::Tensor<> Viajb(4, iajbshape, nsns, dw); 
	{
		int t1_lens[4] = {N, N, N, nvirt};
		int t2_lens[4] = {N, N, nocc, nvirt};
		int t3_lens[4] = {N, nvirt, nocc, nvirt}; 
		CTF::Tensor<> temp1(4, t1_lens, syns, dw);
		CTF::Tensor<> temp2(4, t2_lens, syns, dw); 
		CTF::Tensor<> temp3(4, t3_lens, nsns, dw);	
		
		temp1["uvwb"] = cp_virt["xb"] * aoInts["uvwx"]; 
		temp2["uvjb"] = cp_occ["wj"] * temp1["uvwb"]; 
		temp3["uajb"] = cp_virt["va"] * temp2["uvjb"];
		Viajb["iajb"] = cp_occ["ui"] * temp3["uajb"]; 
	}  
	moInts.push_back(Viajb);
}

double divide_r(double a, double b){
	return a/b;
}

void RPA::fcompute(fInfo& info, bool print) {
	
	Logger& log = focker.getMolecule()->control->log; 
	int xcorder = cmd.get_option<int>("rpax");
	bool sosex = xcorder > 0; 
	bool rpax =  xcorder > 1; 
	int nfrags = info.nocc.size();
	
	// Form MO overlap matrix
	int i_offset = 0; int mu_offset = 0; 
	Matrix sigma(nocc, nocc);
	int f1_nocc, f2_nocc; 
	for (int i = 0; i < nfrags; i++) {
		f1_nocc = info.nocc[i]; 
		int f1_nbfs = f1_nocc + info.nvirt[i]; 
		
		int j_offset = 0; int nu_offset = 0;
		for (int j = 0; j <= i; j++) { 
			f2_nocc = info.nocc[j];
			int f2_nbfs = f2_nocc + info.nvirt[j]; 
			
			sigma.block(i_offset, j_offset, f1_nocc, f2_nocc) = info.T.block(mu_offset, i_offset, f1_nbfs, f1_nocc).transpose() 
				* info.S.block(mu_offset, nu_offset, f1_nbfs, f2_nbfs) * info.T.block(nu_offset, j_offset, f2_nbfs, f2_nocc);

			j_offset += f2_nocc;
			nu_offset += f2_nbfs;
		}
		
		mu_offset += f1_nbfs; 
		i_offset += f1_nocc; 
	}
	sigma = sigma.selfadjointView<Eigen::Lower>();
	sigma = sigma.inverse(); 
	
	// Project virtual orbitals out of occupied subspace
	Matrix P(N, N); 
	i_offset = 0; mu_offset = 0;
	for (int i = 0; i < nfrags; i++) {
		f1_nocc = info.nocc[i];
		int f1_nbfs = f1_nocc + info.nvirt[i];  
		
		int j_offset = 0; int nu_offset = 0;
		for (int j = 0; j <= i; j++) {
			f2_nocc = info.nocc[j];
			int f2_nbfs = f2_nocc + info.nvirt[j]; 
			
			P.block(mu_offset, nu_offset, f1_nbfs, f2_nbfs) = info.T.block(mu_offset, i_offset, f1_nbfs, f1_nocc) 
				* sigma.block(i_offset, j_offset, f1_nocc, f2_nocc) * info.T.block(nu_offset, j_offset, f2_nbfs, f2_nocc).transpose();
			
			j_offset += f2_nocc;
			nu_offset += f2_nbfs;
		}
		
		mu_offset += f1_nbfs; 
		i_offset += f1_nocc; 
	}
	P = P.selfadjointView<Eigen::Lower>();
	
	sigma = Matrix::Identity(N, N) - P * info.S; 
	info.V = sigma * info.V; 
	
	int offset = 0; 
	if(!focker.getMolecule()->control->get_option<bool>("withcore")) {
		for (int i = 0; i < nfrags; i++) 
			offset += info.ncore[i];  
		nocc -= offset; 
		N = nocc + nvirt; 
	}

	int offnocc = nocc + offset; 
	int offN = N + offset; 
	
	int dim = nvirt*nocc; 
	int nvnono = nocc*dim;
	
	int NSNS[4] = {NS, NS, NS, NS}; 
	int iajb_dims[4] = {nocc, nvirt, nocc, nvirt}; 
	
	CTF::Tensor<> A(4, iajb_dims, NSNS, dw); 
	CTF::Tensor<> Ap(4, iajb_dims, NSNS, dw); 
	CTF::Tensor<> B(4, iajb_dims, NSNS, dw); 
	CTF::Tensor<> K(4, iajb_dims, NSNS, dw);   

	double *avals, *apvals, *bvals, *kvals; 
	int64_t asz, apsz, bsz, ksz; 
	avals = A.get_raw_data(&asz); 
	apvals = Ap.get_raw_data(&apsz); 
	bvals = B.get_raw_data(&bsz); 
	kvals = K.get_raw_data(&ksz); 
	
	Matrix F = Matrix::Zero(offN, offN); 
	F.block(0, 0, offnocc, offnocc) = info.T.transpose() * info.F * info.T; 
	F.block(0, offnocc, offnocc, nvirt) = info.T.transpose() * info.F * info.V; 
	F.block(offnocc, 0, nvirt, offnocc) = F.block(0, offnocc, offnocc, nvirt).transpose();
	F.block(offnocc, offnocc, nvirt, nvirt) = info.V.transpose() * info.F * info.V; 
	
	// Compute and transform eris
	auto &shells = info.shells;
	bool density_fitted = cmd.get_option<bool>("df"); 
	int ix1, ix2, ix3; 	
	
	if (print) log.print("Transforming integrals"); 
	if(!density_fitted) {
		std::vector<CTF::Tensor<> > V;
	
		if (cmd.get_option<bool>("longrange")) {
			if(print) log.print("Using long-range potential, with mu = " + std::to_string(cmd.get_option<double>("mu"))); 
			eris(shells, V, info.T, info.V, true);
		} else { 
			eris(shells, V, info.T, info.V); 
		}
		if (print) {
			log.localTime();  
			log.print("\nForming excitation matrices"); 
		}
		
		int64_t sz1, sz2; 
		double *viajb = V[0].get_raw_data(&sz1);
	
		int axis = nocc*(nocc+1);
		axis /= 2;  
		for (int i = 0; i < nocc; i++) {
			for (int a = 0; a < nvirt; a++) {
				ix1 = i*nvirt+a;
			
				for (int j = 0; j < nocc; j++) {
					for (int b = 0; b < nvirt; b++) {
						ix2 = j*nvirt + b; 
						if (ix2 > ix1) continue; 
						ix2 = i + a*nocc + j*dim + b*nvnono; 
						int offix2 = ix2 + offset*(1+dim); 
						avals[ix2] = bvals[ix2] = kvals[ix2] = 2.0 * viajb[offix2]; 
						if (sosex) {
							double xc = viajb[j + a*nocc + i*dim + b*nvnono + offset*(1 + dim)]; 
							bvals[ix2] -= xc; 
							if (rpax) {
								avals[ix2] -= xc;
								kvals[ix2] -= xc; 
							}
						}
						if (i==j) avals[ix2] += F(a+offnocc, b+offnocc); 
						if (a==b) avals[ix2] -= F(i+offset, j+offset); 
						
						ix3 = j + b*nocc + i*dim + a*nvnono; 
						kvals[ix3] = kvals[ix2]; 
						avals[ix3] = avals[ix2];
						bvals[ix3] = bvals[ix2]; 
						
						if (!(i == j && a == b)) apvals[ix2] = apvals[ix3] = avals[ix2]; 
					}
				}
			}
		}
	} else {
		Matrix L; 
		nocc += offset;
		N += offset; 
		df_eris(shells, info.df_shells, info.T, info.V, L); 
		nocc -= offset; 
		N -= offset; 
		Matrix K = L * L.transpose(); 
		L.resize(0, 0); 
		
		if (print) {
			log.localTime();  
			log.print("\nForming excitation matrices"); 
		} 
		
		int ov = offset * nvirt; 
		for (int i = 0; i < nocc; i++) {
			for (int a = 0; a < nvirt; a++) {
				ix1 = i*nvirt+a;
			
				for (int j = 0; j < nocc; j++) {
					for (int b = 0; b < nvirt; b++) {
						ix2 = j*nvirt + b;
						if (ix2 > ix1) continue;
					
						double viajb = K(ix1+ov, ix2+ov); 
						
						ix2 = i + a*nocc + j*dim + b*nvnono; 
						avals[ix2] = bvals[ix2] = kvals[ix2] = 2.0 * viajb; 
						if (sosex) {
							double vjaib = K(j*nvirt+a+ov, i*nvirt+b+ov); 
							bvals[ix2] -= vjaib;
							if (rpax) {
								avals[ix2] -= vjaib; 
								kvals[ix2] -= vjaib; 
							}
						}
						if (i==j) avals[ix2] += F(a+offnocc, b+offnocc); 
						if (a==b) avals[ix2] -= F(i+offset, j+offset); 
						
						ix3 = j + b*nocc + i*dim + a*nvnono; 
						kvals[ix3] = kvals[ix2]; 
						avals[ix3] = avals[ix2];
						bvals[ix3] = bvals[ix2]; 
						
						if (!(i == j && a == b)) {
							apvals[ix2] = apvals[ix3] = avals[ix2]; 
							avals[ix2] = avals[ix3] = 0.0; 
						}
					}
				}
			}
		}
	}
	if (print) log.localTime(); 
	
	if (print) log.print("\nSolving Riccatti equations"); 
	
	CTF::Tensor<> T(4, iajb_dims, NSNS, dw);
	T["iajb"] = A["iaia"] + A["jbjb"]; 
	A["iajb"] = T["iajb"]; 
	CTF::Tensor<> newT(4, iajb_dims, NSNS, dw);
	CTF::Function<> fctr(&divide_r);
	T["iajb"] = -K["iajb"];
	 
	double delta = 1.0;
	int iter = 0; 
	double tnorm = T.norm2(); 
	while (delta > 1e-4 && iter < 30) {
		newT["iajb"] = -Ap["iajb"]; 
		newT["iajb"] -= 0.5 * T["iakc"] * K["kcjb"]; 
		newT["iajb"] = newT["iakc"] * T["kcjb"]; 
		newT["iajb"] += newT["jbia"];
		newT["iajb"] -= K["iajb"];

		newT.contract(1.0, newT, "iajb", A, "iajb", 0.0, "iajb", fctr);

		delta = -tnorm;
		tnorm = newT.norm2();
		delta += tnorm; 
		delta = fabs(delta); 
		
		iter++; 
		T["iajb"] = newT["iajb"]; 
	}
		
	if (print) {
		log.localTime();
		log.print("Decomposing energy contributions...");
	}
	
	energy = 0.0; 
	info.eintra = info.edisp = info.edispexch = info.eionic = info.ebsse = 0.0; 

	std::vector<int> cum_occ, cum_virt; 
	int currnocc = 0; int currvirt = 0;
	for (int i = 0; i < nfrags; i++) {
		currnocc += info.nocc[i] - info.ncore[i]; 
		currvirt += info.nvirt[i];
		cum_occ.push_back(currnocc);
		cum_virt.push_back(currvirt); 
	}
		
	double *tvals = T.get_raw_data(&asz); 
	
	int d1, d2, d3, d4, d; 
	int ifrag, jfrag, afrag, bfrag; 
	double scale = 1.0; 
	for (int i = 0; i < nocc; i++) {	
			
		for (int n = 0; n < nfrags; n++) {
			if(i < cum_occ[n]) { ifrag = n; break; }
		}
		
		for (int j = 0; j <= i; j++) {
			scale = i == j ? 0.5 : 1.0; 
			
			for (int n = 0; n < nfrags; n++) {
				if(j < cum_occ[n]) { jfrag = n; break; }
			}
		
			for (int a = 0; a < nvirt; a++) {
							
				for (int n = 0; n < nfrags; n++) {
					if(a < cum_virt[n]) { afrag = n; break; }
				}
							
				d1 = ifrag == afrag ? 0 : 1; 
				d2 = jfrag == afrag ? 0 : 1; 
				ix1 = i*nvirt + a; 

				for (int b = 0; b < nvirt; b++) {
					
					for (int n = 0; n < nfrags; n++) {
						if(b < cum_virt[n]) { bfrag = n; break; }
					}
								
					d3 = jfrag == bfrag ? 0 : 1;
					d4 = ifrag == bfrag ? 0 : 1;  
								
					d = 8*d1 + 4*d2 + 2*d3 + d4; 
					ix2 = i + a*nocc + j*dim + b*nvnono;
							
					switch(d) {
						case 0: {
							info.eintra += scale * tvals[ix2] * bvals[ix2]; 
							break;
						}
										
						case 5: {
							info.edisp += tvals[ix2] * bvals[ix2]; 
							break;
						}
										
						case 10: {
							info.edispexch += tvals[ix2] * bvals[ix2]; 
							break;
						}
						case 15: {
							info.ebsse += scale * tvals[ix2] * bvals[ix2]; 
							break;
						}
										
						default: {
							info.eionic += scale * tvals[ix2] * bvals[ix2]; 
						}
					}
									
				}
			}
		}

	}
	
	if (rpax) {
		info.edisp *= 0.5;
		info.edispexch *= 0.5;
		info.eionic *= 0.5;
		info.ebsse *= 0.5;
		info.eintra *= 0.5;
	} 
	energy = info.edisp + info.edispexch;
		
	if (print) log.localTime(); 
}

void RPA::df_eris(const std::vector<libint2::Shell>& obs, const std::vector<libint2::Shell>& auxbs, Matrix& T, Matrix& V, Matrix& BmnP) {
	IntegralEngine& integrals = focker.getIntegrals(); 
	
	// Calculate (P|Q) and (P|mn) eris
	Matrix JPQ = integrals.compute_eris_2index(auxbs);  
	Matrix KmnP;
	integrals.compute_eris_3index(obs, auxbs, KmnP);

	// Form inverse J-metric
	EigenSolver es(JPQ); 
	JPQ = es.operatorInverseSqrt(); 
	
	int nobs = integrals.nbasis(obs);
	int nabs = JPQ.rows(); 

	int nvirt = N-nocc; 
	int offset = 0;
	int offN = offset + N; 
	int offnocc = offset + nocc; 
	
	BmnP = Matrix::Zero(nocc*nobs, nabs); 
	int ix; 
	for (int i = 0; i < nocc; i++) {
		for (int nu = 0; nu < nobs; nu++) {
			for (int P = 0; P < nabs; P++) {
				
				ix = i*nobs+nu; 
				for (int mu = 0; mu < nobs; mu++) {
					int X = std::max(mu, nu); 
					int Y = std::min(mu, nu); 
					BmnP(ix, P) +=  T(mu, i) * KmnP((X*(X+1))/2 + Y, P); 
				}
				
			} 
		}
	}
	
	KmnP = Matrix::Zero(nocc*nvirt, nabs); 
	for (int i = 0; i < nocc; i++) {
		for (int a = 0; a < nvirt; a++) {
			for (int P = 0; P < nabs; P++) {
				
				ix = i*nvirt+a; 
				for (int nu = 0; nu < nobs; nu++)
					KmnP(ix, P) += V(nu, a) * BmnP(i*nobs+nu, P); 
				
			}
		}
	}
	
	BmnP.noalias() =  KmnP * JPQ; 
}
