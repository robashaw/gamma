
#include "mp2.hpp"
#include "logger.hpp"
#include "eigen_wrapper.hpp"
#include "molecule.hpp"
#include "integrals.hpp"
#include "error.hpp"
#include <iostream>
#include <thread>
#include "ProgramController.hpp"
#include <libint2.hpp>

// Constructor
MP2::MP2(Fock& _focker) : spinBasis(false), focker(_focker)
{
	N = focker.getDens().rows();
	nocc = focker.getMolecule()->getNel()/2;
	energy = 0.0;
	moInts.assign(N, 0.0);
}

double dividef(double a, double b){
	return a/b;
}

// Integral transformation
void MP2::transformIntegrals(bool withSpin)
{
	if (focker.getMolecule()->control->get_option<bool>("direct")) {
		Error e("MP2TRANS", "Integral direct MP2 not implemented yet.");
		focker.getMolecule()->control->log.error(e);
		nocc = 0;
	} else {
		
		// Multithread
		int nthreads = focker.getMolecule()->control->get_option<int>("nthreads");

		std::vector<Tensor4> moTemps(nthreads);
		std::vector<std::thread> thrds(nthreads);
		std::vector<int> startPoints(nthreads);
		std::vector<int> endPoints(nthreads);

		int spacing = N/nthreads;
		spacing = (double(N)/double(nthreads) - spacing > 0.5 ? spacing+1 : spacing);
		startPoints[0] = 0;
		endPoints[0] = spacing;
	
		for (int i = 0; i < nthreads; i++){
			moTemps[i] = Tensor4(endPoints[i]-startPoints[i], N, N, N, 0.0);
			thrds[i] = std::thread(&MP2::transformThread, *this, startPoints[i], endPoints[i], std::ref(moTemps[i]));
			if (i < nthreads - 1) { 
				startPoints[i+1] = endPoints[i];
				endPoints[i+1] = (i == (nthreads - 2) ? N : endPoints[i] + spacing);
			}
		}
		for (int i = 0; i < nthreads; i++) {
			thrds[i].join();

			// Copy in temporary matrix to moInts
			for (int p = startPoints[i]; p < endPoints[i]; p++){
				for (int q = 0; q < N; q++){
					for (int r = 0; r < N; r++){
						for (int s = 0; s < N; s++){
							moInts(p, q, r, s) = moTemps[i](p-startPoints[i], q, r, s);
						}
					}
				}
			}
		}
		
		focker.getIntegrals().clearTwoInts();
	}
}

void MP2::transformThread(int start, int end, Tensor4& moTemp)
{
	Matrix& C = focker.getCP();
	IntegralEngine& aoInts = focker.getIntegrals();

	int offset = end - start;
	Tensor4 temp1(offset, N, N, N, 0.0);
	Tensor4 temp2(offset, N, N, N, 0.0);
	Tensor4 temp3(offset, N, N, N, 0.0);

	// Transform as four 'quarter-transforms'
	for (int p = start; p < end; p++){
		
		for (int mu = 0; mu < N; mu++){
			for (int a = 0; a < N; a++){
				for (int b = 0; b < N; b++){
					for (int c = 0; c < N; c++)
						temp1(p-start, a, b, c) += C(mu, p)*aoInts.getERI(mu, a, b, c);
				} // b
			} // a
		} // mu
		
		for (int q = 0; q < N; q++){
			for (int nu = 0; nu < N; nu++){
				for (int b = 0; b < N; b++){
					for (int c = 0; c < N; c++)
						temp2(p-start, q, b, c) += C(nu, q)*temp1(p-start, nu, b, c);
				} // b
			} // nu

			for (int r = 0; r < N; r++){
				for (int lam = 0; lam < N; lam++){
					for (int c = 0; c < N; c++)
						temp3(p-start, q, r, c) += C(lam, r)*temp2(p-start, q, lam, c);
				} // lam
				
				for (int s = 0; s < N; s++){
					for (int sig = 0; sig < N; sig++)
						moTemp(p-start, q, r, s) += C(sig, s)*temp3(p-start, q, r, sig);
				} // s
			} // r
		} // q
	} // p
}

void MP2::cctrans() {
	int nvirt = N-nocc; 
	int offset = 0;

	if(!focker.getMolecule()->control->get_option<bool>("withcore")) {
		offset = nocc - focker.getMolecule()->getNValence() / 2; 
		nocc -= offset; 
		N = nocc + nvirt; 
	}

	int offnocc = nocc + offset; 
	int offN = N + offset; 

	focker.getMolecule()->control->log.print("No. of occ. orbitals: " + std::to_string(nocc));
	focker.getMolecule()->control->log.print("No. of virt. orbitals: " + std::to_string(nvirt)); 
	focker.getMolecule()->control->log.flush();
	
	int ao_lens[4] = {offN, offN, offN, offN};
	int sysy[4] = {SY, NS, SY, NS}; 
	int nssy[4] = {NS, NS, SY, NS};
	int nsns[4] = {NS, NS, NS, NS};
	int syns[4] = {SY, NS, NS, NS}; 
	CTF::Tensor<> ao_integrals(4, ao_lens, sysy, dw); 
	int64_t sz, *indices;
	double *values;
	int ctr;
	
	IntegralEngine& aoInts = focker.getIntegrals();
	
	ao_integrals.read_local(&sz, &indices, &values);
	ctr = 0;
	for (int a = 0; a < offN; a++) 
		for (int b = 0; b <= a; b++)
			for (int c = 0; c < offN; c++)
				for (int d = 0; d <= c; d++)
					values[ctr++] = aoInts.getERI(d, c, b, a); 
	ao_integrals.write(sz, indices, values);
	
	focker.getIntegrals().clearTwoInts();
	
	CTF::Matrix<> cp_occ(offN, nocc, NS, dw); 
	CTF::Matrix<> cp_virt(offN, nvirt, NS, dw); 
	Matrix& CP = focker.getCP(); 
	
	ctr = 0; 
	cp_occ.read_local(&sz, &indices, &values);
	for (int i = offset; i < offnocc; i++)
		for (int mu = 0; mu < offN; mu++)
			values[ctr++] = CP(mu, i); 
	cp_occ.write(sz, indices, values); 
	
	ctr = 0; 
	cp_virt.read_local(&sz, &indices, &values);
	for (int a = offnocc; a < offN; a++)
		for (int mu = 0; mu < offN; mu++)
			values[ctr++] = CP(mu, a); 
	cp_virt.write(sz, indices, values); 
	
	spinInts = std::make_shared<Integrals>(nocc, nvirt,  dw);
	Integrals& V = *spinInts; 
	
	// trans ABCD
	{
		int t1_lens[4] = {offN, offN, offN, nvirt}; 
		int t2_lens[4] = {offN, offN, nvirt, nvirt};
		int t3_lens[4] = {offN, nvirt, nvirt, nvirt};
		CTF::Tensor<> temp1(4, t1_lens, syns, dw);
		CTF::Tensor<> temp2(4, t2_lens, sysy, dw); 
		CTF::Tensor<> temp3(4, t3_lens, nssy, dw);	
		
		temp1["uvwd"] = cp_virt["xd"] * ao_integrals["uvwx"]; 
		temp2["uvcd"] = cp_virt["wc"] * temp1["uvwd"]; 
		temp3["ubcd"] = cp_virt["vb"] * temp2["uvcd"];
		V["abcd"] = 0.25 * cp_virt["ua"] * temp3["ubcd"]; 
	} 
	
	// trans ABCI, ABIC, AIBC, IABC
	{
		int t1_lens[4] = {offN, offN, offN, nocc}; 
		int t2_lens[4] = {offN, offN, nvirt, nocc};
		int t3_lens[4] = {offN, nvirt, nvirt, nocc};
		CTF::Tensor<> temp1(4, t1_lens, syns, dw);
		CTF::Tensor<> temp2(4, t2_lens, syns, dw); 
		CTF::Tensor<> temp3(4, t3_lens, nsns, dw);	
		
		temp1["uvwi"] = cp_occ["xi"] * ao_integrals["uvwx"]; 
		temp2["uvci"] = cp_virt["wc"] * temp1["uvwi"]; 
		temp3["ubci"] = cp_virt["vb"] * temp2["uvci"];
		V["abci"] = 0.5 * cp_virt["ua"] * temp3["ubci"]; 
		V["aibc"] = V["bcai"];
		V["iabc"] = V["bcai"];
		V["abic"] = V["abci"]; 
 	} 
	
	// trans ABIJ, IJAB
	{
		int t1_lens[4] = {offN, offN, offN, nocc}; 
		int t2_lens[4] = {offN, offN, nocc, nocc};
		int t3_lens[4] = {offN, nvirt, nocc, nocc};
		CTF::Tensor<> temp1(4, t1_lens, syns, dw);
		CTF::Tensor<> temp2(4, t2_lens, sysy, dw); 
		CTF::Tensor<> temp3(4, t3_lens, nssy, dw);	
		
		temp1["uvwj"] = cp_occ["xj"] * ao_integrals["uvwx"]; 
		temp2["uvij"] = cp_occ["wi"] * temp1["uvwj"]; 
		temp3["ubij"] = cp_virt["vb"] * temp2["uvij"];
		V["abij"] = 0.25 * cp_virt["ua"] * temp3["ubij"]; 
		V["ijab"] = V["abij"];
 	} 
	
	// trans AIBJ, AIJB, IABJ, IAJB
	{
		int t1_lens[4] = {offN, offN, offN, nocc}; 
		int t2_lens[4] = {offN, offN, nvirt, nocc};
		int t3_lens[4] = {offN, nocc, nvirt, nocc};
		CTF::Tensor<> temp1(4, t1_lens, syns, dw);
		CTF::Tensor<> temp2(4, t2_lens, syns, dw); 
		CTF::Tensor<> temp3(4, t3_lens, nsns, dw);	
		
		temp1["uvwj"] = cp_occ["xj"] * ao_integrals["uvwx"]; 
		temp2["uvbj"] = cp_virt["wb"] * temp1["uvwj"]; 
		temp3["uibj"] = cp_occ["vi"] * temp2["uvbj"];
		V["aibj"] = cp_virt["ua"] * temp3["uibj"]; 
		V["aijb"] = V["aibj"];
		V["iabj"] = V["aibj"];
		V["iajb"] = V["aibj"]; 
 	} 
	
	// trans AIJK, IAJK, IJAK, IJKA
	{
		int t1_lens[4] = {offN, offN, offN, nocc}; 
		int t2_lens[4] = {offN, offN, nocc, nocc};
		int t3_lens[4] = {offN, nocc, nocc, nocc};
		CTF::Tensor<> temp1(4, t1_lens, syns, dw);
		CTF::Tensor<> temp2(4, t2_lens, sysy, dw); 
		CTF::Tensor<> temp3(4, t3_lens, nssy, dw);	
		
		temp1["uvwk"] = cp_occ["xk"] * ao_integrals["uvwx"]; 
		temp2["uvjk"] = cp_occ["wj"] * temp1["uvwk"]; 
		temp3["uijk"] = cp_occ["vi"] * temp2["uvjk"];
		V["aijk"] = 0.5 * cp_virt["ua"] * temp3["uijk"]; 
		V["iajk"] = V["aijk"];
		V["ijak"] = V["akij"];
		V["ijka"] = V["akij"]; 
 	} 
	
	// trans IJKL
	{
		int t1_lens[4] = {offN, offN, offN, nocc}; 
		int t2_lens[4] = {offN, offN, nocc, nocc};
		int t3_lens[4] = {offN, nocc, nocc, nocc};
		CTF::Tensor<> temp1(4, t1_lens, syns, dw);
		CTF::Tensor<> temp2(4, t2_lens, sysy, dw); 
		CTF::Tensor<> temp3(4, t3_lens, nssy, dw);	
		
		temp1["uvwl"] = cp_occ["xl"] * ao_integrals["uvwx"]; 
		temp2["uvkl"] = cp_occ["wk"] * temp1["uvwl"]; 
		temp3["ujkl"] = cp_occ["vj"] * temp2["uvkl"];
		V["ijkl"] = 0.25 * cp_occ["ui"] * temp3["ujkl"]; 
 	} 
	
	Vector& eps = focker.getEps();
	amplitudes = std::make_shared<Amplitudes>(nocc, nvirt, dw); 
	Amplitudes& T = *amplitudes;
	
	int dlens[4] = {nocc, nocc, nvirt, nvirt}; 
	CTF::Tensor<> Dijab(4, dlens, sysy, dw); 
	ctr = 0;
	Dijab.read_local(&sz, &indices, &values); 
	for (int b = offnocc; b < offN; b++)
		for (int a = offnocc; a <= b; a++)
			for (int j = offset; j < offnocc; j++)
				for (int i = offset; i <= j; i++) 
					values[ctr++] = eps[i] + eps[j] - eps[a] - eps[b]; 
	Dijab.write(sz, indices, values); 
	
	CTF::Function<> fctr(&dividef);
	T.abij->contract(1.0, *V.iajb, "iajb", Dijab, "ijab", 0.0, "abij", fctr);
	CTF::Tensor<> Aijab(4, dlens, nsns, dw); 
	Aijab["ijab"] = 2.0 * V["iajb"] - V["ibja"]; 
	energy = T["abij"] * Aijab["ijab"];
	
	T.ai->read_local(&sz, &indices, &values); 
	for (int i = 0; i < sz; i++) values[i] = 0.0;
	T.ai->write(sz, indices, values); 
	
	Matrix F = eps.asDiagonal(); 

	V.ab->read_local(&sz, &indices, &values);
	ctr = 0;  
	for (int a = 0; a < nvirt; a++) { 
		for (int b = 0; b <= a; b++) {
			values[ctr] = F(b+offnocc, a+offnocc); 
			ctr++; 
		}
	}
	V.ab->write(sz, indices, values); 

	V.ai->read_local(&sz, &indices, &values);
	ctr = 0;  
	for (int i = 0; i < nocc; i++) { 
		for (int a = 0; a < nvirt; a++) {
			values[ctr] = F(a+offnocc, i+offset); 
			ctr++; 
		}
	}
	V.ai->write(sz, indices, values); 
	V["ia"] = V["ai"]; 

	V.ij->read_local(&sz, &indices, &values);
	ctr = 0;  
	for (int i = 0; i < nocc; i++) { 
		for (int j = 0; j <= i; j++) {
			values[ctr] = F(j+offset, i+offset); 
			ctr++; 
		}
	}
	V.ij->write(sz, indices, values);

}

void MP2::tensormp2(bool print) {
	
	int nvirt = N-nocc; 
	int offset = 0;

	if(!focker.getMolecule()->control->get_option<bool>("withcore")) {
		offset = nocc - focker.getMolecule()->getNValence() / 2; 
		nocc -= offset; 
		N = nocc + nvirt; 
	}

	int offnocc = nocc + offset; 
	int offN = N + offset; 

	if (print) {
		focker.getMolecule()->control->log.print("No. of occ. orbitals: " + std::to_string(nocc));
		focker.getMolecule()->control->log.print("No. of virt. orbitals: " + std::to_string(nvirt)); 
	}
	
	int ao_lens[4] = {offN, offN, offN, offN};
	int mo_lens[4] = {nocc, nvirt, nocc, nvirt};
	int t1_lens[4] = {offN, offN, offN, nvirt}; 
	int t2_lens[4] = {offN, offN, nocc, nvirt};
	int t3_lens[4] = {offN, nvirt, nocc, nvirt};
	int sysy[4] = {SY, NS, SY, NS}; 
	int nssy[4] = {NS, NS, SY, NS};
	int nsns[4] = {NS, NS, NS, NS};
	int syns[4] = {SY, NS, NS, NS}; 
	CTF::Tensor<> ao_integrals(4, ao_lens, sysy, dw); 
	int64_t sz, *i1, *i2, *i3, *i4;
	double *v1, *v2, *v3, *v4;
	int ctr;
	
	IntegralEngine& aoInts = focker.getIntegrals();
	
	ao_integrals.read_local(&sz, &i1, &v1);
	ctr = 0;
	for (int a = 0; a < offN; a++) 
		for (int b = 0; b <= a; b++)
			for (int c = 0; c < offN; c++)
				for (int d = 0; d <= c; d++)
					v1[ctr++] = aoInts.getERI(d, c, b, a); 
	ao_integrals.write(sz, i1, v1);
	
	delete v1; 
	delete i1; 
	
	CTF::Matrix<> cp_occ(offN, nocc, NS, dw); 
	CTF::Matrix<> cp_virt(offN, nvirt, NS, dw); 
	Matrix& CP = focker.getCP(); 
	
	ctr = 0; 
	cp_occ.read_local(&sz, &i2, &v2);
	for (int i = offset; i < offnocc; i++)
		for (int mu = 0; mu < offN; mu++)
			v2[ctr++] = CP(mu, i); 
	cp_occ.write(sz, i2, v2); 
	
	delete i2;
	delete v2; 
	
	ctr = 0; 
	cp_virt.read_local(&sz, &i3, &v3);
	for (int a = offnocc; a < offN; a++)
		for (int mu = 0; mu < offN; mu++)
			v3[ctr++] = CP(mu, a); 
	cp_virt.write(sz, i3, v3); 
	
	delete i3; 
	delete v3; 
	
	CTF::Tensor<> temp1(4, t1_lens, syns, dw);
	CTF::Tensor<> temp2(4, t2_lens, syns, dw); 
	CTF::Tensor<> temp3(4, t3_lens, nsns, dw);
	CTF::Tensor<> Viajb(4, mo_lens, nsns, dw); 
	
	temp1["uvwb"] = cp_virt["xb"] * ao_integrals["uvwx"]; 
	temp2["uvjb"] = cp_occ["wj"] * temp1["uvwb"]; 
	temp3["uajb"] = cp_virt["va"] * temp2["uvjb"];
	Viajb["iajb"] = cp_occ["ui"] * temp3["uajb"]; 
	
	int dlens[4] = {nocc, nocc, nvirt, nvirt}; 
	CTF::Tensor<> Dijab(4, dlens, sysy, dw); 
	
	ctr = 0;
	Vector& eps = focker.getEps(); 
	Dijab.read_local(&sz, &i4, &v4); 
	for (int b = offnocc; b < offN; b++)
		for (int a = offnocc; a <= b; a++)
			for (int j = offset; j < offnocc; j++)
				for (int i = offset; i <= j; i++) 
					v4[ctr++] = eps[i] + eps[j] - eps[a] - eps[b]; 
	Dijab.write(sz, i4, v4); 
	
	delete i4;
	delete v4; 
	
	CTF::Function<> fctr(&dividef);
	CTF::Tensor<> Tijab(4, dlens, nsns, dw);
	Tijab.contract(1.0, Viajb, "iajb", Dijab, "ijab", 0.0, "ijab", fctr);
	CTF::Tensor<> Aijab(4, dlens, nsns, dw); 
	Aijab["ijab"] = 2.0 * Viajb["iajb"] - Viajb["ibja"]; 
	energy = Tijab["ijab"] * Aijab["ijab"];
	 
}

// Determine the MP2 energy
void MP2::calculateEnergy()
{
	energy = 0.0; 
	
	int nvirt = N-nocc; 
	int offset = 0;

	if(!focker.getMolecule()->control->get_option<bool>("withcore")) {
		offset = nocc - focker.getMolecule()->getNValence() / 2; 
		nocc -= offset; 
		N = nocc + nvirt; 
	}

	int offnocc = nocc + offset; 
	int offN = N + offset; 
	
	Vector& eps = focker.getEps(); 
	double resolvent, temp; 
	for (int i = offset; i < offnocc; i++)
		for (int j = offset; j < offnocc; j++)
			for (int a = offnocc; a < offN; a++)
				for (int b = offnocc; b < offN; b++) {
					temp = moInts(i, a, j, b) * (2.0 * moInts(i, a, j, b) - moInts(i, b, j, a)); 
					resolvent = eps[i] + eps[j] - eps[a] - eps[b]; 
					energy += temp / resolvent; 
				}
	
}

void MP2::dfmp2(bool print) {
	
	Logger& log = focker.getMolecule()->control->log; 
	
	IntegralEngine& integrals = focker.getIntegrals(); 
	std::vector<libint2::Shell>& obs = focker.getMolecule()->getBasis().getIntShells();
	std::vector<libint2::Shell>& auxbs = focker.getMolecule()->getBasis().getRIShells();
	
	// Calculate (P|Q) and (P|mn) eris
	
	if (print) log.print("Calculating two-index eris...");
	Matrix JPQ = integrals.compute_eris_2index(auxbs); 
	if (print) log.localTime();  
	
	if (print) log.print("Calculating three-index eris..."); 
	Matrix KmnP;
	integrals.compute_eris_3index(obs, auxbs, KmnP);
	if (print) log.localTime(); 

	// Form inverse J-metric
	JPQ = JPQ.inverse(); 
	
	int nobs = integrals.nbasis(obs);
	int nabs = JPQ.rows(); 
	
	int nvirt = N-nocc; 
	int offset = 0;

	if(!focker.getMolecule()->control->get_option<bool>("withcore")) {
		offset = nocc - focker.getMolecule()->getNValence() / 2; 
		nocc -= offset; 
		N = nocc + nvirt; 
	}
	int offN = offset + N; 
	int offnocc = offset + nocc; 
	
	if (print) {
		focker.getMolecule()->control->log.print("No. of occ. orbitals: " + std::to_string(nocc));
		focker.getMolecule()->control->log.print("No. of virt. orbitals: " + std::to_string(nvirt)); 
	}
	
	if (print) log.print("Transforming integrals to MO basis...");
	Matrix cp_occ = focker.getCP().block(0, offset, offN, nocc); 
	Matrix cp_virt = focker.getCP().block(0, offnocc, offN, nvirt); 
	
	Matrix BmnP = Matrix::Zero(nocc*nobs, nabs); 
	int ix; 
	for (int i = 0; i < nocc; i++) {
		for (int nu = 0; nu < nobs; nu++) {
			for (int P = 0; P < nabs; P++) {
				
				ix = i*nobs+nu; 
				for (int mu = 0; mu < nobs; mu++) {
					int X = std::max(mu, nu);
					int Y = std::min(mu, nu); 
					BmnP(ix, P) += cp_occ(mu, i) * KmnP((X*(X+1))/2+Y, P); 
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
					KmnP(ix, P) += cp_virt(nu, a) * BmnP(i*nobs+nu, P); 
				
			}
		}
	}
	if (print) log.localTime();
	
	BmnP = JPQ * KmnP.transpose(); 
	
	Vector& eps = focker.getEps(); 
	double resolvent, viajb, vibja; 
	int ia, jb, ib, ja; 
	for (int i = 0; i < nocc; i++) {
		for (int a = 0; a < nvirt; a++) {
			ia = i * nvirt + a;
			for (int j = 0; j < nocc; j++) {
				ja = j*nvirt + a; 
				for (int b = 0; b < nvirt; b++) {
					ib = i*nvirt + b;
					jb = j*nvirt + b;
					resolvent = eps[i+offset] + eps[j+offset] - eps[a+offnocc] - eps[b+offnocc]; 
					
					viajb = vibja = 0.0; 
					for (int Q = 0; Q < nabs; Q++) {
						viajb += KmnP(ia, Q) * BmnP(Q, jb);
						vibja += KmnP(ib, Q) * BmnP(Q, ja); 
					}
					
					energy += viajb*(2.0*viajb - vibja) / resolvent;
					
				}
			}
		}
	}
	
}

					
