/*
*
*   PURPOSE: To implement class Fock, a class containing the data and routines
*            needed for Fock space based methods.
*
*   DATE          AUTHOR              CHANGES
*  ===========================================================================
*  14/09/15       Robert Shaw         Original code.
* 
*/

//#ifndef EIGEN_USE_MKL_ALL
//#define EIGEN_USE_MKL_ALL
//#endif

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
#include "sto3g_atomic_density.hpp"
#include "basis.hpp"
#include "basisreader.hpp"
#include <map> 

// Constructor
Fock::Fock(Command& cmd, IntegralEngine& ints, SharedMolecule m) : integrals(ints), molecule(m)
{
	
	// Make the core hamiltonian matrix
	formHCore();

	// Form the orthogonalising matrix and get initial density guess
	try {
		formOrthog();
	} catch (Error e) {
		molecule->control->log.error(e);
	}

	// Retrieve whether this calculation should be done direct, or whether
	// the twoints matrix has been formed, or if the 2e integrals need to be
	// read from file.
	direct = molecule->control->get_option<bool>("direct");
	density_fitted = cmd.get_option<bool>("df"); 
	diis = cmd.get_option<bool>("diis");
	iter = 0;
	nocc = molecule->getNel() / 2; 
	MAX = cmd.get_option<int>("maxdiis");
	precision = cmd.get_option<double>("precision");
	guess = cmd.get_option<std::string>("guess"); 
	
	reset_incremental = started_incremental = false; 
	next_reset = 0.0; 
	rms_error = 1.0; 
	incremental_threshold = 1e-5; 
	
	twoints = false;
	if (!direct){
		Vector ests = integrals.getEstimates();
		if (ests[3] < molecule->control->get_option<double>("memory")) 
			twoints = true;
	}

	fromfile = false;
	if (!twoints && !direct)
		fromfile = true;

}

Fock::Fock(const Fock& other) : integrals(other.integrals), molecule(other.molecule) {
	hcore = other.hcore;
	jkints = other.jkints;
	jints = other.jints;
	kints = other.kints;
	xyK = other.xyK; 
	orthog = other.orthog;
	fockm = other.fockm;
	focka = other.focka;
	fockinc = other.fockinc; 
	CP = other.CP;
	forces = other.forces;
	hessian = other.hessian;
	eps = other.eps;
	focks = other.focks;
	dens = other.dens;
	dens_diff = other.dens_diff; 
	direct = other.direct;
	twoints = other.twoints;
	fromfile = other.fromfile;
	density_fitted = other.density_fitted; 
	guess = other.guess; 
	diis = other.diis;
	nbfs = other.nbfs;
	nocc = other.nocc;
	iter = other.iter;
	MAX = other.MAX; 
	precision = other.precision; 
	
    reset_incremental = other.reset_incremental;
	started_incremental = other.started_incremental; 
   	next_reset = other.next_reset; 
	rms_error = other.rms_error;
	incremental_threshold = other.incremental_threshold; 
	last_reset = other.last_reset; 
	
}

// Form the core hamiltonian matrix
void Fock::formHCore()
{
	// The core Hamiltonian matrix is defined to be
	// the sum of the kinetic and nuclear attraction
	// matrices
	hcore = integrals.getKinetic() + integrals.getNucAttract();
	nbfs = hcore.rows();

}

void Fock::formOrthog()
{

	Matrix& S = integrals.getOverlap();
	
	// Diagonalise the overlap matrix into lambda and U,
	// so that U(T)SU = lambda
	// We can now form S^(-1/2) - the orthogonalising matrix
	EigenSolver es(S);
	Matrix U = es.eigenvectors();
	Vector lambda = es.eigenvalues(); 
	
	molecule->control->log.print("Lowest eigenvalue of the metric: " + std::to_string(lambda[0])); 
  
	orthog = Matrix::Zero(nbfs, nbfs);
	for (int i = 0; i < nbfs; i++)
		orthog(i, i) = 1.0/(std::sqrt(lambda(i)));
  
	// S^-1/2  = U(lambda^-1/2)U(T)
	orthog = U * orthog * U.transpose();
}
  
void Fock::average(Vector &w) {
	if (diis && iter > 2) {
		// Average the fock matrices according to the weights
		focka = Matrix::Zero(nbfs, nbfs);
		int offset = focks.size() - w.size();
		for (int i = offset; i < focks.size(); i++) {
			focka = focka + w[i-offset]*focks[i]; 
		} 
	}
}

// Transform the AO fock matrix to the MO basis 
void Fock::transform(bool first)
{
	if (first) { 
		if (guess == "hcore") {
			// Form the core Fock matrix as (S^-1/2)(T)H(S^-1/2)
			fockm = (orthog.transpose()) * ( hcore * orthog);
			fockinc = hcore; 
		} else {
			// use superposition of atomic densities
			compute_soad_guess();
			fockm = (orthog.transpose()) * ( focka * orthog ); 
		}
	} else {
		// Form the orthogonalised fock matrix
		fockm = orthog.transpose() * (focka * orthog);
	}
}

void Fock::compute_soad_guess() {
	std::vector<Atom> atoms; 
	for (int i = 0; i < molecule->getNAtoms(); i++) atoms.push_back(molecule->getAtom(i));
		
	std::map<int, std::string> names;
	std::vector<libint2::Atom> q(atoms.size()); 
	size_t nao = 0;
	int i = 0;
	for(const auto& atom : atoms) {
		q[i].atomic_number = atom.getCharge(); 
		q[i].x = atom.getX();
		q[i].y = atom.getY();
		q[i].z = atom.getZ(); 
		
		nao += sto3g_num_ao(atom.getCharge()); 
		names[atom.getCharge()] = "STO-3G"; 
	}

	// compute the minimal basis density
	Matrix D = Matrix::Zero(nao, nao);
	size_t ao_offset = 0;  // first AO of this atom
	for (const auto& atom : atoms) {
		const auto Z = atom.getCharge();
		const auto& occvec = sto3g_ao_occupation_vector(Z);
		for(const auto& occ: occvec) {
			D(ao_offset, ao_offset) = occ;
			++ao_offset;
		}
	}
	D *= 0.5; 

	Basis minbasis; 
	BasisReader breader(names); 
	for (int i = 0; i < atoms.size(); i++) {
		Atom &a = atoms[i];
		double pos[3] = { a.getX(), a.getY(), a.getZ() };
		breader.readShellBasis(minbasis, a.getCharge(), pos, i, false);
	}
	
	std::vector<libint2::Shell>& minbs = minbasis.getIntShells(); 
	
	// if basis != minimal basis, map non-representable SOAD guess
	// into the AO basis
	// by diagonalizing a Fock matrix
	focka = hcore;
	focka += compute_2body_fock_general(molecule->getBasis().getIntShells(),
	D, minbs, true /* SOAD_D_is_shelldiagonal */,
	std::numeric_limits<double>::epsilon()  // this is cheap, no reason
		// to be cheaper
		);
	fockinc = focka; 
}	

// Diagonalise the MO fock matrix to get CP and eps
void Fock::diagonalise() 
{
	EigenSolver es(fockm);
	CP = orthog * es.eigenvectors();
	eps = es.eigenvalues();
}

// Construct the density matrix from CP, 
// for nocc number of occupied orbitals
void Fock::makeDens()
{
	// Form the density matrix
	if(dens.rows() != 0) dens_diff = -dens; 
	dens = 2.0 * CP.block(0, 0, nbfs, nocc) * CP.transpose().block(0, 0, nocc, nbfs); 
	if(dens_diff.rows() == 0) dens_diff = dens;
	else dens_diff += dens;  
}

// Make the JK matrix, depending on how two electron integrals are stored/needed

void Fock::simpleAverage(Matrix& D0, double weight)
{
	dens = weight*dens + (1.0-weight)*D0;
}


void Fock::clearDiis() {
	iter = 0;
	focks.clear(); 
}

FockFragment::FockFragment(Command& cmd, IntegralEngine& ints, SharedMolecule m, int _start, int _end) : Fock(cmd, ints, m), start(_start), end(_end) { 
	Sxx = ints.getOverlap();
	incremental_threshold = 0.0; 
}

FockFragment::FockFragment(const FockFragment& other) : Fock(other) {
	Sxx = other.Sxx;
	start = other.start;
	end = other.end;
}

void FockFragment::gensolve()
{
	GeneralizedEigenSolver es(focka, Sxx); 
	CP = es.eigenvectors();
	eps = es.eigenvalues();
	iter++;
}

UnrestrictedFock::UnrestrictedFock(Command& cmd, IntegralEngine& ints, SharedMolecule m) : Fock(cmd, ints, m)
{
	nalpha = molecule->nalpha();
	nbeta = molecule->nbeta(); 
}

UnrestrictedFock::UnrestrictedFock(const UnrestrictedFock& other) : Fock(other)
{
	fock_alpha_ao = other.fock_alpha_ao;
	fock_beta_ao = other.fock_alpha_ao;
	fock_alpha_mo = other.fock_alpha_ao;
	fock_beta_mo = other.fock_alpha_ao;
	dens_alpha = other.dens_alpha;
	dens_beta = other.dens_beta;
	CP_alpha = other.CP_alpha;
	eps_alpha = other.eps_alpha;
	CP_beta = other.CP_beta;
	eps_beta = other.eps_beta;
	kints_alpha = other.kints_alpha;
	kints_beta = other.kints_beta;
	jints_alpha = other.kints_alpha;
	jints_beta = other.kints_beta;
	nalpha = other.nalpha;
	nbeta = other.nbeta; 
}

// Transform the AO fock matrix to the MO basis 
void UnrestrictedFock::transform(bool first)
{
	if (first) { 
		if (guess == "hcore") {
			// Form the core Fock matrix as (S^-1/2)(T)H(S^-1/2)
			fock_alpha_mo = (orthog.transpose()) * ( hcore * orthog);
			fock_beta_mo = fock_alpha_mo; 
		} else {
			compute_soad_guess(); 
			fock_alpha_mo = orthog.transpose() * ( focka * orthog ); 
			fock_beta_mo = fock_alpha_mo; 
		}
	} else {
		// Form the orthogonalised fock matrix
		fock_alpha_mo = orthog.transpose() * (fock_alpha_ao * orthog);
		fock_beta_mo = orthog.transpose() * (fock_beta_ao * orthog);
	}
}

// Diagonalise the MO fock matrix to get CP and eps
void UnrestrictedFock::diagonalise() 
{
	EigenSolver es_alpha(fock_alpha_mo);
	EigenSolver es_beta(fock_beta_mo);
	CP_alpha = orthog * es_alpha.eigenvectors();
	CP_beta = orthog * es_beta.eigenvectors(); 
	eps_alpha = es_alpha.eigenvalues();
	eps_beta = es_beta.eigenvalues();
}

// Construct the density matrix from CP, 
// for nocc number of occupied orbitals
void UnrestrictedFock::makeDens()
{
	// Form the density matrix
	dens_alpha =  CP_alpha.block(0, 0, nbfs, nalpha) * CP_alpha.transpose().block(0, 0, nalpha, nbfs);
	dens_beta = CP_beta.block(0, 0, nbfs, nbeta) * CP_beta.transpose().block(0, 0, nbeta, nbfs); 
}

void UnrestrictedFock::average(Vector &w) {
	if (diis && iter > 2) {
		// Average the fock matrices according to the weights
		fock_alpha_ao = Matrix::Zero(nbfs, nbfs);
		fock_beta_ao = Matrix::Zero(nbfs, nbfs);
		int offset = alpha_focks.size() - w.size();
		for (int i = offset; i < alpha_focks.size(); i++) {
			fock_alpha_ao = fock_alpha_ao + w[i-offset]*alpha_focks[i];
			fock_beta_ao = fock_beta_ao + w[i-offset]*beta_focks[i]; 
		} 
	}
}

void UnrestrictedFock::formJKfile() {
	std::cerr << "File-read not implemented yet!" << std::endl;  
}

void UnrestrictedFock::clearDiis() {
	iter = 0;
	alpha_focks.clear();
	beta_focks.clear();
}

UnrestrictedFockFragment::UnrestrictedFockFragment(Command& cmd, IntegralEngine& ints, SharedMolecule m, int _start, int _end) : UnrestrictedFock(cmd, ints, m), start(_start), end(_end) { 
	Sxx = ints.getOverlap();
}

UnrestrictedFockFragment::UnrestrictedFockFragment(const UnrestrictedFockFragment& other) : UnrestrictedFock(other) {
	Sxx = other.Sxx;
	start = other.start;
	end = other.end;
}

void UnrestrictedFockFragment::gensolve()
{
	GeneralizedEigenSolver es_alpha(fock_alpha_ao, Sxx); 
	CP_alpha = es_alpha.eigenvectors();
	eps_alpha = es_alpha.eigenvalues(); 
	
	GeneralizedEigenSolver es_beta(fock_beta_ao, Sxx);
	CP_beta = es_beta.eigenvectors();
	eps_beta = es_beta.eigenvalues(); 
	
	iter++;
}

