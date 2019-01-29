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

#ifndef FOCKHEADERDEF
#define FOCKHEADERDEF

// Includes
#include "integrals.hpp"
#include "molecule.hpp"
#include "eigen_wrapper.hpp"
#include <vector>
#include <string>

// Forward declarations
class Atom;
class Command; 

/*! \struct Domain
	\brief  Orbital domains for the local exchange procedure.
 */
struct Domain {
	/// The basis function indices where each part of the Domain starts.
	std::vector<int> starts;  
	/// The number of basis functions in each part of the Domain.
	std::vector<int> sizes;
	/// The indices of the fragments on which each part of the Domain is centred.
	std::vector<int> centres;
	
	int totalsize; ///< The total number of basis functions in the Domain. 
	Matrix G; ///< The inverse density-fitting metric for this Domain.
	
	Domain() {} ///< Creates an empty Domain. 
	
	/// Copy constructor
	Domain(const Domain& other) {
		starts = other.starts;
		sizes = other.sizes;
		centres = other.centres;
		totalsize = other.totalsize; 
	}
	
	/// Determines the total size (number of basis functions) of the Domain.
	void sumsizes() {
		totalsize = 0; 
		for (auto s : sizes)
			totalsize += s;
	}
	
	/// Prints details of the Domain to the primary output (for debugging purposes only). 
	void print() {
		std::cout << "STARTS:" << std::endl; 
		for (auto s : starts)
			std::cout << s << ", "; 
		std::cout << std::endl;
		
		std::cout << "SIZES:" << std::endl;
		for (auto s : sizes)
			std::cout << s << ", "; 
		std::cout << std::endl;
		
		std::cout << "CENTRES:" << std::endl;
		for (auto c : centres)
			std::cout << c << ", ";
		std::cout << std::endl << std::endl; 
	}
};

/*! \struct FragmentInfo
 	\brief  Stores all information about a molecular fragment necessary for local density-fitting. 
 */
struct FragmentInfo {
	int occ; 			///< The number of occupied orbitals on this fragment.
	int nbfs; 			///< The number of orbital basis functions on this fragment. 
	int naux;			///< The number of auxiliary basis functions on this fragment.
	int start;			///< The index of the first orbital basis function on this fragment.
	int auxstart;		///< The index of the first auxiliary basis function on this fragment.
	int nshells;		///< The total number of shells in the orbital basis on this fragment.
	int ndfshells;		///< The total number of shells in the auxiliary basis on this fragment.
	
	double radius;		///< The radial extent (in a.u.) of the Basis on this fragment.  
	double mo_thresh;	///< The threshold for including a fragment in the LMO domain. 
	double fit_thresh; 	///< The threshold for including a fragment in the fitting domain. 
	double r_thresh;	///< Threshold distance for connecting two fragments in the fitting domain. 
	
	Vector com; 	 	///< The centre of mass xyz coordinates (in a.u.) of the fragment.
	
	/// Creates an empty FragmentInfo, with default threshold values. 
	FragmentInfo() : occ(0), nbfs(0), naux(0), start(0), auxstart(0), nshells(0), ndfshells(0), radius(0.0), mo_thresh(1e-6), fit_thresh(0.05), r_thresh(15.0) {}
	
	/// Copy constructor
	FragmentInfo(const FragmentInfo& other) {
		occ = other.occ;
		nbfs = other.nbfs;
		naux = other.naux;
		start = other.start;
		auxstart = other.auxstart;
		nshells = other.nshells; 
		ndfshells = other.ndfshells; 
		radius = other.radius;
		mo_thresh = other.mo_thresh;
		fit_thresh = other.fit_thresh;
		r_thresh = other.r_thresh; 
		com = other.com; 
	}
};

/*! \class Fock
	\brief Encapsulates the data and methods shared between all methods that use a Fock matrix.

	All methods currently in GAMMA use Fock matrices in some way. The Fock matrix itself, along with
	references to various integrals, the core Hamiltonian, and orbital domains are stored here, as
	they are the necessary ingredients for building a Fock matrix. 

	The methods are for building Fock matrices in various ways - using stored two-electron integrals 
	(from memory or file), using an integral-direct procedure, density fitting, or local density-fitting.
	This includes gradient and Hessian Fock contributions, and generalised Fock builds necessary for 
	generating superposition of atomic density guesses. 
 */ 
class Fock
{
protected:
  Matrix hcore;		///< The core Hamiltonian matrix. 												
  Matrix jkints;	///< The JK antisymmetrised two-electron integral matrix.			
  Matrix jints;		///< The Coulomb contribution to the Fock matrix.
  Matrix kints;		///< The exchange contribution to the Fock matrix.
  Matrix orthog;	///< The Lowdin symmetric orthogonaliser.
  Matrix fockm;		///< The Fock matrix in the MO basis.
  Matrix focka;		///< The Fock matrix in the AO basis.
  Matrix fockinc; 	///< The incremental Fock matrix.
  Matrix CP;		///< Matrix of MO coefficients (columns are eigenvectors).
  Matrix forces;	///< Nuclear forces (N x 3)
  Matrix hessian;	///< Nuclear Hessian matrix
  Matrix dens;		///< The 1-particle reduced density matrix in the AO basis. 
  Matrix dens_diff; ///< The change in density from the last iteration.
  
  Matrix xyK;				///< Three-centre density fitting ERIs (xy | K)  
  DFBlocks blocked_xyK; 	///< Fragment-blocked three-centre density fitting ERIs 
  SparseMatrix Linv; 		///< The Cholesky decomposition of the inverse density fitting metric.  
  SparseMatrix Linv2; 		///< LL^t of the Cholesky decomposition in Linv. 
   
  Vector eps; 		///< Vector of eigenvalues of the Fock matrix, in ascending order.
  
  std::vector<Matrix> focks; ///< Previous Fock matrices, needed for DIIS extrapolation.
  
  std::vector<Domain> lmo_domains;	///< 
  std::vector<Domain> ao_domains;
  std::vector<Domain> fit_domains;
  
  std::vector<std::vector<int> > Kblocks; 
  
  IntegralEngine& integrals;
  SharedMolecule molecule;
  bool direct, twoints, fromfile, diis, density_fitted;
  double precision; 
  int nbfs, iter, MAX, nocc;
  std::string guess; 
public:
	
    bool reset_incremental, started_incremental; 
    double next_reset, rms_error, incremental_threshold; 
	int last_reset; 
  
  Fock(Command& cmd, IntegralEngine& ints, SharedMolecule m);
  Fock(const Fock& other);
  IntegralEngine& getIntegrals() { return integrals; }
  SharedMolecule getMolecule() { return molecule; }
  Matrix& getHCore() { return hcore; }
  Matrix& getFockAO() { return focka; }
  Matrix& getFockMO() { return fockm; }
  Matrix& getOrthog() { return orthog; }
  Matrix& getCP() { return CP; }
  Vector& getEps() { return eps; }
  Matrix& getJK() { return jkints; }
  Matrix& getJ() { return jints; } 
  Matrix& getK() { return kints; }
  Matrix& getXYK() { return xyK; }
  DFBlocks& getBlockedXYK() { return blocked_xyK; }
  SparseMatrix& getLinv() { return Linv; }
  SparseMatrix& getLinv2() { return Linv2; }
  Matrix& getDens() { return dens; }
  Matrix& getForces() { return forces; }
  Matrix& getHessian() { return hessian; }
  virtual Matrix& getS() { return integrals.getOverlap(); }
  Matrix getHCore() const { return hcore; }
  Matrix getFockAO() const { return focka; }
  Matrix getFockMO() const { return fockm; }
  Matrix getOrthog() const { return orthog; }
  Matrix getCP() const { return CP; }
  Vector getEps() const { return eps; }
  Matrix getJK() const { return jkints; }
  Matrix getJ() const { return jints; } 
  Matrix getK() const { return kints; }
  Matrix getForces() const { return forces; }
  Matrix getHessian() const { return hessian; }
  std::vector<Domain>& getLMODomains() { return lmo_domains; }
  std::vector<Domain>& getAODomains() { return ao_domains; }
  std::vector<Domain>& getFitDomains() { return fit_domains; }
  virtual Matrix getS() const { return integrals.getOverlap(); }
  virtual Matrix getDens() const { return dens; }
  void setDIIS(bool d) { diis = d; } 
  void formHCore();
  void formOrthog();
  virtual void transform(bool first = false);
  virtual void diagonalise();
  virtual void makeJK(Matrix& P, double multiplier = 1.0);
  void formJK(Matrix& P, double multiplier = 1.0);
  void formJKlocal(Matrix& Cocc, const Matrix& sigmainv, Matrix& Pt, std::vector<FragmentInfo>& finfo, bool from_file = false); 
  void formJKdirect(const Matrix& Schwarz, Matrix& P1, double multiplier = 1.0);
  void formJKdf(Matrix& Cocc, double multiplier = 1.0); 
  virtual void formJKfile();
  virtual void makeFock(); 
  virtual void makeFock(Matrix& P, double multiplier = 1.0 );
  virtual void makeDens();
  virtual void average(Vector &w);
  void simpleAverage(Matrix& D0, double weight = 0.5);
  
  virtual void addDiis(); 
  virtual void clearDiis(); 
  
  virtual void compute_forces(const std::vector<Atom> &atoms, int nocc); 
  virtual void compute_hessian(const std::vector<Atom> &atoms, int nocc);
  virtual void compute_hessian_numeric(const std::vector<Atom> &atoms, int nocc, Command &cmd);
  
  Matrix compute_2body_fock_general(
  	const std::vector<libint2::Shell>& obs, const Matrix& D, const std::vector<libint2::Shell>& D_bs,
  bool D_is_sheldiagonal = false,  // set D_is_shelldiagonal if doing SOAD
  double precision = std::numeric_limits<double>::epsilon()  // discard contributions smaller than this
  		); 

  void compute_2body_fock_df(const Matrix& Cocc, Matrix& j, Matrix& k);
  void compute_2body_fock_df_local(Matrix& Cocc, const Matrix& sigmainv, Matrix& Pt, std::vector<FragmentInfo>& finfo, Matrix& j, Matrix& k,
  		std::vector<Domain>& lmod, std::vector<Domain>& aod, std::vector<Domain>& fitd); 
  void compute_2body_fock_df_local_file(Matrix& Cocc, const Matrix& sigmainv, Matrix& Pt, std::vector<FragmentInfo>& finfo, Matrix& j, Matrix& k,
  	    std::vector<Domain>& lmod, std::vector<Domain>& aod, std::vector<Domain>& fitd);
  void build_domains(Matrix& Cocc, Matrix& V, std::vector<FragmentInfo>& finfo, 
  	std::vector<Domain>& lmod, std::vector<Domain>& aod, std::vector<Domain>& fitd);  
  void build_blocked_eris(std::vector<FragmentInfo>& finfo, Matrix& JPQ, Matrix& Pt); 
  
  void compute_soad_guess(); 
  
  virtual Vector compute_xgrad(double fx, Matrix& xhessian, std::vector<int>& activex, Command& cmd); 
  
  template <unsigned deriv_order>
  std::vector<EMatrix> compute_2body_fock_deriv(const std::vector<Atom> &atoms, const EMatrix& D);
};

class UnrestrictedFock : public Fock
{
protected:
	Matrix fock_alpha_ao, fock_beta_ao; 
	Matrix fock_alpha_mo, fock_beta_mo;
	Matrix dens_alpha, dens_beta; 
	Matrix kints_alpha, kints_beta, jints_alpha, jints_beta; 
	Matrix CP_alpha, CP_beta;
	Vector eps_alpha, eps_beta;
	std::vector<Domain> lmo_alpha, ao_alpha, fit_alpha, lmo_beta, ao_beta, fit_beta; 
	std::vector<Matrix> alpha_focks, beta_focks;
	int nalpha, nbeta; 
public:
    UnrestrictedFock(Command& cmd, IntegralEngine& ints, SharedMolecule m);
    UnrestrictedFock(const UnrestrictedFock& other);
	
	Matrix& getFockAlphaAO() { return fock_alpha_ao; }
	Matrix& getFockBetaAO() { return fock_beta_ao; }
	Matrix& getFockAlphaMO() { return fock_alpha_ao; }
	Matrix& getFockBetaMO() { return fock_beta_ao; }
	Matrix& getDensAlpha() { return dens_alpha; }
	Matrix& getDensBeta() { return dens_beta; }
	Matrix& getKAlpha() { return kints_alpha; }
	Matrix& getKBeta() { return kints_beta; }
	Matrix& getJAlpha() { return jints_alpha; }
	Matrix& getJBeta() { return jints_beta; }
	Matrix& getCPAlpha() { return CP_alpha; }
	Matrix& getCPBeta() { return CP_beta; }
	Vector& getEpsAlpha() { return eps_alpha; }
	Vector& getEpsBeta() { return eps_beta; }
	
    virtual void transform(bool first = false);
    virtual void diagonalise();
    virtual void makeJK();
   	void formJK(Matrix& P1, Matrix& P2, double multiplier = 1.0);
    void formJKdirect(const Matrix& Schwarz, Matrix& P1, Matrix& P2, double multiplier = 1.0);
	void formJKdf(Matrix& ca, Matrix& cb, double multiplier = 1.0);
	void formJKlocal(Matrix& ca, Matrix& cb, const Matrix& sa, const Matrix& sb, Matrix& Pa, Matrix& Pb,
		 std::vector<FragmentInfo>& finfoa, std::vector<FragmentInfo>& finfob);   
    virtual void formJKfile();
    virtual void makeFock();
    virtual void makeDens();
    virtual void average(Vector &w);

	virtual void clearDiis(); 
};

class FockFragment : public Fock 
{
private:
	int start, end;
public:
	Matrix Sxx; 
	FockFragment(Command& cmd, IntegralEngine& ints, SharedMolecule m, int start, int end);
	FockFragment(const FockFragment& other);
	Vector buildFock(Matrix& qfq, Matrix& qfp, Matrix& pfp, bool alpha = false); 
	virtual void gensolve();
};


class UnrestrictedFockFragment : public UnrestrictedFock
{
private: 
	int start, end;
public:
	Matrix Sxx;
	UnrestrictedFockFragment(Command& cmd, IntegralEngine& ints, SharedMolecule m, int start, int end);
	UnrestrictedFockFragment(const UnrestrictedFockFragment& other);
	Vector buildFock(Matrix& qfq, Matrix& qfp, Matrix& pfp, bool alpha); 
	virtual void gensolve();
};

#endif
