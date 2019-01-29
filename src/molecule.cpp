/*
*
*     Implements molecule.hpp, defining class Molecule
*
*     DATE         AUTHOR           CHANGES  
*     ==============================================================
*     26/08/15     Robert Shaw      Original code.
*
*/

#include <iostream>
#include "molecule.hpp"
#include "error.hpp"
#include "basisreader.hpp"
#include "ioutil.hpp"
#include <cmath>
#include <libint2.hpp>
#include "ProgramController.hpp"

Atom Molecule::parseGeomLine(std::string line) {
	std::string delimiter = ","; // Define the separation delimiter
	int q;
	double m;
	Vector coords(3);
	double multiplier = (angstrom ? Logger::TOBOHR : 1.0); // Convert to bohr from angstrom?
		
	size_t position = line.find(delimiter); // Find first delimiter
	std::string token = line.substr(0, position); // Tokenise

	// Get rid of whitespace, if any
	token.erase(std::remove(token.begin(), token.end(), ' '), token.end());
	token.erase(std::remove(token.begin(), token.end(), '\t'), token.end());
		
	// Get the atom type and mass
	q = getAtomCharge(token);
	m = getAtomMass(q);

	// Move on to next token
	line.erase(0, position+delimiter.length());
	position = line.find(delimiter);

	if (position != std::string::npos){
		token = line.substr(0, position); 

		// This should now be the x coord
		token.erase(std::remove(token.begin(), token.end(), ' '), token.end());
		token.erase(std::remove(token.begin(), token.end(), '\t'), token.end());
		coords[0] = multiplier*std::stod(token);  // Convert it to double

		// Repeat for y and z
		line.erase(0, position+delimiter.length());
		position = line.find(delimiter);
		if (position != std::string::npos) {
			token = line.substr(0, position);
			token.erase(std::remove(token.begin(), token.end(), ' '), token.end());
			token.erase(std::remove(token.begin(), token.end(), '\t'), token.end());
			coords[1] = multiplier*std::stod(token);
			line.erase(0, position+delimiter.length());
			if (line.length() > 0){
				line.erase(std::remove(line.begin(), line.end(), ' '), line.end());
				line.erase(std::remove(line.begin(), line.end(), '\t'), line.end());
				coords[2] = multiplier*std::stod(line);
	
			} else {
				Error e2("INPUT", "Missing coordinate.");
				control->log.error(e2);
			}
		} else { 
			Error e3("INPUT", "Missing coordinate.");
			control->log.error(e3);
		} 
	} else {
		Error e4("INPUT", "Missing coordinate.");
		control->log.error(e4);
	}					
	
	return Atom(coords, q, m);
}

void Molecule::parseBasis(Construct& sc, std::map<int, std::string>& names) { 

	for (auto& line : sc.content) {
		size_t pos = line.find(',');
		std::string token; 
		
		if (pos != std::string::npos) {
			token = line.substr(0, pos);
			token.erase(std::remove(token.begin(), token.end(), ' '), token.end());
			std::transform(token.begin(), token.end(), token.begin(), ::tolower);
			std::string token2 = line.substr(pos+1, line.length());
			token2.erase(std::remove(token2.begin(), token2.end(), ' '), token2.end());
			if (token == "default")
				names[0] = token2;
			else {
				int q = getAtomCharge(token);
				names[q] = token2;
			}	
		}
	}

}

// An initialisation function, for code reuse purposes
void Molecule::init(Construct& c)
{
	std::vector<Option> options; 
	for (auto& line : c.content)
		options.push_back(Option(line));
	
	charge = 0;
	multiplicity = 1;
	parent = true; 
	angstrom = false;
	std::string xyz_file = ""; 
	for (auto& op : options) {
		if (op.name == "charge") charge = op.get_value<int>();
		else if (op.name == "multiplicity") multiplicity = op.get_value<int>(); 
		else if (op.name == "angstrom") angstrom = true; 
		else if (op.name == "geomxyz") xyz_file = op.get_value<std::string>(); 
	}
	
	if (xyz_file != "") {
		control->log.print("\nReading geometry from " + xyz_file); 
		std::vector<std::string> xyz_coords = read_xyz(xyz_file); 
		angstrom = true; // Assume angstroms
		
		if (xyz_coords.size() == 0) {
			Error xyzerr("XYZ", "Unable to read xyz file.");
			control->log.error(xyzerr); 
		} else {
			std::string name = "geometry"; 
			c.subconstructs.push_back( Construct(xyz_coords, name) );
		}
	} 
	
	bool basis_found = false;
	bool geom_found = false; 
	fragmented = false; 
	has_ecps = false; 
	std::vector<std::vector <int>> frags;
	iVector qs; 
	std::string token; 
	
	for (auto& sc : c.subconstructs) {
		if (sc.name == "basis" || sc.name == "jkbasis" || sc.name == "ribasis") {
			// Parse basis
			if (sc.content.size() > 0) {
				if (sc.name == "basis") {
					basis_found = true; 
					parseBasis(sc, bnames); 
				} else if (sc.name == "jkbasis") 
					parseBasis(sc, jkbnames);
				else if (sc.name == "ribasis")
					parseBasis(sc, ribnames); 
				
			} else {
				Error e("INIT", "Basis specification is empty!"); 
				control->log.error(e); 
			}
				
		} else if (sc.name == "geometry" && !geom_found) {
			// Parse geometry
			natoms = sc.content.size(); 
			qs.resize(natoms);
			if (natoms > 0) {
				
				geom_found = true; 
				atoms = new Atom[natoms];
				nel = -charge; 
				
				int i = 0;
				for (auto& line : sc.content) {
					atoms[i] = parseGeomLine(line); 
					qs[i] = atoms[i].getCharge(); 
					nel += atoms[i].getCharge(); 
					i++;
				}
				
			} else {
				Error e("INIT", "There are no atoms!");
				control->log.error(e); 
				atoms = NULL;
			}
			
		} else if (sc.name == "fragments") {
			// Parse fragments
			if (sc.content.size() > 0) {
				fragmented = true;
				
				for (auto& line : sc.content) {
					size_t pos = line.find(',');
					std::vector<int> f;
					while (pos != std::string::npos) {
						token = line.substr(0, pos);
						f.push_back(std::stoi(token));
						line.erase(0, pos+1);
						pos = line.find(',');
					}
					f.push_back(std::stoi(line));
					if (f.size() > 0) {
						f[0] -= 1;
						frags.push_back(f);
					}
				}
			}
		} else if (sc.name == "ecp"){ 
			// Parse ECPs
			if (sc.content.size() > 0) {
				has_ecps = true; 
				
				for (auto& line : sc.content) {
					size_t pos = line.find(',');
					if (pos != std::string::npos) {
						token = line.substr(0, pos);
						token.erase(std::remove(token.begin(), token.end(), ' '), token.end());
						std::transform(token.begin(), token.end(), token.begin(), ::tolower);
						std::string token2 = line.substr(pos+1, line.length());
						token2.erase(std::remove(token2.begin(), token2.end(), ' '), token2.end());
						int q = getAtomCharge(token);
						bnames[-q] = token2;
					}
				}
			}	
		} else if (sc.name == "frozen") { 
			if (sc.content.size() > 0) {
				for (auto& line : sc.content) {
					size_t pos = line.find(',');
					while (pos != std::string::npos) {
						token = line.substr(0, pos);
						frozen_atoms.push_back(std::stoi(token)-1);
						line.erase(0, pos+1);
						pos = line.find(',');
					}
					frozen_atoms.push_back(std::stoi(line) - 1);
				}
				
				std::sort(frozen_atoms.begin(), frozen_atoms.end(), std::less<int>());
			}
		}else {
			Error e("INIT", "Construct not recognised."); 
			control->log.error(e);
		}
	}

	// Get the basis set     
	if (basis_found && geom_found) {
		bfset = Basis(bnames, has_ecps); 
		
		if (fragmented) {
			for(auto& f : frags) 
				if (f.size() > 3)
					fragments.push_back(std::make_shared<Fragment>(control, &atoms[f[0]], f[1] - f[0], bfset, bnames, jkbnames, ribnames, has_ecps, f[2], f[3]));
		}
	}

}

void Fragment::init(Atom* as, int nat, int q, int mult, bool _has_ecps, const Basis& _bfset)
{
	charge = q;
	parent = false;
	multiplicity = mult;
	atoms = as; 
	natoms = nat;
	has_ecps = _has_ecps; 
	
	if (nat <= 0) {
		Error e("FRAGINIT", "There are no atoms in this fragment!");
		control->log.error(e);
	} else {
		bfset = _bfset; 
		
		nel = -charge;
		for (int i = 0; i < natoms; i++) 
			nel += atoms[i].getCharge();
	}
}

// Constructor
Molecule::Molecule(SharedPC _control, Construct& c, int q) : control(_control)
{
	// Set the log and then initialise
	init(c);
}

Molecule::Molecule(SharedPC _control, int q) : control(_control) { }

// Copy constructor
Molecule::Molecule(const Molecule& other) : control(other.control)
{
	charge = other.charge;
	nel = other.nel;
	multiplicity = other.multiplicity;
	natoms = other.natoms;
	parent = other.parent;
	angstrom = other.angstrom;
	fragmented = other.fragmented;
	has_ecps = other.has_ecps;
	fragments = other.fragments; 
	if (natoms > 0) {
		atoms = new Atom[natoms];
		for (int i = 0; i < natoms; i++) atoms[i] = other.atoms[i]; 
	}
	frozen_atoms = other.frozen_atoms; 
		
	bfset = other.bfset;
	ecpset = other.ecpset;
	bnames = other.bnames;
	jkbnames = other.jkbnames; 
	ribnames = other.ribnames; 
	enuc = other.enuc;
}


// Destructor
Molecule::~Molecule()
{ 
	if(natoms != 0 && parent){
		delete [] atoms;
	}
}

Molecule& Molecule::operator=(const Molecule& other) {
	control = other.control;
	charge = other.charge;
	nel = other.nel;
	multiplicity = other.multiplicity;
	natoms = other.natoms;
	parent = other.parent;
	angstrom = other.angstrom;
	fragmented = other.fragmented;
	has_ecps = other.has_ecps;
	fragments = other.fragments; 
	if (natoms > 0) {
		atoms = new Atom[natoms];
		for (int i = 0; i < natoms; i++) atoms[i] = other.atoms[i]; 
	}
		
	frozen_atoms = other.frozen_atoms; 	
	bfset = other.bfset;
	ecpset = other.ecpset;
	bnames = other.bnames;
	jkbnames = other.jkbnames;
	ribnames = other.ribnames;  
	enuc = other.enuc;
	return *this; 
}

Fragment::Fragment(SharedPC control, Atom* as, int nat, const Basis& _bfset, std::map<int, std::string> _bnames, std::map<int, std::string> _jkbnames, std::map<int, std::string> _ribnames, bool _has_ecps, int q, int mult) : Molecule(control, q)  {
	bnames = _bnames; 
	jkbnames = _jkbnames; 
	ribnames = _ribnames; 
	init(as, nat, q, mult, _has_ecps, _bfset); 
}

Fragment::Fragment(const Fragment& other) : Molecule(other.control, other.charge) {
	bnames = other.bnames; 
	jkbnames = other.jkbnames; 
	ribnames = other.ribnames; 
	init(other.atoms, other.natoms, other.charge, other.multiplicity, other.has_ecps, other.bfset);
}
Fragment::~Fragment() { }

Fragment& Fragment::operator=(const Fragment& other) {
	control = other.control;
	bnames = other.bnames;
	jkbnames = other.jkbnames; 
	ribnames = other.ribnames; 
	init(other.atoms, other.natoms, other.charge, other.multiplicity, other.has_ecps, other.bfset);
	return *this;
}

// Routines


//Get number of core/valence electrons in molecule

int Molecule::getNCore() const 
{
	return nel - getNValence(); 
}

int Molecule::getNValence() const
{
	int nvalence = -charge; 
	
	for (int i = 0; i < natoms; i++)
		nvalence += getAtomValence(atoms[i].getCharge());
	
	return nvalence; 
}

std::vector<int> Molecule::getActiveList() const {
	
	std::vector<int> list; 
	
	int ninactive = frozen_atoms.size(); 
	int currix = ninactive > 0 ? 0 : -1; 
	for (int i = 0; i < natoms; i++) {
		int curr = currix > -1 ? frozen_atoms[currix] : -1; 
		if (i == curr) {
			currix++; 
			currix = currix < ninactive ? currix : -1;  
		} else {
			list.push_back(i); 
		}
	}
	
	return list; 
	
}

// rotate(Matrix U) rotates the coordinate system according
// to the transformation specified by the unitary matrix U
// U needs to be 3x3
void Molecule::rotate(const Matrix& U)
{
	// Check if 3 x 3
	if(U.rows() == 3 && U.cols() == 3){
		for (int i = 0; i < natoms; i++){
			atoms[i].rotate(U); // Do the rotation
		}
		updateBasisPositions();
	} else {
		Error e("ROTATE", "Unsuitable rotation matrix.");
		control->log.error(e);
	}
}

// translate(x, y, z) translates the coordinate system by
// adding the vector (x, y, z) to the coordinate vector
// of each atom.
void Molecule::translate(double x, double y, double z)
{
	for (int i = 0; i < natoms; i++){
		atoms[i].translate(x, y, z);
	}
	updateBasisPositions();
}

// nalpha and nbeta return the number of alpha and beta
// (spin up/spin down) electrons
int Molecule::nalpha() const
{
	// multiplicity = 2s+1; each unpaired alpha electron contributes
	// s = 1/2, therefore s = nunpairedalpha*(1/2), so
  
	return (nel-multiplicity+1)/2 + multiplicity -1;
}

int Molecule::nbeta() const
{
	// The number of beta electrons is (nel - nunpairedalpha)/2
	return (nel - multiplicity + 1)/2;
}

// Calculate the nuclear energy, i.e. the sum of terms
// (Zi * Zj)/Rij for each distinct pair of nucleii, i, j
void Molecule::calcEnuc()
{
	double zi;
	enuc = 0.0;

	// Outer loop over all atoms
	for (int i = 0; i < natoms; i++) {
		zi = atoms[i].getEffectiveCharge();
		// Inner loop over all atoms with index
		// greater than i, to avoid double counting
		for (int j = i+1; j < natoms; j++){
			enuc += (zi*atoms[j].getEffectiveCharge())/dist(i, j);
		}
	}
}
  
// Compute the centre of mass of the molecule
Vector Molecule::com() const 
{
	// Coordinate vector of centre of mass, set to all
	// zeroes for summations
	Vector c = Vector::Zero(3);
  
	// Loop over all atoms
	double sum = 0.0;
	double m;
	for (int i = 0; i < natoms; i++){
		m = atoms[i].getMass();
		c = c + m*atoms[i].getCoords();
		sum += m;
	}

	c = (1.0/sum)*c;
	return c;
}

// Calculate the principal moments of inertia;
// optionally, shift to the inertial coordinate system
// First the elements of the inertia tensor are
// calculated, and then this is diagonalised to 
// give the eigenvalues - the principal moments - 
// which are returned as a 3-vector
Vector Molecule::getInertia(bool shift) 
{
	Vector v(3); // To return the moments in
	Matrix inertia = Matrix::Zero(3, 3); // To contain the elements

	// We must translate to centre of mass coordinates
	Vector c = -1.0*com();
	translate(c(0), c(1), c(2));

	// Loop over all atoms
	double m;
	for (int i = 0; i < natoms; i++){
		m = atoms[i].getMass();
		c = atoms[i].getCoords();

		// The diagonal elements are the mass times the sum
		// of the squares of the other coordinates
		// e.g. I(0, 0) = m(y^2 + z^2)
		inertia(0, 0) += m*(c(1)*c(1) + c(2)*c(2));
		inertia(1, 1) += m*(c(0)*c(0) + c(2)*c(2));
		inertia(2, 2) += m*(c(0)*c(0) + c(1)*c(1));

		// The off diagonal elements are minus the 
		// mass times the two coordinates
		// e.g. I(1, 0) = -m*x*y = I(0 , 1)
		inertia(0, 1) -= m*c(0)*c(1);
		inertia(0, 2) -= m*c(0)*c(2);
		inertia(1, 2) -= m*c(1)*c(2);
	}
	// Finally fill out the symmetric elements
	inertia(1, 0) = inertia(0, 1); 
	inertia(2, 0) = inertia(0, 2); 
	inertia(2, 1) = inertia(1, 2);
    
	// Now we diagonalise, only need eigenvalues
	EigenSolver es(inertia);
	v = es.eigenvalues();
	if(shift)
		rotate(es.eigenvectors());
	
	return v;
}


// Calculate the (Euclidean) distance between the atoms i and j
// dist2 gives the square of this distance
double Molecule::dist2(int i, int j) const
{
	Vector dvec = atoms[j].getCoords() - atoms[i].getCoords(); 
	return dvec.dot(dvec);
}

double Molecule::dist(int i, int j) const {
	return sqrt(dist2(i, j));
}

// Calculate bond, out of plane, and torsional angles
double Molecule::bondAngle(int i, int j, int k) const
{
	// Get the differences of the position vectors
	// and calculate the angle between them
	Vector rij = atoms[j].getCoords();
	Vector rjk = rij;
	rij = rij - atoms[i].getCoords();
	rjk = atoms[k].getCoords() - rjk;
	return angle(rij, rjk); // Vector class friend function
}

double Molecule::oopAngle(int i, int j, int k, int l) const
{
	// The out of plane angle is given by the formula:
	// [(elj x elk).eli]/sin(anglejlk)
	// and sin(anglejlk) = ||ejl x elk||
	// therefore start by forming the unit vectors
	Vector co = atoms[l].getCoords();
	Vector elk = atoms[k].getCoords() - co;
	Vector elj = atoms[j].getCoords() - co;
	Vector eli = atoms[i].getCoords() - co;
  
	co = cross(elj, elk);
	double sinangle = co.norm();
  
	// Normalise
	elk.normalize();
	elj.normalize();
	eli.normalize();
  
	// Do the triple product and divide by the sine
	return (cross(elj, elk)).dot(eli)/sinangle;
}

double Molecule::torsionAngle(int i, int j, int k, int l) const
{
	// The torsional angle is given by the formula:
	// [(eij x ejk).(ejk x ekl)]/[sin(ijk)*sin(jkl)]

	// Form the identity vectors
	Vector jv = atoms[j].getCoords();
	Vector kv = atoms[k].getCoords();
	Vector eij = jv - atoms[i].getCoords();
	Vector ejk = kv - jv;
	Vector ekl = atoms[l].getCoords() - kv;
  
	// Normalise
	eij.normalize();
	ejk.normalize();
	ekl.normalize();

	// Get angles and cross products
	double sinijk = (cross(eij, ejk)).norm();
	double sinjkl = (cross(ejk, ekl)).norm();
	jv = cross(eij, ejk);
	kv = cross(ejk, ekl);
 
	// Return the inner product over the product of angles
	return jv.dot(kv)/(sinijk*sinjkl);
}

// Compute the rotational type and constants of the molecule
// using the getInertia method
std::string Molecule::rType()
{
	std::string rstring = "asymmetric"; // Asymm. is default case
	// Diatomic case
	if(natoms == 2){
		rstring = "diatomic";
	} else {
		// Calculate the principal moments of inertia
		Vector I = getInertia(false);

		// All moments are >= 0, so no need to use fabs
		// Linear if Ic = Ib > Ia = 0
		// Prolate symm. if Ic = Ib > Ia /= 0
		// Oblate symm. if Ic > Ib = Ia /= 0
		// Spherical if Ic = Ib = Ia /= 0
		double CUTOFF = 0.01;

		if ( (I(2)-I(1)) < CUTOFF && I(0) < CUTOFF ){
			rstring = "linear";
		}   else if ( (I(2) - I(1)) < CUTOFF && (I(1)-I(0)) > CUTOFF &&  I(0) >= CUTOFF){
			rstring = "prolate";
		} else if ( (I(1) - I(0)) < CUTOFF && (I(2) - I(1)) > CUTOFF &&  I(0) >= CUTOFF){
			rstring = "oblate";
		} else if ( (I(1) - I(0)) < CUTOFF && (I(2) - I(1)) < CUTOFF){
			rstring = "spherical";
		}
	}
	return rstring;
}

Vector Molecule::rConsts(int units)
{
	// First get the principal moments of inertia
	Vector I = getInertia(false);
  
	double K;
	// Choose units
	if (units == 0){ // reciprocal centimetres
		K = Logger::RTOCM; 
	} else { // GHz
		K = Logger::RTOGHZ; // ~8.39214994 x 10^2 GHz . A^2 . amu
	}
	// Constant Bi = K/Ii 
	for (int i = 0; i < 3; i++){
		if (I(i) < control->get_option<double>("precision")) { I[i] = 0.0; }  
		else { I[i] = K/I(i); }
	}
	return I;
}

void Molecule::buildShellBasis() {
	BasisReader b(bnames);
	BasisReader b2(jkbnames);
	BasisReader b3(ribnames);  
	
	for (int i = 0; i < natoms; i++) {
		Atom &a = atoms[i];
		double pos[3] = { a.getX(), a.getY(), a.getZ() };
		b.readShellBasis(bfset, a.getCharge(), pos, i, 0);
		
		if (jkbnames.size() != 0) 
			b2.readShellBasis(bfset, a.getCharge(), pos, i, 1);
		
		if (ribnames.size() != 0)
			b3.readShellBasis(bfset, a.getCharge(), pos, i, 2);
	}
	
}

void Molecule::buildECPBasis() {
	if (bfset.hasECPS()) {
		BasisReader b(bnames); 
		for (int i = 0; i < natoms; i++) {
			Atom& a = atoms[i];
			auto it = bnames.find(-a.getCharge());
			if (it != bnames.end()) {
				libecpint::ECP newECP = b.readECP(a.getCharge(), ecpset, a.getPos());
				ecpset.addECP(newECP, i);
			}
			atoms[i].setCore(ecpset);
			nel -= atoms[i].getCharge() - atoms[i].getEffectiveCharge();
		}
	}
}

void Molecule::updateBasisPositions() {
	std::vector<libint2::Shell>& shells = bfset.getIntShells();
	std::vector<libint2::Shell>& jkshells = bfset.getJKShells();
	std::vector<libint2::Shell>& rishells = bfset.getRIShells(); 
	
	for (int i = 0; i < natoms; i++) {
		auto x = atoms[i].getX(); auto y = atoms[i].getY(); auto z = atoms[i].getZ();
		for (int j = 0; j < shells.size(); j++) {
			if (bfset.getShellAtom(j) == i) {
				shells[j].O[0] = x;
				shells[j].O[1] = y;
				shells[j].O[2] = z;
			}
		}
		
		for (int j = 0; j < jkshells.size(); j++) {
			if (bfset.getJKShellAtom(j) == i) {
				jkshells[j].O[0] = x;
				jkshells[j].O[1] = y;
				jkshells[j].O[2] = z;
			}
		}
		
		for (int j = 0; j < rishells.size(); j++) {
			if (bfset.getRIShellAtom(j) == i) {
				rishells[j].O[0] = x;
				rishells[j].O[1] = y;
				rishells[j].O[2] = z;
			}
		}
		
		if (bfset.hasECPS()) {
			for (int j = 0; j < ecpset.getN(); j++) {
				if (ecpset.getAtom(j) == i) ecpset.getECP(j).setPos(x, y, z);
			}
		}
	}
}
  


  
  
