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

#ifndef MOLECULEHEADERDEF
#define MOLECULEHEADERDEF

// Includes
#include "basis.hpp"
#include "atom.hpp"
#include <string>
#include <map>
#include <vector>
#include <libecpint/ecp.hpp>
#include "eigen_wrapper.hpp"

// Declare forward dependcies
class Fragment;
class Molecule;
class ProgramController;
struct Construct;

using SharedMolecule = std::shared_ptr<Molecule>; 
using SharedFragment = std::shared_ptr<Fragment>;  

// Begin class definition
class Molecule : public std::enable_shared_from_this<Molecule> 
{
protected:
  libecpint::ECPBasis ecpset;
  Basis bfset;
  Atom* atoms;
  std::vector<int> frozen_atoms; 
  int charge, nel, multiplicity, natoms;
  bool parent, angstrom, fragmented, has_ecps;
  std::vector<SharedFragment> fragments;
  std::map<int, std::string> bnames, jkbnames, ribnames;
  double enuc;
public:
	
  std::shared_ptr<ProgramController> control;
  
  // Constructors and destructor
  void init(Construct& c); // An initialisation function
  Molecule(std::shared_ptr<ProgramController> control, Construct& c, int q = 0); // Need the log for input, q is charge
  Molecule(std::shared_ptr<ProgramController> control, int q); 
  Molecule(const Molecule& other); // Copy constructor
  ~Molecule(); // Deletes the atom array
  
  Atom parseGeomLine(std::string line);
  void parseBasis(Construct& sc, std::map<int, std::string>& names); 
  
  void buildShellBasis();
  void buildECPBasis();
  void updateBasisPositions();
  
  // Accessors
  int getNAtoms() const { return natoms; }
  int getNFrozenAtoms() const { return frozen_atoms.size(); }
  int getNActiveAtoms() const { return natoms - frozen_atoms.size(); }
  int getCharge() const { return charge; }
  int getNel() const { return nel; }
  int getMultiplicity() const { return multiplicity; }
  int getNCore() const; 
  int getNValence() const;
  double getEnuc() const { return enuc; }
  Atom& getAtom(int i) { return atoms[i]; } // Return atom i
  std::vector<int> getActiveList() const; 
  Basis& getBasis() { return bfset; }
  libecpint::ECPBasis& getECPBasis() { return ecpset; }
  libecpint::ECP& getECP(int i) { return ecpset.getECP(i); }
  bool hasECPS() { return has_ecps; }
  std::vector<SharedFragment>& getFragments() { return fragments; }
  
  // Routines
  void rotate(const Matrix& U);
  void translate(double x, double y, double z);
  int nalpha() const;
  int nbeta() const;
  void calcEnuc(); 
  Vector com() const;
  Vector getInertia(bool shift = false);
  double dist(int i, int j) const;
  double dist2(int i, int j) const;
  double bondAngle(int i, int j, int k) const;
  double oopAngle(int i, int j, int k, int l) const;
  double torsionAngle(int i, int j, int k, int l) const;
  std::string rType();
  Vector rConsts(int units);
  
  Molecule& operator=(const Molecule& other); 
};

class Fragment : public Molecule 
{
private:
	std::vector<Atom> frag_atoms; 
public:
	Fragment(std::shared_ptr<ProgramController> control, Atom* as, int nat, const Basis& _bfset, std::map<int, std::string> _bnames, std::map<int, std::string> _jkbnames, std::map<int, std::string> _ribnames, bool _has_ecps = false, int q = 0, int mult = 1); 
	Fragment(const Fragment& other);
	~Fragment();
	void init(Atom* as, int nat, int q, int mult, bool _has_ecps, const Basis& _bfset);
	
	Fragment& operator=(const Fragment& other); 
};

#endif
