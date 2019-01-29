/*
 *
 *     PURPOSE: To define the class Basis representing a
 *              basis set.
 *
 *     class Basis:
 *              owns: bfs - a set of BFs
 *              data: name - the name of the basis set, needed for file io
 *                    charges - a list of atomic numbers corresponding to
 *                            each BF in bfs.
 *                    shells - a list representing which bfs are in which
 *                             shells - i.e. if the first 3 are in one shell,
 *                             the next two in another, etc, it would be
 *                             3, 2, ...
 *                    lnums - a list with the angular quantum number 
 *                            of each shell in shells                 
 *              routines:
 *                    findPosition(q) - find the first position in the bfs 
 *                                      array of an atom of atomic number q
 *                    findShellPosition(q) - find the first position in the 
 *                                           shells array
 *                    getNBFs() - returns the total number of basis functions
 *                    getName() - returns the name
 *                    getCharges() - returns the vector of charges
 *                    getBF(charge, i) - return the ith BF corresponding to
 *                                       an atom of atomic number charge   
 *                    getSize(charge) - return how many basis functions an atom
 *                                      of atomic number charge has
 *                    getShellSize(charge) - return how many shells it has
 *                    getShells(charge) - return a vector of the correct subset
 *                                        of shells
 *                    getLnums(charge) - return vector of subset of lnums
 *                 
 *
 *     DATE        AUTHOR            CHANGES
 *     ====================================================================
 *     27/08/15    Robert Shaw       Original code.
 *
 */

#ifndef BASISHEADERDEF
#define BASISHEADERDEF

// Includes
#include "eigen_wrapper.hpp"
#include <string>
#include <libint2.hpp>
#include <vector>
#include <map>
#include "libecpint/gshell.hpp"

// Forward declarations
class BF;

// Begin class definitions

class Basis
{
private:
  std::map<int, std::string> names;
  std::string name;  
  int maxl, nexps; 
  bool ecps;
  
  std::vector<std::vector <libint2::Shell::Contraction>> raw_contractions;  
  std::vector<libint2::Shell> intShells;
  std::vector<libint2::Shell> jkShells;
  std::vector<libint2::Shell> riShells;  
  std::vector<int> shellAtomList;
  std::vector<int> jkShellAtomList;
  std::vector<int> riShellAtomList;
  
public:
  // Constructors and destructor
  // Note - no copy constructor, as it doesn't really seem necessary
  // Need to specify the name of the basis, n, and a list of the 
  // distinct atoms that are needed (as a vector of atomic numbers)
  Basis() : name("Undefined"), nexps(-1) { } // Default constructor
  Basis(std::map<int, std::string> ns, bool _ecps = false);
  
  // Accessors

  bool hasECPS() const { return ecps; }
  std::string getName(int q) const; 
  std::string getName() const { return name; }
  std::map<int, std::string>& getNames() { return names; }
  int getShellAtom(int i) const { return shellAtomList[i]; }
  int getJKShellAtom(int i) const { return jkShellAtomList[i]; }
  int getRIShellAtom(int i) const { return riShellAtomList[i]; }
  std::vector<libint2::Shell>& getIntShells() { return intShells; }
  std::vector<libint2::Shell>& getJKShells() { return jkShells; }
  std::vector<libint2::Shell>& getRIShells() { return riShells; }
  std::vector<int>& getShellAtomList() { return shellAtomList; }
  std::vector<int>& getJKShellAtomList() { return jkShellAtomList; }
  std::vector<int>& getRIShellAtomList() { return riShellAtomList; }
  
  int getNExps(); 
  int getMaxL() const { return maxl; }
  void setMaxL(int l) { maxl = l; }
  double getExp(int i) const; 
  void setExp(int i, double value); 
  double extent() const; 
  
  void addShell(libecpint::GaussianShell& g, int atom, int type = 0);
  
  // Overloaded operators
  Basis& operator=(const Basis& other);
};

#endif
