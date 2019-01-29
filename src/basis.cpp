/*
 *    PURPOSE: To implement basis.hpp, defining class Basis, representing a basis set.
 *   
 *    DATE            AUTHOR              CHANGES
 *    =====================================================================
 *    27/08/15        Robert Shaw         Original code.
 *
 */

// Includes
#include "basis.hpp"
#include "basisreader.hpp"
#include "ioutil.hpp"
#include <iostream>
#include <array>

// Constructors and destructors
Basis::Basis(std::map<int, std::string> ns, bool _ecps) : names(ns), maxl(0), nexps(-1), ecps(_ecps)
{
	auto it = ns.find(0);
	if (it != ns.end()) name = it->second;
	else name = "sto-3g";
	std::transform(name.begin(), name.end(), name.begin(), ::toupper);
}

std::string Basis::getName(int q) const {
	auto it = names.find(q);
	std::string n = name;
	if (it != names.end()) n = it->second;
	return n; 
}

void Basis::addShell(libecpint::GaussianShell& g, int atom, int type) {
	using libint2::Shell;
	
	std::vector<libint2::Shell::Contraction> contr_arr; 
	libint2::Shell::Contraction newC;
	newC.l = g.l;
	newC.pure = true; 
	for (auto c : g.coeffs) newC.coeff.push_back(c);
	contr_arr.push_back(newC);
	
	std::array<double, 3> O = { {g.center()[0], g.center()[1], g.center()[2]} };

	switch(type) {
		
		case 1: {
			jkShells.push_back(Shell(g.exps, contr_arr, O)); 
			jkShellAtomList.push_back(atom); 
			break;
		}
		
		case 2: {
			riShells.push_back(Shell(g.exps, contr_arr, O));
			riShellAtomList.push_back(atom); 
			break;
		}
		
		default: {
			intShells.push_back(Shell(g.exps, contr_arr, O));
			shellAtomList.push_back(atom);
			raw_contractions.push_back(contr_arr);
		}
	}
}

  
// Overloaded operators
Basis& Basis::operator=(const Basis& other)
{
  // Assign attributes
  names = other.names;
  ecps = other.ecps;
  name = other.name; 
  
  raw_contractions = other.raw_contractions; 
  intShells = other.intShells; 
  jkShells = other.jkShells; 
  riShells = other.riShells;
  shellAtomList = other.shellAtomList; 
  jkShellAtomList = other.jkShellAtomList; 
  riShellAtomList = other.riShellAtomList;

  return *this;
}

int Basis::getNExps() {
	
	if (nexps == -1) {
		nexps = 0; 
		
		for (auto& s : intShells)
			nexps += s.alpha.size();
	}

	return nexps; 
	
}

double Basis::getExp(int i) const {
	int ctr = 0; 
	bool found = false; 
	double ex = 0.0;
	
	for (auto& s : intShells) {
		for (auto x : s.alpha) {
			if (ctr == i) {
				ex = x; 
				found = true;
				break;
			}
			ctr++; 
		}
		if (found) break;
	}
	
	return ex; 
}

void Basis::setExp(int i, double value) {
	if (value > 0) {
		int ctr = 0;
			bool found = false; 
		for (int w = 0; w < intShells.size(); w++) {
			libint2::Shell& s = intShells[w]; 
			
			for (int x = 0; x < s.alpha.size(); x++) {
				if (ctr == i) {
					s.alpha[x] = value; 
					
					libint2::Shell newShell(s.alpha,
						raw_contractions[w],
						s.O);
					
					intShells[w] = newShell; 
					
					found = true;
					break;
				}
				ctr++; 
			}
			if (found) break; 
		}
	}	
}

double Basis::extent() const {
	double minex = -1.0; 
	for(auto& s : intShells)
		for (auto x : s.alpha)
			minex = minex < 0 || x < minex ? x : minex; 
	
	double extent = 2.52 * std::pow(minex, -0.526); 
	
	return extent; 
}

  


