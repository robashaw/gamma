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

#ifndef RPAHEADERDEF
#define RPAHEADERDEF

#include <fock.hpp>
#include "ProgramController.hpp"
#include <vector>
#include "eigen_wrapper.hpp"
#include <libint2.hpp>

struct fInfo {

	std::vector<int> nocc, nvirt, ncore;
	Matrix F, S, T, V; 
	std::vector<libint2::Shell>& shells;
	std::vector<libint2::Shell>& df_shells; 
	
	double eintra, edisp, edispexch, eionic, ebsse;  
	
	fInfo(std::vector<libint2::Shell>& s1, std::vector<libint2::Shell>& s2) : shells(s1), df_shells(s2) {} 

};


class RPA
{
private:
	Command& cmd; 
	Fock& focker; 
	
	int N, nocc, nvirt; 
	double energy; 

public:
	CTF::World& dw;
	
	RPA(Command& c, Fock& f, int _N, int _nocc, CTF::World& _dw);
	
	void eris(const std::vector<libint2::Shell>&, std::vector<CTF::Tensor<> >&, Matrix& T, Matrix& V, bool longrange = false); 
	void df_eris(const std::vector<libint2::Shell>& obs, const std::vector<libint2::Shell>& auxbs, Matrix& T, Matrix& V, Matrix& BmnP); 
	
	void compute(bool print = true); 
	void fcompute(fInfo& info, bool print = false); 
	
	double getEnergy() const { return energy; }
};

#endif 
