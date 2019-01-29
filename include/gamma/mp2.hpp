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

#ifndef MP2HEADERDEF
#define MP2HEADERDEF

#include "tensor4.hpp"
#include "fock.hpp"
#include <ctf.hpp>

class IntegralEngine;

class Integrals {
  public:
  CTF::World * dw;
  
  CTF::Tensor<> * aa;
  CTF::Tensor<> * ii;
  CTF::Tensor<> * ab;
  CTF::Tensor<> * ai;
  CTF::Tensor<> * ia;
  CTF::Tensor<> * ij;
  
  CTF::Tensor<> * abcd;
  CTF::Tensor<> * abci;
  CTF::Tensor<> * abic;
  CTF::Tensor<> * aibc;
  CTF::Tensor<> * iabc;
  CTF::Tensor<> * abij;
  CTF::Tensor<> * aibj; 
  CTF::Tensor<> * aijb;
  CTF::Tensor<> * iajb;
  CTF::Tensor<> * ijab; 
  CTF::Tensor<> * iabj; 
  CTF::Tensor<> * aijk; 
  CTF::Tensor<> * iajk; 
  CTF::Tensor<> * ijak;
  CTF::Tensor<> * ijka;
  CTF::Tensor<> * ijkl;     

  Integrals(int no, int nv, CTF::World &dw_){
    int shapeSYSY[] = {SY,NS,SY,NS};
    int shapeSYNS[] = {SY,NS,NS,NS};
    int shapeNSNS[] = {NS,NS,NS,NS};
    int shapeNSSY[] = {NS,NS,SY,NS};
    int vvvv[]      = {nv,nv,nv,nv};
    int vvvo[]      = {nv,nv,nv,no};
	int vvov[] 		= {nv,nv,no,nv};
    int vovv[]      = {nv,no,nv,nv};
    int ovvv[]		= {no,nv,nv,nv};
	int vvoo[]		= {nv,nv,no,no};
	int vovo[]		= {nv,no,nv,no};
	int ovvo[]		= {no,nv,nv,no};
	int voov[]		= {nv,no,no,nv};
	int ovov[]		= {no,nv,no,nv};
	int oovv[]		= {no,no,nv,nv};
	int vooo[]		= {nv,no,no,no};
	int ovoo[]		= {no,nv,no,no};
	int oovo[]		= {no,no,nv,no};
	int ooov[] 		= {no,no,no,nv};
	int oooo[]		= {no,no,no,no}; 
    
    dw = &dw_;
    
    ab = new CTF_Matrix(nv,nv,SY,dw_,"Vab",1);
    ai = new CTF_Matrix(nv,no,NS,dw_,"Vai",1);
    ia = new CTF_Matrix(no,nv,NS,dw_,"Via",1);
    ij = new CTF_Matrix(no,no,SY,dw_,"Vij",1);

    abcd = new CTF::Tensor<>(4,vvvv,shapeSYSY,dw_,"Vabcd",1);
    abci = new CTF::Tensor<>(4,vvvo,shapeSYNS,dw_,"Vabci",1);
	abic = new CTF::Tensor<>(4,vvov,shapeSYNS,dw_,"Vabic",1);
    aibc = new CTF::Tensor<>(4,vovv,shapeNSSY,dw_,"Vaibc",1);
	iabc = new CTF::Tensor<>(4,ovvv,shapeNSSY,dw_,"Viabc",1);
	abij = new CTF::Tensor<>(4,vvoo,shapeSYSY,dw_,"Vabij",1);
    aibj = new CTF::Tensor<>(4,vovo,shapeNSNS,dw_,"Vaibj",1);
	aijb = new CTF::Tensor<>(4,voov,shapeNSNS,dw_,"Vaijb",1);
	iajb = new CTF::Tensor<>(4,ovov,shapeNSNS,dw_,"Viajb",1);
	ijab = new CTF::Tensor<>(4,oovv,shapeSYSY,dw_,"Vijab",1);
	iabj = new CTF::Tensor<>(4,ovvo,shapeNSNS,dw_,"Viabj",1);
	aijk = new CTF::Tensor<>(4,vooo,shapeNSSY,dw_,"Vaijk",1);
	iajk = new CTF::Tensor<>(4,ovoo,shapeNSSY,dw_,"Viajk",1);
	ijak = new CTF::Tensor<>(4,oovo,shapeSYNS,dw_,"Vijak",1);
	ijka = new CTF::Tensor<>(4,ooov,shapeSYNS,dw_,"Vijka",1);
	ijkl = new CTF::Tensor<>(4,oooo,shapeSYSY,dw_,"Vijkl",1);
	
  }

  ~Integrals(){
    delete ab;
    delete ai;
    delete ia;
    delete ij;
    
    delete abcd;
    delete abci;
	delete abic;
	delete aibc;
	delete iabc;
	delete abij;
	delete aibj;
	delete aijb;
	delete iajb;
	delete ijab;
	delete iabj;
	delete aijk;
	delete iajk;
	delete ijak;
	delete ijka;
	delete ijkl;
  }
  
  CTF::Idx_Tensor operator[](char const * idx_map_){
    int i, lenm, no, nv;
    lenm = strlen(idx_map_);
    char new_idx_map[lenm+1];
    new_idx_map[lenm]='\0';
    no = 0;
    nv = 0;
    for (i=0; i<lenm; i++){
      if (idx_map_[i] >= 'a' && idx_map_[i] <= 'h'){
        new_idx_map[i] = 'a'+nv;
        nv++;
      } else if (idx_map_[i] >= 'i' && idx_map_[i] <= 'n'){
        new_idx_map[i] = 'i'+no;
        no++;
      }
    }
//    printf("indices %s are %s\n",idx_map_,new_idx_map);
    if (0 == strcmp("ab",new_idx_map)) return (*ab)[idx_map_];
    if (0 == strcmp("ai",new_idx_map)) return (*ai)[idx_map_];
    if (0 == strcmp("ia",new_idx_map)) return (*ia)[idx_map_];
    if (0 == strcmp("ij",new_idx_map)) return (*ij)[idx_map_];
    if (0 == strcmp("abcd",new_idx_map)) return (*abcd)[idx_map_];
    if (0 == strcmp("abci",new_idx_map)) return (*abci)[idx_map_];
	if (0 == strcmp("abic",new_idx_map)) return (*abic)[idx_map_];
    if (0 == strcmp("aibc",new_idx_map)) return (*aibc)[idx_map_];
	if (0 == strcmp("iabc",new_idx_map)) return (*iabc)[idx_map_];
    if (0 == strcmp("aibj",new_idx_map)) return (*aibj)[idx_map_];
    if (0 == strcmp("abij",new_idx_map)) return (*abij)[idx_map_];
	if (0 == strcmp("aijb",new_idx_map)) return (*aijb)[idx_map_];
	if (0 == strcmp("iabj",new_idx_map)) return (*iabj)[idx_map_];
	if (0 == strcmp("iajb",new_idx_map)) return (*iajb)[idx_map_];
    if (0 == strcmp("ijab",new_idx_map)) return (*ijab)[idx_map_];
    if (0 == strcmp("aijk",new_idx_map)) return (*aijk)[idx_map_];
	if (0 == strcmp("iajk",new_idx_map)) return (*iajk)[idx_map_];
    if (0 == strcmp("ijak",new_idx_map)) return (*ijak)[idx_map_];
	if (0 == strcmp("ijka",new_idx_map)) return (*ijka)[idx_map_];
    if (0 == strcmp("ijkl",new_idx_map)) return (*ijkl)[idx_map_];
    printf("Invalid integral indices\n");
	std::cout << new_idx_map << std::endl;  
    assert(0);
//shut up compiler
    return (*aa)[idx_map_];
  }
};

class Amplitudes {
  public:
  CTF::Tensor<> * ai;
  CTF::Tensor<> * abij;
  CTF::World * dw;

  Amplitudes(int no, int nv, CTF::World &dw_){
    dw = &dw_;
    int shapeNSNS[] = {NS,NS,NS,NS};
    int vvoo[]      = {nv,nv,no,no};

    ai = new CTF_Matrix(nv,no,NS,dw_,"Tai",1);
    abij = new CTF::Tensor<>(4,vvoo,shapeNSNS,dw_,"Tabij",1);
  }

  ~Amplitudes(){
    delete ai;
    delete abij;
  }

  CTF::Idx_Tensor operator[](char const * idx_map_){
    if (strlen(idx_map_) == 4) return (*abij)[idx_map_];
    else return (*ai)[idx_map_];
  }
};

class MP2
{
private:
	int N, nocc;
	double energy;
	S8EvenTensor4 moInts;
	std::shared_ptr<Integrals> spinInts;
	std::shared_ptr<Amplitudes> amplitudes;  
	bool spinBasis;
	Fock& focker;
public:
	CTF::World dw; 
	
	MP2(Fock& _focker);
	void transformIntegrals(bool withSpin = true);
	void transformThread(int start, int end, Tensor4& moTemp);
	void calculateEnergy();
	
	void tensormp2(bool print = true); 
	void cctrans(); 
	void dfmp2(bool print = true);
	
	double getEnergy() const { return energy; }
	std::shared_ptr<Integrals>& getSpinInts() { 
		return spinInts;
	}
	std::shared_ptr<Amplitudes>& getAmplitudes() {
		return amplitudes;
	}
	S8EvenTensor4& getMOInts() {
		return moInts;
	}
	int getN() const { return N; }
	int getNocc() const { return nocc; }
	Fock& getFock() { return focker; }
};

#endif 
