/* Implementation for error.hpp  
*
*     DATE             AUTHOR                CHANGES                             
*   =======================================================================          
*     14/08/15         Robert Shaw           Original code             
*
*/

#include "diis.hpp"
#include <iostream>
#include <cmath>

/// Random comment
DIISEngine::DIISEngine() { }

void DIISEngine::init(int _maxDiis, bool _useDiis, double _damp) {
	maxDiis = _maxDiis > 0 ? _maxDiis : 8;
	useDiis = _useDiis;
	damping_factor = _damp;
	lastB = Matrix::Zero(1, 1);
}

Vector DIISEngine::compute(std::vector<Vector> &errors) {
	
	Vector weights;
	int iter = errs.size() + 1;
	// Perform DIIS averaging
	// Direct inversion of the iterative subspace
	// Greatly improves convergence and numerical behaviour
	// of the scf iterations.
	if (useDiis) {

		if (iter > maxDiis) errs.erase(errs.begin());
		int length = 0;
		for (auto e : errors) length += e.size();
		Vector newErr(length);
		int i = 0;
		for (auto e : errors) {
			for (int j = 0; j < e.size(); j++) newErr[i++] = e[j];
		}
		errs.push_back(newErr);
		if (iter == 1) lastB(0, 0) = newErr.dot(newErr);
		
		if (iter > 1) {
			
			int lim = (iter < maxDiis ? iter : maxDiis);			
			Matrix B(lim+1, lim+1); // Error norm matrix
			B(lim, lim) = 0.0;

			// The elements of B are <e_i | e_j >
			int start = iter > lim ? 1 : 0;
			for (int i = 0; i < lim-1; i++){
				for (int j = i; j < lim-1; j++){
					B(i, j) = lastB(i+start, j+start);
					B(j, i) = B(i, j);
				}
				B(lim, i) = -1.0;
				B(i, lim) = -1.0;	
			}
			
			for (int i = 0; i < lim; i++) {
				B(lim-1, i) = errs[i].dot(errs[lim-1]);
				B(i, lim-1) = B(lim-1, i);
			}
			B(lim, lim-1) = -1.0;
			B(lim-1, lim) = -1.0;
			lastB = B;
		
			weights = solve(B, 0);
			int rank = weights[weights.size() - 1];
			
			start = 1;
			while (rank < (lim - start + 1) && start < lim) {
				weights = solve(B, start++);
				rank = weights[weights.size() - 1];
			}
			
			weights.conservativeResize(weights.size()-2);
		}

	}

	//weights.print();
	return weights;
}

Vector DIISEngine::solve(Matrix &B, int start) {
	
	int size = B.rows() - start; 
	double norm = 1.0;
	double scale = 1.0 + damping_factor;
	if (B(start, start) > 1e-10) norm = 1.0/B(start,start);

	Matrix newB(size, size);
	for (int i = 0; i < size-1; i++) {
		newB(i, i) = B(i+start, i+start) * scale;
		for (int j = i; j < size-1; j++) {
			newB(i, j) = norm * B(i+start, j+start);
			newB(j, i) = newB(i, j);
		}
		newB(i, size-1) = newB(size-1, i) = -1.0;
	}
	newB(size-1, size-1) = 0.0;
	
	// Solve the linear system of equations for the weights
	Vector w = Vector::Zero(size); 
	w[size-1] = -1.0;
	
	FullLU lu(newB);
	w = lu.solve(w);
	
	Vector weights(size+1);
	for (int i = 0; i < size; i++) weights[i] = w(i);
	weights[size] = lu.rank();
	
	return weights;
}

