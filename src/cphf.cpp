#include "cphf.hpp"
#include <iostream>

Vector cphf_single_solver(Vector& Q, Matrix& V, double tol, int maxiter) {
	
	Vector un = Q; 
	
	int iter = 0; 
	bool converged = false;
	Vector un1;
	double delta; 

	while (!converged && iter < maxiter) {
		Vector un1 = Q + V * un; 
		delta = (un1 - un).norm(); 
		un = un1;
		
		converged = delta < tol; 
		iter++; 
	}
	
	if (!converged) std::cout << "CPHF failed to converge" << std::endl;
	return un; 

}

std::vector<Vector> cphf_group_solver(std::vector<Vector>& Q, Matrix& V, double tol, int maxiter) {
	
	std::vector<Vector> uns; 
	
	for (Vector& q : Q) uns.push_back(cphf_single_solver(q, V, tol, maxiter)); 
	
	return uns; 
	
}
