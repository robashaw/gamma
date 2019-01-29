#include "optimiser.hpp"
#include "logger.hpp"
#include "scf.hpp"
#include "atom.hpp"
#include "integrals.hpp"
#include <vector>
#include <Eigen/Core>
#include <iomanip>

double rfo(Vector& dx, Vector& grad, Matrix& hessian, double alpha, double stepsize) {
	
	Vector step; 
	std::pair<double, double> results = rfo_prime(step, grad, hessian, alpha); 
	Vector m_dy; 
	double m_expect; 
	double step_norm = step.norm(); 
	if (step_norm > stepsize) {
		double v = alpha; 
		int niter = 0; 
		double ndy_last = 0.0; 
		
		double m_ndy = step_norm; 
		m_dy = step; 
		m_expect = results.first; 
		
		while(niter < 100) {
			v += (1-step_norm / stepsize) * step_norm / results.second; 
			results = rfo_prime(step, grad, hessian, v); 
			step_norm = step.norm(); 
			
			if (fabs(step_norm - stepsize)/stepsize < 0.001) {
				m_dy = step;
				m_expect = results.first;
				break; 
			}
			else if (niter > 10 && fabs(ndy_last - step_norm)/step_norm < 0.001) {
				m_dy = step;
				m_expect = results.first;
				break; 	
			}
			
			niter++; 
			ndy_last = step_norm; 
			if (step_norm < m_ndy) {
				m_ndy = step_norm;
				m_dy = step;
				m_expect = results.first; 
			}	
		}
	} else {
		m_dy = step; 
		m_expect = results.first; 
	}
	
	dx = m_dy; 
	return m_expect; 
}

double trust_newton(Vector& dx, Vector& grad, Matrix& hessian, double stepsize) {
	EigenSolver hsolve(hessian); 
	Vector gtilde = hsolve.eigenvectors().transpose() * grad; 
	Vector step = Vector::Zero(hessian.rows());
	for (int i = 0; i < hessian.rows(); i++) {
		double eival = hsolve.eigenvalues()[i]; 
		if (fabs(eival) > 1e-12)
			step[i] = - gtilde[i] / eival; 
	}
	
	dx = hsolve.eigenvectors() * step; 
	
	Matrix ident = Matrix::Identity(hessian.rows(), hessian.cols()); 
	if (dx.norm() > stepsize) {
		double eimin = hsolve.eigenvalues()[0]; 
		double factor = eimin < 0 ? 1.0 : 0.0; 
		double mu = (factor + 0.9)*eimin;
		
		Matrix hshift = hsolve.eigenvalues().asDiagonal();
		for (int i = 0; i < hshift.rows(); i++) {
			double hii = hshift(i, i) - mu;
			hshift(i, i) = 1.0 / (hii * hii); 
		}
		
		double delta = gtilde.transpose() * hshift * gtilde;
		delta = sqrt(delta) - stepsize; 
		
		for (int i = 1; i < 9; i++) {
			double new_mu = (factor + 0.1) * i * eimin; 

			for (int i = 0; i < hshift.rows(); i++) {
				double hii = hsolve.eigenvalues()[i] - new_mu;
				hshift(i, i) = 1.0 / (hii * hii); 
			}
		
			double new_delta = gtilde.transpose() * hshift * gtilde;
			new_delta = sqrt(new_delta) - stepsize; 
			
			if (fabs(new_delta) < fabs(delta)) {
				mu = new_mu;
				delta = new_delta;
			}
		}

		for (int i = 0; i < hessian.rows(); i++) {
			double eival = hsolve.eigenvalues()[i]; 
			if (fabs(eival - mu) > 1e-12)
				step[i] = - gtilde[i] / (eival - mu); 
			else step[i] = 0.0; 
		}
		
	}
	
	dx = hsolve.eigenvectors() * step; 
	double expect = -grad.dot(dx) + dx.transpose() * hessian * dx; 
	
	return expect; 
}

std::pair<double, double> rfo_prime(Vector& dy, Vector &grad, Matrix& hessian, double alpha) {
	
	int ndim = hessian.rows();
	int eival = 2; 
	
	Matrix augH = Matrix::Zero(ndim+1, ndim+1); 
	augH.block(0, 0, ndim, ndim) = hessian; 
	augH.block(ndim, 0, 1, ndim) = -grad.transpose();
	augH.block(0, ndim, ndim, 1) = -grad; 
	
	Matrix B = alpha * Matrix::Identity(ndim+1, ndim+1);
	B(ndim, ndim) = 1.0; 

	GeneralizedEigenSolver solver(augH, B); 
	double lmin = solver.eigenvalues()[eival]; 
	Vector step = solver.eigenvectors().block(0, eival, ndim, 1);
	std::cout << solver.eigenvectors() << std::endl << std::endl;
	double vmin = solver.eigenvectors()(ndim, eival); 
	if (fabs(vmin) > 1e-12)
		step /= vmin;
	
	double nu = alpha * lmin; 
	EigenSolver hsolve(hessian); 
	const Matrix& hvec = hsolve.eigenvectors();
	const Vector& heig = hsolve.eigenvalues(); 
	double dyprime2 = 0.0;
	double dy2 = 0.0; 
	double dotp, denom; 
	for (int i = 0; i < hessian.cols(); i++) {
		dotp = hvec.col(i).dot(grad); 
		dotp *= dotp; 
		denom = heig[i] - nu; 
		dyprime2 += dotp / (denom * denom * denom);
		dy2 += dotp / (denom * denom); 
	}
	double expect = 1 + alpha * step.dot(step); 
	dyprime2 *= 2.0 * lmin / expect; 
	expect *= 0.5 * lmin; 
	double dyprime = 0.5 * dyprime2 / sqrt(dy2); 
	
	dy = step; 
	std::pair<double, double> results = {expect, dyprime}; 
	return results; 
	
}



void RHFOptimiser::optimise() {
	Logger& log = mol->control->log;
	log.title("GEOMETRY OPTIMIZATION"); 
	log.initIterationOpt();
	int MAXITER = cmd.get_option<int>("maxiter");
	double CONVERGE = cmd.get_option<double>("gradconverge");  
	
	Vector dx = Vector::Zero(3*mol->getNActiveAtoms()); 
	
	bool converged = false;
	Vector grad;
	Matrix hessian; 
	calc.iter = 0;
	double trust = cmd.get_option<double>("trust"); 
	double trust_ratio;
	double expect = 0.0;
	double step_norm = 0.0;
	while (!converged && calc.iter < MAXITER) {
		calc(dx, grad, hessian); 
		if (calc.iter > 1) {
		
			trust_ratio = calc.delta_e / expect; 
			
			double stepnorm = dx.norm();
			if (trust_ratio > 0.75 && 1.25*stepnorm > trust) trust *= 2.0;
			else if (trust_ratio < 0.25) trust = 0.25*stepnorm; 
		}
		
		mol->control->log.optg_dump(calc.iter, grad, dx, mol, 
								hessian, trust, calc.delta_e, calc.grad_norm,
								step_norm, calc.energy, expect); 
		
		converged = (calc.grad_norm < CONVERGE) || (fabs(calc.delta_e) < CONVERGE / 1000.0); 
		if (!converged) {
			expect = trust_newton(dx, grad, hessian, trust); 
			step_norm = dx.norm();
		}
	}
	
	if (calc.iter < MAXITER) {
		log.print("Converged geometry:\n");
		for (int i = 0; i <mol->getNAtoms(); i++) log.print(mol->getAtom(i));
		log.result("HF Energy", calc.energy, "Hartree");
		
		if (cmd.get_option<bool>("freq")) frequencies(hessian); 
		
	} else {
		log.result("Geometry optimisation failed to converge.");
		log.print("Last trial geometry\n"); 
		for (int i = 0; i <mol->getNAtoms(); i++) log.print(mol->getAtom(i));
	}
}

void RHFOptimiser::frequencies(Matrix& hessian) {
	
	int natom = hessian.rows(); 
	
	// form mass-weighted hessian
	int row = 0; 
	double mi, mj; 
	Matrix mass_hessian = hessian; 
	
	std::vector<int> activeAtoms = mol->getActiveList(); 
	for (int i : activeAtoms) {
		mi = mol->getAtom(i).getMass(); 
		mi = sqrt(mi); 
		
		for (int xyz = 0; xyz < 3; xyz++) {
			
			int col = 0; 
			for (int j : activeAtoms) {
				mj = mol->getAtom(j).getMass();
				mj = sqrt(mj); 
				
				for (int abc = 0; abc < 3; abc++) {
					mass_hessian(row, col) /= mi * mj; 
					col++;	
				}
			}
			
			row++;
		}
	}
	
	EigenSolver hsolve(mass_hessian); 
	Vector freqs = hsolve.eigenvalues(); 
	Matrix modes = hsolve.eigenvectors(); 
	
	int nimaginary = 0; 
	for (int i = 0; i < hessian.rows(); i++) {
		if (freqs[i] < -1e-4) {
			freqs[i] = -sqrt(-freqs[i]); 
			nimaginary++; 
		} else if (freqs[i] < 1e-4) {
			freqs[i] = 0.0; 
		} else freqs[i] = sqrt(freqs[i]); 
		freqs[i] /= (2.0 * M_PI);  
		freqs[i] *= Logger::TOWAVENUMBERS; 
	}
	
	mol->control->log.frequencies(freqs, modes, cmd.get_option<bool>("modes"));  
	
}

void RHFOptimiser::exponents() { 
	
	int nexp = mol->getBasis().getNExps(); 
	
	Logger& log = mol->control->log;
	log.title("EXPONENT OPTIMIZATION"); 
	log.initIterationOpt();
	int MAXITER = cmd.get_option<int>("maxiter");
	double CONVERGE = cmd.get_option<double>("gradconverge"); 
	bool withmp2 = cmd.get_option<bool>("mp2");  
	
	std::string active = cmd.get_option<std::string>("active"); 
	std::vector<int> activex;
	if (active == "all") {
		for (int i = 0; i < nexp; i++) activex.push_back(i); 
	} else {
		size_t pos = active.find(',');
		std::string token; 
		while (pos != std::string::npos) {
			token = active.substr(0, pos);
			activex.push_back(std::stoi(token)-1);
			active.erase(0, pos+1);
			pos = active.find(',');
		}
		activex.push_back(std::stoi(active) - 1);
		
		nexp = activex.size(); 
	}
	
	int orbital = cmd.get_option<int>("orbital"); 
	bool momap = cmd.get_option<bool>("momap"); 
	Vector mocoeffs; 
	
	bool converged = false;
	Matrix xhessian;
	Vector xgrad; 
	Vector dxx = Vector::Zero(nexp); 
	int iter = 0;
	double trust = cmd.get_option<double>("trust"); 
	double trust_ratio, delta_e, grad_norm;
	double energy = 0.0;
	double expect = 0.0; 
	double step_norm = 0.0;  
	while (!converged && iter < MAXITER) {
		
		IntegralEngine ints(mol, false);
		Fock f(cmd, ints, mol);
		SCF hf(cmd, mol, f); 
		hf.rhf(false); 
		delta_e = hf.getEnergy() - energy;
		energy = hf.getEnergy();
		
		if (withmp2) {
			MP2 mp2obj(f);
			mp2obj.tensormp2(false);
			energy += mp2obj.getEnergy();
			delta_e += mp2obj.getEnergy(); 
		}
		 
		xgrad = f.compute_xgrad(energy, xhessian, activex, cmd); 
		grad_norm = xgrad.norm();
		if (momap) 
			mocoeffs = f.getCP().col(orbital); 
		
		if (iter > 1) {
		
			trust_ratio = delta_e / expect; 
			
			step_norm = dxx.norm();
			if (trust_ratio > 0.75 && 1.25*step_norm > trust) trust *= 2.0;
			else if (trust_ratio < 0.25) trust = 0.25*step_norm; 
		}
		
		mol->control->log.optx_dump(iter++, xgrad, dxx, mol, 
								xhessian, trust, delta_e, grad_norm,
								step_norm, energy, expect, activex); 
		
		converged = (grad_norm < CONVERGE) || (fabs(delta_e) < CONVERGE / 1000.0);  
		if (!converged) {
			expect = trust_newton(dxx, xgrad, xhessian, trust); 
			
			double currexp, newexp;
			int ctr = 0;
			for (int i : activex) {
				currexp = mol->getBasis().getExp(i); 
				newexp = currexp + dxx[ctr++]; 
				newexp = newexp <= 0.0 ? currexp / 10.0 : newexp; 
				mol->getBasis().setExp(i, newexp);  
			}
		}
			
	}
	
	if (calc.iter < MAXITER) {
		log.print("Converged exponents:\n");
		
		for (int i : activex)
			log.print(std::to_string(i) + ": " + std::to_string(mol->getBasis().getExp(i))); 
		
		log.result("HF Energy", energy, "Hartree");
		
		if (momap) {
			std::string mapfile = cmd.get_option<std::string>("mapfile");
			log.mo_map(mocoeffs, mol, cmd.get_option<int>("fineness"), mapfile);
		} 
		
	} else {
		log.result("Exponent optimisation failed to converge.");
		log.print("Last trial exponents\n"); 
		
		for (int i : activex)
			log.print(std::to_string(i) + ": " + std::to_string(mol->getBasis().getExp(i)));
	}

}

double RHFCalculator::operator()(const Vector& dx, Vector& grad, Matrix& hessian) {
		
	int offset = 0;
	int natoms = mol->getNAtoms();
	
	std::vector<int> activeAtoms = mol->getActiveList(); 
	for (int i : activeAtoms) {
		Atom& a = mol->getAtom(i); 
		a.translate(dx[offset], dx[offset+1], dx[offset+2]); 
		offset += 3; 
	}
	mol->updateBasisPositions(); 
	mol->calcEnuc(); 
		
	IntegralEngine ints(mol, false);
	Fock f(cmd, ints, mol);
	SCF hf(cmd, mol, f); 
	hf.rhf(false); 
	energy = hf.getEnergy();
		
	std::vector<Atom> atomlist; 
	for (int i = 0; i < natoms; i++) atomlist.push_back(mol->getAtom(i));
	
	//f.compute_forces(atomlist, mol->getNel()/2);
	f.compute_hessian(atomlist, mol->getNel()/2); 
	Matrix g = f.getForces().transpose(); 
	Vector full_grad = Eigen::Map<Vector>(g.data(), g.cols()*g.rows());
	
	//f.compute_hessian_numeric(atomlist, mol->getNel()/2, cmd); 
	if (natoms == activeAtoms.size()) {
		hessian = f.getHessian(); 
		grad = full_grad; 
	} else {
		Matrix& full_hessian = f.getHessian(); 

		// Restrict to only active atoms
		int nactive = activeAtoms.size(); 
		grad = Vector::Zero(3*nactive);
		Matrix temphessian = Matrix::Zero(3*nactive, 3*natoms);
		int row = 0; 
		for (int i : activeAtoms) {
			grad.segment(row, 3) = full_grad.segment(3*i, 3); 
			temphessian.block(row, 0, 3, 3*natoms) = full_hessian.block(3*i, 0, 3, 3*natoms); 
			row += 3; 
		}
	
		row = 0; 
		hessian = Matrix::Zero(3*nactive, 3*nactive); 
		for (int j : activeAtoms) {
			hessian.block(0, row, 3*nactive, 3) = temphessian.block(0, 3*j, 3*nactive, 3); 
			row += 3; 
		}
	
	}
	
	delta_e = energy - old_e; 
	grad_norm = grad.norm(); 
	old_e = energy; 
	iter++;
	
	return energy; 
	
}
