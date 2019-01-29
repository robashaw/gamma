#include "InternalCoords.hpp"
#include "mathutil.hpp"
#include <cmath>

/*
void Rotator::calcAxis() {
	// Assumes molecule is linear 
	int natoms = atoms.size();
	if (natoms > 1) {
		Vector vy = x0.row(atoms[natoms-1]) - x0.row(atoms[0]).normalized(); 
		
		double theta_x = acos(vy[0]);
		double theta_y = acos(vy[1]);
		double theta_z = acos(vy[2]); 
		
		Vector e0 = Vector::Zero(3); 
		if (theta_x < theta_y && theta_x < theta_z) e0[0] = 1.0;
		else if (theta_y < theta_x && theta_y < theta_z) e0[1] = 1.0;
		else e0[2] = 1.0; 
		
		axis = cross(vy, e0).normalized(); 
	}
}

Vector Rotator::value(const Matrix& xyz) {
	Vector rVal;
	if ((xyz - valXYZ).norm() < 1e-12) rVal = stored_value; 
	else {
		int natoms = atoms.size();

		Matrix subXYZ = Matrix::Zero(natoms, 3);
		Matrix subX0 = subXYZ;
		make_sub_mats(subXYZ, subX0, xyz); 
		
		rVal = get_exp_map(subXYZ, subX0);
		norm = rVal.norm(); 
		valXYZ = xyz; 
		stored_value = rVal;  
	}
	return rVal; 
}

Vector Rotator::make_sub_mats(Matrix& subXYZ, Matrix& subX0, const Matrix& xyz) {
	int natoms = atoms.size();
	
	Vector rVal = Vector::Zero(3); 
	
	for (int i = 0; i < natoms; i++) {
		subXYZ.row(i) = xyz.row(atoms[i]); 
		subX0.row(i) = x0.row(atoms[i]); 
	}
	
	Vector xmean(3), ymean(3);
	for (int i = 0; i < 3; i++) {
		xmean[i] = subXYZ.row(i).sum() / ((double) natoms);
		ymean[i] = subX0.row(i).sum() / ((double) natoms);
	}
	
	if (!linear && is_linear(subXYZ, subX0)) linear = true; 
	
	if (linear) {
		Vector vx = subXYZ.row(natoms-1) - subXYZ.row(0);
		Vector vy = subX0.row(natoms-1) - subX0.row(0); 
		
		calcAxis();
		
		Vector ev = vx.normalized();
		dot2 = ev.transpose() * axis; 
		dot2 *= dot2; 
		
		// Add dummy atom one Bohr from molecular center
		Vector xdum = cross(vx, axis).normalized();
		Vector ydum = cross(vy, axis).normalized();
		
		subXYZ.conservativeResize(subXYZ.rows() + 1, Eigen::NoChange); 
		subX0.conservativeResize(subX0.rows() + 1, Eigen::NoChange); 
		subXYZ.row(natoms) = xdum+xmean;
		subX0.row(natoms) = ydum+ymean; 
		
		xdum = cross(vx, axis);
		double dummy_norm = xdum.norm();
		Vector dxdum = d_cross(vx, axis); 
		Vector dnxdum = d_ncross(vx, axis); 
		rVal = dxdum * dummy_norm - dnxdum.outer(xdum) / (dummy_norm*dummy_norm); 
	}
	
	return rVal; 
}

std::vector<Matrix> Rotator::derivative(const Matrix& xyz) {
	int natoms = atoms.size();
	
	std::vector<Matrix> derivatives; 
	
	if ((xyz - derXYZ).norm() < 1e-12) derivatives = stored_deriv;
	else {
		Matrix subXYZ = Matrix::Zero(natoms, 3);
		Matrix subX0 = subXYZ;
		Vector dexdum = make_sub_mats(subXYZ, subX0, xyz); 
		
		std::vector<Matrix> deriv_raw = get_exp_map_der(subXYZ, subX0); 
		if (linear) {
			deriv_raw[0] -= dexdum * deriv_raw[natoms-1]; 
			for (int i = 0; i < natoms; i++) 
				deriv_raw[i] += Matrix::Identity(3, 3) * deriv_raw[natoms-1] / ((double) natoms);
			deriv_raw[natoms-2] += dexdum * deriv_raw[natoms-1]; 
		}
		
		int i = 0;
		for (int n = 0; n < xyz.rows(); n++) {
			if (n == atoms[i]) {
				derivatives.push_back(deriv_raw[i]);
				i++;
			} else {
				Matrix temp = Matrix::Zero(3, 3);
				derivatives.push_back(temp);
			}
		}
	}
	
	derXYZ = xyz;
	stored_deriv = derivatives; 
	
	return derivatives; 
}

double Angle::value(const Matrix& xyz) {
	Vector v1 = xyz.row(a) - xyz.row(b);
	Vector v2 = xyz.row(c) - xyz.row(b); 
	double n1 = v1.norm();
	double n2 = v2.norm();
	
	double cosx = v1.dot(v2) / (n1 * n2); 
	double value = 0.0;
	
	if (cosx + 1.0 < 1e-12)
		value = M_PI; 
	else if (cosx - 1.0 > 1e-12)
		value = 0.0; 
	else
		value = acos(cosx);  
	
	return value;
}

Vector Angle::normal_vector(const Matrix& xyz) {
	Vector v1 = xyz.row(a) - xyz.row(b);
	Vector v2 = xyz.row(c) - xyz.row(b); 
	return cross(v1, v2).normalized(); 
}

Matrix Angle::derivative(const Matrix& xyz) {
	Matrix derivs = Matrix::Zero(xyz.rows(), xyz.cols()); 
	
	Vector u = xyz.row(a) - xyz.row(b);
	Vector v = xyz.row(c) - xyz.row(b);
	double unorm = u.norm();
	double vnorm = v.norm();
	u /= unorm;
	v /= vnorm;

	double osq3 = 1.0 / sqrt(3.0); 
	Vector v1 = {  osq3, -osq3, osq3 };
	Vector v2 = { -osq3,  osq3, osq3 }; 
	
	Vector w; 
	if ( (u+v).norm() < 1e-10 || (u-v).norm() < 1e-10) {
		// Check if parallel
		if ((u+v1).norm() < 1e-10 || (u-v2).norm() < 1e-10) w = cross(u, v1).normalized();
		else w = cross(u, v2).normalized(); 
	} else {
		w = cross(u, v).normalized(); 
	}
	
	Vector term1 = cross(u, w) / unorm; 
	Vector term2 = cross(w, v) / vnorm; 
	
	derivs.row(a) = term1;
	derivs.row(b) = term2;
	derivs.row(c) = -(term1 + term2); 
	
	return derivs; 
}*/
