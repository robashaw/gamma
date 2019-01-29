/* 
 *    PURPOSE: To implement mathutil.hpp, a library of maths functions needed
 *             by various other programs.
 *
 *    DATE         AUTHOR          CHANGES
 *    =======================================================================
 *    27/08/15     Robert Shaw     Original Code
 *    28/08/15     Robert Shaw     Implemented Boys function.
 *    02/09/15     Robert Shaw     Implemented binom function.
 *    05/09/15     Robert Shaw     Implemented clebsch, sphernorm.
 *    07/09/15     Robert Shaw     Implemented formTransMat.
 *    12/09/15     Robert Shaw     Move vecxmatrix multiplication here.
 */

#include <cmath>
#include "mathutil.hpp"
#include "error.hpp"
#include <algorithm> 

// Factorial and double factorial functions
unsigned long int fact(int i)
{
  unsigned long int rval = (i<1 ? 1 : i);
  for(int j = i-1; j > 0; j--){
    rval *= j;
  }
  return rval;
}

unsigned long int fact2(int i)
{
  unsigned long int rval = (i<1 ? 1 : i);
  for (int j = i-2; j > 0; j-=2){
    rval = rval*j;
  }
  return rval;
}

void factArray(int i, double *values) {
  values[0] = 1.0;
  if ( i > 0 ) {
    values[1] = 1.0;
    for (int j = 2; j <= i; j++) values[j] = values[j-1]*j;
  }				
}

void fact2Array(int i, double *values) {
  values[0] = 1.0;
  if ( i > 0 ) {
    values[1] = 1.0;
    for (int j = 2; j <= i; j++) values[j] = values[j-2]*j;
  }
}

// Calculate the binomial coefficient (n m)T
// given by the formula:
// n! / [(n-m)!*m!]
unsigned long int binom(int n, int m)
{  
  unsigned long int rval;
  if (m >= n || m <= 0 || n <= 0){
    rval = 1;
  } else {
    unsigned long int numerator = n;
    // Set m to be the lesser of m, n-m
    m = (m < n-m ? m : n-m);
    if ( m == 1 ) {
      rval = n;
    } else {
      unsigned long int denominator = m;
      for (int i = m-1; i > 0; i--){
	numerator *= (n-i);
	denominator *= i;
      }
      rval =  numerator/denominator;
    }
  }
  return rval;
}

// Calculate the Clebsch-Gordan Coefficient C^{lm}_{tuv}
// according to the formula in Helgaker, Jorgensen, Olsen,
// Molecular Electronic Structure Theory, chapter 9 pg 338
double clebsch(int l, int m, int t, int u, double v)
{
  double vm = (m < 0 ? 0.5 : 0.0);
  int sign = ( ((int)(t+v-vm))%2 == 0 ? 1 : -1);
  double cval = sign*std::pow(0.25, t);

  // Get the binomial factors
  int mabs = std::abs(m);
  double ltb = binom(l, t);
  double lmtb = binom(l-t, mabs + t);
  double tub = binom(t, u);
  cval = cval*ltb*lmtb*tub*binom(mabs, (int)(2*v));
  return cval;
}

// Calculate the spherical normalisation Nlm
double sphernorm(int l, int m)
{
  int mabs = std::abs(m);
  double nval = std::pow(0.5, mabs)/((double)(fact(l)));
  
  // Get factorials
  double lpmfact = fact(l + mabs);
  double lmmfact = fact(l - mabs);
  double zerom = (m==0 ? 2.0 : 1.0);
  nval = nval*std::sqrt(2*lpmfact*lmmfact/zerom);
  return nval;
}

// Fill in a portion of the cart-to-spher. transformation matrix
// corresponding to quantum numbers l, m. Assumes canonical order.
void formTransMat(Matrix& mat, int row, int col, int l, int m){
  double temp; 
  switch(l){
  case 1: { // p-type, one of 3 identity 3-vectors
    switch(m){
    case -1: { mat(row, col+1) = 1.0; break; } // py
    case 0: { mat(row, col+2) = 1.0; break; }  // pz
    default: { mat(row, col) = 1.0;} // px
    }
    break;
  }
  case 2: { // d-type
    switch(m) {
    case -2: { mat(row, col+4) = 1.0; break; } // dxy
    case -1: { mat(row, col+1) = 1.0; break; } // dyz
    case 0: { // dz2 
      mat(row, col) = 1.0; 
      mat(row, col+3) = mat(row, col+5) = -0.5;
      break;
    }
    case 1: { mat(row, col+2) = 1.0; break; } // dxz
    default: { // d(x^2-y^2)
      temp = std::sqrt(3.0/4.0);
      mat(row, col+5) = temp;
      mat(row, col+3) = -temp;
    }
    }
    break;
  }
  case 3: { // f-type
    switch(m){
    case -3: { 
      mat(row, col+6) = std::sqrt(10)/4.0;
      mat(row, col+8) = -3.0*std::sqrt(2)/4.0;
      break;
    }
    case -2: { mat(row, col+4) = 1.0; break; }
    case -1: { 
      temp = std::sqrt(6.0/5.0);
      mat(row, col+1) = temp;
      mat(row, col+6) = -std::sqrt(6)/4.0;
      mat(row, col+8) = -temp/4.0;
      break;
    }
    case 0: {
      mat(row, col) = mat(row, col+3) = 1.0;
      mat(row, col+5) = -3.0/(2.0*std::sqrt(5));
      break;
    }
    case 1: {
      temp = std::sqrt(6.0/5.0);
      mat(row, col+2) = temp;
      mat(row, col+9) = -std::sqrt(6)/4.0;
      mat(row, col+7) = -temp/4.0;
      break;
    }
    case 2: {
      temp = std::sqrt(3.0/4.0);
      mat(row, col+5) = temp;
      mat(row, col+3) = -temp;
      break;
    }
    default: {
      mat(row, col+9) = std::sqrt(10)/4.0;
      mat(row, col+7) = -3.0*std::sqrt(2)/4.0;
    }
    }
    break;
  }
  case 4: { // g-type
    switch(m) {
    case -4:{
      temp = std::sqrt(5.0/4.0);
      mat(row, col+13) = temp;
      mat(row, col+11) = -temp;
      break;
    }
    case -3:{
      mat(row, col+6) = -1.0*std::sqrt(10.0)/4.0;
      mat(row, col+8) = 0.75*std::sqrt(2.0);
      break;
    }
    case -2:{
      temp = -1.0*std::sqrt(5.0/7.0)/2.0;
      mat(row, col+4) = 3.0*std::sqrt(1.0/7.0);
      mat(row, col+13) = temp;
      mat(row, col+11) = temp;
      break;
    }
    case -1:{
      temp = std::sqrt(2.0/7.0);
      mat(row, col+1) = temp*std::sqrt(5.0);
      mat(row, col+6) = -0.75*temp*std::sqrt(5.0);
      mat(row, col+8) = -0.75*temp;
      break;
    }
    case 0:{
      temp = -3.0*std::sqrt(3.0/35.0);
      mat(row, col) = 1.0;
      mat(row, col+14) = mat(row, col+10) = 3.0/8.0;
      mat(row, col+5) = mat(row, col+3) = temp;
      mat(row, col+12) = -0.25*temp;
      break;
    }
    case 1:{ 
      temp = std::sqrt(2.0/7.0);
      mat(row, col+2) = std::sqrt(5.0)*temp;
      mat(row, col+9) = -0.75*std::sqrt(5.0)*temp;
      mat(row, col+7) = -0.75*temp;
      break;
    }
    case 2:{
      temp = 1.5*std::sqrt(3.0/7.0);
      mat(row, col+5) = temp;
      mat(row, col+3) = -temp;
      temp = 0.25*std::sqrt(5.0);
      mat(row, col+14) = -temp;
      mat(row, col+10) = temp;
      break;
    }
    case 3:{
      temp = 0.25*std::sqrt(2.0);
      mat(row, col+9) = temp*std::sqrt(5.0);
      mat(row, col+7) = -3.0*temp;
      break;
    }
    default:{
      temp = std::sqrt(35.0)/8.0;
      mat(row, col+14) = mat(row, col+10) = temp;
      mat(row, col+12) = -0.75*std::sqrt(3.0);
    }     
    }
    break;
  }
  default: { // s-type
    mat(row, col) = 1.0;
  }
  }
}


// Vector x matrix and matrix x vector- will throw an error if wrong shapes                                                                                         
// Left and right vector x matrix multiplication functions

Vector lmultiply(const Vector& v, const Matrix& mat)
{
  // Assume left multiplication implies transpose                                                                                                                   
  int n = v.size();
  int cols = mat.cols();
  int rows = mat.rows();
  Vector rVec(cols); // Return vector should have dimension cols                                                                                                    
  // For this to work we require n = rows                                                                                                                           
  if (n != rows) {
    throw(Error("VECMATMULT", "Vector and matrix are wrong sizes to multiply."));
  } else { // Do the multiplication                                                                                                                                 
    for (int i = 0; i < cols; i++){
      Vector temp = mat.col(i);
      rVec[i] = v.dot(temp);
    }
  }
  return rVec;
}


Vector rmultiply(const Matrix& mat, const Vector& v)
{
  int n = v.size();
  int cols = mat.cols();
  int rows = mat.rows();
  Vector rVec(rows); // Return vector should have dimension rows                                                                                                    
  // For this to work we require n = cols                                                                                                                           
  if (n != cols) {
    throw(Error("MATVECMULT", "Vector and matrix are wrong sizes to multiply."));
  } else { // Do the multiplication                                                                                                                                 
    for (int i = 0; i < rows; i++){
      Vector temp = mat.row(i);
      rVec[i] = v.dot(temp);
    }
  }
  return rVec;
}

/*
Vector cross(const Vector& v1, const Vector& v2) {
	Eigen::Vector3d u1, u2;
	u1[0] = v1[0]; u2[0] = v2[0];
	u1[1] = v1[1]; u2[1] = v2[1];
	u1[2] = v1[2]; u2[2] = v2[2];
	return u1.cross(u2); 
}

Matrix d_cross(const Vector& v1, const Vector& v2) {
	Matrix dc = Matrix::Zero(3, 3); 
	for (int i = 0; i < 3; i++) {
		Vector ei = Vector::Zero(3);
		ei[i] = 1.0;
		dc.row(i) = cross(ei, v2); 
	}
	return dc;
}

Matrix d_cross_ab(const Vector& v1, const Vector& v2, const Matrix& da, const Matrix& db) {
	int nrow = da.rows(); 
	Matrix dnc = Matrix::Zero(nrow, 3); 
	for (int i = 0; i < nrow; i++)
		dnc.row(i) = cross(v1, db.row(i)) + cross(da.row(i), v2); 
	return dnc;
}

double ncross(const Vector& v1, const Vector& v2) {
	return cross(v1, v2).norm(); 
}

Vector d_ncross(const Vector& v1, const Vector& v2) {
	double nv1v2 = cross(v1, v2).norm();
	Vector result = v1 * v2.dot(v2) - v2 * v1.dot(v2); 
	return result / nv1v2; 
}

Matrix pseudo_inverse(Matrix& mat, double threshold) {
	Eigen::JacobiSVD<Matrix> svd(mat, Eigen::ComputeThinU | Eigen::ComputeThinV); 
	double tolerance = threshold * std::max(mat.cols(), mat.rows()) * svd.singularValues().array().abs()(0); 
	return svd.matrixV() *  (svd.singularValues().array().abs() > tolerance).select(svd.singularValues().array().inverse(), 0).matrix().asDiagonal() * svd.matrixU().adjoint();
}

// Build the F-matrix for constructing the rotation quaternion
Matrix build_F(const Matrix& x, const Matrix& y) {
	Matrix R = x.transpose() * y; 
	Matrix F = Matrix::Zero(4, 4); 
	
	F(0, 0) = R.trace(); 
	F(0, 1) = F(1, 0) = R(1, 2) - R(2, 1); 
	F(0, 2) = F(2, 0) = R(2, 0) - R(0, 2);
	F(0, 3) = F(3, 0) = R(0, 1) - R(1, 0);
	F(1, 1) = R(0, 0) - R(1, 1) - R(2, 2);
	F(1, 2) = F(2, 1) = R(0, 1) + R(1, 0); 
	F(1, 3) = F(3, 1) = R(0, 2) + R(2, 0);
	F(2, 2) = R(1, 1) - R(2, 2) - R(0, 0);
	F(2, 3) = F(3, 2) = R(1, 2) + R(2, 1);
	F(3, 3) = R(2, 2) - R(0, 0) - R(1, 1); 
	
	return F; 
}

bool is_linear(const Matrix& xyz, const Matrix& x0) {
	bool linear = false; 
	
	int natoms = xyz.rows(); 
	
	Vector xmean(3), ymean(3);
	for (int i = 0; i < 3; i++) {
		xmean[i] = xyz.col(i).sum() / ((double) natoms);
		ymean[i] = x0.col(i).sum() / ((double) natoms); 
	}
	
	Matrix x = xyz;
	Matrix y = x0; 
	for (int i = 0; i < natoms; i++) {
		x.row(i) -= xmean;
		y.row(i) -= ymean; 
	}
	
	Matrix F = build_F(x, y);
	EigenSolver es(F);
	Vector& L = es.eigenvalues();
	
	if (L[0]/L[1] < 1.01 && L[0]/L[1] > 0.0) linear = true;
	return linear;
}

// Quaternion that rotates into y, minimizing RMSD
Vector get_quat(const Matrix& x, const Matrix& y) {
	int natoms = xyz.rows(); 
	
	Vector xmean(3), ymean(3);
	for (int i = 0; i < 3; i++) {
		xmean[i] = xyz.col(i).sum() / ((double) natoms);
		ymean[i] = x0.col(i).sum() / ((double) natoms); 
	}
	
	Matrix x = xyz;
	Matrix y = x0; 
	for (int i = 0; i < natoms; i++) {
		x.row(i) -= xmean;
		y.row(i) -= ymean; 
	}
	
	Matrix F = build_F(x, y); 
	EigenSolver es(F);
	
	Vector q = es.eigenvectors().col(0); 
	if (q[0] < 0) q *= -1.0;
	
	return q;
}

Tensor4 get_R_der(const Matrix& x, const Matrix& y) {
	int natoms = x.rows();
	
	Tensor4 ADiffR(natoms, 3, 3, 3, 0.0); 
	for (int u = 0; u < natoms; u++)
		for (int w = 0; w < 3 w++) 
			for (int j = 0; j < 3; j++) ADiffR(u, w, w, j) = y(u, j); 
	
	return ADiffR; 
}

Tensor4 get_F_der(const Matrix& x, const Matrix& y) {
	int natoms = x.rows();
	
	Tensor4 dR = get_R_der(x, y); 
	Tensor4 dF(natoms, 3, 4, 4, 0.0); 
	
	double dr11, dr12, dr13, dr21, dr22, dr23, dr31, dr32, dr33; 
	for (int u = 0; u < natoms; u++) {
		for (int w = 0; w < 3; w++) {
			dr11 = dR(u, w, 0, 0); 
			dr12 = dR(u, w, 0, 1);
			dr13 = dR(u, w, 0, 2);
			dr21 = dR(u, w, 1, 0);
			dr22 = dR(u, w, 1, 1);
			dr23 = dR(u, w, 1, 2);
			dr31 = dR(u, w, 2, 0);
			dr32 = dR(u, w, 2, 1);
			dr33 = dR(u, w, 2, 2);
			
			dF(u, w, 0, 0) = dr11 + dr22 + dr33;
			dF(u, w, 0, 1) = dF(u, w, 1, 0) = dr23 - dr32; 
			dF(u, w, 0, 2) = dF(u, w, 2, 0) = dr31 - dr13;
			dF(u, w, 0, 3) = dF(u, w, 3, 0) = dr12 - dr21;
			dF(u, w, 1, 1) = dr11 - dr22 - dr33;
			dF(u, w, 1, 2) = dF(u, w, 2, 1) = dr12 + dr21;
			dF(u, w, 1, 3) = dF(u, w, 3, 1) = dr13 + dr31; 
			dF(u, w, 2, 2) = dr22 - dr33 - dr11;
			dF(u, w, 2, 3) = dF(u, w, 3, 2) = dr23 + dr32; 
			dF(u, w, 3, 3) = dr33 - dr11 - dr22; 
		}
	}
	return dF; 
}

std::vector<Matrix> get_q_der(const Matrix& x, const Matrix& y) {
	int natoms = xyz.rows(); 
	
	Vector xmean(3), ymean(3);
	for (int i = 0; i < 3; i++) {
		xmean[i] = xyz.col(i).sum() / ((double) natoms);
		ymean[i] = x0.col(i).sum() / ((double) natoms); 
	}
	
	Matrix x = xyz;
	Matrix y = x0; 
	for (int i = 0; i < natoms; i++) {
		x.row(i) -= xmean;
		y.row(i) -= ymean; 
	}
	
	Matrix F = build_F(x, y); 
	Tensor4 dF = get_F_der(x, y); 
	
	EigenSolver es(F);
	
	Vector q = es.eigenvectors().col(0); 
	if (q[0] < 0) q *= -1.0;
	double l = es.eigenvalues()[0]; 
	
	Matrix mat = Matrix::Identity(4, 4) * l - F; 
	Matrix pinv = pseudo_inverse(mat, 1e-6); 
	
	std::vector<Matrix> qder;
	for (int u = 0; u < natoms; u++) {
		Matrix temp1(3, 4);
		for (int w = 0; w < 3; w++) {
			Matrix temp2(4, 4);
			for (int i = 0; i < 4)
				for (int j = 0; j < 4) temp2(i, j) = dF(u, w, i, j); 
			
			temp1.row(w) = pinv * temp2 * q; 
		}
		qder.push_back(temp1);
	}
	
	return qder; 
}

Vector get_exp_map(const Matrix& xyz, const Matrix& x0) {
	Vector q = get_quat(xyz, x0); 
	double q0 = q[0] - 1.0;
	Vector v(3) = { q[1], q[2], q[3] }; 
	if (fabs(q0) < 1e-8) v*=(2.0 - 2.0 * q0/3.0); 
	else {
		q0 += 1.0;
		double f = 2 * acos(q0) / sqrt(1.0 - q0 * q0); 
		v *= f; 
	} 
	return v;
}

std::vector<Matrix> get_exp_map_der(const Matrix& xyz, const Matrix& x0) {
	int natoms = xyz.rows(); 
	
	Vector q1 = get_quat(xyz, x0); 
	Vector q(3) = { q1[1], q1[2], q1[3] }; 
	double q0 = q1[0] - 1.0; 
	
	Vector v; 
	double f, df; 
	if (fabs(q0) < 1e-8) {
		f = 2.0 - 2.0 * q0/3.0;
		df = -2.0/3.0; 
	} else {
		q0 += 1.0;
		double aq = acos(q0);
		double oq2 = sqrt(1.0 - q0 * q0);
		f = 2 * aq / oq2; 
		df = -2.0 / (oq2 * oq2); 
		df += 2.0 * q0 * aq / (oq2 * oq2 * oq2); 
	} 
	v = q * f; 
	
	Matrix dvdq = Matrix::Zero(4, 3); 
	dvdq.row(0) = df * q; 
	for (int i = 0; i < 3; i++) dvdq(i+1, i) = f; 
	
	std::vector<Matrix> dqdx = get_q_der(xyz, x0); 
	std::vector<Matrix> dvdx; 
	for (int i = 0; i < natoms; i++) {
		Matrix temp = Matrix::Zero(3, 3); 
		for (int w = 0; w < 3; w++) {
			Vector dqdx_iw = dqdx[i].row(w); 
			for (int j = 0; j < 4; j++) temp.row(w) += dvdq.row(j) * dqdx[i](w, j); 
		}
		dvdx.push_back(temp); 
	}  
	
	return dvdx; 
}
*/
