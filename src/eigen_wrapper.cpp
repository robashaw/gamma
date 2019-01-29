#include "eigen_wrapper.hpp"
#include <vector>

Vector tensorToVector(const Tensor4& t, Tensor4::TYPE symm) {
	Vector v;
	if (symm == Tensor4::ODD_8 || symm == Tensor4::EVEN_8) {
		// 8-fold symmetry
		int w = t.getW();
		int n = (w*(w+1)*(w+2)*(3*w+1)) / 24;
		v = Vector::Zero(n);
		int ix = 0;
		for (int i = 0; i < w; i++)
			for (int j = 0; j <= i; j++)
				for (int k = 0; k <= i; k++)
					for (int l = 0; l <= k; l++)
						v[ix++] = t(i, j, k, l);
	} else if (symm == Tensor4::ODD_4 || symm == Tensor4::EVEN_4) {
		// 4-fold symmetry
		int w = t.getW(); int x = t.getX();
		int n = (w * x * (w+1) * (x+1)) / 4; 
		v = Vector::Zero(n);
		int ix = 0;
		for (int i = 0; i < w; i++)
			for (int j = 0; j <= i; j++)
				for (int k = 0; k < x; k++)
					for (int l = 0; l <= k; l++)
						v[ix++] = t(i, j, k, l);
	} else {
		// No symmetries
		int w = t.getW(); int x = t.getX();
		int y = t.getY(); int z = t.getZ();
		int n = w * x * y * z;
		v = Vector::Zero(n);
		int ix = 0;
		for (int i = 0; i < w; i++)
			for (int j = 0; j < x; j++)
				for (int k = 0; k < y; k++)
					for (int l = 0; l < z; l++)
						v[ix++] = t(i, j, k, l);
	}
	
	return v;
}

// Get the angle between two vectors
double angle(const Vector& u, const Vector& w)
{
  // Get the magnitudes of the vectors
  double unorm = u.norm();
  double wnorm = w.norm();
  // Get the dot product
  double dprod = u.dot(w);
  // Use the cosine rule
  // but make sure neither is a zero vector
  double rval = 0.0;
  if(dprod > 1E-12){
    rval = std::acos(dprod/(unorm*wnorm));
  }
  return rval;
}

Vector cross(const Vector& u, const Vector& w)
{
    Vector r(3);
    r[0] = u(1)*w(2) - w(1)*u(2);
    r[1] = w(0)*u(2) - w(2)*u(0);
    r[2] = u(0)*w(1) - u(1)*w(0);
    return r;
}
