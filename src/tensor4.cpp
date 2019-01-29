/* 
 *
 *   PURPOSE: To implement class Tensor4
 *
 */

#include "tensor4.hpp"
#include <iostream>
#include <cmath>
#include <algorithm>
#include <cstdio>

Tensor4::Tensor4(int a, int b, int c, int d) : w(a), x(b), y(c), z(d)
{
  // Resize the data matrix
  data.resize(a*b*c*d);
}

Tensor4::Tensor4(int a, int b, int c, int d, double val) : w(a), x(b), y(c), z(d)
{
  // Resize and assign the data matrix
  data.resize(a*b*c*d);
  std::fill(data.begin(), data.end(), val);
}

Tensor4::Tensor4(const Tensor4& other)
{
  w = other.w; x = other.x; y = other.y; z = other.z;
  data.resize(w*x*y*z);
  for (int i = 0; i < w; i++){
    for (int j = 0; j < x; j++){
      for (int k = 0; k < y; k++){
        for (int l = 0; l < z; l++){
          data[i*x*y*z + j*y*z + k*z +l] = other(i, j, k, l);
        }
      } 
    }
  }
}  

int Tensor4::writeToFile(std::string filename) {
	int success = 1; 
	
	FILE* pFile = std::fopen(filename.c_str(), "wb"); 
    if (pFile == NULL) success = -1;
	else {
		std::fwrite(&data[0], sizeof(double), data.size(), pFile); 
		std::fclose(pFile);
		file_name = filename; 
		data_size = data.size();
		data.clear();
	}
	
	return success; 
}

int Tensor4::readFromFile() {
	int success = 1;

	FILE* pFile = std::fopen(file_name.c_str(), "rb"); 
	if (pFile == NULL) success = -1;
	else {
		double* buffer = (double*) malloc (sizeof(double)*data_size);
		size_t result = std::fread(buffer, sizeof(double), data_size, pFile);
		if (result != data_size) success = -2;
		else
			for (int i = 0; i < data_size; i++) data.push_back(buffer[i]); 
		free(buffer);
	}
	
	return success;
}

void Tensor4::resize(int a, int b, int c, int d)
{
  w = a; x = b; y = c; z = d;
  data.resize(a*b*c*d);
}

void Tensor4::assign(int a, int b, int c, int d, double val)
{
  w = a; x = b; y = c; z = d;
  data.resize(a*b*c*d);
  std::fill(data.begin(), data.end(), val);
}

void Tensor4::print() const 
{
  for (int i = 0; i < w; i++){
    for (int j = 0; j < x; j++){
      for (int k = 0; k < y; k++){
	for (int l = 0; l < z; l++){
	  std::cout << i << " " << j << " " << k << " " << l << "   " << data[i*x*y*z + j*y*z + k*z + l] << "\n";
	}
      }
    }
  }
  std::cout << "\n\n";
}

double& Tensor4::operator()(int i, int j, int k, int l)
{
  return data[i*x*y*z + j*y*z + k*z + l];
}

double Tensor4::operator()(int i, int j, int k, int l) const
{
  return data[i*x*y*z + j*y*z + k*z + l];
}

Tensor4& Tensor4::operator=(const Tensor4& other)
{
  resize(other.w, other.x, other.y, other.z);
  for (int i = 0; i < w; i++){
    for (int j = 0; j < x; j++){
      for (int k = 0; k < y; k++){
	for (int l = 0; l < z; l++){
	  data[i*x*y*z + j*y*z + k*z + l] = other(i, j, k, l);
	}
      }
    }
  }
  return *this;
}

Tensor4 Tensor4::operator+() const 
{
	Tensor4 retVal(w, x, y, z);
	for (int i = 0; i < w; i++){
		for (int j = 0; j < x; j++){
			for (int k = 0; k < y; k++){
				for (int l = 0; l < z; l++){
					retVal(i, j, k, l) = data[i*x*y*z + j*y*z + k*z + l];
				}
			}
		}
	}
	return retVal;
}

Tensor4 Tensor4::operator-() const 
{
	Tensor4 retVal(w, x, y, z);
	for (int i = 0; i < w; i++){
		for (int j = 0; j < x; j++){
			for (int k = 0; k < y; k++){
				for (int l = 0; l < z; l++){
					retVal(i, j, k, l) = -data[i*x*y*z + j*y*z + k*z + l];
				}
			}
		}
	}
	return retVal;
}

Tensor4 Tensor4::operator+(const Tensor4& other) const
{
	Tensor4 retVal(w, x, y, z);
	for (int i = 0; i < w; i++){
		for (int j = 0; j < x; j++){
			for (int k = 0; k < y; k++){
				for (int l = 0; l < z; l++){
					retVal(i, j, k, l) = data[i*x*y*z + j*y*z + k*z + l] + other(i, j, k, l);
				}
			}
		}
	}
	return retVal;
}

Tensor4 Tensor4::operator-(const Tensor4& other) const
{
	Tensor4 retVal(w, x, y, z);
	for (int i = 0; i < w; i++){
		for (int j = 0; j < x; j++){
			for (int k = 0; k < y; k++){
				for (int l = 0; l < z; l++){
					retVal(i, j, k, l) = data[i*x*y*z + j*y*z + k*z +l] - other(i, j, k, l);
				}
			}
		}
	}
	return retVal;
}

Tensor4 Tensor4::operator*=(const double& scalar) const {
	Tensor4 retVal(w, x, y, z);
	for (int i = 0; i < w; i++){
		for (int j = 0; j < x; j++){
			for (int k = 0; k < y; k++){
				for (int l = 0; l < z; l++){
					retVal(i, j, k, l) = data[i*x*y*z + j*y*z + k*z +l] * scalar;
				}
			}
		}
	}
	return retVal;
}

// Frobenius norm
double fnorm(const Tensor4& t)
{
	double retval = 0.0;
	for (int i = 0; i < t.w; i++)
		for (int j = 0; j < t.x; j++)
			for (int k = 0; k < t.y; k++)
				for (int l = 0; l < t.z; l++)
					retval += t(i, j, k, l) * t(i, j, k, l);
	
	return std::sqrt(retval);
}


S8EvenTensor4::S8EvenTensor4(int _w)
{
  // Resize the data matrix
  w = _w;
  resize(w);
}

S8EvenTensor4::S8EvenTensor4(int _w, double val)
{
  // Resize and assign the data matrix
	w = _w;
	assign(w, val);
}

S8EvenTensor4::S8EvenTensor4(const S8EvenTensor4& other)
{
  w = other.w;
	long long size = w*(w+1)*(w+2);
	size *= (3*w+1);
	size /= 24; 
  data.resize(size);
  updateMults();
  for (int i = 0; i < w; i++){
    for (int j = 0; j <= i; j++){
      for (int k = 0; k <= i; k++){
        for (int l = 0; l <= k; l++){
          data[imults[i] + j*jkmults[i+1] + jkmults[k] + l] = other(i, j, k, l);
        }
      } 
    }
  }
}  

void S8EvenTensor4::resize(int _w)
{
	w = _w;
  	long long size = w*(w+1)*(w+2);
	size *= (3*w+1);
	size /= 24; 
  	data.resize(size);
	updateMults();
}

void S8EvenTensor4::assign(int _w, double val)
{
	w = _w;
  	long long size = w*(w+1)*(w+2);
	size *= (3*w+1);
	size /= 24; 
  	data.resize(size);
  	std::fill(data.begin(), data.end(), val);
	updateMults();
}

void S8EvenTensor4::print() const 
{
  for (int i = 0; i < w; i++){
    for (int j = 0; j <= i; j++){
      for (int k = 0; k <= i; k++){
	for (int l = 0; l <= k; l++){
	  std::cout << i << " " << j << " " << k << " " << l << "   " << data[imults[i] + j*jkmults[i+1] + jkmults[k] + l] << "\n";
	}
      }
    }
  }
  std::cout << "\n\n";
}

double& S8EvenTensor4::operator()(int i, int j, int k, int l)
{
	int I = i > j ? i : j;
	int J = i > j ? j : i;
	int K = k > l ? k : l;
	int L = k > l ? l : k;
	if (K > I) {
		int tmp = I;
		I = K;
		K = tmp;
		tmp = J;
		J = L;
		L = tmp;
	}
  	return data[imults[I] + J*jkmults[I+1] + jkmults[K] + L];
}

double S8EvenTensor4::operator()(int i, int j, int k, int l) const
{
	int I = i > j ? i : j;
	int J = i > j ? j : i;
	int K = k > l ? k : l;
	int L = k > l ? l : k;
	if (K > I) {
		int tmp = I;
		I = K;
		K = tmp;
		tmp = J;
		J = L;
		L = tmp;
	}
  	return data[imults[I] + J*jkmults[I+1] + jkmults[K] + L];
}

S8EvenTensor4& S8EvenTensor4::operator=(const S8EvenTensor4& other)
{
  resize(other.getW());
  for (int i = 0; i < w; i++){
    for (int j = 0; j <= i; j++){
      for (int k = 0; k <= i; k++){
	for (int l = 0; l <= k; l++){
	  data[imults[i] + j*jkmults[i+1] + jkmults[k] + l] = other(i, j, k, l);
	}
      }
    }
  }
  return *this;
}

S8EvenTensor4 S8EvenTensor4::operator+() const
{
	S8EvenTensor4 retVal(w);
	for (int i = 0; i < w; i++){
		for (int j = 0; j <= i; j++){
			for (int k = 0; k <= i; k++){
				for (int l = 0; l <= k; l++){
					retVal(i, j, k, l) = data[imults[i] + j*jkmults[i+1] + jkmults[k] + l] ;
				}
			}
		}
	}
	return retVal;
}

S8EvenTensor4 S8EvenTensor4::operator-() const
{
	S8EvenTensor4 retVal(w);
	for (int i = 0; i < w; i++){
		for (int j = 0; j <= i; j++){
			for (int k = 0; k <= i; k++){
				for (int l = 0; l <= k; l++){
					retVal(i, j, k, l) = -data[imults[i] + j*jkmults[i+1] + jkmults[k] + l] ;
				}
			}
		}
	}
	return retVal;
}

S8EvenTensor4 S8EvenTensor4::operator+(const S8EvenTensor4& other) const
{
	S8EvenTensor4 retVal(w);
	for (int i = 0; i < w; i++){
		for (int j = 0; j <= i; j++){
			for (int k = 0; k <= i; k++){
				for (int l = 0; l <= k; l++){
					retVal(i, j, k, l) = data[imults[i] + j*jkmults[i+1] + jkmults[k] + l] + other(i, j, k, l);
				}
			}
		}
	}
	return retVal;
}

S8EvenTensor4 S8EvenTensor4::operator-(const S8EvenTensor4& other) const
{
	S8EvenTensor4 retVal(w);
	for (int i = 0; i < w; i++){
		for (int j = 0; j <= i; j++){
			for (int k = 0; k <= i; k++){
				for (int l = 0; l <= k; l++){
					retVal(i, j, k, l) = data[imults[i] + j*jkmults[i+1] + jkmults[k] + l] - other(i, j, k, l);
				}
			}
		}
	}
	return retVal;
}

S8EvenTensor4 S8EvenTensor4::operator*=(const double& scalar) const
{
	S8EvenTensor4 retVal(w);
	for (int i = 0; i < w; i++){
		for (int j = 0; j <= i; j++){
			for (int k = 0; k <= i; k++){
				for (int l = 0; l <= k; l++){
					retVal(i, j, k, l) = data[imults[i] + j*jkmults[i+1] + jkmults[k] + l] * scalar;
				}
			}
		}
	}
	return retVal;
}

void S8EvenTensor4::updateMults() {
	
	imults.clear(); jkmults.clear();
	for (int i = 0; i < w; i++) {
		int val = ( i * (i+1) ) / 2;
		jkmults.push_back(val);
		val *= (i+2) * (3*i + 1);
		imults.push_back(val/12);
	}
	jkmults.push_back( (w * (w+1))/2 );
}

// Frobenius norm
double fnorm(const S8EvenTensor4& t)
{
	double retval = 0.0;
	for (int i = 0; i < t.w; i++)
		for (int j = 0; j <= i; j++)
			for (int k = 0; k <= i; k++)
				for (int l = 0; l <= k; l++)
					retval += t(i, j, k, l) * t(i, j, k, l);
	
	return std::sqrt(8*retval);
}

S8OddTensor4::S8OddTensor4(int N) : S8EvenTensor4(N) { }
S8OddTensor4::S8OddTensor4(int N, double val) : S8EvenTensor4(N, val) { }
S8OddTensor4::S8OddTensor4(const S8OddTensor4& other) : S8EvenTensor4(other) { } 

void S8OddTensor4::set(int i, int j, int k, int l, double val) {
	int parity = 0;
	int I = i;
	int J = j;
	int K = k;
	int L = l;
	if (j > i) { I = j; J = i; parity++; }
	if (l > k) { K = l; L = k; parity++; }
	if (K > I) {
		int tmp = I;
		I = K;
		K = tmp;
		tmp = J;
		J = L;
		L = tmp;
	}
	data[imults[I] + J*jkmults[I+1] + jkmults[K] + L] = (1 - 2*(parity%2)) * val;
}

double S8OddTensor4::operator()(int i, int j, int k, int l) const {
	int parity = 0;
	int I = i;
	int J = j;
	int K = k;
	int L = l;
	if (j > i) { I = j; J = i; parity++; }
	if (l > k) { K = l; L = k; parity++; }
	if (K > I) {
		int tmp = I;
		I = K;
		K = tmp;
		tmp = J;
		J = L;
		L = tmp;
	}
  	return (1 - 2*(parity%2)) * data[imults[I] + J*jkmults[I+1] + jkmults[K] + L];
}

S4OddTensor4::S4OddTensor4(int _w, int _x)
{
  // Resize the data matrix
  resize(_w, _x);
}

S4OddTensor4::S4OddTensor4(int _w, int _x, double val)
{
  // Resize and assign the data matrix
	assign(_w, _x, val);
}

S4OddTensor4::S4OddTensor4(const S4OddTensor4& other)
{
  resize(other.w, other.x);
  for (int i = 0; i < w; i++){
    for (int j = 0; j <= i; j++){
      for (int k = 0; k < x; k++){
        for (int l = 0; l <= k; l++){
          data[(mults[i] + j)*mults[x] + mults[k] + l] = other(i, j, k, l);
        }
      } 
    }
  }
}  

void S4OddTensor4::resize(int _w, int _x)
{
	w = _w;
	x = _x;
	int size = (w * x * (w+1) * (x+1)) / 4;
  	data.resize(size);
	updateMults();
}

void S4OddTensor4::assign(int _w, int _x, double val)
{
	w = _w;
	x = _x;
  	int size = (w * x * (w+1) * (x+1)) / 4;
  	data.resize(size);
  	std::fill(data.begin(), data.end(), val);
	updateMults();
}

void S4OddTensor4::print() const 
{
  for (int i = 0; i < w; i++){
    for (int j = 0; j <= i; j++){
      for (int k = 0; k < x; k++){
	for (int l = 0; l <= k; l++){
	  std::cout << i << " " << j << " " << k << " " << l << "   " << data[(mults[i] + j)*mults[x] + mults[k] + l] << "\n";
	}
      }
    }
  }
  std::cout << "\n\n";
}

void S4OddTensor4::set(int i, int j, int k, int l, double value) {
	int parity = 0;
	int tmp = i;
	if ( j > i ) {
		i = j;
		j = tmp;
		parity++;
	} 
	tmp = k;
	if ( l > k ) {
		k = l;
		l = tmp;
		parity++;
	}
	data[(mults[i] + j)*mults[x] + mults[k] + l] = (1- 2*(parity%2)) * value;
}

double S4OddTensor4::operator()(int i, int j, int k, int l) const
{
	int parity = 0;
	int tmp = i;
	if ( j > i ) {
		i = j;
		j = tmp;
		parity++;
	} 
	tmp = k;
	if ( l > k ) {
		k = l;
		l = tmp;
		parity++;
	}
			
  	return (1 - 2*(parity%2)) * data[(mults[i] + j)*mults[x] + mults[k] + l];
}

S4OddTensor4& S4OddTensor4::operator=(const S4OddTensor4& other)
{
  resize(other.getW(), other.getX());
  for (int i = 0; i < w; i++){
    for (int j = 0; j <= i; j++){
      for (int k = 0; k < x; k++){
	for (int l = 0; l <= k; l++){
	  data[(mults[i] + j)*mults[x] + mults[k] + l] = other(i, j, k, l);
	}
      }
    }
  }
  return *this;
}

S4OddTensor4 S4OddTensor4::operator+() const
{
	S4OddTensor4 retVal(w, x);
	for (int i = 0; i < w; i++){
		for (int j = 0; j <= i; j++){
			for (int k = 0; k < x; k++){
				for (int l = 0; l <= k; l++){
					retVal.set(i, j, k, l, data[(mults[i] + j)*mults[x] + mults[k] + l] );
				}
			}
		}
	}
	return retVal;
}

S4OddTensor4 S4OddTensor4::operator-() const
{
	S4OddTensor4 retVal(w, x);
	for (int i = 0; i < w; i++){
		for (int j = 0; j <= i; j++){
			for (int k = 0; k < x; k++){
				for (int l = 0; l <= k; l++){
					retVal.set(i, j, k, l, -data[(mults[i] + j)*mults[x] + mults[k] + l] );
				}
			}
		}
	}
	return retVal;
}

S4OddTensor4 S4OddTensor4::operator+(const S4OddTensor4& other) const
{
	S4OddTensor4 retVal(w, x);
	for (int i = 0; i < w; i++){
		for (int j = 0; j <= i; j++){
			for (int k = 0; k < x; k++){
				for (int l = 0; l <= k; l++){
					retVal.set(i, j, k, l, data[(mults[i] + j)*mults[x] + mults[k] + l] + other(i, j, k, l) );
				}
			}
		}
	}
	return retVal;
}

S4OddTensor4 S4OddTensor4::operator-(const S4OddTensor4& other) const
{
	S4OddTensor4 retVal(w, x);
	for (int i = 0; i < w; i++){
		for (int j = 0; j <= i; j++){
			for (int k = 0; k < x; k++){
				for (int l = 0; l <= k; l++){
					retVal.set(i, j, k, l, data[(mults[i] + j)*mults[x] + mults[k] + l] - other(i, j, k, l) );
				}
			}
		}
	}
	return retVal;
}

S4OddTensor4 S4OddTensor4::operator*=(const double& scalar) const
{
	S4OddTensor4 retVal(w, x);
	for (int i = 0; i < w; i++){
		for (int j = 0; j <= i; j++){
			for (int k = 0; k < x; k++){
				for (int l = 0; l <= k; l++){
					retVal.set(i, j, k, l, data[(mults[i] + j)*mults[x] + mults[k] + l] * scalar);
				}
			}
		}
	}
	return retVal;
}

void S4OddTensor4::updateMults() {
	
	mults.clear();
	auto max_i = w > x ? w : x;
	for (int i = 0; i <= max_i; i++) {
		int val = ( i * (i+1) ) / 2;
		mults.push_back(val);
	}
	
}

// Frobenius norm
double fnorm(const S4OddTensor4& t)
{
	double retval = 0.0;
	for (int i = 0; i < t.w; i++)
		for (int j = 0; j <= i; j++)
			for (int k = 0; k < t.x; k++)
				for (int l = 0; l <= k; l++)
					retval += t(i, j, k, l) * t(i, j, k, l);
	
	return 2*std::sqrt(retval);
}

S4EvenTensor4::S4EvenTensor4(int n, int m) : S4OddTensor4(n, m) { }
S4EvenTensor4::S4EvenTensor4(int n, int m, double val) : S4OddTensor4(n, m, val) { }
S4EvenTensor4::S4EvenTensor4(const S4EvenTensor4& other) : S4OddTensor4(other) { } 

void S4EvenTensor4::set(int i, int j, int k, int l, double val) {
	int I = i > j ? i : j;
	int J = i > j ? j : i;
	int K = k > l ? k : l;
	int L = k > l ? l : k;
	data[(mults[i] + j)*mults[x] + mults[k] + l] = val;
}

double S4EvenTensor4::operator()(int i, int j, int k, int l) const {
	int I = i > j ? i : j;
	int J = i > j ? j : i;
	int K = k > l ? k : l;
	int L = k > l ? l : k;
	return data[(mults[i] + j)*mults[x] + mults[k] + l];
}
