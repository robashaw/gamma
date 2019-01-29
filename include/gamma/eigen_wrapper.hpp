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

#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <Eigen/Cholesky>
#include <Eigen/SparseCholesky>
#include <Eigen/Eigenvalues>
#include "tensor4.hpp"

/*! \file eigen_wrapper
 	\brief Wrappers and definitions to make using the Eigen library easier. 
 */

using Matrix = Eigen::MatrixXd; 
using SparseMatrix = Eigen::SparseMatrix<double>; 
using Vector = Eigen::VectorXd; 
using iVector = Eigen::VectorXi; 
using EigenSolver = Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd>; 
using GeneralizedEigenSolver = Eigen::GeneralizedSelfAdjointEigenSolver<Eigen::MatrixXd>;
using LLT = Eigen::LLT<Eigen::MatrixXd>; 
using FullLU = Eigen::FullPivLU<Eigen::MatrixXd>;
using SparseLLT = Eigen::SimplicialLLT<Eigen::SparseMatrix<double>>;
using SparseLU = Eigen::SparseLU<Eigen::SparseMatrix<double>>; 

typedef Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> EMatrix;
typedef Eigen::DiagonalMatrix<double, Eigen::Dynamic, Eigen::Dynamic> DiagonalMatrix;

/*! Converts a Tensor4 object into an Eigen vector object. 
	The leading index is taken to be the last, i.e. the Tensor Tijkl will be unrolled as
	T0000, T0001, T0002, ..., T0010, T0011, ... 

	@param t - the Tensor4 object to be converted.
	@param symm - the symmetry class of the Tensor4 object.
	@return an Eigen::VectorXd object of the ravelled Tensor4 object.  
 */
Vector tensorToVector(const Tensor4& t, Tensor4::TYPE symm);

/*! Calculates the angle in radians (from u to v) between two Vector objects.
	@param u - the first Vector.
	@param v - the second Vector.
	@return The angle in radians from u to v. 
 */ 
double angle(const Vector& u, const Vector& w);

/*! Calculates the cross product between two Vector3d objects.
	@param u - the first Vector.
	@param v - the second Vector.
	@return The cross product vector, u x v.
 */
Vector cross(const Vector& u, const Vector& w);
