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

#include "eigen_wrapper.hpp"
#include <fstream>
#include <iostream>

/*! \file eigen_util
	\brief Utility functions for dealing with Eigen objects. 
 */
namespace eigenutil {
	
	/*! Writes any Eigen matrix as a binary stream to file.
		@tparam Matrix - the type of Eigen matrix to be written to file.
		@param filename - the name of the binary file to be created.
		@param matrix - the matrix to be written to file.
	 */
	template<class Matrix>
	void write_binary(const char* filename, const Matrix& matrix){
		std::ofstream out(filename, std::ios::out | std::ios::binary | std::ios::trunc);
		typename Matrix::Index rows=matrix.rows(), cols=matrix.cols();
		out.write((char*) (&rows), sizeof(typename Matrix::Index));
		out.write((char*) (&cols), sizeof(typename Matrix::Index));
		out.write((char*) matrix.data(), rows*cols*sizeof(typename Matrix::Scalar) );
		out.close();
	}
	
	/*! Reads any Eigen matrix from a binary file written with write_binary. 
		@tparam Matrix - the type of the Eigen matrix to be read from file.
		@param filename - the name of the binary file containing the matrix.
		@param matrix - a reference to the matrix into which the matrix on file will be read.
	 */
	template<class Matrix>
	void read_binary(const char* filename, Matrix& matrix){
		std::ifstream in(filename, std::ios::in | std::ios::binary);
		typename Matrix::Index rows=0, cols=0;
		in.read((char*) (&rows),sizeof(typename Matrix::Index));
		in.read((char*) (&cols),sizeof(typename Matrix::Index));
		matrix.resize(rows, cols);
		in.read( (char *) matrix.data() , rows*cols*sizeof(typename Matrix::Scalar) );
		in.close();
	}
}
