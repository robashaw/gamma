/*
 *
 *    PURPOSE: Implements atom.hpp, defining class Atom
 * 
 *    DATE         AUTHOR            CHANGES
 *    =======================================================================
 *    27/08/15     Robert Shaw       Original code.
 *    05/09/15     Robert Shaw       Added function to count number of cgbfs
 *                                   in spherical basis.
 */

// Includes
#include "atom.hpp"
#include <iostream>

//Constructors
Atom::Atom(const Vector& coords, int q, double m)
{
  x = coords(0);
  y = coords(1);
  z = coords(2);
  core = 0;
  pos = new double[3];
  pos[0] = x;
  pos[1] = y;
  pos[2] = z;
  charge = q;
  mass = m;
}

// Copy constructor
Atom::Atom(const Atom& other)
{
  x = other.x; 
  y = other.y; 
  z = other.z; 
  core = other.core;
  pos = new double[3];
  pos[0] = x;
  pos[1] = y;
  pos[2] = z;
  mass = other.mass;
  charge = other.charge;
}

// Coordinate accessor
Vector Atom::getCoords() const
{
  Vector v(3);
  v[0] = x; v[1] = y; v[2] = z;
  return v;
}

void Atom::setCore(libecpint::ECPBasis& ecpset) {
	core = ecpset.getECPCore(charge);	
}	

// Routines

// Translate(dx, dy, dz) translates the coordinates by [dx, dy, dz]
// while rotate(Matrix U) applies the unitary transformation U to the coords
void Atom::rotate(const Matrix& U)
{
  // Brute force is quicker than matrix-vector multiplication for
  // a 3x3 situation such as this
  x = U(0, 0)*x + U(0, 1)*y + U(0, 2)*z;
  y = U(1, 0)*x + U(1, 1)*y + U(1, 2)*z;
  z = U(2, 0)*x + U(2, 1)*y + U(2, 2)*z;
  pos[0] = x; pos[1] = y; pos[2] = z;
}

void Atom::translate(double dx, double dy, double dz)
{
  x += dx; y += dy; z += dz;
  pos[0] = x; pos[1] = y; pos[2] = z;
}

// Overloaded operators

// Overload the assignment operator, =
Atom& Atom::operator=(const Atom& other)
{
  // Assign attributes
  charge = other.charge;
  core = other.core;
  x = other.x; y = other.y; z = other.z;
  pos = new double[3];
  pos[0] = other.pos[0];
  pos[1] = other.pos[1];
  pos[2] = other.pos[2];
  mass = other.mass;
  
  return *this;
}

