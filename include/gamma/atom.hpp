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

#ifndef ATOMHEADERDEF
#define ATOMHEADERDEF

// Includes
#include "eigen_wrapper.hpp"
#include "basis.hpp"
#include <libecpint/ecp.hpp>

/*! \class Atom
	\brief The building block of a Molecule.
 
 	Contains the charge, mass, and position of an atom.
	Also contains the effective charge if ECPs are in use.
	Provides methods for coordinate transformations.
 */
class Atom
{
private:
  int charge; ///< Total atomic charge, not taking into account core.
  int core; ///< Number of frozen core electrons
  
  /// xyz coordinates of the atom in atomic units.
  double x, y, z;
  double mass; ///< Mass of the atom in atomic units.
  double *pos; ///< xyz coordinates of the atom as an array, in atomic units.

  public:
  /// Creates an empty Atom with negative charge positioned at the origin.
  Atom() : charge(-1), core(0) { } 
  
  /*! Creates a new Atom at the given coordinates with the given charge and mass. 
   	  @param coords - a reference to a Vector of the xyz coordinates of the Atom in a.u.
  	  @param q - the charge of the Atom.
  	  @param m - the mass of the Atom in a.u.
   */
  Atom(const Vector& coords, int q, double m);
  
  /// Copy constructor 
  Atom(const Atom& other); // Copy constructor

  /*! Get the total charge of the Atom.
  	@return The total charge of the Atom.
   */
  int getCharge() const { return charge; }
  
  /*! Get the charge of the Atom with core electrons removed.
   @return The charge of the Atom with core electrons removed.
   */
  int getEffectiveCharge() const { return charge - core; }
  
  /*! Sets the number of core electrons for this Atom.
   @param ecpset - the ECP basis that defines how many core electrons each atom type has.
   */
  void setCore(libecpint::ECPBasis& ecpset); 
  
  /// @return The mass in atomic units of the Atom.
  double getMass() const { return mass; }

  /// @return A Vector of the xyz coordinates of the Atom in atomic units.
  Vector getCoords() const; 
  
  /// @return An array of the xyz coordinates of the Atom in atomic units.
  double* getPos() const { return pos; }
  
  double getX() const { return x; } ///< @return The x-coordinate of the Atom in a.u.
  double getY() const { return y; } ///< @return The y-coordinate of the Atom in a.u.
  double getZ() const { return z; } ///< @return the z-coordinate of the Atom in a.u.

  /*! Rotates the Atom according to a unitary transformation matrix.
      @param U - the unitary transformation matrix for the rotation.
   */
  void rotate(const Matrix& U); 
  
  /*! Translates the Atom from its current position.
      @param dx - the change in x-coordinate in a.u.
  	  @param dy - the change in y-coordinate in a.u.
  	  @param dz - the change in z-coordinate in a.u.
   */
  void translate(double dx, double dy, double dz);
  
  /// Overloaded equality operator.
  Atom& operator=(const Atom& other);
};  
  
#endif
