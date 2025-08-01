/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the class library                   */
/*       SoPlex --- the Sequential object-oriented simPlex.                  */
/*                                                                           */
/*  Copyright (c) 1996-2025 Zuse Institute Berlin (ZIB)                      */
/*                                                                           */
/*  Licensed under the Apache License, Version 2.0 (the "License");          */
/*  you may not use this file except in compliance with the License.         */
/*  You may obtain a copy of the License at                                  */
/*                                                                           */
/*      http://www.apache.org/licenses/LICENSE-2.0                           */
/*                                                                           */
/*  Unless required by applicable law or agreed to in writing, software      */
/*  distributed under the License is distributed on an "AS IS" BASIS,        */
/*  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied. */
/*  See the License for the specific language governing permissions and      */
/*  limitations under the License.                                           */
/*                                                                           */
/*  You should have received a copy of the Apache-2.0 license                */
/*  along with SoPlex; see the file LICENSE. If not email to soplex@zib.de.  */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/**@file  slinsolver_rational.h
 * @brief Sparse Linear Solver virtual base class with Rational precision.
 */
#ifndef _SLINSOLVER_RATIONAL_H_
#define _SLINSOLVER_RATIONAL_H_


#include <assert.h>
#include <string.h>

#include "soplex/spxdefines.h"
#include "soplex/svector.h"
#include "soplex/ssvector.h"
#include "soplex/dsvector.h"
#include "soplex/vector.h"
#include "soplex/didxset.h"

namespace soplex
{
/**@brief   Sparse Linear Solver virtual base class with Rational precision.
   @ingroup Algo

   Class SLinSolverRational provides a class for solving sparse linear systems with
   a matrix \f$A\f$ and arbitrary right-hand side vectors. For doing so, the
   matrix must be first #load%ed to an #SLinSolverRational object as an array of
   pointers to the \em column \ref SVectorRational "SVectorsRational" of this matrix.
*/
class SLinSolverRational
{
public:

   //---------------------------------------
   /**@name Types */
   ///@{
   /// status flags of the SLinSolverRational class.
   enum Status
   {
      /** The SLinSolverRational is ready for solving linear systems with the
          loaded matrix */
      OK       = 0,
      /** The loaded matrix allows only for instable solutions to be
          computed */
      INSTABLE = 1,
      /// The loaded matrix is singular.
      SINGULAR = 2,
      /// No matrix has yet been loaded.
      UNLOADED = 4,
      /// An error has occurred.
      ERROR    = 8,
      /// The time limit has been hit
      TIME     = 16
   };
   ///@}

   //---------------------------------------
   /**@name Miscellaneous */
   ///@{
   /// returns the name of the SLinSolverRational.
   virtual const char* getName() const = 0;

   /// returns the Status of the SLinSolverRational.
   virtual Status status() const = 0;

   /// unloads any matrix.
   virtual void clear() = 0;

   /// returns current memory consumption.
   virtual int memory() const = 0;

   /// returns dimension of loaded matrix.
   virtual int dim() const = 0;

   /// loads \p dim column vectors \p vec into the solver.
   /** Initializes SLinSolverRational for the solution of linear systems
       with the matrix consisting of \p dim column vectors given in \p vec.
   */
   virtual Status load(const SVectorRational* vec[], int dim) = 0;

   /// returns a stability number (0: singularity, 1: perfect stability).
   /** Returns a stability parameter between 0 and 1, where 0 indicates
       singularity, while 1 indicates perfect stability.
   */
   virtual Rational stability() const = 0;

   /// returns statistical information in form of a string.
   virtual std::string statistics() const = 0;

   /// Substitute column \p idx with \p subst.
   /** The change method is used to modify the loaded matrix by substituting
       column \p idx with the new vector \p subst. One may also pass the
       optional parameter \p eta to the solution of #solveRight() if
       readily  availabble. This may improve on the performance of the update.
   */
   virtual Status change(int idx, const SVectorRational& subst,
                         const SSVectorRational* eta = nullptr) = 0;

   /// consistency check.
   virtual bool isConsistent() const = 0;

   /// get number of factorizations
   virtual int getFactorCount() const = 0;
   ///@}


   /**@name Solving linear systems
      For solving linear systems with an SLinSolverRational object, it must
      have previously been loaded with the matrix to use.

      Two types of systems can be solved \f$A x = b\f$ and \f$x^T A = b^T\f$.
      Method names related to the first and second type are solveRight() and
      solveLeft(), respectively.

      The methods receive their right hand-side vector \f$b\f$ as a
      \c const parameter, that will hence be unchanged after termination.

      Some methods are available with two parameters for right hand-side
      vectors. Then two system are solved in one method invocation. This
      should generally be faster than solving two systems seperately.

      The result vector(s) are allways given as the first parameter(s). Two
      types of result vectors are supported, VectorRational and SSVectorRational.
   */
   ///@{
   /// Solves \f$Ax=b\f$.
   virtual void solveRight(VectorRational& x, const VectorRational& b) /* const */ = 0;
   /// Solves \f$Ax=b\f$.
   virtual void solveRight(SSVectorRational& x, const SVectorRational& b) /* const */ = 0;

   /** @brief Solves \f$Ax=b\f$.
       Possibly sets up internal data structures suitable for an optimized
       subsequent change() call with \f$b\f$ as entering column.
   */
   virtual void solveRight4update(SSVectorRational& x, const SVectorRational& b) = 0;

   /// Solves \f$Ax=b\f$ and \f$Ay=d\f$.
   virtual void solve2right4update(SSVectorRational& x,
                                   VectorRational& y,
                                   const SVectorRational& b,
                                   SSVectorRational& d) = 0;
   /// Solves \f$Ax=b\f$, \f$Ay=d\f$ and \f$Az=e\f$.
   virtual void solve3right4update(SSVectorRational& x,
                                   VectorRational& y,
                                   VectorRational& z,
                                   const SVectorRational& b,
                                   SSVectorRational& d,
                                   SSVectorRational& e) = 0;
   /// solves \f$x^TA=b^T\f$.
   virtual void solveLeft(VectorRational& x, const VectorRational& b) /* const */ = 0;
   /// solves \f$x^TA=b^T\f$.
   virtual void solveLeft(SSVectorRational& x, const SVectorRational& b) /* const */ = 0;
   /// solves \f$x^TA=b^T\f$ and \f$x^TA=rhs2^T\f$ internally using \f$rhs2\f$.
   virtual void solveLeft(SSVectorRational& x,
                          VectorRational& two,
                          const SVectorRational& b,
                          SSVectorRational& rhs2) /* const */ = 0;
   /// solves \f$x^TA=b^T\f$, \f$y^TA=d^T\f$ and \f$z^TA=e^T\f$
   virtual void solveLeft(SSVectorRational& x, VectorRational& y, VectorRational& z,
                          const SVectorRational& b, SSVectorRational& d, SSVectorRational& e) = 0;
   ///@}


   //---------------------------------------
   /**@name Constructors / Destructors */
   ///@{
   /// default constructor
   SLinSolverRational() = default;
   /// copy constructor
   SLinSolverRational(const SLinSolverRational&) = default;
   /// move constructor
   SLinSolverRational(SLinSolverRational&&) = default;
   /// destructor
   virtual ~SLinSolverRational() = default;
   /// clone function for polymorphism
   virtual SLinSolverRational* clone() const = 0;
   ///@}



};

} // namespace soplex
#endif // _SLINSOLVER_RATIONAL_H_
