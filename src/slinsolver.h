/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the class library                   */
/*       SoPlex --- the Sequential object-oriented simPlex.                  */
/*                                                                           */
/*    Copyright (C) 1997-1999 Roland Wunderling                              */
/*                  1997-2001 Konrad-Zuse-Zentrum                            */
/*                            fuer Informationstechnik Berlin                */
/*                                                                           */
/*  SoPlex is distributed under the terms of the ZIB Academic Licence.       */
/*                                                                           */
/*  You should have received a copy of the ZIB Academic License              */
/*  along with SoPlex; see the file COPYING. If not email to soplex@zib.de.  */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
#pragma ident "@(#) $Id: slinsolver.h,v 1.3 2001/11/07 17:31:20 bzfbleya Exp $"

#ifndef _SLINSOLVER_H_
#define _SLINSOLVER_H_


//@ ----------------------------------------------------------------------------
/*      \Section{Imports}
    Import required system include files
 */
#include <assert.h>


/*  and class header files
 */

#include "svector.h"
#include "ssvector.h"
#include "dsvector.h"
#include "dvector.h"
#include "didxset.h"

namespace soplex
{

//@ ----------------------------------------------------------------------------
/*      \Section{Class Declaration}
 */
/** Sparse Linear Solver virtual base class.
    Class #SLinSolver# provides a class for solving sparse linear systems with
    a matrix $A$ and arbitrary right-hand side vectors. For doing so, the matrix
    must be first #load#ed to an #SLinSolver# object as an array of pointers to
    the {\em column} #SVector#s of this matrix.
 */
class SLinSolver
{
public:
   ///
   enum Status
   {
      /*/ The #SLinSolver# is ready for solving linear systems with the
          loaded matrix
       */
      OK = 0,
      /*/ The loaded matrix allows only for instable solutions to be
          computed
       */
      INSTABLE = 1 ,
      /// The loaded matrix is singular.
      SINGULAR = 2 ,
      /// No matrix has yet been loaded.
      UNLOADED = 4 ,
      /// An error has occurred.
      ERROR = 8
   };

   ///
   virtual Status status() const = 0;

   /// unload any matrix.
   virtual void clear() = 0;

   /// return current memory consumption.
   virtual int memory() const = 0;

   /// return dimension of loaded matrix.
   virtual int dim() const = 0;

   /*/ Initialize #SLinSolver# for the solution of linear systems
       with the matrix consisting of #dim# column vectors given in #vec#
    */
   virtual Status load(const SVector* vec[], int dim) = 0;

   /*/     return a stability parameter between 0 and 1. 0 indicates
           singularity, while 1 indicates perfect stability.
    */
   virtual double stability() const = 0;


   /**@name Solving linear systems
       For solving linear systems with an #SLinSolver# object, it must
       have previously been #load#ed with the matrix to use.

       Two types of systems can be solved $A x = b$ and $x^T A = b^T$. Method
       names related to the first and second type are #solveRight# and
       #solveLeft#, respectively.

       The default methods receive their right hand-side vector $b$ as a
       #const# parameter, that will hence be unchanged after termination.
       Instead methods named #solve2...# are allowed to use the passed
       non-#const# right hand-side vector internally for the solution. This may
       save one copy operation, or may be used for solving different linear
       systems concurrently. Upon termination, this vector is cleared.

       Some methods are available with two parameters for right hand-side
       vectors. Then two system are solved in one method invocation. This
       should generally be faster than solving two systems seperately.

       The result vector(s) are allways given as the first parameter(s). Two
       types of result vectors are supported, #Vector# and #SSVector#.
    */
   //@{
   /// solve $Ax=b$ using $b$ and leaving $b=0$.
   virtual void solve2right(Vector& x, Vector& b) /* const */ = 0;
   /// solve $Ax=b$ using $b$ and leaving $b=0$.
   virtual void solve2right(Vector& x, SSVector& b) /* const */ = 0;
   /// solve $Ax=b$ using $b$ and leaving $b=0$.
   virtual void solve2right(SSVector& x, Vector& b) /* const */ = 0;
   /// solve $Ax=b$ using $b$ and leaving $b=0$.
   virtual void solve2right(SSVector& x, SSVector& b) /* const */ = 0;

   /// solve $Ax=b$.
   virtual void solveRight (Vector& x,
                             const Vector& b) /* const */ = 0;
   /// solve $Ax=b$.
   virtual void solveRight (Vector& x,
                             const SVector& b) /* const */ = 0;
   /// solve $Ax=b$.
   virtual void solveRight (SSVector& x,
                             const Vector& b) /* const */ = 0;
   /// solve $Ax=b$.
   virtual void solveRight (SSVector& x,
                             const SVector& b) /* const */ = 0;

   /** Solve $Ax=b$.
    * This method solves $Ax=b$ thereby (possibly) setting up internal data
    * structures suitable for an optimized subsequent #change# call with
    * $b$ as entering column.
    */
   virtual void solveRight4update(SSVector& x,
                                   const SVector& b) = 0;

   ///
   virtual void solve2right4update(SSVector& x,
                                    Vector& two,
                                    const SVector& b,
                                    SSVector& rhs) = 0;

   /// solve $x^TA=b$ using $b$ and leaveing $b=0$.
   virtual void solve2left(Vector& x, Vector& b) /* const */ = 0;
   /// solve $x^TA=b$ using $b$ and leaveing $b=0$.
   virtual void solve2left(Vector& x, SSVector& b) /* const */ = 0;
   /// solve $x^TA=b$ using $b$ and leaveing $b=0$.
   virtual void solve2left(SSVector& x, Vector& b) /* const */ = 0;
   /// solve $x^TA=b$ using $b$ and leaveing $b=0$.
   virtual void solve2left(SSVector& x, SSVector& b) /* const */ = 0;

   /// solve $x^TA=b$.
   virtual void solveLeft (Vector& x,
                            const SVector& b) /* const */ = 0;
   /// solve $x^TA=b$.
   virtual void solveLeft (Vector& x,
                            const Vector& b) /* const */ = 0;
   /// solve $x^TA=b$.
   virtual void solveLeft (SSVector& x,
                            const SVector& b) /* const */ = 0;
   /// solve $x^TA=b$.
   virtual void solveLeft (SSVector& x,
                            const Vector& b) /* const */ = 0;
   /// solve $x^TA=b$ and $x^TA=rhs2$ internally using $rhs2$.
   virtual void solveLeft (SSVector& x,
                            Vector& two,
                            const SVector& b,
                            SSVector& rhs2) /* const */ = 0;
   //@}


   /** Substitute column #idx# with #subst#.
    *  The #change# method is used to modify the loaded matrix by substituting
    *  column #idx# with the new vector #subst#. One may also pass the optional
    *  parameter #eta# to the solution of #solveRight(subst)# if readily
    *  availabble. This may improve on the performance of the update.
    */
   virtual Status change(int idx, const SVector& subst,
                         const SSVector* eta = 0) = 0;

   ///
   virtual ~SLinSolver()
   {}
   ///

   virtual int isConsistent() const = 0;
};

} // namespace soplex
#endif

//-----------------------------------------------------------------------------
//Emacs Local Variables:
//Emacs mode:c++
//Emacs c-basic-offset:3
//Emacs tab-width:8
//Emacs indent-tabs-mode:nil
//Emacs End:
//-----------------------------------------------------------------------------
