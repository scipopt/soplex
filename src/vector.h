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
#pragma ident "@(#) $Id: vector.h,v 1.1 2001/11/06 16:18:33 bzfkocht Exp $"

#ifndef _VECTOR_H_
#define _VECTOR_H_

#include <assert.h>
#include <string.h>
#include <math.h>
#include <iostream>

#include "vector_c.h"

namespace soplex
{

class SLUFactor;
class SVector;
class SSVector;
class SubSVector;
class IdxSet;

/** Dense Vector.
 
    Class #Vector# provides dense linear algebra vectors. It does not
    provide memory management for the array of values. Instead, the
    constructor requires a pointer to a memory block large enough to
    fit the desired dimension of #double# values.
 
    After construction, the values of a #Vector# can be accessed with
    the subscript operator.  A #Vector# ensures that no access to
    values bejond the specified dimension is done. However, this
    feature can be turned off, if compiled with #-DNDEBUG#.
 
    A #Vector# is distinguished from a simple array of #double#s, by
    providing a set of mathematical operations. Since #Vector# does
    not provide any memory management features, no operations are
    available that would require allocation of temporary memory
    space. Here is a list of provided operations, where the examples
    asume that #Vector a, b# and #double x# have been declared.
 
    \begin{center}
    \begin{tabular}{lll}
        Operation   & Description      &                                  \\
        \hline
        #-=#        & subtraction      & #a -= b#                         \\
        #+=#        & addition         & #a += b#                         \\
        #*#         & scalar product   & #x = a * b#                      \\
        #*=#        & scaling          & #a *= x#                         \\
        #maxAbs()#  & infinity norm    & #a.maxAbs()# == $\|a\|_{\infty}$ \\
        #length()#  & norm             & #a.length()# == $\sqrt(a*a)$ \\
        #length2()# & squared norm     & #a.length2()# == $a*a$       \\
        #multAdd()# & add scaled vector& #a.multAdd(b, x)#                \\
    \end{tabular}
    \end{center}
 
    When using any of these operations, the vector involved must be of
    the same dimension. For #b# also #SVector b# are allowed, if it
    does not contain nonzeros with index greater than the dimension of
    #x#.
*/

class Vector
{
   /* A vector consists of its dimension #dimen# and a pointer to
    * its values #val#. Both are protected, since derived classes may want
    * to do some automatic memory management.
    */
   friend class LP;
   friend Vector& Usolve(Vector&, const SLUFactor&);
   friend Vector& Usolve2(Vector&, Vector&, const SLUFactor&);

protected:
   /// dimension of vector.
   int dimen;

   /** Pointer to values of a vector.
    *  The memory block pointed to by #val# must at least have size
    *  #dimen * sizeof(double)#.  
    */
   double* val;

public:
   /**@name Construction/Destruction
    *  There is no default constructor since the storage for a 
    *  #Vector# must be provides in any case.
    */
   //@{
   /** constructor.
    *  #Vector#s do not provide their own memory for storing values.
    *  Instead, it must be passed a memory block #val# at construction. It
    *  must be large enough to fit at least #dimen# double values.
    */
   //lint -e1712
   Vector(int p_dimen, double *p_val)
      : dimen(p_dimen)
         , val(p_val)
   {
      assert(dimen >= 0);
   }
   /// destructor
   //lint -esym(1540,Vector::val) we do not manage the storage val points to.
   //~Vector() { }

   //@}
   /**@name Query */
   //@{
   /// dimension of vector.
   int dim() const
   {
      return dimen;
   }
   //@}

   /// return #n#-th value.
   double& operator[](int n)
   {
      assert(n >= 0 && n < dim());
      return val[n];
   }

   /// return #n#-th value.
   double operator[](int n) const
   {
      assert(n >= 0 && n < dim());
      return val[n];
   }
   ///
   Vector& operator+=(const Vector& vec);
   ///
   Vector& operator+=(const SVector& vec);
   ///
   Vector& operator+=(const SubSVector& vec);
   ///
   Vector& operator+=(const SSVector& vec);

   ///
   Vector& operator-=(const Vector& vec);
   ///
   Vector& operator-=(const SVector& vec);
   ///
   Vector& operator-=(const SubSVector& vec);
   ///
   Vector& operator-=(const SSVector& vec);

   /// scaling.
   Vector& operator*=(double x);

   ///
   double operator*(const SSVector& v) const;
   ///
   double operator*(const SVector& v) const;
   ///
   double operator*(const SubSVector& v) const;
   /// inner product.
   double operator*(const Vector& v) const
   {
      assert(v.dim() == dim());
      return MultiplyVectorVector(val, dimen, v.val);
   }

   /// infinity norm of a Vector.
   double maxAbs() const;
   /// euclidian norm of a Vector.
   double length() const;
   /// squared norm of a Vector.
   double length2() const;

   ///
   Vector& multAdd(double x, const SVector& vec);
   ///
   Vector& multAdd(double x, const SubSVector& vec);
   ///
   Vector& multAdd(double x, const SSVector& svec);
   /// add scaled vector (#+= x*vec#).
   Vector& multAdd(double x, const Vector& vec)
   {
      assert(vec.dim() == dim());
      Vector_MultAddVector(x, dim(), this->get_ptr(), vec.get_const_ptr());
      return *this;
   }
   /** Conversion to C-style pointer.
       This function serves for using a #Vector# in an old C
       function. It returns a pointer to the first value of the array.
   */
   double* get_ptr()
   {
      return val;
   }
   const double* get_const_ptr() const
   {
      return val;
   }

   /// output operator.
   friend std::ostream& operator<<(std::ostream& s, const Vector& vec);

   /// consistency check.
   int isConsistent() const;

   /// set vector to 0.
   void clear()
   {
      if (dimen)
         memset(val, 0, dimen*sizeof(double));
   }

   /// zero values given by #idx#.
   Vector& clear(const IdxSet& idx);

   ///
   Vector& operator=(const Vector& vec);
   ///
   Vector& operator=(const SVector& vec);
   /** Assignment operator.
       Assigning a #SVector# or #SSVector# to a #Vector# using the operator
       will set all values to 0 except the nonzeros of the right handside
       #SVector#.  This is in constrast to method #assign()#.
   */
   Vector& operator=(const SSVector& vec);

   /// assign values of #sv# only.
   Vector& assign(const SVector& psv);
   /** Assign values of #sv#.
       Assigns all nonzeros of #sv# to a vector. All other values remain
       unchanged.
   */
   Vector& assign(const SSVector& sv);
};

} // namespace soplex
#endif // _VECTOR_H_

//-----------------------------------------------------------------------------
//Emacs Local Variables:
//Emacs c-basic-offset:3
//Emacs tab-width:8
//Emacs indent-tabs-mode:nil
//Emacs End:
//-----------------------------------------------------------------------------
