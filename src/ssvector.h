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
#pragma ident "@(#) $Id: ssvector.h,v 1.6 2001/11/25 14:58:29 bzfkocht Exp $"


#ifndef _SSVECTOR_H_
#define _SSVECTOR_H_

#include <assert.h>

#include "dvector.h"
#include "subsvector.h"
#include "svector.h"
#include "didxset.h"

namespace soplex
{
class SVSet;

/** semi sparse vector.
    This class implements {\bf S}emi {\bf S}parse {\bf Vector}s. Such are
    #DVector#s where the indices of its nonzero elements can be stored in an
    extra #IdxSet#. Only elements with absolute value #> epsilon# are considered
    to be nonzero. Since really storing the nonzeros is not allways convenient,
    an #SSVector# provides two different statuses: setup and not setup.
    An #SSVector# being setup means that the nonzero indices are available,
    otherwise an #SSVector# is just an ordinary #Vector# with an empty #IdxSet#.
 */
class SSVector : protected DVector, protected DIdxSet
{
   SSVector& assign2product1(const SVSet& A, const SSVector& x);
   SSVector& assign2productShort(const SVSet& A, const SSVector& x);
   SSVector& assign2productFull(const SVSet& A, const SSVector& x);

   int setupStatus;

   friend class DVector;
   friend class Vector;
   friend class DSVector;
   friend class SMoPlex;

public:
   ///
   double epsilon;

   /**@name Status of an #SSVector#
       An #SSVector# can be setup or not. In case it is setup, its #IdxSet#
       correctly contains all indices of nonzero elements of the #SSVector#.
       Otherwise, it does not contain any usefull data. Wheter or not an
       #SSVector# is setup can be determined with method #isSetup()#.

       There are three method for directly affecting the setup status of an
       #SSVector#
       \begin{description}
       \item[unSetup]          This method sets the status to ``not setup''.
       \item[setup]            This method initializes the #IdxSet# to the
                               #SSVector#s nonzero indices and sets the status
                               to ``setup''.
       \item[forceSetup]       This method sets the status to ``setup'' without
                               verifying, that the #IdxSet# correctly contains
                               all nonzero indices. It may be used when the
                               nonzero indices have been computed externally.
       \end{description}
    */
   //@{

   /// only used in slufactor.cpp
   double* get_ptr()
   {
      return DVector::get_ptr();
   }
   /// return setup status.
   int isSetup() const
   {
      return setupStatus;
   }

   /// make #SSVector# not setup.
   void unSetup()
   {
      setupStatus = 0;
   }

   /*/ initialize nonzero indices for all elements with absolute values
       #> eps# and set all other elements to 0.
    */
   void setup();

   /// force setup status.
   void forceSetup()
   {
      setupStatus = 1;
   }
   //@}

   /**@name Methods for a setup #SSVectors# */
   //@{
   /// return #n#-th index.
   int index(int n) const
   {
      assert(isSetup());
      return IdxSet::index(n);
   }

   /// return value to #n#-th index.
   double value(int n) const
   {
      assert(isSetup());
      assert(n >= 0 && n < size());
      return val[idx[n]];
   }

   /// find number of index #i# (or -1).
   int number(int i) const
   {
      assert(isSetup());
      return IdxSet::number(i);
   }

   /// number of nonzeros.
   int size() const
   {
      assert(isSetup());
      return IdxSet::size();
   }

   /*/ add nonzero #(i,x)# to #SSVector#
       (no nonzero with index #i# must exist!)
    */
   void add(int i, double x)
   {
      assert(val[i] == 0);
      assert(number(i) < 0);
      addIdx(i);
      val[i] = x;
   }

   /// set #i#-t element to x.
   void setValue(int i, double x);

   /// clear element #i#.
   void clearIdx(int i)
   {
      if (isSetup())
      {
         int n = number(i);
         if (n >= 0)
            remove(n);
      }
      val[i] = 0;
   }

   /// set #n#-th nonzero element to 0 (must exist!).
   void clearNum(int n)
   {
      assert(isSetup());
      assert(index(n) >= 0);
      val[index(n)] = 0;
      remove(n);
   }
   //@}


   /**@name Methods independend on the Status */
   //@{
   /// return #i#-th value.
   double operator[](int i) const
   {
      return val[i];
   }

   /// return array indices.
   const int* indexMem() const
   {
      return IdxSet::indexMem();
   }

   /// return array values.
   const double* values() const
   {
      return val;
   }

   /// return indices.
   const IdxSet& indices() const
   {
      return *this;
   }

   /// return array indices.
   int* altIndexMem()
   {
      unSetup();
      return IdxSet::indexMem();
   }

   /// return array values.
   double* altValues()
   {
      unSetup();
      return val;
   }

   /// return indices.
   IdxSet& altIndices()
   {
      unSetup();
      return *this;
   }
   //@}


   /**@name Mathematical opeations */
   //@{
   ///
   SSVector& operator+=(const Vector& vec);
   ///
   SSVector& operator+=(const SVector& vec);
   ///
   SSVector& operator+=(const SubSVector& vec);
   ///
   SSVector& operator+=(const SSVector& vec);

   ///
   SSVector& operator-=(const Vector& vec);
   ///
   SSVector& operator-=(const SVector& vec);
   ///
   SSVector& operator-=(const SubSVector& vec);
   ///
   SSVector& operator-=(const SSVector& vec);

   ///
   SSVector& operator*=(double x);

   /// add scaled vector (#+= x*vec#).
   SSVector& multAdd(double x, const SSVector& vec);
   ///
   SSVector& multAdd(double x, const SVector& vec);
   ///
   SSVector& multAdd(double x, const SubSVector& vec);
   ///
   SSVector& multAdd(double x, const Vector& vec);

   /// assign #SSVector# to $x^T \cdot A$.
   SSVector& assign2product(const SSVector& x, const SVSet& A);
   /// assign #SSVector# to $A \cdot x$.
   SSVector& assign2product(const SVSet& A, const SSVector& x);
   /// assign #SSVector# to $A \cdot x$ for a setup #x#.
   SSVector& assign2product4setup(const SVSet& A, const SSVector& x);
   /// assign #SSVector# to $A \cdot x$ thereby setting up #x#.
   SSVector& assign2productAndSetup(const SVSet& A, SSVector& x);

   /// infinity norm of a Vector.
   double maxAbs() const;
   /// euclidian norm of a Vector.
   double length() const;
   /// squared norm of a Vector.
   double length2() const;
   //@}

   /**@name Miscellaneous */
   //@{
   ///
   int dim() const
   {
      return dimen;
   }
   /// reset dimension.
   void reDim (int newdim);
   /// set number of nonzeros (thereby unSetup SSVector).
   void setSize(int n)
   {
      unSetup();
      IdxSet::setSize(n);
   }
   /// reset memory consumption.
   void reMem(int newsize);
   /// set to 0.
   void clear ();

   ///
   SSVector& operator=(const SSVector& rhs);

 public:
   /// setup #rhs# vector, if it is not allready.
   void setup_and_assign(SSVector& rhs);
   ///
   SSVector& operator=(const SVector& rhs);
   ///
   SSVector& operator=(const Vector& rhs)
   {
      unSetup();
      Vector::operator=(rhs);
      return *this;
   }

   /// assign only the elements of #rhs#.
   SSVector& assign(const SVector& rhs);

   /// construct nonsetup copy of #vec#.
   SSVector(const Vector& vec, double eps = 1e-16)
      : DVector (vec)
      , DIdxSet (vec.dim() + 1)
      , setupStatus(0)
      , epsilon (eps)
   { }

   ///
   SSVector(int pdim = 0, double peps = 1e-16)
      : DVector (pdim)
      , DIdxSet (pdim + 1)
      , setupStatus(1)
      , epsilon (peps)
   {
      Vector::clear();
   }

   ///
   SSVector(const SSVector& vec)
      : DVector (vec)
      , DIdxSet (vec.dim() + 1)
      , setupStatus(vec.setupStatus)
      , epsilon (vec.epsilon)
   {
      DIdxSet::operator= ( vec );
      //*((DIdxSet*)this) = vec;
   }

   ///
   int isConsistent() const;
   //@}
};


//@ ----------------------------------------------------------------------------

inline Vector& Vector::multAdd(double x, const SSVector& svec)
{
   assert(svec.dim() <= dim());

   if (svec.isSetup())
   {
      const int* idx = svec.indexMem();

      for(int i = 0; i < svec.size(); i++)
         val[idx[i]] += x * svec[idx[i]];
   }
   else
   {
      assert(svec.dim() == dim());

      for(int i = 0; i < dim(); i++)
         val[i] += x * svec.val[i];
   }
   //multAdd(x, static_cast<const Vector&>(svec));

   return *this;
}

inline Vector& Vector::assign(const SSVector& svec)
{
   assert(svec.dim() <= dim());

   if (svec.isSetup())
   {
      const int* idx = svec.indexMem();

      for(int i = svec.size(); i > 0; i--)
      {
         val[*idx] = svec.val[*idx];
         idx++;
      }
   }
   else
      operator= (static_cast<const Vector&>(svec));

   return *this;
}

inline Vector& Vector::operator=(const SSVector& vec)
{
   if (vec.isSetup())
   {
      clear ();
      assign(vec);
   }
   else
      operator= (static_cast<const Vector&>(vec));

   return *this;
}

inline double Vector::operator*(const SSVector& v) const
{
   assert(dim() == v.dim());

   if (v.isSetup())
   {
      const int* idx = v.indexMem();
      double     x   = 0;

      for(int i = v.size(); i > 0; i--)
      {
         x += val[*idx] * v.val[*idx];
         idx++;
      }
      return x;
   }
   else
      return operator*(static_cast<const Vector&>(v));
}
} // namespace soplex
#endif // _SSVECTOR_H_

//-----------------------------------------------------------------------------
//Emacs Local Variables:
//Emacs mode:c++
//Emacs c-basic-offset:3
//Emacs tab-width:8
//Emacs indent-tabs-mode:nil
//Emacs End:
//-----------------------------------------------------------------------------
