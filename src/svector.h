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
#pragma ident "@(#) $Id: svector.h,v 1.6 2001/11/13 21:01:27 bzfkocht Exp $"

#ifndef _SVECTOR_H_ 
#define _SVECTOR_H_

/*      \Section{Imports}
 */
#include <iostream>
#include <stdlib.h>
#include <assert.h>
#include <math.h>

#include "vector.h"

namespace soplex
{

//@ -----------------------------------------------------------------------------
/** Sparse vectors.
 *  Class #SVector# provides packed sparse vectors. Such are a sparse vectors,
 *  with a storage scheme that keeps all data in one contiguous block of memory.
 *  This is best suited for using them for parallel computing on a distributed
 *  memory multiprocessor.
 *
 * #SVector# does not provide any memory management (this will be done by class
 * #DSVector#). This means, that the constructor of #SVector# expects memory
 * where to save the nonzeros. Further, adding nonzeros to an #SVector# may fail
 * if no more memory is available for saving them (see also #DSVector#).
 *
 * When nonzeros are added to an #SVector#, they are appended to the set of
 * nonzeros, i.e. they recieve numbers #size()#, #size()+1# ... . An #SVector#
 * can hold atmost #max()# nonzeros, where #max()# is given in the constructor.
 * When removing nonzeros, the remaining nonzeros are renumbered. However, only
 * the numbers greater than the number of the first removed nonzero are
 * affected.
 *
 * The following mathematical operations are provided by class #SVector#
 * (#SVector a, b, c; double x#): \\
 * \begin{center}
 * \begin{tabular}{lll}
 *     Operation        & Description           & \\
 *       \hline
 *     #-=#             & subtraction           & #a -= b#      \\
 *     #+=#             & addition              & #a += b#      \\
 *     #*#              & skalar product        & #x = a * b#   \\
 *     #*=#             & scaling               & #a *= x#      \\
 *     #maxAbs()#       & infinity norm         & #a.maxAbs()# == $\|a\|_{\infty}$ \\
 *     #length()#       & eucledian norm        & #a.length()# == $\sqrt{a^2}$  \\
 *     #length2()#      & square norm           & #a.length2()# == $a^2$        \\
 * \end{tabular}
 * \end{center}
 * 
 * Operators #+=# and #-=# should be used with caution, since no efficient
 * implementation is available. One should think of assigning the left handside
 * vector to a dense #Vector# first and perform the addition on it. The same
 * applies to the scalar product #*#.
 *
 * There are two numberings of the nonzeros of an #SVector#. First, an #SVector#
 * is supposed to act like a linear algebra #Vector#. An {\em index} reffers to
 * this view of an #SVector#: #operator[]# is provided which return the value of
 * the vector to the given index, i.e. 0 for all indeces not in the set of
 * nonzeros.  The other view of #SVector#s is that of a set of nonzeros. The
 * nonzeros are numbered from 0 to #size()-1#. Methods #index(n)# and #value(n)#
 * allow to access the index and value of the #n#-th nonzero. #n# is reffered to
 * as the {\em number} of a nonzero.
 */
class SVector
{
   friend class Vector;
   friend class SSVector;
public:
   //typedef SVector_Element Element;

   /** Sparse vector nonzero element.
    *  #SVector# keep their nonzeros in an array of #Element#s providing
    *  members for saving the nonzero's index and value.
    */
   struct Element
   {
      /// Value of nonzero element
      double val;
      /// Index of nonzero element
      int idx;
   };


protected:
   /*   \Section{Datastructures}
        An #SVector# keeps its data in an array of #Element#s. The size and maximum
        number of elements allowed is stored in the -1st #Element# in its members
        #idx# and #val# respectively.
   */
   Element *m_elem;

public:
   /**@name Modification */
   //@{
   /// switch n'th with 0'th nonzero.
   void toFront(int n);

   /// append one nonzero #(i,v)#.
   void add(int i, double v)
   {
      int n = size();
      m_elem[n].idx = i;
      m_elem[n].val = v;
      set_size( n + 1 );
      assert(size() <= max());
   }

   /// append nonzeros of #sv#.
   void add(const SVector& sv)
   {
      add(sv.size(), sv.m_elem);
   }

   /// append #n# nonzeros.
   void add(int n, const int i[], const double v[]);

   /// append #n# nonzeros.
   void add(int n, const Element e[]);

   /// remove nonzeros n thru m.
   void remove(int n, int m);

   /// remove n-th nonzero.
   void remove(int n)
   {
      assert(n < size() && n >= 0);
      set_size( size() - 1 );
      m_elem[n] = m_elem[size()];
   }
   /// remove all indices.
   void clear()
   {
      set_size(0);
   }

   /// sort nonzero to increasing indices.
   void sort();
   //@}


   /**@name Inquiery*/
   //@{
   /// number of used indeces.
   int size() const
   {
      return m_elem[ -1].idx;
   }

   /// maximal number indeces.
   int max() const
   {
      return int(m_elem[ -1].val);
   }

   /// maximal index.
   int dim() const;

   /** Number of index #i#.
       Return the number of the first index #i#. If no index #i# is available
       in the #IdxSet#, -1 is returned. Otherwise, #index(number(i)) == i#
       hods.
   */
   int number(int i) const
   {
      int n = size();
      Element* e = &(m_elem[n]);
      while (n--)
      {
         --e;
         if (e->idx == i)
            return n;
      }
      return -1;
   }

   /// get value to index #i#.
   double operator[](int i) const
   {
      int n = number(i);
      if (n >= 0)
         return m_elem[n].val;
      return 0;
   }

   ///
   Element& element(int n)
   {
      assert(n >= 0 && n < max());
      return m_elem[n];
   }

   /// get #n#-th nonzero element.
   Element element(int n) const
   {
      assert(n >= 0 && n < size());
      return m_elem[n];
   }

   ///
   int& index(int n)
   {
      assert(n >= 0 && n < size());
      return m_elem[n].idx;
   }

   /// get index of #n#-th nonzero.
   int index(int n) const
   {
      assert(n >= 0 && n < size());
      return m_elem[n].idx;
   }

   ///
   double& value(int n)
   {
      assert(n >= 0 && n < size());
      return m_elem[n].val;
   }

   /// get value of #n#-th nonzero.
   double value(int n) const
   {
      assert(n >= 0 && n < size());
      return m_elem[n].val;
   }
   //@}


   /**@name Mathematical Operations */
   //@{
   /// infinity norm.
   double maxAbs() const;

   /// eucledian norm.
   double length() const
   {
      return sqrt(length2());
   }

   /// squared eucledian norm.
   double length2() const;

   /// scale with #x#.
   SVector& operator*=(double x);

   /// inner product.
   double operator*(const Vector& w) const
   {
      double x = 0;
      int n = size();
      Element* e = m_elem;

      while (n--)
      {
         x += e->val * w[e->idx];
         e++;
      }
      return x;
   }

   //@}


   /**@name Miscellaneous*/
   //@{
   ///
   friend std::ostream& operator<<(std::ostream& os, const SVector& v);

   ///
   SVector& operator=(const SSVector& sv);
   /// assignment operator.
   SVector& operator=(const SVector& sv);
   ///
   SVector& operator=(const Vector& sv);
   ///
   SVector& assign(const Vector& vec, double eps = 1e-12);

   /// consistency check.
   int isConsistent() const;

   /** default constructor.
       The constructor expects one memory block where to store the nonzero
       elements. This must passed to the constructor, where the {\em number
       of #Element#s} needs that fit into the memory must be given and a
       pointer to the begining of the memory block. Once this memory has
       been passed, it shall not be modified until the #SVector# is no
       longer used. Note, that when a memory block for $n$, say, #Element#s
       has been passed, only $n-1$ are available for actually storing
       nonzeros. The remaining one is used for bookkeeping purposes.
   */
   SVector(int n = 0, Element* p_mem = 0)
   {
      setMem(n, p_mem);
   }

   // Internals.
   Element* mem() const
   {
      return m_elem -1;
   }
//     int& size()
//     {
//        assert(m_elem != 0);
//        return m_elem[ -1].idx;
//     }
   void set_size(int s)
   {
      assert(m_elem != 0);
      m_elem[ -1].idx = s;
   }
      
   void set_max(int m)
   {
      assert(m_elem != 0);
      m_elem[ -1].val = m;
   }
   void setMem(int n, Element* elmem)
   {
      if (n)
      {
         assert(n > 0);
         assert(elmem != 0);
         elmem->val = 0;        // for purify to shut up
         m_elem = &(elmem[1]);
         set_size( 0 );
         set_max ( n - 1 );
      }
      else
         m_elem = 0;
   }
   //@}

};

inline Vector& Vector::multAdd(double x, const SVector& vec)
{
   assert(vec.dim() <= dim());
   Vector_MultAddSVector(val, x, vec.size(), vec.m_elem);
   return *this;
}

} // namespace soplex
#endif  // _SVECTOR_H_

//-----------------------------------------------------------------------------
//Emacs Local Variables:
//Emacs mode:c++
//Emacs c-basic-offset:3
//Emacs tab-width:8
//Emacs indent-tabs-mode:nil
//Emacs End:
//-----------------------------------------------------------------------------
