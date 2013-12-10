/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the class library                   */
/*       SoPlex --- the Sequential object-oriented simPlex.                  */
/*                                                                           */
/*    Copyright (C) 1996-2013 Konrad-Zuse-Zentrum                            */
/*                            fuer Informationstechnik Berlin                */
/*                                                                           */
/*  SoPlex is distributed under the terms of the ZIB Academic Licence.       */
/*                                                                           */
/*  You should have received a copy of the ZIB Academic License              */
/*  along with SoPlex; see the file COPYING. If not email to soplex@zib.de.  */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/**@file  svectorbase.h
 * @brief Sparse vectors.
 */
#ifndef _SVECTORBASE_H_
#define _SVECTORBASE_H_

#include <iostream>
#include <assert.h>
#include <math.h>

namespace soplex
{
template < class R > class VectorBase;
template < class R > class SSVectorBase;

/// Sparse vector nonzero element.
/** SVectorBase keeps its nonzeros in an array of Nonzero%s providing members for saving the index and value.
 */
template < class R >
class Nonzero
{
public:

   R val;        ///< Value of nonzero element.
   int idx;      ///< Index of nonzero element.

   template < class S >
   Nonzero<R>& operator=(const Nonzero<S>& vec)
   {
      val = R(vec.val);
      idx = vec.idx;
      return *this;
   }

   template < class S >
   Nonzero<R>(const Nonzero<S>& vec)
      : val(vec.val)
      , idx(vec.idx)
   {
   }

   Nonzero<R>()
      : val()
      , idx(0)
   {
   }
};

/**@brief   Sparse vectors.
 * @ingroup Algebra
 *
 *  Class SVectorBase provides packed sparse vectors. Such are a sparse vectors, with a storage scheme that keeps all
 *  data in one contiguous block of memory.  This is best suited for using them for parallel computing on a distributed
 *  memory multiprocessor.
 *
 *  SVectorBase does not provide any memory management (this will be done by class DSVectorBase). This means, that the
 *  constructor of SVectorBase expects memory where to save the nonzeros. Further, adding nonzeros to an SVectorBase may
 *  fail if no more memory is available for saving them (see also DSVectorBase).
 *
 *  When nonzeros are added to an SVectorBase, they are appended to the set of nonzeros, i.e., they recieve numbers
 *  size(), size()+1 ... . An SVectorBase can hold atmost max() nonzeros, where max() is given in the constructor.  When
 *  removing nonzeros, the remaining nonzeros are renumbered. However, only the numbers greater than the number of the
 *  first removed nonzero are affected.
 *
 *  The following mathematical operations are provided by class SVectorBase (SVectorBase \p a, \p b, \p c; R \p x):
 *
 *  <TABLE>
 *  <TR><TD>Operation</TD><TD>Description   </TD><TD></TD>&nbsp;</TR>
 *  <TR><TD>\c -=    </TD><TD>subtraction   </TD><TD>\c a \c -= \c b </TD></TR>
 *  <TR><TD>\c +=    </TD><TD>addition      </TD><TD>\c a \c += \c b </TD></TR>
 *  <TR><TD>\c *     </TD><TD>skalar product</TD>
 *      <TD>\c x = \c a \c * \c b </TD></TR>
 *  <TR><TD>\c *=    </TD><TD>scaling       </TD><TD>\c a \c *= \c x </TD></TR>
 *  <TR><TD>maxAbs() </TD><TD>infinity norm </TD>
 *      <TD>\c a.maxAbs() == \f$\|a\|_{\infty}\f$ </TD></TR>
 *  <TR><TD>length() </TD><TD>eucledian norm</TD>
 *      <TD>\c a.length() == \f$\sqrt{a^2}\f$ </TD></TR>
 *  <TR><TD>length2()</TD><TD>square norm   </TD>
 *      <TD>\c a.length2() == \f$a^2\f$ </TD></TR>
 *  </TABLE>
 *
 *  Operators \c += and \c -= should be used with caution, since no efficient implementation is available. One should
 *  think of assigning the left handside vector to a dense VectorBase first and perform the addition on it. The same
 *  applies to the scalar product \c *.
 *
 *  There are two numberings of the nonzeros of an SVectorBase. First, an SVectorBase is supposed to act like a linear
 *  algebra VectorBase. An \em index refers to this view of an SVectorBase: operator[]() is provided which returns the
 *  value at the given index of the vector, i.e., 0 for all indices which are not in the set of nonzeros.  The other view
 *  of SVectorBase%s is that of a set of nonzeros. The nonzeros are numbered from 0 to size()-1.  The methods index(int
 *  n) and value(int n) allow to access the index and value of the \p n 'th nonzero.  \p n is referred to as the \em
 *  number of a nonzero.
 *
 *  @todo SVectorBase should get a new implementation.  There maybe a lot of memory lost due to padding the Nonzero
 *        structure. A better idea seems to be class SVectorBase { int size; int used; int* idx; R* val; }; which for
 *        several reasons could be faster or slower.  If SVectorBase is changed, also DSVectorBase and SVSet have to be
 *        modified.
 */
template < class R >
class SVectorBase
{
   template < class S > friend class SVectorBase;

#if 0 // needed?
   friend class VectorBase;
   friend class SSVector;
   friend std::ostream& operator<<(std::ostream& os, const SVectorBase<R>& v);
#endif

private:

   // ------------------------------------------------------------------------------------------------------------------
   /**@name Data */
   //@{

   Nonzero<R>* m_elem;
   int memsize;
   int memused;

   //@}

public:

   typedef Nonzero<R> Element;

   // ------------------------------------------------------------------------------------------------------------------
   /**@name Access */
   //@{

   /// Number of used indices.
   int size() const
   {
      assert(m_elem != 0 || memused == 0);
      return memused;
   }

   /// Maximal number of indices.
   int max() const
   {
      assert(m_elem != 0 || memused == 0);
      return memsize;
   }

   /// Dimension of the vector defined as maximal index + 1
   int dim() const
   {
      const Nonzero<R>* e = m_elem;
      int d = -1;
      int n = size();

      while( n-- )
      {
         d = (d > e->idx) ? d : e->idx;
         e++;
      }

      return d+1;
   }

   /// Number of index \p i.
   /** @return The number of the first index \p i. If no index \p i is available in the IdxSet, -1 is
    *          returned. Otherwise, index(number(i)) == i holds.
    */
   int number(int i) const
   {
      if( m_elem != 0 )
      {
         int n = size();
         Nonzero<R>* e = &(m_elem[n]);

         while( n-- )
         {
            --e;
            if( e->idx == i )
            {
               assert(index(n) == i);
               return n;
            }
         }
      }

      return -1;
   }

   /// Value to index \p i.
   R operator[](int i) const
   {
      int n = number(i);

      if( n >= 0 )
         return m_elem[n].val;

      return 0;
   }

   /// Reference to the \p n 'th nonzero element.
   Nonzero<R>& element(int n)
   {
      assert(n >= 0);
      assert(n < max());

      return m_elem[n];
   }

   /// The \p n 'th nonzero element.
   const Nonzero<R>& element(int n) const
   {
      assert(n >= 0);
      assert(n < size());

      return m_elem[n];
   }

   /// Reference to index of \p n 'th nonzero.
   int& index(int n)
   {
      assert(n >= 0);
      assert(n < size());

      return m_elem[n].idx;
   }

   /// Index of \p n 'th nonzero.
   int index(int n) const
   {
      assert(n >= 0);
      assert(n < size());

      return m_elem[n].idx;
   }

   /// Reference to value of \p n 'th nonzero.
   R& value(int n)
   {
      assert(n >= 0);
      assert(n < size());

      return m_elem[n].val;
   }

   /// Value of \p n 'th nonzero.
   R value(int n) const
   {
      assert(n >= 0);
      assert(n < size());

      return m_elem[n].val;
   }

   /// Append one nonzero \p (i,v).
   void add(int i, R v)
   {
      assert(m_elem != 0);
      assert(size() < max());

      int n = size();

      m_elem[n].idx = i;
      m_elem[n].val = v;
      set_size( n + 1 );

      assert(size() <= max());
   }

   /// Append nonzeros of \p sv.
   void add(const SVectorBase& sv)
   {
      add(sv.size(), sv.m_elem);
   }

   /// Append \p n nonzeros.
   void add(int n, const int i[], const R v[])
   {
      assert(n + size() <= max());

      if( n <= 0 )
         return;

      Nonzero<R>* e = m_elem + size();

      set_size( size() + n );
      while( n-- )
      {
         e->idx = *i++;
         e->val = *v++;
         e++;
      }
   }

   /// Append \p n nonzeros.
   void add(int n, const Nonzero<R> e[])
   {
      assert(n + size() <= max());

      if( n <= 0 )
         return;

      Nonzero<R>* ee = m_elem + size();

      set_size( size() + n );
      while( n-- )
         *ee++ = *e++;
   }

   /// Remove nonzeros \p n thru \p m.
   void remove(int n, int m)
   {
      assert(n <= m);
      assert(m < size());
      assert(n >= 0);

      ++m;

      int cpy = m - n;
      cpy = (size() - m >= cpy) ? cpy : size() - m;

      Nonzero<R>* e = &m_elem[size() - 1];
      Nonzero<R>* r = &m_elem[n];

      set_size(size() - cpy);
      do
      {
         *r++ = *e--;
      }
      while( --cpy );
   }

   /// Remove \p n 'th nonzero.
   void remove(int n)
   {
      assert(n >= 0);
      assert(n < size());

      set_size(size() - 1);
      m_elem[n] = m_elem[size()];
   }

   /// Remove all indices.
   void clear()
   {
      set_size(0);
   }

   /// Sort nonzeros to increasing indices.
   void sort()
   {
      if( m_elem != 0 )
      {
         Nonzero<R> dummy;
         Nonzero<R>* w;
         Nonzero<R>* l;
         Nonzero<R>* s = &(m_elem[0]);
         Nonzero<R>* e = s + size();

         for( l = s, w = s + 1; w < e; l = w, ++w )
         {
            if( l->idx > w->idx )
            {
               dummy = *w;

               do
               {
                  l[1] = *l;
                  if( l-- == s )
                     break;
               }
               while( l->idx > dummy.idx );

               l[1] = dummy;
            }
         }
      }
   }

   //@}

   // ------------------------------------------------------------------------------------------------------------------
   /**@name Arithmetic operations */
   //@{

   /// Maximum absolute value, i.e., infinity norm.
   R maxAbs() const
   {
      R maxi = 0.0;

      for( int i = size() - 1; i >= 0; --i )
      {
         if( abs(m_elem[i].val) > maxi )
            maxi = abs(m_elem[i].val);
      }

      assert(maxi >= 0.0);

      return maxi;
   }

   /// Minimum absolute value.
   R minAbs() const
   {
      R mini = infinity;

      for( int i = size() - 1; i >= 0; --i )
      {
         if( abs(m_elem[i].val) < mini )
            mini = abs(m_elem[i].val);
      }

      assert(mini >= 0.0);

      return mini;
   }

   /// Floating point approximation of euclidian norm (without any approximation guarantee).
   Real length() const
   {
      return sqrt((Real)length2());
   }

   /// Squared norm.
   R length2() const
   {
      R x = 0;
      int n = size();
      const Nonzero<R>* e = m_elem;

      while( n-- )
      {
         x += e->val * e->val;
         e++;
      }

      return x;
   }

   /// Scaling.
   SVectorBase<R>& operator*=(R x)
   {
      int n = size();
      Nonzero<R>* e = m_elem;

      while( n-- )
      {
         e->val *= x;
         e++;
      }

      return *this;
   }

   /// Inner product.
   R operator*(const VectorBase<R>& w) const;

   //@}

   // ------------------------------------------------------------------------------------------------------------------
   /**@name Constructions, destruction, and assignment */
   //@{

   /// Default constructor.
   /** The constructor expects one memory block where to store the nonzero elements. This must be passed to the
    *  constructor, where the \em number of Nonzero%s needs that fit into the memory must be given and a pointer to the
    *  beginning of the memory block. Once this memory has been passed, it shall not be modified until the SVectorBase
    *  is no longer used.
    */
   explicit SVectorBase<R>(int n = 0, Nonzero<R>* p_mem = 0)
   {
      setMem(n, p_mem);
   }

   /// Assignment operator.
   template < class S >
   SVectorBase<R>& operator=(const VectorBase<S>& vec);

   /// Assignment operator.
   SVectorBase<R>& operator=(const SVectorBase<R>& sv)
   {
      if( this != &sv )
      {
         assert(max() >= sv.size());

         int i = sv.size();
         Nonzero<R>* e = m_elem;
         const Nonzero<R>* s = sv.m_elem;

         while( i-- )
         {
            assert(e != 0);
            *e++ = *s++;
         }

         set_size(sv.size());
      }

      return *this;
   }

   /// Assignment operator.
   template < class S >
   SVectorBase<R>& operator=(const SVectorBase<S>& sv)
   {
      if( this != (SVectorBase<R>*)(&sv) )
      {
         assert(max() >= sv.size());

         int i = sv.size();
         Nonzero<R>* e = m_elem;
         const Nonzero<S>* s = sv.m_elem;

         while( i-- )
         {
            assert(e != 0);
            *e++ = *s++;
         }

         set_size(sv.size());
      }

      return *this;
   }

   /// Assignment operator.
   template < class S >
   SVectorBase<R>& operator=(const SSVectorBase<S>& sv);

   //@}

   // ------------------------------------------------------------------------------------------------------------------
   /**@name Memory */
   //@{

   /// get pointer to internal memory.
   Nonzero<R>* mem() const
   {
      return m_elem;
   }

   /// Set size of the vector.
   void set_size(int s)
   {
      assert(m_elem != 0 || s == 0);
      memused = s;
   }

   /// Set the maximum number of nonzeros in the vector.
   void set_max(int m)
   {
      assert(m_elem != 0 || m == 0);
      memsize = m;
   }

   /// Set the memory area where the nonzeros will be stored.
   void setMem(int n, Nonzero<R>* elmem)
   {
      assert(n >= 0);
      assert(n == 0 || elmem != 0);

      m_elem = elmem;
      set_size(0);
      set_max(n);
   }

   //@}

   // ------------------------------------------------------------------------------------------------------------------
   /**@name Utilities */
   //@{

   /// Consistency check.
   bool isConsistent() const
   {
#ifdef ENABLE_CONSISTENCY_CHECKS
      if( m_elem != 0 )
      {
         const int my_size = size();
         const int my_max = max();

         if( my_size < 0 || my_max < 0 || my_size > my_max )
            return MSGinconsistent("SVectorBase");

         for( int i = 1; i < my_size; ++i )
         {
            for( int j = 0; j < i; ++j )
            {
               // allow trailing zeros
               if( m_elem[i].idx == m_elem[j].idx && m_elem[i].val != 0.0 )
                  return MSGinconsistent("SVectorBase");
            }
         }
      }
#endif

      return true;
   }

   //@}
};

} // namespace soplex
#endif // _SVECTORBASE_H_

//-----------------------------------------------------------------------------
//Emacs Local Variables:
//Emacs mode:c++
//Emacs c-basic-offset:3
//Emacs tab-width:8
//Emacs indent-tabs-mode:nil
//Emacs End:
//-----------------------------------------------------------------------------