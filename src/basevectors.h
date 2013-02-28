/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the class library                   */
/*       SoPlex --- the Sequential object-oriented simPlex.                  */
/*                                                                           */
/*    Copyright (C) 1996      Roland Wunderling                              */
/*                  1996-2013 Konrad-Zuse-Zentrum                            */
/*                            fuer Informationstechnik Berlin                */
/*                                                                           */
/*  SoPlex is distributed under the terms of the ZIB Academic Licence.       */
/*                                                                           */
/*  You should have received a copy of the ZIB Academic License              */
/*  along with SoPlex; see the file COPYING. If not email to soplex@zib.de.  */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/**@file  vectors.h
 * @brief Collection of dense, sparse, and semi-sparse vectors.
 */
#ifndef _VECTORS_H_
#define _VECTORS_H_

#include "spxdefines.h"
#include "rational.h"
#include "vectorbase.h"
#include "dvectorbase.h"
#include "ssvectorbase.h"
#include "svectorbase.h"
#include "dsvectorbase.h"
#include "svsetbase.h"

namespace soplex
{

// ---------------------------------------------------------------------------------------------------------------------
//  Methods of VectorBase
// ---------------------------------------------------------------------------------------------------------------------

/// Assignment operator.
/** Assigning an SVectorBase to a VectorBase using operator=() will set all values to 0 except the nonzeros of \p vec.
 *  This is diffent in method assign().
 */
template < class R >
template < class S >
inline VectorBase<R>& VectorBase<R>::operator=(const SVectorBase<S>& vec)
{
   clear();

   for( int i = 0; i < vec.size(); i++ )
   {
      assert(vec.index(i) < dim());
      val[vec.index(i)] = vec.value(i);
   }

   assert(isConsistent());

   return *this;
}

/// Assign values of \p vec.
/** Assigns all nonzeros of \p vec to the vector.  All other values remain unchanged. */
template < class R >
template < class S >
VectorBase<R>& VectorBase<R>::assign(const SVectorBase<S>& vec)
{
   for( int i = vec.size(); i-- > 0; )
      val[vec.index(i)] = vec.value(i);

   assert(isConsistent());

   return *this;
}

/// Assignment operator.
/** Assigning an SSVectorBase to a VectorBase using operator=() will set all values to 0 except the nonzeros of \p vec.
 *  This is diffent in method assign().
 */
template < class R >
template < class S >
VectorBase<R>& VectorBase<R>::operator=(const SSVectorBase<S>& vec)
{
   if( vec.isSetup() )
   {
      clear();
      assign(vec);
   }
   else
      operator=(static_cast<const VectorBase<R>&>(vec));

   return *this;
}

/// Assign values of \p vec.
/** Assigns all nonzeros of \p vec to the vector.  All other values remain unchanged. */
template < class R >
template < class S >
inline VectorBase<R>& VectorBase<R>::assign(const SSVectorBase<S>& vec)
{
   assert(vec.dim() <= dim());

   if (vec.isSetup())
   {
      const int* idx = vec.indexMem();

      for(int i = vec.size(); i > 0; i--)
      {
         val[*idx] = vec.val[*idx];
         idx++;
      }
   }
   else
      operator=(static_cast<const VectorBase<R>&>(vec));

   return *this;
}

/// Addition.
template < class R >
template < class S >
VectorBase<R>& VectorBase<R>::operator+=(const SVectorBase<S>& vec)
{
   for( int i = 0; i < vec.size(); i++ )
   {
      assert(vec.index(i) >= 0);
      assert(vec.index(i) < dim());
      val[vec.index(i)] += vec.value(i);
   }

   return *this;
}

/// Subtraction.
template < class R >
template < class S >
VectorBase<R>& VectorBase<R>::operator-=(const SVectorBase<S>& vec)
{
   for( int i = 0; i < vec.size() ; i++ )
   {
      assert(vec.index(i) >= 0);
      assert(vec.index(i) < dim());
      val[vec.index(i)] -= vec.value(i);
   }

   return *this;
}

/// Subtraction.
template < class R >
template < class S >
VectorBase<R>& VectorBase<R>::operator-=(const SSVectorBase<S>& vec)
{
   assert(dim() == vec.dim());

   for( int i = 0; i < dim(); i++ )
      val[i] -= vec[i];

   return *this;
}

/// Inner product.
template < class R >
R VectorBase<R>::operator*(const SVectorBase<R>& vec) const
{
   assert(dim() >= vec.dim());

   R x = 0;

   for( int i = vec.size(); i > 0; i-- )
      x += val[vec.index(i)] * vec.value(i);

   return x;
}

/// Inner product.
template < class R >
R VectorBase<R>::operator*(const SSVectorBase<R>& vec) const
{
   assert(dim() == vec.dim());

   if( vec.isSetup() )
   {
      const int* idx = vec.indexMem();

      R x = 0;

      for( int i = vec.size(); i > 0; i-- )
      {
         x += val[*idx] * vec.val[*idx];
         idx++;
      }

      return x;
   }
   else
      return operator*(static_cast<const VectorBase<R>&>(vec));
}

/// Addition of scaled vector.
template < class R >
template < class S, class T >
VectorBase<R>& VectorBase<R>::multAdd(S x, const SVectorBase<T>& vec)
{
   for( int i = 0; i < vec.size(); i++ )
   {
      assert(vec.index(i) < dim());
      val[vec.index(i)] += x * vec.value(i);
   }

   return *this;
}

/// Addition of scaled vector.
template < class R >
template < class S, class T >
inline VectorBase<R>& VectorBase<R>::multAdd(S x, const SSVectorBase<T>& vec)
{
   assert(vec.dim() <= dim());

   if( vec.isSetup() )
   {
      const int* idx = vec.indexMem();

      for( int i = 0; i < vec.size(); i++ )
         val[idx[i]] += x * vec[idx[i]];
   }
   else
   {
      assert(vec.dim() == dim());

      for( int i = 0; i < dim(); i++ )
         val[i] += x * vec.val[i];
   }

   return *this;
}

// ---------------------------------------------------------------------------------------------------------------------
//  Methods of DVectorBase
// ---------------------------------------------------------------------------------------------------------------------

/// Assignment operator.
template < class R >
template < class S >
DVectorBase<R>& DVectorBase<R>::operator=(const SVectorBase<S>& vec)
{
   if( vec.dim() != VectorBase<R>::dim() )
      reDim(vec.dim());

   VectorBase<R>::operator=(vec);

   assert(isConsistent());

   return *this;
}

// ---------------------------------------------------------------------------------------------------------------------
// Methods of SSVectorBase
// ---------------------------------------------------------------------------------------------------------------------

/// Addition.
template < class R >
template < class S >
SSVectorBase<R>& SSVectorBase<R>::operator+=(const SVectorBase<S>& vec)
{
   VectorBase<R>::operator+=(vec);

   if( isSetup() )
   {
      setupStatus = false;
      setup();
   }

   return *this;
}

/// Subtraction.
template < class R >
template < class S >
SSVectorBase<R>& SSVectorBase<R>::operator-=(const SVectorBase<S>& vec)
{
   VectorBase<R>::operator-=(vec);

   if( isSetup() )
   {
      setupStatus = false;
      setup();
   }

   return *this;
}

/// Addition of a scaled vector.
///@todo SSVectorBase::multAdd() should be rewritten without pointer arithmetic.
template < class R >
template < class S, class T >
SSVectorBase<R>& SSVectorBase<R>::multAdd(S xx, const SVectorBase<T>& vec)
{
   if( isSetup() )
   {
      R* v = VectorBase<R>::val;
      R x;
      int adjust = 0;
      int j;

      for( int i = vec.size() - 1; i >= 0; --i )
      {
         j = vec.index(i);

         if( v[j] )
         {
            x = v[j] + xx * vec.value(i);
            if( isNotZero(x, epsilon) )
               v[j] = x;
            else
            {
               adjust = 1;
               v[j] = MARKER;
            }
         }
         else
         {
            x = xx * vec.value(i);
            if( isNotZero(x, epsilon) )
            {
               v[j] = x;
               addIdx(j);
            }
         }
      }

      if( adjust )
      {
         int* iptr = idx;
         int* iiptr = idx;
         int* endptr = idx + num;

         for( ; iptr < endptr; ++iptr )
         {
            x = v[*iptr];
            if( isNotZero(x, epsilon) )
               *iiptr++ = *iptr;
            else
               v[*iptr] = 0;
         }

         num = int(iiptr - idx);
      }
   }
   else
      VectorBase<R>::multAdd(xx, vec);

   assert(isConsistent());

   return *this;
}

/// Assigns \f$x^T \cdot A\f$ to SSVectorBase.
template < class R >
template < class S, class T >
SSVectorBase<R>& SSVectorBase<R>::assign2product(const SSVectorBase<S>& x, const SVSetBase<T>& A)
{
   assert(A.num() == dim());

   R y;

   clear();

   for( int i = dim(); i-- > 0; )
   {
      y = A[i] * x;

      if( isNotZero(y, epsilon) )
      {
         VectorBase<R>::val[i] = y;
         IdxSet::addIdx(i);
      }
   }

   assert(isConsistent());

   return *this;
}

/// Assigns SSVectorBase to \f$A \cdot x\f$ for a setup \p x.
#define shortProductFactor 0.5
template < class R >
template < class S, class T >
SSVectorBase<R>& SSVectorBase<R>::assign2product4setup(const SVSetBase<S>& A, const SSVectorBase<T>& x)
{
   assert(A.num() == x.dim());
   assert(x.isSetup());
   clear();

   if( x.size() == 1 )
   {
      assign2product1(A, x);
      setupStatus = true;
   }
   else if( isSetup() && (double(x.size()) * A.memSize() <= shortProductFactor * dim() * A.num()) )
   {
      assign2productShort(A, x);
      setupStatus = true;
   }
   else
   {
      assign2productFull(A, x);
      setupStatus = false;
   }

   assert(isConsistent());

   return *this;
}

/// Assignment helper.
template < class R >
template < class S, class T >
SSVectorBase<R>& SSVectorBase<R>::assign2product1(const SVSetBase<S>& A, const SSVectorBase<T>& x)
{
   assert(x.isSetup());
   assert(x.size() == 1);

   // get the nonzero value of x and the corresponding vector in A:
   const int nzidx = x.idx[0];
   const T nzval = x.val[nzidx];
   const SVectorBase<S>& Ai = A[nzidx];

   // compute A[nzidx] * nzval:
   if( isZero(nzval, epsilon) || Ai.size() == 0 )
      clear();    // this := zero vector
   else
   {
      num = Ai.size();
      for( register int j = 0; j < num; j++ )
      {
         const Nonzero<S>& Aij = Ai.element(j);
         idx[j] = Aij.idx;
         VectorBase<R>::val[Aij.idx] = nzval * Aij.val;
      }
   }

   assert(isConsistent());

   return *this;
}

/// Assignment helper.
template < class R >
template < class S, class T >
SSVectorBase<R>& SSVectorBase<R>::assign2productShort(const SVSetBase<S>& A, const SSVectorBase<T>& x)
{
   assert(x.isSetup());

   if( x.size() == 0 ) // x can be setup but have size 0 => this := zero vector
   {
      clear();
      return *this;
   }

   // compute x[0] * A[0]
   int curidx = x.idx[0];
   const T x0 = x.val[curidx];
   const SVectorBase<S>& A0 = A[curidx];
   int nonzero_idx = 0;

   num = A0.size();
   if( isZero(x0, epsilon) || num == 0 )
   {
      // A[0] == 0 or x[0] == 0 => this := zero vector
      clear();
   }
   else
   {
      for( register int j = 0; j < num; ++j )
      {
         const Nonzero<S>& elt = A0.element(j);
         const R product = x0 * elt.val;

         // store the value in any case
         idx[nonzero_idx] = elt.idx;
         VectorBase<R>::val[elt.idx] = product;

         // count only non-zero values; not 'isNotZero(product, epsilon)'
         if( product != 0 )
            ++nonzero_idx;
      }
   }

   // Compute the other x[i] * A[i] and add them to the existing vector.
   for( register int i = 1; i < x.size(); ++i )
   {
      curidx = x.idx[i];
      const T xi     = x.val[curidx];
      const SVectorBase<S>& Ai = A[curidx];

      // If A[i] == 0 or x[i] == 0, do nothing.
      if ( isNotZero(xi, epsilon) || Ai.size() == 0 )
      {
         // Compute x[i] * A[i] and add it to the existing vector.
         for( register int j = 0; j < Ai.size(); ++j )
         {
            const Nonzero<S>& elt = Ai.element(j);
            idx[nonzero_idx] = elt.idx;
            R oldval  = VectorBase<R>::val[elt.idx];

            // An old value of exactly 0 means the position is still unused.
            // It will be used now (either by a new nonzero or by a MARKER),
            // so increase the counter. If oldval != 0, we just
            // change an existing NZ-element, so don't increase the counter.
            if( oldval == 0 )
               ++nonzero_idx;

            // Add the current product x[i] * A[i][j]; if oldval was
            // MARKER before, it does not hurt because MARKER is really small.
            oldval += xi * elt.val;

            // If the new value is exactly 0, mark the index as used
            // by setting a value which is nearly 0; otherwise, store
            // the value. Values below epsilon will be removed later.
            if( oldval == 0 )
               VectorBase<R>::val[elt.idx] = MARKER;
            else
               VectorBase<R>::val[elt.idx] = oldval;
         }
      }
   }

   // Clean up by shifting all nonzeros (w.r.t. epsilon) to the front of idx,
   // zeroing all values which are nearly 0, and setting #num# appropriately.
   int nz_counter = 0;

   for( register int i = 0; i < nonzero_idx; ++i )
   {
      curidx = idx[i];

      if( isZero( VectorBase<R>::val[curidx], epsilon ) )
         VectorBase<R>::val[curidx] = 0;
      else
      {
         idx[nz_counter] = curidx;
         ++nz_counter;
      }

      num = nz_counter;
   }

   assert(isConsistent());

   return *this;
}

/// Assignment helper.
template < class R >
template < class S, class T >
SSVectorBase<R>& SSVectorBase<R>::assign2productFull(const SVSetBase<S>& A, const SSVectorBase<T>& x)
{
   assert(x.isSetup());

   if( x.size() == 0 ) // x can be setup but have size 0 => this := zero vector
   {
      clear();
      return *this;
   }

   bool A_is_zero = true;

   for( int i = 0; i < x.size(); ++i )
   {
      const int curidx = x.idx[i];
      const T xi = x.val[curidx];
      const SVectorBase<S>& Ai = A[curidx];

      if( A_is_zero && Ai.size() > 0 )
         A_is_zero = false;

      for( register int j = 0; j < Ai.size(); ++j )
      {
         const Nonzero<S>& elt = Ai.element(j);
         VectorBase<R>::val[elt.idx] += xi * elt.val;
      }
   }

   if( A_is_zero )
      clear(); // case x != 0 but A == 0

   return *this;
}

/// Assigns SSVectorBase to \f$A \cdot x\f$ thereby setting up \p x.
template < class R >
template < class S, class T >
SSVectorBase<R>& SSVectorBase<R>::assign2productAndSetup(const SVSetBase<S>& A, SSVectorBase<T>& x)
{
   if( x.isSetup() )
      return assign2product4setup(A, x);

   if( x.dim() == 0 )
   { // x == 0 => this := zero vector
      clear();
      x.num = 0;
   }
   else
   {
      // x is not setup, so walk through its value vector
      int nzcount = 0;
      int end = x.dim();

      for( int i = 0; i < end; ++i )
      {
         // advance to the next element != 0
         T& xval = x.val[i];

         if( xval != 0 )
         {
            // If x[i] is really nonzero, compute A[i] * x[i] and adapt x.idx,
            // otherwise set x[i] to 0.
            if( isNotZero(xval, epsilon) )
            {
               const SVectorBase<S>& Ai = A[i];
               x.idx[ nzcount++ ] = i;

               for( int j = Ai.size() - 1; j >= 0; --j )
               {
                  const Nonzero<S>& elt = Ai.element(j);
                  VectorBase<R>::val[elt.idx] += xval * elt.val;
               }
            }
            else
               xval = 0;
         }
      }

      x.num = nzcount;
   }

   x.setupStatus = true;
   setupStatus = false;

   assert(isConsistent());

   return *this;
}

/// Assigns only the elements of \p rhs.
template < class R >
template < class S >
SSVectorBase<R>& SSVectorBase<R>::assign(const SVectorBase<S>& rhs)
{
   assert(rhs.dim() <= VectorBase<R>::dim());

   num = 0;

   for( int i = 0; i < rhs.size(); ++i )
   {
      int k = rhs.index(i);
      S v = rhs.value(i);

      if( isZero(v, epsilon) )
         VectorBase<R>::val[k] = 0;
      else
      {
         VectorBase<R>::val[k] = v;
         idx[num++] = k;
      }
   }

   setupStatus = true;

   assert(isConsistent());

   return *this;
}

/// Assignment operator.
template < class R >
template < class S >
SSVectorBase<R>& SSVectorBase<R>::operator=(const SVectorBase<S>& rhs)
{
   clear();

   return assign(rhs);
}

// ---------------------------------------------------------------------------------------------------------------------
//  Methods of SVectorBase
// ---------------------------------------------------------------------------------------------------------------------

/// Assignment operator.
template < class R >
template < class S >
SVectorBase<R>& SVectorBase<R>::operator=(const VectorBase<S>& vec)
{
   int n = 0;
   int i = vec.dim();
   Nonzero<R> *e = m_elem;

   clear();

   while( i-- )
   {
      if( vec[i] )
      {
         assert(n < max());

         e->idx = i;
         e->val = vec[i];
         ++e;
         ++n;
      }
   }

   set_size(n);

   return *this;
}

/// Assignment operator.
template < class R >
template < class S >
SVectorBase<R>& SVectorBase<R>::operator=(const SSVectorBase<S>& sv)
{
   assert(max() >= sv.size());

   set_size(sv.size());

   int i = size();
   Nonzero<R> *e = m_elem;

   while( i-- )
   {
      e->idx = sv.index(i);
      e->val = sv[e->idx];
      ++e;
   }

   return *this;
}

/// Inner product.
template < class R >
R SVectorBase<R>::operator*(const VectorBase<R>& w) const
{
   R x = 0;
   int n = size();
   Nonzero<R>* e = m_elem;

   while( n-- )
   {
      x += e->val * w[e->idx];
      e++;
   }

   return x;
}

// ---------------------------------------------------------------------------------------------------------------------
//  Methods of DSVectorBase
// ---------------------------------------------------------------------------------------------------------------------

/// Copy constructor.
template < class R >
template < class S >
DSVectorBase<R>::DSVectorBase(const VectorBase<S>& vec)
   : theelem(0)
{
   allocMem((vec.dim() < 1) ? 2 : vec.dim());
   *this = vec;

   assert(isConsistent());
}

/// Copy constructor.
template < class R >
template < class S >
DSVectorBase<R>::DSVectorBase(const SSVectorBase<S>& old)
   : theelem(0)
{
   allocMem(old.size() < 1 ? 2 : old.size());
   SVectorBase<R>::operator=(old);

   assert(isConsistent());
}

/// Assignment operator.
template < class R >
template < class S >
DSVectorBase<R>& DSVectorBase<R>::operator=(const VectorBase<S>& vec)
{
   assert(this != (DSVectorBase<R>*)(&vec));

   SVectorBase<R>::clear();
   setMax(vec.dim());
   SVectorBase<R>::operator=(vec);

   assert(isConsistent());

   return *this;
}

/// Assignment operator.
template < class R >
template < class S >
DSVectorBase<R>& DSVectorBase<R>::operator=(const SSVectorBase<S>& vec)
{
   assert(this != &vec);

   SVectorBase<R>::clear();
   makeMem(vec.size());
   SVectorBase<R>::operator=(vec);

   return *this;
}

// ---------------------------------------------------------------------------------------------------------------------
//  Operators
// ---------------------------------------------------------------------------------------------------------------------

/// Output operator.
template < class R >
std::ostream& operator<<(std::ostream& s, const VectorBase<R>& vec)
{
   int i;

   s << '(';

   for( i = 0; i < vec.dim() - 1; ++i )
      s << vec[i] << ", ";

   s << vec[i] << ')';

   return s;
}

/// Negation.
template < class R >
DVectorBase<R> operator-(const VectorBase<R>& vec)
{
   DVectorBase<R> res(vec.dim());

   for( int i = 0; i < res.dim(); ++i )
      res[i] = -vec[i];

   return res;
}

/// Addition.
template < class R >
DVectorBase<R> operator+(const VectorBase<R>& v, const VectorBase<R>& w)
{
   assert(v.dim() == w.dim());

   DVectorBase<R> res(v.dim());

   for( int i = 0; i < res.dim(); ++i )
      res[i] = v[i] + w[i];

   return res;
}

/// Addition.
template < class R >
DVectorBase<R> operator+(const VectorBase<R>& v, const SVectorBase<R>& w)
{
   DVectorBase<R> res(v);

   res += w;

   return res;
}

/// Addition.
template < class R >
DVectorBase<R> operator+(const SVectorBase<R>& v, const VectorBase<R>& w)
{
   return w + v;
}

/// Subtraction.
template < class R >
DVectorBase<R> operator-(const VectorBase<R>& v, const VectorBase<R>& w)
{
   assert(v.dim() == w.dim());

   DVectorBase<R> res(v.dim());

   for( int i = 0; i < res.dim(); ++i )
      res[i] = v[i] - w[i];

   return res;
}

/// Subtraction.
template < class R >
DVectorBase<R> operator-(const VectorBase<R>& v, const SVectorBase<R>& w)
{
   DVectorBase<R> res(v);

   res -= w;

   return res;
}

/// Subtraction.
template < class R >
DVectorBase<R> operator-(const SVectorBase<R>& v, const VectorBase<R>& w)
{
   DVectorBase<R> res(w.dim());

   for( int i = 0; i < res.dim(); ++i )
      res[i] = -w[i];

   res += v;

   return res;
}

/// Scaling.
template < class R >
DVectorBase<R> operator*(const VectorBase<R>& v, R x)
{
   DVectorBase<R> res(v.dim());

   for( int i = 0; i < res.dim(); ++i )
      res[i] = x * v[i];

   return res;
}

/// Scaling.
template < class R >
DVectorBase<R> operator*(R x, const VectorBase<R>& v)
{
   return v * x;
}

/// Input operator.
template < class R >
std::istream& operator>>(std::istream& s, DVectorBase<R>& vec)
{
   char c;
   R val;
   int i = 0;

   while( s.get(c).good() )
   {
      if( c != ' ' && c != '\t' && c != '\n' )
         break;
   }

   if( c != '(' )
      s.putback(c);
   else
   {
      do
      {
         s >> val;

         if( i >= vec.dim() - 1 )
            vec.reDim(i + 16);
         vec[i++] = val;

         while( s.get(c).good() )
         {
            if( c != ' ' && c != '\t' && c != '\n' )
               break;
         }

         if( c != ',' )
         {
            if (c != ')')
               s.putback(c);
            break;
         }
      }
      while( s.good() );
   }

   vec.reDim(i);

   return s;
}

/// Output operator.
template < class R >
std::ostream& operator<<(std::ostream& os, const SVectorBase<R>& v)
{
   for( int i = 0, j = 0; i < v.size(); ++i )
   {
      if( j )
      {
         if( v.value(i) < 0 )
            os << " - " << -v.value(i);
         else
            os << " + " << v.value(i);
      }
      else
         os << v.value(i);

      os << " x" << v.index(i);
      j = 1;

      if( (i + 1) % 4 == 0 )
         os << "\n\t";
   }

   return os;
}

/// Output operator.
template < class R >
std::ostream& operator<<(std::ostream& os, const SVSetBase<R>& s)
{
   for( int i = 0; i < s.num(); i++ )
      os << s[i] << "\n";

   return os;
}

// ---------------------------------------------------------------------------------------------------------------------
//  Explicit instantiations
// ---------------------------------------------------------------------------------------------------------------------

template class VectorBase < Real >;
template class DVectorBase < Real >;
template class SSVectorBase < Real >;
template class SVectorBase < Real >;
template class DSVectorBase < Real >;
template class SVSetBase < Real >;

#ifdef SOPLEX_WITH_GMP
template class VectorBase < Rational >;
template class DVectorBase < Rational >;
template class SVectorBase < Rational >;
template class DSVectorBase < Rational >;
template class SSVectorBase < Rational >;
template class SVSetBase < Rational >;
#endif

}

#endif // _VECTORS_H_

// ---------------------------------------------------------------------------------------------------------------------
//Emacs Local Variables:
//Emacs mode:c++
//Emacs c-basic-offset:3
//Emacs tab-width:8
//Emacs indent-tabs-mode:nil
//Emacs End:
// ---------------------------------------------------------------------------------------------------------------------