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
#pragma ident "@(#) $Id: lpcol.h,v 1.4 2001/11/17 22:15:58 bzfkocht Exp $"

/**@file  lpcol.h
 * @brief LP column.
 */
#ifndef _LPCOL_H_
#define _LPCOL_H_

#include <assert.h>

#include "dsvector.h"

namespace soplex
{
/**@brief   LP column.
   @ingroup Algo

   Class LPCol provides a datatype for storing the column of an LP a the
   form similar to
   \f[
      \begin{array}{rl}
         \hbox{max}  & c^T x         \\
         \hbox{s.t.} & Ax \le b      \\
                     & l \le x \le u
      \end{array}
   \f]
   Hence, an LPCol consists of an objective value, a column DSVector and
   an upper and lower bound to the corresponding variable, which may include
   \f$\pm\infty\f$. However, it depends on the LP code to use, what values are
   actually treated as \f$\infty\f$.
 */
class LPCol
{
private:
   double   up;
   double   low;
   double   object;
   DSVector vec;

public:
   /// get objective value.
   double obj() const
   {
      return object;
   }
   /// access objective value.
   double& obj()
   {
      return object;
   }

   /// get upper bound.
   double upper() const
   {
      return up;
   }
   /// access upper bound.
   double& upper()
   {
      return up;
   }

   /// get lower bound.
   double lower() const
   {
      return low;
   }
   /// access lower bound.
   double& lower()
   {
      return low;
   }

   /// get constraint column vector.
   const SVector& colVector() const
   {
      return vec;
   }
   /// access constraint column vector.
   DSVector& colVector()
   {
      return vec;
   }

   /// copy constructor.
   LPCol(const LPCol& old)
      : up(old.up), low(old.low), object(old.object), vec(old.vec)
   {}

   /// default constructor.
   /** Construct LPCol with a column vector ready for taking \p defDim
    *  nonzeros.
    */
   explicit LPCol(int defDim = 0)
      : up(1e+300), low(0), object(0), vec(defDim)
   {}

   /// initializing constructor.
   /*  Construct LPCol with the given objective value \p obj, a column
    *  %vector \p vec, upper bound \p upper and lower bound \p lower.
    */
   LPCol(double pobj, const SVector& pvector, double pupper, double plower)
      : up(pupper), low(plower), object(pobj), vec(pvector)
   {}

   /// check consistency.
   int isConsistent() const
   {
      return vec.isConsistent();
   }
};
} // namespace soplex
#endif // _LPCOL_H_

//-----------------------------------------------------------------------------
//Emacs Local Variables:
//Emacs mode:c++
//Emacs c-basic-offset:3
//Emacs tab-width:8
//Emacs indent-tabs-mode:nil
//Emacs End:
//-----------------------------------------------------------------------------
