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
#pragma ident "@(#) $Id: lpcol.h,v 1.2 2001/11/06 23:31:02 bzfkocht Exp $"


#ifndef _LPCOL_H_
#define _LPCOL_H_

//@ ----------------------------------------------------------------------------
/*      \Section{Imports}
    Import required system include files
 */
#include <assert.h>


/*  and class header files
 */

#include "dsvector.h"

namespace soplex
{






//@ ----------------------------------------------------------------------------
/* \Section{Class Declaration}
 */
/** LP column.
    Class #LPCol# provides a datatype for storing the column of an LP a the
    form similar to
        \[
        \begin{array}{rl}
            \hbox{max}  & c^T x         \\
            \hbox{s.t.} & Ax \le b      \\
                        & l \le x \le u
        \end{array}
        \]
    Hence, an #LPCol# consists of an objective value, a column #DSVector# and
    an upper and lower bound to the corresponding variable, which may include
    $\pm\infty$. However, it depends on the LP code to use, what values are
    actually treated as $\infty$.
 */
class LPCol
{
private:
   double up, low, object;
   DSVector vec;
public:
   ///
   double obj() const
   {
      return object;
   }
   /// access objective value.
   double& obj()
   {
      return object;
   }

   ///
   double upper() const
   {
      return up;
   }
   /// access upper bound.
   double& upper()
   {
      return up;
   }

   ///
   double lower() const
   {
      return low;
   }
   /// access lower bound.
   double& lower()
   {
      return low;
   }

   ///
   const SVector& colVector() const
   {
      return vec;
   }
   /// access constraint column vector.
   DSVector& colVector()
   {
      return vec;
   }

   ///
   LPCol(const LPCol& old)
      : up(old.up), low(old.low), object(old.object), vec(old.vec)
   { }

   /** Default Constructor.
       Construct #LPCol# with a column vector ready for taking #defDim#
       nonzeros.
    */
   LPCol(int defDim = 0)
      : up(1e+300), low(0), object(0), vec(defDim)
   { }

   /** Initializing Constructor.
       Construct #LPCol# with the given objective value #obj#, a column
       vector #vec#, upper bound #upper# and lower bound #lower#.
    */
   LPCol(double obj, const SVector& vector, double upper, double lower)
      : up(upper), low(lower), object(obj), vec(vector)
   { }

   ///

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
