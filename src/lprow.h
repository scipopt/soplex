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
#pragma ident "@(#) $Id: lprow.h,v 1.1 2001/11/06 16:18:32 bzfkocht Exp $"


#ifndef _LPROW_H_
#define _LPROW_H_

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

/** (In)equality for LPs.
    Class #LPRow# provides constraints for linear programs in the form
    \[
                        l \le a^Tx \le r,
    \]
    where $a$ is a \Ref{DSVector}. $l$ is referred to as {\em left hand side},
    $r$ as {\rm right hand side} and $a$ as {\em row vector} or the constraint
    vector. $l$ and $r$ may also take values $\pm$#infinity#. This static member
    is predefined, but may be overridden to meet the needs of the LP solver to
    be used.
 
    #LPRow#s allow to specify regular inequalities of the form 
    \[
                            a^Tx \sim \alpha,
    \]
    where $\sim$ can take any value of $\{\le, =, \ge\}$, by setting #rhs# and
    #lhs# to the same value or setting one of them to $\infty$.
 
    Since constraints in the regular form occur often, #LPRow#s offers methods
    #type()# and #value()# for retreiving $\sim$ and $\alpha$ of an #LPRow# in
    this form, respectively. Also, a constructor for #LPRow#s given in regular
    form is provided.
 */
class LPRow
{
private:
   double left, right;
   DSVector vec;

public:
   /// values #>= infinity# are treated as $\infty$.
   static double infinity;

   /** (In)Equality of an LP row.
       #LPRow#s may be of one of the above #Type#s. This datatype may be
       used for constructing new #LPRow#s in the regular form.
    */
   enum Type
   {                           /// $a^Tx \le \alpha$.
      LESS_EQUAL,              /// $a^Tx = \alpha$.
      EQUAL,                   /// $a^Tx \ge \alpha$.
      GREATER_EQUAL,           /// $\lambda \le a^Tx \le \rho$.
      RANGE
   };

   /**@name Inquiry */
   //@{
   ///
   Type type() const;
   /** Right-hand side value of (in)equality.
       This method returns $\alpha$ for a #LPRow# in regular form.
       However, #value()# may only be called for #LPRow#s with
       #type() != RANGE#.
    */

   /// set type of (in)equality
   void setType(Type type);

   double value() const;

   ///
   double lhs() const
   {
      return left;
   }
   /// access left hand side value.
   double& lhs()
   {
      return left;
   }

   ///
   double rhs() const
   {
      return right;
   }
   /// access right hand side value.
   double& rhs()
   {
      return right;
   }

   ///
   const SVector& rowVector() const
   {
      return vec;
   }
   /// access constraint rowVector.
   DSVector& rowVector()
   {
      return vec;
   }
   //@}

   /**@name Miscellaneous */
   //@{
   /** Construct #LPRow# with a vector ready to hold #defDim# nonzeros
    */
   LPRow(int defDim = 0)
      : left(0), right(infinity), vec(defDim)
   {}

   ///

   LPRow(const LPRow& row)
      : left(row.left), right(row.right), vec(row.vec)
   {}

   /** Construct #LPRow# with the given left-hand side, right-hand side
       and rowVector.
    */
   LPRow(double lhs, const SVector& rowVector, double rhs)
      : left(lhs), right(rhs), vec(rowVector)
   {}

   /** Construct #LPRow# from passed #rowVector#, #type# and #value#
    */
   LPRow(const SVector& rowVector, Type type, double value);

   /// check consistency.
   int isConsistent() const
   {
      return vec.isConsistent();
   }
   //@}
};

} // namespace soplex
#endif // _LPROW_H_

//-----------------------------------------------------------------------------
//Emacs Local Variables:
//Emacs c-basic-offset:3
//Emacs tab-width:8
//Emacs indent-tabs-mode:nil
//Emacs End:
//-----------------------------------------------------------------------------
