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
#pragma ident "@(#) $Id: updatevector.h,v 1.3 2001/11/07 17:31:26 bzfbleya Exp $"


#ifndef _UPDATEVECTOR_H_
#define _UPDATEVECTOR_H_

#include <assert.h>


#include "dvector.h"
#include "ssvector.h"

namespace soplex
{




extern void UpdateUpdateVector(double*, double, int, const int*, const double*);

//@ ----------------------------------------------------------------------------

/*  \Section{Class Declaration}
    The datastructur of #UpdateVector# is straightforward. It adds to the
    baseclass #DVector# another #SSVector thedelta# for the update vector
    $\delta$ and a #double theval# for the update value $\alpha$.
 */

/** vector with updates.
    In many algorithms vectors are updated in every iteration, by adding a
    multiple of another vector to it, i.e., given a vector $x$, a scalar
    $\alpha$ and another vector $\delta$, the update to $x$ constists of
    substituting it by $x \leftarrow x + \alpha\cdot\delta$.
 
    While the update itself can easily be expressed with methods of #Vector#,
    it is often desirable to save the last update vector $\delta$ and value
    $\alpha$. This is provided by class #UpdateVector#.
 
    #UpdateVector#s are derived from #DVector# and provide additional methods
    for saving and setting the multiplicator $\alpha$ and the update vector
    $\delta$. Further, it allows for efficient sparse updates, by providing an
    #IdxSet# idx containing the nonzero indeces of $\delta$.
 */
class UpdateVector : public DVector
{
   double theval;
   SSVector thedelta;

public:
   ///
   double& value()
   {
      return theval;
   }
   /// update multiplicator $\alpha$.
   double value() const
   {
      return theval;
   }

   ///
   SSVector& delta()
   {
      return thedelta;
   }
   /// update vector $\delta$.
   const SSVector& delta() const
   {
      return thedelta;
   }

   /// nonzero indeces of $\delta$.
   const IdxSet& idx() const
   {
      return thedelta.indices();
   }

   /** Update vector with $\alpha \cdot \delta$.
    *  Add #value() * delta()# to the #UpdateVector#. Only the indeces set
    *  in #idx()# are affected. For all other indeces, #delta()# is asumed
    *  to be 0.
    */
   void update()
   {
      multAdd(theval, thedelta);
   }

   /// clear vector and update vector.
   void clear()
   {
      DVector::clear();
      clearUpdate();
   }

   /// clear $\delta$, $\alpha$ and #idx()# using #idx()#.
   void clearUpdate()
   {
      thedelta.clear();
      theval = 0;
   }


   /// reset dimension.
   void reDim(int newdim)
   {
      DVector::reDim(newdim);
      thedelta.reDim(newdim);
   }

   ///
   UpdateVector& operator=(const DVector& rhs)
   {
      DVector::operator=(rhs);
      return *this;
   }

   ///
   UpdateVector& operator=(const UpdateVector& rhs);

   /// default constructor.
   UpdateVector(int p_dim /*=0*/, double p_eps /*=1e-16*/)
      : DVector (p_dim),
         theval (0),
         thedelta(p_dim, p_eps)
   { }

   ///

   int isConsistent() const;
};


} // namespace soplex
#endif // _UPDATEVECTOR_H_

//-----------------------------------------------------------------------------
//Emacs Local Variables:
//Emacs mode:c++
//Emacs c-basic-offset:3
//Emacs tab-width:8
//Emacs indent-tabs-mode:nil
//Emacs End:
//-----------------------------------------------------------------------------
