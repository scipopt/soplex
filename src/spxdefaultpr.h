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
#pragma ident "@(#) $Id: spxdefaultpr.h,v 1.5 2001/12/25 14:25:55 bzfkocht Exp $"


/**@file  spxdefaultpr.h
 * @brief Default pricer.
 */
#ifndef _SPXDEFAULTPR_H_
#define _SPXDEFAULTPR_H_

#include <assert.h>

#include "spxpricer.h"

namespace soplex
{

/**@brief   Default pricer.
   @ingroup Algo

   Class #SPxDefaultPR is an implementation class for #SPxPricer implementing
   Dantzig's the default pricing strategy, i.e. maximal/minimal reduced cost or
   maximal violated constraint.

   See #SPxPricer for a class documentation.
*/
class SPxDefaultPR : public SPxPricer
{
protected:
   SoPlex* thesolver;
   double theeps;

   ///
   int selectLeaveX(int start, int incr);

   ///
   SoPlex::Id selectEnterX(int start1, int incr1, int start2, int incr2);

public:
   ///
   SoPlex* solver() const
   {
      return thesolver;
   }
   ///
   double epsilon() const
   {
      return theeps;
   }
   ///
   void setEpsilon(double eps)
   {
      theeps = eps;
   }
   ///
   void load(SoPlex* p_solver)
   {
      thesolver = p_solver;
   }
   ///
   void clear()
   {
      thesolver = 0;
   }
   ///
   void setType(SoPlex::Type)
   {}
   ///
   void setRep(SoPlex::Representation)
   {}
   ///
   int selectLeave();
   ///
   void left4(int, SoPlex::Id)
   {}

   ///
   SoPlex::Id selectEnter();
   ///
   void entered4(SoPlex::Id, int)
   {}

   ///
   void addedVecs (int)
   {}
   ///
   void addedCoVecs(int)
   {}


   ///
   void removedVec(int)
   {}
   ///
   void removedVecs(const int[])
   {}
   ///
   void removedCoVec(int)
   {}
   ///
   void removedCoVecs(const int[])
   {}


   ///
   void changeObj(const Vector&)
   {}
   ///
   void changeObj(int, double)
   {}
   ///
   void changeLower(const Vector&)
   {}
   ///
   void changeLower(int, double)
   {}
   ///
   void changeUpper(const Vector&)
   {}
   ///
   void changeUpper(int, double)
   {}
   ///
   void changeLhs(const Vector&)
   {}
   ///
   void changeLhs(int, double)
   {}
   ///
   void changeRhs(const Vector&)
   {}
   ///
   void changeRhs(int, double)
   {}
   ///
   void changeRow(int, const LPRow&)
   {}
   ///
   void changeCol(int, const LPCol&)
   {}
   ///
   void changeElement(int, int, double)
   {}
   ///
   void changeSense(SoPlex::Sense)
   {}
};


} // namespace soplex
#endif // _SPXDEFAULTPRR_H_

//-----------------------------------------------------------------------------
//Emacs Local Variables:
//Emacs mode:c++
//Emacs c-basic-offset:3
//Emacs tab-width:8
//Emacs indent-tabs-mode:nil
//Emacs End:
//-----------------------------------------------------------------------------
