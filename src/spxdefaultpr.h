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
#pragma ident "@(#) $Id: spxdefaultpr.h,v 1.3 2001/11/07 17:31:21 bzfbleya Exp $"


#ifndef _SPXDEFAULTPR_H_
#define _SPXDEFAULTPR_H_

//@ ----------------------------------------------------------------------------
/*      \Section{Imports}
    Import required system include files
 */
#include <assert.h>


/*  and class header files
 */

#include "spxpricer.h"

namespace soplex
{






//@ ----------------------------------------------------------------------------
/* \Section{Class Declaration}
 */

/** default pricer.
    Class #SPxDefaultPR# is an implementation class for #SPxPricer# implementing
    Dantzig's the default pricing strategy, i.e. maximal/minimal reduced cost or
    maximal violated constraint.
 */
class SPxDefaultPR : public SPxPricer
{
protected:
   SoPlex* thesolver;
   double theeps;

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
   int selectLeave(double& bst, int start, int incr);
   ///
   int selectLeave();
   ///
   void left4(int, SoPlex::Id)
   {}

   ///
   SoPlex::Id selectEnter(double& bst, int start1, int incr1, int start2, int incr2);
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
