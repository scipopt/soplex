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
#pragma ident "@(#) $Id: spxdefaultpr.h,v 1.2 2001/11/06 23:31:03 bzfkocht Exp $"


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
   void load(SoPlex* solver)
   {
      thesolver = solver;
   }

   ///
   void clear()
   {
      thesolver = 0;
   }

   ///
   void setType(SoPlex::Type tp)
   {
      (void)tp;
   }
   ///
   void setRep(SoPlex::Representation rep)
   {
      (void)rep;
   }

   ///
   int selectLeave(double& bst, int start, int incr);
   ///
   int selectLeave();
   ///
   void left4(int n, SoPlex::Id id)
   {
      (void)n;
      (void)id;
   }

   ///
   SoPlex::Id selectEnter(double& bst, int start1, int incr1, int start2, int incr2);
   ///
   SoPlex::Id selectEnter();
   ///
   void entered4(SoPlex::Id id, int n)
   {
      (void)n;
      (void)id;
   }


   ///
   void addedVecs (int n)
   {
      (void)n;
   }
   ///
   void addedCoVecs(int n)
   {
      (void)n;
   }


   ///
   void removedVec(int i)
   {
      (void)i;
   }
   ///
   void removedVecs(const int perm[])
   {
      (void)perm;
   }
   ///
   void removedCoVec(int i)
   {
      (void)i;
   }
   ///
   void removedCoVecs(const int perm[])
   {
      (void)perm;
   }


   ///
   void changeObj(const Vector& newObj)
   {
      (void)newObj;
   }
   ///
   void changeObj(int i, double newVal)
   {
      (void)newVal;
      (void)i;
   }
   ///
   void changeLower(const Vector& newLower)
   {
      (void)newLower;
   }
   ///
   void changeLower(int i, double newLower)
   {
      (void)i;
      (void)newLower;
   }
   ///
   void changeUpper(const Vector& newUpper)
   {
      (void)newUpper;
   }
   ///
   void changeUpper(int i, double newUpper)
   {
      (void)i;
      (void)newUpper;
   }
   ///
   void changeLhs(const Vector& newLhs)
   {
      (void)newLhs;
   }
   ///
   void changeLhs(int i, double newLhs)
   {
      (void)i;
      (void)newLhs;
   }
   ///
   void changeRhs(const Vector& newRhs)
   {
      (void)newRhs;
   }
   ///
   void changeRhs(int i, double newRhs)
   {
      (void)i;
      (void)newRhs;
   }
   ///
   void changeRow(int i, const LPRow& newRow)
   {
      (void)i;
      (void)newRow;
   }
   ///
   void changeCol(int i, const LPCol& newCol)
   {
      (void)i;
      (void)newCol;
   }
   ///
   void changeElement(int i, int j, double val)
   {
      (void)i;
      (void)j;
      (void)val;
   }
   ///
   void changeSense(SoPlex::Sense sns)
   {
      (void)sns;
   }
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
