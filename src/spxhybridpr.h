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
#pragma ident "@(#) $Id: spxhybridpr.h,v 1.1 2001/11/06 16:18:33 bzfkocht Exp $"


#ifndef _SPXHYBRIDPR_H_
#define _SPXHYBRIDPR_H_

//@ ----------------------------------------------------------------------------
/*  \Section{Imports}
    Import required system include files ...
 */
#include <assert.h>


/*  ... and class header files
 */

#include "spxpricer.h"
#include "spxdevexpr.h"
#include "spxparmultpr.h"
#include "spxsteeppr.h"

namespace soplex
{






//@ ----------------------------------------------------------------------------
/* \Section{Class Declaration}
 */

/** Hybrid Pricer for SoPlex.
    The hybrid pricer for SoPlex tries to guess the best pricing strategy to
    use for pricing the loaded LP with the loaded algorithm type and basis
    representation. Currently it does so by switching between \Ref{SPxSteepPR},
    \Ref{SPxDevexPR} and \Ref{SPxParMultPR}.
 */
class SPxHybridPR : public SPxPricer
{
   SPxSteepPR steep;
   SPxParMultPR parmult;
   SPxDevexPR devex;

   SPxPricer* thepricer;
   SoPlex* thesolver;
   double theeps;

   double hybridFactor; 

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
   void setEpsilon(double eps);

   ///
   void load(SoPlex* solver);

   ///
   void clear();

   ///
   void setType(SoPlex::Type tp);
   ///
   void setRep(SoPlex::Representation rep);

   ///
   int selectLeave();
   ///
   void left4(int n, SoPlex::Id id);

   ///
   SoPlex::Id selectEnter();
   ///
   void entered4(SoPlex::Id id, int n);

   ///
   void addedVecs (int n);
   ///
   void addedCoVecs(int n);

   ///
   void removedVec(int i);
   ///
   void removedVecs(const int perm[]);
   ///
   void removedCoVec(int i);
   ///
   void removedCoVecs(const int perm[]);


   ///
   void changeObj(const Vector& newObj);
   ///
   void changeObj(int i, double newVal);
   ///
   void changeLower(const Vector& newLower);
   ///
   void changeLower(int i, double newLower);
   ///
   void changeUpper(const Vector& newUpper);
   ///
   void changeUpper(int i, double newUpper);
   ///
   void changeLhs(const Vector& newLhs);
   ///
   void changeLhs(int i, double newLhs);
   ///
   void changeRhs(const Vector& newRhs);
   ///
   void changeRhs(int i, double newRhs);
   ///
   void changeRow(int i, const LPRow& newRow);
   ///
   void changeCol(int i, const LPCol& newCol);
   ///
   void changeElement(int i, int j, double val);
   ///
   void changeSense(SoPlex::Sense sns);
   ///
   int isConsistent() const;

   ///
   SPxHybridPR()
   {
      thesolver = 0;
      // ??? TK20011102 I have no idea what is a reasonable value here.
      hybridFactor = 0.5;
   }
   virtual ~SPxHybridPR()
   {}

};

} // namespace soplex
#endif // _SPXHYBRIDPR_H_

//-----------------------------------------------------------------------------
//Emacs Local Variables:
//Emacs c-basic-offset:3
//Emacs tab-width:8
//Emacs indent-tabs-mode:nil
//Emacs End:
//-----------------------------------------------------------------------------
