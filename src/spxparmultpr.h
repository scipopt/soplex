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
#pragma ident "@(#) $Id: spxparmultpr.h,v 1.3 2001/11/07 17:31:23 bzfbleya Exp $"


#ifndef _SPXPARMULTPR_H_
#define _SPXPARMULTPR_H_

//@ ----------------------------------------------------------------------------
/*      \Section{Imports}
    Import required system include files
 */
#include <assert.h>


/*  and class header files
 */

#include "spxpricer.h"
#include "dataarray.h"
#include "array.h"
#include "ssvector.h"

namespace soplex
{






//@ ----------------------------------------------------------------------------
/* \Section{Class Declaration}
 */
struct SPxParMultPr_Tmp
{
   SoPlex::Id id;
   double test;
};

/** partial multiple pricing.
    Class #SPxParMultPr# is an implementation class for #SPxPricer# implementing
    Dantzig's the default pricing strategy with partial multiple pricing.
    Partial multiple pricing applies to the #ENTER#ing Simplex only. A set of
    #partialSize# eligible pivot indices is selected (partial pricing). In the
    following Simplex iterations pricing is are restricted to these indices
    (multiple pricing) until no more eligable pivots are available. Partial
    multiple pricing significantly reduces the computation time for computing
    the matrix-vector-product in the Simplex algorithm.
 */
class SPxParMultPR : public SPxPricer
{
protected:
   SoPlex* thesolver;
   double theeps;
   DataArray < SPxParMultPr_Tmp > pricSet;
   int multiParts;
   int used;
   int min;
   int last;
   int count;

public:
   /// Set size for partial pricing.
   static int partialSize;

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
   void load(SoPlex* solver);

   ///
   void clear()
   {
      thesolver = 0;
   }

   ///
   void setType(SoPlex::Type tp);
   ///
   void setRep(SoPlex::Representation)
   {}

   ///
   int selectLeave();
   ///
   void left4(int, SoPlex::Id)
   {}

   ///
   SPxLP::Id selectEnter();
   ///
   void entered4(SoPlex::Id id, int n);


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

   ///
   int isConsistent() const;
};


} // namespace soplex
#endif // _SPXPARMULTPRR_H_

//-----------------------------------------------------------------------------
//Emacs Local Variables:
//Emacs mode:c++
//Emacs c-basic-offset:3
//Emacs tab-width:8
//Emacs indent-tabs-mode:nil
//Emacs End:
//-----------------------------------------------------------------------------
