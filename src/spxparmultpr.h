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
#pragma ident "@(#) $Id: spxparmultpr.h,v 1.2 2001/11/06 23:31:04 bzfkocht Exp $"


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
   void setRep(SoPlex::Representation rep)
   {
      (void)rep;
   }

   ///
   int selectLeave();
   ///
   void left4(int n, SoPlex::Id id)
   {
      (void)n;
      (void)id;
   }

   ///
   SPxLP::Id selectEnter();
   ///
   void entered4(SoPlex::Id id, int n);


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
