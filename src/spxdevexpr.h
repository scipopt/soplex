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
#pragma ident "@(#) $Id: spxdevexpr.h,v 1.2 2001/11/06 23:31:03 bzfkocht Exp $"


#ifndef _SPXDEVEXPR_H_
#define _SPXDEVEXPR_H_

//@ ----------------------------------------------------------------------------
/*  \Section{Imports}
    Import required system include files ...
 */
#include <assert.h>


/*  ... and class header files
 */

#include "spxpricer.h"

namespace soplex
{






//@ ----------------------------------------------------------------------------
/* \Section{Class Declaration}
 */

/// Devex Pricer for SoPlex.
/** The Devex Pricer for SoPlex implements an approximate steepest edge pricing,
    that does without solving an extra linear system and computing the scalar
    products.
 */
class SPxDevexPR : public SPxPricer
{
private:
protected:
   double last;           // penalty, selected at last iteration
   DVector penalty;        // vector of pricing penalties
   DVector coPenalty;      // vector of pricing penalties

   SoPlex* thesolver;
   double theeps;

public:
   /// return loaded solver.
   SoPlex* solver() const
   {
      return thesolver;
   }
   /// bound violations up to #epsilon# are tollerated.
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
   void load(SoPlex* base);
   ///
   void clear()
   {
      thesolver = 0;
   }
   ///
   void setType(SoPlex::Type);
   ///
   void setRep(SoPlex::Representation rep);

   ///
   int selectLeave();
protected:
   int selectLeave(double& best, int start = 0, int incr = 1);
public:

   ///
   void left4(int n, SoPlex::Id id);
protected:
   void left4(int n, SoPlex::Id id, int start, int incr);
public:

   ///
   SoPlex::Id selectEnter();
protected:
   SoPlex::Id selectEnter(double& best, int start1 = 0, int incr1 = 1,
                          int start2 = 0, int incr2 = 1);
public:

   ///
   void entered4(SoPlex::Id id, int n);
protected:
   void entered4(SoPlex::Id id, int n, int start1, int incr1, int start2, int incr2);
public:


   /// #n# vectors have been added to loaded LP.
   virtual void addedVecs (int n);
   /// #n# covectors have been added to loaded LP.
   virtual void addedCoVecs(int n);


   /// These  methods are use for implemementing the public remove methods.
   virtual void removedVec(int i);
   ///
   virtual void removedCoVecs(const int perm[]);
   ///
   virtual void removedVecs(const int perm[]);
   ///
   virtual void removedCoVec(int i);


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
#endif // _SPXDEVEXPR_H_

//-----------------------------------------------------------------------------
//Emacs Local Variables:
//Emacs mode:c++
//Emacs c-basic-offset:3
//Emacs tab-width:8
//Emacs indent-tabs-mode:nil
//Emacs End:
//-----------------------------------------------------------------------------
