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
#pragma ident "@(#) $Id: spxsteeppr.h,v 1.4 2001/11/28 16:41:23 bzfpfend Exp $"


/**@file  spxsteeppr.h
 * @brief Steepest edge pricer.
 */
#ifndef _SPXSTEEPPR_H_
#define _SPXSTEEPPR_H_


#include <assert.h>

#include "spxpricer.h"
#include "random.h"

namespace soplex
{

/**@brief   Steepest edge pricer.
   @ingroup Algo
      
   Class #SPxSteepPR implements a steepest edge pricer to be used with
   #SoPlex.
   
   See #SPxPricer for a class documentation.
*/
class SPxSteepPR : public SPxPricer
{
public:
   /// How to setup the direction multipliers.
   /** Possible settings are #EXACT for starting with exactly computed
       values, or #DEFAULT for starting with multipliers set to 1. The
       latter is the default.
   */
   enum Setup {
      EXACT,   ///< starting with exactly computed values
      DEFAULT  ///< starting with multipliers set to 1
   };

protected:
   DVector penalty;                // vector of pricing penalties
   DVector coPenalty;              // vector of pricing penalties

   DVector workVec;                // working vector
   SSVector workRhs;               // working vector

   int lastIdx;
   SoPlex::Id lastId;
   double pi_p;
   double theeps;

   SoPlex* thesolver;

   int prefSetup;
   DataArray < double > coPref; // preference multiplier for selecting as pivot
   DataArray < double > pref;   // preference multiplier for selecting as pivot
   DataArray < double > leavePref;

   void setupPrefs(double mult,
                   double tie, double cotie,
                   double shift, double coshift,
                   int rstart = 0, int cstart = 0,
                   int rend = -1, int cend = -1);

   virtual void setupPrefs(SoPlex::Type);

public:
   /**@todo make setup and accuracy private or protected */
   /// setup type.
   Setup setup;

   /// accuracy for computing steepest directions.
   double accuracy;

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
   void load(SoPlex* base);
   ///
   void clear();
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
   void left4(int n, SoPlex::Id id, int start, int incr);

   ///
   SoPlex::Id selectEnter();
protected:
   SoPlex::Id selectEnter(double& best, int start1 = 0, int incr1 = 1,
                          int start2 = 0, int incr2 = 1);
   SoPlex::Id otherSelectEnter(double& best, int start1 = 0, int incr1 = 1,
                               int start2 = 0, int incr2 = 1);
public:

   ///
   void entered4(SoPlex::Id id, int n);
   void entered4(SoPlex::Id id, int n, int start1, int incr1, int start2, 
                 int incr2);

   ///
   virtual void addedVecs (int n);
   ///
   virtual void addedCoVecs(int n);

   ///
   virtual void removedVec(int i);
   ///
   virtual void removedCoVecs(const int perm[]);
   ///
   virtual void removedVecs(const int perm[]);
   ///
   virtual void removedCoVec(int i);

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

   ///
   SPxSteepPR()
      : workRhs (0, 1e-16)
      , setup (DEFAULT)
      , accuracy(1e-4)
   {}

};

} // namespace soplex
#endif // _SPXSTEEPPR_H_

//-----------------------------------------------------------------------------
//Emacs Local Variables:
//Emacs mode:c++
//Emacs c-basic-offset:3
//Emacs tab-width:8
//Emacs indent-tabs-mode:nil
//Emacs End:
//-----------------------------------------------------------------------------
