/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the class library                   */
/*       SoPlex --- the Sequential object-oriented simPlex.                  */
/*                                                                           */
/*    Copyright (C) 1997-1999 Roland Wunderling                              */
/*                  1997-2002 Konrad-Zuse-Zentrum                            */
/*                            fuer Informationstechnik Berlin                */
/*                                                                           */
/*  SoPlex is distributed under the terms of the ZIB Academic Licence.       */
/*                                                                           */
/*  You should have received a copy of the ZIB Academic License              */
/*  along with SoPlex; see the file COPYING. If not email to soplex@zib.de.  */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
#pragma ident "@(#) $Id: spxsteeppr.h,v 1.18 2005/09/16 12:42:37 bzfhille Exp $"


/**@file  spxsteeppr.h
 * @brief Steepest edge pricer.
 */
#ifndef _SPXSTEEPPR_H_
#define _SPXSTEEPPR_H_


#include <assert.h>

#include "spxdefines.h"
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

private:
   DVector penalty;                // vector of pricing penalties
   DVector coPenalty;              // vector of pricing penalties

   DVector workVec;                // working vector
   SSVector workRhs;               // working vector

   Real pi_p;

   int prefSetup;
   DataArray < Real > coPref; // preference multiplier for selecting as pivot
   DataArray < Real > pref;   // preference multiplier for selecting as pivot
   DataArray < Real > leavePref;

   ///
   void setupPrefsX(Real mult, Real /*tie*/, Real /*cotie*/, Real shift, Real coshift);
   ///
   void setupPrefs(SPxSolver::Type);

public:
   /**@todo make setup and accuracy private or protected */
   /// setup type.
   Setup setup;

   /// accuracy for computing steepest directions.
   Real accuracy;

   ///
   virtual void load(SPxSolver* base);
   ///
   virtual void clear();
   ///
   virtual void setType(SPxSolver::Type);
   ///
   virtual void setRep(SPxSolver::Representation rep);
   ///
   virtual int selectLeave();
   ///
   virtual void left4(int n, SPxId id);
   ///
   virtual SPxId selectEnter();
   ///
   virtual void entered4(SPxId id, int n);
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
#ifndef NO_CONSISTENCY_CHECKS
   ///
   virtual bool isConsistent() const;
#endif

   ///
   SPxSteepPR()
      : SPxPricer("Steep")
      , workRhs (0, 1e-16)
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
