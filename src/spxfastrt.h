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
#pragma ident "@(#) $Id: spxfastrt.h,v 1.1 2001/11/06 16:18:33 bzfkocht Exp $"

#ifndef _SPXFASTRT_H_
#define _SPXFASTRT_H_


//@ ----------------------------------------------------------------------------
/*      \Section{Imports}
    Import required system include files
 */
#include <assert.h>


/*  and class header files
 */

#include "spxratiotester.h"

namespace soplex
{






//@ ----------------------------------------------------------------------------
/* \Section{Class Declaration}
 */

/** fast shifting ratio test.
    Class #SPxFastRT# is an implementation class of #SPxRatioTester# providing
    fast and stable ratio test. Stability is achieved by allowing some
    infeasibility to ensure numerical stability such as the Harris procedure.
    Performance is acchieved by skipping the second phase is the first phase
    allready shows a stable enough pivot.
 */
class SPxFastRT : public SPxRatioTester
{
protected:
   /// minimum stability parameter for stopping after phase 1.
   double minStab;
   /// #|value| < epsilon# is considered 0.
   double epsilon;
   /// currently allowed infeasibility.
   double delta;
   /// initially allowed infeasibility.
   double delta0;

   /// reset tolerances.
   void resetTols();
   /// relax stability requirements.
   void relax();
   /// tighten stability requirements.
   void tighten();

   /** Max phase 1 value.
       Compute the maximum value #val# that could be used for updating #upd#
       such that it would still fullfill the upper and lower bounds #up# and
       #low#, respectively, within #delta#. Return value is the index where the
       minimum value is encounterd. At the same time the maximum absolute value
       of #upd.delta()# is computed and returned in #abs#. Internally all loops
       are started at #start# and incremented by #incr#.
    */
   int maxDelta(double& val,
                 double& abs,
                 UpdateVector& upd,
                 Vector& low,
                 Vector& up,
                 int start,
                 int incr);

   ///
   virtual int maxDelta(double& val,
                         double& abs);

   ///
   virtual SoPlex::Id maxDelta(
      int& nr,
      double& val,
      double& abs);

   /** Min phase 1 value.
       Compute the minimum value #val# that could be used for updating #upd#
       such that it would still fullfill the upper and lower bounds #up# and
       #low#, respectively, within #delta#. Return value is the index where the
       minimum value is encounterd. At the same time the maximum absolute value
       of #upd.delta()# is computed and returned in #abs#. Internally all loops
       are started at #start# and incremented by #incr#.
    */
   int minDelta(double& val,
                 double& abs,
                 UpdateVector& upd,
                 Vector& low,
                 Vector& up,
                 int start,
                 int incr);

   ///
   virtual int minDelta(double& val,
                         double& abs,
                         UpdateVector& upd,
                         Vector& low,
                         Vector& up)
   {
      return minDelta(val, abs, upd, low, up, 0, 1);
   }

   ///
   virtual int minDelta(double& val,
                         double& abs);

   ///
   virtual SoPlex::Id minDelta(
      int& nr,
      double& val,
      double& abs);

   /** Select stable index for maximizing ratio test.
       Select form all update values #val < max# the one with the largest value
       of #upd.delta()# which must be #> stab# and is returned in #stab#. The
       index is returned as well as the corresponding update #val#ue.
       Internally all loops are started at #start# and incremented by #incr#.
    */
   int maxSelect(double& val,
                  double& stab,
                  double& best,
                  double& bestDelta,
                  double max,
                  const UpdateVector& upd,
                  const Vector& low,
                  const Vector& up,
                  int start = 0,
                  int incr = 1);

   virtual int maxSelect(double& val,
                          double& stab,
                          double& bestDelta,
                          double max);

   virtual SoPlex::Id maxSelect(
      int& nr,
      double& val,
      double& stab,
      double& bestDelta,
      double max);

   /** Select stable index for minimizing ratio test.
       Select form all update values #val > max# the one with the largest value
       of #upd.delta()# which must be #> stab# and is returned in #stab#. The
       index is returned as well as the corresponding update #val#ue.
       Internally all loops are started at #start# and incremented by #incr#.
    */
   int minSelect(double& val,
                  double& stab,
                  double& best,
                  double& bestDelta,
                  double max,
                  const UpdateVector& upd,
                  const Vector& low,
                  const Vector& up,
                  int start = 0,
                  int incr = 1);

   virtual int minSelect(double& val,
                          double& stab,
                          double& bestDelta,
                          double max);

   virtual SoPlex::Id minSelect(
      int& nr,
      double& val,
      double& stab,
      double& bestDelta,
      double max);


   ///
   int minReleave(double& sel, int leave, double maxabs);
   /** numerical stability test.
       Test wheater the selected leave index needs to be discarded (and do so)
       and the ratio test is to be recomputed.
    */
   int maxReleave(double& sel, int leave, double maxabs);

   ///
   int minShortLeave(double& sel, int leave, double max, double abs);
   /** test for stop after phase 1.
       Test wheater a shortcut after phase 1 is feasible for the selected leave
       pivot. In this case return the update value in #sel#.
    */
   int maxShortLeave(double& sel, int leave, double max, double abs);

   ///
   virtual int minReenter(double& sel, double max, double maxabs,
                           SoPlex::Id id, int nr);
   /** Numerical stability check.
       Test wheater the selected enter #Id# needs to be discarded (and do so)
       and the ratio test is to be recomputed.
    */
   virtual int maxReenter(double& sel, double max, double maxabs,
                           SoPlex::Id id, int nr);

   /** test for stop after phase 1.
       Test wheater a shortcut after phase 1 is feasible for the selected enter
       pivot. In this case return the update value in #sel#.
    */
   virtual int shortEnter(SoPlex::Id& enterId,
                           int nr,
                           double max,
                           double maxabs
                        );

   SoPlex* thesolver;

public:
   ///
   SoPlex* solver() const
   {
      return thesolver;
   }

   ///
   void load(SoPlex* solver);

   ///
   void clear();

   ///
   int selectLeave(double& val);

   ///
   SoPlex::Id selectEnter(double& val);

   ///
   void setType(SoPlex::Type);
};

} // namespace soplex
#endif // _SPXFASTRT_H_

//-----------------------------------------------------------------------------
//Emacs Local Variables:
//Emacs c-basic-offset:3
//Emacs tab-width:8
//Emacs indent-tabs-mode:nil
//Emacs End:
//-----------------------------------------------------------------------------
