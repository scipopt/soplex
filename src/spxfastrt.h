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
#pragma ident "@(#) $Id: spxfastrt.h,v 1.6 2001/12/28 14:55:13 bzfkocht Exp $"

/**@file  spxfastrt.h
 * @brief Fast shifting ratio test.
 */
#ifndef _SPXFASTRT_H_
#define _SPXFASTRT_H_


#include <assert.h>

#include "spxratiotester.h"

namespace soplex
{

/**@brief   Fast shifting ratio test.
   @ingroup Algo
   
   Class #SPxFastRT is an implementation class of #SPxRatioTester providing
   fast and stable ratio test. Stability is achieved by allowing some
   infeasibility to ensure numerical stability such as the Harris procedure.
   Performance is acchieved by skipping the second phase is the first phase
   allready shows a stable enough pivot.

   See #SPxRatioTester for a class documentation.
*/
class SPxFastRT : public SPxRatioTester
{
protected:
   SoPlex* thesolver;

   /// minimum stability parameter for stopping after phase 1.
   double minStab;
   /// |value| < epsilon is considered 0.
   double epsilon;
   /// currently allowed infeasibility.
   double delta;
   /// initially allowed infeasibility.
   double delta0;

   /// resets tolerances.
   void resetTols();
   /// relaxes stability requirements.
   void relax();
   /// tightens stability requirements.
   void tighten();

   /// Max phase 1 value.
   /** Computes the maximum value \p val that could be used for updating \p upd
       such that it would still fullfill the upper and lower bounds \p up and
       \p low, respectively, within #delta. Return value is the index where the
       minimum value is encounterd. At the same time the maximum absolute value
       of \p upd.delta() is computed and returned in \p abs. Internally all
       loops are started at \p start and incremented by \p incr.
    */
   int maxDelta(double& val,
                double& p_abs,
                UpdateVector& upd,
                Vector& low,
                Vector& up,
                int start,
                int incr);

   ///
   virtual int maxDelta(double& val,
                        double& p_abs);

   ///
   virtual SoPlex::Id maxDelta(int& nr,
                               double& val,
                               double& p_abs);

   /// Min phase 1 value.
   /** Computes the minimum value \p val that could be used for updating \p upd
       such that it would still fullfill the upper and lower bounds \p up and
       \p low, respectively, within #delta. Return value is the index where the
       minimum value is encounterd. At the same time the maximum absolute value
       of \p upd.delta() is computed and returned in \p abs. Internally all
       loops are started at \p start and incremented by \p incr.
   */
   int minDelta(double& val,
                double& p_abs,
                UpdateVector& upd,
                Vector& low,
                Vector& up,
                int start,
                int incr);

   ///
   virtual int minDelta(double& val,
                         double& p_abs,
                         UpdateVector& upd,
                         Vector& low,
                         Vector& up)
   {
      return minDelta(val, p_abs, upd, low, up, 0, 1);
   }

   ///
   virtual int minDelta(double& val,
                        double& p_abs);

   ///
   virtual SoPlex::Id minDelta(int& nr,
                               double& val,
                               double& p_abs);
   
   /// selects stable index for maximizing ratio test.
   /** Selects form all update values \p val < \p max the one with the largest
       value of \p upd.delta() which must be greater than \p stab and is
       returned in \p stab. The index is returned as well as the corresponding
       update value \p val. Internally all loops are started at \p start and
       incremented by \p incr.
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
   
   virtual SoPlex::Id maxSelect(int& nr,
                                double& val,
                                double& stab,
                                double& bestDelta,
                                double max);

   /// selects stable index for minimizing ratio test.
   /** Select form all update values \p val > \p max the one with the largest
       value of \p upd.delta() which must be greater than \p stab and is
       returned in \p stab. The index is returned as well as the corresponding
       update value \p val. Internally all loops are started at \p start and
       incremented by \p incr.
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

   virtual SoPlex::Id minSelect(int& nr,
                                double& val,
                                double& stab,
                                double& bestDelta,
                                double max);


   ///
   int minReleave(double& sel, int leave, double maxabs);
   /// numerical stability tests.
   /** Tests whether the selected leave index needs to be discarded (and do so)
       and the ratio test is to be recomputed.
   */
   int maxReleave(double& sel, int leave, double maxabs);

   ///
   int minShortLeave(double& sel, int leave, double /*max*/, double p_abs);
   /// tests for stop after phase 1.
   /** Tests whether a shortcut after phase 1 is feasible for the 
       selected leave
       pivot. In this case return the update value in \p sel.
   */
   int maxShortLeave(double& sel, int leave, double /*max*/, double p_abs);

   ///
   virtual int minReenter(double& sel, double /*max*/, double maxabs,
                          SoPlex::Id id, int nr);
   /// numerical stability check.
   /** Tests whether the selected enter \p id needs to be discarded (and do so)
       and the ratio test is to be recomputed.
   */
   virtual int maxReenter(double& sel, double /*max*/, double maxabs,
                          SoPlex::Id id, int nr);

   /**@todo the documentation seems to be incorrect. 
            No parameter \p sel exists.
    */
   /// tests for stop after phase 1.
   /** Tests whether a shortcut after phase 1 is feasible for 
       the selected enter
       pivot. In this case return the update value in \p sel.
   */
   virtual int shortEnter(SoPlex::Id& enterId,
                          int nr,
                          double max,
                          double maxabs
                          );

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
//Emacs mode:c++
//Emacs c-basic-offset:3
//Emacs tab-width:8
//Emacs indent-tabs-mode:nil
//Emacs End:
//-----------------------------------------------------------------------------
