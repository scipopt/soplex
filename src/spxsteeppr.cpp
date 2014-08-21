/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the class library                   */
/*       SoPlex --- the Sequential object-oriented simPlex.                  */
/*                                                                           */
/*    Copyright (C) 1996-2014 Konrad-Zuse-Zentrum                            */
/*                            fuer Informationstechnik Berlin                */
/*                                                                           */
/*  SoPlex is distributed under the terms of the ZIB Academic Licence.       */
/*                                                                           */
/*  You should have received a copy of the ZIB Academic License              */
/*  along with SoPlex; see the file COPYING. If not email to soplex@zib.de.  */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

//#define DEBUGGING 1

//TODO may be faster to have a greater zero tolerance for sparse pricing vectors
//     to reduce the number of nonzero entries, e.g. for workVec

#include <assert.h>
#include <iostream>

#include "spxdefines.h"
#include "spxsteeppr.h"
#include "random.h"

#define STEEP_REFINETOL 2.0

namespace soplex
{

// #define EQ_PREF 1000

void SPxSteepPR::clear()
{
   thesolver = 0;
   prefSetup = 0;
}

void SPxSteepPR::load(SPxSolver* base)
{
   thesolver = base;

   if (base)
   {
      workVec.clear();
      workVec.reDim(base->dim());
      workRhs.clear();
      workRhs.reDim(base->dim());

      leavePref.reSize(base->dim());
      coPref.reSize (base->dim());
      pref.reSize (base->coDim());
      prefSetup = 0;
   }
}

void SPxSteepPR::setType(SPxSolver::Type type)
{
   workRhs.setEpsilon(thesolver->epsilon());

   pref.reSize (thesolver->coDim());
   coPref.reSize(thesolver->dim());
   setupPrefs(type);
   setupWeights(type);
   workVec.clear();
   workRhs.clear();
   refined = false;
   if( type == SPxSolver::ENTER && thesolver->hyperPricingEnter )
   {
      if( thesolver->sparsePricingEnter )
      {
         bestPrices.clear();
         bestPrices.setMax(thesolver->dim());
         prices.reMax(thesolver->dim());
      }
      if( thesolver->sparsePricingEnterCo )
      {
         bestPricesCo.clear();
         bestPricesCo.setMax(thesolver->coDim());
         pricesCo.reMax(thesolver->coDim());
      }
   }
   if( type == SPxSolver::LEAVE && thesolver->sparsePricingLeave && thesolver->hyperPricingLeave )
   {
      bestPrices.clear();
      bestPrices.setMax(thesolver->dim());
      prices.reMax(thesolver->dim());
   }
}

void SPxSteepPR::setupWeights(SPxSolver::Type type)
{
   int i;
   if (setup == DEFAULT)
   {
      if (type == SPxSolver::ENTER)
      {
         coPenalty.reDim(thesolver->dim());
         for (i = thesolver->dim() - 1; i >= 0; --i)
            coPenalty[i] = 2;
         penalty.reDim(thesolver->coDim());
         for (i = thesolver->coDim() - 1; i >= 0; --i)
            penalty[i] = 1;
      }
      else
      {
         assert(type == SPxSolver::LEAVE);
         coPenalty.reDim(thesolver->dim());
         for (i = thesolver->dim() - 1; i >= 0; --i)
         {
            const SPxId id = thesolver->basis().baseId(i);
            const int n    = thesolver->number(id);
            assert(n >= 0);
            leavePref[i]   = thesolver->isId(id) ? pref[n] : coPref[n];
            coPenalty[i]   = 1.0;
         }
      }
   }
   else
   {
      MSG_INFO1( spxout << " --- initializing steepest edge multipliers" << std::endl; )

      if (type == SPxSolver::ENTER)
      {
         coPenalty.reDim(thesolver->dim());
         for (i = thesolver->dim() - 1; i >= 0; --i)
            coPenalty[i] = 1;
         penalty.reDim(thesolver->coDim());
         for (i = thesolver->coDim() - 1; i >= 0; --i)
            penalty[i] = 1 + thesolver->vector(i).length2();
      }
      else
      {
         assert(type == SPxSolver::LEAVE);
         coPenalty.reDim(thesolver->dim());
         SSVector tmp(thesolver->dim(), thesolver->epsilon());
         for (i = thesolver->dim() - 1; i >= 0; --i)
         {
            const SPxId id = thesolver->basis().baseId(i);
            const int n    = thesolver->number(id);
            assert(n >= 0);
            leavePref[i]   = thesolver->isId(id) ? pref[n] : coPref[n];
            thesolver->basis().coSolve(tmp, thesolver->unitVector(i));
            coPenalty[i] = tmp.length2();
         }
      }
   }
}

void SPxSteepPR::setupPrefsX(
   Real mult, 
   Real /*tie*/, 
   Real /*cotie*/,
   Real shift, 
   Real coshift)
{
   DataArray<Real>* p;
   DataArray<Real>* cp;
   // Real rtie;
   // Real ctie;
   Real rshift;
   Real cshift;
   int  i;

   if (thesolver->rep() == SPxSolver::COLUMN)
   {
      cp = &pref;
      p  = &coPref;
      // ctie = tie;
      // rtie = cotie;
      cshift = shift;
      rshift = coshift;
   }
   else
   {
      p  = &pref;
      cp = &coPref;
      // rtie = tie;
      // ctie = cotie;
      rshift = shift;
      cshift = coshift;
   }

   //      p[i] += rtie * thesolver->rowVector(i).size() / Real(thesolver->nCols());
   //      p[i] += EQ_PREF * (thesolver->rhs(i) == thesolver->lhs(i));
   //      p[i] += EQ_PREF * (thesolver->rhs(i) >=  infinity
   //                     &&  thesolver->lhs(i) <= -infinity);
   for(i = 0; i < thesolver->nRows(); ++i)
      (*p)[i] = rshift;

   //      cp[i] += ctie * thesolver->colVector(i).size() / Real(thesolver->nRows());
   //      cp[i] += EQ_PREF * (thesolver->upper(i) == thesolver->lower(i));
   //      cp[i] += EQ_PREF * (thesolver->upper(i) >=  infinity
   //                      &&  thesolver->lower(i) <= -infinity);
   for(i = 0; i < thesolver->nCols(); ++i)
      (*cp)[i] = cshift;

   for(i = 0; i < coPref.size(); ++i)
      coPref[i] *= 1.0 - mult * i;

   for(i = 0; i < pref.size(); ++i)
      pref[i] *= 1.0 + mult * i;
}

void SPxSteepPR::setupPrefs(SPxSolver::Type tp)
{
   if (tp != prefSetup)
   {
      Real mult = 1e-8 / Real(1 + thesolver->dim() + thesolver->coDim());

      if (tp == SPxSolver::ENTER)
         setupPrefsX(-mult, -1e-5, -1e-5, 1.0, 1.0);
      else
         setupPrefsX(mult, 1e-5, 1e-5, 1.0, 1.0);

      prefSetup = tp;
   }
}

void SPxSteepPR::setRep(SPxSolver::Representation)
{
   if (workVec.dim() != thesolver->dim())
   {
      DVector tmp = penalty;
      penalty = coPenalty;
      coPenalty = tmp;

      workVec.clear();
      workVec.reDim(thesolver->dim());
   }
}

void SPxSteepPR::left4(int n, SPxId id)
{
   assert(thesolver->type() == SPxSolver::LEAVE);

   //  Update preference multiplier in #leavePref#
   if (thesolver->isId(id))
      leavePref[n] = pref[thesolver->number(id)];
   else if (thesolver->isCoId(id))
      leavePref[n] = coPref[thesolver->number(id)];

   if (id.isValid())
   {
      // Real               delta         = 0.1;   // thesolver->epsilon();
      Real        delta         = 0.1 + 1.0 / thesolver->basis().iteration();
      Real*       coPenalty_ptr = coPenalty.get_ptr();
      const Real* workVec_ptr   = workVec.get_const_ptr();
      const Real* rhoVec        = thesolver->fVec().delta().values();
      Real        rhov_1        = 1.0 / rhoVec[n];
      Real        beta_q        = thesolver->coPvec().delta().length2() * rhov_1 * rhov_1;

      //TK: I gave the 0.5 extra, because I am not sure how hard this assert is.
#ifndef NDEBUG
      if (fabs(rhoVec[n]) < theeps * 0.5)
      {
         MSG_ERROR( spxout << "WSTEEP04: rhoVec = "
                           << rhoVec[n] << " with smaller absolute value than 0.5*theeps = " << 0.5*theeps << std::endl; )
      }
#endif  // NDEBUG

      //  Update #coPenalty# vector
      const IdxSet& rhoIdx = thesolver->fVec().idx();
      int           len    = thesolver->fVec().idx().size();

      for(int i = 0; i < len; ++i)
      {
         int  j = rhoIdx.index(i);
         
         coPenalty_ptr[j] += rhoVec[j] * (beta_q * rhoVec[j] - 2.0 * rhov_1 * workVec_ptr[j]);

         if (coPenalty_ptr[j] < delta)
            coPenalty_ptr[j] = delta; // coPenalty_ptr[j] = delta / (1+delta-x);
         else if (coPenalty_ptr[j] >= infinity)
            coPenalty_ptr[j] = 1.0 / theeps;
      }
      coPenalty_ptr[n] = beta_q;
      //@ coPenalty_ptr[n] = 0.999*beta_q;
      //@ coPenalty_ptr[n] = 1.001*beta_q;
   }
}

int SPxSteepPR::buildBestPriceVectorLeave( Real feastol )
{
   int idx;
   int nsorted;
   Real x;
   const Real* fTest = thesolver->fTest().get_const_ptr();
   const Real* cpen = coPenalty.get_const_ptr();
   const Real* prefPtr = leavePref.get_const_ptr();
   IdxElement price;
   prices.clear();
   bestPrices.clear();

   // construct vector of all prices
   for (int i = thesolver->infeasibilities.size() - 1; i >= 0; --i)
   {
      idx = thesolver->infeasibilities.index(i);
      x = fTest[idx];
      if (x < -feastol)
      {
         if( cpen[idx] < feastol )
            x = x * x / feastol * prefPtr[idx];
         else
            x = x * x / cpen[idx] * prefPtr[idx];
         price.val = x;
         price.idx = idx;
         prices.append(price);
      }
   }
   // set up structures for the quicksort implementation
   compare.elements = prices.get_const_ptr();
   // do a partial sort to move the best ones to the front
   // TODO this can be done more efficiently, since we only need the indices
   nsorted = SPxQuicksortPart(prices.get_ptr(), compare, 0, prices.size(), thesolver->getMaxUpdates());
   // copy indices of best values to bestPrices
   for( int i = 0; i < nsorted; ++i )
   {
      bestPrices.addIdx(prices[i].idx);
      thesolver->isInfeasible[prices[i].idx] = VIOLATED_AND_CHECKED;
   }

   if( nsorted > 0 )
      return prices[0].idx;
   else
      return -1;
}


int SPxSteepPR::selectLeave()
{
   assert(isConsistent());

   int retid;

   if (thesolver->hyperPricingLeave && thesolver->sparsePricingLeave)
   {
      if ( bestPrices.size() < 2 || thesolver->basis().lastUpdate() == 0 )
      {
         // call init method to build up price-vector and return index of largest price
         retid = buildBestPriceVectorLeave(theeps);
      }
      else
         retid = selectLeaveHyper(theeps);
   }
   else if (thesolver->sparsePricingLeave)
      retid = selectLeaveSparse(theeps);
   else
      retid = selectLeaveX(theeps);

   if( retid < 0 && !refined )
   {
      refined = true;
      MSG_INFO3( spxout << "WSTEEP03 trying refinement step..\n"; )
      retid = selectLeaveX(theeps/STEEP_REFINETOL);
   }

   if( retid >= 0 )
   {
      assert( thesolver->coPvec().delta().isConsistent() );
      thesolver->basis().coSolve(thesolver->coPvec().delta(),
                                 thesolver->unitVector(retid));
      assert( thesolver->coPvec().delta().isConsistent() );
      workRhs.setup_and_assign(thesolver->coPvec().delta());
      thesolver->setup4solve(&workVec, &workRhs);
   }

   return retid;
}

int SPxSteepPR::selectLeaveX(Real tol)
{
   const Real* coPenalty_ptr = coPenalty.get_const_ptr();
   const Real* fTest         = thesolver->fTest().get_const_ptr();
   //    const Real* low     = thesolver->lbBound();
   //    const Real* up      = thesolver->ubBound();
   const Real* p             = leavePref.get_const_ptr();

   Real best = -infinity;
   Real x;

   int lastIdx = -1;

   for (int i = thesolver->dim() - 1; i >= 0; --i)
   {
      x = fTest[i];

      if (x < -tol)
      {
         /**@todo this was an assert! is an assertion correct?*/
         // assert(coPenalty_ptr[i] >= theeps);
         if( coPenalty_ptr[i] < tol )
         {
#ifdef ENABLE_ADDITIONAL_CHECKS
            MSG_WARNING( spxout << "WSTEEP02 SPxSteepPR::selectLeaveX(): coPenalty too small ("
                                << coPenalty_ptr[i] << "), assuming epsilon (" << tol << ")!" << std::endl; )
#endif
            x = x * x / tol * p[i];
         }
         else
            x = x * x / coPenalty_ptr[i] * p[i];

         if (x > best)
         {
            best = x;
            lastIdx = i;
         }
      }
   }

   return lastIdx;
}

int SPxSteepPR::selectLeaveSparse(Real tol)
{
   const Real* coPenalty_ptr = coPenalty.get_const_ptr();
   const Real* fTest         = thesolver->fTest().get_const_ptr();
   const Real* p             = leavePref.get_const_ptr();
   Real best = -infinity;
   Real x;
   int lastIdx = -1;
   int idx;

   for (int i = thesolver->infeasibilities.size() - 1; i >= 0; --i)
   {
      idx = thesolver->infeasibilities.index(i);
      x = fTest[idx];

      if (x < -tol)
      {
         if( coPenalty_ptr[idx] < tol )
         {
#ifdef ENABLE_ADDITIONAL_CHECKS
            MSG_WARNING( spxout << "WSTEEP02 SPxSteepPR::selectLeaveSparse(): coPenalty too small ("
                                << coPenalty_ptr[idx] << "), assuming epsilon (" << tol << ")!" << std::endl; )
#endif
            x = x * x / tol * p[idx];
         }
         else
            x = x * x / coPenalty_ptr[idx] * p[idx];

         if (x > best)
         {
            best = x;
            lastIdx = idx;
         }
      }
      else
      {
         thesolver->infeasibilities.remove(i);
         assert(thesolver->isInfeasible[idx] == VIOLATED || thesolver->isInfeasible[idx] == VIOLATED_AND_CHECKED);
         thesolver->isInfeasible[idx] = NOT_VIOLATED;
      }
   }

   return lastIdx;
}

int SPxSteepPR::selectLeaveHyper(Real tol)
{
   const Real* coPen = coPenalty.get_const_ptr();
   const Real* fTest = thesolver->fTest().get_const_ptr();
   const Real* prefPtr  = leavePref.get_const_ptr();

   Real leastBest = infinity;
   Real best = -infinity;
   Real x;

   int bestIdx = -1;
   int idx = 0;

   // find the best price from the short candidate list
   for( int i = bestPrices.size() - 1; i >= 0; --i )
   {
      idx = bestPrices.index(i);
      x = fTest[idx];
      if( x < -tol )
      {
         assert(thesolver->isInfeasible[idx] == VIOLATED || thesolver->isInfeasible[idx] == VIOLATED_AND_CHECKED);
         if( coPen[idx] < -tol )
         {
#ifdef ENABLE_ADDITIONAL_CHECKS
            MSG_WARNING( spxout << "WSTEEP02 SPxSteepPR::selectLeaveSparse(): coPenalty too small ("
                                << coPenalty_ptr[idx] << "), assuming epsilon (" << tol << ")!" << std::endl; )
#endif
            x = x * x / tol * prefPtr[idx];
         }
         else
            x = x * x / coPen[idx] * prefPtr[idx];

         if( x > best )
         {
            best = x;
            bestIdx = idx;
         }
         if( x < leastBest )
            leastBest = x;
      }
      else
      {
         bestPrices.remove(i);
         thesolver->isInfeasible[idx] = NOT_VIOLATED;
      }
   }

   // make sure we do not skip potential candidates due to a high leastBest value
   if( leastBest == infinity )
   {
      assert(bestPrices.size() == 0);
      leastBest = 0;
   }

   // scan the updated indices for a better price
   for( int i = thesolver->updateViols.size() - 1; i >= 0; --i )
   {
      idx = thesolver->updateViols.index(i);
      // is this index a candidate for bestPrices?
      if( thesolver->isInfeasible[idx] == VIOLATED )
      {
         x = fTest[idx];
         assert(x < -tol);
         if( coPen[idx] < -tol )
         {
#ifdef ENABLE_ADDITIONAL_CHECKS
            MSG_WARNING( spxout << "WSTEEP02 SPxSteepPR::selectLeaveSparse(): coPenalty too small ("
                                << coPenalty_ptr[idx] << "), assuming epsilon (" << tol << ")!" << std::endl; )
#endif
            x = x * x / tol * prefPtr[idx];
         }
         else
            x = x * x / coPen[idx] * prefPtr[idx];

         if( x > leastBest )
         {
            if( x > best )
            {
               best = x;
               bestIdx = idx;
            }
            thesolver->isInfeasible[idx] = VIOLATED_AND_CHECKED;
            bestPrices.addIdx(idx);
         }
      }
   }

   return bestIdx;
}

/* Entering Simplex
 */
void SPxSteepPR::entered4(SPxId /* id */, int n)
{
   assert(thesolver->type() == SPxSolver::ENTER);

   if (n >= 0 && n < thesolver->dim())
   {
      Real delta = 2 + 1.0 / thesolver->basis().iteration();
      Real* coPenalty_ptr = coPenalty.get_ptr();
      Real* penalty_ptr = penalty.get_ptr();
      const Real* workVec_ptr = workVec.get_const_ptr();
      const Real* pVec = thesolver->pVec().delta().values();
      const IdxSet& pIdx = thesolver->pVec().idx();
      const Real* coPvec = thesolver->coPvec().delta().values();
      const IdxSet& coPidx = thesolver->coPvec().idx();
      Real xi_p = 1 / thesolver->fVec().delta()[n];
      int i, j;
      Real xi_ip, x;

      assert(thesolver->fVec().delta()[n] > thesolver->epsilon()
              || thesolver->fVec().delta()[n] < -thesolver->epsilon());

      for (j = coPidx.size() - 1; j >= 0; --j)
      {
         i = coPidx.index(j);
         xi_ip = xi_p * coPvec[i];
         x = coPenalty_ptr[i] += xi_ip * (xi_ip * pi_p - 2 * workVec_ptr[i]);
         /*
         if(x < 1)
             coPenalty_ptr[i] = 1 / (2-x);
         */
         if (x < delta)
            coPenalty_ptr[i] = delta;
         // coPenalty_ptr[i] = 1;
         else if (x > infinity)
            coPenalty_ptr[i] = 1 / thesolver->epsilon();
      }

      for (j = pIdx.size() - 1; j >= 0; --j)
      {
         i = pIdx.index(j);
         xi_ip = xi_p * pVec[i];
         x = penalty_ptr[i] += xi_ip * (xi_ip * pi_p - 2.0 * (thesolver->vector(i) * workVec));
         /*
         if(x < 1)
             penalty_ptr[i] = 1 / (2-x);
         */
         if (x < delta)
            penalty_ptr[i] = delta;
         // penalty_ptr[i] = 1;
         else if (x > infinity)
            penalty_ptr[i] = 1.0 / thesolver->epsilon();
      }
   }

   /*@
       if(thesolver->isId(id))
           penalty[   thesolver->number(id) ] *= 1.0001;
       else if(thesolver->isCoId(id))
           coPenalty[ thesolver->number(id) ] *= 1.0001;
   */

}


SPxId SPxSteepPR::buildBestPriceVectorEnterDim( Real& best, Real feastol )
{
   const Real* cp            = coPref.get_const_ptr();
   const Real* coTest        = thesolver->coTest().get_const_ptr();
   const Real* coPenalty_ptr = coPenalty.get_const_ptr();
   int idx;
   int nsorted;
   Real x;
   IdxElement price;

   prices.clear();
   bestPrices.clear();

   // construct vector of all prices
   for( int i = thesolver->infeasibilities.size() - 1; i >= 0; --i )
   {
      idx = thesolver->infeasibilities.index(i);
      x = coTest[idx];
      if ( x < -feastol)
      {
         assert(thesolver->isInfeasible[idx] == VIOLATED || thesolver->isInfeasible[idx] == VIOLATED_AND_CHECKED);
         x = x * x / coPenalty_ptr[idx];
         price.val = x * cp[idx];
         price.idx = idx;
         prices.append(price);
      }
      else
      {
         thesolver->infeasibilities.remove(i);
         thesolver->isInfeasible[idx] = NOT_VIOLATED;
      }
   }
   // set up structures for the quicksort implementation
   compare.elements = prices.get_const_ptr();
   // do a partial sort to move the best ones to the front
   // TODO this can be done more efficiently, since we only need the indices
   nsorted = SPxQuicksortPart(prices.get_ptr(), compare, 0, prices.size(), thesolver->getMaxUpdates());
   // copy indices of best values to bestPrices
   for( int i = 0; i < nsorted; ++i )
   {
      bestPrices.addIdx(prices[i].idx);
      thesolver->isInfeasible[prices[i].idx] = VIOLATED_AND_CHECKED;
   }

   if( nsorted > 0 )
   {
      best = prices[0].val;
      return thesolver->coId(prices[0].idx);
   }
   else
      return SPxId();
}


SPxId SPxSteepPR::buildBestPriceVectorEnterCoDim( Real& best, Real feastol )
{
   const Real* p           = pref.get_const_ptr();
   const Real* test        = thesolver->test().get_const_ptr();
   const Real* penalty_ptr = penalty.get_const_ptr();
   int idx;
   int nsorted;
   Real x;
   IdxElement price;

   pricesCo.clear();
   bestPricesCo.clear();

   // construct vector of all prices
   for( int i = thesolver->infeasibilitiesCo.size() - 1; i >= 0; --i )
   {
      idx = thesolver->infeasibilitiesCo.index(i);
      x = test[idx];
      if ( x < -feastol)
      {
         assert(thesolver->isInfeasibleCo[idx] == VIOLATED || thesolver->isInfeasibleCo[idx] == VIOLATED_AND_CHECKED);
         x = x * x / penalty_ptr[idx];
         price.val = x * p[idx];
         price.idx = idx;
         pricesCo.append(price);
      }
      else
      {
         thesolver->infeasibilitiesCo.remove(i);
         thesolver->isInfeasibleCo[idx] = NOT_VIOLATED;
      }
   }
   // set up structures for the quicksort implementation
   compare.elements = pricesCo.get_const_ptr();
   // do a partial sort to move the best ones to the front
   // TODO this can be done more efficiently, since we only need the indices
   nsorted = SPxQuicksortPart(pricesCo.get_ptr(), compare, 0, pricesCo.size(), thesolver->getMaxUpdates());
   // copy indices of best values to bestPrices
   for( int i = 0; i < nsorted; ++i )
   {
      bestPricesCo.addIdx(pricesCo[i].idx);
      thesolver->isInfeasibleCo[pricesCo[i].idx] = VIOLATED_AND_CHECKED;
   }

   if( nsorted > 0 )
   {
      best = pricesCo[0].val;
      return thesolver->id(pricesCo[0].idx);
   }
   else
      return SPxId();
}


SPxId SPxSteepPR::selectEnter()
{
   assert(thesolver != 0);
   SPxId enterId;

   enterId = selectEnterX(theeps);

   if( !enterId.isValid() && !refined )
   {
      refined = true;
      MSG_INFO3( spxout << "WSTEEP05 trying refinement step..\n"; )
      enterId = selectEnterX(theeps/STEEP_REFINETOL);
   }

   assert(isConsistent());

   if (enterId.isValid())
   {
      SSVector& delta = thesolver->fVec().delta();

      thesolver->basis().solve4update(delta, thesolver->vector(enterId));

      workRhs.setup_and_assign(delta);
      pi_p = 1 + delta.length2();

      thesolver->setup4coSolve(&workVec, &workRhs);
   }
   return enterId;
}

SPxId SPxSteepPR::selectEnterX(Real tol)
{
   SPxId enterId;
   SPxId enterCoId;
   Real best;
   Real bestCo;

   best = -infinity;
   bestCo = -infinity;

   if( thesolver->hyperPricingEnter && !refined )
   {
      if( bestPrices.size() < 2 || thesolver->basis().lastUpdate() == 0 )
         enterCoId = (thesolver->sparsePricingEnter) ? buildBestPriceVectorEnterDim(best, tol) : selectEnterDenseDim(best, tol);
      else
         enterCoId = (thesolver->sparsePricingEnter) ? selectEnterHyperDim(best, tol) : selectEnterDenseDim(best, tol);

      if( bestPricesCo.size() < 2 || thesolver->basis().lastUpdate() == 0 )
         enterId = (thesolver->sparsePricingEnterCo) ? buildBestPriceVectorEnterCoDim(bestCo, tol) : selectEnterDenseCoDim(bestCo, tol);
      else
         enterId = (thesolver->sparsePricingEnterCo) ? selectEnterHyperCoDim(bestCo, tol) : selectEnterDenseCoDim(bestCo, tol);
   }
   else
   {
      enterCoId = (thesolver->sparsePricingEnter && !refined) ? selectEnterSparseDim(best, tol) : selectEnterDenseDim(best, tol);
      enterId = (thesolver->sparsePricingEnterCo && !refined) ? selectEnterSparseCoDim(bestCo, tol) : selectEnterDenseCoDim(bestCo, tol);
   }

   // prefer slack indices to reduce nonzeros in basis matrix
   if( enterCoId.isValid() && (best > SPARSITY_TRADEOFF * bestCo || !enterId.isValid()) )
      return enterCoId;
   else
      return enterId;
}


SPxId SPxSteepPR::selectEnterHyperDim(Real& best, Real tol)
{
   const Real* cp            = coPref.get_const_ptr();
   const Real* coTest        = thesolver->coTest().get_const_ptr();
   const Real* coPenalty_ptr = coPenalty.get_const_ptr();

   Real leastBest = infinity;
   Real x;
   int enterIdx = -1;
   int idx;

   // find the best price from short candidate list
   for( int i = bestPrices.size() - 1; i >= 0; --i )
   {
      idx = bestPrices.index(i);
      x = coTest[idx];
      if( x < -tol )
      {
         x = x * x / coPenalty_ptr[idx];
         x = x * cp[idx];
         if( x > best )
         {
            best = x;
            enterIdx = idx;
         }
         if( x < leastBest )
            leastBest = x;
      }
      else
      {
         bestPrices.remove(i);
         thesolver->isInfeasible[idx] = NOT_VIOLATED;
      }
   }

   // make sure we do not skip potential candidates due to a high leastBest value
   if( leastBest == infinity )
   {
      assert(bestPrices.size() == 0);
      leastBest = 0;
   }

   // scan hte updated indeces for a better price
   for( int i = thesolver->updateViols.size() -1; i >= 0; --i )
   {
      idx = thesolver->updateViols.index(i);
      // only look at indeces that were not checked already
      if( thesolver->isInfeasible[idx] == VIOLATED )
      {
         x = coTest[idx];
         if( x < -tol )
         {
            x = x * x / coPenalty_ptr[idx];
            x = x * cp[idx];
            if( x > leastBest )
            {
               if (x > best)
               {
                  best = x;
                  enterIdx = idx;
               }
               // put index into candidate list
               thesolver->isInfeasible[idx] = VIOLATED_AND_CHECKED;
               bestPrices.addIdx(idx);
            }
         }
         else
         {
            thesolver->isInfeasible[idx] = NOT_VIOLATED;
         }
      }
   }

   if( enterIdx >= 0 )
      return thesolver->coId(enterIdx);
   else
      return SPxId();
}


SPxId SPxSteepPR::selectEnterHyperCoDim(Real& best, Real tol)
{
   const Real* p           = pref.get_const_ptr();
   const Real* test        = thesolver->test().get_const_ptr();
   const Real* penalty_ptr = penalty.get_const_ptr();

   Real leastBest = infinity;
   Real x;
   int enterIdx = -1;
   int idx;

   // find the best price from short candidate list
   for( int i = bestPricesCo.size() - 1; i >= 0; --i )
   {
      idx = bestPricesCo.index(i);
      x = test[idx];
      if( x < -tol )
      {
         x = x * x / penalty_ptr[idx];
         x = x * p[idx];
         if( x > best )
         {
            best = x;
            enterIdx = idx;
         }
         if( x < leastBest )
            leastBest = x;
      }
      else
      {
         bestPricesCo.remove(i);
         thesolver->isInfeasibleCo[idx] = NOT_VIOLATED;
      }
   }

   // make sure we do not skip potential candidates due to a high leastBest value
   if( leastBest == infinity )
   {
      assert(bestPricesCo.size() == 0);
      leastBest = 0;
   }

   // scan the updated indeces for a better price
   for( int i = thesolver->updateViolsCo.size() -1; i >= 0; --i )
   {
      idx = thesolver->updateViolsCo.index(i);
      // only look at indeces that were not checked already
      if( thesolver->isInfeasibleCo[idx] == VIOLATED )
      {
         x = test[idx];
         if( x < -tol )
         {
            x = x * x / penalty_ptr[idx];
            x = x * p[idx];
            if( x > leastBest )
            {
               if (x > best)
               {
                  best = x;
                  enterIdx = idx;
               }
               // put index into candidate list
               thesolver->isInfeasibleCo[idx] = VIOLATED_AND_CHECKED;
               bestPricesCo.addIdx(idx);
            }
         }
         else
         {
            thesolver->isInfeasibleCo[idx] = NOT_VIOLATED;
         }
      }
   }

   if( enterIdx >= 0 )
      return thesolver->id(enterIdx);
   else
      return SPxId();
}


SPxId SPxSteepPR::selectEnterSparseDim(Real& best, Real tol)
{
   SPxId enterId;
   const Real* cp            = coPref.get_const_ptr();
   const Real* coTest        = thesolver->coTest().get_const_ptr();
   const Real* coPenalty_ptr = coPenalty.get_const_ptr();

   int idx;
   Real x;
   Real coPen;
   Real coPrefValue;

   for (int i = thesolver->infeasibilities.size() -1; i >= 0; --i)
   {
      idx = thesolver->infeasibilities.index(i);
      x = coTest[idx];

      if (x < -tol)
      {
         coPen = coPenalty_ptr[idx];
         x = x * x / coPen;
         coPrefValue = cp[idx];
         x = x * coPrefValue;
         // x *= 1 + cp[i];
         if (x > best)
         {
            best = x;
            enterId = thesolver->coId(idx);
         }
      }
      else
      {
         thesolver->infeasibilities.remove(i);
         thesolver->isInfeasible[idx] = NOT_VIOLATED;
      }
   }
   return enterId;
}

SPxId SPxSteepPR::selectEnterSparseCoDim(Real& best, Real tol)
{
   SPxId enterId;
   const Real* p             = pref.get_const_ptr();
   const Real* test          = thesolver->test().get_const_ptr();
   const Real* penalty_ptr   = penalty.get_const_ptr();

   int idx;
   Real x;
   Real pen;
   Real prefValue;

   for (int i = thesolver->infeasibilitiesCo.size() -1; i >= 0; --i)
   {
      idx = thesolver->infeasibilitiesCo.index(i);
      x = test[idx];

      if (x < -tol)
      {
         pen = penalty_ptr[idx];
         x = x * x / pen;
         prefValue = p[idx];
         x = x * prefValue;
         // x *= 1 + p[i];
         if (x > best)
         {
            best   = x;
            enterId = thesolver->id(idx);
         }
      }
      else
      {
         thesolver->infeasibilitiesCo.remove(i);
         thesolver->isInfeasibleCo[idx] = NOT_VIOLATED;
      }
   }
   return enterId;
}

SPxId SPxSteepPR::selectEnterDenseDim(Real& best, Real tol)
{
   SPxId enterId;
   const Real* cp            = coPref.get_const_ptr();
   const Real* coTest        = thesolver->coTest().get_const_ptr();
   const Real* coPenalty_ptr = coPenalty.get_const_ptr();

   Real x;

   for (int i = 0, end = thesolver->dim(); i < end; ++i)
   {
      x = coTest[i];
      if (x < -tol)
      {
         x *= x / coPenalty_ptr[i];
         x *= cp[i];
         // x *= 1 + cp[i];
         if (x > best)
         {
            best = x;
            enterId = thesolver->coId(i);
         }
      }
   }
   return enterId;
}

SPxId SPxSteepPR::selectEnterDenseCoDim(Real& best, Real tol)
{
   SPxId enterId;
   const Real* p             = pref.get_const_ptr();
   const Real* test          = thesolver->test().get_const_ptr();
   const Real* penalty_ptr   = penalty.get_const_ptr();

   Real x;

   for(int i = 0, end = thesolver->coDim(); i < end; ++i)
   {
      x = test[i];
      if (x < -tol)
      {
         x *= x / penalty_ptr[i];
         x *= p[i];
         // x *= 1 + p[i];
         if (x > best)
         {
            best   = x;
            enterId = thesolver->id(i);
         }
      }
   }
   return enterId;
}


void SPxSteepPR::addedVecs(int n)
{
   n = penalty.dim();
   pref.reSize (thesolver->coDim());
   penalty.reDim(thesolver->coDim());

   if (thesolver->type() == SPxSolver::ENTER)
   {
      setupPrefs(thesolver->type());
      for (; n < penalty.dim(); ++n)
         penalty[n] = 2;
   }
   prefSetup = 0;
}

void SPxSteepPR::addedCoVecs(int n)
{
   n = coPenalty.dim();

   leavePref.reSize(thesolver->dim());
   coPref.reSize (thesolver->dim());
   setupPrefs(thesolver->type());

   workVec.reDim (thesolver->dim());
   coPenalty.reDim (thesolver->dim());
   for (; n < coPenalty.dim(); ++n)
      coPenalty[n] = 1;
   prefSetup = 0;
}

void SPxSteepPR::removedVec(int i)
{
   assert(thesolver != 0);
   penalty[i] = penalty[penalty.dim()];
   penalty.reDim(thesolver->coDim());
   prefSetup = 0;
}

void SPxSteepPR::removedVecs(const int perm[])
{
   assert(thesolver != 0);
   if (thesolver->type() == SPxSolver::ENTER)
   {
      int i;
      int j = penalty.dim();
      for (i = 0; i < j; ++i)
         if (perm[i] >= 0)
            penalty[perm[i]] = penalty[i];
   }
   penalty.reDim(thesolver->coDim());
   prefSetup = 0;
}

void SPxSteepPR::removedCoVec(int i)
{
   assert(thesolver != 0);
   coPenalty[i] = coPenalty[coPenalty.dim()];
   coPenalty.reDim(thesolver->dim());
   prefSetup = 0;
}

void SPxSteepPR::removedCoVecs(const int perm[])
{
   assert(thesolver != 0);
   int i;
   int j = coPenalty.dim();
   for (i = 0; i < j; ++i)
      if (perm[i] >= 0)
         coPenalty[perm[i]] = coPenalty[i];
   coPenalty.reDim(thesolver->dim());
   prefSetup = 0;
}

bool SPxSteepPR::isConsistent() const
{
#ifdef ENABLE_CONSISTENCY_CHECKS
   if (thesolver != 0 && thesolver->type() == SPxSolver::LEAVE && setup == EXACT)
   {
      int i;
      SSVector tmp(thesolver->dim(), thesolver->epsilon());
      Real x;
      for (i = thesolver->dim() - 1; i >= 0; --i)
      {
         thesolver->basis().coSolve(tmp, thesolver->unitVector(i));
         x = coPenalty[i] - tmp.length2();
         if (x > thesolver->leavetol() || -x > thesolver->leavetol())
         {
            MSG_ERROR( spxout << "ESTEEP03 x[" << i << "] = " << x << std::endl; )
         }
      }
   }

   if (thesolver != 0 && thesolver->type() == SPxSolver::ENTER)
   {
      int i;
      for (i = thesolver->dim() - 1; i >= 0; --i)
         if (coPenalty[i] < thesolver->epsilon())
            return MSGinconsistent("SPxSteepPR");

      for (i = thesolver->coDim() - 1; i >= 0; --i)
         if (penalty[i] < thesolver->epsilon())
            return MSGinconsistent("SPxSteepPR");
   }
#endif

   return true;
}
} // namespace soplex
