/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the class library                   */
/*       SoPlex --- the Sequential object-oriented simPlex.                  */
/*                                                                           */
/*    Copyright (C) 1996-2019 Konrad-Zuse-Zentrum                            */
/*                            fuer Informationstechnik Berlin                */
/*                                                                           */
/*  SoPlex is distributed under the terms of the ZIB Academic Licence.       */
/*                                                                           */
/*  You should have received a copy of the ZIB Academic License              */
/*  along with SoPlex; see the file COPYING. If not email to soplex@zib.de.  */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

#include <assert.h>
#include <iostream>

#include "soplex/spxdefines.h"
#include "soplex/spxparmultpr.h"

namespace soplex
{

template <>
void SPxParMultPR<Real>::setType(typename SPxSolverBase<Real>::Type tp)
{
   if(tp == SPxSolverBase<Real>::ENTER)
   {
      used = 0;
      this->thesolver->setPricing(SPxSolverBase<Real>::PARTIAL);
   }
   else
   {
      this->thesolver->setPricing(SPxSolverBase<Real>::FULL);
   }

   this->thesolver->weights.reDim(0);
   this->thesolver->coWeights.reDim(0);
   this->thesolver->weightsAreSetup = false;

   last = 0;
   min = partialSize / 2;
}

template <>
void SPxParMultPR<Real>::load(SPxSolverBase<Real>* p_solver)
{
   this->thesolver = p_solver;
   multiParts = (this->thesolver->dim() + this->thesolver->coDim()) / partialSize + 1;
   pricSet.reSize(10 * partialSize);
}

template <>
SPxId SPxParMultPR<Real>::selectEnter()
{
   SPxId id;
   Real x;
   int i;
   int best = -1;
   //    const SPxBasisBase<Real>::Desc& ds   = this->thesolver->basis().desc();

   assert(thesolver != 0);
   int lastlast = -1;

   if(this->thesolver->pricing() == SPxSolverBase<Real>::PARTIAL)
   {
      Real val;
      Real eps = -this->theeps;
      lastlast = last;

      for(i = used - 1; i >= 0; --i)
      {
         int n = this->thesolver->number(pricSet[i].id);

         if(this->thesolver->isId(pricSet[i].id))
         {
            this->thesolver->computePvec(n);
            pricSet[i].test = val = this->thesolver->computeTest(n);
         }
         else
            pricSet[i].test = val = this->thesolver->coTest()[n];

         if(val >= eps)
            pricSet[i] = pricSet[--used];
      }

      while(pricSet.size() - used < partialSize)
      {
         best = 0;

         for(i = 1; i < used; ++i)
         {
            if(pricSet[i].test > pricSet[best].test)
               best = i;
         }

         pricSet[best] = pricSet[--used];
      }

      do
      {
         last = (last + 1) % multiParts;

         for(i = this->thesolver->coDim() - last - 1;
               i >= 0; i -= multiParts)
         {
            this->thesolver->computePvec(i);
            x = this->thesolver->computeTest(i);

            if(x < eps)
            {
               pricSet[used].id = this->thesolver->id(i);
               pricSet[used].test = x;
               used++;
            }
         }

         for(i = this->thesolver->dim() - last - 1;
               i >= 0; i -= multiParts)
         {
            x = this->thesolver->coTest()[i];

            if(x < eps)
            {
               pricSet[used].id = this->thesolver->coId(i);
               pricSet[used].test = x;
               used++;
            }
         }

         assert(used < pricSet.size());
      }
      while(used < min && last != lastlast);

      if(used > 0)
      {
         min = (used + 1);

         if(min < 1)
            min = 1;

         if(min > partialSize)
            min = partialSize;

         best = 0;

         for(i = 1; i < used; ++i)
         {
            if(pricSet[i].test < pricSet[best].test)
               best = i;
         }

         id = pricSet[best].id;
      }

      return id;
   }

   else
   {
      assert(this->thesolver->pricing() == SPxSolverBase<Real>::FULL);
      Real bestx = -this->theeps;

      for(i = this->thesolver->dim() - 1; i >= 0; --i)
      {
         x = this->thesolver->coTest()[i];

         // x *= EQ_PREF * (1 + (ds.coStatus(i) == SPxBasisBase<Real>::Desc::P_FREE
         //                || ds.coStatus(i) == SPxBasisBase<Real>::Desc::D_FREE));
         if(x < bestx)
         {
            id = this->thesolver->coId(i);
            bestx = this->thesolver->coTest()[i];
         }
      }

      for(i = this->thesolver->coDim() - 1; i >= 0; --i)
      {
         x = this->thesolver->test()[i];

         // x *= EQ_PREF * (1 + (ds.status(i) == SPxBasisBase<Real>::Desc::P_FREE
         //                || ds.status(i) == SPxBasisBase<Real>::Desc::D_FREE));
         if(x < bestx)
         {
            id = this->thesolver->id(i);
            bestx = this->thesolver->test()[i];
         }
      }

      return id;
   }
}

template <>
int SPxParMultPR<Real>::selectLeave()
{
   int i, n;
   Real x;
   Real best = -this->theeps;
   //    const Real* up  = this->thesolver->ubBound();
   //    const Real* low = this->thesolver->lbBound();

   assert(thesolver != 0);
   n = -1;

   for(i = this->thesolver->dim() - 1; i >= 0; --i)
   {
      x = this->thesolver->fTest()[i];

      // x *= EQ_PREF * (1 + (up[i] == low[i]));
      if(x < best)
      {
         n = i;
         best = this->thesolver->fTest()[i];
      }
   }

   return n;
}
} // namespace soplex
