/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the class library                   */
/*       SoPlex --- the Sequential object-oriented simPlex.                  */
/*                                                                           */
/*  Copyright (c) 1996-2025 Zuse Institute Berlin (ZIB)                      */
/*                                                                           */
/*  Licensed under the Apache License, Version 2.0 (the "License");          */
/*  you may not use this file except in compliance with the License.         */
/*  You may obtain a copy of the License at                                  */
/*                                                                           */
/*      http://www.apache.org/licenses/LICENSE-2.0                           */
/*                                                                           */
/*  Unless required by applicable law or agreed to in writing, software      */
/*  distributed under the License is distributed on an "AS IS" BASIS,        */
/*  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied. */
/*  See the License for the specific language governing permissions and      */
/*  limitations under the License.                                           */
/*                                                                           */
/*  You should have received a copy of the Apache-2.0 license                */
/*  along with SoPlex; see the file LICENSE. If not email to soplex@zib.de.  */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

#include <assert.h>
#include <stdio.h>

#include "soplex/spxdefines.h"

/*
  Here comes our implementation of the Harris procedure improved by shifting
  bounds. The basic idea is to allow a slight infeasibility within |delta| to
  allow for more freedom when selecting the leaving variable. This freedom
  may then be used for selecting numerically stable variables yielding great
  improvements.

  The algorithms operates in two phases. In a first phase, the maximum value
  |val| is determined, when infeasibility within |delta| is allowed. In the second
  phase, between all variables with value < |val| the one is selected which
  has the largest update Vector component. However, this may not
  always yield an improvement. In that case, we shift the variable towards
  infeasibility and retry. This avoids cycling in the shifted LP.
*/

namespace soplex
{


#define SOPLEX_MINSTAB          1e-5
#define SOPLEX_LOWSTAB          1e-10
#define SOPLEX_TRIES                   2
#define SOPLEX_SHORTVAL         1e-5
#define SOPLEX_DELTA_SHIFT      1e-5
#define SOPLEX_EPSILON          1e-10

template <class R>
void SPxFastRT<R>::resetTols()
{
   epsilon = this->tolerances()->scaleAccordingToEpsilon(SOPLEX_EPSILON);
}

template <class R>
void SPxFastRT<R>::tighten()
{
   R delta_shift = this->tolerances()->scaleAccordingToEpsilon(SOPLEX_DELTA_SHIFT);

   if(fastDelta >= this->delta + delta_shift)
   {
      fastDelta -= delta_shift;

      if(fastDelta > this->tolerances()->scaleAccordingToEpsilon(1e-4))
         fastDelta -= 2 * delta_shift;
   }

   if(minStab < this->tolerances()->scaleAccordingToEpsilon(SOPLEX_MINSTAB))
   {
      minStab /= 0.90;

      if(minStab < this->tolerances()->floatingPointFeastol())
         minStab /= 0.90;
   }
}

template <class R>
void SPxFastRT<R>::relax()
{
   R delta_shift = this->tolerances()->scaleAccordingToEpsilon(SOPLEX_DELTA_SHIFT);
   minStab *= 0.95;
   fastDelta += 3 * delta_shift;
}

template <class R>
R SPxFastRT<R>::minStability(R maxabs)
{
   if(maxabs < 1000.0)
      return minStab;

   return maxabs * minStab / 1000.0;
}

/* The code below implements phase 1 of the ratio test. It differs from the description in the
 * Ph.D. thesis page 57 as follows: It uses \f$\delta_i = d_i - s_i - \delta\f$ if \f$d_i > s_i\f$.
 *
 * This change leads to the following behavior of the source code. Consider the first case (x >
 * epsilon, u < infinity): If u - vec[i] <= 0, vec[i] violates the upper bound. In the Harris ratio
 * test, we would compute (u - vec[i] + delta)/upd[i]. The code computes instead delta/upd[i].
 */
template <class R>
int SPxFastRT<R>::maxDelta(
   R& val,                                /* on return: maximum step length */
   R& maxabs,                             /* on return: maximum absolute value in upd VectorBase<R> */
   UpdateVector<R>& update,
   const VectorBase<R>& lowBound,
   const VectorBase<R>& upBound,
   int start,
   int incr) const
{
   int i, sel;
   R x, y, max;
   R u, l;
   bool leaving = this->m_type == SPxSolverBase<R>::LEAVE;
   bool enterrowrep = !leaving && this->thesolver->theRep == SPxSolverBase<R>::ROW;

   R mabs = maxabs;

   const R* up = upBound.get_const_ptr();
   const R* low = lowBound.get_const_ptr();
   const R* vec = update.get_const_ptr();
   const R* upd = update.delta().values();
   const int* idx = update.delta().indexMem();

   sel = -1;
   max = val;

   if(update.delta().isSetup())
   {
      const int* last = idx + update.delta().size();

      for(idx += start; idx < last; idx += incr)
      {
         i = *idx;

         /* in the dual algorithm, bound flips cannot happen, hence we only consider nonbasic variables */
         if(leaving && ((iscoid && this->thesolver->isCoBasic(i)) || (!iscoid
                        && this->thesolver->isBasic(i))))
            continue;

         if(enterrowrep && this->thesolver->baseId(i).isSPxColId()
               && this->thesolver->desc().colStatus(this->thesolver->number(SPxColId(this->thesolver->baseId(i))))
               == SPxBasisBase<R>::Desc::P_FIXED)
            continue;

         x = upd[i];

         if(x > epsilonZero())
         {
            // @todo check wether mabs should be computed only over bounded vars, i.e., in the if block below
            mabs = (x > mabs) ? x : mabs;
            u = up[i];

            if(u < R(infinity))
            {
               y = u - vec[i];

               if(y <= 0)
                  x = fastDelta / x;
               else
                  x = (y + fastDelta) / x;

               if(x < max)
               {
                  max = x;
                  sel = i;
               }
            }
         }
         else if(x < -epsilonZero())
         {
            // @todo check wether mabs should be computed only over bounded vars, i.e., in the if block below
            mabs = (-x > mabs) ? -x : mabs;
            l = low[i];

            if(l > R(-infinity))
            {
               y = l - vec[i];

               if(y >= 0)
                  x = - fastDelta / x;
               else
                  x = (y - fastDelta) / x;

               if(x < max)
               {
                  max = x;
                  sel = i;
               }
            }
         }
      }
   }
   else
   {
      /* In this case, the indices of the semi-sparse Vector update.delta() are not set up and are filled below. */
      int* l_idx = update.delta().altIndexMem();
      R* uval = update.delta().altValues();
      const R* uend = uval + update.dim();
      assert(uval == upd);

      for(i = 0; uval < uend; ++uval, ++i)
      {
         x = *uval;

         if(x != 0.0)
         {
            if(x >= -epsilonZero() && x <= epsilonZero())
            {
               *uval = 0.0;
               continue;
            }
            else
               *l_idx++ = i;

            /* in the dual algorithm, bound flips cannot happen, hence we only consider nonbasic variables */
            if(leaving && ((iscoid && this->thesolver->isCoBasic(i)) || (!iscoid
                           && this->thesolver->isBasic(i))))
               continue;

            if(enterrowrep && this->thesolver->baseId(i).isSPxColId()
                  && this->thesolver->desc().colStatus(this->thesolver->number(SPxColId(this->thesolver->baseId(i))))
                  == SPxBasisBase<R>::Desc::P_FIXED)
               continue;

            if(x > epsilonZero())
            {
               mabs = (x > mabs) ? x : mabs;
               u = up[i];

               if(u < R(infinity))
               {
                  y = u - vec[i];

                  if(y <= 0)
                     x = fastDelta / x;
                  else
                     x = (y + fastDelta) / x;

                  if(x < max)
                  {
                     max = x;
                     sel = i;
                  }
               }
            }
            else if(x < -epsilonZero())
            {
               mabs = (-x > mabs) ? -x : mabs;
               l = low[i];

               if(l > R(-infinity))
               {
                  y = l - vec[i];

                  if(y >= 0)
                     x = - fastDelta / x;
                  else
                     x = (y - fastDelta) / x;

                  if(x < max)
                  {
                     max = x;
                     sel = i;
                  }
               }
            }
         }
      }

      update.delta().setSize(int(l_idx - update.delta().indexMem()));
      update.delta().forceSetup();
   }

   val = max;
   maxabs = mabs;
   return sel;
}

/* See maxDelta() */
template <class R>
int SPxFastRT<R>::minDelta(
   R& val,
   R& maxabs,
   UpdateVector<R>& update,
   const VectorBase<R>& lowBound,
   const VectorBase<R>& upBound,
   int start,
   int incr) const
{
   int i, sel;
   R x, y, max;
   R u, l;
   bool leaving = (this->m_type == SPxSolverBase<R>::LEAVE);
   bool enterrowrep = !leaving && this->thesolver->theRep == SPxSolverBase<R>::ROW;

   R mabs = maxabs;

   const R* up = upBound.get_const_ptr();
   const R* low = lowBound.get_const_ptr();
   const R* vec = update.get_const_ptr();
   const R* upd = update.delta().values();
   const int* idx = update.delta().indexMem();

   sel = -1;
   max = val;

   if(update.delta().isSetup())
   {
      const int* last = idx + update.delta().size();

      for(idx += start; idx < last; idx += incr)
      {
         i = *idx;
         x = upd[i];

         /* in the dual algorithm, bound flips cannot happen, hence we only consider nonbasic variables */
         if(leaving && ((iscoid && this->thesolver->isCoBasic(i)) || (!iscoid
                        && this->thesolver->isBasic(i))))
            continue;

         if(enterrowrep && this->thesolver->baseId(i).isSPxColId()
               && this->thesolver->desc().colStatus(this->thesolver->number(SPxColId(this->thesolver->baseId(i))))
               == SPxBasisBase<R>::Desc::P_FIXED)
            continue;

         if(x > epsilonZero())
         {
            // @todo check wether mabs should be computed only over bounded vars, i.e., in the if block below
            mabs = (x > mabs) ? x : mabs;
            l = low[i];

            if(l > R(-infinity))
            {
               y = l - vec[i];

               if(y >= 0)
                  x = - fastDelta / x;
               else
                  x = (y - fastDelta) / x;

               if(x > max)
               {
                  max = x;
                  sel = i;
               }
            }
         }
         else if(x < -epsilonZero())
         {
            // @todo check wether mabs should be computed only over bounded vars, i.e., in the if block below
            mabs = (-x > mabs) ? -x : mabs;
            u = up[i];

            if(u < R(infinity))
            {
               y = u - vec[i];

               if(y <= 0)
                  x = fastDelta / x;
               else
                  x = (y + fastDelta) / x;

               if(x > max)
               {
                  max = x;
                  sel = i;
               }
            }
         }
      }
   }
   else
   {
      /* In this case, the indices of the semi-sparse VectorBase<R> update.delta() are not set up and are filled below. */
      int* l_idx = update.delta().altIndexMem();
      R* uval = update.delta().altValues();
      const R* uend = uval + update.dim();
      assert(uval == upd);

      for(i = 0; uval < uend; ++uval, ++i)
      {
         x = *uval;

         if(x != 0.0)
         {
            if(x >= -epsilonZero() && x <= epsilonZero())
            {
               *uval = 0.0;
               continue;
            }
            else
               *l_idx++ = i;

            /* in the dual algorithm, bound flips cannot happen, hence we only consider nonbasic variables */
            if(leaving && ((iscoid && this->thesolver->isCoBasic(i)) || (!iscoid
                           && this->thesolver->isBasic(i))))
               continue;

            if(enterrowrep && this->thesolver->baseId(i).isSPxColId()
                  && this->thesolver->desc().colStatus(this->thesolver->number(SPxColId(this->thesolver->baseId(i))))
                  == SPxBasisBase<R>::Desc::P_FIXED)
               continue;


            if(x > epsilonZero())
            {
               mabs = (x > mabs) ? x : mabs;
               l = low[i];

               if(l > R(-infinity))
               {
                  y = l - vec[i];

                  if(y >= 0)
                     x = - fastDelta / x;
                  else
                     x = (y - fastDelta) / x;

                  if(x > max)
                  {
                     max = x;
                     sel = i;
                  }
               }
            }
            else if(x < -epsilonZero())
            {
               mabs = (-x > mabs) ? -x : mabs;
               u = up[i];

               if(u < R(infinity))
               {
                  y = u - vec[i];

                  if(y <= 0)
                     x = fastDelta / x;
                  else
                     x = (y + fastDelta) / x;

                  if(x > max)
                  {
                     max = x;
                     sel = i;
                  }
               }
            }
         }
      }

      update.delta().setSize(int(l_idx - update.delta().indexMem()));
      update.delta().forceSetup();
   }

   val = max;
   maxabs = mabs;

   return sel;
}

template <class R>
int SPxFastRT<R>::maxDelta(
   R& val,
   R& maxabs)
{
   assert(this->m_type == SPxSolverBase<R>::ENTER);
   return maxDelta(val, maxabs,
                   this->thesolver->fVec(), this->thesolver->lbBound(), this->thesolver->ubBound(), 0, 1);
}

template <class R>
int SPxFastRT<R>::minDelta(
   R& val,
   R& maxabs)
{
   assert(this->m_type == SPxSolverBase<R>::ENTER);
   return minDelta(val, maxabs,
                   this->thesolver->fVec(), this->thesolver->lbBound(), this->thesolver->ubBound(), 0, 1);
}

template <class R>
SPxId SPxFastRT<R>::maxDelta(
   int& nr,
   R& max,                                /* on return: maximum step length */
   R& maxabs)                             /* on return: maximum absolute value in delta VectorBase<R> */
{
   /* The following cause side effects on coPvec and pVec - both changes may be needed later in
      maxSelect(). We can therefore not move the first function after the (indp >= 0) check. */
   iscoid = true;
   int indc = maxDelta(max, maxabs,
                       this->thesolver->coPvec(), this->thesolver->lcBound(), this->thesolver->ucBound(), 0, 1);
   iscoid = false;
   int indp = maxDelta(max, maxabs,
                       this->thesolver->pVec(), this->thesolver->lpBound(), this->thesolver->upBound(), 0, 1);

   if(indp >= 0)
   {
      nr = indp;
      return this->thesolver->id(indp);
   }

   if(indc >= 0)
   {
      nr = indc;
      return this->thesolver->coId(indc);
   }

   nr = -1;
   return SPxId();
}

template <class R>
SPxId SPxFastRT<R>::minDelta(
   int& nr,
   R& max,
   R& maxabs)
{
   /* The following cause side effects on coPvec and pVec - both changes may be needed later in
      minSelect(). We can therefore not move the first function after the (indp >= 0) check. */
   iscoid = true;
   const int indc = minDelta(max, maxabs,
                             this->thesolver->coPvec(), this->thesolver->lcBound(), this->thesolver->ucBound(), 0, 1);
   iscoid = false;
   const int indp = minDelta(max, maxabs,
                             this->thesolver->pVec(), this->thesolver->lpBound(), this->thesolver->upBound(), 0, 1);

   if(indp >= 0)
   {
      nr = indp;
      return this->thesolver->id(indp);
   }

   if(indc >= 0)
   {
      nr = indc;
      return this->thesolver->coId(indc);
   }

   nr = -1;
   return SPxId();
}

/* \p best returns the minimum update value such that the corresponding value of \p upd.delta() is
 * at least \p stab and the update value is smaller than \p max. If no valid update value has been
 * found \p bestDelta returns the slack to the bound corresponding to the index used for \p best. */
template <class R>
int SPxFastRT<R>::minSelect(
   R& val,
   R& stab,
   R& best,
   R& bestDelta,
   R max,
   const UpdateVector<R>& update,
   const VectorBase<R>& lowBound,
   const VectorBase<R>& upBound,
   int start,
   int incr) const
{
   int i;
   R x, y;
   bool leaving = this->m_type == SPxSolverBase<R>::LEAVE;
   bool enterrowrep = !leaving && this->thesolver->theRep == SPxSolverBase<R>::ROW;

   const R* up = upBound.get_const_ptr();
   const R* low = lowBound.get_const_ptr();
   const R* vec = update.get_const_ptr();
   const R* upd = update.delta().values();
   const int* idx = update.delta().indexMem();
   const int* last = idx + update.delta().size();

   int nr = -1;
   int bestNr = -1;

   for(idx += start; idx < last; idx += incr)
   {
      i = *idx;
      x = upd[i];

      // in the dual algorithm, bound flips cannot happen, hence we only consider nonbasic variables
      if(leaving && ((iscoid && this->thesolver->isCoBasic(i)) || (!iscoid
                     && this->thesolver->isBasic(i))))
         continue;

      if(enterrowrep && this->thesolver->baseId(i).isSPxColId()
            && this->thesolver->desc().colStatus(this->thesolver->number(SPxColId(this->thesolver->baseId(i))))
            == SPxBasisBase<R>::Desc::P_FIXED)
         continue;

      if(x > stab)
      {
         y = (low[i] - vec[i]) / x;

         if(y >= max)
         {
            val = y;
            nr = i;
            stab = x;
         }
         else if(y < best)
         {
            best = y;
            bestNr = i;
         }
      }
      else if(x < -stab)
      {
         y = (up[i] - vec[i]) / x;

         if(y >= max)
         {
            val = y;
            nr = i;
            stab = -x;
         }
         else if(y < best)
         {
            best = y;
            bestNr = i;
         }
      }
   }

   if(nr < 0 && bestNr > 0)
   {
      if(upd[bestNr] < 0)
         bestDelta = up[bestNr] - vec[bestNr];
      else
         bestDelta = vec[bestNr] - low[bestNr];
   }

   return nr;
}

/* See minSelect() */
template <class R>
int SPxFastRT<R>::maxSelect(
   R& val,
   R& stab,
   R& best,
   R& bestDelta,
   R max,
   const UpdateVector<R>& update,
   const VectorBase<R>& lowBound,
   const VectorBase<R>& upBound,
   int start,
   int incr) const
{
   int i;
   R x, y;
   bool leaving = this->m_type == SPxSolverBase<R>::LEAVE;
   bool enterrowrep = !leaving && this->thesolver->theRep == SPxSolverBase<R>::ROW;

   const R* up = upBound.get_const_ptr();
   const R* low = lowBound.get_const_ptr();
   const R* vec = update.get_const_ptr();
   const R* upd = update.delta().values();
   const int* idx = update.delta().indexMem();
   const int* last = idx + update.delta().size();

   int nr = -1;
   int bestNr = -1;

   for(idx += start; idx < last; idx += incr)
   {
      i = *idx;
      x = upd[i];

      // in the dual algorithm, bound flips cannot happen, hence we only consider nonbasic variables
      if(leaving && ((iscoid && this->thesolver->isCoBasic(i)) || (!iscoid
                     && this->thesolver->isBasic(i))))
         continue;

      if(enterrowrep && this->thesolver->baseId(i).isSPxColId()
            && this->thesolver->desc().colStatus(this->thesolver->number(SPxColId(this->thesolver->baseId(i))))
            == SPxBasisBase<R>::Desc::P_FIXED)
         continue;

      if(x > stab)
      {
         y = (up[i] - vec[i]) / x;

         if(y <= max)
         {
            val = y;
            nr = i;
            stab = x;
         }
         else if(y > best)
         {
            best = y;
            bestNr = i;
         }
      }
      else if(x < -stab)
      {
         y = (low[i] - vec[i]) / x;

         if(y <= max)
         {
            val = y;
            nr = i;
            stab = -x;
         }
         else if(y > best)
         {
            best = y;
            bestNr = i;
         }
      }
   }

   if(nr < 0 && bestNr > 0)
   {
      if(upd[bestNr] > 0)
         bestDelta = up[bestNr] - vec[bestNr];
      else
         bestDelta = vec[bestNr] - low[bestNr];
   }

   return nr;
}

template <class R>
int SPxFastRT<R>::maxSelect(
   R& val,
   R& stab,
   R& bestDelta,
   R max)
{
   R best = R(-infinity);
   bestDelta = 0.0;
   assert(this->m_type == SPxSolverBase<R>::ENTER);
   return maxSelect(val, stab, best, bestDelta, max,
                    this->thesolver->fVec(), this->thesolver->lbBound(), this->thesolver->ubBound(),  0, 1);
}

template <class R>
SPxId SPxFastRT<R>::maxSelect(
   int& nr,
   R& val,
   R& stab,
   R& bestDelta,
   R max
)
{
   int indp, indc;
   R best = R(-infinity);
   bestDelta = 0.0;
   iscoid = true;
   indc = maxSelect(val, stab, best, bestDelta, max,
                    this->thesolver->coPvec(), this->thesolver->lcBound(), this->thesolver->ucBound(), 0, 1);
   iscoid = false;
   indp = maxSelect(val, stab, best, bestDelta, max,
                    this->thesolver->pVec(), this->thesolver->lpBound(), this->thesolver->upBound(), 0, 1);

   if(indp >= 0)
   {
      nr = indp;
      return this->thesolver->id(indp);
   }

   if(indc >= 0)
   {
      nr = indc;
      return this->thesolver->coId(indc);
   }

   nr = -1;
   return SPxId();
}

template <class R>
int SPxFastRT<R>::minSelect(
   R& val,
   R& stab,
   R& bestDelta,
   R max)
{
   R best = R(infinity);
   bestDelta = 0.0;
   assert(this->m_type == SPxSolverBase<R>::ENTER);
   return minSelect(val, stab, best, bestDelta, max,
                    this->thesolver->fVec(), this->thesolver->lbBound(), this->thesolver->ubBound(), 0, 1);
}

template <class R>
SPxId SPxFastRT<R>::minSelect(
   int& nr,
   R& val,
   R& stab,
   R& bestDelta,
   R max)
{
   R best = R(infinity);
   bestDelta = 0.0;
   iscoid = true;
   int indc = minSelect(val, stab, best, bestDelta, max,
                        this->thesolver->coPvec(), this->thesolver->lcBound(), this->thesolver->ucBound(), 0, 1);
   iscoid = false;
   int indp = minSelect(val, stab, best, bestDelta, max,
                        this->thesolver->pVec(), this->thesolver->lpBound(), this->thesolver->upBound(), 0, 1);

   if(indp >= 0)
   {
      nr = indp;
      return this->thesolver->id(indp);
   }

   if(indc >= 0)
   {
      nr = indc;
      return this->thesolver->coId(indc);
   }

   nr = -1;
   return SPxId();
}

template <class R>
bool SPxFastRT<R>::maxShortLeave(R& sel, int leave, R maxabs)
{
   assert(leave >= 0);
   assert(maxabs >= 0);

   R shortval = this->tolerances()->scaleAccordingToEpsilon(SOPLEX_SHORTVAL);

   sel = this->thesolver->fVec().delta()[leave];

   if(sel > maxabs * shortval)
   {
      sel = (this->thesolver->ubBound()[leave] - this->thesolver->fVec()[leave]) / sel;
      return true;
   }

   if(sel < -maxabs * shortval)
   {
      sel = (this->thesolver->lbBound()[leave] - this->thesolver->fVec()[leave]) / sel;
      return true;
   }

   return false;
}

template <class R>
bool SPxFastRT<R>::minShortLeave(R& sel, int leave, R maxabs)
{
   assert(leave >= 0);
   assert(maxabs >= 0);

   R shortval = this->tolerances()->scaleAccordingToEpsilon(SOPLEX_SHORTVAL);

   sel = this->thesolver->fVec().delta()[leave];

   if(sel > maxabs * shortval)
   {
      sel = (this->thesolver->lbBound()[leave] - this->thesolver->fVec()[leave]) / sel;
      return true;
   }

   if(sel < -maxabs * shortval)
   {
      sel = (this->thesolver->ubBound()[leave] - this->thesolver->fVec()[leave]) / sel;
      return true;
   }

   return false;
}

template <class R>
bool SPxFastRT<R>::maxReLeave(R& sel, int leave, R maxabs, bool polish)
{
   UpdateVector<R>& vec = this->thesolver->fVec();
   VectorBase<R>& low = this->thesolver->lbBound();
   VectorBase<R>& up = this->thesolver->ubBound();

   if(leave < 0)
      return true;

   if(up[leave] > low[leave])
   {
      R x = vec.delta()[leave];

      if(sel < -fastDelta / maxabs)
      {
         sel = 0.0;

         // prevent shifts in polishing mode to avoid a final cleanup step (i.e. simplex type switch)
         if(!polish
               && this->thesolver->dualStatus(this->thesolver->baseId(leave)) != SPxBasisBase<R>::Desc::D_ON_BOTH)
         {
            if(x < 0.0)
               this->thesolver->shiftLBbound(leave, vec[leave]);
            else
               this->thesolver->shiftUBbound(leave, vec[leave]);
         }
      }
   }
   else
   {
      sel = 0.0;

      // prevent shifts in polishing mode to avoid a final cleanup step (i.e. simplex type switch)
      if(!polish)
      {
         this->thesolver->shiftLBbound(leave, vec[leave]);
         this->thesolver->shiftUBbound(leave, vec[leave]);
      }
   }

   return false;
}

template <class R>
bool SPxFastRT<R>::minReLeave(R& sel, int leave, R maxabs, bool polish)
{
   UpdateVector<R>& vec = this->thesolver->fVec();
   VectorBase<R>& low = this->thesolver->lbBound();
   VectorBase<R>& up = this->thesolver->ubBound();

   if(leave < 0)
      return true;

   if(up[leave] > low[leave])
   {
      R x = vec.delta()[leave];

      if(sel > fastDelta / maxabs)
      {
         sel = 0.0;

         if(!polish
               && this->thesolver->dualStatus(this->thesolver->baseId(leave)) != SPxBasisBase<R>::Desc::D_ON_BOTH)
            // prevent shifts in polishing mode to avoid a final cleanup step (i.e. simplex type switch)
         {
            if(x > 0.0)
               this->thesolver->shiftLBbound(leave, vec[leave]);
            else
               this->thesolver->shiftUBbound(leave, vec[leave]);
         }
      }
   }
   else
   {
      sel = 0.0;

      // prevent shifts in polishing mode to avoid a final cleanup step (i.e. simplex type switch)
      if(!polish)
      {
         this->thesolver->shiftLBbound(leave, vec[leave]);
         this->thesolver->shiftUBbound(leave, vec[leave]);
      }
   }

   return false;
}

template <class R>
int SPxFastRT<R>::selectLeave(R& val, R, bool polish)
{
   R maxabs, max, sel;
   R delta_shift = this->tolerances()->scaleAccordingToEpsilon(SOPLEX_DELTA_SHIFT);
   int leave = -1;
   int cnt = 0;

   assert(this->m_type == SPxSolverBase<R>::ENTER);

   // force instable pivot iff true (see explanation in enter.cpp and spxsolve.hpp)
   bool instable = this->solver()->instableEnter;
   R lowstab = this->tolerances()->scaleAccordingToEpsilon(SOPLEX_LOWSTAB);
   assert(!instable || this->solver()->instableEnterId.isValid());

   resetTols();

   if(val > epsilonZero())
   {
      do
      {
         // phase 1:
         max = val;
         maxabs = 0.0;
         leave = maxDelta(max, maxabs);

         assert(leave < 0 || !(this->thesolver->baseId(leave).isSPxColId()) ||
                this->thesolver->desc().colStatus(this->thesolver->number(SPxColId(this->thesolver->baseId(
                         leave)))) != SPxBasisBase<R>::Desc::P_FIXED);

         if(max == val || leave == -1)
         {
            assert(max == val && leave == -1);
            return -1;
         }

         if(!maxShortLeave(sel, leave, maxabs))
         {
            // phase 2:
            R stab, bestDelta;

            stab = 100.0 * minStability(maxabs);

            // force instable pivot iff instable is true (see explanation in enter.hpp and spxsolve.hpp)
            if(instable)
               leave = maxSelect(sel, lowstab, bestDelta, max);
            else
               leave = maxSelect(sel, stab, bestDelta, max);

            if(bestDelta < delta_shift * SOPLEX_TRIES)
               cnt++;
            else
               cnt += SOPLEX_TRIES;
         }

         if(!maxReLeave(sel, leave, maxabs, polish))
            break;

         relax();
      }
      while(cnt < SOPLEX_TRIES);
   }
   else if(val < -epsilonZero())
   {
      do
      {
         max = val;
         maxabs = 0;
         leave = minDelta(max, maxabs);

         assert(leave < 0 || !(this->thesolver->baseId(leave).isSPxColId()) ||
                this->thesolver->desc().colStatus(this->thesolver->number(SPxColId(this->thesolver->baseId(
                         leave)))) != SPxBasisBase<R>::Desc::P_FIXED);

         if(max == val || leave == -1)
         {
            assert(max == val && leave == -1);
            return -1;
         }

         if(!minShortLeave(sel, leave, maxabs))
         {
            // phase 2:
            R stab, bestDelta;

            stab = 100.0 * minStability(maxabs);

            // force instable pivot iff instable is true (see explanation in enter.hpp and spxsolve.hpp)
            if(instable)
               leave = minSelect(sel, lowstab, bestDelta, max);
            else
               leave = minSelect(sel, stab, bestDelta, max);

            assert(leave < 0 || !(this->thesolver->baseId(leave).isSPxColId())
                   || this->thesolver->desc().colStatus(this->thesolver->number(SPxColId(this->thesolver->baseId(
                            leave)))) != SPxBasisBase<R>::Desc::P_FIXED);

            if(bestDelta < delta_shift * SOPLEX_TRIES)
               cnt++;
            else
               cnt += SOPLEX_TRIES;
         }

         if(!minReLeave(sel, leave, maxabs, polish))
            break;

         relax();
      }
      while(cnt < SOPLEX_TRIES);
   }
   else
      return -1;

   SPX_DEBUG(

      if(leave >= 0)
      std::cout
      << "DFSTRT01 "
      << this->thesolver->basis().iteration() << "("
      << std::setprecision(6) << this->thesolver->value() << ","
      << std::setprecision(2) << this->thesolver->basis().stability() << "):"
      << leave << "\t"
      << std::setprecision(4) << sel << " "
      << std::setprecision(4) << this->thesolver->fVec().delta()[leave] << " "
      << std::setprecision(6) << maxabs
      << std::endl;
      else
         std::cout << "DFSTRT02 " << this->thesolver->basis().iteration()
         << ": skipping instable pivot" << std::endl;
      )

         if(polish && leave >= 0)
         {
            assert(this->thesolver->rep() == SPxSolverBase<R>::COLUMN);
            SPxId leaveId = this->thesolver->baseId(leave);

            // decide whether the chosen leave index contributes to the polishing objective
            if(this->thesolver->polishObj == SPxSolverBase<R>::POLISH_INTEGRALITY)
            {
               // only allow (integer) variables to leave the basis
               if(leaveId.isSPxRowId())
                  return -1;
               else if(this->thesolver->integerVariables.size() == this->thesolver->nCols())
               {
                  if(leaveId.isSPxColId() && this->thesolver->integerVariables[this->thesolver->number(leaveId)] == 0)
                     return -1;
               }
            }
            else if(this->thesolver->polishObj == SPxSolverBase<R>::POLISH_FRACTIONALITY)
            {
               // only allow slacks and continuous variables to leave the basis
               if(this->thesolver->integerVariables.size() == this->thesolver->nCols())
               {
                  if(this->thesolver->baseId(leave).isSPxColId()
                        && this->thesolver->integerVariables[this->thesolver->number(leaveId)] == 1)
                     return -1;
               }
               else if(this->thesolver->baseId(leave).isSPxColId())
                  return -1;
            }
         }

   if(leave >= 0 || minStab > 2 * this->solver()->epsilon())
   {
      val = sel;

      if(leave >= 0)
         tighten();
   }

   assert(leave < 0 || !(this->thesolver->baseId(leave).isSPxColId())
          || this->thesolver->desc().colStatus(this->thesolver->number(SPxColId(this->thesolver->baseId(
                   leave)))) != SPxBasisBase<R>::Desc::P_FIXED);

   return leave;
}


template <class R>
bool SPxFastRT<R>::maxReEnter(R& sel,
                              R maxabs,
                              const SPxId& id,
                              int nr,
                              bool polish)
{
   R x, d;
   VectorBase<R>* up;
   VectorBase<R>* low;

   UpdateVector<R>& pvec = this->thesolver->pVec();
   SSVectorBase<R>& pupd = this->thesolver->pVec().delta();
   VectorBase<R>& upb = this->thesolver->upBound();
   VectorBase<R>& lpb = this->thesolver->lpBound();
   UpdateVector<R>& cvec = this->thesolver->coPvec();
   SSVectorBase<R>& cupd = this->thesolver->coPvec().delta();
   VectorBase<R>& ucb = this->thesolver->ucBound();
   VectorBase<R>& lcb = this->thesolver->lcBound();

   if(this->thesolver->isCoId(id))
   {
      if(this->thesolver->isCoBasic(nr))
      {
         cupd.clearIdx(nr);
         return true;
      }

      x = cvec[nr];
      d = cupd[nr];
      up = &ucb;
      low = &lcb;

      if(d < 0.0)
         sel = (lcb[nr] - cvec[nr]) / d;
      else
         sel = (ucb[nr] - cvec[nr]) / d;
   }
   else if(this->thesolver->isId(id))
   {
      pvec[nr] = this->thesolver->vector(nr) * cvec;

      if(this->thesolver->isBasic(nr))
      {
         pupd.clearIdx(nr);
         return true;
      }

      x = pvec[nr];
      d = pupd[nr];
      up = &upb;
      low = &lpb;

      if(d < 0.0)
         sel = (lpb[nr] - pvec[nr]) / d;
      else
         sel = (upb[nr] - pvec[nr]) / d;
   }
   else
      return true;

   if((*up)[nr] != (*low)[nr])
   {
      if(sel < -fastDelta / maxabs)
      {
         sel = 0.0;

         // prevent shifts in polishing mode to avoid a final cleanup step (i.e. simplex type switch)
         if(!polish)
         {
            if(d > 0.0)
            {
               this->thesolver->theShift -= (*up)[nr];
               (*up)[nr] = x;
               this->thesolver->theShift += (*up)[nr];
            }
            else
            {
               this->thesolver->theShift += (*low)[nr];
               (*low)[nr] = x;
               this->thesolver->theShift -= (*low)[nr];
            }
         }
      }
   }
   else
   {
      sel = 0.0;

      // prevent shifts in polishing mode to avoid a final cleanup step (i.e. simplex type switch)
      if(!polish)
      {
         if(x > (*up)[nr])
            this->thesolver->theShift += x - (*up)[nr];
         else
            this->thesolver->theShift += (*low)[nr] - x;

         (*up)[nr] = (*low)[nr] = x;
      }
   }

   return false;
}

template <class R>
bool SPxFastRT<R>::minReEnter(
   R& sel,
   R maxabs,
   const SPxId& id,
   int nr,
   bool polish)
{
   R x, d;
   VectorBase<R>* up;
   VectorBase<R>* low;

   UpdateVector<R>& pvec = this->thesolver->pVec();
   SSVectorBase<R>& pupd = this->thesolver->pVec().delta();
   VectorBase<R>& upb = this->thesolver->upBound();
   VectorBase<R>& lpb = this->thesolver->lpBound();
   UpdateVector<R>& cvec = this->thesolver->coPvec();
   SSVectorBase<R>& cupd = this->thesolver->coPvec().delta();
   VectorBase<R>& ucb = this->thesolver->ucBound();
   VectorBase<R>& lcb = this->thesolver->lcBound();

   if(this->thesolver->isCoId(id))
   {
      if(this->thesolver->isCoBasic(nr))
      {
         cupd.clearIdx(nr);
         return true;
      }

      x = cvec[nr];
      d = cupd[nr];
      up = &ucb;
      low = &lcb;

      if(d > 0.0)
         sel = (this->thesolver->lcBound()[nr] - cvec[nr]) / d;
      else
         sel = (this->thesolver->ucBound()[nr] - cvec[nr]) / d;
   }

   else if(this->thesolver->isId(id))
   {
      pvec[nr] = this->thesolver->vector(nr) * cvec;

      if(this->thesolver->isBasic(nr))
      {
         pupd.clearIdx(nr);
         return true;
      }

      x = pvec[nr];
      d = pupd[nr];
      up = &upb;
      low = &lpb;

      if(d > 0.0)
         sel = (this->thesolver->lpBound()[nr] - pvec[nr]) / d;
      else
         sel = (this->thesolver->upBound()[nr] - pvec[nr]) / d;
   }

   else
      return true;

   if((*up)[nr] != (*low)[nr])
   {
      if(sel > fastDelta / maxabs)
      {
         sel = 0.0;

         // prevent shifts in polishing mode to avoid a final cleanup step (i.e. simplex type switch)
         if(!polish)
         {
            if(d < 0.0)
            {
               this->thesolver->theShift -= (*up)[nr];
               (*up)[nr] = x;
               this->thesolver->theShift += (*up)[nr];
            }
            else
            {
               this->thesolver->theShift += (*low)[nr];
               (*low)[nr] = x;
               this->thesolver->theShift -= (*low)[nr];
            }
         }
      }
   }
   else
   {
      sel = 0.0;

      // prevent shifts in polishing mode to avoid a final cleanup step (i.e. simplex type switch)
      if(!polish)
      {
         if(x > (*up)[nr])
            this->thesolver->theShift += x - (*up)[nr];
         else
            this->thesolver->theShift += (*low)[nr] - x;

         (*up)[nr] = (*low)[nr] = x;
      }
   }

   return false;
}

template <class R>
bool SPxFastRT<R>::shortEnter(
   const SPxId& enterId,
   int nr,
   R max,
   R maxabs) const
{
   R shortval = this->tolerances()->scaleAccordingToEpsilon(SOPLEX_SHORTVAL);

   if(this->thesolver->isCoId(enterId))
   {
      if(max != 0.0)
      {
         R x = this->thesolver->coPvec().delta()[nr];

         if(x < maxabs * shortval && -x < maxabs * shortval)
            return false;
      }

      return true;
   }
   else if(this->thesolver->isId(enterId))
   {
      if(max != 0.0)
      {
         R x = this->thesolver->pVec().delta()[nr];

         if(x < maxabs * shortval && -x < maxabs * shortval)
            return false;
      }

      return true;
   }

   return false;
}

template <class R>
SPxId SPxFastRT<R>::selectEnter(R& val, int, bool polish)
{
   SPxId enterId;
   R max, sel;
   R maxabs = 0.0;
   R delta_shift = this->tolerances()->scaleAccordingToEpsilon(SOPLEX_DELTA_SHIFT);
   int nr;
   int cnt = 0;

   assert(this->m_type == SPxSolverBase<R>::LEAVE);

   resetTols();

   // force instable pivot iff true (see explanation in leave.hpp and spxsolve.hpp)
   bool instable = this->solver()->instableLeave;
   R lowstab = this->tolerances()->scaleAccordingToEpsilon(SOPLEX_LOWSTAB);
   assert(!instable || this->solver()->instableLeaveNum >= 0);

   sel = 0.0;

   if(val > epsilonZero())
   {
      do
      {
         maxabs = 0.0;
         max = val;

         enterId = maxDelta(nr, max, maxabs);

         if(!enterId.isValid())
            return enterId;

         assert(max >= 0.0);
         assert(!enterId.isValid() || !this->solver()->isBasic(enterId));

         if(!shortEnter(enterId, nr, max, maxabs))
         {
            R bestDelta, stab;

            stab = minStability(maxabs);

            // force instable pivot iff instable is true (see explanation in leave.hpp and spxsolve.hpp)
            if(instable)
            {
               enterId = maxSelect(nr, sel, lowstab, bestDelta, max);
            }
            else
            {
               enterId = maxSelect(nr, sel, stab, bestDelta, max);
            }

            if(bestDelta < delta_shift * SOPLEX_TRIES)
               cnt++;
            else
               cnt += SOPLEX_TRIES;
         }

         if(!maxReEnter(sel, maxabs, enterId, nr, polish))
            break;

         relax();
      }
      while(cnt < SOPLEX_TRIES);
   }
   else if(val < -epsilonZero())
   {
      do
      {
         maxabs = 0.0;
         max = val;
         enterId = minDelta(nr, max, maxabs);

         if(!enterId.isValid())
            return enterId;

         assert(max <= 0.0);
         assert(!enterId.isValid() || !this->solver()->isBasic(enterId));

         if(!shortEnter(enterId, nr, max, maxabs))
         {
            R bestDelta, stab;

            stab = minStability(maxabs);

            // force instable pivot iff instable is true (see explanation in leave.hpp and spxsolve.hpp)
            if(instable)
            {
               enterId = minSelect(nr, sel, lowstab, bestDelta, max);
            }
            else
            {
               enterId = minSelect(nr, sel, stab, bestDelta, max);
            }

            if(bestDelta < delta_shift * SOPLEX_TRIES)
               cnt++;
            else
               cnt += SOPLEX_TRIES;
         }

         if(!minReEnter(sel, maxabs, enterId, nr, polish))
            break;

         relax();
      }
      while(cnt < SOPLEX_TRIES);
   }

   SPX_DEBUG(

      if(enterId.isValid())
{
   assert(!enterId.isValid() || !this->solver()->isBasic(enterId));

      R x;

      if(this->thesolver->isCoId(enterId))
         x = this->thesolver->coPvec().delta()[ this->thesolver->number(enterId) ];
      else
         x = this->thesolver->pVec().delta()[ this->thesolver->number(enterId) ];

      std::cout << "DFSTRT03 " << this->thesolver->basis().iteration() << ": "
                << sel << '\t' << x << " (" << maxabs << ")" << std::endl;
   }
   else
      std::cout << "DFSTRT04 " << this->thesolver->basis().iteration()
                << ": skipping instable pivot" << std::endl;
      )

      if(polish && enterId.isValid())
      {
         assert(this->thesolver->rep() == SPxSolverBase<R>::ROW);

         // decide whether the chosen entering index contributes to the polishing objective
         if(this->thesolver->polishObj == SPxSolverBase<R>::POLISH_INTEGRALITY)
         {
            // only allow (integer) variables to enter the basis
            if(enterId.isSPxRowId())
               return SPxId();
            else if(this->thesolver->integerVariables.size() == this->thesolver->nCols()
                    && this->thesolver->integerVariables[this->thesolver->number(enterId)] == 0)
               return SPxId();
         }
         else if(this->thesolver->polishObj == SPxSolverBase<R>::POLISH_FRACTIONALITY)
         {
            // only allow slacks and continuous variables to enter the basis
            if(this->thesolver->integerVariables.size() == this->thesolver->nCols())
            {
               if(enterId.isSPxColId() && this->thesolver->integerVariables[this->thesolver->number(enterId)] == 1)
                  return SPxId();
            }
            else if(enterId.isSPxColId())
               return SPxId();
         }
      }


   if(enterId.isValid() || minStab > 2 * epsilonZero())
   {
      val = sel;

      if(enterId.isValid())
         tighten();
   }

   assert(!enterId.isValid() || !this->solver()->isBasic(enterId));

   return enterId;
}

template <class R>
void SPxFastRT<R>::load(SPxSolverBase<R>* spx)
{
   this->thesolver = spx;
   setType(spx->type());
}

template <class R>
void SPxFastRT<R>::setType(typename SPxSolverBase<R>::Type type)
{
   this->m_type = type;

   minStab = this->tolerances()->scaleAccordingToEpsilon(SOPLEX_MINSTAB);
   fastDelta = this->delta;
}
} // namespace soplex
