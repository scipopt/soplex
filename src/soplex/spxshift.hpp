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
#include <iostream>

#include "soplex/spxdefines.h"
#include "soplex/spxsolver.h"
#include "soplex/spxout.h"

namespace soplex
{
template <class R>
void SPxSolverBase<R>::shiftFvec()
{

   /* the allowed tolerance is (rep() == COLUMN) ? tolerances()->floatingPointFeastol() : tolerances()->floatingPointOpttol() because theFvec is the primal VectorBase<R> in COLUMN
    * and the dual VectorBase<R> in ROW representation; this is equivalent to entertol()
    */
   R minrandom = 10.0 * entertol();
   R maxrandom = 100.0 * entertol();
   R allow = entertol() - epsilon();

   assert(type() == ENTER);
   assert(allow > 0);

   for(int i = dim() - 1; i >= 0; --i)
   {
      if(theUBbound[i] + allow < (*theFvec)[i])
      {
         SPxOut::debug(this, "DSHIFT08 theUBbound[{}] violated by {}", i,
                       (*theFvec)[i] - theUBbound[i] - allow);

         if(theUBbound[i] != theLBbound[i])
         {
            // since minrandom and maxrandom are of the order 10 different,
            // we currently doesn't care about higher precision random
            // numbers. Hence the cast to double.
            shiftUBbound(i, (*theFvec)[i] + random.next((double)minrandom, (double)maxrandom));
         }
         else
         {
            shiftUBbound(i, (*theFvec)[i]);
            theLBbound[i] = theUBbound[i];
         }
      }
      else if((*theFvec)[i] < theLBbound[i] - allow)
      {
         SPxOut::debug(this, "DSHIFT08 theLBbound[{}] violated by {}", i,
                       theLBbound[i] - (*theFvec)[i] - allow);

         if(theUBbound[i] != theLBbound[i])
            shiftLBbound(i, (*theFvec)[i] - random.next((double)minrandom, (double)maxrandom));
         else
         {
            shiftLBbound(i, (*theFvec)[i]);
            theUBbound[i] = theLBbound[i];
         }
      }
   }

#ifndef NDEBUG
   testBounds();
   SPxOut::debug(this, "DSHIFT01 shiftFvec: OK\n");
#endif
}

// -----------------------------------------------------------------

/*
  This methods assumes correctly setup vectors |pVec| and |coPvec| and bound
  vectors for leaving simplex. Then it checks all values of |pVec| and
  |coPvec| to obey these bounds and enlarges them if neccessary.
*/
template <class R>
void SPxSolverBase<R>::shiftPvec()
{

   /* the allowed tolerance is (rep() == ROW) ? tolerances()->floatingPointFeastol() : tolerances()->floatingPointOpttol() because thePvec is the primal VectorBase<R> in ROW and the
    * dual VectorBase<R> in COLUMN representation; this is equivalent to leavetol()
    */
   R minrandom = 10.0 * leavetol();
   R maxrandom = 100.0 * leavetol();
   R allow = leavetol() - epsilon();
   bool tmp;
   int i;

   assert(type() == LEAVE);
   assert(allow > 0.0);

   for(i = dim() - 1; i >= 0; --i)
   {
      tmp = !isBasic(coId(i));

      if((*theCoUbound)[i] + allow <= (*theCoPvec)[i] && tmp)
      {
         if((*theCoUbound)[i] != (*theCoLbound)[i])
            shiftUCbound(i, (*theCoPvec)[i] + random.next((double)minrandom, (double)maxrandom));
         else
         {
            shiftUCbound(i, (*theCoPvec)[i]);
            (*theCoLbound)[i] = (*theCoUbound)[i];
         }
      }
      else if((*theCoLbound)[i] - allow >= (*theCoPvec)[i] && tmp)
      {
         if((*theCoUbound)[i] != (*theCoLbound)[i])
            shiftLCbound(i, (*theCoPvec)[i] - random.next((double)minrandom, (double)maxrandom));
         else
         {
            shiftLCbound(i, (*theCoPvec)[i]);
            (*theCoUbound)[i] = (*theCoLbound)[i];
         }
      }
   }

   for(i = coDim() - 1; i >= 0; --i)
   {
      tmp = !isBasic(id(i));

      if((*theUbound)[i] + allow <= (*thePvec)[i] && tmp)
      {
         if((*theUbound)[i] != (*theLbound)[i])
            shiftUPbound(i, (*thePvec)[i] + random.next((double)minrandom, (double)maxrandom));
         else
         {
            shiftUPbound(i, (*thePvec)[i]);
            (*theLbound)[i] = (*theUbound)[i];
         }
      }
      else if((*theLbound)[i] - allow >= (*thePvec)[i] && tmp)
      {
         if((*theUbound)[i] != (*theLbound)[i])
            shiftLPbound(i, (*thePvec)[i] - random.next((double)minrandom, (double)maxrandom));
         else
         {
            shiftLPbound(i, (*thePvec)[i]);
            (*theUbound)[i] = (*theLbound)[i];
         }
      }
   }

#ifndef NDEBUG
   testBounds();
   SPxOut::debug(this, "DSHIFT02 shiftPvec: OK\n");
#endif
}
// -----------------------------------------------------------------
template <class R>
void SPxSolverBase<R>::perturbMin(
   const UpdateVector<R>& uvec,
   VectorBase<R>& p_low,
   VectorBase<R>& p_up,
   R eps,
   R p_delta,
   int start,
   int incr)
{
   assert(uvec.dim() == p_low.dim());
   assert(uvec.dim() == p_up.dim());

   const R* vec = uvec.get_const_ptr();
   R minrandom = 10.0 * p_delta;
   R maxrandom = 100.0 * p_delta;
   R x, l, u;
   int i;

   if(fullPerturbation)
   {
      eps = p_delta;

      for(i = uvec.dim() - start - 1; i >= 0; i -= incr)
      {
         u = p_up[i];
         l = p_low[i];
         x = vec[i];

         if(LT(u, R(infinity), eps) && NE(l, u, eps) && u <= x + eps)
         {
            p_up[i] = x + random.next((double) minrandom, (double)maxrandom);
            theShift += p_up[i] - u;
         }

         if(GT(l, R(-infinity), eps) && NE(l, u, eps) && l >= x - eps)
         {
            p_low[i] = x - random.next((double)minrandom, (double)maxrandom);
            theShift -= p_low[i] - l;
         }
      }
   }
   else
   {
      const R* upd = uvec.delta().values();
      const IdxSet& idx = uvec.delta().indices();

      for(int j = uvec.delta().size() - start - 1; j >= 0; j -= incr)
      {
         i = idx.index(j);
         x = upd[i];
         u = p_up[i];
         l = p_low[i];

         // do not permute these bounds! c.f. with computeFrhs2() in spxvecs.cpp
         if(this->dualStatus(this->baseId(i)) == SPxBasisBase<R>::Desc::D_ON_BOTH)
         {
            continue;
         }

         if(x < -eps)
         {
            if(LT(u, R(infinity), eps) && NE(l, u, eps) && vec[i] >= u - eps)
            {
               p_up[i] = vec[i] + random.next((double)minrandom, (double)maxrandom);
               theShift += p_up[i] - u;
            }
         }
         else if(x > eps)
         {
            if(GT(l, R(-infinity), eps) && NE(l, u, eps) && vec[i] <= l + eps)
            {
               p_low[i] = vec[i] - random.next((double)minrandom, (double)maxrandom);
               theShift -= p_low[i] - l;
            }
         }
      }
   }
}
// -----------------------------------------------------------------
template <class R>
void SPxSolverBase<R>::perturbMax(
   const UpdateVector<R>& uvec,
   VectorBase<R>& p_low,
   VectorBase<R>& p_up,
   R eps,
   R p_delta,
   int start,
   int incr)
{
   assert(uvec.dim() == p_low.dim());
   assert(uvec.dim() == p_up.dim());

   const R* vec = uvec.get_const_ptr();
   R minrandom = 10.0 * p_delta;
   R maxrandom = 100.0 * p_delta;
   R x, l, u;
   int i;

   if(fullPerturbation)
   {
      eps = p_delta;

      for(i = uvec.dim() - start - 1; i >= 0; i -= incr)
      {
         u = p_up[i];
         l = p_low[i];
         x = vec[i];

         if(LT(u, R(infinity), eps) && NE(l, u, eps) && u <= x + eps)
         {
            p_up[i] = x + random.next((double)minrandom, (double)maxrandom);
            theShift += p_up[i] - u;
         }

         if(GT(l, R(-infinity), eps) && NE(l, u, eps) && l >= x - eps)
         {
            p_low[i] = x - random.next((double)minrandom, (double)maxrandom);
            theShift -= p_low[i] - l;
         }
      }
   }
   else
   {
      const R* upd = uvec.delta().values();
      const IdxSet& idx = uvec.delta().indices();

      for(int j = uvec.delta().size() - start - 1; j >= 0; j -= incr)
      {
         i = idx.index(j);
         x = upd[i];
         u = p_up[i];
         l = p_low[i];

         // do not permute these bounds! c.f. computeFrhs2() in spxvecs.cpp
         if(this->dualStatus(this->baseId(i)) == SPxBasisBase<R>::Desc::D_ON_BOTH)
         {
            continue;
         }

         if(x > eps)
         {
            if(LT(u, R(infinity), eps) && NE(l, u, eps) && vec[i] >= u - eps)
            {
               p_up[i] = vec[i] + random.next((double)minrandom, (double)maxrandom);
               theShift += p_up[i] - u;
            }
         }
         else if(x < -eps)
         {
            if(GT(l, R(-infinity), eps) && NE(l, u, eps) && vec[i] <= l + eps)
            {
               p_low[i] = vec[i] - random.next((double)minrandom, (double)maxrandom);
               theShift -= p_low[i] - l;
            }
         }
      }
   }
}

template <class R>
void SPxSolverBase<R>::perturbMinEnter(void)
{
   SPxOut::debug(this, "DSHIFT03 iteration= {} perturbing {}", this->iteration(), shift());
   fVec().delta().setup();
   perturbMin(fVec(), lbBound(), ubBound(), epsilon(), entertol());
   SPxOut::debug(this, "\t->{}\n", shift());
}


template <class R>
void SPxSolverBase<R>::perturbMaxEnter(void)
{
   SPxOut::debug(this, "DSHIFT04 iteration= {} perturbing {}", this->iteration(), shift());
   fVec().delta().setup();
   perturbMax(fVec(), lbBound(), ubBound(), epsilon(), entertol());
   SPxOut::debug(this, "\t->{}\n", shift());
}


template <class R>
R SPxSolverBase<R>::perturbMin(
   const UpdateVector<R>& uvec,
   VectorBase<R>& p_low,
   VectorBase<R>& p_up,
   R eps,
   R p_delta,
   const typename SPxBasisBase<R>::Desc::Status* stat,
   int start,
   int incr)
{
   assert(uvec.dim() == p_low.dim());
   assert(uvec.dim() == p_up.dim());

   const R* vec = uvec.get_const_ptr();
   R minrandom = 10.0 * p_delta;
   R maxrandom = 100.0 * p_delta;
   R x, l, u;
   int i;
   R l_theShift = 0;

   if(fullPerturbation)
   {
      eps = p_delta;

      for(i = uvec.dim() - start - 1; i >= 0; i -= incr)
      {
         u = p_up[i];
         l = p_low[i];
         x = vec[i];

         if(LT(u, R(infinity), eps) && NE(l, u, eps) && u <= x + eps && rep() * stat[i] < 0)
         {
            p_up[i] = vec[i] + random.next((double)minrandom, (double)maxrandom);
            l_theShift += p_up[i] - u;
         }

         if(GT(l, R(-infinity), eps) && NE(l, u, eps) && l >= x - eps && rep() * stat[i] < 0)
         {
            p_low[i] = vec[i] - random.next((double)minrandom, (double)maxrandom);
            l_theShift -= p_low[i] - l;
         }
      }
   }
   else
   {
      const R* upd = uvec.delta().values();
      const IdxSet& idx = uvec.delta().indices();

      for(int j = uvec.delta().size() - start - 1; j >= 0; j -= incr)
      {
         i = idx.index(j);
         x = upd[i];
         u = p_up[i];
         l = p_low[i];

         if(x < -eps)
         {
            if(LT(u, R(infinity), eps) && NE(l, u, eps) && vec[i] >= u - eps && rep() * stat[i] < 0)
            {
               p_up[i] = vec[i] + random.next((double)minrandom, (double)maxrandom);
               l_theShift += p_up[i] - u;
            }
         }
         else if(x > eps)
         {
            if(GT(l, R(-infinity), eps) && NE(l, u, eps) && vec[i] <= l + eps && rep() * stat[i] < 0)
            {
               p_low[i] = vec[i] - random.next((double)minrandom, (double)maxrandom);
               l_theShift -= p_low[i] - l;
            }
         }
      }
   }

   return l_theShift;
}

template <class R>
R SPxSolverBase<R>::perturbMax(
   const UpdateVector<R>& uvec,
   VectorBase<R>& p_low,
   VectorBase<R>& p_up,
   R eps,
   R p_delta,
   const typename SPxBasisBase<R>::Desc::Status* stat,
   int start,
   int incr)
{
   assert(uvec.dim() == p_low.dim());
   assert(uvec.dim() == p_up.dim());

   const R* vec = uvec.get_const_ptr();
   R minrandom = 10.0 * p_delta;
   R maxrandom = 100.0 * p_delta;
   R x, l, u;
   int i;
   R l_theShift = 0;

   if(fullPerturbation)
   {
      eps = p_delta;

      for(i = uvec.dim() - start - 1; i >= 0; i -= incr)
      {
         u = p_up[i];
         l = p_low[i];
         x = vec[i];

         if(LT(u, R(infinity), eps) && NE(l, u, eps) && u <= x + eps && rep() * stat[i] < 0)
         {
            p_up[i] = vec[i] + random.next((double)minrandom, (double)maxrandom);
            l_theShift += p_up[i] - u;
         }

         if(GT(l, R(-infinity), eps) && NE(l, u, eps) && l >= x - eps && rep() * stat[i] < 0)
         {
            p_low[i] = vec[i] - random.next((double)minrandom, (double)maxrandom);
            l_theShift -= p_low[i] - l;
         }
      }
   }
   else
   {
      const R* upd = uvec.delta().values();
      const IdxSet& idx = uvec.delta().indices();

      for(int j = uvec.delta().size() - start - 1; j >= 0; j -= incr)
      {
         i = idx.index(j);
         x = upd[i];
         u = p_up[i];
         l = p_low[i];

         if(x > eps)
         {
            if(LT(u, R(infinity), eps) && NE(l, u, eps) && vec[i] >= u - eps && rep() * stat[i] < 0)
            {
               p_up[i] = vec[i] + random.next((double)minrandom, (double)maxrandom);
               l_theShift += p_up[i] - u;
            }
         }
         else if(x < -eps)
         {
            if(GT(l, R(-infinity), eps) && NE(l, u, eps) && vec[i] <= l + eps && rep() * stat[i] < 0)
            {
               p_low[i] = vec[i] - random.next((double)minrandom, (double)maxrandom);
               l_theShift -= p_low[i] - l;
            }
         }
      }
   }

   return l_theShift;
}


template <class R>
void SPxSolverBase<R>::perturbMinLeave(void)
{
   SPxOut::debug(this, "DSHIFT05 iteration= {} perturbing {}", this->iteration(), shift());
   pVec().delta().setup();
   coPvec().delta().setup();
   theShift += perturbMin(pVec(), lpBound(), upBound(), epsilon(), leavetol(),
                          this->desc().status(), 0, 1);
   theShift += perturbMin(coPvec(), lcBound(), ucBound(), epsilon(), leavetol(),
                          this->desc().coStatus(), 0, 1);
   SPxOut::debug(this, "\t->{}\n", shift());
}


template <class R>
void SPxSolverBase<R>::perturbMaxLeave(void)
{
   SPxOut::debug(this, "DSHIFT06 iteration= {} perturbing {}", this->iteration(), shift());
   pVec().delta().setup();
   coPvec().delta().setup();
   theShift += perturbMax(pVec(), lpBound(), upBound(), epsilon(), leavetol(),
                          this->desc().status(), 0, 1);
   theShift += perturbMax(coPvec(), lcBound(), ucBound(), epsilon(), leavetol(),
                          this->desc().coStatus(), 0, 1);
   SPxOut::debug(this, "\t->{}\n", shift());
}


template <class R>
void SPxSolverBase<R>::unShift(void)
{
   SPX_MSG_INFO3((*this->spxout), (*this->spxout) << "DSHIFT07 = " << "unshifting ..." << std::endl;);

   if(isInitialized())
   {
      int i;
      R t_up, t_low;
      const typename SPxBasisBase<R>::Desc& ds = this->desc();

      theShift = 0;

      if(type() == ENTER)
      {
         R eps = entertol();

         if(rep() == COLUMN)
         {
            for(i = dim(); i-- > 0;)
            {
               SPxId l_id = this->baseId(i);
               int l_num = this->number(l_id);

               if(l_id.type() == SPxId::ROW_ID)
               {
                  t_up = -this->lhs(l_num);
                  t_low = -this->rhs(l_num);
               }
               else
               {
                  assert(l_id.type() == SPxId::COL_ID);
                  t_up = this->upper(l_num);
                  t_low = this->lower(l_num);
               }

               if(t_up != t_low)
               {
                  if((*theFvec)[i] < t_up + eps)  // check allowed violation
                     theUBbound[i] = t_up; // reset shifted bound to original
                  else if((*theFvec)[i] > t_up)  // shifted bound is required for feasibility
                     theShift += theUBbound[i] - t_up;

                  if((*theFvec)[i] > t_low - eps)  // check allowed violation
                     theLBbound[i] = t_low; // reset shifted bound to original
                  else if((*theFvec)[i] < t_low)  // shifted bound is required for feasibility
                     theShift -= theLBbound[i] - t_low;
               }
               else
               {
                  if(theUBbound[i] > t_up)
                     theShift += theUBbound[i] - t_up;
                  else if(theLBbound[i] < t_low)
                     theShift += t_low - theLBbound[i];
               }
            }

            for(i = this->nRows(); i-- > 0;)
            {
               if(!isBasic(ds.rowStatus(i)))
               {
                  t_up = -this->lhs(i);
                  t_low = -this->rhs(i);

                  if(theURbound[i] > t_up)  // what about t_up == t_low ?
                     theShift += theURbound[i] - t_up;

                  if(t_low > theLRbound[i])  // what about t_up == t_low ?
                     theShift += t_low - theLRbound[i];
               }
            }

            for(i = this->nCols(); i-- > 0;)
            {
               if(!isBasic(ds.colStatus(i)))
               {
                  t_up = this->upper(i);
                  t_low = this->lower(i);

                  if(theUCbound[i] > t_up)  // what about t_up == t_low ?
                     theShift += theUCbound[i] - t_up;

                  if(t_low > theLCbound[i])  // what about t_up == t_low ?
                     theShift += t_low - theLCbound[i];
               }
            }
         }
         else
         {
            assert(rep() == ROW);

            for(i = dim(); i-- > 0;)
            {
               SPxId l_id = this->baseId(i);
               int l_num = this->number(l_id);
               t_up = t_low = 0;

               if(l_id.type() == SPxId::ROW_ID)
                  clearDualBounds(ds.rowStatus(l_num), t_up, t_low);
               else
                  clearDualBounds(ds.colStatus(l_num), t_up, t_low);

               if(theUBbound[i] != theLBbound[i])
               {
                  if(theUBbound[i] > t_up)
                     theShift += theUBbound[i] - t_up;
                  else
                     theShift -= theUBbound[i] - t_up;
               }
               else
               {
                  /* if the basic (primal or dual) variable is fixed (e.g., basis status P_FREE in row representation)
                   * then shiftFvec() and shiftPvec() do not relax the bounds, but shift both, hence they may be outside
                   * of [t_low,t_up] */
                  assert(theLBbound[i] == theUBbound[i] || theUBbound[i] >= t_up);
                  assert(theLBbound[i] == theUBbound[i] || theLBbound[i] <= t_low);

                  if((*theFvec)[i] < t_up - eps)
                     theUBbound[i] = t_up;
                  else if((*theFvec)[i] > t_up)
                     theShift += theUBbound[i] - t_up;

                  if((*theFvec)[i] > t_low + eps)
                     theLBbound[i] = t_low;
                  else if((*theFvec)[i] < t_low)
                     theShift -= theLBbound[i] - t_low;
               }
            }

            for(i = this->nRows(); i-- > 0;)
            {
               if(!isBasic(ds.rowStatus(i)))
               {
                  t_up = t_low = 0;
                  clearDualBounds(ds.rowStatus(i), t_up, t_low);

                  if(theURbound[i] > t_up)  // what about t_up == t_low ?
                     theShift += theURbound[i] - t_up;

                  if(t_low > theLRbound[i])  // what about t_up == t_low ?
                     theShift += t_low - theLRbound[i];
               }
            }

            for(i = this->nCols(); i-- > 0;)
            {
               if(!isBasic(ds.colStatus(i)))
               {
                  t_up = t_low = 0;
                  clearDualBounds(ds.colStatus(i), t_up, t_low);

                  if(theUCbound[i] > t_up)  // what about t_up == t_low ?
                     theShift += theUCbound[i] - t_up;

                  if(t_low > theLCbound[i])  // what about t_up == t_low ?
                     theShift += t_low - theLCbound[i];
               }
            }
         }
      }
      else
      {
         assert(type() == LEAVE);

         R eps = leavetol();

         if(rep() == COLUMN)
         {
            for(i = this->nRows(); i-- > 0;)
            {
               t_up = t_low = this->maxRowObj(i);
               clearDualBounds(ds.rowStatus(i), t_up, t_low);

               if(!isBasic(ds.rowStatus(i)))
               {
                  if((*theCoPvec)[i] < t_up + eps)
                  {
                     theURbound[i] = t_up; // reset bound to original value

                     if(t_up == t_low)
                        theLRbound[i] = t_low; // for fixed rows we change both bounds
                  }
                  else
                     theShift += theURbound[i] - t_up;

                  if((*theCoPvec)[i] > t_low - eps)
                  {
                     theLRbound[i] = t_low; // reset bound to original value

                     if(t_up == t_low)
                        theURbound[i] = t_up; // for fixed rows we change both bounds
                  }
                  else
                     theShift += t_low - theLRbound[i];
               }
               else if(theURbound[i] > t_up)
                  theShift += theURbound[i] - t_up;
               else if(theLRbound[i] < t_low)
                  theShift += t_low - theLRbound[i];
            }

            for(i = this->nCols(); i-- > 0;)
            {
               t_up = t_low = -this->maxObj(i);
               clearDualBounds(ds.colStatus(i), t_low, t_up);

               if(!isBasic(ds.colStatus(i)))
               {
                  if((*thePvec)[i] < -t_up + eps)
                  {
                     theUCbound[i] = -t_up; // reset bound to original value

                     if(t_up == t_low)
                        theLCbound[i] = -t_low; // for fixed variables we change both bounds
                  }
                  else
                     theShift += theUCbound[i] - (-t_up);

                  if((*thePvec)[i] > -t_low - eps)
                  {
                     theLCbound[i] = -t_low; // reset bound to original value

                     if(t_up == t_low)
                        theUCbound[i] = -t_up; // for fixed variables we change both bounds
                  }
                  else
                     theShift += (-t_low) - theLCbound[i];
               }
               else if(theUCbound[i] > -t_up)
                  theShift += theUCbound[i] - (-t_up);
               else if(theLCbound[i] < -t_low)
                  theShift += (-t_low) - theLCbound[i];
            }
         }
         else
         {
            assert(rep() == ROW);

            for(i = this->nRows(); i-- > 0;)
            {
               t_up = this->rhs(i);
               t_low = this->lhs(i);

               if(t_up == t_low)
               {
                  if(theURbound[i] > t_up)
                     theShift += theURbound[i] - t_up;
                  else
                     theShift += t_low - theLRbound[i];
               }
               else if(!isBasic(ds.rowStatus(i)))
               {
                  if((*thePvec)[i] < t_up + eps)
                     theURbound[i] = t_up; // reset bound to original value
                  else
                     theShift += theURbound[i] - t_up;

                  if((*thePvec)[i] > t_low - eps)
                     theLRbound[i] = t_low; // reset bound to original value
                  else
                     theShift += t_low - theLRbound[i];
               }
               else if(theURbound[i] > t_up)
                  theShift += theURbound[i] - t_up;
               else if(theLRbound[i] < t_low)
                  theShift += t_low - theLRbound[i];
            }

            for(i = this->nCols(); i-- > 0;)
            {
               t_up = this->upper(i);
               t_low = this->lower(i);

               if(t_up == t_low)
               {
                  if(theUCbound[i] > t_up)
                     theShift += theUCbound[i] - t_up;
                  else
                     theShift += t_low - theLCbound[i];
               }
               else if(!isBasic(ds.colStatus(i)))
               {
                  if((*theCoPvec)[i] < t_up + eps)
                     theUCbound[i] = t_up; // reset bound to original value
                  else
                     theShift += theUCbound[i] - t_up;

                  if((*theCoPvec)[i] > t_low - eps)
                     theLCbound[i] = t_low; // reset bound to original value
                  else
                     theShift += t_low - theLCbound[i];
               }
               else if(theUCbound[i] > t_up)
                  theShift += theUCbound[i] - t_up;
               else if(theLCbound[i] < t_low)
                  theShift += t_low - theLCbound[i];
            }
         }
      }
   }
}
} // namespace soplex
