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
#pragma ident "@(#) $Id: spxfastrt.cpp,v 1.1 2001/11/06 16:18:32 bzfkocht Exp $"

/*      \Section{Complex Methods}
 */

/*  Import system include files
 */
#include <assert.h>
#include <stdio.h>
#include <iostream>


/*  and class header files
 */
#include "spxfastrt.h"

namespace soplex
{

//@ #define     MINSTAB thesolver->delta()
#define MINSTAB         1e-5
#define TRIES           2
#define SHORT           1e-5
#define DELTA_SHIFT     1e-5
#define EPSILON         1e-10

//@ ----------------------------------------------------------------------------

void SPxFastRT::resetTols()
{
   // epsilon = thesolver->epsilon();
   epsilon = EPSILON;
   /*
       if(thesolver->basis().stability() < 1e-4)
           epsilon *= 1e-4 / thesolver->basis().stability();
       std::cerr << "epsilon = " << epsilon << '\t';
       std::cerr << "delta   = " << delta   << '\t';
       std::cerr << "minStab = " << minStab << std::endl;
    */
}

void SPxFastRT::tighten()
{
   /*
   if((delta > 1.99 * DELTA_SHIFT  &&  thesolver->theShift < 1e-4) ||
       (delta > 1e-4   &&  thesolver->theShift > 1e-4))
    */
   // if(delta > 1.99 * DELTA_SHIFT)
   if (delta >= delta0 + DELTA_SHIFT)
   {
      delta -= DELTA_SHIFT;
      if (delta > 1e-4)
         delta -= 2 * DELTA_SHIFT;
   }

   if (minStab < MINSTAB)
   {
      minStab /= 0.90;
      if (minStab < 1e-6)
         minStab /= 0.90;
   }
   /*
    */
}

void SPxFastRT::relax()
{
   minStab *= 0.95;
   delta += 3 * DELTA_SHIFT;
   // delta   += 2 * (thesolver->theShift > delta) * DELTA_SHIFT;
   //@ std::cerr << '\t' << minStab << '\t' << delta << std::endl;
}


static double minStability(double minStab, double maxabs)
{
   if (maxabs < 1000)
      return minStab;
   return maxabs*minStab / 1000;
}



//@ ----------------------------------------------------------------------------

int SPxFastRT::maxDelta(
   double& val,
   double& abs,
   UpdateVector& update,
   Vector& lowBound,
   Vector& upBound,
   int start,
   int incr
)
{
   int i, sel;
   double x, y, max;
   double u, l;

   double delta = this->delta;
   // double           delta01 = 0.5*delta;
   double delta01 = 0;
   double inf = SPxLP::infinity;
   double mabs = abs;

   double* up = upBound.get_ptr();
   double* low = lowBound.get_ptr();
   const double* vec = update.get_const_ptr();
   const double* upd = update.delta().values();
   const int* idx = update.delta().indexMem();

   sel = -1;
   max = val;

   if (update.delta().isSetup())
   {
      const int* last = idx + update.delta().size();
      for (idx += start; idx < last; idx += incr)
      {
         i = *idx;
         x = upd[i];
         if (x > epsilon)
         {
            mabs = (x > mabs) ? x : mabs;
            u = up[i];
            if (u < inf)
            {
               y = u - vec[i];
               // x = ((1 - (y<=0)) * y + delta) / x;
               x = (y - (y <= 0) * (y + delta01) + delta) / x;
               if (x < max)
               {
                  max = x;
                  sel = i;
               }
            }
         }
         else if (x < -epsilon)
         {
            mabs = (-x > mabs) ? -x : mabs;
            l = low[i];
            if (l > -inf)
            {
               y = l - vec[i];
               // x = ((1 - (y>=0)) * y - delta) / x;
               x = (y - (y >= 0) * (y - delta01) - delta) / x;
               if (x < max)
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
      int* idx = update.delta().altIndexMem();
      double* uval = update.delta().altValues();
      const double* uend = uval + update.dim();

      for (; uval < uend; ++uval)
      {
         if (*uval)
         {
            x = *uval;
            i = uval - upd;
            if (x > epsilon)
            {
               *idx++ = i;
               mabs = (x > mabs) ? x : mabs;
               u = up[i];
               if (u < inf)
               {
                  y = u - vec[i];
                  // x = ((1 - (y<=0)) * y + delta) / x;
                  x = (y - (y <= 0) * (y + delta01) + delta) / x;
                  if (x < max)
                  {
                     max = x;
                     sel = i;
                  }
               }
            }
            else if (x < -epsilon)
            {
               *idx++ = i;
               mabs = (-x > mabs) ? -x : mabs;
               l = low[i];
               if (l > -inf)
               {
                  y = l - vec[i];
                  // x = ((1 - (y>=0)) * y - delta) / x;
                  x = (y - (y >= 0) * (y - delta01) - delta) / x;
                  if (x < max)
                  {
                     max = x;
                     sel = i;
                  }
               }
            }
            else
               *uval = 0;
         }
      }
      update.delta().setSize(idx - update.delta().indexMem());
      update.delta().forceSetup();
   }

   val = max;
   abs = mabs;

   return sel;
}

int SPxFastRT::minDelta(
   double& val,
   double& abs,
   UpdateVector& update,
   Vector& lowBound,
   Vector& upBound,
   int start,
   int incr
)
{
   int i, sel;
   double x, y, max;
   double u, l;

   double delta = this->delta;
   // double           delta01 = 0.5*delta;
   double delta01 = 0;
   double inf = SPxLP::infinity;
   double mabs = abs;

   double* up = upBound.get_ptr();
   double* low = lowBound.get_ptr();
   const double* vec = update.get_const_ptr();
   const double* upd = update.delta().values();
   const int* idx = update.delta().indexMem();

   sel = -1;
   max = val;

   if (update.delta().isSetup())
   {
      const int* last = idx + update.delta().size();
      for (idx += start; idx < last; idx += incr)
      {
         i = *idx;
         x = upd[i];
         if (x > epsilon)
         {
            mabs = (x > mabs) ? x : mabs;
            l = low[i];
            if (l > -inf)
            {
               y = l - vec[i];
               // x = ((1 - (y>=0)) * y - delta) / x;
               x = (y - (y >= 0) * (y - delta01) - delta) / x;
               if (x > max)
               {
                  max = x;
                  sel = i;
               }
            }
         }
         else if (x < -epsilon)
         {
            mabs = (-x > mabs) ? -x : mabs;
            u = up[i];
            if (u < inf)
            {
               y = u - vec[i];
               // x = ((1 - (y<=0)) * y + delta) / x;
               x = (y - (y <= 0) * (y + delta01) + delta) / x;
               if (x > max)
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
      int* idx = update.delta().altIndexMem();
      double* uval = update.delta().altValues();
      const double* uend = uval + update.dim();

      for (; uval < uend; ++uval)
      {
         if (*uval)
         {
            x = *uval;
            i = uval - upd;
            if (x > epsilon)
            {
               *idx++ = i;
               mabs = (x > mabs) ? x : mabs;
               l = low[i];
               if (l > -inf)
               {
                  y = l - vec[i];
                  // x = ((1 - (y>=0)) * y - delta) / x;
                  x = (y - (y >= 0) * (y - delta01) - delta) / x;
                  if (x > max)
                  {
                     max = x;
                     sel = i;
                  }
               }
            }
            else if (x < -epsilon)
            {
               *idx++ = i;
               mabs = (-x > mabs) ? -x : mabs;
               u = up[i];
               if (u < inf)
               {
                  y = u - vec[i];
                  // x = ((1 - (y<=0)) * y + delta) / x;
                  x = (y - (y <= 0) * (y + delta01) + delta) / x;
                  if (x > max)
                  {
                     max = x;
                     sel = i;
                  }
               }
            }
            else
               *uval = 0;
         }
      }
      update.delta().setSize(idx - update.delta().indexMem());
      update.delta().forceSetup();
   }

   val = max;
   abs = mabs;

   return sel;
}

int SPxFastRT::maxDelta(
   double& val,
   double& abs)
{
   return maxDelta(val, abs,
                    thesolver->fVec(), thesolver->lbBound(), thesolver->ubBound(),
                    0, 1);
}

int SPxFastRT::minDelta(
   double& val,
   double& abs)
{
   return minDelta(val, abs,
                    thesolver->fVec(), thesolver->lbBound(), thesolver->ubBound(),
                    0, 1);
}

SoPlex::Id SPxFastRT::maxDelta(
   int& nr,
   double& max,
   double& maxabs)
{
   int indc, indp;
   indc = maxDelta(max, maxabs,
                    thesolver->coPvec(), thesolver->lcBound(), thesolver->ucBound(),
                    0, 1);
   indp = maxDelta(max, maxabs,
                    thesolver->pVec(), thesolver->lpBound(), thesolver->upBound(),
                    0, 1);

   if (indp >= 0)
   {
      nr = indp;
      return thesolver->id(indp);
   }
   if (indc >= 0)
   {
      nr = indc;
      return thesolver->coId(indc);
   }

   SoPlex::Id enterId;
   nr = -1;
   return enterId;
}

SoPlex::Id SPxFastRT::minDelta(
   int& nr,
   double& max,
   double& maxabs)
{
   int indc, indp;
   indc = minDelta(max, maxabs,
                    thesolver->coPvec(), thesolver->lcBound(), thesolver->ucBound(),
                    0, 1);
   indp = minDelta(max, maxabs,
                    thesolver->pVec(), thesolver->lpBound(), thesolver->upBound(),
                    0, 1);

   if (indp >= 0)
   {
      nr = indp;
      return thesolver->id(indp);
   }
   if (indc >= 0)
   {
      nr = indc;
      return thesolver->coId(indc);
   }

   SoPlex::Id enterId;
   nr = -1;
   return enterId;
}


int SPxFastRT::minSelect
(
   double& val,
   double& stab,
   double& best,
   double& bestDelta,
   double max,
   const UpdateVector& update,
   const Vector& lowBound,
   const Vector& upBound,
   int start,
   int incr
)
{
   int i;
   double x, y;

   const double* up = upBound.get_const_ptr();
   const double* low = lowBound.get_const_ptr();
   const double* vec = update.get_const_ptr();
   const double* upd = update.delta().values();
   const int* idx = update.delta().indexMem();
   const int* last = idx + update.delta().size();

   int nr = -1;
   int bestNr = -1;

   for (idx += start; idx < last; idx += incr)
   {
      i = *idx;
      x = upd[i];
      if (x > stab)
      {
         y = (low[i] - vec[i]) / x;
         if (y >= max)
         {
            val = y;
            nr = i;
            stab = x;
         }
         else if (y < best)
         {
            best = y;
            bestNr = i;
         }
      }
      else if (x < -stab)
      {
         y = (up[i] - vec[i]) / x;
         if (y >= max)
         {
            val = y;
            nr = i;
            stab = -x;
         }
         else if (y < best)
         {
            best = y;
            bestNr = i;
         }
      }
   }

   if (nr < 0 && bestNr > 0)
   {
      if (upd[bestNr] < 0)
         bestDelta = up[bestNr] - vec[bestNr];
      else
         bestDelta = vec[bestNr] - low[bestNr];
   }

   return nr;
}

int SPxFastRT::maxSelect
(
   double& val,
   double& stab,
   double& best,
   double& bestDelta,
   double max,
   const UpdateVector& update,
   const Vector& lowBound,
   const Vector& upBound,
   int start,
   int incr
)
{
   int i;
   double x, y;

   const double* up = upBound.get_const_ptr();
   const double* low = lowBound.get_const_ptr();
   const double* vec = update.get_const_ptr();
   const double* upd = update.delta().values();
   const int* idx = update.delta().indexMem();
   const int* last = idx + update.delta().size();

   int nr = -1;
   int bestNr = -1;

   for (idx += start; idx < last; idx += incr)
   {
      i = *idx;
      x = upd[i];
      if (x > stab)
      {
         y = (up[i] - vec[i]) / x;
         if (y <= max)
         {
            val = y;
            nr = i;
            stab = x;
         }
         else if (y > best)
         {
            best = y;
            bestNr = i;
         }
      }
      else if (x < -stab)
      {
         y = (low[i] - vec[i]) / x;
         if (y <= max)
         {
            val = y;
            nr = i;
            stab = -x;
         }
         else if (y > best)
         {
            best = y;
            bestNr = i;
         }
      }
   }

   if (nr < 0 && bestNr > 0)
   {
      if (upd[bestNr] > 0)
         bestDelta = up[bestNr] - vec[bestNr];
      else
         bestDelta = vec[bestNr] - low[bestNr];
   }

   return nr;
}


int SPxFastRT::maxSelect(
   double& val,
   double& stab,
   double& bestDelta,
   double max
)
{
   double best = -SPxLP::infinity;
   bestDelta = 0;
   return maxSelect(val, stab, best, bestDelta, max,
                     thesolver->fVec(), thesolver->lbBound(), thesolver->ubBound(),
                     0, 1);
}

SoPlex::Id SPxFastRT::maxSelect(
   int& nr,
   double& val,
   double& stab,
   double& bestDelta,
   double max
)
{
   int indp, indc;
   double best = -SPxLP::infinity;
   bestDelta = 0;
   indc = maxSelect(val, stab, best, bestDelta, max,
                     thesolver->coPvec(), thesolver->lcBound(), thesolver->ucBound(),
                     0, 1);
   indp = maxSelect(val, stab, best, bestDelta, max,
                     thesolver->pVec(), thesolver->lpBound(), thesolver->upBound(),
                     0, 1);

   if (indp >= 0)
   {
      nr = indp;
      return thesolver->id(indp);
   }
   if (indc >= 0)
   {
      nr = indc;
      return thesolver->coId(indc);
   }


   nr = -1;
   SoPlex::Id enterId;
   return enterId;
}

int SPxFastRT::minSelect(
   double& val,
   double& stab,
   double& bestDelta,
   double max
)
{
   double best = SPxLP::infinity;
   bestDelta = 0;
   return minSelect(val, stab, best, bestDelta, max,
                     thesolver->fVec(), thesolver->lbBound(), thesolver->ubBound(),
                     0, 1);
}

SoPlex::Id SPxFastRT::minSelect(
   int& nr,
   double& val,
   double& stab,
   double& bestDelta,
   double max
)
{
   int indp, indc;
   double best = SPxLP::infinity;
   bestDelta = 0;
   indc = minSelect(val, stab, best, bestDelta, max,
                     thesolver->coPvec(), thesolver->lcBound(), thesolver->ucBound(),
                     0, 1);
   indp = minSelect(val, stab, best, bestDelta, max,
                     thesolver->pVec(), thesolver->lpBound(), thesolver->upBound(),
                     0, 1);

   if (indp >= 0)
   {
      nr = indp;
      return thesolver->id(indp);
   }
   if (indc >= 0)
   {
      nr = indc;
      return thesolver->coId(indc);
   }


   nr = -1;
   SoPlex::Id enterId;
   return enterId;
}


//@ ----------------------------------------------------------------------------
/*
    Here comes our implementation of the Haris procedure improved by shifting
    bounds. The basic idea is to allow a slight infeasibility within |delta| to
    allow for more freedom when selecting the leaveing variable. This freedom
    may than be used for selecting numerical stable variables with great
    improves.
 
    The algorithms operates in two phases. In a first phase, the maximum |val|
    is determined, when in feasibility within |inf| is allowed. In the second
    phase, between all variables with values |< val| the one is selected which
    gives the best step forward in the simplex iteration. However, this may not
    allways yield an improvement. In that case, we shift the variable toward
    infeasibility and retry. This avoids cycling in the shifted LP.
 */
int SPxFastRT::maxShortLeave(double& sel, int leave, double max, double abs)
{
   assert(leave >= 0);
   sel = thesolver->fVec().delta()[leave];
   if (sel > abs*SHORT || -sel > abs*SHORT)
   {
      if (sel > 0)
         sel = (thesolver->ubBound()[leave] - thesolver->fVec()[leave]) / sel;
      else
         sel = (thesolver->lbBound()[leave] - thesolver->fVec()[leave]) / sel;
      //@ std::cerr << '+';
      return 1;
   }
   return 0;
}

int SPxFastRT::minShortLeave(double& sel, int leave, double max, double abs)
{
   assert(leave >= 0);
   sel = thesolver->fVec().delta()[leave];
   if (sel > abs*SHORT || -sel > abs*SHORT)
   {
      if (sel > 0)
         sel = (thesolver->lbBound()[leave] - thesolver->fVec()[leave]) / sel;
      else
         sel = (thesolver->ubBound()[leave] - thesolver->fVec()[leave]) / sel;
      //@ std::cerr << '+';
      return 1;
   }
   return 0;
}

int SPxFastRT::maxReleave(double& sel, int leave, double maxabs)
{
   UpdateVector& vec = thesolver->fVec();
   Vector& low = thesolver->lbBound();
   Vector& up = thesolver->ubBound();

   //    std::cerr << thesolver->theShift << "\t->\t";

   if (leave >= 0)
   {
      if (up[leave] > low[leave])
      {
         double x = vec.delta()[leave];

         if (sel < -delta / maxabs)
         {
            sel = 0;
            if (x < 0)
               thesolver->shiftLBbound(leave, vec[leave]);
            else
               thesolver->shiftUBbound(leave, vec[leave]);
         }
      }
      else
      {
         sel = 0;
         thesolver->shiftLBbound(leave, vec[leave]);
         thesolver->shiftUBbound(leave, vec[leave]);
      }
   }
   else
      return 1;

   return 0;
}

int SPxFastRT::minReleave(double& sel, int leave, double maxabs)
{
   UpdateVector& vec = thesolver->fVec();
   Vector& low = thesolver->lbBound();
   Vector& up = thesolver->ubBound();

   //    std::cerr << thesolver->theShift << "\t-)\t";

   if (leave >= 0)
   {
      if (up[leave] > low[leave])
      {
         double x = vec.delta()[leave];

         if (sel > delta / maxabs)
         {
            if (x > 0)
            {
               thesolver->theShift += low[leave];
               sel = 0;
               low[leave] = vec[leave] + sel * x;
               thesolver->theShift -= low[leave];
            }
            else
            {
               thesolver->theShift -= up[leave];
               sel = 0;
               up[leave] = vec[leave] + sel * x;
               thesolver->theShift += up[leave];
            }
         }
      }
      else
      {
         sel = 0;
         if (vec[leave] < low[leave])
            thesolver->theShift += low[leave] - vec[leave];
         else
            thesolver->theShift += vec[leave] - up[leave];
         low[leave] = up[leave] = vec[leave];
      }
   }
   else
      return 1;

   return 0;
}

int SPxFastRT::selectLeave(double& val)
{
   double maxabs, max, sel;
   int leave = -1;
   int cnt = 0;

   resetTols();

   if (val > epsilon)
   {
      do
      {
         // phase 1:
         max = val;
         maxabs = 0;
         leave = maxDelta(max, maxabs);
         if (max == val)
            return -1;

         if (!maxShortLeave(sel, leave, max, maxabs))
         {
            // phase 2:
            double stab, bestDelta;
            stab = 100 * minStability(minStab, maxabs);
            // std::cerr << '\t' << stab << '\t' << maxabs << std::endl;
            leave = maxSelect(sel, stab, bestDelta, max);
            if (bestDelta < DELTA_SHIFT*TRIES)
               cnt++;
            else
               cnt += TRIES;
         }
         if (!maxReleave(sel, leave, maxabs))
            break;
         relax();
      }
      while (cnt < TRIES);
   }

   else if (val < -epsilon)
   {
      do
      {
         max = val;
         maxabs = 0;
         leave = minDelta(max, maxabs);
         if (max == val)
            return -1;

         if (!minShortLeave(sel, leave, max, maxabs));
         {
            // phase 2:
            double stab, bestDelta;
            stab = 100 * minStability(minStab, maxabs);
            // std::cerr << '\t' << stab << '\t' << maxabs << std::endl;
            leave = minSelect(sel, stab, bestDelta, max);
            if (bestDelta < DELTA_SHIFT*TRIES)
               cnt++;
            else
               cnt += TRIES;
         }
         if (!minReleave(sel, leave, maxabs))
            break;
         relax();
      }
      while (cnt < TRIES);
   }

   else
      return -1;

#ifndef NDEBUG
   if (leave >= 0)
      fprintf(stderr, "%d(%10.6g,%8.2g):%10d\t%10.4g %10.4g %10.6g\n",
               thesolver->basis().iteration(), thesolver->value(),
               thesolver->basis().stability(),
               leave, sel, thesolver->fVec().delta()[leave], maxabs);
   else
      std::cerr << thesolver->basis().iteration() << ": skipping instable pivot\n";
#endif  // NDEBUG

   if (leave >= 0 || minStab > 2*solver()->epsilon())
   {
      val = sel;
      if (leave >= 0)
         tighten();
   }

   return leave;
}

//@ ----------------------------------------------------------------------------
int SPxFastRT::maxReenter(double& sel, double max, double maxabs,
                           SoPlex::Id id, int nr)
{
   double x, d;
   Vector* up;
   Vector* low;

   UpdateVector& pvec = thesolver->pVec();
   SSVector& pupd = thesolver->pVec().delta();
   Vector& upb = thesolver->upBound();
   Vector& lpb = thesolver->lpBound();
   UpdateVector& cvec = thesolver->coPvec();
   SSVector& cupd = thesolver->coPvec().delta();
   Vector& ucb = thesolver->ucBound();
   Vector& lcb = thesolver->lcBound();

   if (thesolver->isCoId(id))
   {
      if (thesolver->isCoBasic(nr))
      {
         cupd.clearIdx(nr);
         return 1;
      }

      x = cvec[nr];
      d = cupd[nr];
      up = &ucb;
      low = &lcb;

      if (d < 0)
         sel = (lcb[nr] - cvec[nr]) / d;
      else
         sel = (ucb[nr] - cvec[nr]) / d;
   }

   else if (thesolver->isId(id))
   {
      pvec[nr] = thesolver->vector(nr) * cvec;
      if (thesolver->isBasic(nr))
      {
         pupd.clearIdx(nr);
         return 1;
      }

      x = pvec[nr];
      d = pupd[nr];
      up = &upb;
      low = &lpb;

      if (d < 0)
         sel = (lpb[nr] - pvec[nr]) / d;
      else
         sel = (upb[nr] - pvec[nr]) / d;
   }
   else
      return 1;

   if ((*up)[nr] != (*low)[nr])
   {
      if (sel < -delta / maxabs)
      {
         if (d > 0)
         {
            thesolver->theShift -= (*up)[nr];
            sel = 0;
            (*up)[nr] = x + sel * d;
            thesolver->theShift += (*up)[nr];
         }
         else
         {
            thesolver->theShift += (*low)[nr];
            sel = 0;
            (*low)[nr] = x + sel * d;
            thesolver->theShift -= (*low)[nr];
         }
      }
   }
   else
   {
      sel = 0;
      if (x > (*up)[nr])
         thesolver->theShift += x - (*up)[nr];
      else
         thesolver->theShift += (*low)[nr] - x;
      (*up)[nr] = (*low)[nr] = x;
   }

   return 0;
}

int SPxFastRT::minReenter(double& sel, double max, double maxabs,
                           SoPlex::Id id, int nr)
{
   double x, d;
   Vector* up;
   Vector* low;

   UpdateVector& pvec = thesolver->pVec();
   SSVector& pupd = thesolver->pVec().delta();
   Vector& upb = thesolver->upBound();
   Vector& lpb = thesolver->lpBound();
   UpdateVector& cvec = thesolver->coPvec();
   SSVector& cupd = thesolver->coPvec().delta();
   Vector& ucb = thesolver->ucBound();
   Vector& lcb = thesolver->lcBound();

   if (thesolver->isCoId(id))
   {
      if (thesolver->isCoBasic(nr))
      {
         cupd.clearIdx(nr);
         return 1;
      }
      x = cvec[nr];
      d = cupd[nr];
      up = &ucb;
      low = &lcb;
      if (d > 0)
         sel = (thesolver->lcBound()[nr] - cvec[nr]) / d;
      else
         sel = (thesolver->ucBound()[nr] - cvec[nr]) / d;
   }

   else if (thesolver->isId(id))
   {
      pvec[nr] = thesolver->vector(nr) * cvec;
      if (thesolver->isBasic(nr))
      {
         pupd.clearIdx(nr);
         return 1;
      }
      x = pvec[nr];
      d = pupd[nr];
      up = &upb;
      low = &lpb;
      if (d > 0)
         sel = (thesolver->lpBound()[nr] - pvec[nr]) / d;
      else
         sel = (thesolver->upBound()[nr] - pvec[nr]) / d;
   }

   else
      return 1;

   if ((*up)[nr] != (*low)[nr])
   {
      if (sel > delta / maxabs)
      {
         if (d < 0)
         {
            thesolver->theShift -= (*up)[nr];
            sel = 0;
            (*up)[nr] = x + sel * d;
            thesolver->theShift += (*up)[nr];
         }
         else
         {
            thesolver->theShift += (*low)[nr];
            sel = 0;
            (*low)[nr] = x + sel * d;
            thesolver->theShift -= (*low)[nr];
         }
      }
   }
   else
   {
      sel = 0;
      if (x > (*up)[nr])
         thesolver->theShift += x - (*up)[nr];
      else
         thesolver->theShift += (*low)[nr] - x;
      (*up)[nr] = (*low)[nr] = x;
   }

   return 0;
}

int SPxFastRT::shortEnter(
   SoPlex::Id& enterId,
   int nr,
   double max,
   double maxabs
)
{
   if (thesolver->isCoId(enterId))
   {
      if (max != 0)
      {
         double x = thesolver->coPvec().delta()[nr];
         if (x < maxabs * SHORT && -x < maxabs * SHORT)
            return 0;
      }
      return 1;
   }

   else if (thesolver->isId(enterId))
   {
      if (max != 0)
      {
         double x = thesolver->pVec().delta()[nr];
         // std::cerr << x << ' ';
         if (x < maxabs * SHORT && -x < maxabs * SHORT)
            return 0;
      }
      return 1;
   }

   return 0;
}

SoPlex::Id SPxFastRT::selectEnter(double& val)
{
   SoPlex::Id enterId;
   double max, sel;
   double maxabs;
   int nr;
   int cnt = 0;

   resetTols();
   sel = 0;

   if (val > epsilon)
   {
      do
      {
         maxabs = 0;
         max = val;

         enterId = maxDelta(nr, max, maxabs);
         if (!enterId.isValid())
            return enterId;
         assert(max >= 0);

         if (!shortEnter(enterId, nr, max, maxabs))
         {
            double bestDelta, stab;
            // stab = minStab;
            stab = minStability(minStab, maxabs);
            enterId = maxSelect(nr, sel, stab, bestDelta, max);
            if (bestDelta < DELTA_SHIFT*TRIES)
               cnt++;
            else
            {
               // std::cerr << "break";
               cnt += TRIES;
            }
         }
         if (!maxReenter(sel, max, maxabs, enterId, nr))
            break;
         // std::cerr << " ++ ";
         relax();
      }
      while (cnt < TRIES);
   }

   else if (val < -epsilon)
   {
      do
      {
         maxabs = 0;
         max = val;
         enterId = minDelta(nr, max, maxabs);
         if (!enterId.isValid())
            return enterId;
         assert(max <= 0);

         if (!shortEnter(enterId, nr, max, maxabs))
         {
            double bestDelta, stab;
            // stab = minStab;
            stab = minStability(minStab, maxabs);
            enterId = minSelect(nr, sel, stab, bestDelta, max);
            if (bestDelta < DELTA_SHIFT*TRIES)
               cnt++;
            else
               cnt += TRIES;
         }
         if (!minReenter(sel, max, maxabs, enterId, nr))
            break;
         relax();
      }
      while (cnt < TRIES);
   }

#ifndef NDEBUG
   if (enterId.isValid())
   {
      double x;
      if (thesolver->isCoId(enterId))
         x = thesolver->coPvec().delta()[ thesolver->number(enterId) ];
      else
         x = thesolver->pVec().delta()[ thesolver->number(enterId) ];
      std::cerr << thesolver->basis().iteration() << ": " << sel << '\t'
      << x << " (" << maxabs << ")\n";
   }
   else
      std::cerr << thesolver->basis().iteration() << ": skipping instable pivot\n";
#endif  // NDEBUG

   if (enterId.isValid() || minStab > 2*epsilon)
   {
      val = sel;
      if (enterId.isValid())
         tighten();
   }

   return enterId;
}

void SPxFastRT::load(SoPlex* spx)
{
   thesolver = spx;
   setType(spx->type());
}

void SPxFastRT::clear()
{
   thesolver = 0;
}

void SPxFastRT::setType(SoPlex::Type tp)
{
   minStab = MINSTAB;
   delta = thesolver->delta();

   /*
       resetTols();
       std::cerr << "epsilon = " << epsilon << '\t';
       std::cerr << "delta   = " << delta   << '\t';
       std::cerr << "minStab = " << minStab << std::endl;
    */

   if (delta > 1e-4)
      delta = 1e-4;

   delta0 = delta;
}
} // namespace soplex

//-----------------------------------------------------------------------------
//Emacs Local Variables:
//Emacs c-basic-offset:3
//Emacs tab-width:8
//Emacs indent-tabs-mode:nil
//Emacs End:
//-----------------------------------------------------------------------------
