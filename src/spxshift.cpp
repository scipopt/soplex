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
#pragma ident "@(#) $Id: spxshift.cpp,v 1.3 2001/11/07 13:48:13 bzfkocht Exp $"

/*      \SubSection{Shifting bounds}
 */

/*  Import system include files
 */
#include <assert.h>
#include <iostream>


/*  and class header files
 */
#include "soplex.h"

namespace soplex
{


//@ ----------------------------------------------------------------------------
void SoPlex::shiftFvec()
{
   // Random   mult(delta(), 100*delta());
   Random mult(10*delta(), 100*delta());
   double allow = delta() - epsilon();

   assert(type() == ENTER);
   assert(allow > 0);
   for (int i = dim() - 1; i >= 0; --i)
   {
      if (theUBbound[i] + allow < (*theFvec)[i])
      {
         if (theUBbound[i] != theLBbound[i])
            shiftUBbound(i, (*theFvec)[i] + double(mult));
         else
         {
            shiftUBbound(i, (*theFvec)[i]);
            theLBbound[i] = theUBbound[i];
         }
      }
      else if ((*theFvec)[i] < theLBbound[i] - allow)
      {
         if (theUBbound[i] != theLBbound[i])
            shiftLBbound(i, (*theFvec)[i] - double(mult));
         else
         {
            shiftLBbound(i, (*theFvec)[i]);
            theUBbound[i] = theLBbound[i];
         }
      }
   }

#ifndef NDEBUG
   testBounds();
   std::cerr << "OK\n\n";
#endif
}

/*
    This methods assumes correctly setup vectors |pVec| and |coPvec| and bound
    vectors for leaving simplex. Then it checks all values of |pVec| and
    |coPvec| to obey these bounds and enlarges them if neccessary.
 */
void SoPlex::shiftPvec()
{
   int i, tmp;
   // Random   mult(delta(), 100*delta());
   Random mult(10*delta(), 100*delta());
   double allow = delta() - epsilon();

   assert(type() == LEAVE);
   for (i = dim() - 1; i >= 0; --i)
   {
      tmp = !isBasic(coId(i));
      if ((*theCoUbound)[i] + allow <= (*theCoPvec)[i] && tmp)
      {
         if ((*theCoUbound)[i] != (*theCoLbound)[i])
            shiftUCbound(i, (*theCoPvec)[i] + double(mult));
         else
         {
            shiftUCbound(i, (*theCoPvec)[i]);
            (*theCoLbound)[i] = (*theCoUbound)[i];
         }
      }
      else if ((*theCoLbound)[i] - allow >= (*theCoPvec)[i] && tmp)
      {
         if ((*theCoUbound)[i] != (*theCoLbound)[i])
            shiftLCbound(i, (*theCoPvec)[i] - double(mult));
         else
         {
            shiftLCbound(i, (*theCoPvec)[i]);
            (*theCoUbound)[i] = (*theCoLbound)[i];
         }
      }
   }

   for (i = coDim() - 1; i >= 0; --i)
   {
      tmp = !isBasic(id(i));
      if ((*theUbound)[i] + allow <= (*thePvec)[i] && tmp)
      {
         if ((*theUbound)[i] != (*theLbound)[i])
            shiftUPbound(i, (*thePvec)[i] + double(mult));
         else
         {
            shiftUPbound(i, (*thePvec)[i]);
            (*theLbound)[i] = (*theUbound)[i];
         }
      }
      else if ((*theLbound)[i] - allow >= (*thePvec)[i] && tmp)
      {
         if ((*theUbound)[i] != (*theLbound)[i])
            shiftLPbound(i, (*thePvec)[i] - double(mult));
         else
         {
            shiftLPbound(i, (*thePvec)[i]);
            (*theUbound)[i] = (*theLbound)[i];
         }
      }
   }

#ifndef NDEBUG
   testBounds();
   std::cerr << "OK\n\n";
#endif
}

void SoPlex::perturbMin
(
   const UpdateVector& uvec,
   Vector& p_low,
   Vector& p_up,
   double eps,
   int start,
   int incr
)
{
   assert(uvec.dim() == p_low.dim());
   assert(uvec.dim() == p_up.dim());

   const double* vec = uvec.get_const_ptr();
   const double* upd = uvec.delta().values();
   const IdxSet& idx = uvec.delta().indices();
   // Random           mult(1*delta(), 100*delta());
   Random mult(10*delta(), 100*delta());
   double x, l, u;
   int i, j;

#ifdef  FULL_SHIFT
   eps = delta();
   for (i = uvec.dim() - start - 1; i >= 0; i -= incr)
   {
      u = p_up[i];
      l = p_low[i];
      if (p_up[i] <= vec[i] + eps)
      {
         p_up[i] = vec[i] + (double)mult;
         theShift += p_up[i] - u;
      }
      if (p_low[i] >= vec[i] - eps)
      {
         p_low[i] = vec[i] - (double)mult;
         theShift -= p_low[i] - l;
      }
   }

#else   // !FULL_SHIFT
for (j = uvec.delta().size() - start - 1; j >= 0; j -= incr)
   {
      i = idx.index(j);
      x = upd[i];
      u = p_up[i];
      l = p_low[i];
      if (x < epsilon())
      {
         if (u != l && vec[i] >= u - eps)
         {
            p_up[i] = vec[i] + (double)mult;
            theShift += p_up[i] - u;
         }
      }
      else if (x > epsilon())
      {
         if (u != l && vec[i] <= l + eps)
         {
            p_low[i] = vec[i] - (double)mult;
            theShift -= p_low[i] - l;
         }
      }
   }
#endif  // !FULL_SHIFT
}

void SoPlex::perturbMax
(
   const UpdateVector& uvec,
   Vector& p_low,
   Vector& p_up,
   double eps,
   int start,
   int incr
)
{
   assert(uvec.dim() == p_low.dim());
   assert(uvec.dim() == p_up.dim());

   const double* vec = uvec.get_const_ptr();
   const double* upd = uvec.delta().values();
   const IdxSet& idx = uvec.delta().indices();
   Random mult(10*delta(), 100*delta());
   double x, l, u;
   int i, j;

#ifdef  FULL_SHIFT
   eps = delta();
   for (i = uvec.dim() - start - 1; i >= 0; i -= incr)
   {
      u = p_up[i];
      l = p_low[i];
      if (p_up[i] <= vec[i] + eps)
      {
         p_up[i] = vec[i] + (double)mult;
         theShift += p_up[i] - u;
      }
      if (p_low[i] >= vec[i] - eps)
      {
         p_low[i] = vec[i] - (double)mult;
         theShift -= p_low[i] - l;
      }
   }

#else   // !FULL_SHIFT
   for (j = uvec.delta().size() - start - 1; j >= 0; j -= incr)
   {
      i = idx.index(j);
      x = upd[i];
      u = p_up[i];
      l = p_low[i];
      if (x > epsilon())
      {
         if (u != l && vec[i] >= u - eps)
         {
            p_up[i] = vec[i] + (double)mult;
            theShift += p_up[i] - u;
         }
      }
      else if (x < epsilon())
      {
         if (u != l && vec[i] <= l + eps)
         {
            p_low[i] = vec[i] - (double)mult;
            theShift -= p_low[i] - l;
         }
      }
   }
#endif  // !FULL_SHIFT
}

void SoPlex::perturbMinEnter(void)
{
   //@ std::cerr << iteration() << ": perturbing " << shift();
   fVec().delta().setup();
   perturbMin(fVec(), lbBound(), ubBound(), epsilon());
   //@ std::cerr << "\t->" << shift() << std::endl;
}


void SoPlex::perturbMaxEnter(void)
{
   //@ std::cerr << iteration() << ": perturbing " << shift();
   fVec().delta().setup();
   perturbMax(fVec(), lbBound(), ubBound(), epsilon());
   //@ std::cerr << "\t->" << shift() << std::endl;
}


double SoPlex::perturbMin
(
   const UpdateVector& uvec,
   Vector& p_low,
   Vector& p_up,
   double eps,
   double delta,
   const SPxBasis::Desc::Status* stat,
   int start,
   int incr
)
{
   assert(uvec.dim() == p_low.dim());
   assert(uvec.dim() == p_up.dim());

   const double* vec = uvec.get_const_ptr();
   const double* upd = uvec.delta().values();
   const IdxSet& idx = uvec.delta().indices();
   Random mult(10*delta, 100*delta);
   double x, l, u;
   int i, j;
   double theShift = 0;

#ifdef  FULL_SHIFT
   eps = delta;
   for (i = uvec.dim() - start - 1; i >= 0; i -= incr)
   {
      u = p_up[i];
      l = p_low[i];
      if (p_up[i] <= vec[i] + eps && rep()*stat[i] < 0)
      {
         p_up[i] = vec[i] + (double)mult;
         theShift += p_up[i] - u;
      }
      if (p_low[i] >= vec[i] - eps && rep()*stat[i] < 0)
      {
         p_low[i] = vec[i] - (double)mult;
         theShift -= p_low[i] - l;
      }
   }

#else   // !FULL_SHIFT
   for (j = uvec.delta().size() - start - 1; j >= 0; j -= incr)
   {
      i = idx.index(j);
      x = upd[i];
      u = p_up[i];
      l = p_low[i];
      if (x < eps)
      {
         if (u != l && vec[i] >= u - eps && rep()*stat[i] < 0)
         {
            p_up[i] = vec[i] + (double)mult;
            theShift += p_up[i] - u;
         }
      }
      else if (x > eps)
      {
         if (u != l && vec[i] <= l + eps && rep()*stat[i] < 0)
         {
            p_low[i] = vec[i] - (double)mult;
            theShift -= p_low[i] - l;
         }
      }
   }
#endif  // !FULL_SHIFT 
   return theShift;
}

double SoPlex::perturbMax
(
   const UpdateVector& uvec,
   Vector& p_low,
   Vector& p_up,
   double eps,
   double delta,
   const SPxBasis::Desc::Status* stat,
   int start,
   int incr
)
{
   assert(uvec.dim() == p_low.dim());
   assert(uvec.dim() == p_up.dim());

   const double* vec = uvec.get_const_ptr();
   const double* upd = uvec.delta().values();
   const IdxSet& idx = uvec.delta().indices();
   Random mult(10*delta, 100*delta);
   double x, l, u;
   int i, j;
   double theShift = 0;

#ifdef  FULL_SHIFT
   eps = delta;
   for (i = uvec.dim() - start - 1; i >= 0; i -= incr)
   {
      u = p_up[i];
      l = p_low[i];
      if (p_up[i] <= vec[i] + eps && rep()*stat[i] < 0)
      {
         p_up[i] = vec[i] + (double)mult;
         theShift += p_up[i] - u;
      }
      if (p_low[i] >= vec[i] - eps && rep()*stat[i] < 0)
      {
         p_low[i] = vec[i] - (double)mult;
         theShift -= p_low[i] - l;
      }
   }

#else   // !FULL_SHIFT
   for (j = uvec.delta().size() - start - 1; j >= 0; j -= incr)
   {
      i = idx.index(j);
      x = upd[i];
      u = p_up[i];
      l = p_low[i];
      if (x > eps)
      {
         if (u != l && vec[i] >= u - eps && rep()*stat[i] < 0)
         {
            p_up[i] = vec[i] + (double)mult;
            theShift += p_up[i] - u;
         }
      }
      else if (x < eps)
      {
         if (u != l && vec[i] <= l + eps && rep()*stat[i] < 0)
         {
            p_low[i] = vec[i] - (double)mult;
            theShift -= p_low[i] - l;
         }
      }
   }
#endif  // !FULL_SHIFT 
   return theShift;
}


void SoPlex::perturbMinLeave(void)
{
   //@ std::cerr << iteration() << ": perturbing " << shift();
   pVec().delta().setup();
   coPvec().delta().setup();
   theShift += perturbMin(pVec(), lpBound(), upBound(), epsilon(), delta(),
                           desc().status(), 0, 1);
   theShift += perturbMin(coPvec(), lcBound(), ucBound(), epsilon(), delta(),
                           desc().coStatus(), 0, 1);
   //@ std::cerr << "\t->" << shift() << std::endl;
}


void SoPlex::perturbMaxLeave(void)
{
   //@ std::cerr << iteration() << ": perturbing " << shift();
   pVec().delta().setup();
   coPvec().delta().setup();
   theShift += perturbMax(pVec(), lpBound(), upBound(), epsilon(), delta(),
                           desc().status(), 0, 1);
   theShift += perturbMax(coPvec(), lcBound(), ucBound(), epsilon(), delta(),
                           desc().coStatus(), 0, 1);
   //@ std::cerr << "\t->" << shift() << std::endl;
}


void SoPlex::unShift(void)
{
   //@ std::cerr << "unshifting ...\n";
   if (isInitialized())
   {
      int i;
      double t_up, t_low;
      double eps = delta();
      const SPxBasis::Desc& ds = desc();

      theShift = 0;
      if (type() == ENTER)
      {
         if (rep() == COLUMN)
         {
            for (i = dim(); i-- > 0;)
            {
               Id id = baseId(i);
               int num = number(id);
               if (id.type() == Id::ROWID)
               {
                  t_up = -lhs(num);
                  t_low = -rhs(num);
               }
               else
               {
                  assert(id.type() == Id::COLID);
                  t_up = upper(num);
                  t_low = lower(num);
               }
               if (t_up != t_low)
               {
                  if ((*theFvec)[i] < t_up - eps)
                     theUBbound[i] = t_up;
                  else if ((*theFvec)[i] > t_up)
                     theShift += theUBbound[i] - t_up;
                  if ((*theFvec)[i] > t_low + eps)
                     theLBbound[i] = t_low;
                  else if ((*theFvec)[i] < t_low)
                     theShift -= theLBbound[i] - t_low;
               }
               else
               {
                  if (theUBbound[i] > t_up)
                     theShift += theUBbound[i] - t_up;
                  else if (theLBbound[i] < t_low)
                     theShift += t_low - theLBbound[i];
               }
            }
            for (i = nRows(); i-- > 0;)
            {
               if (!isBasic(ds.rowStatus(i)))
               {
                  t_up = -lhs(i);
                  t_low = -rhs(i);
                  if (theURbound[i] > t_up)
                     theShift += theURbound[i] - t_up;
                  if (t_low > theLRbound[i])
                     theShift += t_low - theLRbound[i];
               }
            }
            for (i = nCols(); i-- > 0;)
            {
               if (!isBasic(ds.colStatus(i)))
               {
                  t_up = upper(i);
                  t_low = lower(i);
                  if (theUCbound[i] > t_up)
                     theShift += theUCbound[i] - t_up;
                  if (t_low > theLCbound[i])
                     theShift += t_low - theLCbound[i];
               }
            }
         }
         else
         {
            assert(rep() == ROW);
            for (i = dim(); i-- > 0;)
            {
               Id id = baseId(i);
               int num = number(id);
               t_up = t_low = 0;
               if (id.type() == Id::ROWID)
                  clearDualBounds(ds.rowStatus(num), t_up, t_low);
               else
                  clearDualBounds(ds.colStatus(num), t_up, t_low);
               if (theUBbound[i] != theLBbound[i])
               {
                  if (theUBbound[i] > t_up)
                     theShift += theUBbound[i] - t_up;
                  else
                     theShift -= theUBbound[i] - t_up;
               }
               else
               {
                  assert(theUBbound[i] >= t_up);
                  if ((*theFvec)[i] < t_up - eps)
                     theUBbound[i] = t_up;
                  else if ((*theFvec)[i] > t_up)
                     theShift += theUBbound[i] - t_up;
                  assert(theLBbound[i] <= t_low);
                  if ((*theFvec)[i] > t_low + eps)
                     theLBbound[i] = t_low;
                  else if ((*theFvec)[i] < t_low)
                     theShift -= theLBbound[i] - t_low;
               }
            }
            for (i = nRows(); i-- > 0;)
            {
               if (!isBasic(ds.rowStatus(i)))
               {
                  t_up = t_low = 0;
                  clearDualBounds(ds.rowStatus(i), t_up, t_low);
                  if (theURbound[i] > t_up)
                     theShift += theURbound[i] - t_up;
                  if (t_low > theLRbound[i])
                     theShift += t_low - theLRbound[i];
               }
            }
            for (i = nCols(); i-- > 0;)
            {
               if (!isBasic(ds.colStatus(i)))
               {
                  t_up = t_low = 0;
                  clearDualBounds(ds.colStatus(i), t_up, t_low);
                  if (theUCbound[i] > t_up)
                     theShift += theUCbound[i] - t_up;
                  if (t_low > theLCbound[i])
                     theShift += t_low - theLCbound[i];
               }
            }
         }
      }
      else
      {
         assert(type() == LEAVE);
         if (rep() == COLUMN)
         {
            for (i = nRows(); i-- > 0;)
            {
               t_up = t_low = 0;
               clearDualBounds(ds.rowStatus(i), t_up, t_low);
               if (!isBasic(ds.rowStatus(i)))
               {
                  if ((*theCoPvec)[i] < t_up - eps)
                     theURbound[i] = t_up;
                  else
                     theShift += theURbound[i] - t_up;
                  if ((*theCoPvec)[i] > t_low + eps)
                     theLRbound[i] = t_low;
                  else
                     theShift += t_low - theLRbound[i];
               }
               else if (theURbound[i] > t_up)
                  theShift += theURbound[i] - t_up;
               else if (theLRbound[i] < t_low)
                  theShift += t_low - theLRbound[i];
            }
            for (i = nCols(); i-- > 0;)
            {
               t_up = t_low = -maxObj(i);
               clearDualBounds(ds.colStatus(i), t_low, t_up);
               if (!isBasic(ds.colStatus(i)))
               {
                  if ((*thePvec)[i] < -t_up - eps)
                     theUCbound[i] = -t_up;
                  else
                     theShift += theUCbound[i] - (-t_up);
                  if ((*thePvec)[i] > -t_low + eps)
                     theLCbound[i] = -t_low;
                  else
                     theShift += (-t_low) - theLCbound[i];
               }
               else if (theUCbound[i] > -t_up)
                  theShift += theUCbound[i] - (-t_up);
               else if (theLCbound[i] < -t_low)
                  theShift += (-t_low) - theLCbound[i];
            }
         }
         else
         {
            assert(rep() == ROW);
            for (i = nRows(); i-- > 0;)
            {
               t_up = rhs(i);
               t_low = lhs(i);
               if (t_up == t_low)
               {
                  if (theURbound[i] > t_up)
                     theShift += theURbound[i] - t_up;
                  else
                     theShift += t_low - theLRbound[i];
               }
               else
                  if (!isBasic(ds.rowStatus(i)))
                  {
                     if ((*thePvec)[i] < t_up - eps)
                        theURbound[i] = t_up;
                     else
                        theShift += theURbound[i] - t_up;
                     if ((*thePvec)[i] > t_low + eps)
                        theLRbound[i] = t_low;
                     else
                        theShift += t_low - theLRbound[i];
                  }
                  else if (theURbound[i] > t_up)
                     theShift += theURbound[i] - t_up;
                  else if (theLRbound[i] < t_low)
                     theShift += t_low - theLRbound[i];
            }
            for (i = nCols(); i-- > 0;)
            {
               t_up = upper(i);
               t_low = lower(i);
               if (t_up == t_low)
               {
                  if (theUCbound[i] > t_up)
                     theShift += theUCbound[i] - t_up;
                  else
                     theShift += t_low - theLCbound[i];
               }
               else
                  if (!isBasic(ds.colStatus(i)))
                  {
                     if ((*theCoPvec)[i] < t_up - eps)
                        theUCbound[i] = t_up;
                     else
                        theShift += theUCbound[i] - t_up;
                     if ((*theCoPvec)[i] > t_low + eps)
                        theLCbound[i] = t_low;
                     else
                        theShift += t_low - theLCbound[i];
                  }
                  else if (theUCbound[i] > t_up)
                     theShift += theUCbound[i] - t_up;
                  else if (theLCbound[i] < t_low)
                     theShift += t_low - theLCbound[i];
            }
         }
      }
   }
}
} // namespace soplex

//-----------------------------------------------------------------------------
//Emacs Local Variables:
//Emacs mode:c++
//Emacs c-basic-offset:3
//Emacs tab-width:8
//Emacs indent-tabs-mode:nil
//Emacs End:
//-----------------------------------------------------------------------------
