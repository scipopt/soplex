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
#pragma ident "@(#) $Id: spxshift.cpp,v 1.1 2001/11/06 16:18:32 bzfkocht Exp $"

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
   Vector& low,
   Vector& up,
   double eps,
   int start,
   int incr
)
{
   assert(uvec.dim() == low.dim());
   assert(uvec.dim() == up.dim());

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
      u = up[i];
      l = low[i];
      if (up[i] <= vec[i] + eps)
      {
         up[i] = vec[i] + (double)mult;
         theShift += up[i] - u;
      }
      if (low[i] >= vec[i] - eps)
      {
         low[i] = vec[i] - (double)mult;
         theShift -= low[i] - l;
      }
   }

#else   // !FULL_SHIFT
for (j = uvec.delta().size() - start - 1; j >= 0; j -= incr)
   {
      i = idx.index(j);
      x = upd[i];
      u = up[i];
      l = low[i];
      if (x < epsilon())
      {
         if (u != l && vec[i] >= u - eps)
         {
            up[i] = vec[i] + (double)mult;
            theShift += up[i] - u;
         }
      }
      else if (x > epsilon())
      {
         if (u != l && vec[i] <= l + eps)
         {
            low[i] = vec[i] - (double)mult;
            theShift -= low[i] - l;
         }
      }
   }
#endif  // !FULL_SHIFT
}

void SoPlex::perturbMax
(
   const UpdateVector& uvec,
   Vector& low,
   Vector& up,
   double eps,
   int start,
   int incr
)
{
   assert(uvec.dim() == low.dim());
   assert(uvec.dim() == up.dim());

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
      u = up[i];
      l = low[i];
      if (up[i] <= vec[i] + eps)
      {
         up[i] = vec[i] + (double)mult;
         theShift += up[i] - u;
      }
      if (low[i] >= vec[i] - eps)
      {
         low[i] = vec[i] - (double)mult;
         theShift -= low[i] - l;
      }
   }

#else   // !FULL_SHIFT
   for (j = uvec.delta().size() - start - 1; j >= 0; j -= incr)
   {
      i = idx.index(j);
      x = upd[i];
      u = up[i];
      l = low[i];
      if (x > epsilon())
      {
         if (u != l && vec[i] >= u - eps)
         {
            up[i] = vec[i] + (double)mult;
            theShift += up[i] - u;
         }
      }
      else if (x < epsilon())
      {
         if (u != l && vec[i] <= l + eps)
         {
            low[i] = vec[i] - (double)mult;
            theShift -= low[i] - l;
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
   Vector& low,
   Vector& up,
   double eps,
   double delta,
   const SPxBasis::Desc::Status* stat,
   int start,
   int incr
)
{
   assert(uvec.dim() == low.dim());
   assert(uvec.dim() == up.dim());

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
      u = up[i];
      l = low[i];
      if (up[i] <= vec[i] + eps && rep()*stat[i] < 0)
      {
         up[i] = vec[i] + (double)mult;
         theShift += up[i] - u;
      }
      if (low[i] >= vec[i] - eps && rep()*stat[i] < 0)
      {
         low[i] = vec[i] - (double)mult;
         theShift -= low[i] - l;
      }
   }

#else   // !FULL_SHIFT
   for (j = uvec.delta().size() - start - 1; j >= 0; j -= incr)
   {
      i = idx.index(j);
      x = upd[i];
      u = up[i];
      l = low[i];
      if (x < eps)
      {
         if (u != l && vec[i] >= u - eps && rep()*stat[i] < 0)
         {
            up[i] = vec[i] + (double)mult;
            theShift += up[i] - u;
         }
      }
      else if (x > eps)
      {
         if (u != l && vec[i] <= l + eps && rep()*stat[i] < 0)
         {
            low[i] = vec[i] - (double)mult;
            theShift -= low[i] - l;
         }
      }
   }
#endif  // !FULL_SHIFT 
   return theShift;
}

double SoPlex::perturbMax
(
   const UpdateVector& uvec,
   Vector& low,
   Vector& up,
   double eps,
   double delta,
   const SPxBasis::Desc::Status* stat,
   int start,
   int incr
)
{
   assert(uvec.dim() == low.dim());
   assert(uvec.dim() == up.dim());

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
      u = up[i];
      l = low[i];
      if (up[i] <= vec[i] + eps && rep()*stat[i] < 0)
      {
         up[i] = vec[i] + (double)mult;
         theShift += up[i] - u;
      }
      if (low[i] >= vec[i] - eps && rep()*stat[i] < 0)
      {
         low[i] = vec[i] - (double)mult;
         theShift -= low[i] - l;
      }
   }

#else   // !FULL_SHIFT
   for (j = uvec.delta().size() - start - 1; j >= 0; j -= incr)
   {
      i = idx.index(j);
      x = upd[i];
      u = up[i];
      l = low[i];
      if (x > eps)
      {
         if (u != l && vec[i] >= u - eps && rep()*stat[i] < 0)
         {
            up[i] = vec[i] + (double)mult;
            theShift += up[i] - u;
         }
      }
      else if (x < eps)
      {
         if (u != l && vec[i] <= l + eps && rep()*stat[i] < 0)
         {
            low[i] = vec[i] - (double)mult;
            theShift -= low[i] - l;
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
      double up, low;
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
                  up = -lhs(num);
                  low = -rhs(num);
               }
               else
               {
                  assert(id.type() == Id::COLID);
                  up = upper(num);
                  low = lower(num);
               }
               if (up != low)
               {
                  if ((*theFvec)[i] < up - eps)
                     theUBbound[i] = up;
                  else if ((*theFvec)[i] > up)
                     theShift += theUBbound[i] - up;
                  if ((*theFvec)[i] > low + eps)
                     theLBbound[i] = low;
                  else if ((*theFvec)[i] < low)
                     theShift -= theLBbound[i] - low;
               }
               else
               {
                  if (theUBbound[i] > up)
                     theShift += theUBbound[i] - up;
                  else if (theLBbound[i] < low)
                     theShift += low - theLBbound[i];
               }
            }
            for (i = nRows(); i-- > 0;)
            {
               if (!isBasic(ds.rowStatus(i)))
               {
                  up = -lhs(i);
                  low = -rhs(i);
                  if (theURbound[i] > up)
                     theShift += theURbound[i] - up;
                  if (low > theLRbound[i])
                     theShift += low - theLRbound[i];
               }
            }
            for (i = nCols(); i-- > 0;)
            {
               if (!isBasic(ds.colStatus(i)))
               {
                  up = upper(i);
                  low = lower(i);
                  if (theUCbound[i] > up)
                     theShift += theUCbound[i] - up;
                  if (low > theLCbound[i])
                     theShift += low - theLCbound[i];
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
               up = low = 0;
               if (id.type() == Id::ROWID)
                  clearDualBounds(ds.rowStatus(num), up, low);
               else
                  clearDualBounds(ds.colStatus(num), up, low);
               if (theUBbound[i] != theLBbound[i])
               {
                  if (theUBbound[i] > up)
                     theShift += theUBbound[i] - up;
                  else
                     theShift -= theUBbound[i] - up;
               }
               else
               {
                  assert(theUBbound[i] >= up);
                  if ((*theFvec)[i] < up - eps)
                     theUBbound[i] = up;
                  else if ((*theFvec)[i] > up)
                     theShift += theUBbound[i] - up;
                  assert(theLBbound[i] <= low);
                  if ((*theFvec)[i] > low + eps)
                     theLBbound[i] = low;
                  else if ((*theFvec)[i] < low)
                     theShift -= theLBbound[i] - low;
               }
            }
            for (i = nRows(); i-- > 0;)
            {
               if (!isBasic(ds.rowStatus(i)))
               {
                  up = low = 0;
                  clearDualBounds(ds.rowStatus(i), up, low);
                  if (theURbound[i] > up)
                     theShift += theURbound[i] - up;
                  if (low > theLRbound[i])
                     theShift += low - theLRbound[i];
               }
            }
            for (i = nCols(); i-- > 0;)
            {
               if (!isBasic(ds.colStatus(i)))
               {
                  up = low = 0;
                  clearDualBounds(ds.colStatus(i), up, low);
                  if (theUCbound[i] > up)
                     theShift += theUCbound[i] - up;
                  if (low > theLCbound[i])
                     theShift += low - theLCbound[i];
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
               up = low = 0;
               clearDualBounds(ds.rowStatus(i), up, low);
               if (!isBasic(ds.rowStatus(i)))
               {
                  if ((*theCoPvec)[i] < up - eps)
                     theURbound[i] = up;
                  else
                     theShift += theURbound[i] - up;
                  if ((*theCoPvec)[i] > low + eps)
                     theLRbound[i] = low;
                  else
                     theShift += low - theLRbound[i];
               }
               else if (theURbound[i] > up)
                  theShift += theURbound[i] - up;
               else if (theLRbound[i] < low)
                  theShift += low - theLRbound[i];
            }
            for (i = nCols(); i-- > 0;)
            {
               up = low = -maxObj(i);
               clearDualBounds(ds.colStatus(i), low, up);
               if (!isBasic(ds.colStatus(i)))
               {
                  if ((*thePvec)[i] < -up - eps)
                     theUCbound[i] = -up;
                  else
                     theShift += theUCbound[i] - (-up);
                  if ((*thePvec)[i] > -low + eps)
                     theLCbound[i] = -low;
                  else
                     theShift += (-low) - theLCbound[i];
               }
               else if (theUCbound[i] > -up)
                  theShift += theUCbound[i] - (-up);
               else if (theLCbound[i] < -low)
                  theShift += (-low) - theLCbound[i];
            }
         }
         else
         {
            assert(rep() == ROW);
            for (i = nRows(); i-- > 0;)
            {
               up = rhs(i);
               low = lhs(i);
               if (up == low)
               {
                  if (theURbound[i] > up)
                     theShift += theURbound[i] - up;
                  else
                     theShift += low - theLRbound[i];
               }
               else
                  if (!isBasic(ds.rowStatus(i)))
                  {
                     if ((*thePvec)[i] < up - eps)
                        theURbound[i] = up;
                     else
                        theShift += theURbound[i] - up;
                     if ((*thePvec)[i] > low + eps)
                        theLRbound[i] = low;
                     else
                        theShift += low - theLRbound[i];
                  }
                  else if (theURbound[i] > up)
                     theShift += theURbound[i] - up;
                  else if (theLRbound[i] < low)
                     theShift += low - theLRbound[i];
            }
            for (i = nCols(); i-- > 0;)
            {
               up = upper(i);
               low = lower(i);
               if (up == low)
               {
                  if (theUCbound[i] > up)
                     theShift += theUCbound[i] - up;
                  else
                     theShift += low - theLCbound[i];
               }
               else
                  if (!isBasic(ds.colStatus(i)))
                  {
                     if ((*theCoPvec)[i] < up - eps)
                        theUCbound[i] = up;
                     else
                        theShift += theUCbound[i] - up;
                     if ((*theCoPvec)[i] > low + eps)
                        theLCbound[i] = low;
                     else
                        theShift += low - theLCbound[i];
                  }
                  else if (theUCbound[i] > up)
                     theShift += theUCbound[i] - up;
                  else if (theLCbound[i] < low)
                     theShift += low - theLCbound[i];
            }
         }
      }
   }
}
} // namespace soplex

//-----------------------------------------------------------------------------
//Emacs Local Variables:
//Emacs c-basic-offset:3
//Emacs tab-width:8
//Emacs indent-tabs-mode:nil
//Emacs End:
//-----------------------------------------------------------------------------
