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
#pragma ident "@(#) $Id: vector_c.cpp,v 1.3 2001/11/09 13:25:32 bzfpfend Exp $"

#include "vector_c.h"

namespace soplex
{

/*
#define HARD_UNROLL_LOOPS
*/

/*  \SubSection{Methods implemented in pure C}
    Because most C++ compilers don't generate as efficient code as their C
    counterparts, some kernel methods have been implemented in pure C.
    Corresponding member functions of |Vector| are forwarded to theses methods
    by the inline methods.
 */
void Vector_MultAddSSVector(
   double* vec, double x,
   int n, const int *idx,
   const double *val
)
{
#ifdef  HARD_UNROLL_LOOPS
multaddssvec_label:
   switch (n)
   {
   case 16:
      vec[*idx] += x * val[*idx];
      ++idx;
   case 15:
      vec[*idx] += x * val[*idx];
      ++idx;
   case 14:
      vec[*idx] += x * val[*idx];
      ++idx;
   case 13:
      vec[*idx] += x * val[*idx];
      ++idx;
   case 12:
      vec[*idx] += x * val[*idx];
      ++idx;
   case 11:
      vec[*idx] += x * val[*idx];
      ++idx;
   case 10:
      vec[*idx] += x * val[*idx];
      ++idx;
   case 9:
      vec[*idx] += x * val[*idx];
      ++idx;
   case 8:
      vec[*idx] += x * val[*idx];
      ++idx;
   case 7:
      vec[*idx] += x * val[*idx];
      ++idx;
   case 6:
      vec[*idx] += x * val[*idx];
      ++idx;
   case 5:
      vec[*idx] += x * val[*idx];
      ++idx;
   case 4:
      vec[*idx] += x * val[*idx];
      ++idx;
   case 3:
      vec[*idx] += x * val[*idx];
      ++idx;
   case 2:
      vec[*idx] += x * val[*idx];
      ++idx;
   case 1:
      vec[*idx] += x * val[*idx];
      ++idx;
   case 0:
      break;
   default:
      vec[*idx] += x * val[*idx];
      ++idx;
      vec[*idx] += x * val[*idx];
      ++idx;
      vec[*idx] += x * val[*idx];
      ++idx;
      vec[*idx] += x * val[*idx];
      ++idx;
      vec[*idx] += x * val[*idx];
      ++idx;
      vec[*idx] += x * val[*idx];
      ++idx;
      vec[*idx] += x * val[*idx];
      ++idx;
      vec[*idx] += x * val[*idx];
      ++idx;
      vec[*idx] += x * val[*idx];
      ++idx;
      vec[*idx] += x * val[*idx];
      ++idx;
      vec[*idx] += x * val[*idx];
      ++idx;
      vec[*idx] += x * val[*idx];
      ++idx;
      vec[*idx] += x * val[*idx];
      ++idx;
      vec[*idx] += x * val[*idx];
      ++idx;
      vec[*idx] += x * val[*idx];
      ++idx;
      vec[*idx] += x * val[*idx];
      ++idx;
      n -= 16;
      goto multaddssvec_label;
   }
#else   /* HARD_UNROLL_LOOPS */
   while (n-- > 0)
   {
      vec[*idx] += x * val[*idx];
      ++idx;
   }
#endif  /* HARD_UNROLL_LOOPS */
}

#define HARD_UNROLL_LOOPS
void Vector_MultAddSVector(
   double* vec,
   double x,
   int n,
   void* elem
)
{
#ifdef  HARD_UNROLL_LOOPS
   struct Tmp
   {
      double val;
      int idx;
   };
   struct Tmp* e = static_cast<struct Tmp*>(elem);

#ifndef OLD
   struct Tmp* end;
   end = static_cast<struct Tmp*>(elem) + 16 * static_cast<int>(n / 16);
   while (e < end)
   {
      vec[e->idx] += x * e->val;
      ++e;
      vec[e->idx] += x * e->val;
      ++e;
      vec[e->idx] += x * e->val;
      ++e;
      vec[e->idx] += x * e->val;
      ++e;
      vec[e->idx] += x * e->val;
      ++e;
      vec[e->idx] += x * e->val;
      ++e;
      vec[e->idx] += x * e->val;
      ++e;
      vec[e->idx] += x * e->val;
      ++e;
      vec[e->idx] += x * e->val;
      ++e;
      vec[e->idx] += x * e->val;
      ++e;
      vec[e->idx] += x * e->val;
      ++e;
      vec[e->idx] += x * e->val;
      ++e;
      vec[e->idx] += x * e->val;
      ++e;
      vec[e->idx] += x * e->val;
      ++e;
      vec[e->idx] += x * e->val;
      ++e;
      vec[e->idx] += x * e->val;
      ++e;
   }
   end = static_cast<struct Tmp*>(elem) + 4 * static_cast<int>(n / 4);
   while (e < end)
   {
      vec[e->idx] += x * e->val;
      ++e;
      vec[e->idx] += x * e->val;
      ++e;
      vec[e->idx] += x * e->val;
      ++e;
      vec[e->idx] += x * e->val;
      ++e;
   }
   end = static_cast<struct Tmp*>(elem) + n;
   while (e < end)
   {
      vec[e->idx] += x * e->val;
      ++e;
   }

#else   /* OLD */
   do
   {
      if (n < 8)
      {
         if (n < 4)
         {
            if (n < 2)
            {
               if (n < 1)
                  return;
               goto l1;
            }
            else if (n < 3)
               goto l2;
            else
               goto l3;
         }
         else if (n < 6)
         {
            if (n < 5)
               goto l4;
            else
               goto l5;
         }
         else if (n < 7)
            goto l6;
         else
            goto l7;
      }

      vec[e->idx] += x * e->val;
      e++;
      vec[e->idx] += x * e->val;
      e++;
      vec[e->idx] += x * e->val;
      e++;
      vec[e->idx] += x * e->val;
      e++;
      vec[e->idx] += x * e->val;
      e++;
      vec[e->idx] += x * e->val;
      e++;
      vec[e->idx] += x * e->val;
      e++;
      vec[e->idx] += x * e->val;
      e++;
      n -= 8;
   }
   while (n > 0);
   return;

l8:
   vec[e->idx] += x * e->val;
   e++;
l7:
   vec[e->idx] += x * e->val;
   e++;
l6:
   vec[e->idx] += x * e->val;
   e++;
l5:
   vec[e->idx] += x * e->val;
   e++;
l4:
   vec[e->idx] += x * e->val;
   e++;
l3:
   vec[e->idx] += x * e->val;
   e++;
l2:
   vec[e->idx] += x * e->val;
   e++;
l1:
   vec[e->idx] += x * e->val;
   e++;
#endif  /* OLD */

#else   /* HARD_UNROLL_LOOPS */
   struct Tmp
   {
      double val;
      int idx;
   }
   * const end = elem;
   const struct Tmp* e = end + n;
   while (e-- > end)
      vec[e->idx] += x * e->val;
#endif  /* HARD_UNROLL_LOOPS */
}
#undef HARD_UNROLL_LOOPS

void Vector_MultAddVector(
   double x, int n,
   double* v, const double* w)
{
#ifdef  HARD_UNROLL_LOOPS
multaddvec_label:
   switch (n)
   {
   case 16:
      *v++ += x * (*w++);
   case 15:
      *v++ += x * (*w++);
   case 14:
      *v++ += x * (*w++);
   case 13:
      *v++ += x * (*w++);
   case 12:
      *v++ += x * (*w++);
   case 11:
      *v++ += x * (*w++);
   case 10:
      *v++ += x * (*w++);
   case 9:
      *v++ += x * (*w++);
   case 8:
      *v++ += x * (*w++);
   case 7:
      *v++ += x * (*w++);
   case 6:
      *v++ += x * (*w++);
   case 5:
      *v++ += x * (*w++);
   case 4:
      *v++ += x * (*w++);
   case 3:
      *v++ += x * (*w++);
   case 2:
      *v++ += x * (*w++);
   case 1:
      *v++ += x * (*w++);
   case 0:
      break;
   default:
      *v++ += x * (*w++);
      *v++ += x * (*w++);
      *v++ += x * (*w++);
      *v++ += x * (*w++);
      *v++ += x * (*w++);
      *v++ += x * (*w++);
      *v++ += x * (*w++);
      *v++ += x * (*w++);
      *v++ += x * (*w++);
      *v++ += x * (*w++);
      *v++ += x * (*w++);
      *v++ += x * (*w++);
      *v++ += x * (*w++);
      *v++ += x * (*w++);
      *v++ += x * (*w++);
      *v++ += x * (*w++);
      n -= 16;
      goto multaddvec_label;
   }
#else   /* HARD_UNROLL_LOOPS */
while (n-- > 0)
   *v++ += x * (*w++);
#endif  /* HARD_UNROLL_LOOPS */
}

void Vector_Set0toSSVector(
   double*zero, int n,
   const int* idx, const double* val)
{
   while (n-- > 0)
   {
      zero[*idx] = val[*idx];
      ++idx;
   }
}


double MultiplyVectorSSVector(
   const double* dense, int n,
   const int* idx, const double* val)
{
   double x = 0;
   while (n-- > 0)
   {
      x += dense[*idx] * val[*idx];
      ++idx;
   }
   return x;
}

double MultiplyVectorVector(
   const double* v1, int dim,
   const double* v2)
{
   double x = 0;
   while (dim--)
      x += (*v1++) * (*v2++);
   return x;
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
