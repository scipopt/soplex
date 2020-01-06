/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the class library                   */
/*       SoPlex --- the Sequential object-oriented simPlex.                  */
/*                                                                           */
/*    Copyright (C) 1996-2020 Konrad-Zuse-Zentrum                            */
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
#include "soplex/spxsolver.h"

namespace soplex
{

template <class R>
void SPxSolverBase<R>::qualConstraintViolation(R& maxviol, R& sumviol) const
{
   maxviol = 0.0;
   sumviol = 0.0;

   VectorBase<R> solu(this->nCols());

   getPrimalSol(solu);

   for(int row = 0; row < this->nRows(); ++row)
   {
      const SVectorBase<R>& rowvec = this->rowVector(row);

      R val = 0.0;

      for(int col = 0; col < rowvec.size(); ++col)
         val += rowvec.value(col) * solu[rowvec.index(col)];

      R viol = 0.0;

      assert(this->lhs(row) <= this->rhs(row) + 1e-9);

      if(val < this->lhs(row))
         viol = spxAbs(val - this->lhs(row));
      else if(val > this->rhs(row))
         viol = spxAbs(val - this->rhs(row));

      if(viol > maxviol)
         maxviol = viol;

      sumviol += viol;
   }
}

template <class R>
void SPxSolverBase<R>::qualBoundViolation(
   R& maxviol, R& sumviol) const
{
   maxviol = 0.0;
   sumviol = 0.0;

   VectorBase<R> solu(this->nCols());

   getPrimalSol(solu);

   for(int col = 0; col < this->nCols(); ++col)
   {
      assert(this->lower(col) <= this->upper(col) + 1e-9);

      R viol = 0.0;

      if(solu[col] < this->lower(col))
         viol = spxAbs(solu[col] - this->lower(col));
      else if(solu[col] > this->upper(col))
         viol = spxAbs(solu[col] - this->upper(col));

      if(viol > maxviol)
         maxviol = viol;

      sumviol += viol;
   }
}

template <class R>
void SPxSolverBase<R>::qualSlackViolation(R& maxviol, R& sumviol) const
{
   maxviol = 0.0;
   sumviol = 0.0;

   VectorBase<R> solu(this->nCols());
   VectorBase<R> slacks(this->nRows());

   getPrimalSol(solu);
   getSlacks(slacks);

   for(int row = 0; row < this->nRows(); ++row)
   {
      const SVectorBase<R>& rowvec = this->rowVector(row);

      R val = 0.0;

      for(int col = 0; col < rowvec.size(); ++col)
         val += rowvec.value(col) * solu[rowvec.index(col)];

      R viol = spxAbs(val - slacks[row]);

      if(viol > maxviol)
         maxviol = viol;

      sumviol += viol;
   }
}

template <class R>
void SPxSolverBase<R>::qualRedCostViolation(R& maxviol, R& sumviol) const
{
   maxviol = 0.0;
   sumviol = 0.0;

   int i;

   // TODO:   y = c_B * B^-1  => coSolve(y, c_B)
   //         redcost = c_N - yA_N
   // solve system "x = e_i^T * B^-1" to get i'th row of B^-1
   // VectorBase<R> y( this->nRows() );
   // basis().coSolve( x, spx->unitVector( i ) );
   // VectorBase<R> rdcost( this->nCols() );
   if(type() == ENTER)
   {
      for(i = 0; i < dim(); ++i)
      {
         R x = coTest()[i];

         if(x < 0.0)
         {
            sumviol -= x;

            if(x < maxviol)
               maxviol = x;
         }
      }

      for(i = 0; i < coDim(); ++i)
      {
         R x = test()[i];

         if(x < 0.0)
         {
            sumviol -= x;

            if(x < maxviol)
               maxviol = x;
         }
      }
   }
   else
   {
      assert(type() == LEAVE);

      for(i = 0; i < dim(); ++i)
      {
         R x = fTest()[i];

         if(x < 0.0)
         {
            sumviol -= x;

            if(x < maxviol)
               maxviol = x;
         }
      }
   }

   maxviol *= -1;
}

} // namespace soplex
