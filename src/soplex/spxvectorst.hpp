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
#include "soplex/spxvectorst.h"

namespace soplex
{

template <class R>
void SPxVectorST<R>::setupWeights(SPxSolverBase<R>& base)
{
   if(state == PVEC)
   {
      if(vec.dim() != base.nCols())
      {
         SPxWeightST<R>::setupWeights(base);
         return;
      }

      const VectorBase<R>& obj = base.maxObj();
      R eps = base.epsilon();
      R bias = 10000 * eps;
      R x, y;
      int i;

      MSG_DEBUG(std::cout << "DVECST01 colWeight[]: ";)

      for(i = base.nCols(); i--;)
      {
         x = vec[i] - base.SPxLPBase<R>::lower(i);
         y = base.SPxLPBase<R>::upper(i) - vec[i];

         if(x < y)
         {
            this->colWeight[i] = -x - bias * obj[i];
            this->colUp[i] = 0;
         }
         else
         {
            this->colWeight[i] = -y + bias * obj[i];
            this->colUp[i] = 1;
         }

         MSG_DEBUG(std::cout << this->colWeight[i] << " ";)
      }

      MSG_DEBUG(std::cout << std::endl << std::endl;)

      MSG_DEBUG(std::cout << "DVECST02 rowWeight[]: ";)

      for(i = base.nRows(); i--;)
      {
         const SVectorBase<R>& row = base.rowVector(i);
         y = vec * row;
         x = (y - base.lhs(i));
         y = (base.rhs(i) - y);

         if(x < y)
         {
            this->rowWeight[i] = -x - eps * row.size() - bias * (obj * row);
            this->rowRight[i] = 0;
         }
         else
         {
            this->rowWeight[i] = -y - eps * row.size() + bias * (obj * row);
            this->rowRight[i] = 1;
         }

         MSG_DEBUG(std::cout << this->rowWeight[i] << " ";)
      }

      MSG_DEBUG(std::cout << std::endl;)
   }

   else if(state == DVEC)
   {
      if(vec.dim() != base.nRows())
      {
         SPxWeightST<R>::setupWeights(base);
         return;
      }

      R x, y, len;
      int i, j;

      for(i = base.nRows(); i--;)
         this->rowWeight[i] += spxAbs(vec[i]);

      for(i = base.nCols(); i--;)
      {
         const SVectorBase<R>& col = base.colVector(i);

         for(y = len = 0, j = col.size(); j--;)
         {
            x = col.value(j);
            y += vec[col.index(j)] * x;
            len += x * x;
         }

         if(len > 0)
            this->colWeight[i] += spxAbs(y / len - base.maxObj(i));
      }
   }
   else
      SPxWeightST<R>::setupWeights(base);
}
} // namespace soplex
