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

#include <iostream>
#include <assert.h>

#include "soplex.h"
#include "soplex/statistics.h"

namespace soplex
{
/// check scaling of LP
template <class R>
void SoPlexBase<R>::_checkScaling(SPxLPBase<R>* origLP) const
{
   SPX_MSG_INFO1(spxout, spxout << "DEBUG: checking correctness of scaled LP" << std::endl;)
   assert(_realLP->nCols() == origLP->nCols());
   assert(_realLP->nRows() == origLP->nRows());
   assert(_realLP->isScaled());
   assert(!origLP->isScaled());
   bool correct = true;

   SPX_MSG_INFO1(spxout, spxout << "DEBUG: checking rows..." << std::endl;)

   for(int i = 0; i < origLP->nRows(); ++i)
   {
      assert(EQ(origLP->lhs(i), _realLP->lhsUnscaled(i), _solver.epsilon()));
      assert(EQ(origLP->rhs(i), _realLP->rhsUnscaled(i), _solver.epsilon()));

      DSVectorBase<R> row;
      _realLP->getRowVectorUnscaled(i, row);

      assert(origLP->rowVector(i).size() == row.size());

      for(int j = 0; j < row.size(); ++j)
      {
         if(NE(row.value(j), origLP->rowVector(i).value(j), _solver.tolerances()->floatingPointFeastol()))
         {
            SPX_MSG_INFO1(spxout, spxout << "DEBUG: scaling error in row " << i << ", col " << j
                          << ": orig " << origLP->rowVector(i).value(j)
                          << ", unscaled: " << row.value(j) << std::endl;)
            correct = false;
         }
      }
   }

   SPX_MSG_INFO1(spxout, spxout << "DEBUG: checking cols..." << std::endl;)

   for(int i = 0; i < origLP->nCols(); ++i)
   {
      assert(EQ(origLP->lower(i), _realLP->lowerUnscaled(i), _solver.epsilon()));
      assert(EQ(origLP->upper(i), _realLP->upperUnscaled(i), _solver.epsilon()));
      assert(EQ(origLP->obj(i), _realLP->objUnscaled(i), _solver.epsilon()));

      DSVectorBase<R> col;
      _realLP->getColVectorUnscaled(i, col);

      assert(origLP->colVector(i).size() == col.size());

      for(int j = 0; j < col.size(); ++j)
      {
         if(NE(col.value(j), origLP->colVector(i).value(j), _solver.tolerances()->floatingPointFeastol()))
         {
            SPX_MSG_INFO1(spxout, spxout << "DEBUG: scaling error in col " << i << ", row " << j
                          << ": orig " << origLP->colVector(i).value(j)
                          << ", unscaled: " << col.value(j) << std::endl;)
            correct = false;
         }
      }
   }

   if(!correct)
   {
      SPX_MSG_INFO1(spxout, spxout << "DEBUG: scaling check failed" << std::endl;)
   }

   assert(correct);
}


template <class R>
void SoPlexBase<R>::_checkBasisScaling()
{
   if(_status != SPxSolverBase<R>::OPTIMAL)
   {
      SPX_MSG_INFO1(spxout, spxout << "DEBUG: skipping test on non optimal bases\n");
      return;
   }

   assert(&_solver == _realLP);
   VectorBase<R>** binvcol = 0;
   VectorBase<R>** binvrow = 0;
   int* inds = nullptr;
   int basisdim = _solver.rep() == SPxSolverBase<R>::COLUMN ? _solver.nRows() : _solver.nCols();
   spx_alloc(binvcol, basisdim);
   spx_alloc(binvrow, basisdim);
   spx_alloc(inds, basisdim);

   SPX_MSG_INFO1(spxout, spxout << "DEBUG: computing columns of inverted basis matrix\n";)

   // collect columns of the basis inverse
   for(int i = 0; i < basisdim; ++i)
   {
      binvcol[i] = new VectorBase<R>(basisdim);
      binvcol[i]->clear();
      assert(getBasisInverseColReal(i, binvcol[i]->get_ptr(), 0, 0, true));
   }

   SPX_MSG_INFO1(spxout, spxout << "DEBUG: computing rows of inverted basis matrix\n";)

   // collect rows of the basis inverse
   for(int i = 0; i < basisdim; ++i)
   {
      binvrow[i] = new VectorBase<R>(basisdim);
      binvrow[i]->clear();
      assert(getBasisInverseRowReal(i, binvrow[i]->get_ptr(), 0, 0, true));
   }

   SPX_MSG_INFO1(spxout, spxout <<
                 "DEBUG: checking columns for identity after multiplying with basis matrix\n";)

   // multiply with (unscaled) basis matrix and check result (should be unitvecs)
   for(int i = 0; i < basisdim; ++i)
   {
      VectorBase<R> result(*binvcol[i]);
      assert(multBasis(result.get_ptr(), true));
      R sumerror = 0.0;

      for(int j = 0; j < basisdim; ++j)
      {
         R error = 0.0;

         if(j != i)
            error = spxAbs(result[j]);
         else
            error = spxAbs(result[j] - 1.0);

         if(error > _solver.tolerances()->floatingPointFeastol())
            SPX_MSG_INFO1(spxout, spxout << "ERROR: col " << i << " " << j << ", " << result[j] << std::endl);

         sumerror += error;
      }

      assert(LE(sumerror, 0, _solver.tolerances()->floatingPointFeastol()));
   }

   SPX_MSG_INFO1(spxout, spxout <<
                 "DEBUG: checking rows for identity after multiplying with basis matrix\n";)

   for(int i = 0; i < basisdim; ++i)
   {
      VectorBase<R> result(*binvrow[i]);
      assert(multBasisTranspose(result.get_ptr(), true));
      R sumerror = 0.0;

      for(int j = 0; j < basisdim; ++j)
      {
         R error = 0.0;

         if(j != i)
            error = spxAbs(result[j]);
         else
            error = spxAbs(result[j] - 1.0);

         if(error > _solver.tolerances()->floatingPointFeastol())
            SPX_MSG_INFO1(spxout, spxout << "ERROR: row " << i << " " << j << ", " << result[j] << std::endl);

         sumerror += error;
      }

      assert(LE(sumerror, 0, _solver.tolerances()->floatingPointFeastol()));
   }

   if(_solver.isScaled())
   {
      SPX_MSG_INFO1(spxout, spxout << "DEBUG: unscaling LP\n";)

      // unscale LP matrix
      _solver.unscaleLP();

      // invalidate basis factorization
      _solver.loadMatrixVecs();

      VectorBase<R>** binvcol2 = 0;
      VectorBase<R>** binvrow2 = 0;
      spx_alloc(binvcol2, basisdim);
      spx_alloc(binvrow2, basisdim);

      SPX_MSG_INFO1(spxout, spxout << "DEBUG: computing columns of inverted basis matrix again\n";)

      // collect columns of the basis inverse
      for(int i = 0; i < basisdim; ++i)
      {
         binvcol2[i] = new VectorBase<R>(basisdim);
         binvcol2[i]->clear();
         assert(getBasisInverseColReal(i, binvcol2[i]->get_ptr(), 0, 0, false));
      }

      SPX_MSG_INFO1(spxout, spxout << "DEBUG: computing rows of inverted basis matrix again\n";)

      // collect rows of the basis inverse
      for(int i = 0; i < basisdim; ++i)
      {
         binvrow2[i] = new VectorBase<R>(basisdim);
         binvrow2[i]->clear();
         assert(getBasisInverseRowReal(i, binvrow2[i]->get_ptr(), 0, 0, false));
      }

      SPX_MSG_INFO1(spxout, spxout <<
                    "DEBUG: checking rows and columns of scaled/unscaled inverted of basis matrix\n";)

      for(int i = 0; i < basisdim; ++i)
      {
         R sumerror = 0.0;

         for(int j = 0; j < basisdim; ++j)
         {
            if(NE((*binvcol[i])[j], (*binvcol2[i])[j], _solver.tolerances()->floatingPointFeastol()))
            {
               SPX_MSG_INFO1(spxout, spxout << "ERROR: col " << i << " " << j << ", " <<
                             (*binvcol[i])[j] << " " << (*binvcol2[i])[j] << std::endl);
               sumerror += spxAbs((*binvcol[i])[j] - (*binvcol2[i])[j]);
            }

            if(NE((*binvrow[i])[j], (*binvrow2[i])[j], _solver.tolerances()->floatingPointFeastol()))
            {
               SPX_MSG_INFO1(spxout, spxout << "ERROR: row " << i << " " << j << ", " <<
                             (*binvrow[i])[j] << " " << (*binvrow2[i])[j] << std::endl);
               sumerror += spxAbs((*binvrow[i])[j] - (*binvrow2[i])[j]);
            }
         }

         assert(LE(sumerror, 0, _solver.tolerances()->floatingPointFeastol()));
      }

      spx_free(binvcol2);
      spx_free(binvrow2);
   }

   spx_free(inds);
   spx_free(binvrow);
   spx_free(binvcol);
}

} // namespace soplex
