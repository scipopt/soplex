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

/**@file  spxleastsqsc.hpp
 * @brief LP least squares scaling.
 */

#include <assert.h>
#include <cmath>
#include "soplex/spxout.h"
#include "soplex/basevectors.h"
#include "soplex/svsetbase.h"
#include "soplex/svectorbase.h"
#include "soplex/ssvectorbase.h"
#include <array>

namespace soplex
{

/* update scaling VectorBase<R> */
template <class R>
static void updateScale(
   const SSVectorBase<R> vecnnzeroes,
   const SSVectorBase<R> resnvec,
   SSVectorBase<R>& tmpvec,
   SSVectorBase<R>*& psccurr,
   SSVectorBase<R>*& pscprev,
   R qcurr,
   R qprev,
   R eprev1,
   R eprev2,
   R epsilon)
{
   assert(psccurr != nullptr);
   assert(pscprev != nullptr);
   assert(qcurr * qprev != 0.0);

   R fac = -(eprev1 * eprev2);

   SSVectorBase<R>* pssv;

   *pscprev -= *psccurr;

   if(isZero(fac, epsilon))
      (*pscprev).clear();
   else
      *pscprev *= fac;

   *pscprev += tmpvec.assignPWproduct4setup(resnvec, vecnnzeroes);

   *pscprev *= 1.0 / (qcurr * qprev);
   *pscprev += *psccurr;

   /* swap pointers */
   pssv = psccurr;
   psccurr = pscprev;
   pscprev = pssv;
}

/* update scaling VectorBase<R> after main loop */
template <class R>
static void updateScaleFinal(
   const SSVectorBase<R> vecnnzeroes,
   const SSVectorBase<R> resnvec,
   SSVectorBase<R>& tmpvec,
   SSVectorBase<R>*& psccurr,
   SSVectorBase<R>*& pscprev,
   R q,
   R eprev1,
   R eprev2,
   R epsilon)
{
   assert(q != 0);
   assert(psccurr != nullptr);
   assert(pscprev != nullptr);

   R fac = -(eprev1 * eprev2);

   *pscprev -= *psccurr;

   if(isZero(fac, epsilon))
      (*pscprev).clear();
   else
      *pscprev *= fac;

   *pscprev += tmpvec.assignPWproduct4setup(resnvec, vecnnzeroes);
   *pscprev *= 1.0 / q;
   *pscprev += *psccurr;

   psccurr = pscprev;
}

/* update residual VectorBase<R> */
template <class R>
static inline void updateRes(
   const SVSetBase<R> facset,
   const SSVectorBase<R> resvecprev,
   SSVectorBase<R>& resvec,
   SSVectorBase<R>& tmpvec,
   R eprev,
   R qcurr,
   R epsilon)
{
   assert(qcurr != 0.0);

   if(isZero(eprev, epsilon))
      resvec.clear();
   else
      resvec *= eprev;

   int dummy1 = 0;
   int dummy2 = 0;
   tmpvec.assign2product4setup(facset, resvecprev, nullptr, nullptr, dummy1, dummy2);
   tmpvec.setup();
   resvec += tmpvec;

   resvec *= (-1.0 / qcurr);
   resvec.setup();
}


/* initialize constant vectors and matrices */
template <class R>
static void initConstVecs(
   const SVSetBase<R>* vecset,
   SVSetBase<R>& facset,
   SSVectorBase<R>& veclogs,
   SSVectorBase<R>& vecnnzinv,
   R epsilon)
{
   assert(vecset != nullptr);

   const int nvec = vecset->num();

   for(int k = 0; k < nvec; ++k)
   {
      R logsum = 0.0;
      int nnz = 0;
      // get kth row or column of LP
      const SVectorBase<R>& lpvec = (*vecset)[k];
      const int size = lpvec.size();

      for(int i = 0; i < size; ++i)
      {
         const R a = lpvec.value(i);

         if(!isZero(a, epsilon))
         {
            logsum += log2(double(spxAbs(a))); // todo spxLog2?
            nnz++;
         }
      }

      R nnzinv;

      if(nnz > 0)
      {
         nnzinv = 1.0 / nnz;
      }
      else
      {
         /* all-0 entries */
         logsum = 1.0;
         nnzinv = 1.0;
      }

      veclogs.add(k, logsum);
      vecnnzinv.add(k, nnzinv);

      /* create new VectorBase<R> for facset */
      SVectorBase<R>& vecnew = (*(facset.create(nnz)));

      for(int i = 0; i < size; ++i)
      {
         if(!isZero(lpvec.value(i), epsilon))
            vecnew.add(lpvec.index(i), nnzinv);
      }

      vecnew.sort();
   }

   assert(veclogs.isSetup());
   assert(vecnnzinv.isSetup());
}

/* return name of scaler */
static inline const char* makename()
{
   return "Least squares";
}

template <class R>
SPxLeastSqSC<R>::SPxLeastSqSC()
   : SPxScaler<R>(makename(), false, false)
{}

template <class R>
SPxLeastSqSC<R>::SPxLeastSqSC(const SPxLeastSqSC<R>& old)
   : SPxScaler<R>(old), acrcydivisor(old.acrcydivisor), maxrounds(old.maxrounds)
{}

template <class R>
SPxLeastSqSC<R>& SPxLeastSqSC<R>::operator=(const SPxLeastSqSC<R>& rhs)
{
   if(this != &rhs)
   {
      SPxScaler<R>::operator=(rhs);
   }

   return *this;
}

template <class R>
void SPxLeastSqSC<R>::setRealParam(R param, const char* name)
{
   assert(param >= 1.0);
   acrcydivisor = param;
}

template <class R>
void SPxLeastSqSC<R>::setIntParam(int param, const char* name)
{
   assert(param >= 0);
   maxrounds = param;
}

// todo refactor this method. Has no abstraction and is too long
template <class R>
void SPxLeastSqSC<R>::scale(SPxLPBase<R>& lp,  bool persistent)
{
   SPX_MSG_INFO1((*this->spxout), (*this->spxout) << "Least squares LP scaling" <<
                 (persistent ? " (persistent)" : "") << std::endl;)

   this->setup(lp);

   const int nrows = lp.nRows();
   const int ncols = lp.nCols();
   const int lpnnz = lp.nNzos();

   // is contraints matrix empty?
   //todo don't create the scaler in this case!
   if(lpnnz == 0)
   {
      // to keep the invariants, we still need to call this method
      this->applyScaling(lp);

      return;
   }

   assert(nrows > 0 && ncols > 0 && lpnnz > 0);

   /* constant factor matrices;
    * in Curtis-Reid article
    * facnrows equals E^T M^(-1)
    * facncols equals E N^(-1)
    * */
   SVSetBase<R> facnrows(nrows, nrows, 1.1, 1.2);
   SVSetBase<R> facncols(ncols, ncols, 1.1, 1.2);

   /* column scaling factor vectors */
   SSVectorBase<R> colscale1(ncols, this->_tolerances);
   SSVectorBase<R> colscale2(ncols, this->_tolerances);

   /* row scaling factor vectors */
   SSVectorBase<R> rowscale1(nrows, this->_tolerances);
   SSVectorBase<R> rowscale2(nrows, this->_tolerances);

   /* residual vectors */
   SSVectorBase<R> resnrows(nrows, this->_tolerances);
   SSVectorBase<R> resncols(ncols, this->_tolerances);

   /* vectors to store temporary values */
   SSVectorBase<R> tmprows(nrows, this->_tolerances);
   SSVectorBase<R> tmpcols(ncols, this->_tolerances);

   /* vectors storing the row and column sums (respectively) of logarithms of
    *(absolute values of) non-zero elements of left hand matrix of LP
    */
   SSVectorBase<R> rowlogs(nrows, this->_tolerances);
   SSVectorBase<R> collogs(ncols, this->_tolerances);

   /* vectors storing the inverted number of non-zeros in each row and column
    *(respectively) of left hand matrix of LP
    */
   SSVectorBase<R> rownnzinv(nrows, this->_tolerances);
   SSVectorBase<R> colnnzinv(ncols, this->_tolerances);

   /* VectorBase<R> pointers */
   SSVectorBase<R>* csccurr = &colscale1;
   SSVectorBase<R>* cscprev = &colscale2;
   SSVectorBase<R>* rsccurr = &rowscale1;
   SSVectorBase<R>* rscprev = &rowscale2;

   SPX_MSG_INFO2((*this->spxout), (*this->spxout) << "before scaling:"
                 << " min= " << lp.minAbsNzo()
                 << " max= " << lp.maxAbsNzo()
                 << " col-ratio= " << this->maxColRatio(lp)
                 << " row-ratio= " << this->maxRowRatio(lp)
                 << std::endl;)

   /* initialize scalars, vectors and matrices */

   assert(acrcydivisor > 0.0);

   const R smax = lpnnz / acrcydivisor;
   R qcurr = 1.0;
   R qprev = 0.0;

   std::array<R, 3> eprev;
   eprev.fill(0.0);

   initConstVecs(lp.rowSet(), facnrows, rowlogs, rownnzinv, R(this->tolerances()->epsilon()));
   initConstVecs(lp.colSet(), facncols, collogs, colnnzinv, R(this->tolerances()->epsilon()));

   assert(tmprows.isSetup() && tmpcols.isSetup());
   assert(rowscale1.isSetup() && rowscale2.isSetup());
   assert(colscale1.isSetup() && colscale2.isSetup());

   // compute first residual vector r0
   int dummy1 = 0;
   int dummy2 = 0;
   resncols = collogs - tmpcols.assign2product4setup(facnrows, rowlogs, nullptr, nullptr, dummy1,
              dummy2);

   resncols.setup();
   resnrows.setup();

   rowscale1.assignPWproduct4setup(rownnzinv, rowlogs);
   rowscale2 = rowscale1;

   R scurr = resncols * tmpcols.assignPWproduct4setup(colnnzinv, resncols);

   int k;

   /* conjugate gradient loop */
   for(k = 0; k < maxrounds; ++k)
   {
      const R sprev = scurr;

      // termination criterion met?
      if(scurr < smax)
         break;

      // is k even?
      if((k % 2) == 0)
      {
         // not in first iteration?
         if(k != 0)   // true, then update row scaling factor vector
            updateScale(rownnzinv, resnrows, tmprows, rsccurr, rscprev, qcurr, qprev, eprev[1], eprev[2],
                        R(this->tolerances()->epsilon()));

         updateRes(facncols, resncols, resnrows, tmprows, eprev[0], qcurr, R(this->tolerances()->epsilon()));
         scurr = resnrows * tmprows.assignPWproduct4setup(resnrows, rownnzinv);
      }
      else // k is odd
      {
         // update column scaling factor vector
         updateScale(colnnzinv, resncols, tmpcols, csccurr, cscprev, qcurr, qprev, eprev[1], eprev[2],
                     R(this->tolerances()->epsilon()));

         updateRes(facnrows, resnrows, resncols, tmpcols, eprev[0], qcurr, R(this->tolerances()->epsilon()));
         scurr = resncols * tmpcols.assignPWproduct4setup(resncols, colnnzinv);
      }

      // shift eprev entries one to the right
      for(unsigned l = 2; l > 0; --l)
         eprev[l] = eprev[l - 1];

      assert(isNotZero(sprev, R(this->tolerances()->epsilon())));

      eprev[0] = (qcurr * scurr) / sprev;

      const R tmp = qcurr;
      qcurr = 1.0 - eprev[0];
      qprev = tmp;
   }

   if(k > 0 && (k % 2) == 0)
   {
      // update column scaling factor vector
      updateScaleFinal(colnnzinv, resncols, tmpcols, csccurr, cscprev, qprev, eprev[1], eprev[2],
                       R(this->tolerances()->epsilon()));
   }
   else if(k > 0)
   {
      // update row scaling factor vector
      updateScaleFinal(rownnzinv, resnrows, tmprows, rsccurr, rscprev, qprev, eprev[1], eprev[2],
                       R(this->tolerances()->epsilon()));
   }

   /* compute actual scaling factors */

   const SSVectorBase<R>& rowscale = *rsccurr;
   const SSVectorBase<R>& colscale = *csccurr;

   DataArray<int>& colscaleExp = *this->m_activeColscaleExp;
   DataArray<int>& rowscaleExp = *this->m_activeRowscaleExp;

   for(k = 0; k < nrows; ++k)
      rowscaleExp[k] = -int(rowscale[k] + ((rowscale[k] >= 0.0) ? (+0.5) : (-0.5)));

   for(k = 0; k < ncols; ++k)
      colscaleExp[k] = -int(colscale[k] + ((colscale[k] >= 0.0) ? (+0.5) : (-0.5)));

   // scale
   this->applyScaling(lp);

   SPX_MSG_INFO3((*this->spxout), (*this->spxout) << "Row scaling min= " << this->minAbsRowscale()
                 << " max= " << this->maxAbsRowscale()
                 << std::endl
                 << "Col scaling min= " << this->minAbsColscale()
                 << " max= " << this->maxAbsColscale()
                 << std::endl;)

   SPX_MSG_INFO2((*this->spxout), (*this->spxout) << "after scaling: "
                 << " min= " << lp.minAbsNzo(false)
                 << " max= " << lp.maxAbsNzo(false)
                 << " col-ratio= " << this->maxColRatio(lp)
                 << " row-ratio= " << this->maxRowRatio(lp)
                 << std::endl;)
}

} // namespace soplex
