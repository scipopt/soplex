/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the class library                   */
/*       SoPlex --- the Sequential object-oriented simPlex.                  */
/*                                                                           */
/*    Copyright (C) 1996-2016 Konrad-Zuse-Zentrum                            */
/*                            fuer Informationstechnik Berlin                */
/*                                                                           */
/*  SoPlex is distributed under the terms of the ZIB Academic Licence.       */
/*                                                                           */
/*  You should have received a copy of the ZIB Academic License              */
/*  along with SoPlex; see the file COPYING. If not email to soplex@zib.de.  */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/**@file  spxleastsqsc.cpp
 * @brief LP least squares scaling.
 */
#include <cmath>
#include <assert.h>
#include "spxleastsqsc.h"
#include "spxout.h"
#include "basevectors.h"
#include "svsetbase.h"
#include "svectorbase.h"
#include "ssvectorbase.h"

namespace soplex
{

/* update scaling vector */
static inline void updateScale(
   const SSVector vecnnzeroes,
   const SSVector resnvec,
   SSVector& tmpvec,
   SSVector*& psccurr,
   SSVector*& pscprev,
   Real qcurr,
   Real qprev,
   Real eprev1,
   Real eprev2)
{
   assert(psccurr != NULL);
   assert(pscprev != NULL);
   assert(qcurr * qprev != 0.0);

   Real fac = -(eprev1 * eprev2);

   SSVector* pssv;

   *pscprev -= *psccurr;

   if( isZero(fac) )
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

/* update scaling vector after main loop */
static inline void updateScaleFinal(
   const SSVector vecnnzeroes,
   const SSVector resnvec,
   SSVector& tmpvec,
   SSVector*& psccurr,
   SSVector*& pscprev,
   Real q,
   Real eprev1,
   Real eprev2)
{
   assert(q != 0);
   assert(psccurr != NULL);
   assert(pscprev != NULL);

   Real fac = -(eprev1 * eprev2);

   *pscprev -= *psccurr;

   if( isZero(fac) )
      (*pscprev).clear();
   else
      *pscprev *= fac;

   *pscprev += tmpvec.assignPWproduct4setup(resnvec, vecnnzeroes);
   *pscprev *= 1.0 / q;
   *pscprev += *psccurr;

   psccurr = pscprev;
}

/* update residual vector */
static inline void updateRes(
   const SVSet facset,
   const SSVector resvecprev,
   SSVector& resvec,
   SSVector& tmpvec,
   Real eprev,
   Real qcurr)
{
   assert(qcurr != 0);

   if( isZero(eprev) )
      resvec.clear();
   else
      resvec *= eprev;

   tmpvec.assign2product4setup(facset, resvecprev);
   tmpvec.setup();
   resvec += tmpvec;

   resvec *= (-1.0 / qcurr);
   resvec.setup();
}


/* initialize constant vectors and matrices */
static void initConstVecs(
   const SVSet* vecset,
   SVSet& facset,
   SSVector& veclogs,
   SSVector& vecnnzeroes)
{
   assert(vecset != NULL);

   Real a;
   Real x;
   Real sum;
   int l;
   int size;
   int nvec = vecset->num();
   int nnzeros;

   for(int k = 0; k < nvec; ++k )
   {
      sum = 0.0;
      nnzeros = 0;
      const SVector& lpvec = (*vecset)[k];

      size = lpvec.size();

      for( l = 0; l < size; ++l)
      {
         a = spxAbs(lpvec.value(l));

         if( !isZero(a) )
         {
            sum += log2(a);
            nnzeros++;
         }
      }

      if( nnzeros > 0)
      {
         x = (1.0 / nnzeros);
      }
      else
      {
         /* all-0 entries, so assume row is already scaled (all-1) */
         sum = (Real) size;
         x = 1.0 / size;
      }

      veclogs.add(k, sum);

      vecnnzeroes.add(k, x);

      /* create new vector for facset */
      SVector& vecnew = (*(facset.create(nnzeros)));

      for( l = 0; l < size; ++l)
      {
         if( !isZero(lpvec.value(l)) )
            vecnew.add(lpvec.index(l), x);
      }
      vecnew.sort();
   }

   assert(veclogs.isSetup());
   assert(vecnnzeroes.isSetup());
}

/* return name of scaler */
static const char* makename()
{
   return "Least squares";
}

SPxLeastSqSC::SPxLeastSqSC()
   : SPxScaler(makename(), false, false)
{}

SPxLeastSqSC::SPxLeastSqSC(const SPxLeastSqSC& old)
   : SPxScaler(old)
{}

SPxLeastSqSC& SPxLeastSqSC::operator=(const SPxLeastSqSC& rhs)
{
   if(this != &rhs)
   {
      SPxScaler::operator=(rhs);
   }

   return *this;
}

Real SPxLeastSqSC::computeScale(Real /*mini*/, Real maxi) const
{
   return maxi;
}

void SPxLeastSqSC::setRealParam(Real param, const char* name)
{
   assert(param >= 1.0);
   acrcydivisor = param;
}

void SPxLeastSqSC::setIntParam(int param, const char* name)
{
   assert(param >= 0);
   maxrounds = param;
}

void SPxLeastSqSC::scale(SPxLP& lp)
{
   MSG_INFO1( (*spxout), (*spxout) << "Least squares LP scaling" << std::endl; )

   setup(lp);

   Real tmp;
   Real smax;
   Real qcurr;
   Real qprev;
   Real scurr;
   Real sprev;
   Real eprev[3];
   int k;
   int l;
   int nrows = lp.nRows();
   int ncols = lp.nCols();
   int nnzeroes = lp.nNzos();
   int maxscrounds = maxrounds;

   /* constant factor matrices */
   SVSet facnrows(nrows, nrows, 1.1, 1.2);
   SVSet facncols(ncols, ncols, 1.1, 1.2);

   /* column scaling factor vectors */
   SSVector colscale1(ncols);
   SSVector colscale2(ncols);

   /* row scaling factor vectors */
   SSVector rowscale1(nrows);
   SSVector rowscale2(nrows);

   /* residual vectors */
   SSVector resnrows(nrows);
   SSVector resncols(ncols);

   /* vectors to store temporary values */
   SSVector tmprows(nrows);
   SSVector tmpcols(ncols);

   /* vectors storing the row and column sums (respectively) of logarithms of
    *(absolute values of) non-zero elements of left hand matrix of LP
    */
   SSVector rowlogs(nrows);
   SSVector collogs(ncols);

   /* vectors storing the inverted number of non-zeros in each row and column
    *(respectively) of left hand matrix of LP
    */
   SSVector rownnzeroes(nrows);
   SSVector colnnzeroes(ncols);

   /* vector pointers */
   SSVector* csccurr = &colscale1;
   SSVector* cscprev = &colscale2;
   SSVector* rsccurr = &rowscale1;
   SSVector* rscprev = &rowscale2;

   MSG_INFO2( (*spxout), (*spxout) << "LP scaling statistics:"
      << " min= " << lp.minAbsNzo()
      << " max= " << lp.maxAbsNzo()
      << " col-ratio= " << maxColRatio(lp)
      << " row-ratio= " << maxRowRatio(lp)
      << std::endl; )

   /* initialize scalars, vectors and matrices */

   smax = nnzeroes / acrcydivisor;
   qcurr = 1.0;
   qprev = 0.0;

   for(k = 0; k < 3; ++k )
      eprev[k] = 0.0;

   initConstVecs(lp.rowSet(), facnrows, rowlogs, rownnzeroes);
   initConstVecs(lp.colSet(), facncols, collogs, colnnzeroes);

   assert(tmprows.isSetup());
   assert(tmpcols.isSetup());
   assert(rowscale1.isSetup());
   assert(rowscale2.isSetup());
   assert(colscale1.isSetup());
   assert(colscale2.isSetup());

   // compute first residual vector
   resncols = collogs - tmpcols.assign2product4setup(facnrows, rowlogs);

   resncols.setup();
   resnrows.setup();

   rowscale1.assignPWproduct4setup(rownnzeroes, rowlogs);
   rowscale2 = rowscale1;

   scurr = resncols * tmpcols.assignPWproduct4setup(colnnzeroes, resncols);

   /* conjugate gradient loop */
   for( k = 0; k < maxscrounds; ++k )
   {
      sprev = scurr;

      // is k even?
      if( (k % 2) == 0 )
      {
         // not in first iteration?
         if( k != 0 ) // true, then update row scaling factor vector
            updateScale(rownnzeroes, resnrows, tmprows, rsccurr, rscprev, qcurr, qprev, eprev[1], eprev[2]);

         updateRes(facncols, resncols, resnrows, tmprows, eprev[0], qcurr);

         scurr = resnrows * tmprows.assignPWproduct4setup(resnrows, rownnzeroes);
      }
      else // k is odd
      {
         // update column scaling factor vector
         updateScale(colnnzeroes, resncols, tmpcols, csccurr, cscprev, qcurr, qprev, eprev[1], eprev[2]);

         updateRes(facnrows, resnrows, resncols, tmpcols, eprev[0], qcurr);
         scurr = resncols * (tmpcols.assignPWproduct4setup(resncols, colnnzeroes) );
      }

      // shift eprev entries one to the right
      for( l = 2; l > 0; --l)
         eprev[l] = eprev[l - 1];

      eprev[0] = (qcurr * scurr) / sprev;

      tmp = qcurr;
      qcurr = 1.0 - eprev[0];
      qprev = tmp;

      // termination criterion met?
      if( scurr < smax )
         break;
   }

   // is k even?
   if( (k % 2) == 0 )
   {
      // update column scaling factor vector
      updateScaleFinal(colnnzeroes, resncols, tmpcols, csccurr, cscprev, qprev, eprev[1], eprev[2]);
   }
   else // k is odd
   {
      // update row scaling factor vector
      updateScaleFinal(rownnzeroes, resnrows, tmprows, rsccurr, rscprev, qprev, eprev[1], eprev[2]);
   }

   /* compute actual scaling factors */

   SSVector rowscale = *rsccurr;
   SSVector colscale = *csccurr;

   for(k = 0; k < nrows; ++k )
      m_rowscale[k] = pow(2.0, - round(rowscale[k]));

   for(k = 0; k < ncols; ++k )
      m_colscale[k] = pow(2.0, - round(colscale[k]));

   // scale
   applyScaling(lp);

   MSG_INFO3( (*spxout), (*spxout) << "Row scaling min= " << minAbsRowscale()
      << " max= " << maxAbsRowscale()
      << std::endl
      << "\tCol scaling min= " << minAbsColscale()
      << " max= " << maxAbsColscale()
      << std::endl; )

   MSG_INFO2( (*spxout), (*spxout) << "LP scaling statistics:"
      << " min= " << lp.minAbsNzo()
      << " max= " << lp.maxAbsNzo()
      << " col-ratio= " << maxColRatio(lp)
      << " row-ratio= " << maxRowRatio(lp)
      << std::endl; )
}

} // namespace soplex
