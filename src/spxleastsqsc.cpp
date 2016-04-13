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
#include <iostream> // TODO delete
#include "spxleastsqsc.h"
#include "spxout.h"
#include "basevectors.h"
#include "svsetbase.h"
#include "svectorbase.h"
#include "ssvectorbase.h"

#define MAX_SCALINGROUNDS 11
#define TERMINATION_FACTOR 0.001

namespace soplex
{
// colnnzeroes, resncols, tmpcols, csccurr, cscprev, qcurr, qprev, eprev[1], eprev[2]
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
   assert(qcurr * qprev != 0);

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

   /* swap poiners */
   pssv = psccurr;
   psccurr = pscprev;
   pscprev = pssv;
/*
      for(int i = 0; i < 100; ++i )
        std::cout << "psccurr: " << (*psccurr)[i] << " \n";
*/
}


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
   assert(psccurr != NULL);
   assert(pscprev != NULL);
   assert(q != 0);

   Real fac = -(eprev1 * eprev2);

   *pscprev -= *psccurr;

   if( isZero(fac) )
      (*pscprev).clear();
   else
      *pscprev *= fac;

   *pscprev += tmpvec.assignPWproduct4setup(resnvec, vecnnzeroes);
   *pscprev *= 1.0 / (q);
   *pscprev += *psccurr;

   psccurr = pscprev;
}



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

   tmpvec.assign2product4setup(facset, resvecprev); // TODO why setup??
   tmpvec.setup();
   resvec += tmpvec;

   resvec *= (-1.0 / qcurr);
}


/** initialize constant vectors and matrices */
static void initConstVecs(
   const SVSet* vecset,
   SVSet& facset,
   SSVector& veclogs,
   SSVector& vecnnzeroes,
   int nvec)
{
   assert(vecset != NULL);

   Real a;
   Real x;
   Real sum;
   int l;
   int size;
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

      veclogs.add(k, sum);

      assert(nnzeros > 0);

      x = (1.0 / nnzeros);
      vecnnzeroes.add(k, x);

      /* create new vector for facset */
      SVector& vecnew = (*(facset.create(nnzeros)));

      for( l = 0; l < size; ++l)
      {
         if( !isZero(lpvec.value(l)) )
            vecnew.add(lpvec.index(l), x);
      }
   }
   veclogs.setup();
   vecnnzeroes.setup();
}

static const char* makename(bool doBoth)
{
   return doBoth ? "Least squares" : "Least squares";
}

SPxLeastSqSC::SPxLeastSqSC(bool doBoth)
   : SPxScaler(makename(doBoth), false, doBoth)
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

void SPxLeastSqSC::scale(SPxLP& lp)
{

   MSG_INFO1( (*spxout), (*spxout) << "Least squares LP scaling" << std::endl; )

   setup(lp);

   Real tmp;
   Real qcurr;
   Real qprev;
   Real scurr;
   Real sprev;
   Real eprev[3];
   Real termination_constant;
   int k;
   int l;
   int nrows = lp.nRows();
   int ncols = lp.nCols();
   int nnzeroes = lp.nNzos();
   int maxscrounds = MAX_SCALINGROUNDS;

   /* constant factor matrices */
   SVSet facnrows(-1,-1,1.1,1.2); //TODO give memory
   SVSet facncols(-1,-1,1.1,1.2); //TODO give memory

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

#if 0
   MSG_INFO2( (*spxout), (*spxout) << "LP scaling statistics:"
      << " min= " << lp.minAbsNzo()
      << " max= " << lp.maxAbsNzo()
      << " col-ratio= " << colratio
      << " row-ratio= " << rowratio
      << std::endl; )
#endif

   /* initialize */

   std::cout << "nrows: " << nrows << " ncols: " << ncols << "\n";

   qcurr = 1.0;
   qprev = 0.0;
   termination_constant = TERMINATION_FACTOR * nnzeroes;

   for(k = 0; k < 3; ++k )
      eprev[k] = 0.0;

   initConstVecs(lp.rowSet(), facnrows, rowlogs, rownnzeroes, nrows);
   initConstVecs(lp.colSet(), facncols, collogs, colnnzeroes, ncols);

   tmprows.setup();
   tmpcols.setup();
   rowscale1.setup();
   rowscale2.setup();
   colscale1.setup();
   colscale2.setup();

   // compute first residual vector
   resncols = collogs - tmpcols.assign2product4setup(facnrows, rowlogs);

   resncols.setup();
   resnrows.setup();

   rowscale1.assignPWproduct4setup(rownnzeroes, rowlogs);
   rowscale2 = rowscale1;

#if 0
   SSVector a1(3);
   SSVector a2(3);

   SSVector a3(3);

   a1.setValue(0,7);
   a1.setValue(1,7);
   a1.setValue(2,1);

   std::cout << a1 * a2 << " \n";

   a2.setValue(0,1);

std::cout << a1 * a2 << " \n";

a2.setValue(1,3);
  // a2.setValue(2,2);

   std::cout << a1 * a2 << " \n";

   a1.setup();
   a2.setup();
   a3.setup();

   a3.assignPWproduct4setup(a1, a2);

   std::cout << a3[0] << " " << a3[1] << " " << a3[2] <<" \n";
   assert(0);
#endif

   scurr = resncols * tmpcols.assignPWproduct4setup(colnnzeroes, resncols);

   std::cout << "s(0): " << scurr  << " \n";

   resncols.setup();

   /* main loop */
   for(k = 0; k < maxscrounds; ++k )
   {
      sprev = scurr;
      std::cout << "round: " << k << " \n";
      /* is k even? */
      if( (k % 2) == 0 )
      {
         // not in first iteration?
         if( k != 0 )
            updateScale(rownnzeroes, resnrows, tmprows, rsccurr, rscprev, qcurr, qprev, eprev[1], eprev[2]);

         updateRes(facncols, resncols, resnrows, tmprows, eprev[0], qcurr);

         scurr = resnrows * tmprows.assignPWproduct4setup(resnrows, rownnzeroes);
	 std::cout << "s(k+1): " << scurr  << " \n";

      }
      else // k is odd
      {
         updateScale(colnnzeroes, resncols, tmpcols, csccurr, cscprev, qcurr, qprev, eprev[1], eprev[2]);
#if 0
      for(int i = 0; i < ncols; ++i )
         std::cout << "colscale: " << (*csccurr)[i] << " \n";
#endif
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
      if( scurr < termination_constant )
         break;
   }

   /* is k even? */
   if( (k % 2) == 0 )
   {
      /* update column scaling factor vector */
      updateScaleFinal(colnnzeroes, resncols, tmpcols, csccurr, cscprev, qprev, eprev[1], eprev[2]);
   }
   else // k is odd
   {
      /* update row scaling factor vector */
      updateScaleFinal(rownnzeroes, resnrows, tmprows, rsccurr, rscprev, qprev, eprev[1], eprev[2]);
   }

   /* compute actual scaling factors */

   SSVector rowscale = *rsccurr;
   SSVector colscale = *csccurr;

   for(k = 0; k < nrows; k++ )
      m_rowscale[k] = pow(2.0, - round(rowscale[k]));

   for(k = 0; k < ncols; k++ )
      m_colscale[k] = pow(2.0, - round(colscale[k]));

#if 0
    for(int i = 0; i < ncols; ++i )
         std::cout << "FINAL colscale: " << m_colscale[i] << " \n";

    for(int i = 0; i < nrows; ++i )
         std::cout << "FINAL rowscale: " << m_rowscale[i] << " \n";
#endif

   std::cout << "before scaling: min= " << lp.minAbsNzo();
   std::cout << " max= " << lp.maxAbsNzo() << "\n";

   /* scale */
   applyScaling(lp);

   std::cout << "after scaling: min= " << lp.minAbsNzo();
   std::cout << " max= " << lp.maxAbsNzo() << "\n";

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
