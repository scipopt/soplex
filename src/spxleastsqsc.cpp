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


#define MAX_SCALINGROUNDS 20

namespace soplex
{
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

   Real q;
   Real scurr;
   Real sprev;
   Real eprev[4];
   int k;
   int l;
   int nrows = (lp.rowSet())->num();
   int ncols = (lp.colSet())->num();
   int nnzeros;
   int maxscrounds = MAX_SCALINGROUNDS;

   std::cout << "nrows: " << nrows << " ncols: " << ncols;

   SVSet* vecset;

   /* constant factor matrices */
   SVSet facnrows();
   SVSet facncols();

   /* row scaling factor vectors */
   SSVector c1(ncols);
   SSVector c2(ncols);

   /* column scaling factor vectors */
   SSVector p1(nrows);
   SSVector p2(nrows);

   /* residual vectors */
   SSVector res1(nrows);
   SSVector res2(ncols);

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
   SSVector* ccurr;
   SSVector* cprev;
   SSVector* pcurr;
   SSVector* pprev;
   SSVector* rescurr;
   SSVector* resprev;

#if 0
   MSG_INFO2( (*spxout), (*spxout) << "LP scaling statistics:"
                        << " min= " << lp.minAbsNzo()
                        << " max= " << lp.maxAbsNzo()
                        << " col-ratio= " << colratio
                        << " row-ratio= " << rowratio
                        << std::endl; )
#endif
   /* initialize */
   Real a;
   Real sum;

   for( k = 0; k < 3; ++k )
      eprev[k] = 0;

   vecset = lp.rowSet();

   // TODO extra static method
   for(k = 0; k < nrows; ++k )
   {
      sum = 0.0;
      nnzeros = 0;
      SVector& vec = (*vecset)[k];

      assert(vec.size() == ncols);

      for( l = 0; l < ncols; ++l)
      {
         a = spxAbs(vec.value(l));

         if( !isZero(a) )
         {
            sum += log2(a);
            nnzeros++;
         }
      }

      rowlogs.setValue(k, sum);
      rownnzeroes.setValue(k, (1.0 / nnzeros));
   }

   vecset = lp.colSet()

   // TODO extra static method
   for(k = 0; k < ncols; ++k )
   {
      sum = 0.0;
      nnzeros = 0;
      SVector& vec = (*vecset)[k];

      assert(vec.size() == nrows);

      for( l = 0; l < nrows; ++l)
      {
         a = spxAbs(vec.value(l));

         if( !isZero(a) )
         {
            sum += log2(a);
            nnzeros++;
         }
      }

      collogs.setValue(k, sum);
      colnnzeroes.setValue(k, (1.0 / nnzeros));
   }

   for(k = 0; k < nrows; ++k )
      (void) facnrows.create(nrows);

   for(k = 0; k < ncols; ++k )
      (void) facncols.create(ncols);

   /* main loop */
   for( k = 0; k < maxscrounds; ++k )
   {



      /* is k even? */
      if( (k % 2) == 0 )
      {
      }
      else // k is odd
      {
      }

      /* shift eprev entries one to the right */
      for( l = 3; l > 0; --l )
         eprev[l] = eprev[l - 1];

      eprev[0] = (q * scurr) / sprev;
      q = 1 - eprev[0];

      /* termination criterion met? */
      if( )
          break;
   }

    /* is k even? */
    if( (k % 2) == 0 )
    {
       /* update column??? scaling factor vector */
    }
    else // k is odd
    {
       /* update column??? scaling factor vector */
    }

   if (colFirst)
   {
      computeScalingVecs(lp.colSet(), m_rowscale, m_colscale);

      if (m_doBoth)
         computeScalingVecs(lp.rowSet(), m_colscale, m_rowscale);
   }
   else
   {
      computeScalingVecs(lp.rowSet(), m_colscale, m_rowscale);

      if (m_doBoth)
         computeScalingVecs(lp.colSet(), m_rowscale, m_colscale);
   }
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
