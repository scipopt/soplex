/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the class library                   */
/*       SoPlex --- the Sequential object-oriented simPlex.                  */
/*                                                                           */
/*    Copyright (C) 1996-2017 Konrad-Zuse-Zentrum                            */
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

  template <>
  void SPxSolver<Real>::qualConstraintViolation(Real& maxviol, Real& sumviol) const
  {
    maxviol = 0.0;
    sumviol = 0.0;

    DVector solu( this->nCols() );

    getPrimal( solu );

    for( int row = 0; row < this->nRows(); ++row )
      {
        const SVector& rowvec = this->rowVector( row );

        Real val = 0.0;         

        for( int col = 0; col < rowvec.size(); ++col )
          val += rowvec.value( col ) * solu[rowvec.index( col )];

        Real viol = 0.0;

        assert(lhs( row ) <= rhs( row ) + (100 * epsilon()));

        if (val < this->lhs( row )) 
          viol = spxAbs(val - this->lhs( row ));
        else
          if (val > this->rhs( row ))
            viol = spxAbs(val - this->rhs( row ));

        if (viol > maxviol)
          maxviol = viol;

        sumviol += viol;
      }
  }

  template <>
  void SPxSolver<Real>::qualBoundViolation(
                                        Real& maxviol, Real& sumviol) const
  {
    maxviol = 0.0;
    sumviol = 0.0;

    DVector solu( this->nCols() );

    getPrimal( solu );

    for( int col = 0; col < this->nCols(); ++col )
      {
        assert(lower( col ) <= upper( col ) + (100 * epsilon()));

        Real viol = 0.0;

        if (solu[col] < this->lower( col ))
          viol = spxAbs( solu[col] - this->lower( col ));
        else
          if (solu[col] > this->upper( col ))
            viol = spxAbs( solu[col] - this->upper( col ));
         
        if (viol > maxviol)
          maxviol = viol;

        sumviol += viol;
      }
  }

  template <>
  void SPxSolver<Real>::qualSlackViolation(Real& maxviol, Real& sumviol) const
  {
    maxviol = 0.0;
    sumviol = 0.0;

    DVector solu( this->nCols() );
    DVector slacks( this->nRows() );

    getPrimal( solu );
    getSlacks( slacks );

    for( int row = 0; row < this->nRows(); ++row )
      {
        const SVector& rowvec = this->rowVector( row );

        Real val = 0.0;         

        for( int col = 0; col < rowvec.size(); ++col )
          val += rowvec.value( col ) * solu[rowvec.index( col )];

        Real viol = spxAbs(val - slacks[row]);

        if (viol > maxviol)
          maxviol = viol;

        sumviol += viol;
      }
  }

  template <>
  void SPxSolver<Real>::qualRedCostViolation(Real& maxviol, Real& sumviol) const
  {   
    maxviol = 0.0;
    sumviol = 0.0;

    int i;
    // TODO:   y = c_B * B^-1  => coSolve(y, c_B)
    //         redcost = c_N - yA_N 
    // solve system "x = e_i^T * B^-1" to get i'th row of B^-1
    // DVector y( this->nRows() );
    // basis().coSolve( x, spx->unitVector( i ) );
    // DVector rdcost( this->nCols() );
#if 0 // un-const
    if (lastUpdate() > 0)
      factorize();

    computePvec();

    if (type() == ENTER)
      computeTest();
#endif
    if (type() == ENTER)
      {
        for(i = 0; i < dim(); ++i)
          {
            Real x = coTest()[i];
         
            if (x < 0.0)
              {
                sumviol -= x;
            
                if (x < maxviol)
                  maxviol = x;
              }
          }
        for(i = 0; i < coDim(); ++i)
          {
            Real x = test()[i];
         
            if (x < 0.0)
              {
                sumviol -= x;
            
                if (x < maxviol)
                  maxviol = x;
              }
          } 
      }
    else
      {
        assert(type() == LEAVE);

        for(i = 0; i < dim(); ++i)
          {
            Real x = fTest()[i];
         
            if (x < 0.0)
              {
                sumviol -= x;
            
                if (x < maxviol)
                  maxviol = x;
              }
          }
      }
    maxviol *= -1;
  }

} // namespace soplex
