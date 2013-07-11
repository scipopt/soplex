/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the class library                   */
/*       SoPlex --- the Sequential object-oriented simPlex.                  */
/*                                                                           */
/*    Copyright (C) 1996-2013 Konrad-Zuse-Zentrum                            */
/*                            fuer Informationstechnik Berlin                */
/*                                                                           */
/*  SoPlex is distributed under the terms of the ZIB Academic Licence.       */
/*                                                                           */
/*  You should have received a copy of the ZIB Academic License              */
/*  along with SoPlex; see the file COPYING. If not email to soplex@zib.de.  */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

#include <iostream>
#include <assert.h>

#include "soplex2.h"

namespace soplex
{
#if 0
   /// solves rational LP
   SPxSolver::Status SoPlex2::solveRational()
   {
      bool hasUnboundedRay = false;
      bool infeasibilityNotCertified = false;
      bool unboundednessNotCertified = false;
      _statusRational = SPxSolver::UNKNOWN;

      // introduce slack variables to transform inequality constraints into equations
      _transformEqualityRational();

      do
      {
         bool primalFeasible = false;
         bool dualFeasible = false;
         bool infeasible = false;
         bool unbounded = false;
         bool stopped = false;
         bool error = false;

         // solve problem with iterative refinement and recovery mechanism
         _performOptIRStable(primalFeasible, dualFeasible, infeasible, unbounded, stopped, error);

         // case: an unrecoverable error occured
         if( error )
         {
            MSG_INFO1( spxout << "failed during feasibility problem\n" );
            _statusRational = SPxSolver::ERROR;
            break;
         }
         // case: stopped due to some limit
         else if( stopped )
         {
            MSG_INFO1( spxout << "stopped solving\n" );
            _statusRational = SPxSolver::ABORT_TIME;
            break;
         }
         // case: unboundedness detected for the first time
         else if( unbounded && !unboundednessNotCertified )
         {
            _performUnboundedIRStable(hasUnboundedRay, stopped, error);

            if( error )
            {
               MSG_INFO1( spxout << "failed while trying to compute primal unbounded ray\n" );
               _statusRational = SPxSolver::ERROR;
               break;
            }

            assert(!hasUnboundedRay || _solRational.hasPrimalray());
            assert(!_solRational.hasPrimalray() || hasUnboundedRay);

            unboundednessNotCertified = !hasUnboundedRay;

            if( stopped )
            {
               MSG_INFO1( spxout << "stopped solving\n" );
               _statusRational = SPxSolver::ABORT_TIME;
               break;
            }

            _performFeasIRStable(infeasible, stopped, error);

            if( error )
            {
               MSG_INFO1( spxout << "failed during feasibility problem\n" );
               _statusRational = SPxSolver::ERROR;
               break;
            }
            else if( stopped )
            {
               MSG_INFO1( spxout << "stopped solving\n" );
               _statusRational = SPxSolver::ABORT_TIME;
               break;
            }
            else if( infeasible )
            {
               MSG_INFO1( spxout << "proved infeasiblity\n" );
               _statusRational = SPxSolver::INFEASIBLE;
               break;
            }
            else if( hasUnboundedRay )
            {
               MSG_INFO1( spxout << "proved primal unboundedness\n" );
               _statusRational = SPxSolver::UNBOUNDED;
               break;
            }
            else
            {
               MSG_INFO1( spxout << "unboundedness was not confirmed\n" );
               continue;
            }
         }
         // case: infeasibility detected
         else if( infeasible && !infeasibilityNotCertified )
         {
            _performFeasIRStable(infeasible, stopped, error);

            if( error )
            {
               MSG_INFO1( spxout << "failed during feasibility problem\n" );
               _statusRational = SPxSolver::ERROR;
               break;
            }

            infeasibilityNotCertified = !infeasible;

            if( stopped )
            {
               MSG_INFO1( spxout << "stopped solving\n" );
               _statusRational = SPxSolver::ABORT_TIME;
               break;
            }
            else if( infeasible )
            {
               MSG_INFO1( spxout << "proved infeasiblity\n" );
               _statusRational = SPxSolver::INFEASIBLE;
               break;
            }
            else if( hasUnboundedRay )
            {
               MSG_INFO1( spxout << "proved primal unboundedness\n" );
               _statusRational = SPxSolver::UNBOUNDED;
               break;
            }
            else
            {
               MSG_INFO1( spxout << "infeasibility was not confirmed\n" );
               continue;
            }
         }
         else
         {
            MSG_INFO1( spxout << "an unknown error occured\n" );
            continue;
         }
      }
      while( _isSolveStopped() );

      if( _isSolveStopped() )
         _statusRational = SPxSolver::ABORT_TIME;

      // restore original problem
      _untransformEqualityRational();

      return _statusRational;
   }
#endif



   /// introduces slack variables to transform inequality constraints into equations
   void SoPlex2::_transformEqualityRational()
   {
      _slackCols.clear();

      MSG_DEBUG( spxout << "adding slack columns\n" );

      // add artificial slack variables to convert inequality to equality constraints
      for( int i = 0; i < numRowsRational(); i++ )
      {
         if( lhsRational(i) != rhsRational(i) )
         {
            _slackCols.add(0.0, -rhsRational(i), DSVectorRational(UnitVector(i)), -lhsRational(i));
            _rationalLP->changeLhs(i, 0.0);
            _rationalLP->changeRhs(i, 0.0);
         }
      }

      _rationalLP->addCols(_slackCols);

      // adjust basis
      if( _hasBasisRational )
      {
         for( int i = 0; i < _slackCols.num(); i++ )
         {
            int col = numColsRational() - _slackCols.num() + i;
            int row = _slackCols.colVector(i).index(0);

            assert(row >= 0);
            assert(row < numRowsRational());

            switch( _basisStatusRowsRational[row] )
            {
            case SPxSolver::ON_LOWER:
               _basisStatusColsRational.append(SPxSolver::ON_UPPER);
               break;
            case SPxSolver::ON_UPPER:
               _basisStatusColsRational.append(SPxSolver::ON_LOWER);
               break;
            case SPxSolver::BASIC:
            case SPxSolver::FIXED:
            default:
               _basisStatusColsRational.append(_basisStatusRowsRational[row]);
               break;
            }

            _basisStatusRowsRational[row] = SPxSolver::FIXED;
         }
      }
   }



   /// restores original problem
   void SoPlex2::_untransformEqualityRational()
   {
      int numOrigCols = numColsRational() - _slackCols.num();

      // adjust solution
      if( _solRational.hasPrimal() )
      {
         for( int i = 0; i < _slackCols.num(); i++ )
         {
            int col = numOrigCols + i;
            int row = _slackCols.colVector(i).index(0);

            assert(row >= 0);
            assert(row < numRowsRational());

            _solRational._slacks[row] -= _solRational._primal[col];
         }

         _solRational._primal.reDim(numOrigCols);
      }

      if( _solRational.hasPrimalray() )
      {
         _solRational._primalray.reDim(numOrigCols);
      }

      if( _solRational.hasDual() )
      {
         _solRational._redcost.reDim(numOrigCols);
      }

      // adjust basis
      if( _hasBasisRational )
      {
         for( int i = 0; i < _slackCols.num(); i++ )
         {
            int col = numOrigCols + i;
            int row = _slackCols.colVector(i).index(0);

            assert(row >= 0);
            assert(row < numRowsRational());
            assert(_basisStatusRowsRational[row] == SPxSolver::FIXED || _basisStatusRowsRational[row] == SPxSolver::BASIC);

            if( _basisStatusRowsRational[row] == SPxSolver::FIXED )
            {
               switch( _basisStatusColsRational[col] )
               {
               case SPxSolver::ON_LOWER:
                  _basisStatusRowsRational[row] = SPxSolver::ON_UPPER;
                  break;
               case SPxSolver::ON_UPPER:
                  _basisStatusRowsRational[row] = SPxSolver::ON_LOWER;
                  break;
               case SPxSolver::BASIC:
               case SPxSolver::FIXED:
               default:
                  _basisStatusRowsRational[row] = _basisStatusColsRational[col];
                  break;
               }
            }
         }

         _basisStatusColsRational.reSize(numOrigCols);
      }

      // restore sides and remove slack columns
      for( int i = 0; i < _slackCols.num(); i++ )
      {
         int col = numOrigCols + i;
         int row = _slackCols.colVector(i).index(0);

         _rationalLP->changeLhs(row, -upperRational(col));
         _rationalLP->changeRhs(row, -lowerRational(col));
      }

      _rationalLP->removeColRange(numOrigCols, numColsRational() - 1);
   }



   /// solves current problem with iterative refinement and recovery mechanism
   void SoPlex2::_performOptIRStable(bool& primalFeasible, bool& dualFeasible, bool& infeasible, bool& unbounded, bool& stopped, bool& error)
   {
      primalFeasible = false;
      dualFeasible = false;
      infeasible = false;
      unbounded = false;
      stopped = false;
      error = false;
   }



   /// performs iterative refinement on the auxiliary problem for testing unboundedness
   void SoPlex2::_performUnboundedIRStable(bool& hasUnboundedRay, bool& stopped, bool& error)
   {
      hasUnboundedRay = false;
      stopped = false;
      error = false;
   }



   /// performs iterative refinement on the auxiliary problem for testing feasibility
   void SoPlex2::_performFeasIRStable(bool& infeasible, bool& stopped, bool& error)
   {
      infeasible = false;
      stopped = false;
      error = false;
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
