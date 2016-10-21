/*                                                                           */
/*                  This file is part of the class library                   */
/*       SoPlex --- the Sequential object-oriented simPlex.                  */
/*                                                                           */
/*    Copyright (C) 1996-2014 Konrad-Zuse-Zentrum                            */
/*                            fuer Informationstechnik Berlin                */
/*                                                                           */
/*  SoPlex is distributed under the terms of the ZIB Academic Licence.       */
/*                                                                           */
/*  You should have received a copy of the ZIB Academic License              */
/*  along with SoPlex; see the file COPYING. If not email to soplex@zib.de.  */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

#ifndef SOPLEX_LEGACY
#include <iostream>
#include <assert.h>

#include "soplex.h"
#include "statistics.h"
#include "sorter.h"

//#define NO_TOL
#define USE_FEASTOL
//#define NO_TRANSFORM
//#define PERFORM_COMPPROB_CHECK

#define MAX_DEGENCHECK     20    /**< the maximum number of degen checks that are performed before the IDS is abandoned */
#define DEGENCHECK_OFFSET  50    /**< the number of iteration before the degeneracy check is reperformed */
#define SLACKCOEFF         1.0   /**< the coefficient of the slack variable in the incompatible rows. */
#define TIMELIMIT_FRAC     0.5   /**< the fraction of the total time limit given to the setup of the reduced problem */

/* This file contains the private functions for the Improved Dual Simplex (IDS)
 *
 * An important note about the implementation of the IDS is the reliance on the row representation of the basis matrix.
 * The forming of the reduced and complementary problems is involves identifying rows in the row-form of the basis
 * matrix that have zero dual multipliers.
 *
 * Ideally, this work will be extended such that the IDS can be executed using the column-form of the basis matrix. */

namespace soplex
{




   /// solves LP using the improved dual simplex
   void SoPlex::_solveImprovedDualSimplex()
   {
      assert(_solver.rep() == SPxSolver::ROW);
      assert(_solver.type() == SPxSolver::LEAVE);

      bool stop = false;   // flag to indicate that the algorithm must terminate

      // setting the initial status of the reduced problem
      _statistics->redProbStatus = SPxSolver::NO_PROBLEM;

      // start timing
      _statistics->solvingTime->start();

      //@todo Need to check this. Where will the solving of the IDS take place? In the SoPlex::solve, _solvereal or
      //_solverational?
      // remember that last solve was in floating-point
      _lastSolveMode = SOLVEMODE_REAL;

      // setting the current solving mode.
      _currentProb = IDS_ORIG;

      // setting the sense to maximise
      if( intParam(SoPlex::OBJSENSE) == SoPlex::OBJSENSE_MINIMIZE )
      {
         assert(_solver.spxSense() == SPxLPBase<Real>::MINIMIZE);

         _solver.changeObj(_solver.maxObj());
         _solver.changeSense(SPxLPBase<Real>::MAXIMIZE);
      }

      // it is necessary to solve the initial problem to find a starting basis
      _solver.setIdsStatus(SPxSolver::FINDSTARTBASIS);

      // setting the decomposition iteration limit to the parameter setting
      _solver.setDecompIterationLimit(intParam(SoPlex::DECOMP_ITERLIMIT));

      int numDegenCheck = 0;
      Real degeneracyLevel = 0;
      _idsFeasVector.reDim(_solver.nCols());


      bool initSolveFromScratch = true;

      // arrays to store the basis status for the rows and columns at the interruption of the original problem solve.
      DataArray< SPxSolver::VarStatus > basisStatusRows;
      DataArray< SPxSolver::VarStatus > basisStatusCols;
      // since the original LP may have been shifted, the dual multiplier will not be correct for the original LP. This
      // loop will recheck the degeneracy level and compute the proper dual multipliers.
      do
      {
         _idsSimplifyAndSolve(_solver, _slufactor, initSolveFromScratch, initSolveFromScratch);
         initSolveFromScratch = false;
         _solver.basis().solve(_idsFeasVector, _solver.maxObj());
         degeneracyLevel = _solver.getDegeneracyLevel(_idsFeasVector);
         if( _solver.type() == SPxSolver::LEAVE || _solver.status() >= SPxSolver::OPTIMAL
               || _solver.status() == SPxSolver::ABORT_EXDECOMP || numDegenCheck > MAX_DEGENCHECK )
         {
            // returning the sense to minimise
            if( intParam(SoPlex::OBJSENSE) == SoPlex::OBJSENSE_MINIMIZE )
            {
               assert(_solver.spxSense() == SPxLPBase<Real>::MAXIMIZE);

               _solver.changeObj(-(_solver.maxObj()));
               _solver.changeSense(SPxLPBase<Real>::MINIMIZE);
            }

            _solver.setIdsStatus(SPxSolver::DONTFINDSTARTBASIS);

            if( _solver.status() == SPxSolver::UNBOUNDED || _solver.status() == SPxSolver::INFEASIBLE )
               _hasBasis = false;


            // resolving the problem to update the real lp and solve with the correct objective.
            // TODO: With some infeasible problem (e.g. refinery) the dual is violated. Setting fromScratch to true
            // avoids this issue. Need to check what the problem is.
            _idsSimplifyAndSolve(_solver, _slufactor, true, false);

            // retreiving the original problem statistics prior to destroying it.
            getOriginalProblemStatistics();

            // storing the solution from the reduced problem
            _storeSolutionReal();

            // stop timing
            _statistics->solvingTime->stop();
            return;
         }
         else if( _solver.status() == SPxSolver::ABORT_TIME || _solver.status() == SPxSolver::ABORT_ITER
            || _solver.status() == SPxSolver::ABORT_VALUE )
         {
            // at this point, the _idsSimplifyAndSolve does not store the realLP. It stores the _idsLP. As such, it is
            // important to reinstall the _realLP to the _solver.
            if( !_isRealLPLoaded )
            {
               _solver.loadLP(*_idsLP);
               spx_free(_idsLP);
               _idsLP = &_solver;
               _isRealLPLoaded = true;
            }

            // returning the sense to minimise
            if( intParam(SoPlex::OBJSENSE) == SoPlex::OBJSENSE_MINIMIZE )
            {
               assert(_solver.spxSense() == SPxLPBase<Real>::MAXIMIZE);

               _solver.changeObj(-(_solver.maxObj()));
               _solver.changeSense(SPxLPBase<Real>::MINIMIZE);
            }

            _preprocessAndSolveReal(false);

            // storing the solution from the reduced problem
            _storeSolutionReal();

            // stop timing
            _statistics->solvingTime->stop();
            return;
         }

         _solver.setDegenCompOffset(DEGENCHECK_OFFSET);

         numDegenCheck++;
      } while( (degeneracyLevel > 0.9 || degeneracyLevel < 0.1) || !checkBasisDualFeasibility(_idsFeasVector) );

      // decomposition simplex will commence, so the start basis does not need to be found
      _solver.setIdsStatus(SPxSolver::DONTFINDSTARTBASIS);

      // if the original problem has a basis, this will be stored
      if( _hasBasis )
      {
         basisStatusRows.reSize(numRowsReal());
         basisStatusCols.reSize(numColsReal());
         _solver.getBasis(basisStatusRows.get_ptr(), basisStatusCols.get_ptr(), basisStatusRows.size(),
            basisStatusCols.size());
      }

      // updating the algorithm iterations statistic
      _statistics->callsReducedProb++;

      // setting the statistic information
      numIdsIter = 0;
      numRedProbIter = _solver.iterations();
      numCompProbIter = 0;

      // setting the display information
      _idsDisplayLine = 0;

      MSG_INFO1( spxout,
         spxout << "========      Degeneracy Detected       ========" << std::endl
         << std::endl
         << "======== Commencing decomposition solve ========" << std::endl
         );

      MSG_INFO3(spxout, spxout << "Creating the Reduced and Complementary problems." << std::endl );

      // creating copies of the original problem that will be manipulated to form the reduced and complementary
      // problems.
      _createIdsReducedAndComplementaryProblems();

      MSG_INFO3(spxout, spxout << "Forming the Reduced problem." << std::endl );

      // creating the initial reduced problem from the basis information
      _formIdsReducedProblem(stop);

      // setting flags for the decomposition solve
      _hasBasis = false;
      bool hasRedBasis = false;
      bool redProbError = false;
      int algIterCount = 0;


      // a stop will be triggered if the reduced problem exceeded a specified fraction of the time limit
      if( stop )
      {
         MSG_INFO1( spxout, spxout << "==== Error constructing the reduced problem ====" << std::endl );
         redProbError = true;
         _statistics->redProbStatus = SPxSolver::NOT_INIT;
      }

      // setting the verbosity level
      const SPxOut::Verbosity orig_verbosity = spxout.getVerbosity();
      spxout.setVerbosity( SPxOut::ERROR );

      // the main solving loop of the decomposition simplex.
      // This loop solves the Reduced problem, and if the problem is feasible, the complementary problem is solved.
      while( !stop )
      {
         // setting the current solving mode.
         _currentProb = IDS_RED;

         // solve the reduced problem
         MSG_INFO2(spxout,
            spxout << std::endl
            << "=========================" << std::endl
            << "Solving: Reduced Problem." << std::endl
            << "=========================" << std::endl
            << std::endl );

         _hasBasis = hasRedBasis;
         // solving the reduced problem
         _idsSimplifyAndSolve(_solver, _slufactor, !algIterCount, !algIterCount);

         stop = idsTerminate(realParam(SoPlex::TIMELIMIT));  // checking whether the algorithm should terminate

         // updating the algorithm iterations statistics
         _statistics->callsReducedProb++;
         _statistics->redProbStatus = _solver.status();

         assert(_isRealLPLoaded);
         hasRedBasis = _hasBasis;
         _hasBasis = false;

         // checking the status of the reduced problem
         // It is expected that infeasibility and unboundedness will be found in the initialisation solve. If this is
         // not the case and the reduced problem is infeasible or unbounded the decomposition simplex will terminate.
         // If the decomposition simplex terminates, then the original problem will be solved using the stored basis.
         if( _solver.status() != SPxSolver::OPTIMAL )
         {
            if( _solver.status() == SPxSolver::UNBOUNDED )
               MSG_INFO2(spxout, spxout << "Unbounded reduced problem." << std::endl );
            if( _solver.status() == SPxSolver::INFEASIBLE )
               MSG_INFO2(spxout, spxout << "Infeasible reduced problem." << std::endl );

            MSG_INFO2(spxout, spxout << "Reduced problem status: " << _solver.status() << std::endl );

            redProbError = true;
            break;
         }

         // printing display line
         printIdsDisplayLine(_solver, orig_verbosity, !algIterCount, !algIterCount);

         // get the dual solutions from the reduced problem
         DVector reducedLPDualVector(_solver.nRows());
         DVector reducedLPRedcostVector(_solver.nCols());
         _solver.getDual(reducedLPDualVector);
         _solver.getRedCost(reducedLPRedcostVector);


         // Using the solution to the reduced problem:
         // In the first iteration of the algorithm create the complementary problem.
         // In the subsequent iterations of the algorithm update the complementary problem.
         if( algIterCount == 0 )
            _formIdsComplementaryProblem();
         else
         {
            if( boolParam(SoPlex::USECOMPDUAL) )
               _updateIdsComplementaryDualProblem(false);
            else
               _updateIdsComplementaryPrimalProblem(false);
         }

         // setting the current solving mode.
         _currentProb = IDS_COMP;

         // solve the complementary problem
         MSG_INFO2(spxout,
            spxout << std::endl
            << "=========================" << std::endl
            << "Solving: Complementary Problem." << std::endl
            << "=========================" << std::endl
            << std::endl );

         _idsSimplifyAndSolve(_compSolver, _compSlufactor, true, true);
         assert(_isRealLPLoaded);

         MSG_INFO3(spxout, spxout << "Iteration " << algIterCount
            << "Objective Value: " << std::setprecision(10) << _compSolver.objValue()
            << std::endl );


         // Check whether the complementary problem is solved with a non-negative objective function, is infeasible or
         // unbounded. If true, then stop the algorithm.
         if( GE(_compSolver.objValue(), 0.0, 1e-20)
            || _compSolver.status() == SPxSolver::INFEASIBLE
            || _compSolver.status() == SPxSolver::UNBOUNDED )
         {
            _statistics->compProbStatus = _compSolver.status();
            _statistics->finalCompObj = _compSolver.objValue();
            if( _compSolver.status() == SPxSolver::UNBOUNDED )
               MSG_INFO2(spxout, spxout << "Unbounded complementary problem." << std::endl );
            if( _compSolver.status() == SPxSolver::INFEASIBLE )
               MSG_INFO2(spxout, spxout << "Infeasible complementary problem." << std::endl );
            stop = true;
         }

         if( !stop )
         {
            if( !boolParam(SoPlex::EXPLICITVIOL) )
            {
               // get the primal solutions from the complementary problem
               DVector compLPPrimalVector(_compSolver.nCols());
               _compSolver.getPrimal(compLPPrimalVector);

               // get the dual solutions from the complementary problem
               DVector compLPDualVector(_compSolver.nRows());
               _compSolver.getDual(compLPDualVector);

               // updating the reduced problem
               _updateIdsReducedProblem(_compSolver.objValue(), reducedLPDualVector, reducedLPRedcostVector,
                  compLPPrimalVector, compLPDualVector);
            }
            else
            {
               // getting the primal solutions from the reduced problem
               DVector reducedLPPrimalVector(_solver.nCols());
               _solver.getPrimal(reducedLPPrimalVector);

               // checking the optimality of the reduced problem solution with the original problem
               _checkOriginalProblemOptimality(reducedLPPrimalVector, false);

               // updating the reduced problem based upon the violated rows and constraints found during the optimality
               // check
               _updateIdsReducedProblemViol(false);
            }
         }
         // if the complementary problem is infeasible or unbounded, it is possible that the algorithm can continue.
         // a check of the original problem is required to determine whether there are any violations.
         else if( _compSolver.status() == SPxSolver::INFEASIBLE
            || _compSolver.status() == SPxSolver::UNBOUNDED )
         {
            // getting the primal vector from the reduced problem
            DVector reducedLPPrimalVector(_solver.nCols());
            _solver.getPrimal(reducedLPPrimalVector);

            // checking the optimality of the reduced problem solution with the original problem
            _checkOriginalProblemOptimality(reducedLPPrimalVector, false);

            // if there are any violated rows or bounds then stop is reset and the algorithm continues.
            if( _nIdsViolBounds > 0 || _nIdsViolRows > 0 )
               stop = false;
            if( _nIdsViolBounds == 0 && _nIdsViolRows == 0 )
               stop = true;

            // updating the reduced problem with the original problem violated rows
            if( !stop )
               _updateIdsReducedProblemViol(true);
         }

         numIdsIter++;
         algIterCount++;
      }

      // if there is an error in solving the reduced problem, i.e. infeasible or unbounded, then we leave the
      // decomposition solve and resolve the original. Infeasibility should be dealt with in the original problem.
      if( !redProbError )
      {
#ifndef NDEBUG
         // computing the solution for the original variables
         DVector reducedLPPrimalVector(_solver.nCols());
         _solver.getPrimal(reducedLPPrimalVector);

         // checking the optimality of the reduced problem solution with the original problem
         _checkOriginalProblemOptimality(reducedLPPrimalVector, true);

         // resetting the verbosity level
         spxout.setVerbosity( orig_verbosity );
#endif

         // This is an additional check to ensure that the complementary problem is solving correctly.
#ifdef PERFORM_COMPPROB_CHECK
         // solving the complementary problem with the original objective function
         if( boolParam(SoPlex::USECOMPDUAL) )
            _setComplementaryDualOriginalObjective();
         else
            _setComplementaryPrimalOriginalObjective();

         // for clean up
         // the real lp has to be set to not loaded.
         //_isRealLPLoaded = false;
         // the basis has to be set to false
         _hasBasis = false;

         if( boolParam(SoPlex::USECOMPDUAL) )
         {
            SPxLPIds compDualLP;
            _compSolver.buildDualProblem(compDualLP, _idsPrimalRowIDs.get_ptr(), _idsPrimalColIDs.get_ptr(),
                  _idsDualRowIDs.get_ptr(), _idsDualColIDs.get_ptr(), &_nPrimalRows, &_nPrimalCols, &_nDualRows, &_nDualCols);

            _compSolver.loadLP(compDualLP);
         }

         _idsSimplifyAndSolve(_compSolver, _compSlufactor, true, true);

         _solReal._hasPrimal = true;
         _hasSolReal = true;
         // get the primal solutions from the reduced problem
         DVector testPrimalVector(_compSolver.nCols());
         _compSolver.getPrimal(testPrimalVector);
         _solReal._primal.reDim(_compSolver.nCols());
         _solReal._primal = testPrimalVector;

         Real maxviol = 0;
         Real sumviol = 0;

         if( getIdsBoundViolation(maxviol, sumviol) )
            MSG_INFO1(spxout, spxout << "Bound violation - "
               << "Max: "<< std::setprecision(20) << maxviol << " "
               << "Sum: "<< sumviol << std::endl );

         if( getIdsRowViolation(maxviol, sumviol) )
            MSG_INFO1(spxout, spxout << "Row violation - "
               << "Max: "<< std::setprecision(21) << maxviol << " "
               << "Sum: "<< sumviol << std::endl );

         MSG_INFO1(spxout, spxout << "Objective Value: " << _compSolver.objValue() << std::endl );
#endif
      }

      // resetting the verbosity level
      spxout.setVerbosity( orig_verbosity );

      MSG_INFO1( spxout,
         spxout << "========  Decomposition solve completed ========" << std::endl
         << std::endl
         << "========   Resolving original problem   ========" << std::endl
         );

      // if there is a reduced problem error in the first iteration the complementary problme has not been
      // set up. In this case no memory has been allocated for _idsCompProbColIDsIdx and _fixedOrigVars.
      if( !redProbError || algIterCount > 0 )
      {
         spx_free(_idsCompProbColIDsIdx);
         spx_free(_fixedOrigVars);
      }
      spx_free(_idsViolatedRows);
      spx_free(_idsViolatedBounds);
      spx_free(_idsReducedProbCols);
      spx_free(_idsReducedProbRows);

      // retreiving the original problem statistics prior to destroying it.
      getOriginalProblemStatistics();

      // setting the reduced problem statistics
      _statistics->numRedProbRows = numIncludedRows;
      _statistics->numRedProbCols = _solver.nCols();


      // printing display line for resolve of problem
      _solver.printDisplayLine(false, true);

      // if there is an error solving the reduced problem the LP must be solved from scratch
      if (redProbError)
      {
         MSG_INFO1( spxout, spxout << "=== Reduced problem error - Solving original ===" << std::endl );

         // the solver is loaded with the realLP to solve the problem from scratch
         _solver.loadLP(*_realLP);
         spx_free(_realLP);
         _realLP = &_solver;
         _isRealLPLoaded = true;

         // returning the sense to minimise
         if( intParam(SoPlex::OBJSENSE) == SoPlex::OBJSENSE_MINIMIZE )
         {
            assert(_solver.spxSense() == SPxLPBase<Real>::MAXIMIZE);

            _solver.changeObj(-(_solver.maxObj()));
            _solver.changeSense(SPxLPBase<Real>::MINIMIZE);

            // Need to add commands to multiply the objective solution values by -1
         }

         //_solver.setBasis(basisStatusRows.get_const_ptr(), basisStatusCols.get_const_ptr());
         _preprocessAndSolveReal(true);
      }
      else
      {
         // resolving the problem to update the real lp and solve with the correct objective.
         // the realLP is updated with the current solver. This is to resolve the problem to get the correct solution
         // values
         _realLP->~SPxLPReal();
         spx_free(_realLP);
         _realLP = &_solver;
         _isRealLPLoaded = true;

         // returning the sense to minimise
         if( intParam(SoPlex::OBJSENSE) == SoPlex::OBJSENSE_MINIMIZE )
         {
            assert(_solver.spxSense() == SPxLPBase<Real>::MAXIMIZE);

            _solver.changeObj(-(_solver.maxObj()));
            _solver.changeSense(SPxLPBase<Real>::MINIMIZE);

            // Need to add commands to multiply the objective solution values by -1
         }

         _preprocessAndSolveReal(false);
      }

      // storing the solution from the reduced problem
      _storeSolutionReal();

      // stop timing
      _statistics->solvingTime->stop();
   }



   /// creating copies of the original problem that will be manipulated to form the reduced and complementary problems
   void SoPlex::_createIdsReducedAndComplementaryProblems()
   {
      // the reduced problem is formed from the current problem
      // So, we copy the _solver to the _realLP and work on the _solver
      // NOTE: there is no need to preprocess because we always have a starting basis.
      _realLP = 0;
      spx_alloc(_realLP);
      _realLP = new (_realLP) SPxLPIds(_solver);

      // allocating memory for the reduced problem rows and cols flag array
      _idsReducedProbRows = 0;
      spx_alloc(_idsReducedProbRows, numRowsReal());
      _idsReducedProbCols = 0;
      spx_alloc(_idsReducedProbCols, numColsReal());

      // the complementary problem is formulated with all incompatible rows and those from the reduced problem that have
      // a positive reduced cost.
      _compSolver = _solver;
      _compSolver.setOutstream(spxout);
      _compSolver.setSolver(&_compSlufactor);

      // allocating memory for the violated bounds and rows arrays
      _idsViolatedBounds = 0;
      _idsViolatedRows = 0;
      spx_alloc(_idsViolatedBounds, numColsReal());
      spx_alloc(_idsViolatedRows, numRowsReal());
      _nIdsViolBounds = 0;
      _nIdsViolRows = 0;
   }



   /// forms the reduced problem
   void SoPlex::_formIdsReducedProblem(bool& stop)
   {
      MSG_INFO3( spxout, spxout << "Forming the reduced problem" << std::endl );
      int* nonposind = 0;
      int* compatind = 0;
      int* rowsforremoval = 0;
      int* colsforremoval = 0;
      int nnonposind = 0;
      int ncompatind = 0;

      assert(_solver.nCols() == numColsReal());
      assert(_solver.nRows() == numRowsReal());

      // capturing the basis used for the transformation
      _idsTransBasis = _solver.basis();

      // setting row counter to zero
      numIncludedRows = 0;

      // the _idsLP is used as a helper LP object. This is used so that _realLP can still hold the original problem.
      _idsLP = 0;
      spx_alloc(_idsLP);
      _idsLP = new (_idsLP) SPxLPIds(_solver);

      // retreiving the basis information
      _basisStatusRows.reSize(numRowsReal());
      _basisStatusCols.reSize(numColsReal());
      _solver.getBasis(_basisStatusRows.get_ptr(), _basisStatusCols.get_ptr());

      // get the indices of the rows with positive dual multipliers and columns with positive reduced costs.
      spx_alloc(nonposind, numColsReal());
      spx_alloc(colsforremoval, numColsReal());
      if( !stop )
         _getZeroDualMultiplierIndices(_idsFeasVector, nonposind, colsforremoval, &nnonposind, stop);

      // get the compatible columns from the constraint matrix w.r.t the current basis matrix
      MSG_INFO3(spxout, spxout << "Computing the compatible columns" << std::endl
        << "Solving time: " << solveTime() << std::endl );

      spx_alloc(compatind, _solver.nRows());
      spx_alloc(rowsforremoval, _solver.nRows());
      if( !stop )
         _getCompatibleColumns(_idsFeasVector, nonposind, compatind, rowsforremoval, colsforremoval, nnonposind,
            &ncompatind, true, stop);

      int* compatboundcons = 0;
      int ncompatboundcons = 0;
      spx_alloc(compatboundcons, numColsReal());

      LPRowSet boundcons;

      // identifying the compatible bound constraints
      if( !stop )
         _getCompatibleBoundCons(boundcons, compatboundcons, nonposind, &ncompatboundcons, nnonposind, stop);

      // delete rows and columns from the LP to form the reduced problem
      MSG_INFO3(spxout, spxout << "Deleting rows and columns to form the reduced problem" << std::endl
         << "Solving time: " << solveTime() << std::endl );

      // allocating memory to add bound constraints
      SPxRowId* addedrowids = 0;
      spx_alloc(addedrowids, ncompatboundcons);

      // computing the reduced problem obj coefficient vector and updating the problem
      if( !stop )
      {
         _computeReducedProbObjCoeff(stop);

         _solver.loadLP(*_idsLP);

         // the colsforremoval are the columns with a zero reduced cost.
         // the rowsforremoval are the rows identified as incompatible.
         _solver.removeRows(rowsforremoval);
         // adding the rows for the compatible bound constraints
         _solver.addRows(addedrowids, boundcons);

         for( int i = 0; i < ncompatboundcons; i++ )
            _idsReducedProbColRowIDs[compatboundcons[i]] = addedrowids[i];
      }


      // freeing allocated memory
      spx_free(addedrowids);
      spx_free(compatboundcons);
      spx_free(rowsforremoval);
      spx_free(compatind);
      spx_free(colsforremoval);
      spx_free(nonposind);

      _idsLP->~SPxLPIds();
      spx_free(_idsLP);
   }



   /// forms the complementary problem
   void SoPlex::_formIdsComplementaryProblem()
   {
      // delete rows and columns from the LP to form the reduced problem
      _nElimPrimalRows = 0;
      _idsElimPrimalRowIDs.reSize(_realLP->nRows());   // the number of eliminated rows is less than the number of rows
                                                      // in the reduced problem
      _idsCompProbColIDsIdx = 0;
      spx_alloc(_idsCompProbColIDsIdx, _realLP->nCols());

      _fixedOrigVars = 0;
      spx_alloc(_fixedOrigVars, _realLP->nCols());
      for (int i = 0; i < _realLP->nCols(); i++)
      {
         _idsCompProbColIDsIdx[i] = -1;
         _fixedOrigVars[i] = -2; // this must be set to a value differet from -1, 0, 1 to ensure the
                                 // _updateComplementaryFixedPrimalVars function updates all variables in the problem.
      }

      int naddedrows = 0;
      DataArray< SPxRowId > rangedRowIds(numRowsReal());
      _deleteAndUpdateRowsComplementaryProblem(rangedRowIds.get_ptr(), naddedrows);

      if( boolParam(SoPlex::USECOMPDUAL) )   // if we use the dual formulation of the complementary problem, we must
                                             // perform the transformation
      {
         // initialising the arrays to store the row id's from the primal and the col id's from the dual
         _idsPrimalRowIDs.reSize(3*_compSolver.nRows());
         _idsPrimalColIDs.reSize(3*_compSolver.nCols());
         _idsDualRowIDs.reSize(3*_compSolver.nCols());
         _idsDualColIDs.reSize(3*_compSolver.nRows());
         _nPrimalRows = 0;
         _nPrimalCols = 0;
         _nDualRows = 0;
         _nDualCols = 0;

         // convert complementary problem to dual problem
         SPxLPIds compDualLP;
         _compSolver.buildDualProblem(compDualLP, _idsPrimalRowIDs.get_ptr(), _idsPrimalColIDs.get_ptr(),
               _idsDualRowIDs.get_ptr(), _idsDualColIDs.get_ptr(), &_nPrimalRows, &_nPrimalCols, &_nDualRows, &_nDualCols);

         // setting the id index array for the complementary problem
         for( int i = 0; i < _nPrimalRows; i++ )
         {
            // a primal row may result in two dual columns. So the index array points to the first of the indicies for the
            // primal row.
            if( i + 1 < _nPrimalRows &&
                  _realLP->number(SPxRowId(_idsPrimalRowIDs[i])) == _realLP->number(SPxRowId(_idsPrimalRowIDs[i + 1])) )
               i++;
         }

         for( int i = 0; i < _nPrimalCols; i++ )
         {
            if( _idsPrimalColIDs[i].getIdx() != _compSlackColId.getIdx() )
               _idsCompProbColIDsIdx[_realLP->number(_idsPrimalColIDs[i])] = i;
         }

         // retrieving the dual row id for the complementary slack column
         // there should be a one to one relationship between the number of primal columns and the number of dual rows.
         // hence, it should be possible to equate the dual row id to the related primal column.
         assert(_nPrimalCols == _nDualRows);
         assert(_compSolver.nCols() == compDualLP.nRows());
         _compSlackDualRowId = compDualLP.rId(_compSolver.number(_compSlackColId));

         _compSolver.loadLP(compDualLP);

         _compSolver.init();

         _updateComplementaryDualSlackColCoeff();

         // initalising the array for the dual columns representing the fixed variables in the complementary problem.
         _idsFixedVarDualIDs.reSize(_realLP->nCols());
         _idsVarBoundDualIDs.reSize(2*_realLP->nCols());


         // updating the constraints for the complementary problem
         _updateIdsComplementaryDualProblem(false);
      }
      else  // when using the primal of the complementary problem, no transformation of the problem is required.
      {
         // initialising the arrays to store the row and column id's from the original problem the col id's from the dual
         // if a row or column is not included in the complementary problem, the ID is set to INVALID in the
         // _idsPrimalRowIDs or _idsPrimalColIDs respectively.
         _idsPrimalRowIDs.reSize(_compSolver.nRows()*2);
         _idsPrimalColIDs.reSize(_compSolver.nCols());
         _idsCompPrimalRowIDs.reSize(_compSolver.nRows()*2);
         _idsCompPrimalColIDs.reSize(_compSolver.nCols());
         _nPrimalRows = 0;
         _nPrimalCols = 0;

         // at this point in time the complementary problem contains all rows and columns from the original problem
         int addedrangedrows = 0;
         _nCompPrimalRows = 0;
         _nCompPrimalCols = 0;
         for( int i = 0; i < _realLP->nRows(); i++ )
         {
            assert(_nCompPrimalRows < _compSolver.nRows());
            _idsPrimalRowIDs[_nPrimalRows] = _realLP->rId(i);
            _idsCompPrimalRowIDs[_nCompPrimalRows] = _compSolver.rId(i);
            _nPrimalRows++;
            _nCompPrimalRows++;

            if( _realLP->rowType(i) == LPRowBase<Real>::RANGE || _realLP->rowType(i) == LPRowBase<Real>::EQUAL )
            {
               assert(EQ(_compSolver.lhs(rangedRowIds[addedrangedrows]), _realLP->lhs(i)));
               assert(LE(_compSolver.rhs(rangedRowIds[addedrangedrows]), infinity));
               assert(LE(_compSolver.lhs(i), -infinity));
               assert(LT(_compSolver.rhs(i), infinity));

               _idsPrimalRowIDs[_nPrimalRows] = _realLP->rId(i);
               _idsCompPrimalRowIDs[_nCompPrimalRows] = rangedRowIds[addedrangedrows];
               _nPrimalRows++;
               _nCompPrimalRows++;
               addedrangedrows++;
            }
         }

         for( int i = 0; i < _compSolver.nCols(); i++ )
         {

            if( _compSolver.cId(i).getIdx() != _compSlackColId.getIdx() )
            {
               _idsPrimalColIDs[i] = _realLP->cId(i);
               _idsCompPrimalColIDs[i] = _compSolver.cId(i);
               _nPrimalCols++;
               _nCompPrimalCols++;

               _idsCompProbColIDsIdx[_realLP->number(_idsPrimalColIDs[i])] = i;
            }
         }

         // updating the constraints for the complementary problem
         _updateIdsComplementaryPrimalProblem(false);
      }
   }



   /// simplifies the problem and solves
   void SoPlex::_idsSimplifyAndSolve(SPxSolver& solver, SLUFactor& sluFactor, bool fromScratch, bool applyPreprocessing)
   {
      if( realParam(SoPlex::TIMELIMIT) < realParam(SoPlex::INFTY) )
         solver.setTerminationTime(realParam(SoPlex::TIMELIMIT) - _statistics->solvingTime->time());

      solver.changeObjOffset(realParam(SoPlex::OBJ_OFFSET));
      _statistics->preprocessingTime->start();

      SPxSimplifier::Result result = SPxSimplifier::OKAY;

      /* delete starting basis if solving from scratch */
      if ( fromScratch )
      {
         try
         {
            solver.reLoad();
         }
         catch( SPxException E )
         {
            MSG_ERROR( spxout << "Caught exception <" << E.what() << "> during simplify and solve.\n" );
         }
      }
      //assert(!fromScratch || solver.status() == SPxSolver::NO_PROBLEM);

      if( applyPreprocessing )
      {
         _enableSimplifierAndScaler();
         solver.setTerminationValue(realParam(SoPlex::INFTY));
      }
      else
      {
         _disableSimplifierAndScaler();
         ///@todo implement for both objective senses
         solver.setTerminationValue(solver.spxSense() == SPxLPBase<Real>::MINIMIZE
               ? realParam(SoPlex::OBJLIMIT_UPPER) : realParam(SoPlex::INFTY));
      }

      /* store original lp */
      applyPreprocessing = (_scaler != NULL || _simplifier != NULL);

      if( _isRealLPLoaded )
      {
         //assert(_idsLP == &solver); // there should be an assert here, but I don't know what to use.

         // preprocessing is always applied to the LP in the solver; hence we have to create a copy of the original LP
         // if preprocessing is turned on
         if( applyPreprocessing )
         {
            _idsLP = 0;
            spx_alloc(_idsLP);
            _idsLP = new (_idsLP) SPxLPIds(solver);
            _isRealLPLoaded = false;
         }
         else
            _idsLP = &solver;
      }
      else
      {
         assert(_idsLP != &solver);

         // ensure that the solver has the original problem
         solver.loadLP(*_idsLP);

         // load basis if available
         if( _hasBasis )
         {
            assert(_basisStatusRows.size() == solver.nRows());
            assert(_basisStatusCols.size() == solver.nCols());

            ///@todo this should not fail even if the basis is invalid (wrong dimension or wrong number of basic
            ///      entries); fix either in SPxSolver or in SPxBasis
            solver.setBasis(_basisStatusRows.get_const_ptr(), _basisStatusCols.get_const_ptr());
         }

         // if there is no preprocessing, then the original and the transformed problem are identical and it is more
         // memory-efficient to keep only the problem in the solver
         if( !applyPreprocessing )
         {
            _idsLP->~SPxLPIds();
            spx_free(_idsLP);
            _idsLP = &solver;
            _isRealLPLoaded = true;
         }
         else
         {
            _idsLP = 0;
            spx_alloc(_idsLP);
            _idsLP = new (_idsLP) SPxLPIds(solver);
         }
      }

      // assert that we have two problems if and only if we apply preprocessing
      assert(_idsLP == &solver || applyPreprocessing);
      assert(_idsLP != &solver || !applyPreprocessing);

      // apply problem simplification
      if( _simplifier != 0 )
      {
         result = _simplifier->simplify(solver, realParam(SoPlex::EPSILON_ZERO), realParam(SoPlex::FEASTOL),
               realParam(SoPlex::OPTTOL));
         solver.changeObjOffset(_simplifier->getObjoffset() + realParam(SoPlex::OBJ_OFFSET));
      }

      _statistics->preprocessingTime->stop();

      // run the simplex method if problem has not been solved by the simplifier
      if( result == SPxSimplifier::OKAY )
      {
         if( _scaler != 0 )
            _scaler->scale(solver);

         bool _hadBasis = _hasBasis;

         _statistics->simplexTime->start();
         try
         {
            solver.solve();
         }
         catch( const SPxException& E )
         {
            MSG_ERROR( std::cerr << "Caught exception <" << E.what() << "> while solving real LP.\n" );
            _status = SPxSolver::ERROR;
         }
         catch( ... )
         {
            MSG_ERROR( std::cerr << "Caught unknown exception while solving real LP.\n" );
            _status = SPxSolver::ERROR;
         }
         _statistics->simplexTime->stop();

         // record statistics
         // only record the main statistics for the original problem and reduced problem.
         // the complementary problem is a pivoting rule, so these will be held in other statistics.
         // statistics on the number of iterations for each problem is stored in individual variables.
         if( _currentProb == IDS_ORIG || _currentProb == IDS_RED )
         {
            _statistics->iterations += solver.iterations();
            _statistics->iterationsPrimal += solver.primalIterations();
            _statistics->iterationsFromBasis += _hadBasis ? solver.iterations() : 0;
            _statistics->boundflips += solver.boundFlips();
            _statistics->luFactorizationTimeReal += sluFactor.getFactorTime();
            _statistics->luSolveTimeReal += sluFactor.getSolveTime();
            _statistics->luFactorizationsReal += sluFactor.getFactorCount();
            _statistics->luSolvesReal += sluFactor.getSolveCount();
            sluFactor.resetCounters();

            _statistics->degenPivotsPrimal += solver.primalDegeneratePivots();
            _statistics->degenPivotsDual += solver.dualDegeneratePivots();

            _statistics->sumDualDegen += _solver.sumDualDegeneracy();
            _statistics->sumPrimalDegen += _solver.sumPrimalDegeneracy();

            if( _currentProb == IDS_ORIG )
               _statistics->iterationsInit += solver.iterations();
            else
               _statistics->iterationsRedProb += solver.iterations();
         }

         if( _currentProb == IDS_COMP )
            _statistics->iterationsCompProb += solver.iterations();

      }

      // check the result and run again without preprocessing if necessary
      _evaluateSolutionIDS(solver, sluFactor, result);
   }



   /// loads original problem into solver and solves again after it has been solved to optimality with preprocessing
   void SoPlex::_idsResolveWithoutPreprocessing(SPxSolver& solver, SLUFactor& sluFactor, SPxSimplifier::Result result)
   {
      assert(_simplifier != 0 || _scaler != 0);
      assert(result == SPxSimplifier::VANISHED
         || (result == SPxSimplifier::OKAY
            && (solver.status() == SPxSolver::OPTIMAL
               || solver.status() == SPxSolver::ABORT_DECOMP
               || solver.status() == SPxSolver::ABORT_EXDECOMP )));

      // if simplifier is active and LP is solved in presolving or to optimality, then we unsimplify to get the basis
      if( _simplifier != 0 )
      {
         assert(!_simplifier->isUnsimplified());

         bool vanished = result == SPxSimplifier::VANISHED;

         // get solution vectors for transformed problem
         DVectorReal primal(vanished ? 0 : solver.nCols());
         DVectorReal slacks(vanished ? 0 : solver.nRows());
         DVectorReal dual(vanished ? 0 : solver.nRows());
         DVectorReal redCost(vanished ? 0 : solver.nCols());

         assert(!_isRealLPLoaded);
         _basisStatusRows.reSize(_idsLP->nRows());
         _basisStatusCols.reSize(_idsLP->nCols());
         assert(vanished || _basisStatusRows.size() >= solver.nRows());
         assert(vanished || _basisStatusCols.size() >= solver.nCols());

         if( !vanished )
         {
            assert(solver.status() == SPxSolver::OPTIMAL || solver.status() == SPxSolver::ABORT_DECOMP || solver.status() == SPxSolver::ABORT_EXDECOMP);

            solver.getPrimal(primal);
            solver.getSlacks(slacks);
            solver.getDual(dual);
            solver.getRedCost(redCost);

            // unscale vectors
            if( _scaler != 0 )
            {
               _scaler->unscalePrimal(primal);
               _scaler->unscaleSlacks(slacks);
               _scaler->unscaleDual(dual);
               _scaler->unscaleRedCost(redCost);
            }

            // get basis of transformed problem
            solver.getBasis(_basisStatusRows.get_ptr(), _basisStatusCols.get_ptr());
         }

         try
         {
            _simplifier->unsimplify(primal, dual, slacks, redCost, _basisStatusRows.get_ptr(), _basisStatusCols.get_ptr());
            _simplifier->getBasis(_basisStatusRows.get_ptr(), _basisStatusCols.get_ptr());
            _hasBasis = true;
         }
         catch( SPxException E )
         {
            MSG_ERROR( spxout << "Caught exception <" << E.what() << "> during unsimplification. Resolving without simplifier and scaler.\n" );
         }
         catch( ... )
         {
            MSG_ERROR( spxout << "Caught unknown exception during unsimplification. Resolving without simplifier and scaler.\n" );
            _status = SPxSolver::ERROR;
         }
      }
      // if the original problem is not in the solver because of scaling, we also need to store the basis
      else if( _scaler != 0 )
      {
         _basisStatusRows.reSize(numRowsReal());
         _basisStatusCols.reSize(numColsReal());
         assert(_basisStatusRows.size() == solver.nRows());
         assert(_basisStatusCols.size() == solver.nCols());

         solver.getBasis(_basisStatusRows.get_ptr(), _basisStatusCols.get_ptr());
         _hasBasis = true;
      }

      // resolve the original problem
      _idsSimplifyAndSolve(solver, sluFactor, false, false);
      return;
   }


   /// updates the reduced problem with additional rows using the solution to the complementary problem
   void SoPlex::_updateIdsReducedProblem(Real objValue, DVector dualVector, DVector redcostVector,
         DVector compPrimalVector, DVector compDualVector)
   {
      Real feastol = realParam(SoPlex::FEASTOL);

      Real maxDualRatio = infinity;

      bool usecompdual = boolParam(SoPlex::USECOMPDUAL);

      for( int i = 0; i < _nPrimalRows; i++ )
      {
         Real reducedProbDual = 0;
         Real compProbPrimal = 0;
         Real dualRatio = 0;
         int rownumber = _realLP->number(SPxRowId(_idsPrimalRowIDs[i]));
         if( _idsReducedProbRows[rownumber] && _solver.isBasic(_idsReducedProbRowIDs[rownumber]) )
         {
            int solverRowNum = _solver.number(_idsReducedProbRowIDs[rownumber]);

            // retreiving the reduced problem dual solutions and the complementary problem primal solutions
            reducedProbDual = dualVector[solverRowNum]; // this is y
            //if( _solver.rowType(solverRowNum) == LPRowBase<Real>::GREATER_EQUAL )
               //reducedProbDual *= -1;
            if( usecompdual )
               compProbPrimal = compPrimalVector[_compSolver.number(SPxColId(_idsDualColIDs[i]))]; // this is u
            else
               compProbPrimal = compDualVector[_compSolver.number(SPxRowId(_idsCompPrimalRowIDs[i]))];

            // the variable in the basis is degenerate.
            if( EQ(reducedProbDual, 0.0, feastol) )
            {
               MSG_WARNING( spxout,
                 spxout << "WIMDSM01: reduced problem dual value is very close to zero." << std::endl; );
               continue;
            }

            // the translation of the complementary primal problem to the dual some rows resulted in two columns.
            if( usecompdual && i < _nPrimalRows - 1 &&
                  _realLP->number(SPxRowId(_idsPrimalRowIDs[i])) == _realLP->number(SPxRowId(_idsPrimalRowIDs[i + 1])) )
            {
               i++;
               // @todo make sure that this is just a simple sum
               compProbPrimal += compPrimalVector[_compSolver.number(SPxColId(_idsDualColIDs[i]))]; // this is u
            }



            // updating the ratio
            SoPlex::DualSign varSign = getExpectedDualVariableSign(solverRowNum);
            if( varSign == SoPlex::IS_FREE || (varSign == SoPlex::IS_POS && LE(compProbPrimal, 0, feastol)) ||
                  (varSign == SoPlex::IS_NEG && GE(compProbPrimal, 0, feastol)) )
            {
               dualRatio = infinity;
            }
            else
            {
               dualRatio = reducedProbDual/compProbPrimal;
            }

            if( LT(dualRatio, maxDualRatio, feastol) )
               maxDualRatio = dualRatio;
         }
         else
         {
            if( usecompdual )
               compProbPrimal = compPrimalVector[_compSolver.number(SPxColId(_idsDualColIDs[i]))]; // this is v
            else
               compProbPrimal = compDualVector[_compSolver.number(SPxRowId(_idsCompPrimalRowIDs[i]))];
         }

      }

      // setting the maxDualRatio to a maximum of 1
      if( maxDualRatio > 1.0 )
         maxDualRatio = 1.0;

      DVector compProbRedcost(_compSolver.nCols());   // the reduced costs of the complementary problem

      // Retrieving the reduced costs for each original problem row.
      _compSolver.getRedCost(compProbRedcost);

      LPRowSet updaterows;

      DataArray<RowViolation> violatedrows;
      int nviolatedrows = 0;
      int* newrowidx = 0;
      int nnewrowidx = 0;
      spx_alloc(newrowidx, _nPrimalRows);

      violatedrows.reSize(_nPrimalRows);

      bool ratioTest = true;
      //ratioTest = false;
      for( int i = 0; i < _nPrimalRows; i++ )
      {
         LPRowReal origlprow;
         DSVectorBase<Real> rowtoaddVec(_realLP->nCols());
         Real compProbPrimal = 0;
         Real compRowRedcost = 0;
         int rownumber = _realLP->number(SPxRowId(_idsPrimalRowIDs[i]));
         if( !_idsReducedProbRows[_realLP->number(SPxRowId(_idsPrimalRowIDs[i]))] )
         {
            int compRowNumber;
            if( usecompdual )
            {
               compRowNumber = _compSolver.number(_idsDualColIDs[i]);
               // retreiving the complementary problem primal solutions
               compProbPrimal = compPrimalVector[compRowNumber]; // this is v
               compRowRedcost = compProbRedcost[compRowNumber];
            }
            else
            {
               compRowNumber = _compSolver.number(_idsCompPrimalRowIDs[i]);
               // retreiving the complementary problem primal solutions
               compProbPrimal = compDualVector[compRowNumber]; // this is v
            }

            // the translation of the complementary primal problem to the dual some rows resulted in two columns.
            if( usecompdual && i < _nPrimalRows - 1 &&
                  _realLP->number(SPxRowId(_idsPrimalRowIDs[i])) == _realLP->number(SPxRowId(_idsPrimalRowIDs[i + 1])) )
            {
               i++;
               compRowNumber = _compSolver.number(SPxColId(_idsDualColIDs[i]));
               // @todo make sure that this is just a simple sum
               compProbPrimal += compPrimalVector[compRowNumber]; // this is v
               compRowRedcost += compProbRedcost[compRowNumber];
            }

            SoPlex::DualSign varSign = getOrigProbDualVariableSign(rownumber);

            // add row to the reduced problem if the computed dual is of the correct sign for a feasible dual solution
            if( ratioTest && ((varSign == SoPlex::IS_FREE && !isZero(compProbPrimal, feastol)) ||
                  (varSign == SoPlex::IS_POS && GT(compProbPrimal*maxDualRatio, 0, feastol)) ||
                  (varSign == SoPlex::IS_NEG && LT(compProbPrimal*maxDualRatio, 0, feastol)))
                  /*|| (//isZero(objValue, 1e-10) &&
                   * _compSolver.basis().desc().rowStatus(compRowNum) > SPxBasis::Desc::D_FREE &&
                     _compSolver.basis().desc().rowStatus(compRowNum) < SPxBasis::Desc::D_UNDEFINED)*/ )
            {
               //this set of statements are required to add rows to the reduced problem. This will add a row for every
               //violated constraint. I only want to add a row for the most violated constraint. That is why the row
               //adding functionality is put outside the for loop.
               if( !_idsReducedProbRows[rownumber] )
               {
                  numIncludedRows++;
                  assert(numIncludedRows <= _realLP->nRows());
               }
               violatedrows[nviolatedrows].idx = rownumber;
               violatedrows[nviolatedrows].violation = spxAbs(compProbPrimal*maxDualRatio);
               nviolatedrows++;
            }
         }
      }


      // sorting the violated rows by the absolute violation.
      // Only a predefined number of rows will be added to the reduced problem
      RowViolationCompare compare;
      compare.entry = violatedrows.get_const_ptr();

      // if no rows are identified by the pricing rule, we add rows based upon the constraint violations
      if( !ratioTest || nviolatedrows == 0 )
      {
         _findViolatedRows(objValue, violatedrows, nviolatedrows);
      }


      int sorted = 0;
      int sortsize = MINIMUM(intParam(SoPlex::DECOMP_MAXADDEDROWS), nviolatedrows);

      // only sorting if the sort size is less than the number of violated rows.
      if( sortsize > 0 && sortsize < nviolatedrows )
         sorted = SPxQuicksortPart(violatedrows.get_ptr(), compare, sorted + 1, nviolatedrows, sortsize);

      // adding the violated rows.
      for( int i = 0; i < sortsize; i++ )
      {
         updaterows.add(_transformedRows.lhs(violatedrows[i].idx), _transformedRows.rowVector(violatedrows[i].idx),
            _transformedRows.rhs(violatedrows[i].idx));

         _idsReducedProbRows[violatedrows[i].idx] = true;
         newrowidx[nnewrowidx] = violatedrows[i].idx;
         nnewrowidx++;
      }

      SPxRowId* addedrowids = 0;
      spx_alloc(addedrowids, nnewrowidx);
      _solver.addRows(addedrowids, updaterows);

      for( int i = 0; i < nnewrowidx; i++ )
         _idsReducedProbRowIDs[newrowidx[i]] = addedrowids[i];

      // freeing allocated memory
      spx_free(addedrowids);
      spx_free(newrowidx);
   }



   /// update the reduced problem with additional columns and rows based upon the violated original bounds and rows
   void SoPlex::_updateIdsReducedProblemViol(bool allrows)
   {
#ifdef NO_TOL
      Real feastol = 0.0;
#else
#ifdef USE_FEASTOL
      Real feastol = realParam(SoPlex::FEASTOL);
#else
      Real feastol = realParam(SoPlex::EPSILON_ZERO);
#endif
#endif
      LPRowSet updaterows;

      int* newrowidx = 0;
      int nnewrowidx = 0;
      spx_alloc(newrowidx, _nPrimalRows);

      int rowNumber;
      int bestrow = -1;
      Real bestrownorm = infinity;
      Real percenttoadd = 1;

      int nrowstoadd = 1;  // adding the most violated row
      //if( allrows )
         nrowstoadd = _nIdsViolRows;   // adding all violated rows

      SSVector y(_solver.nCols());
      y.unSetup();

      for( int i = 0; i < nrowstoadd; i++ )
      {
         rowNumber = _idsViolatedRows[i];

         if( !allrows )
         {
            Real norm = 0;

            // the rhs of this calculation are the rows of the constraint matrix
            // so we are solving y B = A_{i,.}
            try
            {
               _solver.basis().solve(y, _solver.vector(rowNumber));
            }
            catch( SPxException E )
            {
               MSG_ERROR( spxout << "Caught exception <" << E.what() << "> while computing compatability.\n" );
            }


            if( y.isSetup() )
            {
               for( int j = 0; j < y.size(); j++ )
               {
                  if( isZero(_solver.fVec()[i], feastol) )
                     norm += spxAbs(y.value(j))*spxAbs(y.value(j));
               }
            }
            else
            {
               for( int j = 0; j < numColsReal(); j++ )
               {
                  if( isZero(_solver.fVec()[i], feastol) )
                     norm += spxAbs(y[j])*spxAbs(y[j]);
               }
            }


            norm = spxSqrt(norm);
            if( LT(norm, bestrownorm) )
            {
               bestrow = rowNumber;
               bestrownorm = norm;
            }

            if( isZero(norm, feastol) && LT(nnewrowidx/Real(numRowsReal()), percenttoadd) )
            {
               updaterows.add(_transformedRows.lhs(rowNumber), _transformedRows.rowVector(rowNumber),
                  _transformedRows.rhs(rowNumber));

               _idsReducedProbRows[rowNumber] = true;
               newrowidx[nnewrowidx] = rowNumber;
               nnewrowidx++;
            }

         }
         else
         {
            updaterows.add(_transformedRows.lhs(rowNumber), _transformedRows.rowVector(rowNumber),
                  _transformedRows.rhs(rowNumber));

            _idsReducedProbRows[rowNumber] = true;
            newrowidx[nnewrowidx] = rowNumber;
            nnewrowidx++;
         }
      }

      if( nnewrowidx == 0 )
      {
         for( int i = 0; i < nrowstoadd; i++ )
         {
            rowNumber = _idsViolatedRows[i];

            updaterows.add(_transformedRows.lhs(rowNumber), _transformedRows.rowVector(rowNumber),
                  _transformedRows.rhs(rowNumber));

            _idsReducedProbRows[rowNumber] = true;
            newrowidx[nnewrowidx] = rowNumber;
            nnewrowidx++;
         }
      }

      if( !allrows && bestrow >= 0 )
      {
         updaterows.add(_transformedRows.lhs(bestrow), _transformedRows.rowVector(bestrow),
            _transformedRows.rhs(bestrow));

         _idsReducedProbRows[bestrow] = true;
         newrowidx[nnewrowidx] = bestrow;
         nnewrowidx++;
      }


      SPxRowId* addedrowids = 0;
      spx_alloc(addedrowids, nnewrowidx);
      _solver.addRows(addedrowids, updaterows);

      for( int i = 0; i < nnewrowidx; i++ )
         _idsReducedProbRowIDs[newrowidx[i]] = addedrowids[i];

      // freeing allocated memory
      spx_free(addedrowids);
      spx_free(newrowidx);
   }



   /// builds the update rows with those violated in the complmentary problem
   // A row is violated in the constraint matrix Ax <= b, if b - A_{i}x < 0
   // To aid the computation, all violations are translated to <= constraints
   void SoPlex::_findViolatedRows(Real compObjValue, DataArray<RowViolation>& violatedrows, int& nviolatedrows)
   {
      Real feastol = realParam(SoPlex::FEASTOL);
      DVector compProbRedcost(_compSolver.nCols());   // the reduced costs of the complementary problem
      DVector compProbPrimal(_compSolver.nCols());    // the primal vector of the complementary problem
      DVector compProbActivity(_compSolver.nRows());  // the activity vector of the complementary problem
      Real compProbSlackVal = 0;

      if( boolParam(SoPlex::USECOMPDUAL) )
      {
         // Retrieving the slacks for each row.
         _compSolver.getRedCost(compProbRedcost);
      }
      else
      {
         // Retrieving the primal vector for the complementary problem
         _compSolver.getPrimal(compProbPrimal);
         _compSolver.computePrimalActivity(compProbPrimal, compProbActivity);

         compProbSlackVal = compProbPrimal[_compSolver.number(_compSlackColId)];
      }

      for( int i = 0; i < _nPrimalRows; i++ )
      {
         LPRowReal origlprow;
         DSVectorBase<Real> rowtoaddVec(_realLP->nCols());
         Real compProbViol = 0;
         Real compSlackCoeff = 0;
         int rownumber = _realLP->number(SPxRowId(_idsPrimalRowIDs[i]));
         int comprownum = _compSolver.number(SPxRowId(_idsPrimalRowIDs[i]));

         if( !_idsReducedProbRows[rownumber] )
         {
            // retreiving the violation of the complementary problem primal constraints
            if( boolParam(SoPlex::USECOMPDUAL) )
            {
               compSlackCoeff = getCompSlackVarCoeff(i);
               compProbViol = compProbRedcost[_compSolver.number(SPxColId(_idsDualColIDs[i]))]; // this is b - Ax
               // subtracting the slack variable value
               compProbViol += compObjValue*compSlackCoeff; // must add on the slack variable value.
               compProbViol *= compSlackCoeff;  // translating the violation to a <= constraint
            }
            else
            {
               Real viol = _compSolver.rhs(comprownum) - (compProbActivity[comprownum] + compProbSlackVal);
               if( viol < 0.0 )
                  compProbViol = viol;

               viol = (compProbActivity[comprownum] - compProbSlackVal) - _compSolver.lhs(comprownum);
               if( viol < 0.0 )
                  compProbViol = viol;

            }

            // NOTE: if the row was originally a ranged constraint, we are only interest in one of the inequalities.
            // If one inequality of the range violates the bounds, then we will add the row.

            // the translation of the complementary primal problem to the dual some rows resulted in two columns.
            if( boolParam(SoPlex::USECOMPDUAL) && i < _nPrimalRows - 1 &&
                  _realLP->number(SPxRowId(_idsPrimalRowIDs[i])) == _realLP->number(SPxRowId(_idsPrimalRowIDs[i + 1])) )
            {
               i++;
               compSlackCoeff = getCompSlackVarCoeff(i);
               Real tempViol = compProbRedcost[_compSolver.number(SPxColId(_idsDualColIDs[i]))]; // this is b - Ax
               tempViol += compObjValue*compSlackCoeff;
               tempViol *= compSlackCoeff;

               // if the other side of the range constraint has a larger violation, then this is used for the
               // computation.
               if( tempViol < compProbViol )
                  compProbViol = tempViol;
            }


            if( LT(compProbViol, 0, feastol) )
            {
               if( !_idsReducedProbRows[rownumber] )
               {
                  numIncludedRows++;
                  assert(numIncludedRows <= _realLP->nRows());
               }
               violatedrows[nviolatedrows].idx = rownumber;
               violatedrows[nviolatedrows].violation = spxAbs(compProbViol);
               nviolatedrows++;
            }
         }
      }
   }



   /// identifies the columns of the row-form basis that correspond to rows with zero dual multipliers.
   // This function assumes that the basis is in the row form.
   // @todo extend this to the case when the basis is in the column form.
   void SoPlex::_getZeroDualMultiplierIndices(Vector feasVector, int* nonposind, int* colsforremoval,
      int* nnonposind, bool& stop)
   {
      assert(_solver.rep() == SPxSolver::ROW);

#ifdef NO_TOL
      Real feastol = 0.0;
#else
#ifdef USE_FEASTOL
      Real feastol = realParam(SoPlex::FEASTOL);
#else
      Real feastol = realParam(SoPlex::EPSILON_ZERO);
#endif
#endif

      bool delCol;

      _idsReducedProbColIDs.reSize(numColsReal());
      *nnonposind = 0;

      // iterating over all columns in the basis matrix
      // this identifies the basis indices and the indicies that are positive.
      for( int i = 0; i < _solver.nCols(); ++i )
      {
         _idsReducedProbCols[i] = true;
         _idsReducedProbColIDs[i].inValidate();
         colsforremoval[i] = i;
         delCol = false;
         if( _solver.basis().baseId(i).isSPxRowId() ) // find the row id's for rows in the basis
         {
            //@todo need to check this regarding min and max problems
            if( isZero(feasVector[i], feastol) )
            {
               nonposind[*nnonposind] = i;
               (*nnonposind)++;

               // NOTE: commenting out the delCol flag at this point. The colsforremoval array should indicate the
               // columns that have a zero reduced cost. Hence, the delCol flag should only be set in the isSPxColId
               // branch of the if statement.
               //delCol = true;
            }
         }
         else if( _solver.basis().baseId(i).isSPxColId() )  // get the column id's for the columns in the basis
         {
            if( isZero(feasVector[i], feastol) )
            {
               nonposind[*nnonposind] = i;
               (*nnonposind)++;

               delCol = true;
            }
         }

         // setting an array to identify the columns to be removed from the LP to form the reduced problem
         if( delCol )
         {
            //colsforremoval[_solver.number(_solver.basis().baseId(i))] = -1;
            //_idsReducedProbCols[_solver.number(_solver.basis().baseId(i))] = false;
            colsforremoval[i] = -1;
            _idsReducedProbCols[i] = false;
         }
         else if( _solver.basis().baseId(i).isSPxColId() )
         {
#if 0
            colsforremoval[_solver.number(_solver.basis().baseId(i))] = _solver.number(_solver.basis().baseId(i));
            _idsReducedProbColIDs[_solver.number(_solver.basis().baseId(i))] = SPxColId(_solver.basis().baseId(i));
#endif
#if 1
            if( _solver.basis().baseId(i).isSPxColId() )
               _idsReducedProbColIDs[_solver.number(_solver.basis().baseId(i))] = SPxColId(_solver.basis().baseId(i));
            else
               _idsReducedProbCols[_solver.number(_solver.basis().baseId(i))] = false;
#endif
         }
      }

      stop = idsTerminate(realParam(SoPlex::TIMELIMIT)*TIMELIMIT_FRAC);
   }



   /// retrieves the compatible columns from the constraint matrix
   // This function also updates the constraint matrix of the reduced problem. It is efficient to perform this in the
   // following function because the required linear algebra has been performed.
   void SoPlex::_getCompatibleColumns(Vector feasVector, int* nonposind, int* compatind, int* rowsforremoval,
      int* colsforremoval, int nnonposind, int* ncompatind, bool formRedProb, bool& stop)
   {
#ifdef NO_TOL
      Real feastol = 0.0;
#else
#ifdef USE_FEASTOL
      Real feastol = realParam(SoPlex::FEASTOL);
#else
      Real feastol = realParam(SoPlex::EPSILON_ZERO);
#endif
#endif

      bool compatible;
      SSVector y(_solver.nCols());
      y.unSetup();

      *ncompatind  = 0;

#ifndef NDEBUG
      int bind = 0;
      bool* activerows = 0;
      spx_alloc(activerows, numRowsReal());

      for( int i = 0; i < numRowsReal(); ++i )
         activerows[i] = false;

      for( int i = 0; i < numColsReal(); ++i )
      {
         if( _solver.basis().baseId(i).isSPxRowId() ) // find the row id's for rows in the basis
         {
            if( !isZero(feasVector[i], feastol) )
            {
               bind = _realLP->number(SPxRowId(_solver.basis().baseId(i))); // getting the corresponding row
                                                                                   // for the original LP.
               assert(bind >= 0 && bind < numRowsReal());
               activerows[bind] = true;
            }
         }
      }
#endif

      if( formRedProb )
      {
         _idsReducedProbRowIDs.reSize(_solver.nRows());
         _idsReducedProbColRowIDs.reSize(_solver.nRows());
      }

      for( int i = 0; i < numRowsReal(); ++i )
      {
         rowsforremoval[i] = i;
         if( formRedProb )
            _idsReducedProbRows[i] = true;

         numIncludedRows++;

         // the rhs of this calculation are the rows of the constraint matrix
         // so we are solving y B = A_{i,.}
         try
         {
            _solver.basis().solve(y, _solver.vector(i));
         }
         catch( SPxException E )
         {
            MSG_ERROR( spxout << "Caught exception <" << E.what() << "> while computing compatability.\n" );
         }

         compatible = true;
         // a compatible row is given by zeros in all columns related to the nonpositive indices
         for( int j = 0; j < nnonposind; ++j )
         {
            // @todo really need to check this part of the code. Run through this with Ambros or Matthias.
            // 09.02.2014 getting a tolerance issue with this check. Don't know how to fix it.
            if( !isZero(y[nonposind[j]], feastol) )
            {
               compatible = false;
               break;
            }
         }

         // checking that the active rows are compatible
         assert(!activerows[i] || compatible);

         // changing the matrix coefficients
         DSVector newRowVector;
         LPRowReal rowtoupdate;

         if( y.isSetup() )
         {
            for( int j = 0; j < y.size(); j++ )
               newRowVector.add(y.index(j), y.value(j));
         }
         else
         {
            for( int j = 0; j < numColsReal(); j++ )
            {
               if( !isZero(y[j], feastol) )
                  newRowVector.add(j, y[j]);
            }
         }

         // transforming the original problem rows
         _solver.getRow(i, rowtoupdate);

#ifndef NO_TRANSFORM
         rowtoupdate.setRowVector(newRowVector);
#endif

         if( formRedProb )
            _transformedRows.add(rowtoupdate);


         // Making all equality constraints compatible, i.e. they are included in the reduced problem
         if( EQ(rowtoupdate.lhs(), rowtoupdate.rhs()) )
            compatible = true;

         if( compatible )
         {
            compatind[*ncompatind] = i;
            (*ncompatind)++;

            if( formRedProb )
            {
               _idsReducedProbRowIDs[i] = _solver.rowId(i);

               // updating the compatible row
               _idsLP->changeRow(i, rowtoupdate);
            }

#if 0
            // checking the columns contained in the removed row to ensure that all are removed.
            LPRowReal lprow;
            _realLP->getRow(i, lprow);
            for( int j = 0; j < lprow.rowVector().size(); j++ )
            {
               // if a column not marked for removal exists in the row, the row remains.
               if( colsforremoval[lprow.rowVector().index(j)] < 0 )
               {
                  colsforremoval[lprow.rowVector().index(j)] = lprow.rowVector().index(j);
                  _idsReducedProbCols[lprow.rowVector().index(j)] = true;
                  _idsReducedProbColIDs[lprow.rowVector().index(j)] = SPxColId(_realLP->cId(lprow.rowVector().index(j)));
               }
            }
#endif
         }
         else
         {
            // setting an array to identify the rows to be removed from the LP to form the reduced problem
            rowsforremoval[i] = -1;
            numIncludedRows--;

            if( formRedProb )
               _idsReducedProbRows[i] = false;
#if 0
            // checking the columns contained in the removed row to ensure that all are removed.
            LPRowReal lprow;
            _realLP->getRow(i, lprow);
            for( int j = 0; j < lprow.rowVector().size(); j++ )
            {
               // if a column not marked for removal exists in the row, the row remains.
               if( colsforremoval[lprow.rowVector().index(j)] >= 0 )
               {
                  rowsforremoval[i] = i;
                  _idsReducedProbRows[i] = true;

                  compatind[*ncompatind] = i;
                  (*ncompatind)++;

                  _idsReducedProbRowIDs[i] = _solver.rowId(i);
                  break;
               }
            }
#endif
         }

         // determine whether the reduced problem setup should be terminated
         stop = idsTerminate(realParam(SoPlex::TIMELIMIT)*TIMELIMIT_FRAC);

         if( stop )
            break;
      }
      assert(numIncludedRows <= _solver.nRows());

#ifndef NDEBUG
      spx_free(activerows);
#endif
   }



   /// computes the reduced problem objective coefficients
   void SoPlex::_computeReducedProbObjCoeff(bool& stop)
   {
#ifdef NO_TOL
      Real feastol = 0.0;
#else
#ifdef USE_FEASTOL
      Real feastol = realParam(SoPlex::FEASTOL);
#else
      Real feastol = realParam(SoPlex::EPSILON_ZERO);
#endif
#endif

      SSVector y(numColsReal());
      y.unSetup();

      // the rhs of this calculation is the original objective coefficient vector
      // so we are solving y B = c
      try
      {
         _solver.basis().solve(y, _solver.maxObj());
      }
      catch( SPxException E )
      {
         MSG_ERROR( spxout << "Caught exception <" << E.what() << "> while computing compatability.\n" );
      }

      _transformedObj.reDim(numColsReal());
      if( y.isSetup() )
      {
         int ycount = 0;
         for( int i = 0; i < numColsReal(); i++ )
         {
            if( ycount < y.size() && i == y.index(ycount) )
            {
               _transformedObj[i] = y.value(ycount);
               ycount++;
            }
            else
               _transformedObj[i] = 0.0;
         }
      }
      else
      {
         for( int i = 0; i < numColsReal(); i++ )
         {
            if( isZero(y[i], feastol) )
               _transformedObj[i] = 0.0;
            else
               _transformedObj[i] = y[i];
         }
      }

      // setting the updated objective vector
#ifndef NO_TRANSFORM
      _idsLP->changeObj(_transformedObj);
#endif

      // determine whether the reduced problem setup should be terminated
      stop = idsTerminate(realParam(SoPlex::TIMELIMIT)*TIMELIMIT_FRAC);
   }



   /// computes the compatible bound constraints and adds them to the reduced problem
   // NOTE: No columns are getting removed from the reduced problem. Only the bound constraints are being removed.
   // So in the reduced problem, all variables are free unless the bound constraints are selected as compatible.
   void SoPlex::_getCompatibleBoundCons(LPRowSet& boundcons, int* compatboundcons, int* nonposind,
         int* ncompatboundcons, int nnonposind, bool& stop)
   {
#ifdef NO_TOL
      Real feastol = 0.0;
#else
#ifdef USE_FEASTOL
      Real feastol = realParam(SoPlex::FEASTOL);
#else
      Real feastol = realParam(SoPlex::EPSILON_ZERO);
#endif
#endif

      bool compatible;
      SSVector y(numColsReal());
      y.unSetup();

      _idsReducedProbColRowIDs.reSize(numColsReal());


      // identifying the compatible bound constraints
      for( int i = 0; i < numColsReal(); i++ )
      {
         _idsReducedProbColRowIDs[i].inValidate();

         // setting all variables to free variables.
         // Bound constraints will only be added to the variables with compatible bound constraints.
         _idsLP->changeUpper(i, infinity);
         _idsLP->changeLower(i, -infinity);

         // the rhs of this calculation are the unit vectors of the bound constraints
         // so we are solving y B = I_{i,.}

         try
         {
            _solver.basis().solve(y, _solver.unitVector(i));
         }
         catch( SPxException E )
         {
            MSG_ERROR( spxout << "Caught exception <" << E.what() << "> while computing compatability.\n" );
         }

         compatible = true;
         // a compatible row is given by zeros in all columns related to the nonpositive indices
         for( int j = 0; j < nnonposind; ++j )
         {
            // @todo really need to check this part of the code. Run through this with Ambros or Matthias.
            if( !isZero(y[nonposind[j]]) )
            {
               compatible = false;
               break;
            }
         }

         // changing the matrix coefficients
         DSVector newRowVector;
         LPRowReal rowtoupdate;

#ifndef NO_TRANSFORM
         if( y.isSetup() )
         {
            for( int j = 0; j < y.size(); j++ )
               newRowVector.add(y.index(j), y.value(j));
         }
         else
         {
            for( int j = 0; j < numColsReal(); j++ )
            {
               if( !isZero(y[j], feastol) )
               {
                  newRowVector.add(j, y[j]);
               }
            }
         }
#else
         newRowVector.add(i, 1.0);
#endif

         // will probably need to map this set of rows.
         _transformedRows.add(_solver.lower(i), newRowVector, _solver.upper(i));

         // variable bounds are compatible in the current implementation.
         compatible = true;

         // if the bound constraint is compatible, it remains in the reduced problem.
         // Otherwise the bound is removed and the variables are free.
         if( compatible )
         {
            Real lhs = -infinity;
            Real rhs = infinity;
            if( GT(_solver.lower(i), -infinity) )
               lhs = _solver.lower(i);

            if( LT(_solver.upper(i), infinity) )
               rhs = _solver.upper(i);

            if( GT(lhs, -infinity) || LT(rhs, infinity) )
            {
               compatboundcons[(*ncompatboundcons)] = i;
               (*ncompatboundcons)++;

               boundcons.add(lhs, newRowVector, rhs);
            }


         }

         // determine whether the reduced problem setup should be terminated
         stop = idsTerminate(realParam(SoPlex::TIMELIMIT)*TIMELIMIT_FRAC);

         if( stop )
            break;
      }

   }



   /// computes the rows to remove from the complementary problem
   void SoPlex::_getRowsForRemovalComplementaryProblem(int* nonposind, int* bind, int* rowsforremoval,
         int* nrowsforremoval, int nnonposind)
   {
      *nrowsforremoval = 0;

      for( int i = 0; i < nnonposind; i++ )
      {
         if( bind[nonposind[i]] < 0 )
         {
            rowsforremoval[*nrowsforremoval] = -1 - bind[nonposind[i]];
            (*nrowsforremoval)++;
         }
      }
   }



   /// removing rows from the complementary problem.
   // the rows that are removed from idsCompLP are the rows from the reduced problem that have a non-positive dual
   // multiplier in the optimal solution.
   void SoPlex::_deleteAndUpdateRowsComplementaryProblem(SPxRowId rangedRowIds[], int& naddedrows)
   {
      DSVector slackColCoeff;

      // setting the objective coefficients of the original variables to zero
      DVector newObjCoeff(numColsReal());
      for( int i = 0; i < numColsReal(); i++ )
      {
         _compSolver.changeBounds(_realLP->cId(i), Real(-infinity), Real(infinity));
         newObjCoeff[i] = 0;
      }

      _compSolver.changeObj(newObjCoeff);

      // adding the slack column to the complementary problem
      SPxColId* addedcolid = 0;
      if( boolParam(SoPlex::USECOMPDUAL) )
      {
         spx_alloc(addedcolid, 1);
         LPColSetReal compSlackCol;
         compSlackCol.add(1.0, -infinity, slackColCoeff, infinity);
         _compSolver.addCols(addedcolid, compSlackCol);
         _compSlackColId = addedcolid[0];
      }
      else
      {
         LPRowSet addrangedrows;    // the row set of ranged and equality rows that must be added to the complementary problem.
         naddedrows = 0;
         // finding all of the ranged and equality rows and creating two <= constraints.
         for( int i = 0; i < numRowsReal(); i++ )
         {
            if( _realLP->rowType(i) == LPRowBase<Real>::RANGE || _realLP->rowType(i) == LPRowBase<Real>::EQUAL )
            {
               assert(GT(_compSolver.lhs(i), -infinity) && LT(_compSolver.rhs(i), infinity));
               assert(_compSolver.rowType(i) == LPRowBase<Real>::RANGE || _compSolver.rowType(i) == LPRowBase<Real>::EQUAL);

               _compSolver.changeLhs(i, -infinity);
               addrangedrows.add(_realLP->lhs(i), _realLP->rowVector(i), infinity);
               naddedrows++;
            }
         }

         // adding the rows for the ranged rows to make <= conatraints
         SPxRowId* addedrowid = 0;
         spx_alloc(addedrowid, naddedrows);
         _compSolver.addRows(addedrowid, addrangedrows);

         // collecting the row ids
         for( int i = 0; i < naddedrows; i++ )
            rangedRowIds[i] = addedrowid[i];

         spx_free(addedrowid);


         // adding the slack column
         spx_alloc(addedcolid, 1);
         LPColSetReal compSlackCol;
         compSlackCol.add(-1.0, 0.0, slackColCoeff, infinity);

         _compSolver.addCols(addedcolid, compSlackCol);
         _compSlackColId = addedcolid[0];
      }

      // freeing allocated memory
      spx_free(addedcolid);
   }



   /// update the dual complementary problem with additional columns and rows
   // Given the solution to the updated reduced problem, the complementary problem will be updated with modifications to
   // the constraints and the removal of variables
   void SoPlex::_updateIdsComplementaryDualProblem(bool origObj)
   {
      Real feastol = realParam(SoPlex::FEASTOL);

      int prevNumCols = _compSolver.nCols(); // number of columns in the previous formulation of the complementary prob
      int prevNumRows = _compSolver.nRows();
      int prevPrimalRowIds = _nPrimalRows;
      int prevDualColIds = _nDualCols;

      LPColSetBase<Real> addElimCols(_nElimPrimalRows);  // columns previously eliminated from the
                                                         // complementary problem that must be added
      int numElimColsAdded = 0;
      int currnumrows = prevNumRows;
      // looping over all rows from the original LP that were eliminated during the formation of the complementary
      // problem. The eliminated rows will be added if they are basic in the reduced problem.

      for( int i = 0; i < _nElimPrimalRows; i++ )
      {
         int rowNumber = _realLP->number(_idsElimPrimalRowIDs[i]);

         int solverRowNum = _solver.number(_idsReducedProbRowIDs[rowNumber]);
         assert(solverRowNum >= 0 && solverRowNum < _solver.nRows());

         // checking for the basic rows in the reduced problem
         if( _solver.basis().desc().rowStatus(solverRowNum) == SPxBasis::Desc::P_ON_UPPER ||
           _solver.basis().desc().rowStatus(solverRowNum) == SPxBasis::Desc::P_ON_LOWER ||
           _solver.basis().desc().rowStatus(solverRowNum) == SPxBasis::Desc::P_FIXED ||
           _solver.basis().desc().rowStatus(solverRowNum) == SPxBasis::Desc::D_FREE ||
           (_solver.basis().desc().rowStatus(solverRowNum) == SPxBasis::Desc::D_ON_LOWER &&
           LE(_solver.rhs(solverRowNum) - _solver.pVec()[solverRowNum], 0.0, feastol)) ||
           (_solver.basis().desc().rowStatus(solverRowNum) == SPxBasis::Desc::D_ON_UPPER &&
           LE(_solver.pVec()[solverRowNum] - _solver.lhs(solverRowNum), 0.0, feastol)) )
         {
            LPRowReal origlprow;
            DSVectorBase<Real> coltoaddVec(_realLP->nCols());

            LPRowSet additionalrows;
            int nnewrows = 0;

            _realLP->getRow(rowNumber, origlprow);
            for( int j = 0; j < origlprow.rowVector().size(); j++ )
            {
               // the column of the new row may not exist in the current complementary problem.
               // if the column does not exist, then it is necessary to create the column.
               int colNumber = origlprow.rowVector().index(j);
               if( _idsCompProbColIDsIdx[colNumber] == -1 )
               {
                  assert(!_idsReducedProbColIDs[colNumber].isValid());
                  _idsPrimalColIDs[_nPrimalCols] = _realLP->cId(colNumber);
                  _idsCompProbColIDsIdx[colNumber] = _nPrimalCols;
                  _fixedOrigVars[colNumber] = -2;
                  _nPrimalCols++;

                  // all columns for the complementary problem are converted to unrestricted.
                  additionalrows.create(1, _realLP->maxObj(colNumber), _realLP->maxObj(colNumber));
                  nnewrows++;

                  coltoaddVec.add(currnumrows, origlprow.rowVector().value(j));
                  currnumrows++;
               }
               else
                  coltoaddVec.add(_compSolver.number(_idsDualRowIDs[_idsCompProbColIDsIdx[colNumber]]),
                        origlprow.rowVector().value(j));
            }

            SPxRowId* addedrowids = 0;
            spx_alloc(addedrowids, nnewrows);
            _compSolver.addRows(addedrowids, additionalrows);

            for( int j = 0; j < nnewrows; j++ )
            {
               _idsDualRowIDs[_nDualRows] = addedrowids[j];
               _nDualRows++;
            }

            spx_free(addedrowids);



            if( _solver.basis().desc().rowStatus(solverRowNum) == SPxBasis::Desc::P_ON_UPPER ||
                  _solver.basis().desc().rowStatus(solverRowNum) == SPxBasis::Desc::P_FIXED ||
                  _solver.basis().desc().rowStatus(solverRowNum) == SPxBasis::Desc::D_FREE ||
                 (_solver.basis().desc().rowStatus(solverRowNum) == SPxBasis::Desc::D_ON_LOWER &&
                 LE(_solver.rhs(solverRowNum) - _solver.pVec()[solverRowNum], 0.0, feastol)) )
            {
               assert(LT(_realLP->rhs(_idsElimPrimalRowIDs[i]), infinity));
               addElimCols.add(_realLP->rhs(_idsElimPrimalRowIDs[i]), -infinity, coltoaddVec, infinity);

               if( _nPrimalRows >= _idsPrimalRowIDs.size() )
               {
                  _idsPrimalRowIDs.reSize(_nPrimalRows*2);
                  _idsDualColIDs.reSize(_nPrimalRows*2);
               }

               _idsPrimalRowIDs[_nPrimalRows] = _idsElimPrimalRowIDs[i];
               _nPrimalRows++;

               _idsElimPrimalRowIDs.remove(i);
               _nElimPrimalRows--;
               i--;

               numElimColsAdded++;
            }
            else if( _solver.basis().desc().rowStatus(solverRowNum) == SPxBasis::Desc::P_ON_LOWER ||
              (_solver.basis().desc().rowStatus(solverRowNum) == SPxBasis::Desc::D_ON_UPPER &&
              LE(_solver.pVec()[solverRowNum] - _solver.lhs(solverRowNum), 0.0, feastol)) )
            {
               // this assert should stay, but there is an issue with the status and the dual vector
               //assert(LT(dualVector[_solver.number(_idsReducedProbRowIDs[rowNumber])], 0.0));
               assert(GT(_realLP->lhs(_idsElimPrimalRowIDs[i]), -infinity));
               addElimCols.add(_realLP->lhs(_idsElimPrimalRowIDs[i]), -infinity, coltoaddVec, infinity);

               _idsPrimalRowIDs[_nPrimalRows] = _idsElimPrimalRowIDs[i];
               _nPrimalRows++;

               _idsElimPrimalRowIDs.remove(i);
               _nElimPrimalRows--;
               i--;

               numElimColsAdded++;
            }
         }
      }

      MSG_INFO2(spxout, spxout << "Number of eliminated columns added to complementary problem: "
         << numElimColsAdded << std::endl );

      // updating the _idsDualColIDs with the additional columns from the eliminated rows.
      _compSolver.addCols(addElimCols);
      for( int i = prevNumCols; i < _compSolver.nCols(); i++ )
      {
         _idsDualColIDs[prevDualColIds + i - prevNumCols] = _compSolver.colId(i);
         _nDualCols++;
      }

      assert(_nDualCols == _nPrimalRows);

      // looping over all rows from the original problem that were originally contained in the complementary problem.
      // The basic rows will be set as free variables, the non-basic rows will be eliminated from the complementary
      // problem.
      DSVector slackRowCoeff(_compSolver.nCols());

      int* colsforremoval = 0;
      int ncolsforremoval = 0;
      spx_alloc(colsforremoval, prevPrimalRowIds);
      for( int i = 0; i < prevPrimalRowIds; i++ )
      {
         int rowNumber = _realLP->number(_idsPrimalRowIDs[i]);
         // this loop runs over all rows previously in the complementary problem. If rows are added to the reduced
         // problem, they will be transfered from the incompatible set to the compatible set in the following if
         // statement.
         if( _idsReducedProbRows[rowNumber] )
         {
            // rows added to the reduced problem may have been equality constriants. The equality constraints from the
            // original problem are converted into <= and >= constraints. Upon adding these constraints to the reduced
            // problem, only a single dual column is needed in the complementary problem. Hence, one of the dual columns
            // is removed.
            //
            // 22.06.2015 Testing keeping all constraints in the complementary problem.
            // This requires a dual column to be fixed to zero for the range and equality rows.
            bool incrementI = false;
            if( i + 1 < prevPrimalRowIds
                  && _realLP->number(_idsPrimalRowIDs[i]) == _realLP->number(_idsPrimalRowIDs[i+1]) )
            {
               assert(_idsPrimalRowIDs[i].idx == _idsPrimalRowIDs[i+1].idx);

#if 0 // 22.06.2015
               if( _realLP->rowType(_idsPrimalRowIDs[i]) == LPRowBase<Real>::RANGE )
               {
                  _compSolver.changeObj(_idsDualColIDs[i + 1], 0.0);
                  _compSolver.changeBounds(_idsDualColIDs[i + 1], 0.0, 0.0);
               }
               incrementI = true;
#else
               colsforremoval[ncolsforremoval] = _compSolver.number(SPxColId(_idsDualColIDs[i + 1]));
               ncolsforremoval++;

               _idsPrimalRowIDs.remove(i + 1);
               _nPrimalRows--;
               _idsDualColIDs.remove(i + 1);
               _nDualCols--;

               prevPrimalRowIds--;
#endif
            }
            //assert(i + 1 == prevPrimalRowIds || _idsPrimalRowIDs[i].idx != _idsPrimalRowIDs[i+1].idx);

            int solverRowNum = _solver.number(_idsReducedProbRowIDs[rowNumber]);
            assert(solverRowNum >= 0 && solverRowNum < _solver.nRows());
            if( _solver.basis().desc().rowStatus(solverRowNum) == SPxBasis::Desc::P_ON_UPPER ||
                  _solver.basis().desc().rowStatus(solverRowNum) == SPxBasis::Desc::P_FIXED ||
                  _solver.basis().desc().rowStatus(solverRowNum) == SPxBasis::Desc::D_FREE ||
                 (_solver.basis().desc().rowStatus(solverRowNum) == SPxBasis::Desc::D_ON_LOWER &&
                 LE(_solver.rhs(solverRowNum) - _solver.pVec()[solverRowNum], 0.0, feastol)) )
            {
               //assert(GT(dualVector[solverRowNum], 0.0));
               _compSolver.changeObj(_idsDualColIDs[i], _realLP->rhs(SPxRowId(_idsPrimalRowIDs[i])));
               _compSolver.changeBounds(_idsDualColIDs[i], -infinity, infinity);
            }
            else if( _solver.basis().desc().rowStatus(solverRowNum) == SPxBasis::Desc::P_ON_LOWER ||
              (_solver.basis().desc().rowStatus(solverRowNum) == SPxBasis::Desc::D_ON_UPPER &&
              LE(_solver.pVec()[solverRowNum] - _solver.lhs(solverRowNum), 0.0, feastol)) )
            {
               //assert(LT(dualVector[solverRowNum], 0.0));
               _compSolver.changeObj(_idsDualColIDs[i], _realLP->lhs(SPxRowId(_idsPrimalRowIDs[i])));
               _compSolver.changeBounds(_idsDualColIDs[i], -infinity, infinity);
            }
            else //if ( _solver.basis().desc().rowStatus(solverRowNum) != SPxBasis::Desc::D_FREE )
            {
               //assert(isZero(dualVector[solverRowNum], 0.0));

               // 22.06.2015 Testing keeping all rows in the complementary problem
#if 0
               switch( _realLP->rowType(_idsPrimalRowIDs[i]) )
               {
                  case LPRowBase<Real>::RANGE:
                     assert(_realLP->number(SPxColId(_idsPrimalRowIDs[i])) ==
                           _realLP->number(SPxColId(_idsPrimalRowIDs[i+1])));
                     assert(_compSolver.obj(_compSolver.number(SPxColId(_idsDualColIDs[i]))) !=
                           _compSolver.obj(_compSolver.number(SPxColId(_idsDualColIDs[i + 1]))));

                     _compSolver.changeObj(_idsDualColIDs[i], _realLP->rhs(_idsPrimalRowIDs[i]));
                     _compSolver.changeBounds(_idsDualColIDs[i], 0.0, infinity);
                     _compSolver.changeObj(_idsDualColIDs[i + 1], _realLP->lhs(_idsPrimalRowIDs[i]));
                     _compSolver.changeBounds(_idsDualColIDs[i + 1], -infinity, 0.0);

                     i++;
                     break;
                  case LPRowBase<Real>::GREATER_EQUAL:
                     _compSolver.changeObj(_idsDualColIDs[i], _realLP->lhs(_idsPrimalRowIDs[i]));
                     _compSolver.changeBounds(_idsDualColIDs[i], -infinity, 0.0);
                     break;
                  case LPRowBase<Real>::LESS_EQUAL:
                     _compSolver.changeObj(_idsDualColIDs[i], _realLP->rhs(_idsPrimalRowIDs[i]));
                     _compSolver.changeBounds(_idsDualColIDs[i], 0.0, infinity);
                     break;
                  default:
                     throw SPxInternalCodeException("XIDSSL01 This should never happen.");
               }

#else // 22.06.2015 testing keeping all rows in the complementary problem
               colsforremoval[ncolsforremoval] = _compSolver.number(SPxColId(_idsDualColIDs[i]));
               ncolsforremoval++;

               if( _nElimPrimalRows >= _idsElimPrimalRowIDs.size() )
                  _idsElimPrimalRowIDs.reSize(_realLP->nRows());

               _idsElimPrimalRowIDs[_nElimPrimalRows] = _idsPrimalRowIDs[i];
               _nElimPrimalRows++;
               _idsPrimalRowIDs.remove(i);
               _nPrimalRows--;
               _idsDualColIDs.remove(i);
               _nDualCols--;

               i--;
               prevPrimalRowIds--;
#endif
            }

            if( incrementI )
               i++;
         }
         else
         {
            switch( _realLP->rowType(_idsPrimalRowIDs[i]) )
            {
               case LPRowBase<Real>::RANGE:
                  assert(_realLP->number(SPxColId(_idsPrimalRowIDs[i])) ==
                        _realLP->number(SPxColId(_idsPrimalRowIDs[i+1])));
                  assert(_compSolver.obj(_compSolver.number(SPxColId(_idsDualColIDs[i]))) !=
                        _compSolver.obj(_compSolver.number(SPxColId(_idsDualColIDs[i + 1]))));
                  if( _compSolver.obj(_compSolver.number(SPxColId(_idsDualColIDs[i]))) <
                        _compSolver.obj(_compSolver.number(SPxColId(_idsDualColIDs[i + 1]))))
                  {
                     slackRowCoeff.add(_compSolver.number(SPxColId(_idsDualColIDs[i])), -SLACKCOEFF);
                     slackRowCoeff.add(_compSolver.number(SPxColId(_idsDualColIDs[i + 1])), SLACKCOEFF);
                  }
                  else
                  {
                     slackRowCoeff.add(_compSolver.number(SPxColId(_idsDualColIDs[i])), SLACKCOEFF);
                     slackRowCoeff.add(_compSolver.number(SPxColId(_idsDualColIDs[i + 1])), -SLACKCOEFF);
                  }


                  i++;
                  break;
               case LPRowBase<Real>::EQUAL:
                  assert(_realLP->number(SPxColId(_idsPrimalRowIDs[i])) ==
                        _realLP->number(SPxColId(_idsPrimalRowIDs[i+1])));

                  slackRowCoeff.add(_compSolver.number(SPxColId(_idsDualColIDs[i])), SLACKCOEFF);
                  slackRowCoeff.add(_compSolver.number(SPxColId(_idsDualColIDs[i + 1])), SLACKCOEFF);

                  i++;
                  break;
               case LPRowBase<Real>::GREATER_EQUAL:
                  slackRowCoeff.add(_compSolver.number(SPxColId(_idsDualColIDs[i])), -SLACKCOEFF);
                  break;
               case LPRowBase<Real>::LESS_EQUAL:
                  slackRowCoeff.add(_compSolver.number(SPxColId(_idsDualColIDs[i])), SLACKCOEFF);
                  break;
               default:
                  throw SPxInternalCodeException("XIDSSL01 This should never happen.");
            }

            if( origObj )
            {
               int numRemove = 1;
               int removeCount = 0;
               if( _realLP->number(SPxColId(_idsPrimalRowIDs[i])) ==
                        _realLP->number(SPxColId(_idsPrimalRowIDs[i+1])) )
                  numRemove++;

               do
               {
                  colsforremoval[ncolsforremoval] = _compSolver.number(SPxColId(_idsDualColIDs[i]));
                  ncolsforremoval++;

                  if( _nElimPrimalRows >= _idsElimPrimalRowIDs.size() )
                     _idsElimPrimalRowIDs.reSize(_realLP->nRows());

                  _idsElimPrimalRowIDs[_nElimPrimalRows] = _idsPrimalRowIDs[i];
                  _nElimPrimalRows++;
                  _idsPrimalRowIDs.remove(i);
                  _nPrimalRows--;
                  _idsDualColIDs.remove(i);
                  _nDualCols--;

                  i--;
                  prevPrimalRowIds--;

                  removeCount++;
               } while( removeCount < numRemove );
            }
         }
      }

      // updating the slack column in the complementary problem
      Real lhs = 1.0;
      Real rhs = 1.0;
      LPRowBase<Real> compSlackRow(lhs, slackRowCoeff, rhs);
      _compSolver.changeRow(_compSlackDualRowId, compSlackRow);


      // if the original objective is used, then all dual columns related to primal rows not in the reduced problem are
      // removed from the complementary problem.
      // As a result, the slack row becomes an empty row.
      int* perm = 0;
      spx_alloc(perm, _compSolver.nCols() + numElimColsAdded);
      _compSolver.removeCols(colsforremoval, ncolsforremoval, perm);


      // updating the dual columns to represent the fixed primal variables.
      int* currFixedVars = 0;
      spx_alloc(currFixedVars, _realLP->nCols());
      _identifyComplementaryDualFixedPrimalVars(currFixedVars);
      _removeComplementaryDualFixedPrimalVars(currFixedVars);
      _updateComplementaryDualFixedPrimalVars(currFixedVars);


      // freeing allocated memory
      spx_free(currFixedVars);
      spx_free(perm);
      spx_free(colsforremoval);
   }



   /// update the primal complementary problem with additional columns and rows
   // Given the solution to the updated reduced problem, the complementary problem will be updated with modifications to
   // the constraints and the removal of variables
   void SoPlex::_updateIdsComplementaryPrimalProblem(bool origObj)
   {
      Real feastol = realParam(SoPlex::FEASTOL);

      int prevNumRows = _compSolver.nRows();
      int prevPrimalRowIds = _nPrimalRows;

      assert(_nPrimalRows == _nCompPrimalRows);

      LPRowSetBase<Real> addElimRows(_nElimPrimalRows);  // rows previously eliminated from the
                                                         // complementary problem that must be added
      int numElimRowsAdded = 0;
      // looping over all rows from the original LP that were eliminated during the formation of the complementary
      // problem. The eliminated rows will be added if they are basic in the reduced problem.

      for( int i = 0; i < _nElimPrimalRows; i++ )
      {
         int rowNumber = _realLP->number(_idsElimPrimalRowIDs[i]);

         int solverRowNum = _solver.number(_idsReducedProbRowIDs[rowNumber]);
         assert(solverRowNum >= 0 && solverRowNum < _solver.nRows());

         // checking the rows that are basic in the reduced problem that should be added to the complementary problem
         if( _solver.basis().desc().rowStatus(solverRowNum) == SPxBasis::Desc::P_ON_UPPER ||
           _solver.basis().desc().rowStatus(solverRowNum) == SPxBasis::Desc::P_ON_LOWER ||
           _solver.basis().desc().rowStatus(solverRowNum) == SPxBasis::Desc::P_FIXED ||
           _solver.basis().desc().rowStatus(solverRowNum) == SPxBasis::Desc::D_FREE ||
           (_solver.basis().desc().rowStatus(solverRowNum) == SPxBasis::Desc::D_ON_LOWER &&
           EQ(_solver.rhs(solverRowNum) - _solver.pVec()[solverRowNum], 0.0, feastol)) ||
           (_solver.basis().desc().rowStatus(solverRowNum) == SPxBasis::Desc::D_ON_UPPER &&
           EQ(_solver.pVec()[solverRowNum] - _solver.lhs(solverRowNum), 0.0, feastol)) )
         {
            LPRowReal origlprow;
            _realLP->getRow(rowNumber, origlprow);

            // NOTE: 11.02.2016 I am assuming that all columns from the original problem are contained in the
            // complementary problem. Will need to check this. Since nothing is happenning in the
            // _deleteAndUpdateRowsComplementaryProblem function, I am feeling confident that all columns remain.


            if( _solver.basis().desc().rowStatus(solverRowNum) == SPxBasis::Desc::P_ON_UPPER ||
                  _solver.basis().desc().rowStatus(solverRowNum) == SPxBasis::Desc::P_FIXED ||
                  _solver.basis().desc().rowStatus(solverRowNum) == SPxBasis::Desc::D_FREE ||
                 (_solver.basis().desc().rowStatus(solverRowNum) == SPxBasis::Desc::D_ON_LOWER &&
                 EQ(_solver.rhs(solverRowNum) - _solver.pVec()[solverRowNum], 0.0, feastol)) )
            {
               assert(LT(_realLP->rhs(_idsElimPrimalRowIDs[i]), infinity));

               if( _nPrimalRows >= _idsPrimalRowIDs.size() )
               {
                  _idsPrimalRowIDs.reSize(_nPrimalRows*2);
                  _idsCompPrimalRowIDs.reSize(_nPrimalRows*2);
               }

               addElimRows.add(_realLP->rhs(_idsElimPrimalRowIDs[i]), origlprow.rowVector(),
                  _realLP->rhs(_idsElimPrimalRowIDs[i]));

               _idsPrimalRowIDs[_nPrimalRows] = _idsElimPrimalRowIDs[i];
               _nPrimalRows++;

               _idsElimPrimalRowIDs.remove(i);
               _nElimPrimalRows--;
               i--;

               numElimRowsAdded++;
            }
            else if( _solver.basis().desc().rowStatus(solverRowNum) == SPxBasis::Desc::P_ON_LOWER ||
              (_solver.basis().desc().rowStatus(solverRowNum) == SPxBasis::Desc::D_ON_UPPER &&
              EQ(_solver.pVec()[solverRowNum] - _solver.lhs(solverRowNum), 0.0, feastol)) )
            {
               assert(GT(_realLP->lhs(_idsElimPrimalRowIDs[i]), -infinity));

               if( _nPrimalRows >= _idsPrimalRowIDs.size() )
               {
                  _idsPrimalRowIDs.reSize(_nPrimalRows*2);
                  _idsCompPrimalRowIDs.reSize(_nPrimalRows*2);
               }

               addElimRows.add(_realLP->lhs(_idsElimPrimalRowIDs[i]), origlprow.rowVector(),
                  _realLP->lhs(_idsElimPrimalRowIDs[i]));

               _idsPrimalRowIDs[_nPrimalRows] = _idsElimPrimalRowIDs[i];
               _nPrimalRows++;

               _idsElimPrimalRowIDs.remove(i);
               _nElimPrimalRows--;
               i--;

               numElimRowsAdded++;
            }
         }
      }

      MSG_INFO2(spxout, spxout << "Number of eliminated rows added to the complementary problem: "
         << numElimRowsAdded << std::endl );

      // adding the eliminated rows to the complementary problem.
      _compSolver.addRows(addElimRows);
      for( int i = prevNumRows; i < _compSolver.nRows(); i++ )
      {
         _idsCompPrimalRowIDs[prevPrimalRowIds + i - prevNumRows] = _compSolver.rowId(i);
         _nCompPrimalRows++;
      }

      assert(_nPrimalRows == _nCompPrimalRows);


      // looping over all rows from the original problem that were originally contained in the complementary problem.
      // The basic rows will be set as equalities, the non-basic rows will be eliminated from the complementary
      // problem.
      DSVector slackColCoeff(_compSolver.nRows());

      int* rowsforremoval = 0;
      int nrowsforremoval = 0;
      spx_alloc(rowsforremoval, prevPrimalRowIds);
      for( int i = 0; i < prevPrimalRowIds; i++ )
      {
         int rowNumber = _realLP->number(_idsPrimalRowIDs[i]);
         // this loop runs over all rows previously in the complementary problem. If rows are added to the reduced
         // problem, they will be transfered from the incompatible set to the compatible set in the following if
         // statement.
         if( _idsReducedProbRows[rowNumber] )
         {
            // rows added to the reduced problem may have been equality constriants. The equality constraints from the
            // original problem are converted into <= and >= constraints. Upon adding these constraints to the reduced
            // problem, only a single dual column is needed in the complementary problem. Hence, one of the dual columns
            // is removed.
            //

            int solverRowNum = _solver.number(_idsReducedProbRowIDs[rowNumber]);
            assert(solverRowNum >= 0 && solverRowNum < _solver.nRows());
            if( _solver.basis().desc().rowStatus(solverRowNum) == SPxBasis::Desc::P_ON_UPPER ||
                  _solver.basis().desc().rowStatus(solverRowNum) == SPxBasis::Desc::P_FIXED ||
                  _solver.basis().desc().rowStatus(solverRowNum) == SPxBasis::Desc::D_FREE ||
                 (_solver.basis().desc().rowStatus(solverRowNum) == SPxBasis::Desc::D_ON_LOWER &&
                 EQ(_solver.rhs(solverRowNum) - _solver.pVec()[solverRowNum], 0.0, feastol)) )
            {
               //assert(GT(dualVector[solverRowNum], 0.0));
               _compSolver.changeLhs(_idsCompPrimalRowIDs[i], _realLP->rhs(SPxRowId(_idsPrimalRowIDs[i])));
               // need to also update the RHS because a ranged row could have previously been fixed to LOWER
               _compSolver.changeRhs(_idsCompPrimalRowIDs[i], _realLP->rhs(SPxRowId(_idsPrimalRowIDs[i])));
            }
            else if( _solver.basis().desc().rowStatus(solverRowNum) == SPxBasis::Desc::P_ON_LOWER ||
              (_solver.basis().desc().rowStatus(solverRowNum) == SPxBasis::Desc::D_ON_UPPER &&
              EQ(_solver.pVec()[solverRowNum] - _solver.lhs(solverRowNum), 0.0, feastol)) )
            {
               //assert(LT(dualVector[solverRowNum], 0.0));
               _compSolver.changeRhs(_idsCompPrimalRowIDs[i], _realLP->lhs(SPxRowId(_idsPrimalRowIDs[i])));
               // need to also update the LHS because a ranged row could have previously been fixed to UPPER
               _compSolver.changeLhs(_idsCompPrimalRowIDs[i], _realLP->lhs(SPxRowId(_idsPrimalRowIDs[i])));
            }
            else //if ( _solver.basis().desc().rowStatus(solverRowNum) != SPxBasis::Desc::D_FREE )
            {
               //assert(isZero(dualVector[solverRowNum], 0.0));
               rowsforremoval[nrowsforremoval] = _compSolver.number(SPxRowId(_idsCompPrimalRowIDs[i]));
               nrowsforremoval++;

               if( _nElimPrimalRows >= _idsElimPrimalRowIDs.size() )
                  _idsElimPrimalRowIDs.reSize(_realLP->nRows());

               _idsElimPrimalRowIDs[_nElimPrimalRows] = _idsPrimalRowIDs[i];
               _nElimPrimalRows++;
               _idsPrimalRowIDs.remove(i);
               _nPrimalRows--;
               _idsCompPrimalRowIDs.remove(i);
               _nCompPrimalRows--;

               i--;
               prevPrimalRowIds--;
            }
         }
         else
         {
            switch( _compSolver.rowType(_idsCompPrimalRowIDs[i]) )
            {
               case LPRowBase<Real>::RANGE:
                  assert(false);
                  break;
               case LPRowBase<Real>::EQUAL:
                  assert(false);
                  break;
               case LPRowBase<Real>::LESS_EQUAL:
                  slackColCoeff.add(_compSolver.number(SPxRowId(_idsCompPrimalRowIDs[i])), -SLACKCOEFF);
                  break;
               case LPRowBase<Real>::GREATER_EQUAL:
                  slackColCoeff.add(_compSolver.number(SPxRowId(_idsCompPrimalRowIDs[i])), SLACKCOEFF);
                  break;
               default:
                  throw SPxInternalCodeException("XIDSSL01 This should never happen.");
            }

            if( origObj )
            {
               rowsforremoval[nrowsforremoval] = _compSolver.number(SPxRowId(_idsCompPrimalRowIDs[i]));
               nrowsforremoval++;

               if( _nElimPrimalRows >= _idsElimPrimalRowIDs.size() )
                  _idsElimPrimalRowIDs.reSize(_realLP->nRows());

               _idsElimPrimalRowIDs[_nElimPrimalRows] = _idsPrimalRowIDs[i];
               _nElimPrimalRows++;
               _idsPrimalRowIDs.remove(i);
               _nPrimalRows--;
               _idsCompPrimalRowIDs.remove(i);
               _nCompPrimalRows--;

               i--;
               prevPrimalRowIds--;
            }
         }
      }

      // updating the slack column in the complementary problem
      LPColBase<Real> compSlackCol(-1, slackColCoeff, infinity, 0.0);
      _compSolver.changeCol(_compSlackColId, compSlackCol);


      // if the original objective is used, then all complementary rows related to primal rows not in the reduced problem are
      // removed from the complementary problem.
      // As a result, the slack column becomes an empty column.
      int* perm = 0;
      spx_alloc(perm, _compSolver.nRows() + numElimRowsAdded);
      _compSolver.removeRows(rowsforremoval, nrowsforremoval, perm);

      // updating the dual columns to represent the fixed primal variables.
      int* currFixedVars = 0;
      spx_alloc(currFixedVars, _realLP->nCols());
      _identifyComplementaryPrimalFixedPrimalVars(currFixedVars);
      _updateComplementaryPrimalFixedPrimalVars(currFixedVars);


      // freeing allocated memory
      spx_free(currFixedVars);
      spx_free(perm);
      spx_free(rowsforremoval);
   }



   /// checking the optimality of the original problem.
   // this function is called if the complementary problem is solved with a non-negative objective value. This implies
   // that the rows currently included in the reduced problem are sufficient to identify the optimal solution to the
   // original problem.
   void SoPlex::_checkOriginalProblemOptimality(Vector primalVector, bool printViol)
   {
      SSVector x(_solver.nCols());
      x.unSetup();

      _idsTransBasis.coSolve(x, primalVector);

      if( printViol )
      {
         MSG_INFO1(spxout, spxout << std::endl
           << "Checking consistency between the reduced problem and the original problem." << std::endl );
      }


      Real redObjVal = 0;
      Real objectiveVal = 0;
      for( int i = 0; i < _solver.nCols(); i++ )
      {
         redObjVal += _solver.maxObj(i)*primalVector[i];
         objectiveVal += _realLP->maxObj(i)*x[i];
      }

      if( printViol )
      {
         MSG_INFO1(spxout, spxout << "Reduced Problem Objective Value: " << redObjVal << std::endl
           << "Original Problem Objective Value: " << std::endl );
      }

      _solReal._hasPrimal = true;
      _hasSolReal = true;
      // get the primal solutions from the reduced problem
      _solReal._primal.reDim(_solver.nCols());
      _solReal._primal = x;

      Real maxviol = 0;
      Real sumviol = 0;

      if( getIdsBoundViolation(maxviol, sumviol) )
      {
         if( printViol )
            MSG_INFO1(spxout, spxout << "Bound violation - "
               << "Max violation: " << maxviol << " Sum violation: " << sumviol << std::endl );
      }

      _statistics->totalBoundViol = sumviol;
      _statistics->maxBoundViol = maxviol;

      if( getIdsRowViolation(maxviol, sumviol) )
      {
         if( printViol )
            MSG_INFO1(spxout, spxout << "Row violation - "
               << "Max violation: " << maxviol << " Sum violation: " << sumviol << std::endl );
      }

      _statistics->totalRowViol = sumviol;
      _statistics->maxRowViol = maxviol;

      if( printViol )
         MSG_INFO1(spxout, spxout << std::endl );
   }



   /// updating the slack column coefficients to adjust for equality constraints
   void SoPlex::_updateComplementaryDualSlackColCoeff()
   {
      // the slack column for the equality constraints is not handled correctly in the dual conversion. Hence, it is
      // necessary to change the equality coefficients of the dual row related to the slack column.
      for( int i = 0; i < _nPrimalRows; i++ )
      {
         int rowNumber = _realLP->number(SPxRowId(_idsPrimalRowIDs[i]));
         if( !_idsReducedProbRows[rowNumber] )
         {
            if( _realLP->rowType(_idsPrimalRowIDs[i]) == LPRowBase<Real>::EQUAL )
            {
               assert(_realLP->lhs(_idsPrimalRowIDs[i]) == _realLP->rhs(_idsPrimalRowIDs[i]));
               _compSolver.changeLower(_idsDualColIDs[i], 0.0);   // setting the lower bound of the dual column to zero.

               LPColBase<Real> addEqualityCol(-1.0*_realLP->rhs(_idsPrimalRowIDs[i]),
                     -1.0*_compSolver.colVector(_idsDualColIDs[i]), infinity, 0.0);    // adding a new column to the dual

               SPxColId newDualCol;
               _compSolver.addCol(newDualCol, addEqualityCol);

               // inserting the row and col ids for the added column. This is to be next to the original column that has
               // been duplicated.
               _idsPrimalRowIDs.insert(i + 1, 1, _idsPrimalRowIDs[i]);
               _idsDualColIDs.insert(i + 1, 1, newDualCol);
               assert(_realLP->number(_idsPrimalRowIDs[i]) == _realLP->number(_idsPrimalRowIDs[i + 1]));

               i++;
               _nPrimalRows++;
               _nDualCols++;
            }
         }
      }
   }



   /// identify the dual columns related to the fixed variables
   void SoPlex::_identifyComplementaryDualFixedPrimalVars(int* currFixedVars)
   {
      Real feastol = realParam(SoPlex::FEASTOL);

      int numFixedVar = 0;
      for( int i = 0; i < _realLP->nCols(); i++ )
      {
         currFixedVars[i] = 0;
         if( !_idsReducedProbColRowIDs[i].isValid() )
            continue;

         int rowNumber = _solver.number(_idsReducedProbColRowIDs[i]);
         if( _idsReducedProbColRowIDs[i].isValid() )
         {
            if( _solver.basis().desc().rowStatus(rowNumber) == SPxBasis::Desc::P_ON_UPPER ||
               _solver.basis().desc().rowStatus(rowNumber) == SPxBasis::Desc::P_ON_LOWER ||
               _solver.basis().desc().rowStatus(rowNumber) == SPxBasis::Desc::P_FIXED ||
               _solver.basis().desc().rowStatus(rowNumber) == SPxBasis::Desc::D_FREE )
            {
               // setting the value of the _fixedOrigVars array to indicate which variables are at their bounds.
               currFixedVars[i] = getOrigVarFixedDirection(i);

               numFixedVar++;
            }
            else
            {
               // the dual flags do not imply anything about the primal status of the rows.
               if( _solver.basis().desc().rowStatus(rowNumber) == SPxBasis::Desc::D_ON_LOWER &&
                EQ(_solver.rhs(rowNumber) - _solver.pVec()[rowNumber], 0.0, feastol) )
                  currFixedVars[i] = 1;
               else if(  _solver.basis().desc().rowStatus(rowNumber) == SPxBasis::Desc::D_ON_UPPER &&
                EQ(_solver.pVec()[rowNumber] - _solver.lhs(rowNumber), 0.0, feastol) )
                  currFixedVars[i] = -1;
            }
         }
      }

      MSG_INFO3(spxout, spxout << "Number of fixed primal variables in the complementary (dual) problem: "
         << numFixedVar << std::endl );
   }



   /// removing the dual columns related to the fixed variables
   void SoPlex::_removeComplementaryDualFixedPrimalVars(int* currFixedVars)
   {
      SPxColId tempId;
      int ncolsforremoval = 0;
      int* colsforremoval = 0;
      spx_alloc(colsforremoval, _realLP->nCols()*2);

      tempId.inValidate();
      for( int i = 0; i < _realLP->nCols(); i++ )
      {
         assert(_idsCompProbColIDsIdx[i] != -1);   // this should be true in the current implementation
         if( _idsCompProbColIDsIdx[i] != -1 && _fixedOrigVars[i] != -2 )//&& _fixedOrigVars[i] != currFixedVars[i] )
         {
            if( _fixedOrigVars[i] != 0 )
            {
               assert(_compSolver.number(SPxColId(_idsFixedVarDualIDs[i])) >= 0);
               assert(_fixedOrigVars[i] == -1 || _fixedOrigVars[i] == 1);
               assert(_idsFixedVarDualIDs[i].isValid());

               colsforremoval[ncolsforremoval] = _compSolver.number(SPxColId(_idsFixedVarDualIDs[i]));
               ncolsforremoval++;

               _idsFixedVarDualIDs[i] = tempId;
            }
            else //if( false && !_idsReducedProbColRowIDs[i].isValid() ) // we want to remove all valid columns
               // in the current implementation, the only columns not included in the reduced problem are free columns.
            {
               assert((LE(_realLP->lower(i), -infinity) && GE(_realLP->upper(i), infinity)) ||
                     _compSolver.number(SPxColId(_idsVarBoundDualIDs[i*2])) >= 0);
               int varcount = 0;
               if( GT(_realLP->lower(i), -infinity) )
               {
                  colsforremoval[ncolsforremoval] = _compSolver.number(SPxColId(_idsVarBoundDualIDs[i*2 + varcount]));
                  ncolsforremoval++;

                  _idsVarBoundDualIDs[i*2 + varcount] = tempId;
                  varcount++;
               }

               if( LT(_realLP->upper(i), infinity) )
               {
                  colsforremoval[ncolsforremoval] = _compSolver.number(SPxColId(_idsVarBoundDualIDs[i*2 + varcount]));
                  ncolsforremoval++;

                  _idsVarBoundDualIDs[i*2 + varcount] = tempId;
               }

            }
         }
      }

      int* perm = 0;
      spx_alloc(perm, _compSolver.nCols());
      _compSolver.removeCols(colsforremoval, ncolsforremoval, perm);

      // freeing allocated memory
      spx_free(perm);
      spx_free(colsforremoval);
   }



   /// updating the dual columns related to the fixed primal variables.
   void SoPlex::_updateComplementaryDualFixedPrimalVars(int* currFixedVars)
   {
      DSVectorBase<Real> col(1);
      LPColSetBase<Real> boundConsCols;
      LPColSetBase<Real> fixedVarsDualCols(_nPrimalCols);
      int numFixedVars = 0;
      // the solution to the reduced problem results in a number of variables at their bounds. If such variables exist
      // it is necessary to include a dual column to the complementary problem related to a variable fixing. This is
      // equivalent to the tight constraints being converted to equality constraints.
      int numBoundConsCols = 0;
      int* boundConsColsAdded = 0;
      spx_alloc(boundConsColsAdded, _realLP->nCols());
      // NOTE: this loop only goes over the primal columns that are included in the complementary problem, i.e. the
      // columns from the original problem.
      // 29.04.15 in the current implementation, all bound constraints are included in the reduced problem. So, all
      // variables (columns) are included in the reduced problem.
      for( int i = 0; i < _realLP->nCols(); i++ )
      {
         boundConsColsAdded[i] = 0;
         assert(_idsCompProbColIDsIdx[i] != -1);
         if( _idsCompProbColIDsIdx[i] != -1 )
         {
            int idIndex = _idsCompProbColIDsIdx[i];
            assert(_compSolver.number(SPxRowId(_idsDualRowIDs[idIndex])) >= 0);
            col.add(_compSolver.number(SPxRowId(_idsDualRowIDs[idIndex])), 1.0);
            if( currFixedVars[i] != 0 )
            {
               assert(currFixedVars[i] == -1 || currFixedVars[i] == 1);

               assert(_realLP->lower(i) == _solver.lhs(_idsReducedProbColRowIDs[i]));
               assert(_realLP->upper(i) == _solver.rhs(_idsReducedProbColRowIDs[i]));
               Real colObjCoeff = 0;
               if( currFixedVars[i] == -1 )
                  colObjCoeff = _solver.lhs(_idsReducedProbColRowIDs[i]);
               else
                  colObjCoeff = _solver.rhs(_idsReducedProbColRowIDs[i]);

               fixedVarsDualCols.add(colObjCoeff, -infinity, col, infinity);
               numFixedVars++;
            }
            // 09.02.15 I think that the else should only be entered if the column does not exist in the reduced
            // prob. I have tested by just leaving this as an else (without the if), but I think that this is wrong.
            //else if( !_idsReducedProbColRowIDs[i].isValid() )

            // NOTE: in the current implementation all columns, except the free columns, are included int the reduced
            // problem. There in no need to include the variable bounds in the complementary problem.
            else //if( false && _fixedOrigVars[i] == -2 )
            {
               bool isRedProbCol = _idsReducedProbColRowIDs[i].isValid();
               // 29.04.15 in the current implementation only free variables are not included in the reduced problem
               if( GT(_realLP->lower(i), -infinity) )
               {
                  if( !isRedProbCol )
                     col.add(_compSolver.number(SPxRowId(_compSlackDualRowId)), -SLACKCOEFF);
                  boundConsCols.add(_realLP->lower(i), Real(-infinity), col, 0.0);

                  if( !isRedProbCol )
                     col.remove(col.size() - 1);

                  boundConsColsAdded[i]++;
                  numBoundConsCols++;
               }

               if( LT(_realLP->upper(i), infinity) )
               {
                  if( !isRedProbCol )
                     col.add(_compSolver.number(SPxRowId(_compSlackDualRowId)), SLACKCOEFF);
                  boundConsCols.add(_realLP->upper(i), 0.0, col, Real(infinity));

                  if( !isRedProbCol )
                     col.remove(col.size() - 1);

                  boundConsColsAdded[i]++;
                  numBoundConsCols++;
               }
            }
            col.clear();
            _fixedOrigVars[i] = currFixedVars[i];
         }
      }

      // adding the fixed var dual columns to the complementary problem
      SPxColId* addedcolids = 0;
      spx_alloc(addedcolids, numFixedVars);
      _compSolver.addCols(addedcolids, fixedVarsDualCols);

      SPxColId tempId;
      int addedcolcount = 0;

      tempId.inValidate();
      for( int i = 0; i < _realLP->nCols(); i++ )
      {
         if( _fixedOrigVars[i] != 0 )
         {
            assert(_fixedOrigVars[i] == -1 || _fixedOrigVars[i] == 1);
            _idsFixedVarDualIDs[i] = addedcolids[addedcolcount];
            addedcolcount++;
         }
         else
            _idsFixedVarDualIDs[i] = tempId;
      }

      // adding the bound cons dual columns to the complementary problem
      SPxColId* addedbndcolids = 0;
      spx_alloc(addedbndcolids, numBoundConsCols);
      _compSolver.addCols(addedbndcolids, boundConsCols);

      addedcolcount = 0;
      for( int i = 0; i < _realLP->nCols(); i++ )
      {
         if( boundConsColsAdded[i] > 0 )
         {
            for( int j = 0; j < boundConsColsAdded[i]; j++ )
            {
               _idsVarBoundDualIDs[i*2 + j] = addedbndcolids[addedcolcount];
               addedcolcount++;
            }
         }

         switch( boundConsColsAdded[i] )
         {
            case 0:
               _idsVarBoundDualIDs[i*2] = tempId;
            case 1:
               _idsVarBoundDualIDs[i*2 + 1] = tempId;
               break;
         }
      }

      // freeing allocated memory
      spx_free(addedbndcolids);
      spx_free(addedcolids);
      spx_free(boundConsColsAdded);
   }


   /// identify the dual columns related to the fixed variables
   void SoPlex::_identifyComplementaryPrimalFixedPrimalVars(int* currFixedVars)
   {
      int numFixedVar = 0;
      for( int i = 0; i < _realLP->nCols(); i++ )
      {
         currFixedVars[i] = 0;
         if( !_idsReducedProbColRowIDs[i].isValid() )
            continue;

         int rowNumber = _solver.number(_idsReducedProbColRowIDs[i]);
         if( _idsReducedProbColRowIDs[i].isValid() &&
              (_solver.basis().desc().rowStatus(rowNumber) == SPxBasis::Desc::P_ON_UPPER ||
               _solver.basis().desc().rowStatus(rowNumber) == SPxBasis::Desc::P_ON_LOWER ||
               _solver.basis().desc().rowStatus(rowNumber) == SPxBasis::Desc::P_FIXED) )
         {
            // setting the value of the _fixedOrigVars array to indicate which variables are at their bounds.
            currFixedVars[i] = getOrigVarFixedDirection(i);

            numFixedVar++;
         }
      }

      MSG_INFO3(spxout, spxout << "Number of fixed primal variables in the complementary (primal) problem: "
         << numFixedVar << std::endl );
   }



   /// updating the dual columns related to the fixed primal variables.
   void SoPlex::_updateComplementaryPrimalFixedPrimalVars(int* currFixedVars)
   {
      int numFixedVars = 0;
      // NOTE: this loop only goes over the primal columns that are included in the complementary problem, i.e. the
      // columns from the original problem.
      // 29.04.15 in the current implementation, all bound constraints are included in the reduced problem. So, all
      // variables (columns) are included in the reduced problem.
      for( int i = 0; i < _nCompPrimalCols; i++ )
      {
         int colNumber = _compSolver.number(SPxColId(_idsCompPrimalColIDs[i]));
         if( _fixedOrigVars[colNumber] != currFixedVars[colNumber] )
         {
            if( currFixedVars[colNumber] != 0 )
            {
               assert(currFixedVars[colNumber] == -1 || currFixedVars[colNumber] == 1);

               if( currFixedVars[colNumber] == -1 )
                  _compSolver.changeBounds(colNumber, _realLP->lower(SPxColId(_idsPrimalColIDs[i])),
                     _realLP->lower(SPxColId(_idsPrimalColIDs[i])));
               else
                  _compSolver.changeBounds(colNumber, _realLP->upper(SPxColId(_idsPrimalColIDs[i])),
                     _realLP->upper(SPxColId(_idsPrimalColIDs[i])));

               numFixedVars++;
            }
            else
            {
               _compSolver.changeBounds(colNumber, -infinity, infinity);
            }
         }

         _fixedOrigVars[colNumber] = currFixedVars[colNumber];
      }
   }



   /// updating the complementary dual problem with the original objective function
   void SoPlex::_setComplementaryDualOriginalObjective()
   {
      for( int i = 0; i < _realLP->nCols(); i++ )
      {
         assert(_idsCompProbColIDsIdx[i] != -1);   // this should be true in the current implementation
         int idIndex = _idsCompProbColIDsIdx[i];
         int compRowNumber = _compSolver.number(_idsDualRowIDs[idIndex]);
         if( _idsCompProbColIDsIdx[i] != -1 )
         {
            // In the dual conversion, when a variable has a non-standard bound it is converted to a free variable.
            if( LE(_realLP->lower(i), -infinity) && GE(_realLP->upper(i), infinity) )
            {
               // unrestricted variable
               _compSolver.changeLhs(compRowNumber, _realLP->obj(i));
               _compSolver.changeRhs(compRowNumber, _realLP->obj(i));
               assert(LE(_compSolver.lhs(compRowNumber), _compSolver.rhs(compRowNumber)));
            }
            else if( LE(_realLP->lower(i), -infinity) )
            {
               // variable with a finite upper bound
               _compSolver.changeRhs(compRowNumber, _realLP->obj(i));
               if( isZero(_realLP->upper(i)) )
                  _compSolver.changeLhs(compRowNumber, -infinity);
               else
                  _compSolver.changeLhs(compRowNumber, _realLP->obj(i));
            }
            else if( GE(_realLP->upper(i), infinity) )
            {
               // variable with a finite lower bound
               _compSolver.changeLhs(compRowNumber, _realLP->obj(i));
               if( isZero(_realLP->upper(i)) )
                  _compSolver.changeRhs(compRowNumber, infinity);
               else
                  _compSolver.changeRhs(compRowNumber, _realLP->obj(i));
            }
            else if( NE(_realLP->lower(i), _realLP->upper(i)) )
            {
               // variable with a finite upper and lower bound
               if( isZero(_realLP->upper(i)) )
               {
                  _compSolver.changeLhs(compRowNumber, _realLP->obj(i));
                  _compSolver.changeRhs(compRowNumber, infinity);
               }
               else if( isZero(_realLP->upper(i)) )
               {
                  _compSolver.changeLhs(compRowNumber, -infinity);
                  _compSolver.changeRhs(compRowNumber, _realLP->obj(i));
               }
               else
               {
                  _compSolver.changeLhs(compRowNumber, _realLP->obj(i));
                  _compSolver.changeRhs(compRowNumber, _realLP->obj(i));
               }
            }
            else
            {
               // fixed variable
               _compSolver.changeLhs(compRowNumber, _realLP->obj(i));
               _compSolver.changeRhs(compRowNumber, _realLP->obj(i));
            }
         }
      }

      // removing the complementary problem slack column dual row
      _compSolver.removeRow(_compSlackDualRowId);
   }



   /// updating the complementary primal problem with the original objective function
   void SoPlex::_setComplementaryPrimalOriginalObjective()
   {
      // the comp solver has not removed any columns. Only the slack variables have been added.
      assert(_realLP->nCols() == _compSolver.nCols() - 1);
      for( int i = 0; i < _realLP->nCols(); i++ )
      {
         int colNumber = _realLP->number(_idsPrimalColIDs[i]);
         int compColNumber = _compSolver.number(_idsCompPrimalColIDs[i]);
         _compSolver.changeObj(compColNumber, _realLP->maxObj(colNumber));
      }

      // removing the complementary problem slack column dual row
      _compSolver.removeCol(_compSlackColId);
   }



   /// determining which bound the primal variables will be fixed to.
   int SoPlex::getOrigVarFixedDirection(int colNum)
   {
      if( !_idsReducedProbColRowIDs[colNum].isValid() )
         return 0;

      int rowNumber = _solver.number(_idsReducedProbColRowIDs[colNum]);
      // setting the value of the _fixedOrigVars array to indicate which variables are at their bounds.
      if( _solver.basis().desc().rowStatus(rowNumber) == SPxBasis::Desc::P_ON_UPPER ||
       _solver.basis().desc().rowStatus(rowNumber) == SPxBasis::Desc::P_FIXED ||
       _solver.basis().desc().rowStatus(rowNumber) == SPxBasis::Desc::D_FREE )
      {
         assert(_solver.rhs(rowNumber) < infinity);
         return 1;
      }
      else if( _solver.basis().desc().rowStatus(rowNumber) == SPxBasis::Desc::P_ON_LOWER )
      {
         assert(_solver.lhs(rowNumber) > -infinity);
         return -1;
      }

      return 0;
   }



   // @todo update this function and related comments. It has only been hacked together.
   /// checks result of the solving process and solves again without preprocessing if necessary
   // @todo need to evaluate the solution to ensure that it is solved to optimality and then we are able to perform the
   // next steps in the algorithm.
   void SoPlex::_evaluateSolutionIDS(SPxSolver& solver, SLUFactor& sluFactor, SPxSimplifier::Result result)
   {
      SPxSolver::Status solverStat = SPxSolver::UNKNOWN;
      if( result == SPxSimplifier::INFEASIBLE )
         solverStat = SPxSolver::INFEASIBLE;
      else if( result == SPxSimplifier::DUAL_INFEASIBLE )
         solverStat = SPxSolver::INForUNBD;
      else if( result == SPxSimplifier::UNBOUNDED )
         solverStat = SPxSolver::UNBOUNDED;
      else if( result == SPxSimplifier::VANISHED )
         solverStat = SPxSolver::OPTIMAL;
      else if( result == SPxSimplifier::OKAY )
         solverStat = solver.status();

      // updating the status of SoPlex if the problem solved is the reduced problem.
      if( _currentProb == IDS_ORIG || _currentProb == IDS_RED )
         _status = solverStat;

      // process result
      switch( solverStat )
      {
         case SPxSolver::OPTIMAL:
            if( !_isRealLPLoaded )
            {
               solver.changeObjOffset(realParam(SoPlex::OBJ_OFFSET));
               _idsResolveWithoutPreprocessing(solver, sluFactor, result);
               // Need to solve the complementary problem
               return;
            }
            else
               _hasBasis = true;

            break;

         case SPxSolver::UNBOUNDED:
         case SPxSolver::INFEASIBLE:
         case SPxSolver::INForUNBD:
            // in case of infeasibility or unboundedness, we currently can not unsimplify, but have to solve the original LP again
            if( !_isRealLPLoaded )
            {
               solver.changeObjOffset(realParam(SoPlex::OBJ_OFFSET));
               _idsSimplifyAndSolve(solver, sluFactor, false, false);
               return;
            }
            else
               _hasBasis = (solver.basis().status() > SPxBasis::NO_PROBLEM);
            break;

         case SPxSolver::ABORT_DECOMP:
         case SPxSolver::ABORT_EXDECOMP:
            // in the initialisation of the decomposition simplex, we want to keep the current basis.
            if( !_isRealLPLoaded )
            {
               solver.changeObjOffset(realParam(SoPlex::OBJ_OFFSET));
               _idsResolveWithoutPreprocessing(solver, sluFactor, result);
               //_idsSimplifyAndSolve(solver, sluFactor, false, false);
               return;
            }
            else
               _hasBasis = (solver.basis().status() > SPxBasis::NO_PROBLEM);
            break;

         case SPxSolver::SINGULAR:
         case SPxSolver::ABORT_CYCLING:
            // in case of infeasibility or unboundedness, we currently can not unsimplify, but have to solve the original LP again
            if( !_isRealLPLoaded )
            {
               solver.changeObjOffset(realParam(SoPlex::OBJ_OFFSET));
               _idsSimplifyAndSolve(solver, sluFactor, false, false);
               return;
            }
            break;


            // else fallthrough
         case SPxSolver::ABORT_TIME:
         case SPxSolver::ABORT_ITER:
         case SPxSolver::ABORT_VALUE:
         case SPxSolver::REGULAR:
         case SPxSolver::RUNNING:
            // store regular basis if there is no simplifier and the original problem is not in the solver because of
            // scaling; non-optimal bases should currently not be unsimplified
            if( _simplifier == 0 && solver.basis().status() > SPxBasis::NO_PROBLEM )
            {
               _basisStatusRows.reSize(_idsLP->nRows());
               _basisStatusCols.reSize(_idsLP->nCols());
               assert(_basisStatusRows.size() == solver.nRows());
               assert(_basisStatusCols.size() == solver.nCols());

               solver.getBasis(_basisStatusRows.get_ptr(), _basisStatusCols.get_ptr());
               _hasBasis = true;
            }
            else
               _hasBasis = false;
            break;

         default:
            _hasBasis = false;
            break;
      }
   }




   /// checks the dual feasibility of the current basis
   bool SoPlex::checkBasisDualFeasibility(Vector feasVec)
   {
      assert(_solver.rep() == SPxSolver::ROW);
      assert(_solver.spxSense() == SPxLPBase<Real>::MAXIMIZE);

      Real feastol = realParam(SoPlex::FEASTOL);
      // iterating through all elements in the basis to determine whether the dual feasibility is satisfied.
      for( int i = 0; i < _solver.nCols(); i++ )
      {
         if( _solver.basis().baseId(i).isSPxRowId() ) // find the row id's for rows in the basis
         {
            int rownumber = _solver.number(_solver.basis().baseId(i));
            if( _solver.basis().desc().rowStatus(rownumber) != SPxBasis::Desc::P_ON_UPPER &&
                  _solver.basis().desc().rowStatus(rownumber) != SPxBasis::Desc::P_FIXED )
            {
               if( GT(feasVec[i], 0, feastol) )
                  return false;
            }

            if( _solver.basis().desc().rowStatus(rownumber) != SPxBasis::Desc::P_ON_LOWER &&
                  _solver.basis().desc().rowStatus(rownumber) != SPxBasis::Desc::P_FIXED )
            {
               if( LT(feasVec[i], 0, feastol) )
                  return false;
            }

         }
         else if( _solver.basis().baseId(i).isSPxColId() )  // get the column id's for the columns in the basis
         {
            int colnumber = _solver.number(_solver.basis().baseId(i));
            if( _solver.basis().desc().colStatus(colnumber) != SPxBasis::Desc::P_ON_UPPER &&
                  _solver.basis().desc().colStatus(colnumber) != SPxBasis::Desc::P_FIXED )
            {
               if( GT(feasVec[i], 0, feastol) )
                  return false;
            }

            if( _solver.basis().desc().colStatus(colnumber) != SPxBasis::Desc::P_ON_LOWER &&
                  _solver.basis().desc().colStatus(colnumber) != SPxBasis::Desc::P_FIXED )
            {
               if( LT(feasVec[i], 0, feastol) )
                  return false;
            }
         }
      }

      return true;
   }



   /// returns the expected sign of the dual variables for the reduced problem
   SoPlex::DualSign SoPlex::getExpectedDualVariableSign(int rowNumber)
   {
      if( _solver.isRowBasic(rowNumber) )
      {
         if( _solver.basis().desc().rowStatus(rowNumber) != SPxBasis::Desc::P_ON_UPPER &&
               _solver.basis().desc().rowStatus(rowNumber) != SPxBasis::Desc::P_FIXED )
            return SoPlex::IS_NEG;

         if( _solver.basis().desc().rowStatus(rowNumber) != SPxBasis::Desc::P_ON_LOWER &&
               _solver.basis().desc().rowStatus(rowNumber) != SPxBasis::Desc::P_FIXED )
            return SoPlex::IS_POS;
      }

      return SoPlex::IS_FREE;
   }



   /// returns the expected sign of the dual variables for the original problem
   SoPlex::DualSign SoPlex::getOrigProbDualVariableSign(int rowNumber)
   {
      if( _realLP->rowType(rowNumber) == LPRowBase<Real>::LESS_EQUAL )
         return IS_POS;

      if( _realLP->rowType(rowNumber) == LPRowBase<Real>::GREATER_EQUAL )
         return IS_NEG;

      if( _realLP->rowType(rowNumber) == LPRowBase<Real>::EQUAL )
         return IS_FREE;

      // this final statement needs to be checked.
      if( _realLP->rowType(rowNumber) == LPRowBase<Real>::RANGE )
         return IS_FREE;

      return IS_FREE;
   }


   /// print display line of flying table
   void SoPlex::printIdsDisplayLine(SPxSolver& solver, const SPxOut::Verbosity origVerb, bool force, bool forceHead)
   {
      // setting the verbosity level
      const SPxOut::Verbosity currVerb = spxout.getVerbosity();
      spxout.setVerbosity( origVerb );

      int displayFreq = solver.getDisplayFreq();

      MSG_INFO1( spxout,
         if( forceHead || (_idsDisplayLine % (displayFreq*30) == 0) )
         {
            spxout << "type |   time |   iters | red iter | alg iter |     rows |     cols |  shift   |    value\n";
         }
         if( force || (_idsDisplayLine % displayFreq == 0) )
         {
            Real currentTime = _statistics->solvingTime->time();
            (solver.type() == SPxSolver::LEAVE) ? spxout << "  L  |" : spxout << "  E  |";
            spxout << std::fixed << std::setw(7) << std::setprecision(1) << currentTime << " |";
            spxout << std::scientific << std::setprecision(2);
            spxout << std::setw(8) << _statistics->iterations << " | ";
            spxout << std::scientific << std::setprecision(2);
            spxout << std::setw(8) << solver.iterations() << " | ";
            spxout << std::scientific << std::setprecision(2);
            spxout << std::setw(8) << _statistics->callsReducedProb << " | ";
            spxout << std::scientific << std::setprecision(2);
            spxout << std::setw(8) << numIncludedRows << " | ";
            spxout << std::scientific << std::setprecision(2);
            spxout << std::setw(8) << solver.nCols() << " | "
            << solver.shift() << " | "
            << std::setprecision(8) << solver.value() + solver.objOffset()
            << std::endl;

         }
         _idsDisplayLine++;
      );

      spxout.setVerbosity( currVerb );
   }



   /// stores the problem statistics of the original problem
   void SoPlex::getOriginalProblemStatistics()
   {
      numProbRows = _realLP->nRows();
      numProbCols = _realLP->nCols();
      numNonzeros = _realLP->nNzos();
      minAbsNonzero = _realLP->minAbsNzo();
      maxAbsNonzero = _realLP->maxAbsNzo();

      origCountLower = 0;
      origCountUpper = 0;
      origCountBoxed = 0;
      origCountFreeCol = 0;

      origCountLhs = 0;
      origCountRhs = 0;
      origCountRanged = 0;
      origCountFreeRow = 0;

      for( int i = 0; i < _realLP->nCols(); i++ )
      {
         bool hasLower = false;
         bool hasUpper = false;

         if( _realLP->lower(i) > -infinity )
         {
            origCountLower++;
            hasLower = true;
         }

         if( _realLP->upper(i) < infinity )
         {
            origCountUpper++;
            hasUpper = true;
         }

         if( hasUpper && hasLower )
            origCountBoxed++;

         if( !hasUpper && !hasLower )
            origCountFreeCol++;
      }

      for( int i = 0; i < _realLP->nRows(); i++)
      {
         bool hasRhs = false;
         bool hasLhs = false;

         if( _realLP->lhs(i) > -infinity )
         {
            origCountLhs++;
            hasLhs = true;
         }

         if( _realLP->rhs(i) < infinity )
         {
            origCountRhs++;
            hasRhs = true;
         }

         if( hasRhs && hasLhs )
            origCountRanged++;

         if( !hasRhs && !hasLhs )
            origCountFreeRow++;
      }
   }


   void SoPlex::printOriginalProblemStatistics(std::ostream& os)
   {
      os << "  Columns           : " << numProbCols << "\n"
         << "              boxed : " << origCountBoxed << "\n"
         << "        lower bound : " << origCountLower << "\n"
         << "        upper bound : " << origCountUpper << "\n"
         << "               free : " << origCountFreeCol << "\n"
         << "  Rows              : " << numProbRows << "\n"
         << "             ranged : " << origCountRanged << "\n"
         << "                lhs : " << origCountLhs << "\n"
         << "                rhs : " << origCountRhs << "\n"
         << "               free : " << origCountFreeRow << "\n"
         << "  Nonzeros          : " << numNonzeros << "\n"
         << "         per column : " << Real(numNonzeros) / Real(numProbCols) << "\n"
         << "            per row : " << Real(numNonzeros) / Real(numProbRows) << "\n"
         << "           sparsity : " << Real(numNonzeros) / Real(numProbCols) / Real(numProbRows) << "\n"
         << "    min. abs. value : " << Real(minAbsNonzero) << "\n"
         << "    max. abs. value : " << Real(maxAbsNonzero) << "\n";
   }



   /// gets the coefficient of the slack variable in the primal complementary problem
   Real SoPlex::getCompSlackVarCoeff(int primalRowNum)
   {
      int indDir = 1;
      switch( _realLP->rowType(_idsPrimalRowIDs[primalRowNum]) )
      {
         // NOTE: check the sign of the slackCoeff for the Range constraints. This will depend on the method of
         // dual conversion.
         case LPRowBase<Real>::RANGE:
            assert((primalRowNum < _nPrimalRows - 1 && _realLP->number(SPxColId(_idsPrimalRowIDs[primalRowNum])) ==
               _realLP->number(SPxColId(_idsPrimalRowIDs[primalRowNum+1]))) ||
               (primalRowNum > 0 && _realLP->number(SPxColId(_idsPrimalRowIDs[primalRowNum-1])) ==
               _realLP->number(SPxColId(_idsPrimalRowIDs[primalRowNum]))));

            // determine with primalRowNum and primalRowNum+1 or primalRowNum-1 and primalRowNum have the same row id.
            if( _realLP->number(SPxColId(_idsPrimalRowIDs[primalRowNum-1])) ==
               _realLP->number(SPxColId(_idsPrimalRowIDs[primalRowNum])) )
               indDir = -1;

            if( _compSolver.obj(_compSolver.number(SPxColId(_idsDualColIDs[primalRowNum]))) <
               _compSolver.obj(_compSolver.number(SPxColId(_idsDualColIDs[primalRowNum + indDir]))))
               return -SLACKCOEFF;
            else
               return SLACKCOEFF;

            break;

         case LPRowBase<Real>::GREATER_EQUAL:
            return -SLACKCOEFF;
            break;
         case LPRowBase<Real>::EQUAL:
            assert((primalRowNum < _nPrimalRows - 1 && _realLP->number(SPxColId(_idsPrimalRowIDs[primalRowNum])) ==
               _realLP->number(SPxColId(_idsPrimalRowIDs[primalRowNum+1]))) ||
               (primalRowNum > 0 && _realLP->number(SPxColId(_idsPrimalRowIDs[primalRowNum-1])) ==
               _realLP->number(SPxColId(_idsPrimalRowIDs[primalRowNum]))));
         case LPRowBase<Real>::LESS_EQUAL:
            return SLACKCOEFF;
            break;
         default:
            throw SPxInternalCodeException("XIDSSL01 This should never happen.");
      }

      return 0;
   }



   /// gets violation of bounds; returns true on success
   bool SoPlex::getIdsBoundViolation(Real& maxviol, Real& sumviol)
   {
      Real feastol = realParam(SoPlex::FEASTOL);

      VectorReal& primal = _solReal._primal;
      assert(primal.dim() == _realLP->nCols());

      _nIdsViolBounds = 0;

      maxviol = 0.0;
      sumviol = 0.0;

      bool isViol = false;
      bool isMaxViol = false;

      for( int i = _realLP->nCols() - 1; i >= 0; i-- )
      {
         Real viol = _realLP->lower(i) - primal[i];

         isViol = false;
         isMaxViol = false;

         if( viol > 0.0 )
         {
            sumviol += viol;
            if( viol > maxviol )
            {
               maxviol = viol;
               isMaxViol = true;
            }
         }

         if( GT(viol, 0.0, feastol) )
            isViol = true;

         viol = primal[i] - _realLP->upper(i);
         if( viol > 0.0 )
         {
            sumviol += viol;
            if( viol > maxviol )
            {
               maxviol = viol;
               isMaxViol = true;
            }
         }

         if( GT(viol, 0.0, feastol) )
            isViol = true;

         if( isViol )
         {

            // updating the violated bounds list
            if( isMaxViol )
            {
               _idsViolatedBounds[_nIdsViolBounds] = _idsViolatedBounds[0];
               _idsViolatedBounds[0] = i;
            }
            else
               _idsViolatedBounds[_nIdsViolBounds] = i;

            _nIdsViolBounds++;
         }
      }

      return true;
   }



   /// gets violation of constraints; returns true on success
   bool SoPlex::getIdsRowViolation(Real& maxviol, Real& sumviol)
   {
      Real feastol = realParam(SoPlex::FEASTOL);

      VectorReal& primal = _solReal._primal;
      assert(primal.dim() == _realLP->nCols());

      DVectorReal activity(_realLP->nRows());
      _realLP->computePrimalActivity(primal, activity);

      _nIdsViolRows = 0;

      maxviol = 0.0;
      sumviol = 0.0;

      bool isViol = false;
      bool isMaxViol = false;

      for( int i = _realLP->nRows() - 1; i >= 0; i-- )
      {
         Real currviol = 0.0;

         isViol = false;
         isMaxViol = false;

         Real viol = _realLP->lhs(i) - activity[i];
         if( viol > 0.0 )
         {
            sumviol += viol;
            if( viol > maxviol )
            {
               maxviol = viol;
               isMaxViol = true;
            }

            if( viol > currviol )
               currviol = viol;
         }

         if( GT(viol, 0.0, feastol) )
            isViol = true;

         viol = activity[i] - _realLP->rhs(i);
         if( viol > 0.0 )
         {
            sumviol += viol;
            if( viol > maxviol )
            {
               maxviol = viol;
               isMaxViol = true;
            }

            if( viol > currviol )
               currviol = viol;
         }

         if( GT(viol, 0.0, feastol) )
            isViol = true;

         if( isViol )
         {
            // updating the violated rows list
            if( isMaxViol )
            {
               _idsViolatedRows[_nIdsViolRows] = _idsViolatedRows[0];
               _idsViolatedRows[0] = i;
            }
            else
               _idsViolatedRows[_nIdsViolRows] = i;

            _nIdsViolRows++;
         }
      }

      return true;
   }



   /// function call to terminate the decomposition simplex
   bool SoPlex::idsTerminate(Real timeLimit)
   {
      Real maxTime = timeLimit;

      // check if a time limit is actually set
      if( maxTime < 0 || maxTime >= infinity )
         return false;

      Real currentTime = _statistics->solvingTime->time();
      if ( currentTime >= maxTime )
      {
         MSG_INFO2( spxout, spxout << " --- timelimit (" << _solver.getMaxTime()
            << ") reached" << std::endl; )
         _solver.setSolverStatus(SPxSolver::ABORT_TIME);
         return true;
      }

      return false;
   }



   /// function to retrieve the original problem row basis status from the reduced and complementary problems
   void SoPlex::getOriginalProblemBasisRowStatus(DataArray< int >& degenerateRowNums,
      DataArray< SPxSolver::VarStatus >& degenerateRowStatus, int& nDegenerateRows, int& nNonBasicRows)
   {
      Real feastol = realParam(SoPlex::FEASTOL);
      int basicRow = 0;

      // looping over all rows not contained in the complementary problem
      for( int i = 0; i < _nElimPrimalRows; i++ )
      {
         int rowNumber = _realLP->number(_idsElimPrimalRowIDs[i]);

         int solverRowNum = _solver.number(_idsReducedProbRowIDs[rowNumber]);
         assert(solverRowNum >= 0 && solverRowNum < _solver.nRows());
         if( _solver.basis().desc().rowStatus(solverRowNum) == SPxBasis::Desc::P_ON_UPPER ) /*||
            (_solver.basis().desc().rowStatus(solverRowNum) == SPxBasis::Desc::D_ON_LOWER &&
             EQ(_solver.rhs(solverRowNum) - _solver.pVec()[solverRowNum], 0.0, feastol)) )*/
         {
            _basisStatusRows[rowNumber] = SPxSolver::ON_UPPER;
         }
         else if( _solver.basis().desc().rowStatus(solverRowNum) == SPxBasis::Desc::P_ON_LOWER ) /*||
            (_solver.basis().desc().rowStatus(solverRowNum) == SPxBasis::Desc::D_ON_UPPER &&
             EQ(_solver.pVec()[solverRowNum] - _solver.lhs(solverRowNum), 0.0, feastol)) )*/
         {
            _basisStatusRows[rowNumber] = SPxSolver::ON_LOWER;
         }
         else if( _solver.basis().desc().rowStatus(solverRowNum) == SPxBasis::Desc::P_FIXED )
         {
            _basisStatusRows[rowNumber] = SPxSolver::FIXED;
         }
         else if( _solver.basis().desc().rowStatus(solverRowNum) == SPxBasis::Desc::P_FREE )
         {
            _basisStatusRows[rowNumber] = SPxSolver::ZERO;
         }
         else
         {
            _basisStatusRows[rowNumber] = SPxSolver::BASIC;
            basicRow++;
         }
      }

      for( int i = 0; i < _nPrimalRows; i++ )
      {
         int rowNumber = _realLP->number(_idsPrimalRowIDs[i]);
         // this loop runs over all rows previously in the complementary problem.
         if( _idsReducedProbRows[rowNumber] )
         {
            int solverRowNum = _solver.number(_idsReducedProbRowIDs[rowNumber]);
            assert(solverRowNum >= 0 && solverRowNum < _solver.nRows());
            if( _solver.basis().desc().rowStatus(solverRowNum) == SPxBasis::Desc::P_ON_UPPER )
            {
               _basisStatusRows[rowNumber] = SPxSolver::ON_UPPER;
            }
            else if( _solver.basis().desc().rowStatus(solverRowNum) == SPxBasis::Desc::D_ON_LOWER &&
                EQ(_solver.rhs(solverRowNum) - _solver.pVec()[solverRowNum], 0.0, feastol) )
            {
               // if a row is non basic, but is at its bound then the row number and status is stored
               degenerateRowNums[nDegenerateRows] = rowNumber;
               degenerateRowStatus[nDegenerateRows] = SPxSolver::ON_UPPER;
               nDegenerateRows++;
            }
            else if( _solver.basis().desc().rowStatus(solverRowNum) == SPxBasis::Desc::P_ON_LOWER )
            {
               _basisStatusRows[rowNumber] = SPxSolver::ON_LOWER;
            }
            else if( _solver.basis().desc().rowStatus(solverRowNum) == SPxBasis::Desc::D_ON_UPPER &&
                EQ(_solver.pVec()[solverRowNum] - _solver.lhs(solverRowNum), 0.0, feastol) )
            {
               // if a row is non basic, but is at its bound then the row number and status is stored
               degenerateRowNums[nDegenerateRows] = rowNumber;
               degenerateRowStatus[nDegenerateRows] = SPxSolver::ON_LOWER;
               nDegenerateRows++;
            }
            else if( _solver.basis().desc().rowStatus(solverRowNum) == SPxBasis::Desc::P_FIXED )
            {
               _basisStatusRows[rowNumber] = SPxSolver::FIXED;
            }
            else if( _solver.basis().desc().rowStatus(solverRowNum) == SPxBasis::Desc::P_FREE )
            {
               _basisStatusRows[rowNumber] = SPxSolver::ZERO;
            }
            else
            {
               _basisStatusRows[rowNumber] = SPxSolver::BASIC;
               basicRow++;
            }
         }
         else
         {
            _basisStatusRows[rowNumber] = SPxSolver::BASIC;
            basicRow++;
            switch( _realLP->rowType(_idsPrimalRowIDs[i]) )
            {
               case LPRowBase<Real>::RANGE:
                  //assert(_realLP->number(SPxColId(_idsPrimalRowIDs[i])) ==
                        //_realLP->number(SPxColId(_idsPrimalRowIDs[i+1])));
                  //assert(_compSolver.obj(_compSolver.number(SPxColId(_idsDualColIDs[i]))) !=
                        //_compSolver.obj(_compSolver.number(SPxColId(_idsDualColIDs[i + 1]))));

                  //// NOTE: This is not correct at present. Need to modify the inequalities. Look at the dual conversion
                  //// function to get this correct.
                  //if( _compSolver.obj(_compSolver.number(SPxColId(_idsDualColIDs[i]))) <
                     //_compSolver.obj(_compSolver.number(SPxColId(_idsDualColIDs[i + 1]))))
                  //{
                     //if( _compSolver.basis().desc().rowStatus(_compSolver.number(SPxColId(_idsDualColIDs[i]))) ==
                        //SPxBasis::Desc::D_ON_LOWER )
                     //{
                        //_basisStatusRows[rowNumber] = SPxSolver::ON_LOWER;
                     //}
                     //else if( _compSolver.basis().desc().rowStatus(_compSolver.number(SPxColId(_idsDualColIDs[i + 1]))) ==
                        //SPxBasis::Desc::D_ON_UPPER )
                     //{
                        //_basisStatusRows[rowNumber] = SPxSolver::ON_UPPER;
                     //}
                     //else
                     //{
                        //_basisStatusRows[rowNumber] = SPxSolver::BASIC;
                        //basicRow++;
                     //}
                  //}
                  //else
                  //{
                     //if( _compSolver.basis().desc().rowStatus(_compSolver.number(SPxColId(_idsDualColIDs[i]))) ==
                        //SPxBasis::Desc::D_ON_UPPER )
                     //{
                        //_basisStatusRows[rowNumber] = SPxSolver::ON_UPPER;
                     //}
                     //else if( _compSolver.basis().desc().rowStatus(_compSolver.number(SPxColId(_idsDualColIDs[i + 1]))) ==
                        //SPxBasis::Desc::D_ON_LOWER )
                     //{
                        //_basisStatusRows[rowNumber] = SPxSolver::ON_LOWER;
                     //}
                     //else
                     //{
                        //_basisStatusRows[rowNumber] = SPxSolver::BASIC;
                        //basicRow++;
                     //}
                  //}

                  i++;
                  break;
               case LPRowBase<Real>::EQUAL:
                  //assert(_realLP->number(SPxColId(_idsPrimalRowIDs[i])) ==
                        //_realLP->number(SPxColId(_idsPrimalRowIDs[i+1])));

                  //_basisStatusRows[rowNumber] = SPxSolver::FIXED;

                  i++;
                  break;
               case LPRowBase<Real>::GREATER_EQUAL:
                  //if( _compSolver.basis().desc().rowStatus(_compSolver.number(SPxColId(_idsDualColIDs[i]))) ==
                     //SPxBasis::Desc::D_ON_LOWER )
                  //{
                     //_basisStatusRows[rowNumber] = SPxSolver::ON_LOWER;
                  //}
                  //else
                  //{
                     //_basisStatusRows[rowNumber] = SPxSolver::BASIC;
                     //basicRow++;
                  //}

                  break;
               case LPRowBase<Real>::LESS_EQUAL:
                  //if( _compSolver.basis().desc().rowStatus(_compSolver.number(SPxColId(_idsDualColIDs[i]))) ==
                     //SPxBasis::Desc::D_ON_UPPER )
                  //{
                     //_basisStatusRows[rowNumber] = SPxSolver::ON_UPPER;
                  //}
                  //else
                  //{
                     //_basisStatusRows[rowNumber] = SPxSolver::BASIC;
                     //basicRow++;
                  //}

                  break;
               default:
                  throw SPxInternalCodeException("XIDSSL01 This should never happen.");
            }
         }
      }

      nNonBasicRows = _realLP->nRows() - basicRow - nDegenerateRows;
      MSG_INFO2(spxout, spxout << "Number of non-basic rows: " << nNonBasicRows << " (from "
         << _realLP->nRows() << ")" << std::endl );
   }


   /// function to retrieve the column status for the original problem basis from the reduced and complementary problems
   // all columns are currently in the reduced problem, so it is only necessary to retrieve the status from there.
   // However, a transformation of the columns has been made, so it is only possible to retrieve the status from the
   // variable bound constraints.
   void SoPlex::getOriginalProblemBasisColStatus(int& nNonBasicCols)
   {
      Real feastol = realParam(SoPlex::FEASTOL);
      int basicCol = 0;
      int numZeroDual = 0;

      // looping over all variable bound constraints
      for( int i = 0; i < _realLP->nCols(); i++ )
      {
         if( !_idsReducedProbColRowIDs[i].isValid() )
            continue;

         int rowNumber = _solver.number(_idsReducedProbColRowIDs[i]);
         if( _idsReducedProbColRowIDs[i].isValid() )
         {
            if( _solver.basis().desc().rowStatus(rowNumber) == SPxBasis::Desc::P_ON_UPPER ) /*||
              (_solver.basis().desc().rowStatus(rowNumber) == SPxBasis::Desc::D_ON_LOWER &&
                EQ(_solver.rhs(rowNumber) - _solver.pVec()[rowNumber], 0.0, feastol)) )*/
            {
               _basisStatusCols[i] = SPxSolver::ON_UPPER;
            }
            else if( _solver.basis().desc().rowStatus(rowNumber) == SPxBasis::Desc::P_ON_LOWER ) /*||
              (_solver.basis().desc().rowStatus(rowNumber) == SPxBasis::Desc::D_ON_UPPER &&
                EQ(_solver.pVec()[rowNumber] - _solver.lhs(rowNumber), 0.0, feastol)) )*/
            {
               _basisStatusCols[i] = SPxSolver::ON_LOWER;
            }
            else if( _solver.basis().desc().rowStatus(rowNumber) == SPxBasis::Desc::P_FIXED )
            {
               _basisStatusCols[i] = SPxSolver::FIXED;
            }
            else if( _solver.basis().desc().rowStatus(rowNumber) == SPxBasis::Desc::P_FREE )
            {
               _basisStatusCols[i] = SPxSolver::ZERO;
            }
            else
            {
               _basisStatusCols[i] = SPxSolver::BASIC;
               basicCol++;
            }
         }

         if( EQ(_solver.fVec()[i], 0.0, feastol) )
            numZeroDual++;
      }

      nNonBasicCols = _realLP->nCols() - basicCol;
      MSG_INFO2(spxout, spxout << "Number of non-basic columns: "
         << nNonBasicCols << " (from " << _realLP->nCols() << ")" << std::endl
         << "Number of zero dual columns: " << numZeroDual << " (from " << _realLP->nCols() << ")" << std::endl );
   }



   /// function to build a basis for the original problem as given by the solution to the reduced problem
   void SoPlex::_writeOriginalProblemBasis(const char* filename, NameSet* rowNames, NameSet* colNames, bool cpxFormat)
   {
      int numrows = _realLP->nRows();
      int numcols = _realLP->nCols();
      int nNonBasicRows = 0;
      int nNonBasicCols = 0;
      int nDegenerateRows = 0;
      DataArray< int > degenerateRowNums;   // array to store the orig prob row number of degenerate rows
      DataArray< SPxSolver::VarStatus > degenerateRowStatus;  // possible status for the degenerate rows in the orig prob

      // setting the basis status rows and columns to size of the original problem
      _basisStatusRows.reSize(numrows);
      _basisStatusCols.reSize(numcols);

      degenerateRowNums.reSize(numrows);
      degenerateRowStatus.reSize(numrows);

      // setting the row and column status for the original problem
      getOriginalProblemBasisRowStatus(degenerateRowNums, degenerateRowStatus, nDegenerateRows, nNonBasicRows);
      getOriginalProblemBasisColStatus(nNonBasicCols);

      // checking the consistency of the non-basic rows and columns for printing out the basis.
      // it is necessary to have numcols of non-basic rows and columns.
      // if there are not enought non-basic rows and columns, then the degenerate rows are set as non-basic.
      // all degenerate rows that are not changed to non-basic must be set to basic.
      assert(nDegenerateRows + nNonBasicRows + nNonBasicCols >= numcols);
      int degenRowsSetNonBasic = 0;
      for( int i = 0; i < nDegenerateRows; i++ )
      {
         if( nNonBasicRows + nNonBasicCols + degenRowsSetNonBasic < numcols )
         {
            _basisStatusRows[degenerateRowNums[i]] = degenerateRowStatus[i];
            degenRowsSetNonBasic++;
         }
         else
            _basisStatusRows[degenerateRowNums[i]] = SPxSolver::BASIC;
      }

      // writing the basis file for the original problem
      bool wasRealLPLoaded = _isRealLPLoaded;
      _isRealLPLoaded = false;
      writeBasisFile(filename, rowNames, colNames, cpxFormat);
      _isRealLPLoaded = wasRealLPLoaded;
   }



   /// function call to terminate the decomposition simplex
   void SoPlex::_getIdsPrimalVector(Vector& primal)
   {
   }

} // namespace soplex
#endif
