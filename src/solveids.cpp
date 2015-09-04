/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
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

#define SIMPLEDEBUG
//#define DEBUGGING
//#define SOLVEIDS_DEBUG
//#define WRITE_LP
//#define EXTRA_PRINT

#define SLACKCOEFF   1.0   /**< the coefficient of the slack variable in the incompatible rows. */

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
      int stopCount = 0;

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

      //#ifdef WRITE_LP
      printf("Writing the original lp to a file\n");
      _solver.writeFile("original.lp");
      //#endif

      // it is necessary to solve the initial problem to find a starting basis
      _solver.setIdsStatus(SPxSolver::FINDSTARTBASIS);

      Real degeneracyLevel = 0;
      _idsFeasVector.reDim(_solver.nCols());
      // since the original LP may have been shifted, the dual multiplier will not be correct for the original LP. This
      // loop will recheck the degeneracy level and compute the proper dual multipliers.
      do
      {
         _idsSimplifyAndSolve(_solver, _slufactor, false, false);
         _solver.basis().solve(_idsFeasVector, _solver.maxObj());
         degeneracyLevel = _solver.getDegeneracyLevel(_idsFeasVector);
         if( _solver.type() == SPxSolver::LEAVE || _solver.status() >= SPxSolver::OPTIMAL
               || _solver.getIdsStatus() == SPxSolver::DONTFINDSTARTBASIS )
         {
            // returning the sense to minimise
            if( intParam(SoPlex::OBJSENSE) == SoPlex::OBJSENSE_MINIMIZE )
            {
               assert(_solver.spxSense() == SPxLPBase<Real>::MAXIMIZE);

               _solver.changeObj(-(_solver.maxObj()));
               _solver.changeSense(SPxLPBase<Real>::MINIMIZE);
            }

            _solver.setIdsStatus(SPxSolver::DONTFINDSTARTBASIS);

            // resolving the problem to update the real lp and solve with the correct objective.
            _idsSimplifyAndSolve(_solver, _slufactor, false, false);

            // storing the solution from the reduced problem
            _storeSolutionReal();

            // stop timing
            _statistics->solvingTime->stop();
            return;
         }
      } while( (degeneracyLevel > 0.9 || degeneracyLevel < 0.1) || !checkBasisDualFeasibility(_idsFeasVector) );

      // updating the algorithm iterations statistic
      _statistics->callsReducedProb++;

      _solver.setIdsStatus(SPxSolver::DONTFINDSTARTBASIS);
      _hasBasis = true; // this is probably wrong. Will need to set this properly.

      // setting the statistic information
      numIdsIter = 0;
      numRedProbIter = _solver.iterations();
      numCompProbIter = 0;

      // setting the display information
      _idsDisplayLine = 0;

      MSG_INFO1( spxout,
         spxout << "========      Degeneracy Detected       ========" << std::endl;
         spxout << std::endl;
         spxout << "======== Commencing decomposition solve ========" << std::endl;
         );

      // setting the verbosity level
      const SPxOut::Verbosity orig_verbosity = spxout.getVerbosity();
      spxout.setVerbosity( SPxOut::ERROR );

      // creating copies of the original problem that will be manipulated to form the reduced and complementary
      // problems.
      _createIdsReducedAndComplementaryProblems();

      // creating the initial reduced problem from the basis information
      _formIdsReducedProblem();


      _hasBasis = false;
      bool hasRedBasis = false;

      int algIterCount = 0;
      while( !stop || stopCount <= 0 )
      {
         //#ifdef WRITE_LP
      printf("Writing the reduced lp to a file\n");
      _solver.writeFile("reduced.lp");
      //#endif

#ifdef SOLVEIDS_DEBUG
         printf("Reduced Prob Rows: ");
         for( int j = 0; j < numRowsReal(); j++ )
         {
            if( _idsReducedProbRows[j] )
               printf("%d ", j);
         }
         printf("\n");
#endif

         // setting the current solving mode.
         _currentProb = IDS_RED;

         // solve the reduced problem
#ifdef SOLVEIDS_DEBUG
         printf("\n");
         printf("=========================\n");
         printf("Solving: Reduced Problem.\n");
         printf("=========================\n");
         printf("\n");
#endif
         _hasBasis = hasRedBasis;
         _idsSimplifyAndSolve(_solver, _slufactor, !algIterCount, !algIterCount);

         stop = idsTerminate();  // checking whether the algorithm should terminate

         // printing display line
#ifdef SIMPLEDEBUG
         printIdsDisplayLine(_solver, orig_verbosity, true, !algIterCount);
#else
         printIdsDisplayLine(_solver, orig_verbosity, !algIterCount, !algIterCount);
#endif

         // updating the algorithm iterations statistics
         _statistics->callsReducedProb++;


         assert(_isRealLPLoaded);
         hasRedBasis = _hasBasis;
         _hasBasis = false;

         if( _solver.status() > SPxSolver::OPTIMAL )
            break;

         _idsFeasVector.reDim(_solver.nCols());
         _solver.basis().solve(_idsFeasVector, _solver.maxObj());

         // get the dual solutions from the reduced problem
         DVector reducedLPDualVector(_solver.nRows());
         DVector reducedLPRedcostVector(_solver.nCols());
         _solver.getDual(reducedLPDualVector);
         _solver.getRedCost(reducedLPRedcostVector);

         if( algIterCount == 0 )
            _formIdsComplementaryProblem();  // create the complementary problem using
                                             // the solution to the reduced problem
         else
            _updateIdsComplementaryProblem(false);

#ifdef WRITE_LP
         printf("Writing the complementary lp to a file\n");
         _compSolver.writeFile("complement_dual.lp");
#endif
         // setting the current solving mode.
         _currentProb = IDS_COMP;

         // solve the complementary problem
#ifdef SOLVEIDS_DEBUG
         printf("\n");
         printf("===============================\n");
         printf("Solving: Complementary Problem.\n");
         printf("===============================\n");
         printf("\n");
#endif
         _idsSimplifyAndSolve(_compSolver, _compSlufactor, true, true);
         assert(_isRealLPLoaded);

         //#ifdef SIMPLEDEBUG
         printf("Iteration %d Objective Value: %f\n", algIterCount, _compSolver.objValue());
         //#endif

         //SPxLPIds compPrimalLP;
         //_compSolver.buildDualProblem(compPrimalLP);

         //char buffer[50];
         //sprintf(buffer, "lpfiles/complement_%d.lp", algIterCount);
         //printf("Writing the complementary primal lp to a file\n");
         //compPrimalLP.writeFile(buffer);

         // check the optimality of the original problem with the objective value of the complementary problem
         if( GE(_compSolver.objValue(), 0.0, 1e-10) || _compSolver.status() == SPxSolver::UNBOUNDED )
         {
            if( _compSolver.status() == SPxSolver::UNBOUNDED )
               printf("Unbounded complementary problem.\n");
            stop = true;
         }

         if( stop )
            stopCount++;

         if( !stop )
         {
            // get the primal solutions from the complementary problem
            DVector compLPPrimalVector(_compSolver.nCols());
            _compSolver.getPrimal(compLPPrimalVector);

            _updateIdsReducedProblem(_compSolver.objValue(), reducedLPDualVector, reducedLPRedcostVector,
                  compLPPrimalVector);
         }

         numIdsIter++;
         algIterCount++;
      }

      // computing the solution for the original variables
      DVector reducedLPPrimalVector(_solver.nCols());
      _solver.getPrimal(reducedLPPrimalVector);
      _checkOriginalProblemOptimality(reducedLPPrimalVector);

      // writing the original problem basis
      _writeOriginalProblemBasis("orig_basis.bas", NULL, NULL, false);

      // resetting the verbosity level
      spxout.setVerbosity( orig_verbosity );

#if 0
      // removing all rows not in the reduced problem
      _updateIdsComplementaryProblem(true);
#endif
      // solving the complementary problem with the original objective function
      _setComplementaryOriginalObjective();

      // for clean up
      // the real lp has to be set to not loaded.
      //_isRealLPLoaded = false;
      // the basis has to be set to false
      _hasBasis = false;

      SPxLPIds compDualLP;
      _compSolver.buildDualProblem(compDualLP, _idsPrimalRowIDs.get_ptr(), _idsPrimalColIDs.get_ptr(),
            _idsDualRowIDs.get_ptr(), _idsDualColIDs.get_ptr(), &_nPrimalRows, &_nPrimalCols, &_nDualRows, &_nDualCols);

      _compSolver.loadLP(compDualLP);
      //#ifdef WRITE_LP
      printf("Writing the complementary lp to a file\n");
      _compSolver.writeFile("complement.lp");

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
         printf("Bound violation: %.20f %.20f\n", maxviol, sumviol);

      if( getIdsRowViolation(maxviol, sumviol) )
         printf("Row violation: %.20f %.20f\n", maxviol, sumviol);

      printf("Objective Value: %f\n", _compSolver.objValue());
      //#endif

      // returning the sense to minimise
      if( intParam(SoPlex::OBJSENSE) == SoPlex::OBJSENSE_MINIMIZE )
      {
         assert(_solver.spxSense() == SPxLPBase<Real>::MAXIMIZE);

         _solver.changeObj(-(_solver.maxObj()));
         _solver.changeSense(SPxLPBase<Real>::MINIMIZE);

         // Need to add commands to multiply the objective solution values by -1
      }

      // resetting the verbosity level
      spxout.setVerbosity( orig_verbosity );

      MSG_INFO1( spxout,
         spxout << "========  Decomposition solve completed ========" << std::endl;
         spxout << std::endl;
         spxout << "========   Resolving original problem   ========" << std::endl;
         );

      spx_free(_idsCompProbColIDsIdx);
      spx_free(_fixedOrigVars);
      spx_free(_idsReducedProbCols);
      spx_free(_idsReducedProbRows);

      // retreiving the original problem statistics prior to destroying it.
      getOriginalProblemStatistics();

      // setting the reduced problem statistics
      _statistics->numRedProbRows = numIncludedRows;
      _statistics->numRedProbCols = _solver.nCols();

      // resolving the problem to update the real lp and solve with the correct objective.
      _realLP->~SPxLPReal();
      spx_free(_realLP);
      _realLP = &_solver;
      _isRealLPLoaded = true;

      // printing display line for resolve of problem
      _solver.printDisplayLine(false, true);

      _preprocessAndSolveReal(false);

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
      //_compSolver.reLoad();
      _compSolver.setOutstream(spxout);
      _compSolver.setSolver(&_compSlufactor);
      //_compSolver.loadLP(*_realLP);
      //_compSolver.loadBasis(_solver.basis().desc());

      //_idsSimplifyAndSolve(_compSolver, _compSlufactor, true, true);
   }



   /// forms the reduced problem
   void SoPlex::_formIdsReducedProblem()
   {
#ifdef SOLVEIDS_DEBUG
      printf("Forming the reduced problem\n");
#endif
      int* bind = 0;
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

      _idsLP = 0;
      spx_alloc(_idsLP);
      _idsLP = new (_idsLP) SPxLPIds(_solver);

      // retreiving the basis information
      _basisStatusRows.reSize(numRowsReal());
      _basisStatusCols.reSize(numColsReal());
      _solver.getBasis(_basisStatusRows.get_ptr(), _basisStatusCols.get_ptr());

      // get thie indices of the rows with positive dual multipliers and columns with positive reduced costs.
      spx_alloc(bind, numColsReal());

      spx_alloc(nonposind, numColsReal());
      spx_alloc(colsforremoval, numColsReal());
      _getNonPositiveDualMultiplierInds(_idsFeasVector, nonposind, bind, colsforremoval, &nnonposind);

      // get the compatible columns from the constraint matrix w.r.t the current basis matrix
#ifdef SOLVEIDS_DEBUG
      printf("Computing the compatible columns\n");
      printf("Solving time: %f\n", solveTime());
#endif
      spx_alloc(compatind, _solver.nRows());
      spx_alloc(rowsforremoval, _solver.nRows());
      _getCompatibleColumns(nonposind, compatind, rowsforremoval, colsforremoval, nnonposind, &ncompatind, true);

      int* compatboundcons = 0;
      int ncompatboundcons = 0;
      spx_alloc(compatboundcons, numColsReal());

      LPRowSet boundcons;

      // identifying the compatible bound constraints
      _getCompatibleBoundCons(boundcons, compatboundcons, nonposind, &ncompatboundcons, nnonposind);

      // delete rows and columns from the LP to form the reduced problem
#ifdef SOLVEIDS_DEBUG
      printf("Deleting rows and columns to form the reduced problem\n");
      printf("Solving time: %f\n", solveTime());
#endif

      // computing the reduced problem obj coefficient vector and updating the problem
      _computeReducedProbObjCoeff();

#ifdef SOLVEIDS_DEBUG
      printf("Nonposinds obj coefficients\n");
      for( int i = 0; i < nnonposind; i++ )
      {
         if( !isZero(_idsLP->maxObj(nonposind[i])) )
            printf("nonposind[%d]: %d, coeff: %f\n", i, nonposind[i], _idsLP->maxObj(nonposind[i]));
      }
#endif

      _solver.loadLP(*_idsLP);

      // the colsforremoval are the columns with a zero reduced cost.
      // the rowsforremoval are the rows identified as incompatible.
      _solver.removeRows(rowsforremoval);

      // adding the rows for the compatible bound constraints
      SPxRowId* addedrowids = 0;
      spx_alloc(addedrowids, ncompatboundcons);
      _solver.addRows(addedrowids, boundcons);

      for( int i = 0; i < ncompatboundcons; i++ )
         _idsReducedProbColRowIDs[compatboundcons[i]] = addedrowids[i];


      // freeing allocated memory
      spx_free(addedrowids);
      spx_free(compatboundcons);
      spx_free(rowsforremoval);
      spx_free(compatind);
      spx_free(colsforremoval);
      spx_free(nonposind);
      spx_free(bind);

      _idsLP->~SPxLPReal();
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

      _deleteAndUpdateRowsComplementaryProblem();

      // initialising the arrays to store the row id's from the primal and the col id's from the dual
      _idsPrimalRowIDs.reSize(3*_compSolver.nRows());
      _idsPrimalColIDs.reSize(3*_compSolver.nCols());
      _idsDualRowIDs.reSize(3*_compSolver.nCols());
      _idsDualColIDs.reSize(3*_compSolver.nRows());
      _nPrimalRows = 0;
      _nPrimalCols = 0;
      _nDualRows = 0;
      _nDualCols = 0;

#ifdef WRITE_LP
      printf("Writing the complementary primal lp to a file\n");
      _compSolver.writeFile("comp_primal.lp");
#endif

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

      _updateComplementarySlackColCoeff();

      // initalising the array for the dual columns representing the fixed variables in the complementary problem.
      _idsFixedVarDualIDs.reSize(_realLP->nCols());
      _idsVarBoundDualIDs.reSize(2*_realLP->nCols());


      // updating the constraints for the complementary problem
      _updateIdsComplementaryProblem(false);


#ifdef WRITE_LP
      printf("Writing the complementary primal lp to a file\n");
      _compSolver.writeFile("comp_dual.lp");
#endif
   }



   /// simplifies the problem and solves
   void SoPlex::_idsSimplifyAndSolve(SPxSolver& solver, SLUFactor& sluFactor, bool fromScratch, bool applyPreprocessing)
   {
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
         solver.changeObjOffset(_simplifier->getObjoffset());
      }

      _statistics->preprocessingTime->stop();

      // run the simplex method if problem has not been solved by the simplifier
      if( result == SPxSimplifier::OKAY )
      {
         if( _scaler != 0 )
            _scaler->scale(solver);

         bool _hadBasis = _hasBasis;

         _statistics->simplexTime->start();
         solver.solve();
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
         || (result == SPxSimplifier::OKAY && solver.status() == SPxSolver::OPTIMAL));

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
            assert(solver.status() == SPxSolver::OPTIMAL);

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
         DVector compPrimalVector)
   {
      Real feastol = realParam(SoPlex::FEASTOL);

      Real maxDualRatio = infinity;

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
            compProbPrimal = compPrimalVector[_compSolver.number(SPxColId(_idsDualColIDs[i]))]; // this is u

            // the variable in the basis is degenerate.
            if( EQ(reducedProbDual, 0.0, feastol) )
            {
               MSG_WARNING( spxout,
                 spxout << "WIMDSM01: reduced problem dual value is very close to zero." << std::endl; );
               continue;
            }

#ifdef SOLVEIDS_DEBUG
            printf("%d redProbNum: %d compProbNum: %d origProbNum: %d. y = %.10f u = %.20f, rowtype: %d, DualSign: %d\n",
                  i, solverRowNum, _compSolver.number(SPxColId(_idsDualColIDs[i])),
                  _realLP->number(SPxRowId(_idsPrimalRowIDs[i])), reducedProbDual, compProbPrimal,
                  _solver.rowType(solverRowNum), getExpectedDualVariableSign(solverRowNum));
#endif

            // the translation of the complementary primal problem to the dual some rows resulted in two columns.
            if( i < _nPrimalRows - 1 &&
                  _realLP->number(SPxRowId(_idsPrimalRowIDs[i])) == _realLP->number(SPxRowId(_idsPrimalRowIDs[i + 1])) )
            {
               i++;
               // @todo make sure that this is just a simple sum
               compProbPrimal += compPrimalVector[_compSolver.number(SPxColId(_idsDualColIDs[i]))]; // this is u
            }



            // updating the ratio
#if 0
            if( (LE(reducedProbDual, 0, 1e-10) && LE(compProbPrimal, 0, 1e-10)) ||
                  (GE(reducedProbDual, 0, 1e-10) && GE(compProbPrimal, 0, 1e-10)) )
#endif
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
            compProbPrimal = compPrimalVector[_compSolver.number(SPxColId(_idsDualColIDs[i]))]; // this is v
#ifdef SOLVEIDS_DEBUG
            printf("%d compProbNum: %d origProbNum: %d. v = %.20f\n", i, _compSolver.number(SPxColId(_idsDualColIDs[i])),
                  _realLP->number(SPxRowId(_idsPrimalRowIDs[i])), compProbPrimal);
#endif
         }

      }

#ifdef SOLVEIDS_DEBUG
      printf("maxDualRatio: %f\n", maxDualRatio);
#endif

#if 1
      for( int i = 0; i < _nPrimalCols; i++ )
      {
         Real reducedProbDual = 0;
         Real compProbPrimal = 0;
         Real dualRatio = 0;
         int colNumber = _realLP->number(SPxColId(_idsPrimalColIDs[i])); // Need to check this. It should probably be
                                                                         // just colNumber = i.
         assert(colNumber >= 0 || (colNumber == -1 && _idsPrimalColIDs[i].getIdx() == _compSlackColId.getIdx()));
         if( colNumber >= 0 && _fixedOrigVars[colNumber] != 0 )
         {
            assert(_fixedOrigVars[colNumber] == -1 || _fixedOrigVars[colNumber] == 1);
            assert(_idsReducedProbColRowIDs[colNumber].isValid());
            assert(_compSolver.number(_idsFixedVarDualIDs[colNumber]) >= 0);

            int solverRowNum = _solver.number(_idsReducedProbColRowIDs[colNumber]);

            reducedProbDual = dualVector[solverRowNum]; // this is y
            compProbPrimal = compPrimalVector[_compSolver.number(_idsFixedVarDualIDs[colNumber])]; // this is u

            // the variable in the basis is degenerate.
            if( EQ(reducedProbDual, 0.0, feastol) )
            {
               MSG_WARNING( spxout,
                 spxout << "WIMDSM02: reduced problem dual value is very close to zero." << std::endl; );
               continue;
            }

#ifdef SOLVEIDS_DEBUG
            printf("%d %d %d: y = %.20f u = %.20f, DualSign: %d\n", i, colNumber, _fixedOrigVars[colNumber],
                  reducedProbDual, compProbPrimal, getExpectedDualVariableSign(solverRowNum));
#endif

            assert(_realLP->number(SPxColId(_idsPrimalColIDs[i])) != _realLP->number(SPxColId(_idsPrimalColIDs[i + 1])));

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
      }

#ifdef SOLVEIDS_DEBUG
      printf("maxDualRatio: %f\n", maxDualRatio);
#endif

#endif

      // setting the maxDualRatio to a maximum of 1
      if( maxDualRatio > 1.0 )
         maxDualRatio = 1.0;

      DVector compProbRedcost(_compSolver.nCols());   // the reduced costs of the complementary problem

      // Retrieving the reduced costs for each row.
      _compSolver.getRedCost(compProbRedcost);

      LPRowSet updaterows;

      int* newrowidx = 0;
      int nnewrowidx = 0;
      double largestViolation = 0;
      int bestrow = -1;
      spx_alloc(newrowidx, _nPrimalRows);

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
            int compRowNumber = _compSolver.number(SPxColId(_idsDualColIDs[i]));
            // retreiving the complementary problem primal solutions
            compProbPrimal = compPrimalVector[compRowNumber]; // this is v
            compRowRedcost = compProbRedcost[compRowNumber];

            // the translation of the complementary primal problem to the dual some rows resulted in two columns.
            if( i < _nPrimalRows - 1 &&
                  _realLP->number(SPxRowId(_idsPrimalRowIDs[i])) == _realLP->number(SPxRowId(_idsPrimalRowIDs[i + 1])) )
            {
               i++;
               compRowNumber = _compSolver.number(SPxColId(_idsDualColIDs[i]));
               // @todo make sure that this is just a simple sum
               compProbPrimal += compPrimalVector[compRowNumber]; // this is v
               compRowRedcost += compProbRedcost[compRowNumber];
            }

            SoPlex::DualSign varSign = getOrigProbDualVariableSign(rownumber);

#ifdef SOLVEIDS_DEBUG
            printf("i: %d, rownumber: %d, compProbPrimal: %f, dualSign: %d\n", i, rownumber, compProbPrimal, varSign);
#endif
            // add row to the reduced problem if the computed dual is of the correct sign for a feasible dual solution
            if( ratioTest && ((varSign == SoPlex::IS_FREE && !isZero(compProbPrimal, feastol)) ||
                  (varSign == SoPlex::IS_POS && GT(compProbPrimal*maxDualRatio, 0, feastol)) ||
                  (varSign == SoPlex::IS_NEG && LT(compProbPrimal*maxDualRatio, 0, feastol)))
                  /*|| (//isZero(objValue, 1e-10) &&
                   * _compSolver.basis().desc().rowStatus(compRowNum) > SPxBasis::Desc::D_FREE &&
                     _compSolver.basis().desc().rowStatus(compRowNum) < SPxBasis::Desc::D_UNDEFINED)*/ )
            {
#ifdef SOLVEIDS_DEBUG
               printf("Adding rows!!!\n");
#endif
               //this set of statements are required to add rows to the reduced problem. This will add a row for every
               //violated constraint. I only want to add a row for the most violated constraint. That is why the row
               //adding functionality is put outside the for loop.
#if 1
               updaterows.add(_transformedRows.lhs(rownumber), _transformedRows.rowVector(rownumber),
                     _transformedRows.rhs(rownumber));

               if( !_idsReducedProbRows[rownumber] )
               {
                  numIncludedRows++;
                  assert(numIncludedRows <= _realLP->nRows());
               }
               _idsReducedProbRows[rownumber] = true;
               newrowidx[nnewrowidx] = rownumber;
               nnewrowidx++;
#else
               if( GT(abs(compProbPrimal*maxDualRatio), largestViolation) )
               {
                  largestViolation = abs(compProbPrimal*maxDualRatio);
                  bestrow = rownumber;
               }
#endif
            }
         }
      }

#if 0
      // only adding the most violated row.
      updaterows.add(_transformedRows.lhs(bestrow), _transformedRows.rowVector(bestrow),
            _transformedRows.rhs(bestrow));

      _idsReducedProbRows[bestrow] = true;
      newrowidx[nnewrowidx] = bestrow;
      nnewrowidx++;
#endif

#if 0
      for( int i = 0; i < _nPrimalCols; i++ )
      {
         LPRowReal origlprow;
         DSVectorBase<Real> rowtoaddVec(_realLP->nCols());
         Real compProbPrimal = 0;
         Real compRowRedcost = 0;
         int rownumber = _realLP->number(SPxRowId(_idsPrimalRowIDs[i]));
         if( !_idsReducedProbColRowIDs[i].isValid() )
         {
            compRowRedcost = compProbRedcost[compRowNumber];
         }
      }
#endif

      // if no rows are identified by the pricing rule, we add rows based upon the constraint violations
      if( !ratioTest || nnewrowidx == 0 )
      {
         _findViolatedRows(objValue, updaterows, newrowidx, nnewrowidx);
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



#define LARGEST_VIOL
   /// builds the update rows with those violated in the complmentary problem
   // A row is violated in the constraint matrix Ax <= b, if b - A_{i}x < 0
   // To aid the computation, all violations are translated to <= constraints
   void SoPlex::_findViolatedRows(Real compObjValue, LPRowSet& updaterows, int* newrowidx, int& nnewrowidx)
   {
      Real feastol = realParam(SoPlex::FEASTOL);
      DVector compProbRedcost(_compSolver.nCols());   // the reduced costs of the complementary problem

      // Retrieving the slacks for each row.
      _compSolver.getRedCost(compProbRedcost);

#ifdef LARGEST_VIOL
      double largestViolation = 0;
      int bestrow = -1;
#endif
      for( int i = 0; i < _nPrimalRows; i++ )
      {
         LPRowReal origlprow;
         DSVectorBase<Real> rowtoaddVec(_realLP->nCols());
         Real compProbViol = 0;
         Real compSlackCoeff = 0;
         int rownumber = _realLP->number(SPxRowId(_idsPrimalRowIDs[i]));
         if( !_idsReducedProbRows[rownumber] )
         {
            // retreiving the violation of the complementary problem primal constraints
            compSlackCoeff = getCompSlackVarCoeff(i);
            compProbViol = compProbRedcost[_compSolver.number(SPxColId(_idsDualColIDs[i]))]; // this is b - Ax
            // subtracting the slack variable value
            compProbViol += compObjValue*compSlackCoeff; // must add on the slack variable value.
            compProbViol *= compSlackCoeff;  // translating the violation to a <= constraint

            // NOTE: if the row was originally a ranged constraint, we are only interest in one of the inequalities.
            // If one inequality of the range violates the bounds, then we will add the row.

            // the translation of the complementary primal problem to the dual some rows resulted in two columns.
            if( i < _nPrimalRows - 1 &&
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

            if( GE(compObjValue, 0.0, feastol) )
            {
               printf("Rowtype: %d, Constraints violation: %.20f\n", _realLP->rowType(rownumber), compProbViol);
               //if( _realLP->rowType(rownumber) == LPRowBase<Real>::LESS_EQUAL )
                  //compProbViol = -1;
            }

#ifdef SOLVEIDS_DEBUG
            printf("i: %d, rownumber: %d, compProbViol: %f\n", i, rownumber, compProbViol);
#endif
            if( LT(compProbViol, 0, feastol) )
            {
#ifdef SOLVEIDS_DEBUG
               printf("Adding rows!!!\n");
#endif
#ifndef LARGEST_VIOL
               updaterows.add(_transformedRows.lhs(rownumber), _transformedRows.rowVector(rownumber),
                     _transformedRows.rhs(rownumber));

               if( !_idsReducedProbRows[rownumber] )
               {
                  numIncludedRows++;
                  assert(numIncludedRows <= _realLP->nRows());
               }
               _idsReducedProbRows[rownumber] = true;
               newrowidx[nnewrowidx] = rownumber;
               nnewrowidx++;
#else
               if( LT(compProbViol, largestViolation) )
               {
                  largestViolation = compProbViol;
                  bestrow = rownumber;
               }
#endif
            }
         }
      }
#ifdef LARGEST_VIOL
      // only adding the most violated row.
      updaterows.add(_transformedRows.lhs(bestrow), _transformedRows.rowVector(bestrow),
            _transformedRows.rhs(bestrow));

      assert(!_idsReducedProbRows[bestrow]);
      if( !_idsReducedProbRows[bestrow] )
      {
         numIncludedRows++;
         assert(numIncludedRows <= _realLP->nRows());
      }
      _idsReducedProbRows[bestrow] = true;
      newrowidx[nnewrowidx] = bestrow;
      nnewrowidx++;
#endif
   }



   // This function assumes that the basis is in the row form.
   // @todo extend this to the case when the basis is in the column form.
   //
   // NOTE: Changing "nonposind[*nnonposind] = bind[i]" to "nonposind[*nnonposind] = i"
   void SoPlex::_getNonPositiveDualMultiplierInds(Vector feasVector, int* nonposind, int* bind, int* colsforremoval,
         int* nnonposind)
   {
      assert(_solver.rep() == SPxSolver::ROW);

      Real feastol = realParam(SoPlex::FEASTOL);

      bool delCol;

      _idsReducedProbColIDs.reSize(numColsReal());
      *nnonposind = 0;

      // iterating over all columns in the basis matrix
      // this identifies the basis indices and the indicies that are positive.
      for( int i = 0; i < _solver.nCols(); ++i ) // @todo Check the use of numColsReal for the reduced problem.
      {
#ifdef SOLVEIDS_DEBUG
         printf("%d %d %f ", i, _solver.number(_solver.basis().baseId(i)), feasVector[i]);
#endif
         _idsReducedProbCols[i] = true;
         _idsReducedProbColIDs[i].inValidate();
         colsforremoval[i] = i;
         delCol = false;
         // @todo I have questions about my implementation of this function. I don't think that I am getting the right
         // information. Additionally, in getCompatibleColumns the information may not be used correctly.
         if( _solver.basis().baseId(i).isSPxRowId() ) // find the row id's for rows in the basis
         {
#ifdef SOLVEIDS_DEBUG
            printf("status: %d, Row\n", int(_solver.basis().desc().rowStatus(_solver.number(_solver.basis().baseId(i)))));
#endif
            bind[i] = -1 - _realLP->number(SPxRowId(_solver.basis().baseId(i))); // getting the corresponding row
                                                                                   // for the original LP.

            //@todo need to check this regarding min and max problems
            if( isZero(feasVector[i], feastol) )
            {
               nonposind[*nnonposind] = i;
               (*nnonposind)++;

               // NOTE: commenting out the delCol flag at this point. The colsforremoval array should indicate the
               // columns that have a zero reduced cost. Hence, the delCol flag should only be set in the isSPxColId
               // branch of the if statement.
               delCol = true;
            }
         }
         else if( _solver.basis().baseId(i).isSPxColId() )  // get the column id's for the columns in the basis
         {
#ifdef SOLVEIDS_DEBUG
            printf("status: %d, Col\n", int(_solver.basis().desc().colStatus(_solver.number(_solver.basis().baseId(i)))));
#endif
            bind[i] = _realLP->number(SPxColId(_solver.basis().baseId(i)));

            //int colnumber = _solver.number(_solver.basis().baseId(i));
            if( isZero(feasVector[i], feastol)/* &&
                   _solver.basis().desc().colStatus(colnumber) != SPxBasis::Desc::P_FIXED*/ )
            {
               //nonposind[*nnonposind] = _solver.number(_solver.basis().baseId(i));
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
#ifdef SOLVEIDS_DEBUG
      printf("\n");
#endif
   }



   /// retrieves the compatible columns from the constraint matrix
   // This function also updates the constraint matrix of the reduced problem. It is efficient to perform this in the
   // following function because the required linear algebra has been performed.
   void SoPlex::_getCompatibleColumns(int* nonposind, int* compatind, int* rowsforremoval, int* colsforremoval,
         int nnonposind, int* ncompatind, bool formRedProb)
   {
      Real feastol = realParam(SoPlex::FEASTOL);

      bool compatible;
      SSVector y(_solver.nCols());
      y.unSetup();

      *ncompatind  = 0;

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

         rowtoupdate.setRowVector(newRowVector);

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
      }
      assert(numIncludedRows <= _solver.nRows());
   }



   /// computes the reduced problem objective coefficients
   void SoPlex::_computeReducedProbObjCoeff()
   {
      Real feastol = realParam(SoPlex::EPSILON_ZERO);

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
      _idsLP->changeObj(_transformedObj);
   }



   /// computes the compatible bound constraints and adds them to the reduced problem
   // NOTE: No columns are getting removed from the reduced problem. Only the bound constraints are being removed.
   // So in the reduced problem, all variables are free unless the bound constraints are selected as compatible.
   void SoPlex::_getCompatibleBoundCons(LPRowSet& boundcons, int* compatboundcons, int* nonposind,
         int* ncompatboundcons, int nnonposind)
   {
      Real feastol = realParam(SoPlex::EPSILON_ZERO);

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
   void SoPlex::_deleteAndUpdateRowsComplementaryProblem()
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
      spx_alloc(addedcolid, 1);
      LPColSetReal compSlackCol;
      compSlackCol.add(1.0, -infinity, slackColCoeff, infinity);
      _compSolver.addCols(addedcolid, compSlackCol);
      _compSlackColId = addedcolid[0];

      // freeing allocated memory
      spx_free(addedcolid);
   }



   /// update the dual complementary problem with additional columns and rows
   // Given the solution to the updated reduced problem, the complementary problem will be updated with modifications to
   // the constraints and the removal of variables
   void SoPlex::_updateIdsComplementaryProblem(bool origObj)
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
         if( _realLP->rowType(_idsElimPrimalRowIDs[i]) == LPRowBase<Real>::EQUAL )
            printf("Equality constraint eliminated.\n");

         /*
          * =================================================
          * This is not correct!!! The addition of rows does not take into account the different columns contained in
          * the complementary problem compared to the original problem.
          * MUST FIX!!!!!
          * 07.01.2014 - I think that this is fixed. Still need to check.
          * =================================================
          */
         // Originally this if was surrounded with another statement checking for basis rows. I think that this is
         // unnecessary and the check should be related to the status of the variables.
         int solverRowNum = _solver.number(_idsReducedProbRowIDs[rowNumber]);
         assert(solverRowNum >= 0 && solverRowNum < _solver.nRows());
         if( _solver.basis().desc().rowStatus(solverRowNum) == SPxBasis::Desc::P_ON_UPPER ||
           _solver.basis().desc().rowStatus(solverRowNum) == SPxBasis::Desc::P_ON_LOWER ||
           _solver.basis().desc().rowStatus(solverRowNum) == SPxBasis::Desc::P_FIXED ||
           _solver.basis().desc().rowStatus(solverRowNum) == SPxBasis::Desc::D_FREE ||
           (_solver.basis().desc().rowStatus(solverRowNum) == SPxBasis::Desc::D_ON_LOWER &&
           EQ(_solver.rhs(solverRowNum) - _solver.pVec()[solverRowNum], 0.0, feastol)) ||
           (_solver.basis().desc().rowStatus(solverRowNum) == SPxBasis::Desc::D_ON_UPPER &&
           EQ(_solver.pVec()[solverRowNum] - _solver.lhs(solverRowNum), 0.0, feastol)) )
         {
            // this assert should stay, but there is an issue with the status and the dual vector
            //assert(isNotZero(dualVector[_solver.number(_idsReducedProbRowIDs[rowNumber])]));
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
                  printf("Adding column to complementary problem: %d\n", colNumber);
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
                 EQ(_solver.rhs(solverRowNum) - _solver.pVec()[solverRowNum], 0.0, feastol)) )
            {
               // this assert should stay, but there is an issue with the status and the dual vector
               //assert(GT(dualVector[_solver.number(_idsReducedProbRowIDs[rowNumber])], 0.0));
               // NOTE: This will probably need the SPxRowId passed to the function to get the new row id.
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
              EQ(_solver.pVec()[solverRowNum] - _solver.lhs(solverRowNum), 0.0, feastol)) )
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

#ifdef SOLVEIDS_DEBUG
      printf("Num added cols: %d\n", numElimColsAdded);
#endif

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

#if 1 // 22.06.2015
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
                 EQ(_solver.rhs(solverRowNum) - _solver.pVec()[solverRowNum], 0.0, feastol)) )
            {
               //assert(GT(dualVector[solverRowNum], 0.0));
               _compSolver.changeObj(_idsDualColIDs[i], _realLP->rhs(SPxRowId(_idsPrimalRowIDs[i])));
               _compSolver.changeBounds(_idsDualColIDs[i], -infinity, infinity);
            }
            else if( _solver.basis().desc().rowStatus(solverRowNum) == SPxBasis::Desc::P_ON_LOWER ||
              (_solver.basis().desc().rowStatus(solverRowNum) == SPxBasis::Desc::D_ON_UPPER &&
              EQ(_solver.pVec()[solverRowNum] - _solver.lhs(solverRowNum), 0.0, feastol)) )
            {
               //assert(LT(dualVector[solverRowNum], 0.0));
               _compSolver.changeObj(_idsDualColIDs[i], _realLP->lhs(SPxRowId(_idsPrimalRowIDs[i])));
               _compSolver.changeBounds(_idsDualColIDs[i], -infinity, infinity);
            }
            else //if ( _solver.basis().desc().rowStatus(solverRowNum) != SPxBasis::Desc::D_FREE )
            {
               //assert(isZero(dualVector[solverRowNum], 0.0));

               // 22.06.2015 Testing keeping all rows in the complementary problem
#if 1
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
                     throw SPxInternalCodeException("XLPFRD01 This should never happen.");
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
                  throw SPxInternalCodeException("XLPFRD01 This should never happen.");
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


      int* perm = 0;
      spx_alloc(perm, _compSolver.nCols() + numElimColsAdded);
      _compSolver.removeCols(colsforremoval, ncolsforremoval, perm);

      // updating the dual columns to represent the fixed primal variables.
      int* currFixedVars = 0;
      spx_alloc(currFixedVars, _realLP->nCols());
      _identifyComplementaryFixedPrimalVars(currFixedVars);
      _removeComplementaryFixedPrimalVars(currFixedVars);
#ifdef WRITE_LP
      printf("Writing the complementary problem after column removal\n");
      _compSolver.writeFile("comp_colremove.lp");
#endif
      _updateComplementaryFixedPrimalVars(currFixedVars);


      // freeing allocated memory
      spx_free(currFixedVars);
      spx_free(perm);
      spx_free(colsforremoval);
   }



   /// checking the optimality of the original problem.
   // this function is called if the complementary problem is solved with a non-negative objective value. This implies
   // that the rows currently included in the reduced problem are sufficient to identify the optimal solution to the
   // original problem.
   void SoPlex::_checkOriginalProblemOptimality(Vector primalVector)
   {
      SSVector x(_solver.nCols());
      x.unSetup();

      _idsTransBasis.coSolve(x, primalVector);

      Real redObjVal = 0;
      Real objectiveVal = 0;
      for( int i = 0; i < _solver.nCols(); i++ )
      {
         printf("Variable %d: Red - %f Orig - %f\n", i, primalVector[i], x[i]);
         redObjVal += _solver.maxObj(i)*primalVector[i];
         objectiveVal += _realLP->maxObj(i)*x[i];
      }

      printf("Reduced Problem Objective Value: %f\n", redObjVal);
      printf("Original Problem Objective Value: %f\n", objectiveVal);
   }



   /// updating the slack column coefficients to adjust for equality constraints
   void SoPlex::_updateComplementarySlackColCoeff()
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
   void SoPlex::_identifyComplementaryFixedPrimalVars(int* currFixedVars)
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

#ifdef SOLVEIDS_DEBUG
      printf("Number of fixed variables: %d\n", numFixedVar);
#endif
   }



   /// removing the dual columns related to the fixed variables
   void SoPlex::_removeComplementaryFixedPrimalVars(int* currFixedVars)
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
   void SoPlex::_updateComplementaryFixedPrimalVars(int* currFixedVars)
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


   /// updating the complementary problem with the original objective function
   void SoPlex::_setComplementaryOriginalObjective()
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
               solver.changeObjOffset(0.0);
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
               solver.changeObjOffset(0.0);
               _idsSimplifyAndSolve(solver, sluFactor, false, false);
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
               solver.changeObjOffset(0.0);
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
            throw SPxInternalCodeException("XLPFRD01 This should never happen.");
      }

      return 0;
   }



   /// gets violation of bounds; returns true on success
   bool SoPlex::getIdsBoundViolation(Real& maxviol, Real& sumviol)
   {
      VectorReal& primal = _solReal._primal;
      assert(primal.dim() == _realLP->nCols());

      maxviol = 0.0;
      sumviol = 0.0;

      bool isViol = false;

      for( int i = _realLP->nCols() - 1; i >= 0; i-- )
      {
         Real viol = _realLP->lower(i) - primal[i];

         isViol = false;

         if( viol > 0.0 )
         {
            sumviol += viol;
            if( viol > maxviol )
               maxviol = viol;

            isViol = true;
         }
         //else if( isZero(viol) )
            //printf("Col %d fixed to lower. Redprob fixed: %d\n", i, _fixedOrigVars[i]);

         viol = primal[i] - _realLP->upper(i);
         if( viol > 0.0 )
         {
            sumviol += viol;
            if( viol > maxviol )
               maxviol = viol;

            isViol = true;
         }
         //else if( isZero(viol) )
            //printf("Col %d fixed to upper. Redprob fixed: %d\n", i, _fixedOrigVars[i]);

         if( isViol )
         {
            printf("Violated column: %d. Bounds (%f, %f). Primal: %f\n", i, _realLP->lower(i), _realLP->upper(i),
               primal[i]);
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
      maxviol = 0.0;
      sumviol = 0.0;

      bool isViol = false;

      for( int i = _realLP->nRows() - 1; i >= 0; i-- )
      {
         int solverRowNum = _solver.number(_idsReducedProbRowIDs[i]);

         isViol = false;

         Real viol = _realLP->lhs(i) - activity[i];
         if( viol > 0.0 )
         {
            sumviol += viol;
            if( viol > maxviol )
               maxviol = viol;
         }

         if( GT(viol, 0.0, feastol) )
            isViol = true;

         viol = activity[i] - _realLP->rhs(i);
         if( viol > 0.0 )
         {
            sumviol += viol;
            if( viol > maxviol )
               maxviol = viol;
         }

         //if( GT(viol, 0.0, feastol) )
         if( viol > 0.0 )
            isViol = true;

         if( isViol )
         {
            if( _idsReducedProbRows[i] &&
               (_solver.basis().desc().rowStatus(solverRowNum) == SPxBasis::Desc::P_ON_UPPER ||
               _solver.basis().desc().rowStatus(solverRowNum) == SPxBasis::Desc::P_FIXED ||
               _solver.basis().desc().rowStatus(solverRowNum) == SPxBasis::Desc::P_ON_LOWER ||
               _solver.basis().desc().rowStatus(solverRowNum) == SPxBasis::Desc::D_FREE) )
               printf("I");
            printf("(%d, %d, %d, %d) %.20f ", _solver.basis().desc().rowStatus(solverRowNum), _idsReducedProbRows[i], i,
               solverRowNum, viol);
            printf("Bounds = Lower (%f, %f) Upper (%f, %f)\n", _realLP->lhs(i), _solver.lhs(solverRowNum),
               _realLP->rhs(i), _solver.rhs(solverRowNum));

            if( false && _idsReducedProbRows[i] )
            {
               for( int j = 0; j < _realLP->rowVector(i).size(); j++ )
                  printf("(%d, %f) ", _realLP->rowVector(i).index(j), _realLP->rowVector(i).value(j));
               printf("\n");

               for( int j = 0; j < _solver.rowVector(solverRowNum).size(); j++ )
               {
                  int rowNumber = _solver.number(_idsReducedProbColRowIDs[_solver.rowVector(solverRowNum).index(j)]);
                  printf("(%d, %f, %d) ", _solver.rowVector(solverRowNum).index(j),
                     _solver.rowVector(solverRowNum).value(j), _solver.basis().desc().rowStatus(rowNumber));

                  for( int k = 0; k < _solver.rowVector(rowNumber).size(); k++ )
                     printf("(%d, %f) ", _solver.rowVector(rowNumber).index(k),
                        _solver.rowVector(rowNumber).value(k));

               }
               printf("\n");
            }
         }
      }
      printf("\n");

      return true;
   }



   /// function call to terminate the decomposition simplex
   bool SoPlex::idsTerminate()
   {
      Real maxTime = realParam(SoPlex::TIMELIMIT);

      // check if a time limit is actually set
      if( maxTime < 0 || maxTime >= infinity )
         return false;

      Real currentTime = _statistics->solvingTime->time();
      if ( currentTime >= maxTime )
      {
         MSG_INFO2( spxout, spxout << " --- timelimit (" << _solver.getMaxTime()
            << ") reached" << std::endl; )
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
                  throw SPxInternalCodeException("XLPFRD01 This should never happen.");
            }
         }
      }

      nNonBasicRows = _realLP->nRows() - basicRow - nDegenerateRows;
      printf("Number of non-basic rows: %d (from %d)\n", nNonBasicRows, _realLP->nRows());
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
      printf("Number of non-basic columns: %d (from %d)\n", nNonBasicCols, _realLP->nCols());
      printf("Number of zero dual columns: %d (from %d)\n", numZeroDual, _realLP->nCols());
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

      for( int i = 0; i < _solver.nCols(); i++ )
      {
         int solverNum = _solver.number(_solver.basis().baseId(i));
         if( _solver.basis().baseId(i).isSPxRowId() )
            printf("%d: Row %d, Status: %d\n", i, solverNum, _solver.basis().desc().rowStatus(solverNum));
         else
            printf("%d: Col %d, Status: %d\n", i, solverNum, _solver.basis().desc().colStatus(solverNum));

      }

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
