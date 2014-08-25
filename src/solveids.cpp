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

/* This file contains the private functions for the Improved Dual Simplex (IDS)
 *
 * An important note about the implementation of the IDS is the reliance on the row representation of the basis matrix.
 * The forming of the reduced and complementary problems is involves identifying rows in the row-form of the basis
 * matrix that have zero dual multipliers.
 *
 * Ideally, this work will be extended such that the IDS can be executed using the column-form of the basis matrix. */

//#define SOLVEIDS_DEBUG

namespace soplex
{
   /// solves LP using the improved dual simplex
   void SoPlex::_solveImprovedDualSimplex()
   {
      assert(_solver.rep() == SPxSolver::ROW);
      assert(_isRealLPLoaded == true);

      DVector reducedLPDualVector(_solver.nRows());
      DVector compLPPrimalVector(_solver.nCols());

      // start timing
      _statistics->solvingTime.start();

      // it is necessary to solve the initial problem to find a starting basis
      _solver.setIdsStatus(SPxSolver::FINDSTARTBASIS);
      //while (!happy)
      //{
         _solver.solve();
         //DVector x(_solver.nCols());
         //_solver.basis.solve(x, _solver.maxObj());
         //if()
            //happy = true;
      //}
      _solver.setIdsStatus(SPxSolver::DONTFINDSTARTBASIS);
      _hasBasis = true; // this is probably wrong. Will need to set this properly.

#ifdef SOLVEIDS_DEBUG
      printf("Writing the original lp to a file\n");
      _solver.writeState("original");
#endif

      //@todo Need to check this. Where will the solving of the IDS take place? In the SoPlex::solve, _solvereal or
      //_solverational?
      // remember that last solve was in floating-point
      _lastSolveMode = SOLVEMODE_REAL;

      // creating copies of the original problem that will be manipulated to form the reduced and complementary
      // problems.
      _createIdsReducedAndComplementaryProblems();

      // creating the initial reduced problem from the basis information
      _formIdsReducedProblem();
      _isRealLPLoaded = false;

#ifdef SOLVEIDS_DEBUG
      printf("Writing the reduced lp to a file\n");
      _solver.writeState("reduced");
#endif
      for( int i = 0; i < 2; i++ )
      {
         // solve the reduced problem
         _solver.solve();
         // get the dual solutions from the reduced problem
         reducedLPDualVector.reDim(_solver.nRows());
         _solver.getDual(reducedLPDualVector);

         // create the complementary problem using the solution to the reduced problem
         _formIdsComplementaryProblem();

#ifdef SOLVEIDS_DEBUG
         printf("Writing the complementary lp to a file\n");
         _compSolver.writeState("complement");
#endif
         // solve the complementary problem
         _compSolver.solve();
         // get the primal solutions from the complementary problem
         compLPPrimalVector.reDim(_compSolver.nCols());
         _compSolver.getPrimal(compLPPrimalVector);

//#ifdef SOLVEIDS_DEBUG
      //printf("Writing redLP to a file\n");
      //_solver.writeState("redLP");
//#endif

         //_solver.solve(); // Solving Original Reduced (redLP) problem
         //_solver.loadLP(*_idsRedLP);
//#ifdef SOLVEIDS_DEBUG
      //printf("Writing _idsRedLP to a file\n");
      //_solver.writeState("_idsRedLP");
//#endif
         //_solver.solve(); // Solving reduced (*_idsRedLP)  problem

         _updateIdsReducedProblem(reducedLPDualVector, compLPPrimalVector);
#ifdef SOLVEIDS_DEBUG
      printf("Writing updated _idsRedLP to a file\n");
      _solver.writeState("_idsRedLP_U");
#endif
         _solver.solve(); // Solving the updated reduced problem (*_idsRedLP)
      }

      // stop timing
      _statistics->solvingTime.stop();
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

      // allocating memory for the reduced problem rows flag array
      _idsReducedProbRows = 0;
      spx_alloc(_idsReducedProbRows, numRowsReal());

      // the complementary problem is formulated with all incompatible rows and those from the reduced problem that have
      // a positive reduced cost.
      _compSolver = _solver;
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

      DVector dualSolutions(numRowsReal());
      DVector reducedCosts(numColsReal());
      DVector primalVector(numColsReal());
      _solver.getDual(dualSolutions);
      _solver.getRedCost(reducedCosts);
      _solver.getPrimal(primalVector);
      for( int i = 0; i < numRowsReal(); i++ )
         printf("%f ", dualSolutions[i]);
      printf("\n");
      for( int i = 0; i < numColsReal(); i++ )
         printf("%f ", reducedCosts[i]);
      printf("\n");
      for( int i = 0; i < numColsReal(); i++ )
         printf("%f ", primalVector[i]);
      printf("\n");

      // retreiving the basis information
      _basisStatusRows.reSize(numRowsReal());
      _basisStatusCols.reSize(numColsReal());
      _solver.getBasis(_basisStatusRows.get_ptr(), _basisStatusCols.get_ptr());

      // get thie indices of the rows with positive dual multipliers and columns with positive reduced costs.
      spx_alloc(bind, numColsReal());

      getBasisInd(bind);
      printf("bind: ");
      for( int i = 0; i < numRowsReal(); i++ )
         printf("%d ", bind[i]);
      printf("\n");


      spx_alloc(nonposind, numColsReal());
      spx_alloc(colsforremoval, numColsReal());
      _getNonPositiveDualMultiplierInds(nonposind, bind, colsforremoval, &nnonposind);
      for( int i = 0; i < numColsReal(); i++ )
         printf("%d ", bind[i]);
      printf("\n");

      for( int i = 0; i < numColsReal(); i++ )
         printf("%d %d\n", i, colsforremoval[i]);

      // get the compatible columns from the constraint matrix w.r.t the current basis matrix
#ifdef SOLVEIDS_DEBUG
      printf("Computing the compatible columns\n");
      printf("Solving time: %f\n", solveTime());
#endif
      spx_alloc(compatind, numRowsReal());
      spx_alloc(rowsforremoval, numRowsReal());
      _getCompatibleColumns(nonposind, compatind, rowsforremoval, nnonposind, &ncompatind);

      // delete rows and columns from the LP to form the reduced problem
#ifdef SOLVEIDS_DEBUG
      printf("Deleting rows and columns to form the reduced problem\n");
      printf("Solving time: %f\n", solveTime());
#endif
      for( int i = 0; i < numRowsReal(); i++ )
         printf("%d %d\n", i, rowsforremoval[i]);

      // the colsforremoval are the columns with a zero reduced cost.
      // the rowsforremoval are the rows identified as incompatible.
      _deleteRowsAndColumnsReducedProblem(colsforremoval, rowsforremoval);

      for( int i = 0; i < numRowsReal(); i++ )
         printf("%d %d\n", i, rowsforremoval[i]);

      printf("%f\n", objValueReal());

      // freeing allocated memory
      spx_free(rowsforremoval);
      spx_free(compatind);
      spx_free(colsforremoval);
      spx_free(nonposind);
      spx_free(bind);
   }

   /// forms the complementary problem
   void SoPlex::_formIdsComplementaryProblem()
   {
      int* bind = 0;
      int* nonposind = 0;
      int* colsforremoval = 0;
      int* rowsforremoval = 0;
      int nnonposind = 0;
      int nrowsforremoval = 0;

      // get the indices of the rows with positive dual multipliers and columns with positive reduced costs.
      spx_alloc(bind, _compSolver.nCols());
      spx_alloc(nonposind, _compSolver.nCols());
      spx_alloc(colsforremoval, _compSolver.nCols());
      _getNonPositiveDualMultiplierInds(nonposind, bind, colsforremoval, &nnonposind);

      for( int i = 0; i < _compSolver.nCols(); i++ )
         printf("%d %d\n", i, colsforremoval[i]);
      printf("\n");


      // get the rows for removal from the nonposind and the bind
      spx_alloc(rowsforremoval, _compSolver.nRows());
      _getRowsForRemovalComplementaryProblem(nonposind, bind, rowsforremoval, &nrowsforremoval, nnonposind);

      for( int i = 0; i < nnonposind; i++ )
         printf("%d %d\n", i, nonposind[i]);
      printf("\n");

      for( int i = 0; i < nrowsforremoval; i++ )
         printf("%d %d\n", i, rowsforremoval[i]);
      printf("\n");

      // delete rows and columns from the LP to form the reduced problem
      _deleteAndUpdateRowsComplementaryProblem(rowsforremoval, nrowsforremoval);

      // initialising the arrays to store the row id's from the primal and the col id's from the dual
      _idsPrimalRowIDs.reSize(2*_compSolver.nRows());
      _idsDualColIDs.reSize(2*_compSolver.nRows());
      _nPrimalRows = 0;
      _nDualCols = 0;

      // convert complementary problem to dual problem
      SPxLPIds compDualLP;
      _compSolver.buildDualProblem(&compDualLP, _idsPrimalRowIDs.get_ptr(), _idsDualColIDs.get_ptr(), &_nPrimalRows,
            &_nDualCols);

      _idsPrimalRowIDs.reSize(_nPrimalRows);
      _idsDualColIDs.reSize(_nDualCols);

      _compSolver.loadLP(compDualLP);

      // freeing allocated memory
      spx_free(colsforremoval);
      spx_free(nonposind);
      spx_free(bind);
   }

   /// updates the reduced problem with additional rows using the solution to the complementary problem
   void SoPlex::_updateIdsReducedProblem(DVector dualVector, DVector compPrimalVector)
   {
      LPRowSet updaterows;

      Real maxDualRatio = infinity;

      for( int i = 0; i < _nPrimalRows; i++ )
      {
         Real reducedProbDual = 0;
         Real compProbPrimal = 0;
         Real dualRatio = 0;
         if( _idsReducedProbRows[_realLP->number(SPxRowId(_idsPrimalRowIDs[i]))] )
         {
            // retreiving the reduced problem dual solutions and the complementary problem primal solutions
            reducedProbDual = dualVector[_solver.number(SPxRowId(_idsPrimalRowIDs[i]))]; // this is y
            compProbPrimal = compPrimalVector[_compSolver.number(SPxColId(_idsDualColIDs[i]))]; // this is u

            // the translation of the complementary primal problem to the dual some rows resulted in two columns.
            if( i < _nPrimalRows - 1 &&
                  _realLP->number(SPxRowId(_idsPrimalRowIDs[i])) == _realLP->number(SPxRowId(_idsPrimalRowIDs[i + 1])) )
            {
               i++;
               // @todo make sure that this is just a simple sum
               compProbPrimal += compPrimalVector[_compSolver.number(SPxColId(_idsDualColIDs[i]))]; // this is u
            }

            // updating the ratio
            if( compProbPrimal <= 0 )
               dualRatio = infinity;
            else
               dualRatio = reducedProbDual/compProbPrimal;

            if( dualRatio < maxDualRatio )
               maxDualRatio = dualRatio;
         }
      }

      for( int i = 0; i < _nPrimalRows; i++ )
      {
         LPRowReal rowtoadd;
         Real compProbPrimal = 0;
         if( !_idsReducedProbRows[_realLP->number(SPxRowId(_idsPrimalRowIDs[i]))] )
         {
            // retreiving the complementary problem primal solutions
            compProbPrimal = compPrimalVector[_compSolver.number(SPxColId(_idsDualColIDs[i]))]; // this is u

            // the translation of the complementary primal problem to the dual some rows resulted in two columns.
            if( i < _nPrimalRows - 1 &&
                  _realLP->number(SPxRowId(_idsPrimalRowIDs[i])) == _realLP->number(SPxRowId(_idsPrimalRowIDs[i + 1])) )
            {
               i++;
               // @todo make sure that this is just a simple sum
               compProbPrimal += compPrimalVector[_compSolver.number(SPxColId(_idsDualColIDs[i]))]; // this is u
            }

            // add row to the reduced problem the computed dual is positive
            if( compProbPrimal*maxDualRatio > 0 )
            {
               _realLP->getRow(_realLP->number(SPxRowId(_idsPrimalRowIDs[i])), rowtoadd);
               updaterows.add(rowtoadd);
               _idsReducedProbRows[_realLP->number(SPxRowId(_idsPrimalRowIDs[i]))] = true;
            }
         }
      }

      _solver.addRows(updaterows);
   }


   // This function assumes that the basis is in the row form.
   // @todo extend this to the case when the basis is in the column form.
   //
   // NOTE: Changing "nonposind[*nnonposind] = bind[i]" to "nonposind[*nnonposind] = i"
   void SoPlex::_getNonPositiveDualMultiplierInds(int* nonposind, int* bind, int* colsforremoval, int* nnonposind) const
   {
      assert(_solver.rep() == SPxSolver::ROW);
      bool delCol;

      *nnonposind = 0;

      // iterating over all columns in the basis matrix
      // this identifies the basis indices and the indicies that are positive.
      for( int i = 0; i < _solver.nCols(); ++i ) // @todo Check the use of numColsReal for the reduced problem.
      {
         printf("%.20f ", _solver.fVec()[i]);
         delCol = false;
         // @todo I have questions about my implementation of this function. I don't think that I am getting the right
         // information. Additionally, in getCompatibleColumns the information may not be used correctly.
         if( _solver.basis().baseId(i).isSPxRowId() ) // find the row id's for rows in the basis
         {
            bind[i] = -1 - _realLP->number(SPxRowId(_solver.basis().baseId(i))); // getting the corresponding row
                                                                                   // for the original LP.

            //@todo need to check this regarding min and max problems
            if( _solver.fVec()[i] > -1*_solver.opttol()
                  && _solver.fVec()[i] < _solver.opttol() )
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
            bind[i] = _realLP->number(SPxColId(_solver.basis().baseId(i)));

            if ( _solver.spxSense() == SPxLP::MINIMIZE )
            {
               //@todo need to check this regarding min and max problems
               if( _solver.fVec()[i] > -1*_solver.opttol()
                     && _solver.fVec()[i] < _solver.opttol() )
               {
                  nonposind[*nnonposind] = i;
                  (*nnonposind)++;

                  delCol = true;
               }
            }
            else
            {
               //@todo need to check this regarding min and max problems
               if( _solver.fVec()[i] > -1*_solver.opttol()
                     && _solver.fVec()[i] < _solver.opttol() )
               {
                  nonposind[*nnonposind] = i;
                  (*nnonposind)++;

                  delCol = true;
               }
            }
         }

         // setting an array to identify the columns to be removed from the LP to form the reduced problem
         if( delCol )
            colsforremoval[i] = -1;
         else
            colsforremoval[i] = i;
      }
      printf("\n");
   }


   /// retrieves the compatible columns from the constraint matrix
   void SoPlex::_getCompatibleColumns(int* nonposind, int* compatind, int* rowsforremoval, int nnonposind,
         int* ncompatind)
   {
      bool compatible;
      SSVector y(numColsReal());

      *ncompatind  = 0;

      for( int i = 0; i < numRowsReal(); ++i )
      {
         rowsforremoval[i] = i;
         _idsReducedProbRows[i] = true;

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
            if( y[nonposind[j]] != 0 )
            {
               compatible = false;
               break;
            }
         }

         if( compatible )
         {
            compatind[*ncompatind] = i;
            (*ncompatind)++;
         }
         else
         {
            // setting an array to identify the rows to be removed from the LP to form the reduced problem
            rowsforremoval[i] = -1;
            _idsReducedProbRows[i] = false;
         }
      }
   }

   // @todo need to put in a check similar to _isRealLPLoaded for the ids LP.
   void SoPlex::_deleteRowsAndColumnsReducedProblem(int* colsforremoval, int* rowsforremoval)
   {
      assert(_basisStatusRows.size() == numRowsReal());
      _solver.removeRows(rowsforremoval);
      _solver.removeCols(colsforremoval);


      // removing rows from the solver LP
      // @todo not sure whether the first if statment in this function is valid.
      if( _isRealLPLoaded )
      {
         _hasBasis = (_solver.basis().status() > SPxBasis::NO_PROBLEM);
      }
      else if( _hasBasis )
      {
         for( int i = numRowsReal() - 1; i >= 0 && _hasBasis; i-- )
         {
            printf("%d %d %d\n", rowsforremoval[i], rowsforremoval[rowsforremoval[i]], _basisStatusRows[i]);
            if( rowsforremoval[i] < 0 && _basisStatusRows[i] != SPxSolver::BASIC )
               _hasBasis = false;
            else if( rowsforremoval[i] >= 0 && rowsforremoval[i] != i )
            {
               assert(rowsforremoval[i] < numRowsReal());
               assert(rowsforremoval[rowsforremoval[i]] < 0);

               _basisStatusRows[rowsforremoval[i]] = _basisStatusRows[i];
            }
         }

         if( _hasBasis )
            _basisStatusRows.reSize(numRowsReal());

         // removing columns from the solver LP
         if( _hasBasis )
         {
            for( int i = numColsReal() - 1; i >= 0 && _hasBasis; i-- )
            {
               if( colsforremoval[i] < 0 && _basisStatusCols[i] == SPxSolver::BASIC )
                  _hasBasis = false;
               else if( colsforremoval[i] >= 0 && colsforremoval[i] != i )
               {
                  assert(colsforremoval[i] < numColsReal());
                  assert(colsforremoval[colsforremoval[i]] < 0);

                  _basisStatusCols[colsforremoval[i]] = _basisStatusCols[i];
               }
            }

            if( _hasBasis )
               _basisStatusCols.reSize(numColsReal());
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
   void SoPlex::_deleteAndUpdateRowsComplementaryProblem(int* rowsforremoval, int nrowsforremoval)
   {
      int* perm = 0;
      spx_alloc(perm, numRowsReal());

      _compSolver.removeRows(rowsforremoval, nrowsforremoval, perm);

      for( int i = 0; i < numRowsReal(); i++ )
      {
         if( _idsReducedProbRows[i] && perm[i] >= 0 )
         {
            // it is assumed that a ranged constraint will be at it lower bound for a minimisation problem and its upper
            // bound for a maximisation problem.
            // @todo this is incorrect. The setting of the bounds will depend on the type of constraint. It is possible
            // that a <= constraint can be at its bound in a minimisation problem, hence it will require the equality to
            // be set to the rhs.
            // @todo ask Ambros about the constraint type. Is it set by the user or Soplex, and does it get reset when
            // constraints are changed.
            if (_solver.spxSense() == SPxLP::MINIMIZE)
               _changeRhsReal(i, lhsReal(i));
            else
               _changeLhsReal(i, rhsReal(i));
         }
      }
   }

   // @todo update this function and related comments. It has only been hacked together.
   /// checks result of the solving process and solves again without preprocessing if necessary
   // @todo need to evaluate the solution to ensure that it is solved to optimality and then we are able to perform the
   // next steps in the algorithm.
   void SoPlex::_evaluateSolutionIDS()
   {
      // process result
      switch( _status )
      {
      case SPxSolver::OPTIMAL:
         if( !_isRealLPLoaded )
         {
            MSG_INFO1( spxout << " --- updating the basis partitioning" << std::endl; )
            // Need to solve the complementary problem
            return;
         }
         else
            _hasBasis = true;
         break;

      case SPxSolver::UNBOUNDED:
      case SPxSolver::INFEASIBLE:
      case SPxSolver::INForUNBD:
      case SPxSolver::SINGULAR:
      case SPxSolver::ABORT_CYCLING:
         // when solving the reduced problem it should not be possible that the problem is infeasible or unbounded. This
         // result is a check to make sure that the partitioning is correct.
         if( !_isRealLPLoaded )
         {
            // in this situation the original problem should be solved again.
            _solver.changeObjOffset(0.0); // NOTE: What is this????
            // build in some recourse to quit the IDS if the reduced problem is infeasible or unbounded
            return;
         }
         else
            _hasBasis = (_solver.basis().status() > SPxBasis::NO_PROBLEM);
         break;

         // else fallthrough
      case SPxSolver::ABORT_TIME:
      case SPxSolver::ABORT_ITER:
      case SPxSolver::ABORT_VALUE:
      case SPxSolver::REGULAR:
      case SPxSolver::RUNNING:
         if( _isRealLPLoaded )
            _hasBasis = (_solver.basis().status() > SPxBasis::NO_PROBLEM);
         // store regular basis if there is no simplifier and the original problem is not in the solver because of
         // scaling; non-optimal bases should currently not be unsimplified
         else if( _simplifier == 0 && _solver.basis().status() > SPxBasis::NO_PROBLEM )
         {
            _basisStatusRows.reSize(numRowsReal());
            _basisStatusCols.reSize(numColsReal());
            assert(_basisStatusRows.size() == _solver.nRows());
            assert(_basisStatusCols.size() == _solver.nCols());

            _solver.getBasis(_basisStatusRows.get_ptr(), _basisStatusCols.get_ptr());
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

} // namespace soplex
#endif
