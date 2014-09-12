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

#define SOLVEIDS_DEBUG
//#define EXTRA_PRINT

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

      // start timing
      _statistics->solvingTime.start();

      // setting the sense to maximise
      if( intParam(SoPlex::OBJSENSE) == SoPlex::OBJSENSE_MINIMIZE )
      {
         assert(_solver.spxSense() == SPxLPBase<Real>::MINIMIZE);

         _solver.changeObj(_solver.maxObj());
         _solver.changeSense(SPxLPBase<Real>::MAXIMIZE);
      }

#ifdef SOLVEIDS_DEBUG
      printf("Writing the original lp to a file\n");
      _solver.writeFile("original.lp");
#endif

      // it is necessary to solve the initial problem to find a starting basis
      _solver.setIdsStatus(SPxSolver::FINDSTARTBASIS);

      SPxSolver::Type startingType = _solver.type();
      Real degeneracyLevel = 0;
      _idsFeasVector.reDim(_solver.nCols());
      // since the original LP may have been shifted, the dual multiplier will not be correct for the original LP. This
      // loop will recheck the degeneracy level and compute the proper dual multipliers.
      do
      {
         _solver.solve();
         _solver.basis().solve(_idsFeasVector, _solver.maxObj());
         degeneracyLevel = _solver.getDegeneracyLevel(_idsFeasVector);
         if( _solver.type() != startingType )
         {
            // returning the sense to minimise
            if( intParam(SoPlex::OBJSENSE) == SoPlex::OBJSENSE_MINIMIZE )
            {
               assert(_solver.spxSense() == SPxLPBase<Real>::MAXIMIZE);

               _solver.changeObj(-(_solver.maxObj()));
               _solver.changeSense(SPxLPBase<Real>::MINIMIZE);
            }

            _solveReal();
            return;
         }
         printf("Degeneracy Level: %f\n", degeneracyLevel);
      } while( degeneracyLevel > 0.8 || degeneracyLevel < 0.1 );

      _solver.setIdsStatus(SPxSolver::DONTFINDSTARTBASIS);
      _hasBasis = true; // this is probably wrong. Will need to set this properly.


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

      for( int i = 0; i < 5; i++ )
      {
#ifdef SOLVEIDS_DEBUG
      printf("Writing the reduced lp to a file\n");
      _solver.writeFile("reduced.lp");
#endif
         printf("Reduced Prob Rows: ");
         for( int j = 0; j < numRowsReal(); j++ )
         {
            if( _idsReducedProbRows[j] )
               printf("%d ", j);
         }
         printf("\n");

         // solve the reduced problem
         _solver.solve();

         _idsFeasVector.reDim(_solver.nCols());
         _solver.basis().solve(_idsFeasVector, _solver.maxObj());

         // get the dual solutions from the reduced problem
         reducedLPDualVector.reDim(_solver.nRows());
         _solver.getDual(reducedLPDualVector);

         for( int j = 0; j < _solver.nRows(); j++ )
         {
            if( !isZero(reducedLPDualVector[j]) )
               printf("%d %f\n", j, reducedLPDualVector[j]);
         }

         if( i == 0 )
            _formIdsComplementaryProblem();  // create the complementary problem using
                                             // the solution to the reduced problem
         else
            _updateIdsComplementaryProblem(reducedLPDualVector);
#ifdef SOLVEIDS_DEBUG
         printf("Writing the complementary lp to a file\n");
         _compSolver.writeFile("complement.lp");
#endif
         // solve the complementary problem
         _compSolver.solve();
#ifdef SOLVEIDS_DEBUG
         printf("Objective Value: %f\n", _compSolver.objValue());
#endif
         // check the optimality of the original problem with the objective value of the complementary problem
         //if( _compSolver.objValue() >= 0 )


         // get the primal solutions from the complementary problem
         DVector compLPPrimalVector(_compSolver.nCols());
         _compSolver.getPrimal(compLPPrimalVector);

         _updateIdsReducedProblem(reducedLPDualVector, compLPPrimalVector);
      }

      // returning the sense to minimise
      if( intParam(SoPlex::OBJSENSE) == SoPlex::OBJSENSE_MINIMIZE )
      {
         assert(_solver.spxSense() == SPxLPBase<Real>::MAXIMIZE);

         _solver.changeObj(-(_solver.maxObj()));
         _solver.changeSense(SPxLPBase<Real>::MINIMIZE);

         // Need to add commands to multiple the objective solution values by -1
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


#ifdef EXTRA_PRINT
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
         printf("%d %f\n", i, reducedCosts[i]);
      printf("\n");
      for( int i = 0; i < numColsReal(); i++ )
         printf("%f ", primalVector[i]);
      printf("\n");
#endif

      // retreiving the basis information
      _basisStatusRows.reSize(numRowsReal());
      _basisStatusCols.reSize(numColsReal());
      _solver.getBasis(_basisStatusRows.get_ptr(), _basisStatusCols.get_ptr());

      // get thie indices of the rows with positive dual multipliers and columns with positive reduced costs.
      spx_alloc(bind, numColsReal());

#ifdef EXTRA_PRINT
      getBasisInd(bind);
      printf("bind: ");
      for( int i = 0; i < numRowsReal(); i++ )
         printf("%d ", bind[i]);
      printf("\n");
#endif


      spx_alloc(nonposind, numColsReal());
      spx_alloc(colsforremoval, numColsReal());
      _getNonPositiveDualMultiplierInds(_idsFeasVector, nonposind, bind, colsforremoval, &nnonposind);

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

#ifdef EXTRA_PRINT
      for( int i = 0; i < numRowsReal(); i++ )
         printf("%d %d\n", i, rowsforremoval[i]);
#endif

      // the colsforremoval are the columns with a zero reduced cost.
      // the rowsforremoval are the rows identified as incompatible.
      _deleteRowsAndColumnsReducedProblem(colsforremoval, rowsforremoval);

      for( int i = 0; i < numRowsReal(); i++ )
         printf("%d %d\n", i, rowsforremoval[i]);

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
      // delete rows and columns from the LP to form the reduced problem
      _nElimPrimalRows = 0;
      _idsElimPrimalRowIDs.reSize(_solver.nRows());   // the number of eliminated rows is less than the number of rows
                                                      // in the reduced problem
      _deleteAndUpdateRowsComplementaryProblem();

      // initialising the arrays to store the row id's from the primal and the col id's from the dual
      _idsPrimalRowIDs.reSize(2*_compSolver.nRows());
      _idsDualColIDs.reSize(2*_compSolver.nRows());
      _nPrimalRows = 0;
      _nDualCols = 0;

      // convert complementary problem to dual problem
      printf("Writing the complementary primal lp to a file\n");
      _compSolver.writeFile("comp_primal.lp");

      SPxLPIds compDualLP;
      _compSolver.buildDualProblem(compDualLP, _idsPrimalRowIDs.get_ptr(), _idsDualColIDs.get_ptr(), &_nPrimalRows,
            &_nDualCols);

      // retrieving the dual row id for the complementary slack column
      // there should be a one to one relationship between the number of primal columns and the number of dual rows.
      // hence, it should be possible to equate the dual row id to the related primal column.
      assert(_compSolver.nCols() == compDualLP.nRows());
      _compSlackDualRowId = compDualLP.rId(_compSolver.number(_compSlackColId));

      _compSolver.loadLP(compDualLP);

      _updateComplementarySlackColCoeff();

      printf("Writing the complementary primal lp to a file\n");
      _compSolver.writeFile("comp_dual.lp");

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

            printf("%d %d %f %f\n", i, _realLP->number(SPxRowId(_idsPrimalRowIDs[i])), reducedProbDual, compProbPrimal);

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

      int* addedrowidx = 0;
      int naddedrowidx = 0;
      spx_alloc(addedrowidx, _nPrimalRows);
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
               //updaterows.add(_idsReducedProbRowIDs[_realLP->number(SPxRowId(_idsPrimalRowIDs[i]))], rowtoadd);
               updaterows.add(rowtoadd);

               _idsReducedProbRows[_realLP->number(SPxRowId(_idsPrimalRowIDs[i]))] = true;
               addedrowidx[naddedrowidx] = _realLP->number(SPxRowId(_idsPrimalRowIDs[i]));
               naddedrowidx++;

               SVector addedRow = _realLP->rowVector(_idsPrimalRowIDs[i]);
               for( int j = 0; j < addedRow.size(); j++ )
                  printf("Row Vector: %d %f\n", addedRow.index(j), addedRow.value(j));
            }
         }
      }

      SPxRowId* addedids = 0;
      spx_alloc(addedids, naddedrowidx);
      _solver.addRows(addedids, updaterows);

      for( int i = 0; i < naddedrowidx; i++ )
         _idsReducedProbRowIDs[addedrowidx[i]] = addedids[i];

      spx_free(addedids);
      spx_free(addedrowidx);
   }



   // This function assumes that the basis is in the row form.
   // @todo extend this to the case when the basis is in the column form.
   //
   // NOTE: Changing "nonposind[*nnonposind] = bind[i]" to "nonposind[*nnonposind] = i"
   void SoPlex::_getNonPositiveDualMultiplierInds(Vector feasVector, int* nonposind, int* bind, int* colsforremoval,
         int* nnonposind) const
   {
      assert(_solver.rep() == SPxSolver::ROW);
      bool delCol;

      *nnonposind = 0;

      // iterating over all columns in the basis matrix
      // this identifies the basis indices and the indicies that are positive.
      for( int i = 0; i < _solver.nCols(); ++i ) // @todo Check the use of numColsReal for the reduced problem.
      {
         printf("%d %f ", i, feasVector[i]);
         delCol = false;
         // @todo I have questions about my implementation of this function. I don't think that I am getting the right
         // information. Additionally, in getCompatibleColumns the information may not be used correctly.
         if( _solver.basis().baseId(i).isSPxRowId() ) // find the row id's for rows in the basis
         {
            printf("Row\n");
            bind[i] = -1 - _realLP->number(SPxRowId(_solver.basis().baseId(i))); // getting the corresponding row
                                                                                   // for the original LP.

            //@todo need to check this regarding min and max problems
            if( isZero(feasVector[i]) )
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
            printf("Col\n");
            bind[i] = _realLP->number(SPxColId(_solver.basis().baseId(i)));

            if( isZero(feasVector[i]) )
            {
               nonposind[*nnonposind] = i;
               (*nnonposind)++;

               delCol = true;
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

      _idsReducedProbRowIDs.reSize(numRowsReal());

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

            _idsReducedProbRowIDs[i] = _solver.rowId(i);
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
      // @todo not sure whether the first if statement in this function is valid.
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
   void SoPlex::_deleteAndUpdateRowsComplementaryProblem()
   {
      int nrowsforremoval = 0;
      int* rowsforremoval = 0;
      DSVector slackColCoeff;
      Real slackCoeff = 1.0;

      spx_alloc(rowsforremoval, _solver.nRows());
      for( int i = 0; i < _solver.nCols(); ++i ) // @todo Check the use of numColsReal for the reduced problem.
      {
         // it is assumed that a ranged constraint will be at it lower bound for a minimisation problem and its upper
         // bound for a maximisation problem.
         // @todo this is incorrect. The setting of the bounds will depend on the type of constraint. It is possible
         // that a <= constraint can be at its bound in a minimisation problem, hence it will require the equality to
         // be set to the rhs.
         // @todo ask Ambros about the constraint type. Is it set by the user or Soplex, and does it get reset when
         // constraints are changed.
         if( _solver.basis().baseId(i).isSPxRowId() ) // find the row id's for rows in the basis
         {
            printf("%d %f\n", _solver.number(SPxRowId(_solver.basis().baseId(i))),
                  _solver.fVec()[i]);
            if( isZero(_solver.fVec()[i]) )
            {
               _idsElimPrimalRowIDs[_nElimPrimalRows] = SPxRowId(_solver.basis().baseId(i));
               _nElimPrimalRows++;
               rowsforremoval[nrowsforremoval] = _realLP->number(SPxRowId(_solver.basis().baseId(i)));
               nrowsforremoval++;
            }
            else if( GT(_solver.fVec()[i], 0.0) )
               _compSolver.changeLhs(SPxRowId(_solver.basis().baseId(i)),
                     _solver.rhs(SPxRowId(_solver.basis().baseId(i))));
            else
               _compSolver.changeRhs(SPxRowId(_solver.basis().baseId(i)),
                     _solver.lhs(SPxRowId(_solver.basis().baseId(i))));
         }
      }

      for( int i = 0; i < numRowsReal(); i++ )
      {
         if( !_idsReducedProbRows[i] )
         {
            switch( _realLP->rowType(i) )
            {
               case LPRowBase<Real>::RANGE:
               case LPRowBase<Real>::GREATER_EQUAL:
                  slackColCoeff.add(i, -slackCoeff);
                  break;
               case LPRowBase<Real>::EQUAL:
               case LPRowBase<Real>::LESS_EQUAL:
                  slackColCoeff.add(i, slackCoeff);
                  break;
               default:
                  throw SPxInternalCodeException("XLPFRD01 This should never happen.");
            }
         }
         else if( !_solver.isBasic(_realLP->rId(i)) )
         {
            _idsElimPrimalRowIDs[_nElimPrimalRows] = _realLP->rId(i);
            _nElimPrimalRows++;
            rowsforremoval[nrowsforremoval] = i;
            nrowsforremoval++;
         }
      }

      // setting the objective coefficients of the original variables to zero
      LPRowSetBase<Real> boundCons;
      DSVectorBase<Real> col(1);
      DVector newObjCoeff(numColsReal());
      int addedRowsCounter = 0;
      for( int i = 0; i < numColsReal(); i++ )
      {
         col.add(i, 1.0);
         if( GT(_realLP->lower(i), -infinity) && !isZero(_realLP->lower(i)) )
         {
            if( isZero(_realLP->upper(i)) )
            {
               boundCons.add(_realLP->lower(i), col, Real(infinity));
               slackColCoeff.add(numRowsReal() + addedRowsCounter, -slackCoeff);
               addedRowsCounter++;

               _compSolver.changeBounds(_realLP->cId(i), Real(-infinity), 0.0);
            }
            else
            {
               boundCons.add(_realLP->lower(i), col, Real(infinity));
               slackColCoeff.add(numRowsReal() + addedRowsCounter, -slackCoeff);
               addedRowsCounter++;

               boundCons.add(Real(-infinity), col, _realLP->upper(i));
               slackColCoeff.add(numRowsReal() + addedRowsCounter, slackCoeff);
               addedRowsCounter++;

               _compSolver.changeBounds(_realLP->cId(i), Real(-infinity), Real(infinity));
            }
         }
         else if( LT(_realLP->upper(i), infinity) && !isZero(_realLP->upper(i)) )
         {
            if( isZero(_realLP->lower(i)) )
            {
               boundCons.add(Real(-infinity), col, _realLP->upper(i));
               slackColCoeff.add(numRowsReal() + addedRowsCounter, slackCoeff);
               addedRowsCounter++;

               _compSolver.changeBounds(_realLP->cId(i), 0.0, Real(infinity));
            }
            else
            {
               boundCons.add(_realLP->lower(i), col, Real(infinity));
               slackColCoeff.add(numRowsReal() + addedRowsCounter, -slackCoeff);
               addedRowsCounter++;

               boundCons.add(Real(-infinity), col, _realLP->upper(i));
               slackColCoeff.add(numRowsReal() + addedRowsCounter, slackCoeff);
               addedRowsCounter++;

               _compSolver.changeBounds(_realLP->cId(i), Real(-infinity), Real(infinity));
            }
         }
         col.clear();
         newObjCoeff[i] = 0;
      }

      _compSolver.addRows(boundCons);
      _compSolver.changeObj(newObjCoeff);

      //adding the slack column to the complementary problem
      SPxColId* addedcolid = 0;
      spx_alloc(addedcolid, 1);
      LPColSetReal compSlackCol;
      compSlackCol.add(1.0, -infinity, slackColCoeff, infinity);
      _compSolver.addCols(addedcolid, compSlackCol);
      //_compSlackColId = _compSolver.colId(_compSolver.nCols() - 1);
      _compSlackColId = addedcolid[0];

      int* perm = 0;
      spx_alloc(perm, numRowsReal() + addedRowsCounter);
      _compSolver.removeRows(rowsforremoval, nrowsforremoval, perm);

      // freeing allocated memory
      spx_free(perm);
      spx_free(addedcolid);
   }



   /// update the dual complementary problem with additional columns and rows
   // Given the solution to the updated reduced problem, the complementary problem will be updated with modifications to
   // the constraints and the removal of variables
   void SoPlex::_updateIdsComplementaryProblem(DVector dualVector)
   {
      int prevNumCols = _compSolver.nCols(); // number of columns in the previous formulation of the complementary prob
      int prevPrimalRowIds = _nPrimalRows;
      int prevDualColIds = _nDualCols;

      for( int i = 0; i < _nPrimalRows; i++ )
         printf("{%d,%d}(%d) ", _realLP->number(_idsPrimalRowIDs[i]), _solver.number(_idsPrimalRowIDs[i]),
               _idsReducedProbRows[_realLP->number(_idsPrimalRowIDs[i])]);
      printf("\n");

      LPColSetBase<Real> addElimCols(_nElimPrimalRows);  // columns previously eliminated from the
                                                         // complementary problem that must be added
      int numElimColsAdded = 0;
      // looping over all rows from the original LP that were eliminated during the formation of the complementary
      // problem. The eliminated rows will be added if they are basic in the reduced problem.
      for( int i = 0; i < _nElimPrimalRows; i++ )
      {
         if( _solver.isBasic(_idsElimPrimalRowIDs[i]) )
         {
            if( GT(dualVector[_solver.number(SPxRowId(_idsElimPrimalRowIDs[i]))], 0.0) )
            {
               // NOTE: This will probably need the SPxRowId passed to the function to get the new row id.
               addElimCols.add(_realLP->rhs(_idsElimPrimalRowIDs[i]), -infinity,
                     _realLP->rowVector(_idsElimPrimalRowIDs[i]), infinity);

               _idsPrimalRowIDs[_nPrimalRows] = _idsElimPrimalRowIDs[i];
               _nPrimalRows++;

               _idsElimPrimalRowIDs.remove(i);
               _nElimPrimalRows--;
               i--;

               numElimColsAdded++;
            }
            else if( LT(dualVector[_solver.number(SPxRowId(_idsElimPrimalRowIDs[i]))], 0.0) )
            {
               addElimCols.add(_realLP->lhs(_idsElimPrimalRowIDs[i]), -infinity,
                     _realLP->rowVector(_idsElimPrimalRowIDs[i]), infinity);

               _idsPrimalRowIDs[_nPrimalRows] = _idsElimPrimalRowIDs[i];
               _nPrimalRows++;

               _idsElimPrimalRowIDs.remove(i);
               _nElimPrimalRows--;
               i--;

               numElimColsAdded++;
            }
         }
      }

      printf("Num added cols: %d\n", numElimColsAdded);

      // updating the _idsDualColIDs with the additional columns from the eliminated rows.
      _compSolver.addCols(addElimCols);
      for( int i = prevNumCols; i < _compSolver.nCols(); i++ )
         _idsDualColIDs[prevDualColIds + i - prevNumCols] = _compSolver.colId(i);

      _nDualCols = _nPrimalRows;

      for( int i = 0; i < _nPrimalRows; i++ )
         printf("%d ", _solver.number(_idsPrimalRowIDs[i]));
      printf("\n");

      // looping over all rows from the original problem that were originally contained in the complementary problem.
      // The basic rows will be set as free variables, the non-basic rows will be eliminated from the complementary
      // problem.
      DSVector slackRowCoeff(_compSolver.nCols());
      Real slackCoeff = 1.0;

      int* colsforremoval = 0;
      int ncolsforremoval = 0;
      spx_alloc(colsforremoval, prevPrimalRowIds);
      for( int i = 0; i < prevPrimalRowIds; i++ )
      {
         int rowNumber = _realLP->number(_idsPrimalRowIDs[i]);
         if( _idsReducedProbRows[rowNumber] )
         {
            // rows added to the reduced problem may have been equality constriants. The equality constraints from the
            // original problem are converted into <= and >= constraints. Upon adding these constraints to the reduced
            // problem, only a single dual column is needed in the complementary problem. Hence, one of the dual columns
            // is removed.
            if( _realLP->number(_idsPrimalRowIDs[i]) == _realLP->number(_idsPrimalRowIDs[i+1]) )
            {
               colsforremoval[ncolsforremoval] = _compSolver.number(SPxColId(_idsDualColIDs[i + 1]));
               printf("Remove Col: %d\n", _compSolver.number(SPxColId(_idsDualColIDs[i + 1])));
               ncolsforremoval++;

               _idsPrimalRowIDs.remove(i + 1);
               _nPrimalRows--;
               _idsDualColIDs.remove(i + 1);
               _nDualCols--;

               prevPrimalRowIds--;
            }

            if( _solver.isBasic(_idsReducedProbRowIDs[rowNumber]) )
            {
               if( GT(dualVector[_solver.number(_idsReducedProbRowIDs[rowNumber])], 0.0) )
               {
                  _compSolver.changeObj(_idsDualColIDs[i], _realLP->rhs(SPxRowId(_idsPrimalRowIDs[i])));
                  _compSolver.changeBounds(_idsDualColIDs[i], -infinity, infinity);
               }
               else if( LT(dualVector[_solver.number(_idsReducedProbRowIDs[rowNumber])], 0.0) )
               {
                  _compSolver.changeObj(_idsDualColIDs[i], _realLP->lhs(SPxRowId(_idsPrimalRowIDs[i])));
                  _compSolver.changeBounds(_idsDualColIDs[i], -infinity, infinity);
               }
               else
               {
                  assert(isZero(dualVector[_solver.number(_idsReducedProbRowIDs[rowNumber])]));
                  colsforremoval[ncolsforremoval] = _compSolver.number(SPxColId(_idsDualColIDs[i]));
                  printf("Remove Col: %d\n", _compSolver.number(SPxColId(_idsDualColIDs[i])));
                  ncolsforremoval++;

                  _idsElimPrimalRowIDs[_nElimPrimalRows] = _idsPrimalRowIDs[i];
                  _nElimPrimalRows++;
                  _idsPrimalRowIDs.remove(i);
                  _nPrimalRows--;
                  _idsDualColIDs.remove(i);
                  _nDualCols--;
               }
            }
         }
         else
         {
            switch( _realLP->rowType(_idsPrimalRowIDs[i]) )
            {
               case LPRowBase<Real>::RANGE:
                  assert(_realLP->number(SPxColId(_idsPrimalRowIDs[i])) ==
                        _realLP->number(SPxColId(_idsPrimalRowIDs[i+1])));
                  if( _compSolver.obj(_compSolver.number(SPxColId(_idsDualColIDs[i]))) <
                        _compSolver.obj(_compSolver.number(SPxColId(_idsDualColIDs[i + 1]))))
                  {
                     slackRowCoeff.add(_compSolver.number(SPxColId(_idsDualColIDs[i])), -slackCoeff);
                     slackRowCoeff.add(_compSolver.number(SPxColId(_idsDualColIDs[i + 1])), slackCoeff);
                  }
                  else
                  {
                     slackRowCoeff.add(_compSolver.number(SPxColId(_idsDualColIDs[i])), slackCoeff);
                     slackRowCoeff.add(_compSolver.number(SPxColId(_idsDualColIDs[i + 1])), -slackCoeff);
                  }
                  i++;
                  break;
               case LPRowBase<Real>::EQUAL:
                  assert(_realLP->number(SPxColId(_idsPrimalRowIDs[i])) ==
                        _realLP->number(SPxColId(_idsPrimalRowIDs[i+1])));

                  slackRowCoeff.add(_compSolver.number(SPxColId(_idsDualColIDs[i])), slackCoeff);
                  slackRowCoeff.add(_compSolver.number(SPxColId(_idsDualColIDs[i + 1])), slackCoeff);

                  i++;
                  break;
               case LPRowBase<Real>::GREATER_EQUAL:
                  slackRowCoeff.add(_compSolver.number(SPxColId(_idsDualColIDs[i])), -slackCoeff);
                  break;
               case LPRowBase<Real>::LESS_EQUAL:
                  slackRowCoeff.add(_compSolver.number(SPxColId(_idsDualColIDs[i])), slackCoeff);
                  break;
               default:
                  throw SPxInternalCodeException("XLPFRD01 This should never happen.");
            }
         }
      }

      // updating the slack column in the complementary problem
      LPRowBase<Real> compSlackRow(1.0, slackRowCoeff, 1.0);
      _compSolver.changeRow(_compSlackDualRowId, compSlackRow);


      int* perm = 0;
      spx_alloc(perm, _compSolver.nCols() + numElimColsAdded);
      _compSolver.removeCols(colsforremoval, ncolsforremoval, perm);

      // freeing allocated memory
      spx_free(perm);
   }



   /// checking the optimality of the original problem.
   // this function is called if the complementary problem is solved with a non-negative objective value. This implies
   // that the rows currently included in the reduced problem are sufficient to identify the optimal solution to the
   // original problem.
   void SoPlex::_checkOriginalProblemOptimality()
   {

   }



   /// updating the slack column coefficients to adjust for equality constraints
   void SoPlex::_updateComplementarySlackColCoeff()
   {
      Real slackCoeff = 1.0;

      // the slack column for the equality constraints is not handled correctly in the dual conversion. Hence, it is
      // necessary to change the equality coefficients of the dual row related to the slack column.
      for( int i = 0; i < _nPrimalRows; i++ )
      {
         int rowNumber = _realLP->number(SPxRowId(_idsPrimalRowIDs[i]));
         if( !_idsReducedProbRows[rowNumber] && _realLP->rowType(_idsPrimalRowIDs[i]) == LPRowBase<Real>::EQUAL )
         {
            assert(_realLP->lhs(_idsPrimalRowIDs[i]) == _realLP->rhs(_idsPrimalRowIDs[i]));
            _compSolver.changeLower(_idsDualColIDs[i], 0.0);

            LPColBase<Real> addEqualityCol(_realLP->rhs(_idsPrimalRowIDs[i]),
                  -1.0*_compSolver.colVector(_idsDualColIDs[i]), infinity, 0.0);
            SPxColId newDualCol;
            _compSolver.addCol(newDualCol, addEqualityCol);

            _idsPrimalRowIDs.insert(i + 1, 1, _idsPrimalRowIDs[i]);
            _idsDualColIDs.insert(i + 1, 1, newDualCol);

            _compSolver.changeElement(_compSlackDualRowId, _idsDualColIDs[i], slackCoeff);
            _compSolver.changeElement(_compSlackDualRowId, _idsDualColIDs[i + 1], slackCoeff);

            i++;
            _nPrimalRows++;
            _nDualCols++;
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
