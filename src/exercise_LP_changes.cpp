/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the class library                   */
/*       SoPlex --- the Sequential object-oriented simPlex.                  */
/*                                                                           */
/*    Copyright (C) 1997-1999 Roland Wunderling                              */
/*                  1997-2002 Konrad-Zuse-Zentrum                            */
/*                            fuer Informationstechnik Berlin                */
/*                                                                           */
/*  SoPlex is distributed under the terms of the ZIB Academic Licence.       */
/*                                                                           */
/*  You should have received a copy of the ZIB Academic License              */
/*  along with SoPlex; see the file COPYING. If not email to soplex@zib.de.  */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
#pragma ident "@(#) $Id: exercise_LP_changes.cpp,v 1.5 2005/11/09 13:53:49 bzforlow Exp $"

#include <assert.h>
#include <math.h>
#include <string.h>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>

#include "spxdefines.h"
#include "soplex.h"
#include "spxsolver.h"

#include "timer.h"
#include "spxpricer.h"
#include "spxdefaultpr.h"
#include "spxparmultpr.h"
#include "spxdevexpr.h"
#include "spxhybridpr.h"
#include "spxsteeppr.h"
#include "spxweightpr.h"
#include "spxratiotester.h"
#include "spxharrisrt.h"
#include "spxdefaultrt.h"
#include "spxfastrt.h"
#include "spxsimplifier.h"
#include "spxintervalsm.h"
#include "spxaggregatesm.h"
#include "spxredundantsm.h"
#include "spxgeneralsm.h"
#include "spxscaler.h"
#include "spxequilisc.h"
#include "spxgeometsc.h"
#include "spxsumst.h"
#include "spxweightst.h"
#include "spxvectorst.h"
#include "slufactor.h"
#include "spxout.h"

using namespace soplex;

//------------------------------------------------------------------------
//    class TestSolver
//------------------------------------------------------------------------

/**@class TestSolver
    @brief A simple derived class from SpxSolver, used as object under test.
 */
class TestSolver : public SPxSolver
{
public:

   //------------------------------------
   /**@name Default parameter values
      Parameters used for LP solving. They mainly correspond to the defaults of
      the main binary, except for the ones indicated by "// *".
      Parts of the default settings are also realized in the constructor below.
   */
   //@{
   static const SLUFactor::UpdateType update = SLUFactor::FOREST_TOMLIN;
   static const Real delta = DEFAULT_BND_VIOL;
   static const Real timelimit = -1.0;
   static const Real epsilon = DEFAULT_EPS_ZERO;
   static const Real epsilon_factor = DEFAULT_EPS_FACTOR;
   static const Real epsilon_update = DEFAULT_EPS_UPDATE;
   static const int verbose = 1;
   static const int precision = 12;
   //@}

private:
   //------------------------------------
   /**@name Data */
   //@{
   SLUFactor _solver;              ///< sparse LU factorization
   SPxSteepPR _pricer;             ///< steepest edge pricer
   SPxFastRT _ratiotester;         ///< Harris fast ratio tester
   //@}

public:

   //------------------------------------
   /**@name Construction / destruction */
   //@{
   /// Default constructor.
   explicit
   TestSolver( const SPxSolver::Type type_ = SPxSolver::LEAVE, 
               const SPxSolver::Representation representation_ = SPxSolver::COLUMN )
      : SPxSolver( type_, 
                   representation_ )
   {
      setDelta( delta  );
      setTerminationTime( timelimit );

      Param::setEpsilon( epsilon );
      Param::setEpsilonFactorization( epsilon_factor );
      Param::setEpsilonUpdate( epsilon_update );
      Param::setVerbose( verbose );

      setPricer( &_pricer );
      setTester( &_ratiotester );
      setSolver( &_solver );
      _solver.setUtype( update );
      // no starter, no simplifier, no scaler

      assert( isConsistent() );
   }
   //@}
};

//
// Define static members of "TestSolver".
//
const SLUFactor::UpdateType TestSolver::update;
const Real TestSolver::delta;
const Real TestSolver::timelimit;
const Real TestSolver::epsilon;
const Real TestSolver::epsilon_factor;
const Real TestSolver::epsilon_update;
const int TestSolver::verbose;
const int TestSolver::precision;


//------------------------------------------------------------------------
//    class ChangeExerciser
//------------------------------------------------------------------------

/** 
    Class implementing the tests for the LP-changing interface.
 */
class ChangeExerciser
{
   /**
      Precision used in (relative) equality checks for differing but logically 
      equivalent LP solutions.
   */
   static const Real epsilon_solution_equal = 1e-9;

public:

   //------------------------------------
   /**@name Construction / destruction */
   //@{
   /// Default constructor.
   explicit
   ChangeExerciser( const std::string& instance_name )
      : _asserts_failed( 0 )
      , _instance_name( instance_name )
   {}
   //@}

public:

   //------------------------------------
   /**@name Test methods */
   //@{
   void test_add_delete_row();
   void test_add_delete_rows();
   void test_add_delete_col();
   void test_add_delete_cols();

   void test_change_obj();
   void test_change_lower();
   void test_change_upper();
   void test_change_bounds();
   void test_change_lhs();
   void test_change_rhs();
   void test_change_range();
   void test_change_sense();

   /// Reset test statistics.
   void reset()
   {
      _asserts_failed = 0;
   }

   /// Number of failed asserts.
   long asserts_failed() const { return _asserts_failed; }
   //@}

private:

   //------------------------------------
   /**@name Testing support */
   //@{
   ///
   void _assert( const std::string& description, const bool condition )
   {
      if ( !condition )
         {
            ++_asserts_failed;
            std::cout << "check '" << description << "' failed" << std::endl;
         }
   }
   ///
   void _assert_EQrel( const std::string& description, const Real ref, const Real val )
   {
      if ( !EQrel( ref, val, epsilon_solution_equal ) )
         {
            ++_asserts_failed;
            std::cout << "check '" << description << "' failed: expected " 
                      << ref << ", got " << val << std::endl;
         }
   }

   ///
   TestSolver* _prepare_Solver() const 
   {
      TestSolver* work_ptr = new TestSolver;
      work_ptr->readFile( _instance_name.c_str(), 0, 0 );
      return work_ptr;
   }

   ///
   long _asserts_failed;
   ///
   std::string _instance_name;
   //@}   
};

//
// Define static members of "ChangeExerciser".
//
const Real ChangeExerciser::epsilon_solution_equal;



//------------------------------------------------------------------------
//    main program
//------------------------------------------------------------------------

int main( int argc, 
          const char* const argv[] )
{
   const char* usage =
   "<test_suite>\n\n"
   "<test_suite> is supposed to be file containing names of LPs in either MPS or LPF format\n\n"
   ;

   if ( argc != 2 )
      {
         std::cout << "usage: " << argv[0] << " " << usage << std::endl;
         exit(0);
      }

   const char* suite_name = argv[1];

   std::cout.setf( std::ios::scientific | std::ios::showpoint );
   std::cerr.setf( std::ios::scientific | std::ios::showpoint );

   std::ifstream test_suite( suite_name );
   if ( !test_suite )
   {
      MSG_ERROR( spxout << "error while reading test suite file \""
                        << suite_name << "\"" << std::endl; )
      exit(1);
   }
   
   //
   // Do the actual work.
   //
   long total_asserts_failed = 0;

   while ( test_suite )
      {
         std::string filename;
         std::getline( test_suite, filename );

         if ( filename == "" ) break;
         std::cout << filename << ": " << std::endl;

         // Double-check that instance exists and is readable.
         TestSolver work;

         if ( !work.readFile( filename.c_str(), 0, 0 ) )
            {
               MSG_ERROR( spxout << "error while reading file \"" 
                                 << filename << "\"" << std::endl; )
               break;
            }

         // Do testing.
         ChangeExerciser tester( filename );
         tester.test_add_delete_row();
         tester.test_add_delete_rows();
         tester.test_add_delete_col();
         tester.test_add_delete_cols();

         tester.test_change_obj();
         tester.test_change_lower();
         tester.test_change_upper();
         tester.test_change_bounds();
         tester.test_change_lhs();
         tester.test_change_rhs();
         tester.test_change_range();
         tester.test_change_sense();

         std::cout << tester.asserts_failed() << " asserts failed.\n" << std::endl;
         total_asserts_failed += tester.asserts_failed();
      }

   std::cout << "Total number of failed asserts: " << total_asserts_failed << std::endl;

   return 0;
}


/**
   Deletes and re-adds a single row.
*/
void ChangeExerciser::test_add_delete_row()
{
   TestSolver* work_ptr = _prepare_Solver();

   work_ptr->solve();
   const Real original_obj = work_ptr->objValue();

   for( int row_idx = 0; row_idx != work_ptr->nRows(); ++row_idx )
      {
         const SVector& row_vector = work_ptr->rowVector( row_idx );
         const int non_zeroes = row_vector.size();
  
         LPRow current_row( non_zeroes );

         // Access rows alternatingly by index and ID.
         if ( row_idx % 2 == 0 )
            {
               work_ptr->getRow( row_idx, current_row );
               work_ptr->removeRow( row_idx );
               work_ptr->solve();
               work_ptr->addRow( current_row );
            }
         else
            {
               const SPxRowId& row_ID = work_ptr->rId( row_idx );
               work_ptr->getRow( row_ID, current_row );               
               work_ptr->removeRow( row_ID );
               work_ptr->solve();

               SPxRowId new_row_ID;
               work_ptr->addRow( new_row_ID, current_row );
            }
         work_ptr->solve();

         std::ostringstream description;
         description << "remove+add row " << row_idx;
         _assert_EQrel( description.str(), original_obj, work_ptr->objValue() );
      }

   delete work_ptr;
}


/**
   Deletes and re-adds sets of rows.
*/
void ChangeExerciser::test_add_delete_rows()
{
   TestSolver* work_ptr = _prepare_Solver();
   work_ptr->solve();
   const int orig_rows = work_ptr->nRows();
   const Real original_obj = work_ptr->objValue();

   // First test: Remove the middle third of the rows and reinsert it again.
   LPRowSet row_set;
   const int start = work_ptr->nRows() / 3;
   const int end = ( 2 * work_ptr->nRows() ) / 3;

   work_ptr->getRows( start, end, row_set );

   work_ptr->removeRowRange( start, end );
   work_ptr->solve();

   // Use addRows()-version returning row IDs.
   SPxRowId* row_IDs = new SPxRowId[ row_set.num() ];
   work_ptr->addRows( row_IDs, row_set );
   work_ptr->solve();
   delete row_IDs; row_IDs = 0;

   _assert( "remove+add rows: row number matches", orig_rows == work_ptr->nRows() );

   std::ostringstream description;
   description << "remove+add range of rows " << start << " to " << end;
   _assert_EQrel( description.str(), original_obj, work_ptr->objValue() );

   delete work_ptr; work_ptr = 0;

   // Second test: Remove a (more or less) random third of the rows and reinsert it.
   //      Removal is done once based on the absolute indices of the rows.
   work_ptr = _prepare_Solver();
   srand( 42 );

   int nr_deleted_rows = work_ptr->nRows() / 3;
   int* nums = new int[ nr_deleted_rows ];

   // We need to record used rows to avoid duplicates in row_set below.
   bool* used = new bool[ work_ptr->nRows() ];
   for ( int row_idx = 0; row_idx < work_ptr->nRows(); ++row_idx )
      {
         used[ row_idx ] = false;
      }

   // Mark rows for deletion and remember them.
   row_set.clear();
   for ( int i = 0; i < nr_deleted_rows; ++i )
      {
         int random_row_idx = rand() % work_ptr->nRows();
         while ( used[ random_row_idx ] )
            {
               random_row_idx = rand() % work_ptr->nRows();
            }

         nums[ i ] = random_row_idx;
         used[ random_row_idx ] = true;

         LPRow current_row( work_ptr->rowVector( random_row_idx ).size() );
         work_ptr->getRow( random_row_idx, current_row );
         row_set.add( current_row );
      }
   delete used; used = 0;

   work_ptr->removeRows( nums, nr_deleted_rows );
   work_ptr->solve();

   work_ptr->addRows( row_set );
   work_ptr->solve();

   _assert( "remove+add rows (idx): row number matches", orig_rows == work_ptr->nRows() );
   _assert_EQrel( "remove+add (idx) random set of rows", original_obj, work_ptr->objValue() );

   delete nums; nums = 0;
   delete work_ptr; work_ptr = 0;

   // Third test: As before, but removal based based on SPxRowIds.
   work_ptr = _prepare_Solver();
   nr_deleted_rows = work_ptr->nRows() / 3;
   SPxRowId* IDs = new SPxRowId[ nr_deleted_rows ];

   used = new bool[ work_ptr->nRows() ];
   for ( int row_idx = 0; row_idx < work_ptr->nRows(); ++row_idx )
      {
         used[ row_idx ] = false;
      }

   // Mark rows for deletion and remember them.
   row_set.clear();
   for ( int i = 0; i < nr_deleted_rows; ++i )
      {
         int random_row_idx = rand() % work_ptr->nRows();
         while ( used[ random_row_idx ] )
            {
               random_row_idx = rand() % work_ptr->nRows();
            }

         IDs[ i ] = work_ptr->rId( random_row_idx );
         used[ random_row_idx ] = true;

         LPRow current_row( work_ptr->rowVector( random_row_idx ).size() );
         work_ptr->getRow( random_row_idx, current_row );
         row_set.add( current_row );
      }
   delete used; used = 0;

   work_ptr->removeRows( IDs, nr_deleted_rows );
   work_ptr->solve();

   work_ptr->addRows( row_set );
   work_ptr->solve();

   _assert( "remove+add rows (ID): row number matches", orig_rows == work_ptr->nRows() );
   _assert_EQrel( "remove+add (ID) random set of rows", original_obj, work_ptr->objValue() );

   delete IDs; IDs = 0;
   delete work_ptr; work_ptr = 0;
}


//
// The following two routines are the analogue of the above for columns.
// In fact they have been obtained by a replace and it is a good idea to keep 
// them in sync.
//

/**
   Deletes and re-adds a single col.
*/
void ChangeExerciser::test_add_delete_col()
{
   TestSolver* work_ptr = _prepare_Solver();

   work_ptr->solve();
   const Real original_obj = work_ptr->objValue();

   for( int col_idx = 0; col_idx != work_ptr->nCols(); ++col_idx )
      {
         const SVector& col_vector = work_ptr->colVector( col_idx );
         const int non_zeroes = col_vector.size();
  
         LPCol current_col( non_zeroes );

         // Access cols alternatingly by index and ID.
         if ( col_idx % 2 == 0 )
            {
               work_ptr->getCol( col_idx, current_col );
               work_ptr->removeCol( col_idx );
               work_ptr->solve();
               work_ptr->addCol( current_col );
            }
         else
            {
               const SPxColId& col_ID = work_ptr->cId( col_idx );
               work_ptr->getCol( col_ID, current_col );

               work_ptr->removeCol( col_ID );
               work_ptr->solve();

               SPxColId new_col_ID;
               work_ptr->addCol( new_col_ID, current_col );
            }
         work_ptr->solve();

         std::ostringstream description;
         description << "remove+add col " << col_idx;
         _assert_EQrel( description.str(), original_obj, work_ptr->objValue() );
      }

   delete work_ptr;
}


/**
   Deletes and re-adds sets of cols.
*/
void ChangeExerciser::test_add_delete_cols()
{
   TestSolver* work_ptr = _prepare_Solver();
   work_ptr->solve();
   const int orig_cols = work_ptr->nCols();
   const Real original_obj = work_ptr->objValue();

   // First test: Remove the middle third of the cols and reinsert it again.
   LPColSet col_set;
   const int start = work_ptr->nCols() / 3;
   const int end = ( 2 * work_ptr->nCols() ) / 3;

   work_ptr->getCols( start, end, col_set );

   work_ptr->removeColRange( start, end );
   work_ptr->solve();

   // Use addCols()-version returning col IDs.
   SPxColId* col_IDs = new SPxColId[ col_set.num() ];
   work_ptr->addCols( col_IDs, col_set );
   work_ptr->solve();
   delete col_IDs; col_IDs = 0;

   _assert( "remove+add cols: col number matches", orig_cols == work_ptr->nCols() );

   std::ostringstream description;
   description << "remove+add range of cols " << start << " to " << end;
   _assert_EQrel( description.str(), original_obj, work_ptr->objValue() );

   delete work_ptr; work_ptr = 0;

   // Second test: Remove a (more or less) random third of the cols and reinsert it.
   //      Removal is done once based on the absolute indices of the cols.
   work_ptr = _prepare_Solver();
   srand( 42 );

   int nr_deleted_cols = work_ptr->nCols() / 3;
   int* nums = new int[ nr_deleted_cols ];

   // We need to record used cols to avoid duplicates in col_set below.
   bool* used = new bool[ work_ptr->nCols() ];
   for ( int col_idx = 0; col_idx < work_ptr->nCols(); ++col_idx )
      {
         used[ col_idx ] = false;
      }

   // Mark cols for deletion and remember them.
   col_set.clear();
   for ( int i = 0; i < nr_deleted_cols; ++i )
      {
         int random_col_idx = rand() % work_ptr->nCols();
         while ( used[ random_col_idx ] )
            {
               random_col_idx = rand() % work_ptr->nCols();
            }

         nums[ i ] = random_col_idx;
         used[ random_col_idx ] = true;

         LPCol current_col( work_ptr->colVector( random_col_idx ).size() );
         work_ptr->getCol( random_col_idx, current_col );
         col_set.add( current_col );
      }
   delete used; used = 0;

   work_ptr->removeCols( nums, nr_deleted_cols );
   work_ptr->solve();

   work_ptr->addCols( col_set );
   work_ptr->solve();

   _assert( "remove+add cols (idx): col number matches", orig_cols == work_ptr->nCols() );
   _assert_EQrel( "remove+add (idx) random set of cols", original_obj, work_ptr->objValue() );

   delete nums; nums = 0;
   delete work_ptr; work_ptr = 0;

   // Third test: As before, but removal based based on SPxColIds.
   work_ptr = _prepare_Solver();
   nr_deleted_cols = work_ptr->nCols() / 3;
   SPxColId* IDs = new SPxColId[ nr_deleted_cols ];

   used = new bool[ work_ptr->nCols() ];
   for ( int col_idx = 0; col_idx < work_ptr->nCols(); ++col_idx )
      {
         used[ col_idx ] = false;
      }

   // Mark cols for deletion and remember them.
   col_set.clear();
   for ( int i = 0; i < nr_deleted_cols; ++i )
      {
         int random_col_idx = rand() % work_ptr->nCols();
         while ( used[ random_col_idx ] )
            {
               random_col_idx = rand() % work_ptr->nCols();
            }

         IDs[ i ] = work_ptr->cId( random_col_idx );
         used[ random_col_idx ] = true;

         LPCol current_col( work_ptr->colVector( random_col_idx ).size() );
         work_ptr->getCol( random_col_idx, current_col );
         col_set.add( current_col );
      }
   delete used; used = 0;

   work_ptr->removeCols( IDs, nr_deleted_cols );
   work_ptr->solve();

   work_ptr->addCols( col_set );
   work_ptr->solve();

   _assert( "remove+add cols (ID): col number matches", orig_cols == work_ptr->nCols() );
   _assert_EQrel( "remove+add (ID) random set of cols", original_obj, work_ptr->objValue() );

   delete IDs; IDs = 0;
   delete work_ptr; work_ptr = 0;
}


/**
   Changes objective function.
*/
void ChangeExerciser::test_change_obj()
{
   TestSolver* work_ptr = _prepare_Solver();
   work_ptr->solve();
   const Real original_obj = work_ptr->objValue();

   // First test: Multiply every objective value by 2.0. using changeObj by index
   for ( int col_idx = 0; col_idx != work_ptr->nCols(); ++col_idx )
      {
         Real coeff = work_ptr->obj( col_idx );
         coeff *= 2.0;
         // use changeObj(int, Real)
         work_ptr->changeObj( col_idx, coeff );
         _assert( "objective entry changed by index", work_ptr->obj( col_idx ) == coeff );        
      }

   work_ptr->solve();
   _assert_EQrel( "change objective entry by index", 2.0 * original_obj, work_ptr->objValue() );

   // Second test: Divide every objective value by 2.0, using changeObj by ID
   for ( int col_idx = 0; col_idx != work_ptr->nCols(); ++col_idx )
      {
         const SPxColId& col_ID = work_ptr->cId( col_idx );

         Real coeff = work_ptr->obj( col_ID );
         coeff = coeff/2.0;

         // use changeObj(int, Real)
         work_ptr->changeObj( col_ID, coeff );
         _assert( "objective entry changed by ID", work_ptr->obj( col_ID ) == coeff );        
      }
   
   work_ptr->solve();
   _assert_EQrel( "change objective entry by ID", original_obj, work_ptr->objValue() );

   // Third test: Multiply every objective value by 2.0. using changeObj by Vector
   Real* v = new Real[work_ptr->nCols()];
   Vector original(work_ptr->nCols(), v);

   work_ptr->getObj(original);

   original*=2.0;

   work_ptr->changeObj(original);

   work_ptr->solve();
   _assert_EQrel( "change objective entry by vector", 2.0 * original_obj, work_ptr->objValue() );

   original*=0.5;
   work_ptr->changeObj(original);
   work_ptr->solve();   
   _assert_EQrel( "change objective entry by vector", original_obj, work_ptr->objValue() );

   delete v; v = 0;
  
   delete work_ptr; work_ptr = 0;
}


/**
   Changes lower bounds
*/
void ChangeExerciser::test_change_lower()
{
   TestSolver* work_ptr = _prepare_Solver();
   work_ptr->solve();
   const Real original_obj = work_ptr->objValue();

   // First test : Mupltiply lower vector by 1.0
   Vector original_lower = work_ptr->lower();
   
   work_ptr->changeLower(original_lower*=1.0);

   work_ptr->solve();
   _assert_EQrel( "change lower by vector", original_obj, work_ptr->objValue() );

   // Second test:  get solution vector and compute slacks
   Real* val = new Real[work_ptr->nCols()];
   Vector solution(work_ptr->nCols(), val);
   work_ptr->getPrimal(solution);

   for(int col_idx = 0; col_idx != work_ptr->nCols(); ++col_idx){
      Real lower = work_ptr->lower(col_idx);
      
      // slack is positive
      if(solution[col_idx] > lower){
         // get ID
         const SPxColId& col_ID = work_ptr->cId( col_idx );
         work_ptr->changeLower(col_ID, solution[col_idx]);
         _assert( "lower value changed", work_ptr->lower( col_idx ) == solution[col_idx] );
      }
   }

   work_ptr->solve();
   _assert_EQrel( "change lower", original_obj, work_ptr->objValue() );

   delete val; val=0;
   delete work_ptr; work_ptr = 0;
}


/**
   Changes upper bounds
*/
void ChangeExerciser::test_change_upper()
{
   TestSolver* work_ptr = _prepare_Solver();
   work_ptr->solve();
   const Real original_obj = work_ptr->objValue();

   // First test : Mupltiply lower vector by 1.0
   Vector original_upper = work_ptr->upper();
   
   work_ptr->changeUpper(original_upper*=1.0);

   work_ptr->solve();
   _assert_EQrel( "change upper", original_obj, work_ptr->objValue() );

   // Second test: get solution vector and compute slacks
   Real* val = new Real[work_ptr->nCols()];
   Vector solution(work_ptr->nCols(), val);
   work_ptr->getPrimal(solution);

   for(int col_idx = 0; col_idx != work_ptr->nCols(); ++col_idx){
      Real upper = work_ptr->upper(col_idx);
      
      // slack is bigger than 0
      if(upper > solution[col_idx]){
         // get ID
         const SPxColId& col_ID = work_ptr->cId( col_idx );
         work_ptr->changeUpper(col_ID, solution[col_idx]);
         _assert( "upper value changed", work_ptr->upper( col_idx ) == solution[col_idx] );
      }
   }

   work_ptr->solve();
   _assert_EQrel( "change upper", original_obj, work_ptr->objValue() );

   delete work_ptr; work_ptr = 0;
}


/**
   Changes bounds
*/
void ChangeExerciser::test_change_bounds()
{
   TestSolver* work_ptr = _prepare_Solver();
   work_ptr->solve();
   const Real original_obj = work_ptr->objValue();

   // First test: multiply the bounds as vectors by 1.0
   Vector original_lower = work_ptr->lower();
   Vector original_upper = work_ptr->upper();
   work_ptr->changeBounds(original_lower*=1.0, original_upper*=1.0);
   work_ptr->solve();
   _assert_EQrel( "change bounds by vector", original_obj, work_ptr->objValue() );

   // Second test: use positive slacks to change the bounds (see test_change_upper, test_change_lower)
   Real* val = new Real[work_ptr->nCols()];
   Vector solution(work_ptr->nCols(), val);
   work_ptr->getPrimal(solution);

   for(int col_idx = 0; col_idx != work_ptr->nCols(); ++col_idx){
      Real upper = work_ptr->upper(col_idx);
      Real lower = work_ptr->lower(col_idx);
  
      // slack is bigger than 0
      if(upper > solution[col_idx])
         upper = solution[col_idx];
      
      if(solution[col_idx] > lower)
         lower = solution[col_idx];
      
      // get ID
      const SPxColId& col_ID = work_ptr->cId( col_idx );
      work_ptr->changeBounds(col_ID, lower, upper);
      _assert( "upper bound changed", work_ptr->upper( col_idx ) == upper );
      _assert( "lower bound changed", work_ptr->lower( col_idx ) == lower );
   }

   work_ptr->solve();
   _assert_EQrel( "change bounds", original_obj, work_ptr->objValue() );

   delete val; val = 0;
   delete work_ptr; work_ptr = 0;
}


/**
   Changes left hand side vector
*/
void ChangeExerciser::test_change_lhs()
{
   TestSolver* work_ptr = _prepare_Solver();
   work_ptr->solve();
   const Real original_obj = work_ptr->objValue();    

   // First test: change lhs by vector

   // store lhs vector
   Vector lhs = work_ptr->lhs();
   
   // change lhs to original
   work_ptr->changeLhs(lhs);
   work_ptr->solve();
   _assert_EQrel("lhs changed by vector", original_obj, work_ptr->objValue());

   // Second test: change lhs by indices
   for(int row_idx = 0; row_idx < work_ptr->nRows(); ++row_idx){
      if(row_idx%2==0.0){
         Real lhs_val = work_ptr->lhs(row_idx);
         // get ID
         const SPxRowId& row_ID = work_ptr->rId( row_idx );
         work_ptr->changeLhs(row_ID, lhs_val);
      }      
   }
   work_ptr->solve();
   _assert_EQrel("lhs changed by ID", original_obj, work_ptr->objValue());

   delete work_ptr; work_ptr=0;  
}


/**
   Changes right hand side vector
*/
void ChangeExerciser::test_change_rhs()
{
   TestSolver* work_ptr = _prepare_Solver();
   work_ptr->solve();
   const Real original_obj = work_ptr->objValue();    

   // First test: change rhs by vector

   // store rhs vector
   Vector rhs = work_ptr->rhs();
   
   // change rhs to original
   work_ptr->changeRhs(rhs);
   work_ptr->solve();
   _assert_EQrel("rhs changed by vector", original_obj, work_ptr->objValue());

   // Second test: change lhs by indices
   for(int row_idx = 0; row_idx < work_ptr->nRows(); ++row_idx){
      if(row_idx%2==0.0){
         // get ID
         const SPxRowId& row_ID = work_ptr->rId( row_idx );
         Real rhs_val = work_ptr->rhs(row_idx);
         work_ptr->changeRhs(row_ID, rhs_val);
      }      
   }
   work_ptr->solve();
   _assert_EQrel("rhs changed", original_obj, work_ptr->objValue());

   delete work_ptr; work_ptr=0;
}


/**
   Changes range
*/
void ChangeExerciser::test_change_range()
{
   TestSolver* work_ptr = _prepare_Solver();
   work_ptr->solve();
   const Real original_obj = work_ptr->objValue();    

   // First test: change range by vector
   Vector rhs = work_ptr->rhs();
   Vector lhs = work_ptr->lhs();
   for(int row_idx = 0; row_idx < work_ptr->nRows(); ++row_idx){
      rhs[row_idx]*=1.0;
      lhs[row_idx]*=1.0;
   }
   work_ptr->changeRange(lhs, rhs);
   
   work_ptr->solve();
   _assert_EQrel("range changed by vector", original_obj, work_ptr->objValue());

   // Second test: change range by indices
   for(int row_idx = 0; row_idx < work_ptr->nRows(); ++row_idx){
      if(row_idx%2 == 0){
         Real newr = lhs[row_idx]*(-1.0); 
         Real newl = rhs[row_idx]*(-1.0);
         // get ID
         const SPxRowId& row_ID = work_ptr->rId( row_idx );
         work_ptr->changeRange(row_ID, newl, newr);
         _assert("range value changed", work_ptr->lhs(row_idx)==newl);
         _assert("range value changed", work_ptr->rhs(row_idx)==newr);
      }
   }
   work_ptr->solve();

   // reverse the change
   for(int row_idx=0; row_idx < work_ptr->nRows(); ++row_idx){
      if(row_idx%2 == 0){
         Real newr = lhs[row_idx]*(-1.0); 
         Real newl = rhs[row_idx]*(-1.0);
         // get ID
         const SPxRowId& row_ID = work_ptr->rId( row_idx );
         work_ptr->changeRange(row_ID, newl, newr);
         _assert("range value changed", work_ptr->lhs(row_idx)==newl);
         _assert("range value changed", work_ptr->rhs(row_idx)==newr);
      }
   }

   work_ptr->solve();
   _assert_EQrel("range changed", original_obj, work_ptr->objValue());

   delete work_ptr; work_ptr=0;
}


/**
   Changes sense
*/
void ChangeExerciser::test_change_sense()
{
   TestSolver* work_ptr = _prepare_Solver();
   work_ptr->solve();
   const Real original_obj = work_ptr->objValue();    
   
   // change the sense
   if(work_ptr->sense()==SPxLP::MINIMIZE){    
      SPxLP::SPxSense max = SPxLP::MAXIMIZE;
      work_ptr->changeSense(max);
   }
   else if(work_ptr->sense()==SPxLP::MAXIMIZE){     
      SPxLP::SPxSense min = SPxLP::MINIMIZE;
      work_ptr->changeSense(min);
   }
   
   // change objective function (multiply by -1) 
   Real* v = new Real[work_ptr->nCols()];
   Vector original(work_ptr->nCols(), v);

   work_ptr->getObj(original);   
   original*=-1.0;
   work_ptr->changeObj(original);
   
   work_ptr->solve();
   
   // check failure
   _assert_EQrel( "change sense", original_obj*-1.0, work_ptr->objValue() );
   
   delete work_ptr; work_ptr=0;
   delete v; v=0;
}

//-----------------------------------------------------------------------------
//Emacs Local Variables:
//Emacs mode:c++
//Emacs c-basic-offset:3
//Emacs tab-width:8
//Emacs indent-tabs-mode:nil
//Emacs End:
//-----------------------------------------------------------------------------
