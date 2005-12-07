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
#pragma ident "@(#) $Id: exercise_LP_changes.cpp,v 1.8 2005/12/07 18:03:24 bzfhille Exp $"

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
   static const Real delta;
   static const Real timelimit;
   static const Real epsilon;
   static const Real epsilon_factor;
   static const Real epsilon_update;
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
   TestSolver( const SPxSolver::Type type_, 
               const SPxSolver::Representation representation_ )
      : SPxSolver( type_, 
                   representation_ )
   {
      setDelta( delta  );
      setTerminationTime( timelimit );

      Param::setEpsilon( epsilon );
      Param::setEpsilonFactorization( epsilon_factor );
      Param::setEpsilonUpdate( epsilon_update );

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
const Real TestSolver::delta = DEFAULT_BND_VIOL;
const Real TestSolver::timelimit = -1.0;
const Real TestSolver::epsilon = DEFAULT_EPS_ZERO;
const Real TestSolver::epsilon_factor = DEFAULT_EPS_FACTOR;
const Real TestSolver::epsilon_update = DEFAULT_EPS_UPDATE;
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
   static const Real epsilon_solution_equal;

public:

   //------------------------------------
   /**@name Construction / destruction */
   //@{
   /// Default constructor.
   explicit
   ChangeExerciser( const std::string& instance_name,
                    const SPxSolver::Type type,
                    const SPxSolver::Representation representation )
      : _asserts_failed( 0 )
      , _instance_name( instance_name )
      , _type( type )
      , _representation( representation )
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
   void test_change_row();
   void test_change_col();
   void test_change_element();
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
      TestSolver* work_ptr = new TestSolver( _type, _representation );
      work_ptr->readFile( _instance_name.c_str(), 0, 0 );
      return work_ptr;
   }

   ///
   long _asserts_failed;
   ///
   std::string _instance_name;
   ///
   const SPxSolver::Type _type;
   ///
   const SPxSolver::Representation _representation;
   //@}   
};

//
// Define static members of "ChangeExerciser".
//
const Real ChangeExerciser::epsilon_solution_equal = 1e-9;


//------------------------------------------------------------------------
//    main program
//------------------------------------------------------------------------

/// Global verbosity indicator.
bool verbose = false;

/// Prints message only if verbose mode is on.
#define MSG( message ) { if ( verbose ) { message } }


/**
   Wrapper to run all tests for  fixed representation and algorithm.
   Returns the number of failed asserts.
*/
long run_tests( const std::string& filename, 
                const SPxSolver::Type type,
                const SPxSolver::Representation representation )
{
   ChangeExerciser tester( filename,
                           type,
                           representation );

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
   tester.test_change_row();
   tester.test_change_col();
   tester.test_change_element();
   tester.test_change_sense();
   
   std::cout << representation << ", " << type << ": " << tester.asserts_failed() 
             << " asserts failed." << std::endl;

   return tester.asserts_failed();
}


int main( int argc, 
          const char* const argv[] )
{
   const char* usage =
      "[-vLevel] <test_suite>\n\n"
      "<test_suite> is supposed to be file containing names of LPs in either MPS or LPF format\n\n"
      "-v        activate verbose mode\n"
   ;

   // Cope with command line.
   if ( argc < 2 || argc > 3 )
      {
         std::cout << "usage: " << argv[0] << " " << usage << std::endl;
         exit(0);
      }

   int optidx = 1;
   for (optidx = 1; optidx < argc; optidx++)
   {
      if (*argv[optidx] != '-')
         break;

      switch(argv[optidx][1])
      {
      case 'v' :
         verbose = true;
         break;
      default :
         std::cout << "usage: " << argv[0] << " " << usage << std::endl;
         exit(0);
      }
   }

   const char* suite_name = argv[ optidx ];

   std::ifstream test_suite( suite_name );
   if ( !test_suite )
   {
      std::cerr << "error while reading test suite file \""
                << suite_name << "\"" << std::endl;
      exit(1);
   }

   // Setup.
   std::cout.setf( std::ios::scientific | std::ios::showpoint );
   std::cerr.setf( std::ios::scientific | std::ios::showpoint );

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
         TestSolver work( SPxSolver::ENTER, SPxSolver::COLUMN );

         if ( !work.readFile( filename.c_str(), 0, 0 ) )
            {
               std::cout << "error while reading file \"" 
                         << filename << "\"" << std::endl;
               break;
            }

         // Do testing.
         total_asserts_failed += run_tests( filename, SPxSolver::LEAVE, SPxSolver::COLUMN );
         total_asserts_failed += run_tests( filename, SPxSolver::LEAVE, SPxSolver::ROW );
         total_asserts_failed += run_tests( filename, SPxSolver::ENTER, SPxSolver::COLUMN );
         total_asserts_failed += run_tests( filename, SPxSolver::ENTER, SPxSolver::ROW );

         std::cout << std::endl;
      }

   std::cout << "Total number of failed asserts: " << total_asserts_failed << std::endl;

   return 0;
}

//------------------------------------------------------------------------
//    ChangeExerciser class
//------------------------------------------------------------------------

/**
   Deletes and re-adds a single row.
*/
void ChangeExerciser::test_add_delete_row()
{
   MSG( std::cout << "Testing addRow() / removeRow()" << std::endl; );

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
   MSG( std::cout << "Testing addRows() / removeRows()" << std::endl; );

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
   delete[] row_IDs; row_IDs = 0;

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
   delete[] used; used = 0;

   work_ptr->removeRows( nums, nr_deleted_rows );
   work_ptr->solve();

   work_ptr->addRows( row_set );
   work_ptr->solve();

   _assert( "remove+add rows (idx): row number matches", orig_rows == work_ptr->nRows() );
   _assert_EQrel( "remove+add (idx) random set of rows", original_obj, work_ptr->objValue() );

   delete[] nums; nums = 0;
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
   delete[] used; used = 0;

   work_ptr->removeRows( IDs, nr_deleted_rows );
   work_ptr->solve();

   work_ptr->addRows( row_set );
   work_ptr->solve();

   _assert( "remove+add rows (ID): row number matches", orig_rows == work_ptr->nRows() );
   _assert_EQrel( "remove+add (ID) random set of rows", original_obj, work_ptr->objValue() );

   delete[] IDs; IDs = 0;
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
   MSG( std::cout << "Testing addCol() / removeCol()" << std::endl; );

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
   MSG( std::cout << "Testing addCols() / removeCols()" << std::endl; );

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
   delete[] used; used = 0;

   work_ptr->removeCols( nums, nr_deleted_cols );
   work_ptr->solve();

   work_ptr->addCols( col_set );
   work_ptr->solve();

   _assert( "remove+add cols (idx): col number matches", orig_cols == work_ptr->nCols() );
   _assert_EQrel( "remove+add (idx) random set of cols", original_obj, work_ptr->objValue() );

   delete[] nums; nums = 0;
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
   delete[] used; used = 0;

   work_ptr->removeCols( IDs, nr_deleted_cols );
   work_ptr->solve();

   work_ptr->addCols( col_set );
   work_ptr->solve();

   _assert( "remove+add cols (ID): col number matches", orig_cols == work_ptr->nCols() );
   _assert_EQrel( "remove+add (ID) random set of cols", original_obj, work_ptr->objValue() );

   delete[] IDs; IDs = 0;
   delete work_ptr; work_ptr = 0;
}


/**
   Changes objective function.
*/
void ChangeExerciser::test_change_obj()
{
   MSG( std::cout << "Testing changeObj()" << std::endl; );
   
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

   delete[] v; v = 0;
   delete work_ptr; work_ptr = 0;
}


/**
   Changes lower bounds
*/
void ChangeExerciser::test_change_lower()
{
   MSG( std::cout << "Testing changeLower()" << std::endl; );

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

   delete[] val; val=0;
   delete work_ptr; work_ptr = 0;
}


/**
   Changes upper bounds
*/
void ChangeExerciser::test_change_upper()
{
   MSG( std::cout << "Testing changeUpper()" << std::endl; );

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

   delete[] val; val=0;
   delete work_ptr; work_ptr = 0;
}


/**
   Changes bounds
*/
void ChangeExerciser::test_change_bounds()
{
   MSG( std::cout << "Testing changeBounds()" << std::endl; );

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

   delete[] val; val = 0;
   delete work_ptr; work_ptr = 0;
}


/**
   Changes left hand side vector
*/
void ChangeExerciser::test_change_lhs()
{
   MSG( std::cout << "Testing changeLhs()" << std::endl; );

   TestSolver* work_ptr = _prepare_Solver();
   work_ptr->solve();
   Real original_obj = work_ptr->objValue();    

   // First test: change lhs by vector

   // store lhs vector
   Vector lhs = work_ptr->lhs();
   
   // change lhs to original
   work_ptr->changeLhs( lhs );
   work_ptr->solve();
   _assert_EQrel( "lhs changed by vector", original_obj, work_ptr->objValue() );

   delete work_ptr; work_ptr = 0;

   // Second test: set new lhs to solution*rowVector  
   work_ptr = _prepare_Solver();
   work_ptr->solve();
   original_obj = work_ptr->objValue();    

   Real* val = new Real[ work_ptr->nCols() ];
   Vector solution( work_ptr->nCols(), val );
   work_ptr->getPrimal( solution );

   for (int row_idx = 0; row_idx < work_ptr->nRows(); ++row_idx)
      {
         const SPxRowId& row_ID = work_ptr->rId( row_idx );
                     
         const Real row_prod = solution * work_ptr->rowVector( row_ID );

         // The new LHS must not be larger than the RHS.
         const Real new_lhs = std::min< Real >(  row_prod,
                                                 work_ptr->rhs(row_idx) );

         work_ptr->changeLhs( row_ID, new_lhs );
         
         _assert( "lhs updated", work_ptr->lhs( row_idx ) == new_lhs );
         _assert( "lhs validity ", work_ptr->lhs( row_idx ) <= work_ptr->rhs( row_idx ) );
      }
   
   work_ptr->solve();
   _assert_EQrel( "lhs changed by ID", original_obj, work_ptr->objValue() );
   
   delete[] val; val = 0;
   delete work_ptr; work_ptr = 0;
}


/**
   Changes right hand side vector
*/
void ChangeExerciser::test_change_rhs()
{
   MSG( std::cout << "Testing changeRhs()" << std::endl; );

   TestSolver* work_ptr = _prepare_Solver();
   work_ptr->solve();
   Real original_obj = work_ptr->objValue();    

   // First test: change rhs by vector

   // store rhs vector
   Vector rhs = work_ptr->rhs();
   
   // change rhs to original
   work_ptr->changeRhs( rhs );
   work_ptr->solve();
   _assert_EQrel( "rhs changed by vector", original_obj, work_ptr->objValue() );

   delete work_ptr; work_ptr = 0;

   // Second test: set new rhs to solution*rowVector  
   work_ptr = _prepare_Solver();
   work_ptr->solve();
   original_obj = work_ptr->objValue();    

   Real* val = new Real[ work_ptr->nCols() ];
   Vector solution( work_ptr->nCols(), val );
   work_ptr->getPrimal( solution );

   for (int row_idx = 0; row_idx < work_ptr->nRows(); ++row_idx)
      {
         const SPxRowId& row_ID = work_ptr->rId( row_idx );
                     
         const Real row_prod = solution * work_ptr->rowVector( row_ID );

         // The new RHS must not be smaller than the LHS.
         const Real new_rhs = std::max< Real >(  row_prod,
                                                 work_ptr->lhs(row_idx) );

         work_ptr->changeRhs( row_ID, new_rhs );
         
         _assert( "rhs updated", work_ptr->rhs( row_idx ) == new_rhs );
         _assert( "rhs validity ", work_ptr->rhs( row_idx ) <= work_ptr->rhs( row_idx ) );
      }
   
   work_ptr->solve();
   _assert_EQrel( "rhs changed by ID", original_obj, work_ptr->objValue() );
   
   delete[] val; val = 0;
   delete work_ptr; work_ptr = 0;
}


/**
   Changes range
*/
void ChangeExerciser::test_change_range()
{
   MSG( std::cout << "Testing changeRange()" << std::endl; );
   
   TestSolver* work_ptr = _prepare_Solver();
   work_ptr->solve();
   Real original_obj = work_ptr->objValue();    

   // First test: change range by vector
   Vector rhs = work_ptr->rhs();
   Vector lhs = work_ptr->lhs();
   for (int row_idx = 0; row_idx < work_ptr->nRows(); ++row_idx) 
   {
      rhs[ row_idx ] *= 2.0;
      lhs[ row_idx ] *= 2.0;
   }
   work_ptr->changeRange( lhs, rhs );
   
   work_ptr->solve();
   _assert_EQrel( "range changed by vector", 2.0 * original_obj, work_ptr->objValue() );
   
   delete work_ptr; work_ptr=0;

   // Second test: change range by indices
   work_ptr = _prepare_Solver();
   work_ptr->solve();
   original_obj = work_ptr->objValue();

   for (int row_idx = 0; row_idx < work_ptr->nRows(); ++row_idx) 
   {
      const Real newl = -1.0 * rhs[row_idx];
      const Real newr = -1.0 * lhs[row_idx]; 

      work_ptr->changeRange( work_ptr->rId( row_idx ), newl, newr );
      _assert( "lhs range value changed", work_ptr->lhs( row_idx ) == newl );
      _assert( "rhs range value changed", work_ptr->rhs( row_idx ) == newr );
   }
   work_ptr->solve();

   // reverse the change
   for (int row_idx=0; row_idx < work_ptr->nRows(); ++row_idx)
   {
      const Real newl = -1.0 * rhs[row_idx];
      const Real newr = -1.0 * lhs[row_idx]; 

      work_ptr->changeRange( work_ptr->rId( row_idx ), newl, newr );
      _assert( "lhs range value changed", work_ptr->lhs( row_idx ) == newl );
      _assert( "rhs range value changed", work_ptr->rhs( row_idx ) == newr );
   }

   work_ptr->solve();
   _assert_EQrel( "range changed", original_obj, work_ptr->objValue() );

   delete work_ptr; work_ptr = 0;

   // Third test: set both range values to solution*rowVector
   work_ptr = _prepare_Solver();
   work_ptr->solve();
   original_obj = work_ptr->objValue();
   
   Real* val = new Real[ work_ptr->nCols() ];
   Vector solution( work_ptr->nCols(), val );
   work_ptr->getPrimal( solution );
   
   for(int row_idx = 0; row_idx < work_ptr->nRows(); ++row_idx)
   {
      const SPxRowId& row_ID = work_ptr->rId( row_idx );
      Real row_prod = solution * (work_ptr->rowVector( row_idx ) );
      work_ptr->changeRange( row_ID, row_prod, row_prod );
   }
   work_ptr->solve();
   _assert_EQrel( "range changed", original_obj, work_ptr->objValue() );

   delete[] val; val = 0;
   delete work_ptr; work_ptr=0;
}


/**
   Changes row
*/
void ChangeExerciser::test_change_row()
{
   MSG( std::cout << "Testing changeRow()" << std::endl; );

   //
   // First test: Change rows to non-trivial multiples.
   //
   TestSolver* work_ptr = _prepare_Solver();
   work_ptr->solve();
   Real original_obj = work_ptr->objValue();    
   
   for (int row_idx = 0; row_idx < work_ptr->nRows(); ++row_idx) {
      Real change_coeff = 1.0;

      if ( row_idx % 2 == 0 )
         {
            change_coeff = 2.0;
         }
      else 
         {
            change_coeff = -1.0;
         }

      SVector work_vec( work_ptr->rowVector( row_idx ) );
            
      Real lhs = work_ptr->lhs( row_idx );
      Real rhs = work_ptr->rhs( row_idx );
            
      lhs *= change_coeff;
      rhs *= change_coeff;
      work_vec *= change_coeff;

      if ( change_coeff < 0.0 )
         {
            Real tmp = lhs;
            lhs = rhs;
            rhs = tmp;
         }
      assert( lhs <= rhs );

      const LPRow row( lhs, work_vec, rhs);
            
      // get ID
      const SPxRowId& row_ID = work_ptr->rId( row_idx );
      work_ptr->changeRow( row_ID, row );
            
      work_ptr->solve();
   }

   work_ptr->solve();
   _assert_EQrel( "rows changed", original_obj, work_ptr->objValue() );
   delete work_ptr; work_ptr = 0;

   //
   // Second test: Exchange rows pairwise.
   //
   work_ptr = _prepare_Solver();
   work_ptr->solve();
   original_obj = work_ptr->objValue();

   const int nPairs = work_ptr->nRows()/2;  
   for (int pair_idx = 0; pair_idx < nPairs; ++pair_idx) {
      const int first_idx = 2 * pair_idx; 
      const int second_idx = first_idx + 1;
      
      // Get rows.
      LPRow f_row( work_ptr->rowVector( first_idx ).size() );
      work_ptr->getRow( first_idx, f_row );
      
      LPRow s_row( work_ptr->rowVector( second_idx ).size() );
      work_ptr->getRow( second_idx, s_row);

      // Exchange rows and solve.
      work_ptr->changeRow( first_idx, s_row );
      work_ptr->changeRow( second_idx, f_row ); 
      work_ptr->solve();

      _assert_EQrel( "rows exchanged pairwise", original_obj, work_ptr->objValue() );
   }

   delete work_ptr; work_ptr = 0;
}


/**
   Changes col
*/
void ChangeExerciser::test_change_col()
{
   MSG( std::cout << "Testing changeCol()" << std::endl; );
   
   TestSolver* work_ptr = _prepare_Solver();
   work_ptr->solve();
   Real original_obj = work_ptr->objValue();    
   
   //
   // First test: Multiply every column, lhs, and rhs by a non-trivial number.
   //             Set lower and upper bounds for each value closer to optimal solution.
   // 
   const Real change_coeff = 2.0;
   Real* val = new Real[ work_ptr->nCols() ];
   Vector solution( work_ptr->nCols(), val );
   work_ptr->getPrimal( solution );

   Vector lhs = work_ptr->lhs(); lhs *= change_coeff;
   Vector rhs = work_ptr->rhs(); rhs *= change_coeff;
   work_ptr->changeRange( lhs, rhs );

   for (int col_idx = 0; col_idx < work_ptr->nCols(); ++col_idx) 
   {
      LPCol current_col( work_ptr->colVector( col_idx ).size() );
      work_ptr->getCol( col_idx, current_col );

      const Real new_lb = ( current_col.lower() == -infinity ) ? -infinity :
         std::max< Real >( current_col.lower(),
                           0.5 * current_col.lower() + 0.5 * solution[col_idx] );
      const Real new_ub = ( current_col.upper() == infinity ) ? infinity :
         std::min< Real >( current_col.upper(),
                           0.5 * current_col.upper() + 0.5 * solution[col_idx] );

      SVector new_col = current_col.colVector();
      new_col *= change_coeff;

      LPCol changed_col( current_col.obj(),
                         new_col, 
                         new_ub,
                         new_lb );

      work_ptr->changeCol( work_ptr->cId( col_idx ), changed_col );
   }

   work_ptr->solve();
   _assert_EQrel( "cols changed by changing bounds", original_obj, work_ptr->objValue() );

   delete val; val = 0;
   delete work_ptr; work_ptr = 0;

   //
   // Second test: Exchange columns pairwise.
   //
   work_ptr = _prepare_Solver();
   work_ptr->solve();
   original_obj = work_ptr->objValue();

   const int nPairs = work_ptr->nCols()/2;  
   for (int pair_idx = 0; pair_idx < nPairs; ++pair_idx)
   {
      const int first_idx = 2 * pair_idx; 
      const int second_idx = first_idx + 1;
      
      LPCol f_col( work_ptr->colVector( first_idx ).size() );
      work_ptr->getCol( first_idx, f_col );

      LPCol s_col( work_ptr->colVector( second_idx ).size() );
      work_ptr->getCol( second_idx, s_col );

      work_ptr->changeCol( first_idx, s_col );
      work_ptr->changeCol( second_idx, f_col );
      work_ptr->solve();

      _assert_EQrel( "cols changed by exchanging columns", original_obj, work_ptr->objValue() );
   }

   delete work_ptr; work_ptr = 0;
}


/**
   Changes single elements of the matrix.
*/
void ChangeExerciser::test_change_element()
{ 
   MSG( std::cout << "Testing changeElement()" << std::endl; );

   TestSolver* work_ptr = _prepare_Solver();
   work_ptr->solve();
   const Real original_obj = work_ptr->objValue();    
   
   // Test: Exchange each pair of rows elementwise
   int nPairs = work_ptr->nRows()/2;  
   for (int pair_idx = 0; pair_idx < nPairs; ++pair_idx)
   {
      int first_idx = 2 * pair_idx; 
      int second_idx = first_idx + 1;

      const SVector& first_row = work_ptr->rowVector( first_idx );
      const SVector& second_row = work_ptr->rowVector( second_idx );

      const SPxRowId& first_row_ID = work_ptr->rId( first_idx );
      const SPxRowId& second_row_ID = work_ptr->rId( second_idx );
      
      for (int col_idx=0; col_idx<work_ptr->nCols(); ++col_idx)
      {
         Real new_f = second_row[col_idx];
         Real new_s = first_row[col_idx];
         const SPxColId& col_ID = work_ptr->cId( col_idx );
         
         work_ptr->changeElement( first_row_ID, col_ID, new_f );
         work_ptr->changeElement( second_row_ID, col_ID, new_s );

         const SVector& new_first_row = work_ptr->rowVector( first_idx );
         const SVector& new_second_row = work_ptr->rowVector( second_idx );
         _assert( "elements changed",  new_first_row[ col_idx ] == new_f );
         _assert( "elements changed",  new_second_row[ col_idx ] == new_s );
      }

      // flip left/right hand sides of the rows
      Real new_first_lhs = work_ptr->lhs(second_idx);
      work_ptr->changeLhs( second_idx, work_ptr->lhs(first_idx) );
      work_ptr->changeLhs( first_idx, new_first_lhs );

      Real new_first_rhs = work_ptr->rhs(second_idx);
      work_ptr->changeRhs( second_idx, work_ptr->rhs(first_idx) );
      work_ptr->changeRhs( first_idx, new_first_rhs );
      
      work_ptr->solve();
      _assert_EQrel( "elements changed", original_obj, work_ptr->objValue() );
   }

   delete work_ptr; work_ptr = 0;
}


/**
   Changes sense
*/
void ChangeExerciser::test_change_sense()
{
   MSG( std::cout << "Testing changeSense()" << std::endl; );

   TestSolver* work_ptr = _prepare_Solver();
   work_ptr->solve();
   const Real original_obj = work_ptr->objValue();    
   
   // Flip sense.
   if ( work_ptr->sense() == SPxLP::MINIMIZE )
   {
      work_ptr->changeSense( SPxLP::MAXIMIZE );
   }
   else 
   {     
      assert( work_ptr->sense() == SPxLP::MAXIMIZE );
      work_ptr->changeSense( SPxLP::MINIMIZE );
   }
   
   // Invert objective function.
   Real* val = new Real[ work_ptr->nCols() ];
   Vector original( work_ptr->nCols(), val );

   work_ptr->getObj( original );
   original *= -1.0;
   work_ptr->changeObj( original );
   work_ptr->solve();

   _assert_EQrel( "change sense", original_obj, -work_ptr->objValue() );
   
   delete[] val; val = 0;
   delete work_ptr; work_ptr = 0;
}


//-----------------------------------------------------------------------------
//Emacs Local Variables:
//Emacs mode:c++
//Emacs c-basic-offset:3
//Emacs tab-width:8
//Emacs indent-tabs-mode:nil
//Emacs End:
//-----------------------------------------------------------------------------
