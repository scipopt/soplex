/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the class library                   */
/*       SoPlex --- the Sequential object-oriented simPlex.                  */
/*                                                                           */
/*  Copyright 1996-2022 Zuse Institute Berlin                                */
/*                                                                           */
/*  Licensed under the Apache License, Version 2.0 (the "License");          */
/*  you may not use this file except in compliance with the License.         */
/*  You may obtain a copy of the License at                                  */
/*                                                                           */
/*      http://www.apache.org/licenses/LICENSE-2.0                           */
/*                                                                           */
/*  Unless required by applicable law or agreed to in writing, software      */
/*  distributed under the License is distributed on an "AS IS" BASIS,        */
/*  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied. */
/*  See the License for the specific language governing permissions and      */
/*  limitations under the License.                                           */
/*                                                                           */
/*  You should have received a copy of the Apache-2.0 license                */
/*  along with SoPlex; see the file LICENSE. If not email to soplex@zib.de.  */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

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
#include "spxdantzigpr.h"
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
   ///@{
   static const SLUFactor<Real>::UpdateType update = SLUFactor<Real>::FOREST_TOMLIN;
   static const Real delta;
   static const Real timelimit;
   static const Real epsilon;
   static const Real epsilon_factor;
   static const Real epsilon_update;
   static const int verbose = SPxOut::ERROR;
   static const int precision = 12;
   ///@}

private:
   //------------------------------------
   /**@name Data */
   ///@{
   SLUFactor<Real> _solver;              ///< sparse LU factorization
   SPxSteepPR _pricer;             ///< steepest edge pricer
   SPxFastRT _ratiotester;         ///< Harris fast ratio tester
   ///@}

public:

   //------------------------------------
   /**@name Construction / destruction */
   ///@{
   /// Default constructor.
   explicit
   TestSolver( const SPxSolver::Type type_,
               const SPxSolver::Representation representation_ )
      : SPxSolver( type_,
                   representation_ )
   {
      setDelta( delta  );
      setTerminationTime( timelimit );

      _tolerances->setEpsilon( epsilon );
      _tolerances->setEpsilonFactorization( epsilon_factor );
      _tolerances->setEpsilonUpdate( epsilon_update );
      _tolerances->setVerbose( verbose );

      setPricer( &_pricer );
      setTester( &_ratiotester );
      setBasisSolver( &_solver );
      _solver.setUtype( update );
      // no starter, no simplifier, no scaler

      assert( isConsistent() );
   }
   ///@}
};

//
// Define static members of "TestSolver".
//
const SLUFactor<Real>::UpdateType TestSolver::update;
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

    All tests are designed such that the number of calls to solve() does not
    depend on the size of the LP. However, one can increase the "resolution"
    of the tests at the cost of a higher running time by increasing
    \p calls_solve.
 */
class ChangeExerciser
{
   /**
      Precision used in (relative) equality checks for differing but logically
      equivalent LP solutions.
   */
   static const Real epsilon_solution_equal;

   /// Number of calls to solve() in a loop where something was changed.
public:
   static int calls_solve;

public:

   //------------------------------------
   /**@name Construction / destruction */
   ///@{
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
   ///@}

public:

   //------------------------------------
   /**@name Test methods */
   ///@{
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
   ///@}

private:

   //------------------------------------
   /**@name Testing support */
   ///@{
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
   ///@}
};

//
// Define static members of "ChangeExerciser".
//
const Real ChangeExerciser::epsilon_solution_equal = 1e-6;
int ChangeExerciser::calls_solve = 10;


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
      "[-v] [-s Solves] <test_suite>\n\n"
      "<test_suite> is supposed to be file containing names of LPs in either MPS or LPF format\n\n"
      "-v        activate verbose mode\n"
      "-sSolves  give number of calls to solve() per iteration (the more the slower and the more exhaustive testing), default: 10\n"
   ;

   // Cope with command line.
   if ( argc < 2 || argc > 5 )
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
      case 's' :
         ++optidx;
         if ( optidx < argc )
         {
            ChangeExerciser::calls_solve = atoi( argv[ optidx ] );
            if ( ChangeExerciser::calls_solve <= 0 )
            {
               std::cout << "usage error: invalid parameter for -s" << std::endl;
               exit( 1 );
            }
         }
         else
         {
            std::cout << "usage error: -s needs parameter" << std::endl;
            exit( 1 );
         }
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

   Timer timer;
   timer.start();

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
   timer.stop();

   std::cout << "Total number of failed asserts: " << total_asserts_failed << std::endl;
   std::cout << "Total running time: " << timer.userTime() << std::endl;


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
   Real original_obj = work_ptr->objValue();

   // First test: Remove and readd single rows.
   int solve_interval = std::max< int >( 1, work_ptr->nRows() / calls_solve );

   for( int row_idx = 0; row_idx < work_ptr->nRows(); row_idx += solve_interval )
   {
      const SVector& row_vector = work_ptr->rowVector( row_idx );
      LPRow current_row( row_vector.size() );

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

   delete work_ptr; work_ptr = 0;

   // Second test: Remove and readd blocks of rows.
   work_ptr = _prepare_Solver();
   work_ptr->solve();
   original_obj = work_ptr->objValue();

   Array< LPRow > rows( solve_interval );
   Array< SPxRowId > row_IDs( work_ptr->nRows() );

   // First remember row IDs in order to avoid confusion after rows have been deleted.
   for( int row_idx = 0; row_idx < work_ptr->nRows(); ++row_idx )
   {
      row_IDs[ row_idx ] = work_ptr->rId( row_idx );
   }

   for ( int block = 0; block < calls_solve + 1; ++block )
   {
      int first_block_row_idx = block * solve_interval;
      int last_block_row_idx = std::min( work_ptr->nRows() - 1,
                                         (block+1) * solve_interval - 1 );

      // Save-and-delete loop.
      for( int row_idx = first_block_row_idx; row_idx <= last_block_row_idx; ++row_idx )
      {
         const SPxRowId& row_ID = row_IDs[ row_idx ];
         work_ptr->getRow( row_ID, rows[ row_idx - first_block_row_idx ]  );
         work_ptr->removeRow( row_ID );
      }
      work_ptr->solve();

      // Add loop (in reverse order).
      for( int row_idx = last_block_row_idx; row_idx >= first_block_row_idx; --row_idx )
      {
         work_ptr->addRow( rows[ row_idx - first_block_row_idx ]  );
      }

      work_ptr->solve();
      std::ostringstream description;
      description << "remove+add row block" << block;
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
   Real original_obj = work_ptr->objValue();

   // First test: Remove and readd single cols.
   int solve_interval = std::max< int >( 1, work_ptr->nCols() / calls_solve );

   for( int col_idx = 0; col_idx < work_ptr->nCols(); col_idx += solve_interval )
   {
      const SVector& col_vector = work_ptr->colVector( col_idx );
      LPCol current_col( col_vector.size() );

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

   delete work_ptr; work_ptr = 0;

   // Second test: Remove and readd blocks of cols.
   work_ptr = _prepare_Solver();
   work_ptr->solve();
   original_obj = work_ptr->objValue();

   Array< LPCol > cols( solve_interval );
   Array< SPxColId > col_IDs( work_ptr->nCols() );

   // First remember col IDs in order to avoid confusion after cols have been deleted.
   for( int col_idx = 0; col_idx < work_ptr->nCols(); ++col_idx )
   {
      col_IDs[ col_idx ] = work_ptr->cId( col_idx );
   }

   for ( int block = 0; block < calls_solve + 1; ++block )
   {
      int first_block_col_idx = block * solve_interval;
      int last_block_col_idx = std::min( work_ptr->nCols() - 1,
                                         (block+1) * solve_interval - 1 );

      // Save-and-delete loop.
      for( int col_idx = first_block_col_idx; col_idx <= last_block_col_idx; ++col_idx )
      {
         const SPxColId& col_ID = col_IDs[ col_idx ];
         work_ptr->getCol( col_ID, cols[ col_idx - first_block_col_idx ]  );
         work_ptr->removeCol( col_ID );
      }
      work_ptr->solve();

      // Add loop (in reverse order).
      for( int col_idx = last_block_col_idx; col_idx >= first_block_col_idx; --col_idx )
      {
         work_ptr->addCol( cols[ col_idx - first_block_col_idx ]  );
      }

      work_ptr->solve();
      std::ostringstream description;
      description << "remove+add col block" << block;
      _assert_EQrel( description.str(), original_obj, work_ptr->objValue() );
   }

   delete work_ptr; work_ptr = 0;

   // Third test: Remove and add single columns. This is separate from the first test to reproduce
   // a fail with netlib/scfxm1.mps.gz.
   work_ptr = _prepare_Solver();
   work_ptr->solve();
   original_obj = work_ptr->objValue();

   if ( work_ptr->nCols() >= 7 )
   {
      for( int col_idx = 0; col_idx < 8; ++col_idx )
      {
         const SVector& col_vector = work_ptr->colVector( col_idx );
         LPCol current_col( col_vector.size() );

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
   }

   delete work_ptr; work_ptr = 0;
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
   delete[] col_IDs; col_IDs = 0;

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
   work_ptr->getPrimalSol(solution);

   for (int col_idx = 0; col_idx != work_ptr->nCols(); ++col_idx)
   {
      if ( work_ptr->lower( col_idx ) < solution[col_idx] )
      {
         const Real new_lower = std::min< Real >( solution[col_idx], work_ptr->upper( col_idx) );

         const SPxColId& col_ID = work_ptr->cId( col_idx );
         work_ptr->changeLower( col_ID, new_lower );
         _assert( "lower value changed", EQ( work_ptr->lower( col_idx ), new_lower ) );
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
   work_ptr->getPrimalSol(solution);

   for(int col_idx = 0; col_idx != work_ptr->nCols(); ++col_idx)
   {
      if( work_ptr->upper( col_idx ) > solution[col_idx] )
      {
         const Real new_upper = std::max< Real >( solution[col_idx], work_ptr->lower( col_idx ) );

         const SPxColId& col_ID = work_ptr->cId( col_idx );
         work_ptr->changeUpper( col_ID, new_upper );
         _assert( "upper value changed", EQ( work_ptr->upper( col_idx ), new_upper ) );
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
   work_ptr->getPrimalSol(solution);

   for (int col_idx = 0; col_idx != work_ptr->nCols(); ++col_idx)
   {
      Real upper = work_ptr->upper(col_idx);
      Real lower = work_ptr->lower(col_idx);

      // slack is bigger than 0
      if(upper > solution[col_idx])
         upper = std::max< Real >( solution[col_idx], work_ptr->lower( col_idx ) );

      if(solution[col_idx] > lower)
         lower = std::min< Real >( solution[col_idx], work_ptr->upper( col_idx) );

      // get ID
      const SPxColId& col_ID = work_ptr->cId( col_idx );
      work_ptr->changeBounds( col_ID, lower, upper );

      _assert( "upper bound changed", EQ( work_ptr->upper( col_idx ), upper ) );
      _assert( "lower bound changed", EQ( work_ptr->lower( col_idx ), lower ) );
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
   work_ptr->getPrimalSol( solution );

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
   work_ptr->getPrimalSol( solution );

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
   Real* val = new Real[ work_ptr->nCols() ];

   Vector solution( work_ptr->nCols(), val );
   work_ptr->getPrimalSol( solution );

   Vector new_rhs = work_ptr->rhs();
   Vector new_lhs = work_ptr->lhs();

   for (int row_idx = 0; row_idx < work_ptr->nRows(); ++row_idx)
   {
      const Real row_prod = solution * (work_ptr->rowVector( row_idx ) );

      new_rhs[ row_idx ] = std::max< Real >( work_ptr->lhs( row_idx ), row_prod );
      new_lhs[ row_idx ] = std::min< Real >( work_ptr->rhs( row_idx ), row_prod );
   }
   work_ptr->changeRange( new_lhs, new_rhs );

   for (int row_idx = 0; row_idx < work_ptr->nRows(); ++row_idx)
   {
      _assert( "range changed by vector: lhs", EQ( new_lhs[ row_idx ], work_ptr->lhs( row_idx ) ) );
      _assert( "range changed by vector: rhs", EQ( new_rhs[ row_idx ], work_ptr->rhs( row_idx ) ) );
   }

   work_ptr->solve();
   _assert_EQrel( "range changed by vector: objective", original_obj, work_ptr->objValue() );

   delete[] val; val = 0;
   delete work_ptr; work_ptr = 0;

#if 0
   // Second test: change range by indices
   work_ptr = _prepare_Solver();
   work_ptr->solve();
   original_obj = work_ptr->objValue();

   Vector rhs = work_ptr->rhs();
   Vector lhs = work_ptr->lhs();

   for (int row_idx = 0; row_idx < work_ptr->nRows(); ++row_idx)
   {
      const Real newl = 0.5 * lhs[ row_idx ];
      const Real newr = 2.0 * rhs[ row_idx ];

      work_ptr->changeRange( work_ptr->rId( row_idx ), newl, newr );
      _assert( "lhs range value changed", work_ptr->lhs( row_idx ) == newl );
      _assert( "rhs range value changed", work_ptr->rhs( row_idx ) == newr );
   }
   work_ptr->solve();

   // reverse the change
   for (int row_idx=0; row_idx < work_ptr->nRows(); ++row_idx)
   {
      const Real newl = lhs[ row_idx ];
      const Real newr = rhs[ row_idx ];

      work_ptr->changeRange( work_ptr->rId( row_idx ), newl, newr );
      _assert( "lhs range value changed", work_ptr->lhs( row_idx ) == newl );
      _assert( "rhs range value changed", work_ptr->rhs( row_idx ) == newr );
   }
   work_ptr->solve();

   _assert_EQrel( "range changed and restored", original_obj, work_ptr->objValue() );

   delete work_ptr; work_ptr = 0;
#endif

   // Third test: set both range values to solution*rowVector
   work_ptr = _prepare_Solver();
   work_ptr->solve();
   original_obj = work_ptr->objValue();

   val = new Real[ work_ptr->nCols() ];
   Vector solution3( work_ptr->nCols(), val );
   work_ptr->getPrimalSol( solution3 );

   for(int row_idx = 0; row_idx < work_ptr->nRows(); ++row_idx)
   {
      const SPxRowId& row_ID = work_ptr->rId( row_idx );
      Real row_prod = solution3 * (work_ptr->rowVector( row_idx ) );
      work_ptr->changeRange( row_ID, row_prod, row_prod );
   }
   work_ptr->solve();
   _assert_EQrel( "range changed to optimal bounds", original_obj, work_ptr->objValue() );

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

   int solve_interval = std::max< int >( 1, work_ptr->nRows() / calls_solve );

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

      if ( row_idx % solve_interval == 0 )
      {
         work_ptr->solve();
         _assert_EQrel( "rows changed", original_obj, work_ptr->objValue() );
      }
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
   solve_interval = std::max< int >( 1, nPairs / calls_solve );

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

      if ( pair_idx % solve_interval == 0 )
      {
         work_ptr->solve();
         _assert_EQrel( "rows exchanged pairwise", original_obj, work_ptr->objValue() );
      }
   }

   work_ptr->solve();
   _assert_EQrel( "final rows exchanged pairwise", original_obj, work_ptr->objValue() );

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
   work_ptr->getPrimalSol( solution );

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

   delete[] val; val = 0;
   delete work_ptr; work_ptr = 0;

   //
   // Second test: Exchange columns pairwise.
   //
   work_ptr = _prepare_Solver();
   work_ptr->solve();
   original_obj = work_ptr->objValue();

   const int nPairs = work_ptr->nCols()/2;
   int solve_interval = std::max< int >( 1, nPairs / calls_solve );

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

      if ( pair_idx % solve_interval == 0 )
      {
         work_ptr->solve();
         _assert_EQrel( "cols changed by exchanging columns", original_obj, work_ptr->objValue() );
      }
   }

   work_ptr->solve();
   _assert_EQrel( "final cols changed by exchanging columns", original_obj, work_ptr->objValue() );

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
   int solve_interval = std::max< int >( 1, nPairs / calls_solve );

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

      if ( pair_idx % solve_interval == 0 )
      {
         work_ptr->solve();
         _assert_EQrel( "elements changed", original_obj, work_ptr->objValue() );
      }
   }

   work_ptr->solve();
   _assert_EQrel( "final elements changed", original_obj, work_ptr->objValue() );

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
