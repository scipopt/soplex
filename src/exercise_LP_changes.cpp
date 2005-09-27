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
#pragma ident "@(#) $Id: exercise_LP_changes.cpp,v 1.1 2005/09/27 13:28:08 bzfhille Exp $"

#include <assert.h>
#include <math.h>
#include <string.h>
#include <iostream>
#include <iomanip>
#include <fstream>

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

/** 
    A simple derived class from #SoPlex, used as object under test.
 */
class TestSolver : public SoPlex
{
   /**
      Parameters used for LP solving. They mainly correspond to the defaults of
      the main binary, except for the ones indicated by "// *".
      Part of the default settings are also realized in the ctor below.
   */
public:
   static const SLUFactor::UpdateType update = SLUFactor::FOREST_TOMLIN;
   static const Real delta = DEFAULT_BND_VIOL;
   static const Real timelimit = -1.0;
   static const Real epsilon = DEFAULT_EPS_ZERO;
   static const Real epsilon_factor = DEFAULT_EPS_FACTOR;
   static const Real epsilon_update = DEFAULT_EPS_UPDATE;
   static const int verbose = 1;
   static const int precision = 12;

private:
   SPxPricer* _pricer;
   SPxRatioTester* _ratiotester;
   SPxScaler* _prescaler;
   SPxScaler* _postscaler;
   SPxSimplifier* _simplifier;
   SPxStarter* _starter;

public:
   /// Default constructor.
   TestSolver( const SPxSolver::Type type_ = SPxSolver::LEAVE, 
               const SPxSolver::Representation representation_ = SPxSolver::COLUMN )
      : SoPlex( type_, representation_ )
      , _pricer( 0 )
      , _ratiotester( 0 )
      , _prescaler( 0 )
      , _postscaler( 0 )
      , _simplifier( 0 )
      , _starter( 0 )
   {
      setUtype( update );
      setDelta( delta  );
      setTerminationTime( timelimit );

      Param::setEpsilon( epsilon );
      Param::setEpsilonFactorization( epsilon_factor );
      Param::setEpsilonUpdate( epsilon_update );
      Param::setVerbose( verbose );

      // Use default sub-algorithms as in main soplex binary.
      _pricer = new SPxSteepPR;
      _ratiotester = new SPxFastRT;
      _prescaler = new SPxEquiliSC( representation_ == SPxSolver::COLUMN, true );
      _postscaler = new SPxGeometSC( representation_ == SPxSolver::COLUMN, 1 );
      _simplifier = 0; // * no simplifier
      _starter = 0;

      setPricer( _pricer );
      setTester( _ratiotester );
      setPreScaler( _prescaler );
      setPostScaler( _postscaler );
      setSimplifier( _simplifier );
      setStarter( _starter );

      assert( isConsistent() );
   }

   //------------------------------------------------------------------------
   /// virtual destructor
   virtual ~TestSolver()
   {
      delete _pricer;
      delete _ratiotester;
      delete _prescaler;
      delete _postscaler;
      delete _simplifier;
      delete _starter;
   }
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
   /// Default constructor.
   ChangeExerciser( const std::string& instance_name )
      : _asserts_failed( 0 )
      , _instance_name( instance_name )
   {}

   //@{ Test methods     
public:
   void test_add_delete_row();

   /// Reset test statistics.
   void reset()
   {
      _asserts_failed = 0;
   }

   /// Number of failed asserts.
   long asserts_failed() const { return _asserts_failed; }
   //@}

   //@{ Testing support 
private:
   void _assert_EQrel( const std::string& description, const Real a, const Real b )
   {
      if ( !EQrel( a, b, epsilon_solution_equal ) )
         {
            ++_asserts_failed;
            std::cout << "check '" << description << "' failed" << std::endl;
         }
   }

   TestSolver* _prepare_Solver() const 
   {
      TestSolver* work_ptr = new TestSolver;
      work_ptr->readFile( _instance_name.c_str(), 0, 0 );
      return work_ptr;
   }

   long _asserts_failed;
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

         std::cout << tester.asserts_failed() << " asserts failed.\n" << std::endl;
      }

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
            }
         else
            {
               const SPxRowId& row_ID = work_ptr->rId( row_idx );
               work_ptr->getRow( row_ID, current_row );               
            }
         
         work_ptr->removeRow( row_idx );
         work_ptr->solve();
         
         work_ptr->addRow( current_row );
         work_ptr->solve();

         char description[ 20 ];
         sprintf( description, "remove+add row %i", row_idx );
         _assert_EQrel( description, original_obj, work_ptr->objValue() );
      }

   delete work_ptr;
}


/**
   Deletes and re-adds sets of rows.
*/




//-----------------------------------------------------------------------------
//Emacs Local Variables:
//Emacs mode:c++
//Emacs c-basic-offset:3
//Emacs tab-width:8
//Emacs indent-tabs-mode:nil
//Emacs End:
//-----------------------------------------------------------------------------
