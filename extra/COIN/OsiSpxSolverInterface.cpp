//  LAST EDIT: Tue Nov 27 16:30:04 CET 2001 by Laszlo Ladanyi
//-----------------------------------------------------------------------------
// name:     OSI Interface for SOPLEX
// author:   Tobias Pfender
//           Konrad-Zuse-Zentrum Berlin (Germany)
//           email: pfender@zib.de
// date:     01/16/2002
//-----------------------------------------------------------------------------
// Copyright (C) 2002, Tobias Pfender, International Business Machines
// Corporation and others.  All Rights Reserved.
#ifdef COIN_USE_SPX
#if defined(_MSC_VER)
// Turn off compiler warning about long names
#  pragma warning(disable:4786)
#endif

#include <iostream>
#include <cassert>
#include <string>
#include <numeric>

#include "CoinError.hpp"

#include "OsiSpxSolverInterface.hpp"
#include "OsiRowCut.hpp"
#include "OsiColCut.hpp"
#include "OsiPackedMatrix.hpp"
#include "OsiWarmStartBasis.hpp"

//#############################################################################
// A couple of helper functions
//#############################################################################

inline void freeCacheDouble( double*& ptr )
{
  if( ptr != NULL )
    {
      delete [] ptr;
      ptr = NULL;
    }
}

inline void freeCacheChar( char*& ptr )
{
  if( ptr != NULL )
    {
      delete [] ptr;
      ptr = NULL;
    }
}

inline void freeCacheMatrix( OsiPackedMatrix*& ptr )
{
  if( ptr != NULL )
    {
      delete ptr;
      ptr = NULL;
    }
}

inline void throwSPXerror( std::string error, std::string osimethod )
{
  std::cout << "ERROR: " << error << " (" << osimethod << 
    " in OsiSpxSolverInterface)" << std::endl;
  throw CoinError( error.c_str(), osimethod.c_str(), "OsiSpxSolverInterface" );
}

//#############################################################################
// Solve methods
//#############################################################################

void OsiSpxSolverInterface::initialSolve()
{
  if( spxsolver_.type() != soplex::SoPlex::ENTER )
    spxsolver_.setType( soplex::SoPlex::ENTER );   // switch to primal simplex
  if( spxsolver_.rep() != soplex::SoPlex::COLUMN )
    spxsolver_.setRep( soplex::SoPlex::COLUMN );
  spxsolver_.solve();
}
//-----------------------------------------------------------------------------
void OsiSpxSolverInterface::resolve()
{
  if( spxsolver_.type() != soplex::SoPlex::LEAVE )
    spxsolver_.setType( soplex::SoPlex::LEAVE );   // switch to dual simplex
  if( spxsolver_.rep() != soplex::SoPlex::COLUMN )
    spxsolver_.setRep( soplex::SoPlex::COLUMN );
  spxsolver_.solve();
}
//-----------------------------------------------------------------------------
void OsiSpxSolverInterface::branchAndBound()
{
  throwSPXerror( "SOPLEX does not provide an internal branch and bound procedure",
		 "branchAndBound" );
}

//#############################################################################
// Parameter related methods
//#############################################################################

bool
OsiSpxSolverInterface::setIntParam(OsiIntParam key, int value)
{
  bool retval = false;
  switch (key)
    {
    case OsiMaxNumIteration:
      spxsolver_.setTerminationIter( value );
      retval = true;
      break;
    case OsiMaxNumIterationHotStart:
      if( value >= 0 )
	{
	  hotStartMaxIteration_ = value;
	  retval = true;
	}
      else
	retval = false;
      break;
    case OsiLastIntParam:
      retval = false;
      break;
    }
  return retval;
}

//-----------------------------------------------------------------------------

bool
OsiSpxSolverInterface::setDblParam(OsiDblParam key, double value)
{
  bool retval = false;
  switch (key)
    {
    case OsiDualObjectiveLimit:
      // SOPLEX doesn't support different termination values for primal and dual simplex
      // at the moment, setting a termination value is not supported in SOPLEX
      // spxsolver_.setTerminationValue( value );
      // retval = true;
      retval = false;
      break;
    case OsiPrimalObjectiveLimit:
      // SOPLEX doesn't support different termination values for primal and dual simplex
      // at the moment, setting a termination value is not supported in SOPLEX
      // spxsolver_.setTerminationValue( value );
      // retval = true;
      retval = false;
      break;
    case OsiDualTolerance:
      // ??? Is delta the correct SOPLEX equivalent of dual tolerance?
      // SOPLEX doesn't support different deltas for primal and dual simplex
      spxsolver_.setDelta( value );
      retval = true;
      break;
    case OsiPrimalTolerance:
      // ??? Is delta the correct SOPLEX equivalent of primal tolerance?
      // SOPLEX doesn't support different deltas for primal and dual simplex
      spxsolver_.setDelta( value );
      retval = true;
      break;
    case OsiObjOffset:
      retval = OsiSolverInterface::setDblParam(key,value);
      break;
    case OsiLastDblParam:
      retval = false;
      break;
    }
  return retval;
}

//-----------------------------------------------------------------------------

bool
OsiSpxSolverInterface::getIntParam(OsiIntParam key, int& value) const
{
  bool retval = false;
  switch (key)
    {
    case OsiMaxNumIteration:
      value = spxsolver_.terminationIter();
      retval = true;
      break;
    case OsiMaxNumIterationHotStart:
      value = hotStartMaxIteration_;
      retval = true;
      break;
    case OsiLastIntParam:
      retval = false;
      break;
    }
  return retval;
}

//-----------------------------------------------------------------------------

bool
OsiSpxSolverInterface::getDblParam(OsiDblParam key, double& value) const
{
  bool retval = false;
  switch (key) 
    {
    case OsiDualObjectiveLimit:
      // at the moment, setting a termination value is not supported in SOPLEX
      // value = spxsolver_.terminationValue();
      // retval = true;
      retval = false;
      break;
    case OsiPrimalObjectiveLimit:
      // at the moment, setting a termination value is not supported in SOPLEX
      // value = spxsolver_.terminationValue();
      // retval = true;
      retval = false;
      break;
    case OsiDualTolerance:
      value = spxsolver_.delta();
      retval = true;
      break;
    case OsiPrimalTolerance:
      value = spxsolver_.delta();
      retval = true;
      break;
    case OsiObjOffset:
      retval = OsiSolverInterface::getDblParam(key, value);
      break;
    case OsiLastDblParam:
      retval = false;
      break;
    }
  return retval;
}

//#############################################################################
// Methods returning info on how the solution process terminated
//#############################################################################

bool OsiSpxSolverInterface::isAbandoned() const
{
  int stat = spxsolver_.status();

  return ( stat == soplex::SoPlex::SINGULAR ||
	   stat == soplex::SoPlex::ERROR    );
}

bool OsiSpxSolverInterface::isProvenOptimal() const
{
  int stat = spxsolver_.status();

  return ( stat == soplex::SoPlex::OPTIMAL );
}

bool OsiSpxSolverInterface::isProvenPrimalInfeasible() const
{
  int stat = spxsolver_.status();

  return ( stat == soplex::SoPlex::INFEASIBLE );
}

bool OsiSpxSolverInterface::isProvenDualInfeasible() const
{
  int stat = spxsolver_.status();

  return ( stat == soplex::SoPlex::UNBOUNDED );
}

bool OsiSpxSolverInterface::isPrimalObjectiveLimitReached() const
{
  return ( spxsolver_.status() == soplex::SoPlex::ABORT_VALUE );
}

bool OsiSpxSolverInterface::isDualObjectiveLimitReached() const
{
  return ( spxsolver_.status() == soplex::SoPlex::ABORT_VALUE );
}

bool OsiSpxSolverInterface::isIterationLimitReached() const
{
  return ( spxsolver_.status() == soplex::SoPlex::ABORT_ITER );
}

//#############################################################################
// WarmStart related methods
//#############################################################################

OsiWarmStart* OsiSpxSolverInterface::getWarmStart() const
{
  OsiWarmStartBasis* ws = NULL;
  int numcols = getNumCols();
  int numrows = getNumRows();
  int i;

  ws = new OsiWarmStartBasis();
  ws->setSize( numcols, numrows );

  for( i = 0; i < numrows; ++i )
    {
      switch( spxsolver_.getBasisRowStatus( i ) )
	{
	case soplex::SoPlex::BASIC:
	  ws->setArtifStatus( i, OsiWarmStartBasis::basic );
	  break;	  
	case soplex::SoPlex::FIXED:
	case soplex::SoPlex::ON_LOWER:
	  ws->setArtifStatus( i, OsiWarmStartBasis::atLowerBound );
	  break;
	case soplex::SoPlex::ON_UPPER:
	  ws->setArtifStatus( i, OsiWarmStartBasis::atUpperBound );
	  break;
	case soplex::SoPlex::ZERO:
	  ws->setArtifStatus( i, OsiWarmStartBasis::isFree );
	  break;
	default:
	  throwSPXerror( "invalid row status", "getWarmStart" );
	  break;
	}
    }

  for( i = 0; i < numcols; ++i )
    {
      switch( spxsolver_.getBasisColStatus( i ) )
	{
	case soplex::SoPlex::BASIC:
	  ws->setStructStatus( i, OsiWarmStartBasis::basic );
	  break;	  
	case soplex::SoPlex::FIXED:
	case soplex::SoPlex::ON_LOWER:
	  ws->setStructStatus( i, OsiWarmStartBasis::atLowerBound );
	  break;
	case soplex::SoPlex::ON_UPPER:
	  ws->setStructStatus( i, OsiWarmStartBasis::atUpperBound );
	  break;
	case soplex::SoPlex::ZERO:
	  ws->setStructStatus( i, OsiWarmStartBasis::isFree );
	  break;
	default:
	  throwSPXerror( "invalid column status", "getWarmStart" );
	  break;
	}
    }

  return ws;
}

//-----------------------------------------------------------------------------

bool OsiSpxSolverInterface::setWarmStart(const OsiWarmStart* warmstart)
{
  const OsiWarmStartBasis* ws = dynamic_cast<const OsiWarmStartBasis*>(warmstart);
  int numcols, numrows, i;
  soplex::SoPlex::VarStatus *cstat, *rstat;
  bool retval = false;

  if( !ws )
    return false;

  numcols = ws->getNumStructural();
  numrows = ws->getNumArtificial();
  
  if( numcols != getNumCols() || numrows != getNumRows() )
    return false;

  cstat = new soplex::SoPlex::VarStatus[numcols];
  rstat = new soplex::SoPlex::VarStatus[numrows];

  for( i = 0; i < numrows; ++i )
    {
      switch( ws->getArtifStatus( i ) )
	{
	case OsiWarmStartBasis::basic:
	  rstat[i] = soplex::SoPlex::BASIC;
	  break;
	case OsiWarmStartBasis::atLowerBound:
	  rstat[i] = soplex::SoPlex::ON_LOWER;
	  break;
	case OsiWarmStartBasis::atUpperBound:
	  rstat[i] = soplex::SoPlex::ON_UPPER;
	  break;
	case OsiWarmStartBasis::isFree:
	  rstat[i] = soplex::SoPlex::ZERO;
	  break;
	default:  // unknown row status
	  retval = false;
	  goto TERMINATE;
	}
    }
  for( i = 0; i < numcols; ++i )
    {
      switch( ws->getStructStatus( i ) )
	{
	case OsiWarmStartBasis::basic:
	  cstat[i] = soplex::SoPlex::BASIC;
	  break;
	case OsiWarmStartBasis::atLowerBound:
	  cstat[i] = soplex::SoPlex::ON_LOWER;
	  break;
	case OsiWarmStartBasis::atUpperBound:
	  cstat[i] = soplex::SoPlex::ON_UPPER;
	  break;
	case OsiWarmStartBasis::isFree:
	  cstat[i] = soplex::SoPlex::ZERO;
	  break;
	default:  // unknown column status
	  retval = false;
	  goto TERMINATE;
	}
    }

  spxsolver_.setBasis( rstat, cstat );
  retval = true;

 TERMINATE:
  delete[] cstat;
  delete[] rstat;
  return retval;
}

//#############################################################################
// Hotstart related methods (primarily used in strong branching)
//#############################################################################

void OsiSpxSolverInterface::markHotStart()
{
  int numcols, numrows;

  numcols = getNumCols();
  numrows = getNumRows();
  if( numcols > hotStartCStatSize_ )
    {
      delete[] hotStartCStat_;
      hotStartCStatSize_ = static_cast<int>( 1.2 * static_cast<double>( numcols ) ); // get some extra space for future hot starts
      hotStartCStat_ = new soplex::SoPlex::VarStatus[hotStartCStatSize_];
    }
  if( numrows > hotStartRStatSize_ )
    {
      delete[] hotStartRStat_;
      hotStartRStatSize_ = static_cast<int>( 1.2 * static_cast<double>( numrows ) ); // get some extra space for future hot starts
      hotStartRStat_ = new soplex::SoPlex::VarStatus[hotStartRStatSize_];
    }
  spxsolver_.getBasis( hotStartRStat_, hotStartCStat_ );
}

void OsiSpxSolverInterface::solveFromHotStart()
{
  int maxiter;

  assert( getNumCols() <= hotStartCStatSize_ );
  assert( getNumRows() <= hotStartRStatSize_ );

  spxsolver_.setBasis( hotStartRStat_, hotStartCStat_ );

  maxiter = spxsolver_.terminationIter();
  spxsolver_.setTerminationIter( hotStartMaxIteration_ );

  resolve();

  spxsolver_.setTerminationIter( maxiter );
}

void OsiSpxSolverInterface::unmarkHotStart()
{
  // be lazy with deallocating memory and do nothing here, deallocate memory in the destructor
}

//#############################################################################
// Problem information methods (original data)
//#############################################################################

//------------------------------------------------------------------
// Get number of rows, columns, elements, ...
//------------------------------------------------------------------
int OsiSpxSolverInterface::getNumCols() const
{
  return spxsolver_.nCols();
}
int OsiSpxSolverInterface::getNumRows() const
{
  return spxsolver_.nRows();
}
int OsiSpxSolverInterface::getNumElements() const
{
  return spxsolver_.nNzos();
}

//------------------------------------------------------------------
// Get pointer to rim vectors
//------------------------------------------------------------------  

const double * OsiSpxSolverInterface::getColLower() const
{
  return spxsolver_.lower().get_const_ptr();
}
//------------------------------------------------------------------
const double * OsiSpxSolverInterface::getColUpper() const
{
  return spxsolver_.upper().get_const_ptr();
}
//------------------------------------------------------------------
const char * OsiSpxSolverInterface::getRowSense() const
{
  if ( rowsense_ == NULL )
    {
      // rowsense is determined with rhs, so invoke rhs
      getRightHandSide();
      assert( rowsense_ != NULL || getNumRows() == 0 );
    }
  return rowsense_;
}
//------------------------------------------------------------------
const double * OsiSpxSolverInterface::getRightHandSide() const
{
  if ( rhs_ == NULL )
    {
      int nrows = getNumRows();
      if( nrows > 0 ) 
	{
	  int    row;

	  assert( rowrange_ == NULL );
	  assert( rowsense_ == NULL );

	  rhs_      = new double[nrows];
	  rowrange_ = new double[nrows];
	  rowsense_ = new char[nrows];
	  
	  for( row = 0; row < nrows; ++row )
	    convertBoundToSense( spxsolver_.lhs( row ), spxsolver_.rhs( row ),
				 rowsense_[row], rhs_[row], rowrange_[row] );
	}
    }
  return rhs_;
}
//------------------------------------------------------------------
const double * OsiSpxSolverInterface::getRowRange() const
{
  if ( rowrange_==NULL ) 
    {
      // rowrange is determined with rhs, so invoke rhs
      getRightHandSide();
      assert( rowrange_ != NULL || getNumRows() == 0 );
    }
  return rowrange_;
}
//------------------------------------------------------------------
const double * OsiSpxSolverInterface::getRowLower() const
{
  return spxsolver_.lhs().get_const_ptr();
}
//------------------------------------------------------------------
const double * OsiSpxSolverInterface::getRowUpper() const
{  
  return spxsolver_.rhs().get_const_ptr();
}
//------------------------------------------------------------------
const double * OsiSpxSolverInterface::getObjCoefficients() const
{
  if( obj_ == NULL )
    {
      obj_ = new soplex::DVector( getNumCols() );
      spxsolver_.getObj( *obj_ );
    }
  return obj_->get_const_ptr();
}
//------------------------------------------------------------------
double OsiSpxSolverInterface::getObjSense() const
{
  switch( spxsolver_.spxSense() )
    {
    case soplex::SPxLP::MINIMIZE:
      return +1.0;
    case soplex::SPxLP::MAXIMIZE:
      return -1.0;
    default:
      throwSPXerror( "invalid optimization sense", "getObjSense" );
      return 0.0;
    }
}

//------------------------------------------------------------------
// Return information on integrality
//------------------------------------------------------------------

bool OsiSpxSolverInterface::isContinuous( int colNumber ) const
{
  return( spxintvars_.number( colNumber ) < 0 );
}

//------------------------------------------------------------------
// Row and column copies of the matrix ...
//------------------------------------------------------------------

const OsiPackedMatrix * OsiSpxSolverInterface::getMatrixByRow() const
{
  if( matrixByRow_ == NULL )
    {
      int    nrows     = getNumRows();
      int    ncols     = getNumCols();
      int    nelems    = getNumElements();
      double *elements = new double [nelems];
      int    *indices  = new int    [nelems];
      int    *starts   = new int    [nrows+1];
      int    *len      = new int    [nrows];
      int    row, i, elem;

      elem = 0;
      for( row = 0; row < nrows; ++row )
	{
	  const soplex::SVector& rowvec = spxsolver_.rowVector( row );
	  starts[row] = elem;
	  len   [row] = rowvec.size();
	  for( i = 0; i < len[row]; ++i, ++elem )
	    {
	      assert( elem < nelems );
	      elements[elem] = rowvec.value( i );
	      indices [elem] = rowvec.index( i );
	    }
	}
      starts[nrows] = elem;
      assert( elem == nelems );

      matrixByRow_ = new OsiPackedMatrix();
      matrixByRow_->assignMatrix( false /* not column ordered */,
				  ncols, nrows, nelems,
				  elements, indices, starts, len );      
    }
  return matrixByRow_;
} 

//------------------------------------------------------------------

const OsiPackedMatrix * OsiSpxSolverInterface::getMatrixByCol() const
{
  if( matrixByCol_ == NULL )
    {
      int    nrows     = getNumRows();
      int    ncols     = getNumCols();
      int    nelems    = getNumElements();
      double *elements = new double [nelems];
      int    *indices  = new int    [nelems];
      int    *starts   = new int    [ncols+1];
      int    *len      = new int    [ncols];
      int    col, i, elem;

      elem = 0;
      for( col = 0; col < ncols; ++col )
	{
	  const soplex::SVector& colvec = spxsolver_.colVector( col );
	  starts[col] = elem;
	  len   [col] = colvec.size();
	  for( i = 0; i < len[col]; ++i, ++elem )
	    {
	      assert( elem < nelems );
	      elements[elem] = colvec.value( i );
	      indices [elem] = colvec.index( i );
	    }
	}
      starts[ncols] = elem;
      assert( elem == nelems );

      matrixByCol_ = new OsiPackedMatrix();
      matrixByCol_->assignMatrix( true /* column ordered */,
				  nrows, ncols, nelems,
				  elements, indices, starts, len );      
    }
  return matrixByCol_;
} 

//------------------------------------------------------------------
// Get solver's value for infinity
//------------------------------------------------------------------
double OsiSpxSolverInterface::getInfinity() const
{
  return soplex::infinity;
}

//#############################################################################
// Problem information methods (results)
//#############################################################################

// *FIXME*: what should be done if a certain vector doesn't exist???

const double * OsiSpxSolverInterface::getColSolution() const
{
  if( colsol_ == NULL )
    {
      int ncols = getNumCols();
      if( ncols > 0 )
	{
	  colsol_ = new soplex::DVector( ncols );
	  if( isProvenOptimal() )
	    spxsolver_.getPrimal( *colsol_ );
	  else
	    colsol_->clear();
	}
      else
	return NULL;
    }
  return colsol_->get_const_ptr();
}
//------------------------------------------------------------------
const double * OsiSpxSolverInterface::getRowPrice() const
{
  if( rowsol_ == NULL )
    {
      int nrows = getNumRows();
      if( nrows > 0 )
	{
	  rowsol_ = new soplex::DVector( nrows );
	  if( isProvenOptimal() )
	    spxsolver_.getDual( *rowsol_ );
	  else
	    rowsol_->clear();
	}
      else
	return NULL;
    }
  return rowsol_->get_const_ptr();
}
//------------------------------------------------------------------
const double * OsiSpxSolverInterface::getReducedCost() const
{
  if( redcost_ == NULL )
    {
      int ncols = getNumCols();
      if( ncols > 0 )
	{
	  redcost_ = new soplex::DVector( ncols );
	  if( isProvenOptimal() )
	    spxsolver_.getRdCost( *redcost_ );
	  else
	    redcost_->clear();
	}
      else
	return NULL;
    }
  return redcost_->get_const_ptr();
}
//------------------------------------------------------------------
const double * OsiSpxSolverInterface::getRowActivity() const
{
  if( rowact_ == NULL )
    {
      int nrows = getNumRows();
      if( nrows > 0 )
	{
	  rowact_ = new soplex::DVector( nrows );
	  if( isProvenOptimal() )
	    spxsolver_.getSlacks( *rowact_ );
	  else
	    rowact_->clear();
	}
      else
	return NULL;
    }
  return rowact_->get_const_ptr();
}
//------------------------------------------------------------------
double OsiSpxSolverInterface::getObjValue() const
{
  double objval;

  switch( spxsolver_.status() )
    {
    case soplex::SoPlex::OPTIMAL:
    case soplex::SoPlex::UNBOUNDED:
    case soplex::SoPlex::INFEASIBLE:
      objval = spxsolver_.value();
      break;
    default:
      objval = 0.0;
      break;
    }

  // Adjust objective function value by constant term in objective function
  double objOffset;
  getDblParam(OsiObjOffset,objOffset);
  objval = objval - objOffset;

  return objval;
}
//------------------------------------------------------------------
int OsiSpxSolverInterface::getIterationCount() const
{
  return spxsolver_.iterations();
}
//------------------------------------------------------------------
std::vector<double*> OsiSpxSolverInterface::getDualRays(int maxNumRays) const
{
  // *FIXME* : must write the method
  throw CoinError("method is not yet written", "getDualRays",
		  "OsiSpxSolverInterface");
  return std::vector<double*>();
}
//------------------------------------------------------------------
std::vector<double*> OsiSpxSolverInterface::getPrimalRays(int maxNumRays) const
{
  // *FIXME* : must write the method
  throw CoinError("method is not yet written", "getPrimalRays",
		  "OsiSpxSolverInterface");
  return std::vector<double*>();
}

//#############################################################################
// Problem modifying methods (rim vectors)
//#############################################################################

void OsiSpxSolverInterface::setObjCoeff( int elementIndex, double elementValue )
{
  spxsolver_.changeObj( elementIndex, elementValue );
  freeCachedData( OsiSpxSolverInterface::FREECACHED_COLUMN );
}

void OsiSpxSolverInterface::setColLower(int elementIndex, double elementValue)
{
  spxsolver_.changeLower( elementIndex, elementValue );
  freeCachedData( OsiSpxSolverInterface::FREECACHED_COLUMN );
}
//-----------------------------------------------------------------------------
void OsiSpxSolverInterface::setColUpper(int elementIndex, double elementValue)
{  
  spxsolver_.changeUpper( elementIndex, elementValue );
  freeCachedData( OsiSpxSolverInterface::FREECACHED_COLUMN );
} 
//-----------------------------------------------------------------------------
void OsiSpxSolverInterface::setColBounds( int elementIndex, double lower, double upper )
{
  spxsolver_.changeBounds( elementIndex, lower, upper );
  freeCachedData( OsiSpxSolverInterface::FREECACHED_COLUMN );
}
//-----------------------------------------------------------------------------
void
OsiSpxSolverInterface::setRowLower( int i, double elementValue )
{
  spxsolver_.changeLhs( i, elementValue );
  freeCachedData( OsiSpxSolverInterface::FREECACHED_ROW );
}
//-----------------------------------------------------------------------------
void
OsiSpxSolverInterface::setRowUpper( int i, double elementValue )
{
  spxsolver_.changeRhs( i, elementValue );
  freeCachedData( OsiSpxSolverInterface::FREECACHED_ROW );
}
//-----------------------------------------------------------------------------
void
OsiSpxSolverInterface::setRowBounds( int elementIndex, double lower, double upper )
{
  spxsolver_.changeRange( elementIndex, lower, upper );
  freeCachedData( OsiSpxSolverInterface::FREECACHED_ROW );
}
//-----------------------------------------------------------------------------
void
OsiSpxSolverInterface::setRowType(int i, char sense, double rightHandSide,
				  double range)
{
  double lower, upper;

  convertSenseToBound( sense, rightHandSide, range, lower, upper );
  setRowBounds( i, lower, upper );
}
//#############################################################################
void
OsiSpxSolverInterface::setContinuous(int index)
{
  int pos = spxintvars_.number( index );
  if( pos >= 0 )
    {
      spxintvars_.remove( pos );
      freeCachedData( OsiSpxSolverInterface::FREECACHED_COLUMN );
    }
}
//-----------------------------------------------------------------------------
void
OsiSpxSolverInterface::setInteger(int index)
{
  int pos = spxintvars_.number( index );
  if( pos < 0 )
    {
      spxintvars_.addIdx( index );
      freeCachedData( OsiSpxSolverInterface::FREECACHED_COLUMN );
    }
}
//#############################################################################

void OsiSpxSolverInterface::setObjSense(double s) 
{
  if( s != getObjSense() )
    {
      if( s == +1.0 )
	spxsolver_.changeSense( soplex::SPxLP::MINIMIZE );
      else
	spxsolver_.changeSense( soplex::SPxLP::MAXIMIZE );
      freeCachedData( OsiSpxSolverInterface::FREECACHED_RESULTS );
    }
}
 
//-----------------------------------------------------------------------------

void OsiSpxSolverInterface::setColSolution(const double * cs) 
{
  int col;
  int ncols = getNumCols();

  if( colsol_ != NULL )
    delete colsol_;

  if( ncols > 0 && cs != NULL )
    {
      colsol_ = new soplex::DVector( ncols );
      for( col = 0; col < ncols; ++col )
	(*colsol_)[col] = cs[col];
    }
  else
    colsol_ = NULL;
}

//-----------------------------------------------------------------------------

void OsiSpxSolverInterface::setRowPrice(const double * rs) 
{
  int row;
  int nrows = getNumRows();

  if( rowsol_ != NULL )
    delete rowsol_;

  if( nrows > 0 && rs != NULL )
    {
      rowsol_ = new soplex::DVector( nrows );
      for( row = 0; row < nrows; ++row )
	(*rowsol_)[row] = rs[row];
    }
  else
    rowsol_ = NULL;
}

//#############################################################################
// Problem modifying methods (matrix)
//#############################################################################
void 
OsiSpxSolverInterface::addCol(const OsiPackedVectorBase& vec,
			      const double collb, const double colub,   
			      const double obj)
{
  soplex::DSVector colvec;

  colvec.add( vec.getNumElements(), vec.getIndices(), vec.getElements() );
  spxsolver_.addCol( soplex::LPCol( obj, colvec, colub, collb ) );
  freeCachedData( OsiSpxSolverInterface::KEEPCACHED_ROW );
}
//-----------------------------------------------------------------------------
void 
OsiSpxSolverInterface::deleteCols(const int num, const int * columnIndices)
{
  spxsolver_.removeCols( const_cast<int*>(columnIndices), num );
  freeCachedData( OsiSpxSolverInterface::KEEPCACHED_ROW );
}
//-----------------------------------------------------------------------------
void 
OsiSpxSolverInterface::addRow(const OsiPackedVectorBase& vec,
			      const double rowlb, const double rowub)
{
  soplex::DSVector rowvec;

  rowvec.add( vec.getNumElements(), vec.getIndices(), vec.getElements() );
  spxsolver_.addRow( soplex::LPRow( rowlb, rowvec, rowub ) );
  freeCachedData( OsiSpxSolverInterface::KEEPCACHED_COLUMN );
}
//-----------------------------------------------------------------------------
void 
OsiSpxSolverInterface::addRow(const OsiPackedVectorBase& vec,
			      const char rowsen, const double rowrhs,   
			      const double rowrng)
{
  double rowlb, rowub;

  convertSenseToBound( rowsen, rowrhs, rowrng, rowlb, rowub );
  addRow( vec, rowlb, rowub );
}
//-----------------------------------------------------------------------------
void 
OsiSpxSolverInterface::deleteRows(const int num, const int * rowIndices)
{
  spxsolver_.removeRows( const_cast<int*>(rowIndices), num );
  freeCachedData( OsiSpxSolverInterface::KEEPCACHED_COLUMN );
}

//#############################################################################
// Methods to input a problem
//#############################################################################

void
OsiSpxSolverInterface::loadProblem( const OsiPackedMatrix& matrix,
				    const double* collb, const double* colub,
				    const double* obj,
				    const double* rowlb, const double* rowub )
{
  int          ncols   = matrix.getNumCols();
  int          nrows   = matrix.getNumRows();
  const int    *length = matrix.getVectorLengths();
  const int    *start  = matrix.getVectorStarts();
  const double *elem   = matrix.getElements();
  const int    *index  = matrix.getIndices();
  double       *thecollb, *thecolub, *theobj, *therowlb, *therowub;
  
  // create defaults if parameter is NULL
  if( collb == NULL )
    {
      thecollb = new double[ncols];
      CoinFillN( thecollb, ncols, 0.0 );
    }
  else
    thecollb = const_cast<double*>(collb);
  if( colub == NULL )
    {
      thecolub = new double[ncols];
      CoinFillN( thecolub, ncols, getInfinity() );
    }
  else
    thecolub = const_cast<double*>(colub);
  if( obj == NULL )
    {
      theobj = new double[ncols];
      CoinFillN( theobj, ncols, 0.0 );
    }
  else
    theobj = const_cast<double*>(obj);
  if( rowlb == NULL )
    {
      therowlb = new double[nrows];
      CoinFillN( therowlb, nrows, -getInfinity() );
    }
  else
    therowlb = const_cast<double*>(rowlb);
  if( rowub == NULL )
    {
      therowub = new double[nrows];
      CoinFillN( therowub, nrows, +getInfinity() );
    }
  else
    therowub = const_cast<double*>(rowub);

  // copy problem into spxsolver_
  spxsolver_.clear();
  spxintvars_.clear();
  freeCachedData( OsiSpxSolverInterface::KEEPCACHED_NONE );

  if( matrix.isColOrdered() )
    {
      int row, col, pos;
      soplex::LPRowSet rowset( nrows, 0 );
      soplex::DSVector rowvec;
      soplex::LPColSet colset( ncols, matrix.getNumElements() );
      soplex::DSVector colvec;

      /* insert empty rows */
      rowvec.clear();
      for( row = 0; row < nrows; ++row )
         rowset.add( therowlb[row], rowvec, therowub[row] );
      spxsolver_.addRows( rowset );

      /* create columns */
      for( col = 0; col < ncols; ++col )
	{
	  pos = start[col];
	  colvec.clear();
	  colvec.add( length[col], &(index[pos]), &(elem[pos]) );
	  colset.add( theobj[col], thecollb[col], colvec, thecolub[col] );
	}
      
      spxsolver_.addCols( colset );
      // spxsolver_.changeRange( soplex::Vector( nrows, therowlb ), soplex::Vector( nrows, therowub ) );
    }
  else
    {
      int row, col, pos;
      soplex::LPRowSet rowset( nrows, matrix.getNumElements() );
      soplex::DSVector rowvec;
      soplex::LPColSet colset( ncols, 0 );
      soplex::DSVector colvec;

      /* insert empty columns */
      colvec.clear();
      for( col = 0; col < ncols; ++col )
         colset.add( theobj[col], thecollb[col], colvec, thecolub[col] );
      spxsolver_.addCols( colset );

      /* create rows */
      for( row = 0; row < nrows; ++row )
	{
	  pos = start[row];
	  rowvec.clear();
	  rowvec.add( length[row], &(index[pos]), &(elem[pos]) );
	  rowset.add( therowlb[row], rowvec, therowub[row] );
	}
      
      spxsolver_.addRows( rowset );
      // spxsolver_.changeObj( soplex::Vector( ncols, theobj ) );
      // spxsolver_.changeBounds( soplex::Vector( ncols, thecollb ), soplex::Vector( ncols, thecolub ) );
    }

  // switch sense to minimization problem
  spxsolver_.changeSense( soplex::SoPlex::MINIMIZE );

  // delete default arrays if neccessary
  if( collb == NULL )
    delete[] thecollb;
  if( colub == NULL )
    delete[] thecolub;
  if( obj == NULL )
    delete[] theobj;
  if( rowlb == NULL )
    delete[] therowlb;
  if( rowub == NULL )
    delete[] therowub;
}
			    
//-----------------------------------------------------------------------------

void
OsiSpxSolverInterface::assignProblem( OsiPackedMatrix*& matrix,
				      double*& collb, double*& colub,
				      double*& obj,
				      double*& rowlb, double*& rowub )
{
  loadProblem( *matrix, collb, colub, obj, rowlb, rowub );
  delete matrix;   matrix = 0;
  delete[] collb;  collb = 0;
  delete[] colub;  colub = 0;
  delete[] obj;    obj = 0;
  delete[] rowlb;  rowlb = 0;
  delete[] rowub;  rowub = 0;
}

//-----------------------------------------------------------------------------

void
OsiSpxSolverInterface::loadProblem( const OsiPackedMatrix& matrix,
				    const double* collb, const double* colub,
				    const double* obj,
				    const char* rowsen, const double* rowrhs,
				    const double* rowrng )
{
  int     nrows = matrix.getNumRows();
  double* rowlb = new double[nrows];
  double* rowub = new double[nrows];
  int     row;
  char    *therowsen;
  double  *therowrhs, *therowrng;
  
  if( rowsen == NULL )
    {
      therowsen = new char[nrows];
      CoinFillN( therowsen, nrows, 'G' );
    }
  else
    therowsen = const_cast<char*>(rowsen);
  if( rowrhs == NULL )
    {
      therowrhs = new double[nrows];
      CoinFillN( therowrhs, nrows, 0.0 );
    }
  else
    therowrhs = const_cast<double*>(rowrhs);
  if( rowrng == NULL )
    {
      therowrng = new double[nrows];
      CoinFillN( therowrng, nrows, 0.0 );
    }
  else
    therowrng = const_cast<double*>(rowrng);

  for( row = 0; row < nrows; ++row )
    convertSenseToBound( therowsen[row], therowrhs[row], therowrng[row], rowlb[row], rowub[row] );
  
  loadProblem( matrix, collb, colub, obj, rowlb, rowub );

  if( rowsen == NULL )
    delete[] therowsen;
  if( rowrhs == NULL )
    delete[] therowrhs;
  if( rowrng == NULL )
    delete[] therowrng;

  delete[] rowlb;
  delete[] rowub;
}
   
//-----------------------------------------------------------------------------

void
OsiSpxSolverInterface::assignProblem( OsiPackedMatrix*& matrix,
				      double*& collb, double*& colub,
				      double*& obj,
				      char*& rowsen, double*& rowrhs,
				      double*& rowrng )
{
   loadProblem( *matrix, collb, colub, obj, rowsen, rowrhs, rowrng );
   delete matrix;   matrix = 0;
   delete[] collb;  collb = 0;
   delete[] colub;  colub = 0;
   delete[] obj;    obj = 0;
   delete[] rowsen; rowsen = 0;
   delete[] rowrhs; rowrhs = 0;
   delete[] rowrng; rowrng = 0;
}

//-----------------------------------------------------------------------------

void
OsiSpxSolverInterface::loadProblem(const int numcols, const int numrows,
				   const int* start, const int* index,
				   const double* value,
				   const double* collb, const double* colub,   
				   const double* obj,
				   const double* rowlb, const double* rowub )
{
  int col, pos;
  soplex::LPColSet colset( numcols, start[numcols] );
  soplex::DSVector colvec;
  
  spxsolver_.clear();
  spxintvars_.clear();
  freeCachedData( OsiSpxSolverInterface::KEEPCACHED_NONE );

  for( col = 0; col < numcols; ++col )
    {
      pos = start[col];
      colvec.clear();
      colvec.add( start[col+1] - pos, &(index[pos]), &(value[pos]) );
      colset.add( obj[col], collb[col], colvec, colub[col] );
    }
  
  spxsolver_.addCols( colset );
  spxsolver_.changeRange( soplex::Vector( numrows, const_cast<double*>(rowlb) ), soplex::Vector( numrows, const_cast<double*>(rowub) ) );
}

//-----------------------------------------------------------------------------

void
OsiSpxSolverInterface::loadProblem(const int numcols, const int numrows,
				   const int* start, const int* index,
				   const double* value,
				   const double* collb, const double* colub,   
				   const double* obj,
				   const char* rowsen, const double* rowrhs,
				   const double* rowrng )
{
  double* rowlb = new double[numrows];
  double* rowub = new double[numrows];
  int     row;

  for( row = 0; row < numrows; ++row )
    convertSenseToBound( rowsen[row], rowrhs[row], rowrng[row], rowlb[row], rowub[row] );
  
  loadProblem( numcols, numrows, start, index, value, collb, colub, obj, rowlb, rowub );

  delete[] rowlb;
  delete[] rowub;
}
 
//-----------------------------------------------------------------------------
// Read mps files
//-----------------------------------------------------------------------------
int OsiSpxSolverInterface::readMps( const char * filename,
				    const char * extension )
{
#if 0
  std::string f(filename);
  std::string e(extension);
  std::string fullname = f + "." + e;

  spxsolver_.clear();
  spxintvars_.clear();
  if( !spxsolver_.readFile( fullname.c_str(), NULL, NULL, &spxintvars_ ) )
    throwSPXerror( "error reading file <" + fullname + ">", "readMps" );
#endif
  // just call base class method
  return OsiSolverInterface::readMps(filename,extension);
}

//-----------------------------------------------------------------------------
// Write mps files
//-----------------------------------------------------------------------------
void OsiSpxSolverInterface::writeMps( const char * filename,
				      const char * extension ) const
{
  std::string f(filename);
  std::string e(extension);
  std::string fullname = f + ".lp"; // SOPLEX cannot write MPS files yet! + "." + e;
  spxsolver_.dumpFile( fullname.c_str() );
  std::cout << "WARNING: LP file <" << fullname << "> created instead of "
	    << "MPS file <" << f << "." << e << "> !" << std::endl;
}

//#############################################################################
// Constructors, destructors clone and assignment
//#############################################################################

//-------------------------------------------------------------------
// Default Constructor 
//-------------------------------------------------------------------
OsiSpxSolverInterface::OsiSpxSolverInterface ()
  : OsiSolverInterface(),
    spxsolver_(soplex::SoPlex::ENTER, soplex::SoPlex::COLUMN), // default is primal simplex algorithm
    spxintvars_(),
    hotStartCStat_(NULL),
    hotStartCStatSize_(0),
    hotStartRStat_(NULL),
    hotStartRStatSize_(0),
    hotStartMaxIteration_(1000000), // ??? default iteration limit for strong branching is large
    obj_(NULL),
    rowsense_(NULL),
    rhs_(NULL),
    rowrange_(NULL),
    colsol_(NULL),
    rowsol_(NULL),
    redcost_(NULL),
    rowact_(NULL),
    matrixByRow_(NULL),
    matrixByCol_(NULL)
{
#ifndef NDEBUG
  soplex::Param::setVerbose( 3 );
#else
  soplex::Param::setVerbose( 0 );
#endif
}


//----------------------------------------------------------------
// Clone
//----------------------------------------------------------------
OsiSolverInterface * OsiSpxSolverInterface::clone(bool copyData) const
{
  return( new OsiSpxSolverInterface( *this ) );
}

//-------------------------------------------------------------------
// Copy constructor 
//-------------------------------------------------------------------
OsiSpxSolverInterface::OsiSpxSolverInterface( const OsiSpxSolverInterface & source )
  : OsiSolverInterface(source),
    spxsolver_(source.spxsolver_.type(), source.spxsolver_.rep()),
    spxintvars_(source.spxintvars_),
    hotStartCStat_(NULL),
    hotStartCStatSize_(0),
    hotStartRStat_(NULL),
    hotStartRStatSize_(0),
    hotStartMaxIteration_(source.hotStartMaxIteration_),
    obj_(NULL),
    rowsense_(NULL),
    rhs_(NULL),
    rowrange_(NULL),
    colsol_(NULL),
    rowsol_(NULL),
    redcost_(NULL),
    rowact_(NULL),
    matrixByRow_(NULL),
    matrixByCol_(NULL)
{
  spxsolver_.loadLP( source.spxsolver_ );
  if( source.spxsolver_.basis().status() != soplex::SPxBasis::NO_PROBLEM )
     spxsolver_.loadBasis( source.spxsolver_.basis().desc() );
  spxsolver_.setTerminationTime ( source.spxsolver_.terminationTime() );
  spxsolver_.setTerminationIter ( source.spxsolver_.terminationIter() );
  // spxsolver_.setTerminationValue( source.spxsolver_.terminationValue() ); // ???
  spxsolver_.setPricing( source.spxsolver_.pricing() );
  spxsolver_.setDelta  ( source.spxsolver_.delta()   );
  setColSolution(source.getColSolution());
  setRowPrice(source.getRowPrice());
}


//-------------------------------------------------------------------
// Destructor 
//-------------------------------------------------------------------
OsiSpxSolverInterface::~OsiSpxSolverInterface ()
{
  freeAllMemory();
}

//----------------------------------------------------------------
// Assignment operator 
//-------------------------------------------------------------------
OsiSpxSolverInterface& OsiSpxSolverInterface::operator=( const OsiSpxSolverInterface& source )
{
  if (this != &source)
    {    
      freeAllMemory();

      OsiSolverInterface::operator=( source );
      spxintvars_ = source.spxintvars_;
      
      spxsolver_.setRep        ( source.spxsolver_.rep() );
      spxsolver_.setType       ( source.spxsolver_.type() );
      spxsolver_.loadLP        ( source.spxsolver_ );
      spxsolver_.loadBasis     ( source.spxsolver_.basis().desc() );
      spxsolver_.setTerminationTime ( source.spxsolver_.terminationTime() );
      spxsolver_.setTerminationIter ( source.spxsolver_.terminationIter() );
      // spxsolver_.setTerminationValue( source.spxsolver_.terminationValue() );
      spxsolver_.setPricing    ( source.spxsolver_.pricing() );
      spxsolver_.setDelta      ( source.spxsolver_.delta()   );
      setColSolution(source.getColSolution());
      setRowPrice(source.getRowPrice());
      hotStartMaxIteration_ = source.hotStartMaxIteration_;
    }
  return *this;
}

//#############################################################################
// Applying cuts
//#############################################################################

void OsiSpxSolverInterface::applyColCut( const OsiColCut & cc )
{
  const double * soplexColLB = getColLower();
  const double * soplexColUB = getColUpper();
  const OsiPackedVector & lbs = cc.lbs();
  const OsiPackedVector & ubs = cc.ubs();
  int i;

  for( i = 0; i < lbs.getNumElements(); ++i ) 
    if ( lbs.getElements()[i] > soplexColLB[lbs.getIndices()[i]] )
      setColLower( lbs.getIndices()[i], lbs.getElements()[i] );
  for( i = 0; i < ubs.getNumElements(); ++i )
    if ( ubs.getElements()[i] < soplexColUB[ubs.getIndices()[i]] )
      setColUpper( ubs.getIndices()[i], ubs.getElements()[i] );
}

//-----------------------------------------------------------------------------

void OsiSpxSolverInterface::applyRowCut( const OsiRowCut & rowCut )
{
  addRow( rowCut.row(), rowCut.lb(), rowCut.ub() );
}

//#############################################################################
// Private methods (non-static and static)
//#############################################################################

//-------------------------------------------------------------------
/// free cached vectors

void OsiSpxSolverInterface::freeCachedColRim()
{
  delete obj_;
  obj_ = NULL;
}

void OsiSpxSolverInterface::freeCachedRowRim()
{
  freeCacheChar  ( rowsense_ );
  freeCacheDouble( rhs_ );
  freeCacheDouble( rowrange_ );
  assert( rowsense_ == NULL ); 
  assert( rhs_      == NULL ); 
  assert( rowrange_ == NULL ); 
 }

void OsiSpxSolverInterface::freeCachedMatrix()
{
  freeCacheMatrix( matrixByRow_ );
  freeCacheMatrix( matrixByCol_ );
  assert( matrixByRow_ == NULL ); 
  assert( matrixByCol_ == NULL ); 
}

void OsiSpxSolverInterface::freeCachedResults()
{
  delete colsol_;
  delete rowsol_;
  delete redcost_;
  delete rowact_;
  colsol_  = NULL;
  rowsol_  = NULL;
  redcost_ = NULL;
  rowact_  = NULL;
}


void OsiSpxSolverInterface::freeCachedData( int keepCached )
{
  if( !(keepCached & OsiSpxSolverInterface::KEEPCACHED_COLUMN) )
    freeCachedColRim();
  if( !(keepCached & OsiSpxSolverInterface::KEEPCACHED_ROW) )
    freeCachedRowRim();
  if( !(keepCached & OsiSpxSolverInterface::KEEPCACHED_MATRIX) )
    freeCachedMatrix();
  if( !(keepCached & OsiSpxSolverInterface::KEEPCACHED_RESULTS) )
    freeCachedResults();
}

void OsiSpxSolverInterface::freeAllMemory()
{
  freeCachedData();
  if( hotStartCStat_ != NULL )
    delete[] hotStartCStat_;
  if( hotStartRStat_ != NULL )
    delete[] hotStartRStat_;
  hotStartCStat_     = NULL;
  hotStartCStatSize_ = 0;
  hotStartRStat_     = NULL;
  hotStartRStatSize_ = 0;
  spxsolver_.clear();
  spxintvars_.clear();
}

#endif
