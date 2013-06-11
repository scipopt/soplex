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

/**@file  soplex2.h
 * @brief preconfigured \ref soplex::SoPlex "SoPlex" LP-solver.
 */
#ifndef _SOPLEX2_H_
#define _SOPLEX2_H_

#include <string>

///@todo SoPlex2 should also have an spxout object to avoid using a global one
#include "rational.h"
#include "spxsolver.h"
#include "slufactor.h"

///@todo try to move to cpp file by forward declaration
#include "spxsimplifier.h"
#include "spxmainsm.h"

#include "spxscaler.h"
#include "spxequilisc.h"
#include "spxgeometsc.h"

#include "spxstarter.h"
#include "spxweightst.h"
#include "spxsumst.h"
#include "spxvectorst.h"

#include "spxpricer.h"
#include "spxdantzigpr.h"
#include "spxparmultpr.h"
#include "spxdevexpr.h"
#include "spxsteeppr.h"
#include "spxsteepexpr.h"
#include "spxhybridpr.h"

#include "spxratiotester.h"
#include "spxdefaultrt.h"
#include "spxharrisrt.h"
#include "spxfastrt.h"
#include "spxboundflippingrt.h"


///@todo maximum line length in LP and MPS reader/writer should be 6553? (Dan)
///@todo implement interface to rational LP, including rational basis (Ambros)

///@todo draw flow chart of main solving loop and performInfeasibilityIR() (Dan)
///@todo solution structure (primal, dual, basis, maxviolation) and record and return "best" solutions found during IR (Ambros)
///@todo implement main IR loop for primal and dual feasible case with fail otherwise (Ambros)
///@todo implement statistical info (time, factor time, iters, ...) since last call to solveReal() or solveRational() (Ambros?)
///@todo implement performInfeasibilityIR (Ambros?)
///@todo extend IR loop to infeasible case (Dan?)
///@todo extend IR loop to unbounded case (Dan?)

///@todo interface rational reconstruction code for rational vectors
///@todo integrate rational reconstruction into IR loop
///@todo templatize SPxSolver and necessary components (SLUFactor, pricer, ratiotester)
///@todo integrate rational SPxSolver and distinguish between original and transformed rational LP
///@todo rational scalers
///@todo rational simplifier


namespace soplex
{

/**@class SoPlex2
 * @brief   Preconfigured SoPlex LP-solver.
 * @ingroup Algo
 */
class SoPlex2
{
public:

   //**@name Construction and destruction */
   //@{

   /// default constructor
   SoPlex2();

   /// assignment operator
   SoPlex2& operator=(const SoPlex2& rhs);

   /// copy constructor
   SoPlex2(const SoPlex2& rhs);

   /// destructor
   ~SoPlex2();

   //@}


   //**@name Access of the real LP */
   //@{

   /// returns number of rows
   int numRowsReal() const;

   /// returns number of columns
   int numColsReal() const;

   /// returns number of nonzeros
   int numNonzerosReal() const;

   /// returns smallest non-zero element in absolute value
   Real minAbsNonzeroReal() const;

   /// returns biggest non-zero element in absolute value
   Real maxAbsNonzeroReal() const;

   /// returns row identifier for row \p i
   SPxRowId rowIdReal(int i) const;

   /// returns column identifier for column \p i
   SPxColId colIdReal(int i) const;

   /// returns index of the row with identifier \p id
   int idxReal(const SPxRowId& id) const;

   /// returns index of the column with identifier \p id
   int idxReal(const SPxColId& id) const;

   /// returns index of the row or column with identifier \p id
   int idxReal(const SPxId& id) const;

   /// gets row \p i
   void getRowReal(int i, LPRowReal& lprow) const;

   /// gets row with identifier \p id
   void getRowReal(const SPxRowId& id, LPRowReal& lprow) const;

   /// gets rows \p start, ..., \p end.
   void getRowsReal(int start, int end, LPRowSetReal& lprowset) const;

   /// returns vector of row \p i
   const SVectorReal& rowVectorReal(int i) const;

   /// returns vector of row with identifier \p id
   const SVectorReal& rowVectorReal(const SPxRowId& id) const;

   /// returns right-hand side vector
   const VectorReal& rhsReal() const;

   /// returns right-hand side of row \p i
   Real rhsReal(int i) const;

   /// returns right-hand side of row with identifier \p id
   Real rhsReal(const SPxRowId& id) const;

   /// returns left-hand side vector
   const VectorReal& lhsReal() const;

   /// returns left-hand side of row \p i
   Real lhsReal(int i) const;

   /// returns left-hand side of row with identifier \p id
   Real lhsReal(const SPxRowId& id) const;

   /// returns inequality type of row \p i
   LPRowReal::Type rowTypeReal(int i) const;

   /// returns inequality type of row with identifier \p id
   LPRowReal::Type rowTypeReal(const SPxRowId& id) const;

   /// gets column \p i
   void getColReal(int i, LPColReal& lpcol) const;

   /// gets column with identifier \p id.
   void getColReal(const SPxColId& id, LPColReal& lpcol) const;

   /// gets columns \p start, ..., \p end
   void getColsReal(int start, int end, LPColSetReal& lpcolset) const;

   /// returns vector of column \p i
   const SVectorReal& colVectorReal(int i) const;

   /// returns vector of column with identifier \p id
   const SVectorReal& colVectorReal(const SPxColId& id) const;

   /// returns upper bound vector
   const VectorReal& upperReal() const;

   /// returns upper bound of column \p i
   Real upperReal(int i) const;

   /// returns upper bound of column with identifier \p id
   Real upperReal(const SPxColId& id) const;

   /// returns lower bound vector
   const VectorReal& lowerReal() const;

   /// returns lower bound of column \p i
   Real lowerReal(int i) const;

   /// returns lower bound of column with identifier \p id
   Real lowerReal(const SPxColId& id) const;

   /// gets objective function vector
   void getObjReal(VectorReal& obj) const;

   /// returns objective value of column \p i
   Real objReal(int i) const;

   /// returns objective value of column with identifier \p id
   Real objReal(const SPxColId& id) const;

   /// returns objective function vector after transformation to a maximization problem; since this is how it is stored
   /// internally, this is generally faster
   const VectorReal& maxObjReal() const;

   /// returns objective value of column \p i after transformation to a maximization problem; since this is how it is
   /// stored internally, this is generally faster
   Real maxObjReal(int i) const;

   /// returns objective value of column with identifier \p id after transformation to a maximization problem; since
   /// this is how it is stored internally, this is generally faster
   Real maxObjReal(const SPxColId& id) const;

   //@}


   //**@name Access of the rational LP */
   //@{

   /// returns number of rows
   int numRowsRational() const;

   /// returns number of columns
   int numColsRational() const;

   /// returns number of nonzeros
   int numNonzerosRational() const;

   /// returns smallest non-zero element in absolute value
   Rational minAbsNonzeroRational() const;

   /// returns biggest non-zero element in absolute value
   Rational maxAbsNonzeroRational() const;

   /// returns row identifier for row \p i
   SPxRowId rowIdRational(int i) const;

   /// returns column identifier for column \p i
   SPxColId colIdRational(int i) const;

   /// returns index of the row with identifier \p id
   int idxRational(const SPxRowId& id) const;

   /// returns index of the column with identifier \p id
   int idxRational(const SPxColId& id) const;

   /// returns index of the row or column with identifier \p id
   int idxRational(const SPxId& id) const;

   /// gets row \p i
   void getRowRational(int i, LPRowRational& lprow) const;

   /// gets row with identifier \p id
   void getRowRational(const SPxRowId& id, LPRowRational& lprow) const;

   /// gets rows \p start, ..., \p end.
   void getRowsRational(int start, int end, LPRowSetRational& lprowset) const;

   /// returns vector of row \p i
   const SVectorRational& rowVectorRational(int i) const;

   /// returns vector of row with identifier \p id
   const SVectorRational& rowVectorRational(const SPxRowId& id) const;

   /// returns right-hand side vector
   const VectorRational& rhsRational() const;

   /// returns right-hand side of row \p i
   Rational rhsRational(int i) const;

   /// returns right-hand side of row with identifier \p id
   Rational rhsRational(const SPxRowId& id) const;

   /// returns left-hand side vector
   const VectorRational& lhsRational() const;

   /// returns left-hand side of row \p i
   Rational lhsRational(int i) const;

   /// returns left-hand side of row with identifier \p id
   Rational lhsRational(const SPxRowId& id) const;

   /// returns inequality type of row \p i
   LPRowRational::Type rowTypeRational(int i) const;

   /// returns inequality type of row with identifier \p id
   LPRowRational::Type rowTypeRational(const SPxRowId& id) const;

   /// gets column \p i
   void getColRational(int i, LPColRational& lpcol) const;

   /// gets column with identifier \p id.
   void getColRational(const SPxColId& id, LPColRational& lpcol) const;

   /// gets columns \p start, ..., \p end
   void getColsRational(int start, int end, LPColSetRational& lpcolset) const;

   /// returns vector of column \p i
   const SVectorRational& colVectorRational(int i) const;

   /// returns vector of column with identifier \p id
   const SVectorRational& colVectorRational(const SPxColId& id) const;

   /// returns upper bound vector
   const VectorRational& upperRational() const;

   /// returns upper bound of column \p i
   Rational upperRational(int i) const;

   /// returns upper bound of column with identifier \p id
   Rational upperRational(const SPxColId& id) const;

   /// returns lower bound vector
   const VectorRational& lowerRational() const;

   /// returns lower bound of column \p i
   Rational lowerRational(int i) const;

   /// returns lower bound of column with identifier \p id
   Rational lowerRational(const SPxColId& id) const;

   /// gets objective function vector
   void getObjRational(VectorRational& obj) const;

   /// returns objective value of column \p i
   Rational objRational(int i) const;

   /// returns objective value of column with identifier \p id
   Rational objRational(const SPxColId& id) const;

   /// returns objective function vector after transformation to a maximization problem; since this is how it is stored
   /// internally, this is generally faster
   const VectorRational& maxObjRational() const;

   /// returns objective value of column \p i after transformation to a maximization problem; since this is how it is
   /// stored internally, this is generally faster
   Rational maxObjRational(int i) const;

   /// returns objective value of column with identifier \p id after transformation to a maximization problem; since
   /// this is how it is stored internally, this is generally faster
   Rational maxObjRational(const SPxColId& id) const;

   //@}


   //**@name Modification of the real LP */
   //@{

   /// adds a single row
   void addRowReal(const LPRowReal& lprow);

   /// adds a single row and gets its \p id
   void addRowReal(SPxRowId& id, const LPRowReal& lprow);

   /// adds multiple rows
   void addRowsReal(const LPRowSetReal& lprowset);

   /// adds multiple rows and gets an array of their \p id 's
   void addRowsReal(SPxRowId id[], const LPRowSetReal& lprowset);

   /// adds a single column
   void addColReal(const LPCol& lpcol);

   /// adds a single column and gets its \p id
   void addColReal(SPxColId& id, const LPCol& lpcol);

   /// adds multiple columns
   void addColsReal(const LPColSetReal& lpcolset);

   /// adds multiple columns and gets an array of their \p id 's
   void addColsReal(SPxColId id[], const LPColSetReal& lpcolset);

   /// replaces row \p i with \p lprow
   void changeRowReal(int i, const LPRowReal& lprow);

   /// replaces row with identifier \p id with \p lprow
   void changeRowReal(SPxRowId id, const LPRowReal& lprow);

   /// changes left-hand side vector for constraints to \p lhs
   void changeLhsReal(const VectorReal& lhs);

   /// changes left-hand side of row \p i to \p lhs
   void changeLhsReal(int i, Real lhs);

   /// changes left-hand side of row with identifier \p id to \p lhs
   void changeLhsReal(SPxRowId id, Real lhs);

   /// changes right-hand side vector to \p rhs
   void changeRhsReal(const VectorReal& rhs);

   /// changes right-hand side of row \p i to \p rhs
   void changeRhsReal(int i, Real rhs);

   /// changes right-hand of row with identifier \p id to \p rhs
   void changeRhsReal(SPxRowId id, Real rhs);

   /// changes left- and right-hand side vectors
   void changeRangeReal(const VectorReal& lhs, const VectorReal& rhs);

   /// changes left- and right-hand side of row \p i
   void changeRangeReal(int i, Real lhs, Real rhs);

   /// changes left- and right-hand side of row with identifier \p id
   void changeRangeReal(SPxRowId id, Real lhs, Real rhs);

   /// replaces column \p i with \p lpcol
   void changeColReal(int i, const LPColReal& lpcol);

   /// replaces column with identifier \p id with \p lpcol
   void changeColReal(SPxColId id, const LPColReal& lpcol);

   /// changes vector of lower bounds to \p lower
   void changeLowerReal(const VectorReal& lower);

   /// changes lower bound of column i to \p lower
   void changeLowerReal(int i, Real lower);

   /// changes lower bound of column with identifier \p id to \p lower
   void changeLowerReal(SPxColId id, Real lower);

   /// changes vector of upper bounds to \p upper
   void changeUpperReal(const VectorReal& upper);

   /// changes \p i 'th upper bound to \p upper
   void changeUpperReal(int i, Real upper);

   /// changes upper bound of column with identifier \p id to \p upper
   void changeUpperReal(SPxColId id, Real upper);

   /// changes vectors of column bounds to \p lower and \p upper
   void changeBoundsReal(const VectorReal& lower, const VectorReal& upper);

   /// changes bounds of column \p i to \p lower and \p upper
   void changeBoundsReal(int i, Real lower, Real upper);

   /// changes bounds of column with identifier \p id to \p lower and \p upper
   void changeBoundsReal(SPxColId id, Real lower, Real upper);

   /// changes objective function vector to \p obj
   void changeObjReal(const VectorReal& obj);

   /// changes objective coefficient of column i to \p obj
   void changeObjReal(int i, Real obj);

   /// changes objective coefficient of column with identifier \p id to \p obj
   void changeObjReal(SPxColId id, Real obj);

   /// changes matrix entry in row \p i and column \p j to \p val
   void changeElementReal(int i, int j, Real val);

   /// changes matrix entry identified by (\p rowid, \p colid) to \p val
   void changeElementReal(SPxRowId rowid, SPxColId colid, Real val);

   /// removes row \p i
   void removeRowReal(int i);

   /// removes row with identifier \p id
   void removeRowReal(SPxRowId id);

   /// removes all rows with an index \p i such that \p perm[i] < 0; upon completion, \p perm[i] >= 0 indicates the
   /// new index where row \p i has been moved to; note that \p perm must point to an array of size at least
   /// #numRowsReal()
   void removeRowsReal(int perm[]);

   /// remove all rows with identifier in array \p id of size \p n; an array \p perm of size #numRowsReal() may be
   /// passed as buffer memory
   void removeRowsReal(SPxRowId id[], int n, int perm[] = 0);

   /// remove all rows with indices in array \p idx of size \p n; an array \p perm of size #numRowsReal() may be passed
   /// as buffer memory
   void removeRowsReal(int idx[], int n, int perm[] = 0);

   /// removes rows \p start to \p end including both; an array \p perm of size #numRowsReal() may be passed as buffer
   /// memory
   void removeRowRangeReal(int start, int end, int perm[] = 0);

   /// removes column i
   void removeColReal(int i);

   /// removes column with identifier \p id
   void removeColReal(SPxColId id);

   /// removes all columns with an index \p i such that \p perm[i] < 0; upon completion, \p perm[i] >= 0 indicates the
   /// new index where column \p i has been moved to; note that \p perm must point to an array of size at least
   /// #numColsReal()
   void removeColsReal(int perm[]);

   /// remove all columns with identifier in array \p id of size \p n; an array \p perm of size #numColsReal() may be
   /// passed as buffer memory
   void removeColsReal(SPxColId id[], int n, int perm[] = 0);

   /// remove all columns with indices in array \p idx of size \p n; an array \p perm of size #numColsReal() may be
   /// passed as buffer memory
   void removeColsReal(int idx[], int n, int perm[] = 0);

   /// removes columns \p start to \p end including both; an array \p perm of size #numColsReal() may be passed as
   /// buffer memory
   void removeColRangeReal(int start, int end, int perm[] = 0);

   /// clears the LP
   void clearLPReal();

   //@}


   //**@name Modification of the rational LP */
   //@{

   //@}


   //**@name Solving and solution query for the real LP */
   //@{

   /// solves real LP
   SPxSolver::Status solveReal();

   /// returns the current status
   SPxSolver::Status statusReal() const;

   /// returns the objective value if a primal solution is available
   Real objValueReal() const;

   /// is a primal feasible solution available?
   bool hasPrimalReal() const;

   /// gets the primal solution vector if available; returns true on success
   bool getPrimalReal(VectorReal& vector) const;

   /// gets the vector of slack values if available; returns true on success
   bool getSlacksReal(VectorReal& vector) const;

   /// gets the primal ray if LP is unbounded; returns true on success
   bool getPrimalrayReal(VectorReal& vector) const;

   /// is a dual feasible solution available?
   bool hasDualReal() const;

   /// gets the dual solution vector if available; returns true on success
   bool getDualReal(VectorReal& vector) const;

   /// gets the vector of reduced cost values if available; returns true on success
   bool getRedcostReal(VectorReal& vector) const;

   /// gets the Farkas proof if LP is infeasible; returns true on success
   bool getDualfarkasReal(VectorReal& vector) const;

   /// gets violation of bounds by given primal solution
   void getBoundViolationReal(VectorReal& primal, Real& maxviol, Real& sumviol) const;

   /// gets violation of constraints by given primal solution
   void getConstraintViolationReal(VectorReal& primal, Real& maxviol, Real& sumviol) const;

   //@}


   //**@name Solving and solution query for the rational LP */
   //@{

   /// solves rational LP
   SPxSolver::Status solveRational();

   //@}


   //**@name Basis information for the real LP */
   //@{

   /// is an advanced starting basis available?
   bool hasBasisReal() const;

   /// returns basis status for a single row
   SPxSolver::VarStatus basisRowStatusReal(int row) const;

   /// returns basis status for a single row
   SPxSolver::VarStatus basisRowStatusReal(const SPxRowId& id) const;

   /// returns basis status for a single column
   SPxSolver::VarStatus basisColStatusReal(int col) const;

   /// returns basis status for a single column
   SPxSolver::VarStatus basisColStatusReal(const SPxColId& id) const;

   /// gets current basis and returns solver status
   void getBasisReal(SPxSolver::VarStatus rows[], SPxSolver::VarStatus cols[]) const;

   /// sets starting basis via arrays of statuses
   void setBasisReal(SPxSolver::VarStatus rows[], SPxSolver::VarStatus cols[]);

   /// clears starting basis
   void clearBasisReal();

   //@}

#if 0
   //**@name Statistical information */
   //@{

   /// time spent in factorizations
   Real factorTime() const;

   /// number of factorizations performed
   int factorCount() const;

   /// time spent in solves
   Real luSolveTime() const;

   /// number of solves performed
   int luSolveCount() const;

   /// number of iterations since last call to solve
   int numIterations() const;

   /// statistical information in form of a string
   std::string statisticString() const;

   //@}
#endif


   /// name of starter
   const char* getStarterName();

   /// name of simplifier
   const char* getSimplifierName();

   /// name of scaling method before simplifier
   const char* getFirstScalerName();

   /// name of scaling method after simplifier
   const char* getSecondScalerName();

   /// name of currently loaded pricer
   const char* getPricerName();

   /// name of currently loaded ratiotester
   const char* getRatiotesterName();

   //**@name I/O for the real LP */
   //@{

   /// reads real LP in LP or MPS format from file and returns true on success; gets row names, column names, and
   /// integer variables if desired
   bool readFileReal(const char* filename, NameSet* rowNames = 0, NameSet* colNames = 0, DIdxSet* intVars = 0);

   /// writes real LP to file; LP or MPS format is chosen from the extension in \p filename; if \p rowNames and \p
   /// colNames are \c NULL, default names are used; if \p intVars is not \c NULL, the variables contained in it are
   /// marked as integer
   void writeFileReal(const char* filename, const NameSet* rowNames = 0, const NameSet* colNames = 0, const DIdxSet* intvars = 0) const;

   /// reads basis information from \p filename and returns true on success; if \p rowNames and \p colNames are \c NULL,
   /// default names are assumed
   bool readBasisFileReal(const char* filename, const NameSet* rowNames = 0, const NameSet* colNames = 0);

   /// writes basis information to \p filename; if \p rowNames and \p colNames are \c NULL, default names are used
   void writeBasisFileReal(const char* filename, const NameSet* rowNames, const NameSet* colNames);

#if 0
   /// writes internal LP, basis information, and parameter settings; if \p rowNames and \p colNames are \c NULL,
   /// default names are used
   void writeStateReal(const char* filename, const NameSet* rowNames = 0, const NameSet* colNames = 0) const;
#endif

   //@}


   //**@name I/O for the rational LP */
   //@{

   /// reads rational LP in LP or MPS format from file and returns true on success; gets row names, column names, and
   /// integer variables if desired
   bool readFileRational(const char* filename, NameSet* rowNames = 0, NameSet* colNames = 0, DIdxSet* intVars = 0);

   /// writes rational LP to file; LP or MPS format is chosen from the extension in \p filename; if \p rowNames and \p
   /// colNames are \c NULL, default names are used; if \p intVars is not \c NULL, the variables contained in it are
   /// marked as integer
   void writeFileRational(const char* filename, const NameSet* rowNames = 0, const NameSet* colNames = 0, const DIdxSet* intvars = 0) const;

   /// reads basis information from \p filename and returns true on success; if \p rowNames and \p colNames are \c NULL,
   /// default names are assumed
   bool readBasisFileRational(const char* filename, const NameSet* rowNames = 0, const NameSet* colNames = 0);

   /// writes basis information to \p filename; if \p rowNames and \p colNames are \c NULL, default names are used
   void writeBasisFileRational(const char* filename, const NameSet* rowNames, const NameSet* colNames);

#if 0
   /// writes internal LP, basis information, and parameter settings; if \p rowNames and \p colNames are \c NULL,
   /// default names are used
   void writeStateRational(const char* filename, const NameSet* rowNames = 0, const NameSet* colNames = 0) const;
#endif

   //@}


   //**@name Parameters */
   //@{

   /// boolean parameters
   typedef enum
   {
      /// should partial pricing be used?
      PARTIAL_PRICING = 0,

      /// number of boolean parameters
      BOOLPARAM_COUNT = 1
   } BoolParam;

   /// integer parameters
   typedef enum
   {
      /// objective sense
      OBJSENSE = 0,

      /// type of computational form, i.e., column or row representation
      REPRESENTATION = 1,

      /// type of algorithm, i.e., enter or leave
      ALGORITHM = 2,

      /// type of LU update
      FACTOR_UPDATE_TYPE = 3,

      /// maximum number of updates before fresh factorization
      FACTOR_UPDATE_MAX = 4,

      /// iteration limit (-1 if unlimited)
      ITERLIMIT = 5,

      /// display frequency
      DISPLAY_FREQ = 6,

      /// verbosity level
      VERBOSITY = 7,

      /// type of simplifier
      SIMPLIFIER = 8,

      /// type of scaler applied before simplification
      SCALER_BEFORE_SIMPLIFIER = 9,

      /// type of scaler applied after simplification
      SCALER_AFTER_SIMPLIFIER = 10,

      /// type of starter used to create crash basis
      STARTER = 11,

      /// type of pricer
      PRICER = 12,

      /// type of ratio test
      RATIOTESTER = 13,

      /// number of integer parameters
      INTPARAM_COUNT = 14
   } IntParam;

   /// values for parameter OBJSENSE
   enum
   {
      /// minimization
      OBJSENSE_MINIMIZE = -1,

      /// maximization
      OBJSENSE_MAXIMIZE = 1
   };

   /// values for parameter REPRESENTATION
   enum
   {
      /// column representation Ax - s = 0, lower <= x <= upper, lhs <= s <= rhs
      REPRESENTATION_COLUMN = 0,

      /// row representation (lower,lhs) <= (x,Ax) <= (upper,rhs)
      REPRESENTATION_ROW = 1
   };

   /// values for parameter ALGORITHM
   enum
   {
      /// entering algorithm, i.e., primal simplex for column and dual simplex for row representation
      ALGORITHM_ENTER = 0,

      /// leaving algorithm, i.e., dual simplex for column and primal simplex for row representation
      ALGORITHM_LEAVE = 1
   };

   /// values for parameter FACTOR_UPDATE_TYPE
   enum
   {
      /// product form update
      FACTOR_UPDATE_TYPE_ETA = 0,

      /// Forrest-Tomlin type update
      FACTOR_UPDATE_TYPE_FT = 1
   };

   /// values for parameter VERBOSITY
   enum
   {
      /// only error output
      VERBOSITY_ERROR = 0,

      /// only error and warning output
      VERBOSITY_WARNING = 1,

      /// only error, warning, and debug output
      VERBOSITY_DEBUG = 2,

      /// standard verbosity level
      VERBOSITY_NORMAL = 3,

      /// high verbosity level
      VERBOSITY_HIGH = 4,

      /// full verbosity level
      VERBOSITY_FULL = 5
   };

   /// values for parameter SIMPLIFIER
   enum
   {
      /// no simplifier
      SIMPLIFIER_OFF = 0,

      /// automatic choice
      SIMPLIFIER_AUTO = 1
   };

   /// values for parameter SCALER
   enum
   {
      /// no scaler
      SCALER_OFF = 0,

      /// equilibrium scaling on rows or columns
      SCALER_UNIEQUI = 1,

      /// equilibrium scaling on rows and columns
      SCALER_BIEQUI = 2,

      /// geometric mean scaling on rows and columns, max 1 round
      SCALER_GEO1 = 3,

      /// geometric mean scaling on rows and columns, max 8 rounds
      SCALER_GEO8 = 4,
   };

   /// values for parameter STARTER
   enum
   {
      /// slack basis
      STARTER_OFF = 0,

      /// greedy crash basis weighted by objective, bounds, and sides
      STARTER_WEIGHT = 1,

      /// crash basis from a greedy solution
      STARTER_SUM = 2,

      /// generic solution-based crash basis
      STARTER_VECTOR = 3
   };

   ///@todo is the order different than usual on purpose?
   /// values for parameter PRICER
   enum
   {
      /// automatic pricer
      PRICER_AUTO = 0,

      /// Dantzig pricer
      PRICER_DANTZIG = 1,

      /// partial multiple pricer based on Dantzig pricing
      PRICER_PARMULT = 2,

      /// devex pricer
      PRICER_DEVEX = 3,

      /// steepest edge pricer with initialization to unit norms
      PRICER_QUICKSTEEP = 4,

      /// steepest edge pricer with exact initialization of norms
      PRICER_STEEP = 5,

      /// hyprid pricer choosing between quicksteep and partial multiple pricer
      PRICER_HYBRID = 6
   };

   /// values for parameter RATIOTESTER
   enum
   {
      /// textbook ratio test without stabilization
      RATIOTESTER_TEXTBOOK = 0,

      /// standard Harris ratio test
      RATIOTESTER_HARRIS = 1,

      /// modified Harris ratio test
      RATIOTESTER_FAST = 2,

      /// bound flipping ratio test for long steps in the dual simplex
      RATIOTESTER_BOUNDFLIPPING = 3
   };

   /// real parameters
   typedef enum
   {
      /// general zero tolerance
      EPSILON_ZERO = 0,

      /// zero tolerance used in factorization
      EPSILON_FACTORIZATION = 1,

      /// zero tolerance used in factorization update
      EPSILON_UPDATE = 2,

      /// infinity threshold
      INFTY = 3,

      /// time limit in seconds (INFTY if unlimited)
      TIMELIMIT = 4,

      /// lower limit on objective value
      OBJLIMIT_LOWER = 5,

      /// upper limit on objective value
      OBJLIMIT_UPPER = 6,

      /// threshold for activating iterative refinement
      IRTHRESHOLD = 7,

      /// number of real parameters
      REALPARAM_COUNT = 8
   } RealParam;

   /// rational parameters
   typedef enum
   {
      /// primal feasibility tolerance
      FEASTOL,

      /// dual feasibility tolerance
      OPTTOL,

      /// number of rational parameters
      RATIONALPARAM_COUNT = 2
   } RationalParam;

   /// class of parameter settings
   class Settings;

   /// returns boolean parameter value
   bool boolParam(const BoolParam param) const;

   /// returns integer parameter value
   int intParam(const IntParam param) const;

   /// returns real parameter value
   Real realParam(const RealParam param) const;

   /// returns rational parameter value
   Rational rationalParam(const RationalParam param) const;

   /// returns current parameter settings
   const Settings& settings() const;

   /// sets boolean parameter value; returns true on success
   bool setBoolParam(const BoolParam param, const bool value, const bool quiet = false, const bool init = false);

   /// sets integer parameter value; returns true on success
   bool setIntParam(const IntParam param, const int value, const bool quiet = false, const bool init = false);

   /// sets real parameter value; returns true on success
   bool setRealParam(const RealParam param, const Real value, const bool quiet = false, const bool init = false);

   /// sets rational parameter value; returns true on success
   bool setRationalParam(const RationalParam param, const Rational value, const bool quiet = false, const bool init = false);

   /// sets parameter settings; returns true on success
   bool setSettings(const Settings& settings, const bool quiet = false, const bool init = false);

   //@}


private:
   Settings* _currentSettings;

   SPxSolver _solver;
   SLUFactor _slufactor;
   SPxMainSM _simplifierMainSM;
   SPxEquiliSC _scalerUniequi;
   SPxEquiliSC _scalerBiequi;
   SPxGeometSC _scalerGeo1;
   SPxGeometSC _scalerGeo8;
   SPxWeightST _starterWeight;
   SPxSumST _starterSum;
   SPxVectorST _starterVector;
   SPxDantzigPR _pricerDantzig;
   SPxParMultPR _pricerParMult;
   SPxDevexPR _pricerDevex;
   SPxSteepPR _pricerQuickSteep;
   SPxSteepExPR _pricerSteep;
   SPxHybridPR _pricerHybrid;
   SPxDefaultRT _ratiotesterTextbook;
   SPxHarrisRT _ratiotesterHarris;
   SPxFastRT _ratiotesterFast;
   SPxBoundFlippingRT _ratiotesterBoundFlipping;

   DataArray< SPxSolver::VarStatus > _basisStatusRowsReal;
   DataArray< SPxSolver::VarStatus > _basisStatusColsReal;

   SPxSimplifier* _simplifier;
   SPxScaler* _firstScaler;
   SPxScaler* _secondScaler;
   SPxStarter* _starter;

   SPxLPReal* _realLP;
   SPxLPRational* _rationalLP;

   bool _isRealLPLoaded;
   bool _hasBasisReal;
   bool _hasBasisRational;

   /// checks consistency
   bool _isConsistent() const;

   /// creates a permutation for removing rows/columns from an array of IDs
   void _idToPerm(SPxId* id, int idSize, int* perm, int permSize) const;

   /// creates a permutation for removing rows/columns from an array of indices
   void _idxToPerm(int* idx, int idxSize, int* perm, int permSize) const;

   /// creates a permutation for removing rows/columns from a range of indices
   void _rangeToPerm(int start, int end, int* perm, int permSize) const;
};
} // namespace soplex
#endif // _SOPLEX2_H_

//-----------------------------------------------------------------------------
//Emacs Local Variables:
//Emacs mode:c++
//Emacs c-basic-offset:3
//Emacs tab-width:8
//Emacs indent-tabs-mode:nil
//Emacs End:
//-----------------------------------------------------------------------------
