/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the class library                   */
/*       SoPlex --- the Sequential object-oriented simPlex.                  */
/*                                                                           */
/*    Copyright (C) 1996-2015 Konrad-Zuse-Zentrum                            */
/*                            fuer Informationstechnik Berlin                */
/*                                                                           */
/*  SoPlex is distributed under the terms of the ZIB Academic Licence.       */
/*                                                                           */
/*  You should have received a copy of the ZIB Academic License              */
/*  along with SoPlex; see the file COPYING. If not email to soplex@zib.de.  */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/**@file  soplex.h
 * @brief Preconfigured SoPlex LP solver
 */

#ifndef _SOPLEX_H_
#define _SOPLEX_H_

#ifndef SOPLEX_LEGACY
#include <string.h>

///@todo SoPlex should also have an spxout object to avoid using a global one
#include "spxdefines.h"
#include "basevectors.h"
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
#include "spxautopr.h"
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

#include "sol.h"

///@todo implement automatic rep switch, based on row/col dim
///@todo introduce status codes for SoPlex, especially for rational solving

///@todo record and return "best" solutions found during IR (Ambros)
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

/**@class SoPlex
 * @brief   Preconfigured SoPlex LP-solver.
 * @ingroup Algo
 */
class SoPlex
{
public:

   //**@name Construction and destruction */
   //@{

   /// default constructor
   SoPlex();

   /// assignment operator
   SoPlex& operator=(const SoPlex& rhs);

   /// copy constructor
   SoPlex(const SoPlex& rhs);

   /// destructor
   virtual ~SoPlex();

   //@}


   /// message handler
//   SPxOut spxout;


   //**@name Access to the real LP */
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

   /// gets row \p i
   void getRowReal(int i, LPRowReal& lprow) const;

   /// gets rows \p start, ..., \p end.
   void getRowsReal(int start, int end, LPRowSetReal& lprowset) const;

   /// returns vector of row \p i
   const SVectorReal& rowVectorReal(int i) const;

   /// returns right-hand side vector
   const VectorReal& rhsReal() const;

   /// returns right-hand side of row \p i
   Real rhsReal(int i) const;

   /// returns left-hand side vector
   const VectorReal& lhsReal() const;

   /// returns left-hand side of row \p i
   Real lhsReal(int i) const;

   /// returns inequality type of row \p i
   LPRowReal::Type rowTypeReal(int i) const;

   /// gets column \p i
   void getColReal(int i, LPColReal& lpcol) const;

   /// gets columns \p start, ..., \p end
   void getColsReal(int start, int end, LPColSetReal& lpcolset) const;

   /// returns vector of column \p i
   const SVectorReal& colVectorReal(int i) const;

   /// returns upper bound vector
   const VectorReal& upperReal() const;

   /// returns upper bound of column \p i
   Real upperReal(int i) const;

   /// returns lower bound vector
   const VectorReal& lowerReal() const;

   /// returns lower bound of column \p i
   Real lowerReal(int i) const;

   /// gets objective function vector
   void getObjReal(VectorReal& obj) const;

   /// returns objective value of column \p i
   Real objReal(int i) const;

   /// returns objective function vector after transformation to a maximization problem; since this is how it is stored
   /// internally, this is generally faster
   const VectorReal& maxObjReal() const;

   /// returns objective value of column \p i after transformation to a maximization problem; since this is how it is
   /// stored internally, this is generally faster
   Real maxObjReal(int i) const;

   /// gets number of available dual norms
   void getNdualNorms(int& nnormsRow, int& nnormsCol) const;

   /// gets steepest edge norms and returns false if they are not available
   bool getDualNorms(int& nnormsRow, int& nnormsCol, Real* norms) const;

   /// sets steepest edge norms and returns false if that's not possible
   bool setDualNorms(int nnormsRow, int nnormsCol, Real* norms);

   //@}


   //**@name Access to the rational LP */
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

   /// gets row \p i
   void getRowRational(int i, LPRowRational& lprow) const;

   /// gets rows \p start, ..., \p end.
   void getRowsRational(int start, int end, LPRowSetRational& lprowset) const;

   /// returns vector of row \p i
   const SVectorRational& rowVectorRational(int i) const;

   /// returns right-hand side vector
   const VectorRational& rhsRational() const;

   /// returns right-hand side of row \p i
   const Rational& rhsRational(int i) const;

   /// returns left-hand side vector
   const VectorRational& lhsRational() const;

   /// returns left-hand side of row \p i
   const Rational& lhsRational(int i) const;

   /// returns inequality type of row \p i
   LPRowRational::Type rowTypeRational(int i) const;

   /// gets column \p i
   void getColRational(int i, LPColRational& lpcol) const;

   /// gets columns \p start, ..., \p end
   void getColsRational(int start, int end, LPColSetRational& lpcolset) const;

   /// returns vector of column \p i
   const SVectorRational& colVectorRational(int i) const;

   /// returns upper bound vector
   const VectorRational& upperRational() const;

   /// returns upper bound of column \p i
   const Rational& upperRational(int i) const;

   /// returns lower bound vector
   const VectorRational& lowerRational() const;

   /// returns lower bound of column \p i
   const Rational& lowerRational(int i) const;

   /// gets objective function vector
   void getObjRational(VectorRational& obj) const;

   /// gets objective value of column \p i
   void getObjRational(int i, Rational& obj) const;

   /// returns objective value of column \p i
   Rational objRational(int i) const;

   /// returns objective function vector after transformation to a maximization problem; since this is how it is stored
   /// internally, this is generally faster
   const VectorRational& maxObjRational() const;

   /// returns objective value of column \p i after transformation to a maximization problem; since this is how it is
   /// stored internally, this is generally faster
   const Rational& maxObjRational(int i) const;

   //@}


   //**@name Modification of the real LP */
   //@{

   /// adds a single row
   void addRowReal(const LPRowReal& lprow);

   /// adds multiple rows
   void addRowsReal(const LPRowSetReal& lprowset);

   /// adds a single column
   void addColReal(const LPCol& lpcol);

   /// adds multiple columns
   void addColsReal(const LPColSetReal& lpcolset);

   /// replaces row \p i with \p lprow
   void changeRowReal(int i, const LPRowReal& lprow);

   /// changes left-hand side vector for constraints to \p lhs
   void changeLhsReal(const VectorReal& lhs);

   /// changes left-hand side of row \p i to \p lhs
   void changeLhsReal(int i, const Real& lhs);

   /// changes right-hand side vector to \p rhs
   void changeRhsReal(const VectorReal& rhs);

   /// changes right-hand side of row \p i to \p rhs
   void changeRhsReal(int i, const Real& rhs);

   /// changes left- and right-hand side vectors
   void changeRangeReal(const VectorReal& lhs, const VectorReal& rhs);

   /// changes left- and right-hand side of row \p i
   void changeRangeReal(int i, const Real& lhs, const Real& rhs);

   /// replaces column \p i with \p lpcol
   void changeColReal(int i, const LPColReal& lpcol);

   /// changes vector of lower bounds to \p lower
   void changeLowerReal(const VectorReal& lower);

   /// changes lower bound of column i to \p lower
   void changeLowerReal(int i, const Real& lower);

   /// changes vector of upper bounds to \p upper
   void changeUpperReal(const VectorReal& upper);

   /// changes \p i 'th upper bound to \p upper
   void changeUpperReal(int i, const Real& upper);

   /// changes vectors of column bounds to \p lower and \p upper
   void changeBoundsReal(const VectorReal& lower, const VectorReal& upper);

   /// changes bounds of column \p i to \p lower and \p upper
   void changeBoundsReal(int i, const Real& lower, const Real& upper);

   /// changes objective function vector to \p obj
   void changeObjReal(const VectorReal& obj);

   /// changes objective coefficient of column i to \p obj
   void changeObjReal(int i, const Real& obj);

   /// changes matrix entry in row \p i and column \p j to \p val
   void changeElementReal(int i, int j, const Real& val);

   /// removes row \p i
   void removeRowReal(int i);

   /// removes all rows with an index \p i such that \p perm[i] < 0; upon completion, \p perm[i] >= 0 indicates the
   /// new index where row \p i has been moved to; note that \p perm must point to an array of size at least
   /// #numRowsReal()
   void removeRowsReal(int perm[]);

   /// remove all rows with indices in array \p idx of size \p n; an array \p perm of size #numRowsReal() may be passed
   /// as buffer memory
   void removeRowsReal(int idx[], int n, int perm[] = 0);

   /// removes rows \p start to \p end including both; an array \p perm of size #numRowsReal() may be passed as buffer
   /// memory
   void removeRowRangeReal(int start, int end, int perm[] = 0);

   /// removes column i
   void removeColReal(int i);

   /// removes all columns with an index \p i such that \p perm[i] < 0; upon completion, \p perm[i] >= 0 indicates the
   /// new index where column \p i has been moved to; note that \p perm must point to an array of size at least
   /// #numColsReal()
   void removeColsReal(int perm[]);

   /// remove all columns with indices in array \p idx of size \p n; an array \p perm of size #numColsReal() may be
   /// passed as buffer memory
   void removeColsReal(int idx[], int n, int perm[] = 0);

   /// removes columns \p start to \p end including both; an array \p perm of size #numColsReal() may be passed as
   /// buffer memory
   void removeColRangeReal(int start, int end, int perm[] = 0);

   /// clears the LP
   void clearLPReal();

   /// synchronizes real LP with rational LP, i.e., copies (rounded) rational LP into real LP, if sync mode is manual
   void syncLPReal();

   //@}


   //**@name Modification of the rational LP */
   //@{

   /// adds a single row
   void addRowRational(const LPRowRational& lprow);

#ifdef SOPLEX_WITH_GMP
   /// adds a single row
   void addRowRational(const mpq_t* lhs, const mpq_t* rowValues, const int* rowIndices, const int rowSize, const mpq_t* rhs);

   /// adds a set of rows
   void addRowsRational(const mpq_t* lhs, const mpq_t* rowValues, const int* rowIndices, const int* rowStarts, const int* rowLengths, const int numRows, const int numValues, const mpq_t* rhs);
#endif

   /// adds multiple rows
   void addRowsRational(const LPRowSetRational& lprowset);

   /// adds a single column
   void addColRational(const LPColRational& lpcol);

#ifdef SOPLEX_WITH_GMP
   /// adds a single column
   void addColRational(const mpq_t* obj, const mpq_t* lower, const mpq_t* colValues, const int* colIndices, const int colSize, const mpq_t* upper);

   /// adds a set of columns
   void addColsRational(const mpq_t* obj, const mpq_t* lower, const mpq_t* colValues, const int* colIndices, const int* colStarts, const int* colLengths, const int numCols, const int numValues, const mpq_t* upper);
#endif

   /// adds multiple columns
   void addColsRational(const LPColSetRational& lpcolset);

   /// replaces row \p i with \p lprow
   void changeRowRational(int i, const LPRowRational& lprow);

   /// changes left-hand side vector for constraints to \p lhs
   void changeLhsRational(const VectorRational& lhs);

   /// changes left-hand side of row \p i to \p lhs
   void changeLhsRational(int i, const Rational& lhs);

#ifdef SOPLEX_WITH_GMP
   /// changes left-hand side of row \p i to \p lhs
   void changeLhsRational(int i, const mpq_t* lhs);
#endif

   /// changes right-hand side vector to \p rhs
   void changeRhsRational(const VectorRational& rhs);

#ifdef SOPLEX_WITH_GMP
   /// changes right-hand side vector to \p rhs
   void changeRhsRational(const mpq_t* rhs, int rhsSize);
#endif

   /// changes right-hand side of row \p i to \p rhs
   void changeRhsRational(int i, const Rational& rhs);

   /// changes left- and right-hand side vectors
   void changeRangeRational(const VectorRational& lhs, const VectorRational& rhs);

   /// changes left- and right-hand side of row \p i
   void changeRangeRational(int i, const Rational& lhs, const Rational& rhs);

#ifdef SOPLEX_WITH_GMP
   /// changes left- and right-hand side of row \p i
   void changeRangeRational(int i, const mpq_t* lhs, const mpq_t* rhs);
#endif

   /// replaces column \p i with \p lpcol
   void changeColRational(int i, const LPColRational& lpcol);

   /// changes vector of lower bounds to \p lower
   void changeLowerRational(const VectorRational& lower);

   /// changes lower bound of column i to \p lower
   void changeLowerRational(int i, const Rational& lower);

#ifdef SOPLEX_WITH_GMP
   /// changes lower bound of column i to \p lower
   void changeLowerRational(int i, const mpq_t* lower);
#endif

   /// changes vector of upper bounds to \p upper
   void changeUpperRational(const VectorRational& upper);

   /// changes \p i 'th upper bound to \p upper
   void changeUpperRational(int i, const Rational& upper);

#ifdef SOPLEX_WITH_GMP
   /// changes upper bound of column i to \p upper
   void changeUpperRational(int i, const mpq_t* upper);
#endif

   /// changes vectors of column bounds to \p lower and \p upper
   void changeBoundsRational(const VectorRational& lower, const VectorRational& upper);

   /// changes bounds of column \p i to \p lower and \p upper
   void changeBoundsRational(int i, const Rational& lower, const Rational& upper);

#ifdef SOPLEX_WITH_GMP
   /// changes bounds of column \p i to \p lower and \p upper
   void changeBoundsRational(int i, const mpq_t* lower, const mpq_t* upper);
#endif

   /// changes objective function vector to \p obj
   void changeObjRational(const VectorRational& obj);

   /// changes objective coefficient of column i to \p obj
   void changeObjRational(int i, const Rational& obj);

#ifdef SOPLEX_WITH_GMP
   /// changes objective coefficient of column i to \p obj
   void changeObjRational(int i, const mpq_t* obj);
#endif

   /// changes matrix entry in row \p i and column \p j to \p val
   void changeElementRational(int i, int j, const Rational& val);

#ifdef SOPLEX_WITH_GMP
   /// changes matrix entry in row \p i and column \p j to \p val
   void changeElementRational(int i, int j, const mpq_t* val);
#endif

   /// removes row \p i
   void removeRowRational(int i);

   /// removes all rows with an index \p i such that \p perm[i] < 0; upon completion, \p perm[i] >= 0 indicates the new
   /// index where row \p i has been moved to; note that \p perm must point to an array of size at least
   /// #numRowsRational()
   void removeRowsRational(int perm[]);

   /// remove all rows with indices in array \p idx of size \p n; an array \p perm of size #numRowsRational() may be
   /// passed as buffer memory
   void removeRowsRational(int idx[], int n, int perm[] = 0);

   /// removes rows \p start to \p end including both; an array \p perm of size #numRowsRational() may be passed as
   /// buffer memory
   void removeRowRangeRational(int start, int end, int perm[] = 0);

   /// removes column i
   void removeColRational(int i);

   /// removes all columns with an index \p i such that \p perm[i] < 0; upon completion, \p perm[i] >= 0 indicates the
   /// new index where column \p i has been moved to; note that \p perm must point to an array of size at least
   /// #numColsRational()
   void removeColsRational(int perm[]);

   /// remove all columns with indices in array \p idx of size \p n; an array \p perm of size #numColsRational() may be
   /// passed as buffer memory
   void removeColsRational(int idx[], int n, int perm[] = 0);

   /// removes columns \p start to \p end including both; an array \p perm of size #numColsRational() may be passed as
   /// buffer memory
   void removeColRangeRational(int start, int end, int perm[] = 0);

   /// clears the LP
   void clearLPRational();

   /// synchronizes rational LP with real LP, i.e., copies real LP to rational LP, if sync mode is manual
   void syncLPRational();

   //@}


   //**@name Solving and general solution query */
   //@{

   /// solves the LP
   SPxSolver::Status solve();

   /// returns the current solver status
   SPxSolver::Status status() const;

   /// is a primal feasible solution available?
   bool hasPrimal() const;

   /// is a primal unbounded ray available?
   bool hasPrimalRay() const;

   /// is a dual feasible solution available?
   bool hasDual() const;

   /// is Farkas proof of infeasibility available?
   bool hasDualFarkas() const;

   //@}


   //**@name Query for the real solution data */
   //@{

   /// returns the objective value if a primal solution is available
   Real objValueReal();

   /// gets the primal solution vector if available; returns true on success
   bool getPrimalReal(VectorReal& vector);

   /// gets the vector of slack values if available; returns true on success
   bool getSlacksReal(VectorReal& vector);

   /// gets the primal ray if available; returns true on success
   bool getPrimalRayReal(VectorReal& vector);

   /// gets the dual solution vector if available; returns true on success
   bool getDualReal(VectorReal& vector);

   /// gets the vector of reduced cost values if available; returns true on success
   bool getRedCostReal(VectorReal& vector);

   /// gets the Farkas proof if available; returns true on success
   bool getDualFarkasReal(VectorReal& vector);

   /// gets violation of bounds; returns true on success
   bool getBoundViolationReal(Real& maxviol, Real& sumviol);

   /// gets violation of constraints; returns true on success
   bool getRowViolationReal(Real& maxviol, Real& sumviol);

   /// gets violation of reduced costs; returns true on success
   bool getRedCostViolationReal(Real& maxviol, Real& sumviol);

   /// gets violation of dual multipliers; returns true on success
   bool getDualViolationReal(Real& maxviol, Real& sumviol);

   //@}


   //**@name Query for the rational solution data */
   //@{

   /// returns the objective value if a primal solution is available
   Rational objValueRational();

   /// gets the primal solution vector if available; returns true on success
   bool getPrimalRational(VectorRational& vector);

   /// gets the vector of slack values if available; returns true on success
   bool getSlacksRational(VectorRational& vector);

   /// gets the primal ray if LP is unbounded; returns true on success
   bool getPrimalRayRational(VectorRational& vector);

   /// gets the dual solution vector if available; returns true on success
   bool getDualRational(VectorRational& vector);

   /// gets the vector of reduced cost values if available; returns true on success
   bool getRedCostRational(VectorRational& vector);

   /// gets the Farkas proof if LP is infeasible; returns true on success
   bool getDualFarkasRational(VectorRational& vector);

   /// gets violation of bounds; returns true on success
   bool getBoundViolationRational(Rational& maxviol, Rational& sumviol);

   /// gets violation of constraints; returns true on success
   bool getRowViolationRational(Rational& maxviol, Rational& sumviol);

   /// gets violation of reduced costs; returns true on success
   bool getRedCostViolationRational(Rational& maxviol, Rational& sumviol);

   /// gets violation of dual multipliers; returns true on success
   bool getDualViolationRational(Rational& maxviol, Rational& sumviol);

#ifdef SOPLEX_WITH_GMP
   /// gets the primal solution vector if available; returns true on success
   bool getPrimalRational(mpq_t* vector, const int size);

   /// gets the vector of slack values if available; returns true on success
   bool getSlacksRational(mpq_t* vector, const int size);

   /// gets the primal ray if LP is unbounded; returns true on success
   bool getPrimalRayRational(mpq_t* vector, const int size);

   /// gets the dual solution vector if available; returns true on success
   bool getDualRational(mpq_t* vector, const int size);

   /// gets the vector of reduced cost values if available; returns true on success
   bool getRedCostRational(mpq_t* vector, const int size);

   /// gets the Farkas proof if LP is infeasible; returns true on success
   bool getDualFarkasRational(mpq_t* vector, const int size);
#endif

   /// get size of primal solution
   int totalSizePrimalRational(const int base = 2);

   /// get size of dual solution
   int totalSizeDualRational(const int base = 2);

   /// get size of least common multiple of denominators in primal solution
   int dlcmSizePrimalRational(const int base = 2);

   /// get size of least common multiple of denominators in dual solution
   int dlcmSizeDualRational(const int base = 2);

   /// get size of largest denominator in primal solution
   int dmaxSizePrimalRational(const int base = 2);

   /// get size of largest denominator in dual solution
   int dmaxSizeDualRational(const int base = 2);

   //@}


   //**@name Access and modification of basis information */
   //@{

   /// is an advanced starting basis available?
   bool hasBasis() const;

   /// returns the current basis status
   SPxBasis::SPxStatus basisStatus() const;

   /// returns basis status for a single row
   SPxSolver::VarStatus basisRowStatus(int row) const;

   /// returns basis status for a single column
   SPxSolver::VarStatus basisColStatus(int col) const;

   /// gets current basis via arrays of statuses
   void getBasis(SPxSolver::VarStatus rows[], SPxSolver::VarStatus cols[]) const;

   /// gets the indices of the basic columns and rows; basic column n gives value n, basic row m gives value -1-m
   void getBasisInd(int* bind) const;

   /// computes an estimated condition number for the current basis matrix using the power method; returns true on success
   bool getEstimatedCondition(Real& condition);

   /// computes the exact condition number for the current basis matrix using the power method; returns true on success
   bool getExactCondition(Real& condition);

   /// computes row r of basis inverse; returns true on success
   bool getBasisInverseRowReal(int r, Real* coef, int* inds = NULL, int* ninds = NULL);

   /// computes column c of basis inverse; returns true on success
   bool getBasisInverseColReal(int c, Real* coef, int* inds = NULL, int* ninds = NULL);

   /// computes dense solution of basis matrix B * sol = rhs; returns true on success
   bool getBasisInverseTimesVecReal(Real* rhs, Real* sol);

   /// sets starting basis via arrays of statuses
   void setBasis(SPxSolver::VarStatus rows[], SPxSolver::VarStatus cols[]);

   /// clears starting basis
   void clearBasis();

   //@}


   //**@name Statistical information */
   //@{

   /// number of iterations since last call to solve
   int numIterations() const;

   /// time spent in last call to solve
   Real solveTime() const;

   /// statistical information in form of a string
   std::string statisticString() const;

   /// name of starter
   const char* getStarterName();

   /// name of simplifier
   const char* getSimplifierName();

   /// name of scaling method
   const char* getScalerName();

   /// name of currently loaded pricer
   const char* getPricerName();

   /// name of currently loaded ratiotester
   const char* getRatiotesterName();

   //@}


   //**@name File I/O */
   //@{

   /// reads LP file in LP or MPS format according to READMODE parameter; gets row names, column names, and
   /// integer variables if desired; returns true on success
   bool readFile(const char* filename, NameSet* rowNames = 0, NameSet* colNames = 0, DIdxSet* intVars = 0);

   /// writes real LP to file; LP or MPS format is chosen from the extension in \p filename; if \p rowNames and \p
   /// colNames are \c NULL, default names are used; if \p intVars is not \c NULL, the variables contained in it are
   /// marked as integer; returns true on success
   bool writeFileReal(const char* filename, const NameSet* rowNames = 0, const NameSet* colNames = 0, const DIdxSet* intvars = 0) const;

   /// writes rational LP to file; LP or MPS format is chosen from the extension in \p filename; if \p rowNames and \p
   /// colNames are \c NULL, default names are used; if \p intVars is not \c NULL, the variables contained in it are
   /// marked as integer; returns true on success
   bool writeFileRational(const char* filename, const NameSet* rowNames = 0, const NameSet* colNames = 0, const DIdxSet* intvars = 0) const;

   /// reads basis information from \p filename and returns true on success; if \p rowNames and \p colNames are \c NULL,
   /// default names are assumed; returns true on success
   bool readBasisFile(const char* filename, const NameSet* rowNames = 0, const NameSet* colNames = 0);

   /// writes basis information to \p filename; if \p rowNames and \p colNames are \c NULL, default names are used;
   /// returns true on success
   bool writeBasisFile(const char* filename, const NameSet* rowNames = 0, const NameSet* colNames = 0, const bool cpxFormat = false) const;

   /// writes internal LP, basis information, and parameter settings; if \p rowNames and \p colNames are \c NULL,
   /// default names are used
   void writeStateReal(const char* filename, const NameSet* rowNames = 0, const NameSet* colNames = 0, const bool cpxFormat = false) const;

   /// writes internal LP, basis information, and parameter settings; if \p rowNames and \p colNames are \c NULL,
   /// default names are used
   void writeStateRational(const char* filename, const NameSet* rowNames = 0, const NameSet* colNames = 0, const bool cpxFormat = false) const;

   //@}


   //**@name Parameters */
   //@{

   /// boolean parameters
   typedef enum
   {
      /// should lifting be used to reduce range of nonzero matrix coefficients?
      LIFTING = 0,

      /// should LP be transformed to equality form before a rational solve?
      EQTRANS = 1,

      /// should dual infeasibility be tested in order to try to return a dual solution even if primal infeasible?
      TESTDUALINF = 2,

      /// should a rational factorization be performed after iterative refinement?
      RATFAC = 3,

      /// should cycling solutions be accepted during iterative refinement?
      ACCEPTCYCLING = 4,

      /// apply rational reconstruction after each iterative refinement?
      RATREC = 5,

      /// round scaling factors for iterative refinement to powers of two?
      POWERSCALING = 6,

      /// continue iterative refinement with exact basic solution if not optimal?
      RATFACJUMP = 7,

      /// should feasibility be tested with relaxed bounds and sides?
      FEASRELAX = 8,

      /// use bound flipping also for row representation?
      ROWBOUNDFLIPS = 9,

      /// number of boolean parameters
      BOOLPARAM_COUNT = 10
   } BoolParam;

   /// integer parameters
   typedef enum
   {
      /// objective sense
      OBJSENSE = 0,

      /// type of computational form, i.e., column or row representation
      REPRESENTATION = 1,

      /// type of algorithm, i.e., primal or dual
      ALGORITHM = 2,

      /// type of LU update
      FACTOR_UPDATE_TYPE = 3,

      /// maximum number of updates without fresh factorization
      FACTOR_UPDATE_MAX = 4,

      /// iteration limit (-1 if unlimited)
      ITERLIMIT = 5,

      /// refinement limit (-1 if unlimited)
      REFLIMIT = 6,

      /// stalling refinement limit (-1 if unlimited)
      STALLREFLIMIT = 7,

      /// display frequency
      DISPLAYFREQ = 8,

      /// verbosity level
      VERBOSITY = 9,

      /// type of simplifier
      SIMPLIFIER = 10,

      /// type of scaler
      SCALER = 11,

      /// type of starter used to create crash basis
      STARTER = 12,

      /// type of pricer
      PRICER = 13,

      /// type of ratio test
      RATIOTESTER = 14,

      /// mode for synchronizing real and rational LP
      SYNCMODE = 15,

      /// mode for reading LP files
      READMODE = 16,

      /// mode for iterative refinement strategy
      SOLVEMODE = 17,

      /// mode for a posteriori feasibility checks
      CHECKMODE = 18,

      /// type of timer
      TIMER = 19,

      /// mode for hyper sparse pricing
      HYPER_PRICING = 20,

      /// minimum number of stalling refinements since last pivot to trigger rational factorization
      RATFAC_MINSTALLS = 21,

      /// number of integer parameters
      INTPARAM_COUNT = 22
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
      /// automatic choice according to number of rows and columns
      REPRESENTATION_AUTO = 0,

      /// column representation Ax - s = 0, lower <= x <= upper, lhs <= s <= rhs
      REPRESENTATION_COLUMN = 1,

      /// row representation (lower,lhs) <= (x,Ax) <= (upper,rhs)
      REPRESENTATION_ROW = 2
   };

   /// values for parameter ALGORITHM
   enum
   {
      /// primal simplex algorithm, i.e., entering for column and leaving for row representation
      ALGORITHM_PRIMAL = 0,

      /// dual simplex algorithm, i.e., leaving for column and entering for row representation
      ALGORITHM_DUAL = 1
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
      SCALER_GEO8 = 4
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
      PRICER_STEEP = 5
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

   /// values for parameter SYNCMODE
   enum
   {
      /// store only real LP
      SYNCMODE_ONLYREAL = 0,

      /// automatic sync of real and rational LP
      SYNCMODE_AUTO = 1,

      /// user sync of real and rational LP
      SYNCMODE_MANUAL = 2
   };

   /// values for parameter READMODE
   enum
   {
      /// standard floating-point parsing
      READMODE_REAL = 0,

      /// rational parsing
      READMODE_RATIONAL = 1
   };

   /// values for parameter SOLVEMODE
   enum
   {
      /// apply standard floating-point algorithm
      SOLVEMODE_REAL = 0,

      /// decide depending on tolerances whether to apply iterative refinement
      SOLVEMODE_AUTO = 1,

      /// force iterative refinement
      SOLVEMODE_RATIONAL = 2
   };

   /// values for parameter CHECKMODE
   enum
   {
      /// floating-point check
      CHECKMODE_REAL = 0,

      /// decide according to READMODE
      CHECKMODE_AUTO = 1,

      /// rational check
      CHECKMODE_RATIONAL = 2
   };

   /// values for parameter TIMER
   enum
   {
      /// disable timing
      TIMER_OFF = 0,

      /// cpu or user time
      TIMER_CPU = 1,

      /// wallclock time
      TIMER_WALLCLOCK = 2
   };

   /// values for parameter HYPER_PRICING
   enum
   {
      /// never
      HYPER_PRICING_OFF = 0,

      /// decide according to problem size
      HYPER_PRICING_AUTO = 1,

      /// always
      HYPER_PRICING_ON = 2
   };

   /// real parameters
   typedef enum
   {
      /// primal feasibility tolerance
      FEASTOL = 0,

      /// dual feasibility tolerance
      OPTTOL = 1,

      /// general zero tolerance
      EPSILON_ZERO = 2,

      /// zero tolerance used in factorization
      EPSILON_FACTORIZATION = 3,

      /// zero tolerance used in update of the factorization
      EPSILON_UPDATE = 4,

      /// pivot zero tolerance used in factorization
      EPSILON_PIVOT = 5,

      /// infinity threshold
      INFTY = 6,

      /// time limit in seconds (INFTY if unlimited)
      TIMELIMIT = 7,

      /// lower limit on objective value
      OBJLIMIT_LOWER = 8,

      /// upper limit on objective value
      OBJLIMIT_UPPER = 9,

      /// working tolerance for feasibility in floating-point solver during iterative refinement
      FPFEASTOL = 10,

      /// working tolerance for optimality in floating-point solver during iterative refinement
      FPOPTTOL = 11,

      /// maximum increase of scaling factors between refinements
      MAXSCALEINCR = 12,

      /// lower threshold in lifting (nonzero matrix coefficients with smaller absolute value will be reformulated)
      LIFTMINVAL = 13,

      /// upper threshold in lifting (nonzero matrix coefficients with larger absolute value will be reformulated)
      LIFTMAXVAL = 14,

      /// sparse pricing threshold (\#violations < dimension * SPARSITY_THRESHOLD activates sparse pricing)
      SPARSITY_THRESHOLD = 15,

      /// threshold on number of rows vs. number of columns for switching from column to row representations in auto mode
      REPRESENTATION_SWITCH = 16,

      /// geometric frequency at which to apply rational reconstruction
      RATREC_FREQ = 17,

      /// minimal reduction (sum of removed rows/cols) to continue simplification
      MINRED = 18,

      /// number of real parameters
      REALPARAM_COUNT = 19
   } RealParam;

#ifdef SOPLEX_WITH_RATIONALPARAM
   /// rational parameters
   typedef enum
   {
      /// number of rational parameters
      RATIONALPARAM_COUNT = 0
   } RationalParam;
#endif

   /// class of parameter settings
   class Settings;

   mutable SPxOut spxout;

   /// returns boolean parameter value
   bool boolParam(const BoolParam param) const;

   /// returns integer parameter value
   int intParam(const IntParam param) const;

   /// returns real parameter value
   Real realParam(const RealParam param) const;

#ifdef SOPLEX_WITH_RATIONALPARAM
   /// returns rational parameter value
   Rational rationalParam(const RationalParam param) const;
#endif

   /// returns current parameter settings
   const Settings& settings() const;

   /// sets boolean parameter value; returns true on success
   bool setBoolParam(const BoolParam param, const bool value, const bool quiet = false, const bool init = false);

   /// sets integer parameter value; returns true on success
   bool setIntParam(const IntParam param, const int value, const bool quiet = false, const bool init = false);

   /// sets real parameter value; returns true on success
   bool setRealParam(const RealParam param, const Real value, const bool quiet = false, const bool init = false);

#ifdef SOPLEX_WITH_RATIONALPARAM
   /// sets rational parameter value; returns true on success
   bool setRationalParam(const RationalParam param, const Rational value, const bool quiet = false, const bool init = false);
#endif

   /// sets parameter settings; returns true on success
   bool setSettings(const Settings& newSettings, const bool quiet = false, const bool init = false);

   /// print non-default parameter values
   void printUserSettings();

   /// writes settings file; returns true on success
   bool saveSettingsFile(const char* filename, const bool onlyChanged = false) const;

   /// reads settings file; returns true on success
   bool loadSettingsFile(const char* filename);

   /// parses one setting string and returns true on success; note that string is modified
   bool parseSettingsString(char* line);

   //@}


   //**@name Statistics */
   //@{

   /// prints solution statistics
   void printSolutionStatistics(std::ostream& os);

   /// prints statistics on solving process
   void printSolvingStatistics(std::ostream& os);

   /// prints short statistics
   void printShortStatistics(std::ostream& os);

   /// prints complete statistics
   void printStatistics(std::ostream& os);

   /// prints status
   void printStatus(std::ostream& os, SPxSolver::Status status);

   //@}


   //**@name Miscellaneous */
   //@{

   /// prints version and compilation options
   void printVersion() const;

   /// checks if real LP and rational LP are in sync; dimensions will always be compared,
   /// vector and matrix values only if the respective parameter is set to true.
   /// If quiet is set to true the function will only display which vectors are different.
   bool areLPsInSync(const bool checkVecVals = true, const bool checkMatVals = false, const bool quiet = false) const;

   //@}

private:

   //**@name Statistics on solving process */
   //@{

   /// class of statistics
   class Statistics;

   /// statistics since last call to solveReal() or solveRational()
   Statistics* _statistics;

   //@}


   //**@name Parameter settings */
   //@{

   Settings* _currentSettings;

   Rational _rationalPosInfty;
   Rational _rationalNegInfty;
   Rational _rationalFeastol;
   Rational _rationalOpttol;
   Rational _rationalMaxscaleincr;

   //@}


   //**@name Data for the real LP */
   //@{

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
   SPxAutoPR _pricerAuto;
   SPxDantzigPR _pricerDantzig;
   SPxParMultPR _pricerParMult;
   SPxDevexPR _pricerDevex;
   SPxSteepPR _pricerQuickSteep;
   SPxSteepExPR _pricerSteep;
   SPxDefaultRT _ratiotesterTextbook;
   SPxHarrisRT _ratiotesterHarris;
   SPxFastRT _ratiotesterFast;
   SPxBoundFlippingRT _ratiotesterBoundFlipping;

   SPxLPReal* _realLP;
   SPxSimplifier* _simplifier;
   SPxScaler* _scaler;
   SPxStarter* _starter;

   bool _isRealLPLoaded;

   DVectorReal _manualLower;
   DVectorReal _manualUpper;
   DVectorReal _manualLhs;
   DVectorReal _manualRhs;
   DVectorReal _manualObj;
   SPxLPReal _manualRealLP;

   //@}


   //**@name Data for the rational LP */
   //@{

   SPxLPRational* _rationalLP;

   LPColSetRational _slackCols;
   DVectorRational _unboundedLower;
   DVectorRational _unboundedUpper;
   DVectorRational _unboundedLhs;
   DVectorRational _unboundedRhs;
   DSVectorRational _tauColVector;
   DVectorRational _feasObj;
   DVectorRational _feasLhs;
   DVectorRational _feasRhs;
   DVectorRational _feasLower;
   DVectorRational _feasUpper;
   DVectorRational _modLower;
   DVectorRational _modUpper;
   DVectorRational _modLhs;
   DVectorRational _modRhs;
   DVectorRational _modObj;
   DSVectorRational _primalDualDiff;
   DataArray< SPxSolver::VarStatus > _storedBasisStatusRows;
   DataArray< SPxSolver::VarStatus > _storedBasisStatusCols;
   DataArray< UnitVectorRational* > _unitMatrixRational;
   bool _storedBasis;
   int _beforeLiftRows;
   int _beforeLiftCols;

   /// type of bounds and sides
   typedef enum
   {
      /// both bounds are infinite
      RANGETYPE_FREE = 0,

      /// lower bound is finite, upper bound is infinite
      RANGETYPE_LOWER = 1,

      /// upper bound is finite, lower bound is infinite
      RANGETYPE_UPPER = 2,

      /// lower and upper bound finite, but different
      RANGETYPE_BOXED = 3,

      /// lower bound equals upper bound
      RANGETYPE_FIXED = 4
   } RangeType;

   DataArray< RangeType > _colTypes;
   DataArray< RangeType > _rowTypes;

   //@}


   //**@name Solution data */
   //@{

   SPxSolver::Status _status;
   int _lastSolveMode;

   DataArray< SPxSolver::VarStatus > _basisStatusRows;
   DataArray< SPxSolver::VarStatus > _basisStatusCols;

   SolReal _solReal;
   SolRational _solRational;
   SolRational _workSol;

   bool _hasBasis;
   bool _hasSolReal;
   bool _hasSolRational;

   //@}


   //**@name Constant helper methods */
   //@{

   /// extends sparse vector to hold newmax entries if and only if it holds no more free entries
   void _ensureDSVectorRationalMemory(DSVectorRational& vec, const int newmax) const;

   /// creates a permutation for removing rows/columns from an array of indices
   void _idxToPerm(int* idx, int idxSize, int* perm, int permSize) const;

   /// creates a permutation for removing rows/columns from a range of indices
   void _rangeToPerm(int start, int end, int* perm, int permSize) const;

   /// checks consistency
   bool _isConsistent() const;

   /// should solving process be stopped?
   bool _isSolveStopped() const;

   /// determines RangeType from real bounds
   RangeType _rangeTypeReal(const Real& lower, const Real& upper) const;

   /// determines RangeType from rational bounds
   RangeType _rangeTypeRational(const Rational& lower, const Rational& upper) const;

   /// switches RANGETYPE_LOWER to RANGETYPE_UPPER and vice versa
   RangeType _switchRangeType(const RangeType& rangeType) const;

   /// checks whether RangeType corresponds to finite lower bound
   bool _lowerFinite(const RangeType& rangeType) const;

   /// checks whether RangeType corresponds to finite upper bound
   bool _upperFinite(const RangeType& rangeType) const;

   //@}


   //**@name Non-constant helper methods */
   //@{

   /// adds a single row to the real LP and adjusts basis
   void _addRowReal(const LPRowReal& lprow);

   /// adds a single row to the real LP and adjusts basis
   void _addRowReal(Real lhs, const SVectorReal& lprow, Real rhs);

   /// adds multiple rows to the real LP and adjusts basis
   void _addRowsReal(const LPRowSetReal& lprowset);

   /// adds a single column to the real LP and adjusts basis
   void _addColReal(const LPColReal& lpcol);

   /// adds a single column to the real LP and adjusts basis
   void _addColReal(Real obj, Real lower, const SVectorReal& lpcol, Real upper);

   /// adds multiple columns to the real LP and adjusts basis
   void _addColsReal(const LPColSetReal& lpcolset);

   /// replaces row \p i with \p lprow and adjusts basis
   void _changeRowReal(int i, const LPRowReal& lprow);

   /// changes left-hand side vector for constraints to \p lhs and adjusts basis
   void _changeLhsReal(const VectorReal& lhs);

   /// changes left-hand side of row \p i to \p lhs and adjusts basis
   void _changeLhsReal(int i, const Real& lhs);

   /// changes right-hand side vector to \p rhs and adjusts basis
   void _changeRhsReal(const VectorReal& rhs);

   /// changes right-hand side of row \p i to \p rhs and adjusts basis
   void _changeRhsReal(int i, const Real& rhs);

   /// changes left- and right-hand side vectors and adjusts basis
   void _changeRangeReal(const VectorReal& lhs, const VectorReal& rhs);

   /// changes left- and right-hand side of row \p i and adjusts basis
   void _changeRangeReal(int i, const Real& lhs, const Real& rhs);

   /// replaces column \p i with \p lpcol and adjusts basis
   void _changeColReal(int i, const LPColReal& lpcol);

   /// changes vector of lower bounds to \p lower and adjusts basis
   void _changeLowerReal(const VectorReal& lower);

   /// changes lower bound of column i to \p lower and adjusts basis
   void _changeLowerReal(int i, const Real& lower);

   /// changes vector of upper bounds to \p upper and adjusts basis
   void _changeUpperReal(const VectorReal& upper);

   /// changes \p i 'th upper bound to \p upper and adjusts basis
   void _changeUpperReal(int i, const Real& upper);

   /// changes vectors of column bounds to \p lower and \p upper and adjusts basis
   void _changeBoundsReal(const VectorReal& lower, const VectorReal& upper);

   /// changes bounds of column \p i to \p lower and \p upper and adjusts basis
   void _changeBoundsReal(int i, const Real& lower, const Real& upper);

   /// changes matrix entry in row \p i and column \p j to \p val and adjusts basis
   void _changeElementReal(int i, int j, const Real& val);

   /// removes row \p i and adjusts basis
   void _removeRowReal(int i);

   /// removes all rows with an index \p i such that \p perm[i] < 0; upon completion, \p perm[i] >= 0 indicates the
   /// new index where row \p i has been moved to; note that \p perm must point to an array of size at least
   /// #numRowsReal()
   void _removeRowsReal(int perm[]);

   /// remove all rows with indices in array \p idx of size \p n; an array \p perm of size #numRowsReal() may be passed
   /// as buffer memory
   void _removeRowsReal(int idx[], int n, int perm[]);

   /// removes rows \p start to \p end including both; an array \p perm of size #numRowsReal() may be passed as buffer
   /// memory
   void _removeRowRangeReal(int start, int end, int perm[]);

   /// removes column i
   void _removeColReal(int i);

   /// removes all columns with an index \p i such that \p perm[i] < 0; upon completion, \p perm[i] >= 0 indicates the
   /// new index where column \p i has been moved to; note that \p perm must point to an array of size at least
   /// #numColsReal()
   void _removeColsReal(int perm[]);

   /// remove all columns with indices in array \p idx of size \p n; an array \p perm of size #numColsReal() may be
   /// passed as buffer memory
   void _removeColsReal(int idx[], int n, int perm[]);

   /// removes columns \p start to \p end including both; an array \p perm of size #numColsReal() may be passed as
   /// buffer memory
   void _removeColRangeReal(int start, int end, int perm[]);

   /// invalidates solution
   void _invalidateSolution();

   /// enables simplifier and scaler according to current parameters
   void _enableSimplifierAndScaler();

   /// disables simplifier and scaler
   void _disableSimplifierAndScaler();

   /// ensures that the rational LP is available; performs no sync
   void _ensureRationalLP();

   /// ensures that the real LP and the basis are loaded in the solver; performs no sync
   void _ensureRealLPLoaded();

   /// call floating-point solver and update statistics on iterations etc.
   void _solveRealLPAndRecordStatistics();

   /// reads real LP in LP or MPS format from file and returns true on success; gets row names, column names, and
   /// integer variables if desired
   bool _readFileReal(const char* filename, NameSet* rowNames = 0, NameSet* colNames = 0, DIdxSet* intVars = 0);

   /// reads rational LP in LP or MPS format from file and returns true on success; gets row names, column names, and
   /// integer variables if desired
   bool _readFileRational(const char* filename, NameSet* rowNames = 0, NameSet* colNames = 0, DIdxSet* intVars = 0);

   /// recomputes range types from scratch using real LP
   void _recomputeRangeTypesReal();

   /// recomputes range types from scratch using rational LP
   void _recomputeRangeTypesRational();

   /// synchronizes real LP with rational LP, i.e., copies (rounded) rational LP into real LP, without looking at the sync mode
   void _syncLPReal(bool time = true);

   /// synchronizes rational LP with real LP, i.e., copies real LP to rational LP, without looking at the sync mode
   void _syncLPRational(bool time = true);

   /// synchronizes real solution with rational solution, i.e., copies real solution to rational solution
   void _syncRealSolution();

   /// synchronizes rational solution with real solution, i.e., copies (rounded) rational solution to real solution
   void _syncRationalSolution();

   /// returns pointer to a constant unit vector available until destruction of the SoPlex class
   const UnitVectorRational* _unitVectorRational(const int i);

   /// parses one line in a settings file and returns true on success; note that the string is modified
   bool _parseSettingsLine(char* line, const int lineNumber);

   //@}


   //**@name Private solving methods implemented in solverational.cpp */
   //@{

   /// solves rational LP
   void _solveRational();

   /// solves current problem with iterative refinement and recovery mechanism
   void _performOptIRStable(SolRational& sol, bool acceptUnbounded, bool acceptInfeasible, int minRounds, bool& primalFeasible, bool& dualFeasible, bool& infeasible, bool& unbounded, bool& stopped, bool& error);

   /// performs iterative refinement on the auxiliary problem for testing unboundedness
   void _performUnboundedIRStable(SolRational& sol, bool& hasUnboundedRay, bool& stopped, bool& error);

   /// performs iterative refinement on the auxiliary problem for testing feasibility
   void _performFeasIRStable(SolRational& sol, bool& withDualFarkas, bool& stopped, bool& error);

   /// reduces matrix coefficient in absolute value by the lifting procedure of Thiele et al. 2013
   void _lift();

   /// undoes lifting
   void _project(SolRational& sol);

   /// store basis
   void _storeBasis();

   /// restore basis
   void _restoreBasis();

   /// stores objective, bounds, and sides of real LP
   void _storeLPReal();

   /// restores objective, bounds, and sides of real LP
   void _restoreLPReal();

   /// introduces slack variables to transform inequality constraints into equations for both rational and real LP,
   /// which should be in sync
   void _transformEquality();

   /// undoes transformation to equality form
   void _untransformEquality(SolRational& sol);

   /// transforms LP to unboundedness problem by moving the objective function to the constraints, changing right-hand
   /// side and bounds to zero, and adding an auxiliary variable for the decrease in the objective function
   void _transformUnbounded();

   /// undoes transformation to unboundedness problem
   void _untransformUnbounded(SolRational& sol, bool unbounded);

   /// transforms LP to feasibility problem by removing the objective function, shifting variables, and homogenizing the
   /// right-hand side
   void _transformFeasibility();

   /// undoes transformation to feasibility problem
   void _untransformFeasibility(SolRational& sol, bool infeasible);

   /** computes radius of infeasibility box implied by an approximate Farkas' proof

     Given constraints of the form \f$ lhs <= Ax <= rhs \f$, a farkas proof y should satisfy \f$ y^T A = 0 \f$ and
     \f$ y_+^T lhs - y_-^T rhs > 0 \f$, where \f$ y_+, y_- \f$ denote the positive and negative parts of \f$ y \f$.
     If \f$ y \f$ is approximate, it may not satisfy \f$ y^T A = 0 \f$ exactly, but the proof is still valid as long
     as the following holds for all potentially feasible \f$ x \f$:

     \f[
         y^T Ax < (y_+^T lhs - y_-^T rhs)              (*)
     \f]

     we may therefore calculate \f$ y^T A \f$ and \f$ y_+^T lhs - y_-^T rhs \f$ exactly and check if the upper and lower
     bounds on \f$ x \f$ imply that all feasible \f$ x \f$ satisfy (*), and if not then compute bounds on \f$ x \f$ to
     guarantee (*).  The simplest way to do this is to compute

     \f[
     B = (y_+^T lhs - y_-^T rhs) / \sum_i(|(y^T A)_i|)
     \f]

     noting that if every component of \f$ x \f$ has \f$ |x_i| < B \f$, then (*) holds.

     \f$ B \f$ can be increased by iteratively including variable bounds smaller than \f$ B \f$.  The speed of this
     method can be further improved by using interval arithmetic for all computations.  For related information see
     Sec. 4 of Neumaier and Shcherbina, Mathematical Programming A, 2004.

     Set transformed to true if this method is called after _transformFeasibility().
     */
   void _computeInfeasBox(SolRational& sol, bool transformed);

   /// solves real LP during iterative refinement
   SPxSolver::Status _solveRealForRational(bool fromscratch, VectorReal& primal, VectorReal& dual,
                                           DataArray< SPxSolver::VarStatus >& basisStatusRows,
                                           DataArray< SPxSolver::VarStatus >& basisStatusCols, bool& returnedBasis);

   /// solves real LP with recovery mechanism
   SPxSolver::Status _solveRealStable(bool acceptUnbounded, bool acceptInfeasible, VectorReal& primal, VectorReal& dual,
                                      DataArray< SPxSolver::VarStatus >& basisStatusRows,
                                      DataArray< SPxSolver::VarStatus >& basisStatusCols, bool& returnedBasis, const bool forceNoSimplifier = false);

   /// factorizes rational basis matrix in column representation
   void _factorizeColumnRational(SolRational& sol, DataArray< SPxSolver::VarStatus >& basisStatusRows, DataArray< SPxSolver::VarStatus >& basisStatusCols, bool& stopped, bool& error, bool& optimal);

   /// attempts rational reconstruction of primal-dual solution
   bool _reconstructSolutionRational(SolRational& sol, DataArray< SPxSolver::VarStatus >& basisStatusRows, DataArray< SPxSolver::VarStatus >& basisStatusCols, bool& stopped, bool& error, const Rational& denomBoundSquared);
   //@}


   //**@name Private solving methods implemented in solvereal.cpp */
   //@{

   /// solves real LP
   void _solveReal();

   /// checks result of the solving process and solves again without preprocessing if necessary
   void _evaluateSolutionReal(SPxSimplifier::Result simplificationStatus);

   /// solves real LP with/without preprocessing
   void _preprocessAndSolveReal(bool applyPreprocessing);

   /// loads original problem into solver and solves again after it has been solved to optimality with preprocessing
   void _resolveWithoutPreprocessing(SPxSimplifier::Result simplificationStatus);

   /// stores solution of the real LP; before calling this, the real LP must be loaded in the solver and solved (again)
   void _storeSolutionReal();

   //@}
};
}
#else
#include "soplexlegacy.h"

namespace soplex
{
   typedef SoPlexLegacy SoPlex;
}
#endif
#endif // _SOPLEX_H_
