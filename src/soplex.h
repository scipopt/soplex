/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the class library                   */
/*       SoPlex --- the Sequential object-oriented simPlex.                  */
/*                                                                           */
/*  Copyright (c) 1996-2023 Zuse Institute Berlin (ZIB)                      */
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

/**@file  soplex.h
 * @brief Preconfigured SoPlex LP solver
 */

#ifndef _SOPLEX_H_
#define _SOPLEX_H_

#include <string.h>

#include "soplex/spxgithash.h"
#include "soplex/spxdefines.h"
#include "soplex/basevectors.h"
#include "soplex/spxsolver.h"
#include "soplex/slufactor.h"
#include "soplex/slufactor_rational.h"

///@todo try to move to cpp file by forward declaration
#include "soplex/spxsimplifier.h"
#include "soplex/spxmainsm.h"

#include "soplex/spxscaler.h"
#include "soplex/spxequilisc.h"
#include "soplex/spxleastsqsc.h"
#include "soplex/spxgeometsc.h"

#include "soplex/spxstarter.h"
#include "soplex/spxweightst.h"
#include "soplex/spxsumst.h"
#include "soplex/spxvectorst.h"

#include "soplex/spxpricer.h"
#include "soplex/spxautopr.h"
#include "soplex/spxdantzigpr.h"
#include "soplex/spxparmultpr.h"
#include "soplex/spxdevexpr.h"
#include "soplex/spxsteeppr.h"
#include "soplex/spxsteepexpr.h"
#include "soplex/spxhybridpr.h"

#include "soplex/spxratiotester.h"
#include "soplex/spxdefaultrt.h"
#include "soplex/spxharrisrt.h"
#include "soplex/spxfastrt.h"
#include "soplex/spxboundflippingrt.h"

#include "soplex/solbase.h"
#include "soplex/sol.h"

#include "soplex/spxlpbase.h"

#include "soplex/spxpapilo.h"

#ifdef SOPLEX_WITH_GMP
#include <gmp.h>
#endif

#ifdef SOPLEX_WITH_BOOST
#ifdef SOPLEX_WITH_MPFR
// For multiple precision
#include <boost/multiprecision/mpfr.hpp>
#endif
#ifdef SOPLEX_WITH_CPPMPF
#include <boost/multiprecision/cpp_dec_float.hpp>
#endif

// An alias for boost multiprecision
namespace mpf = boost::multiprecision;
#endif

#define SOPLEX_DEFAULT_RANDOM_SEED   0   // used to suppress output when the seed was not changed

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
///@todo integrate rational SPxSolver and distinguish between original and transformed rational LP
///@todo rational scalers
///@todo rational simplifier

namespace soplex
{

/**@class SoPlex
 * @brief   Preconfigured SoPlex LP-solver.
 * @ingroup Algo
 */
template <class R>
class SoPlexBase
{
public:

   ///@name Construction and destruction
   ///@{

   /// default constructor
   SoPlexBase();

   /// assignment operator
   SoPlexBase<R>& operator=(const SoPlexBase<R>& rhs);

   /// copy constructor
   SoPlexBase(const SoPlexBase<R>& rhs);

   /// destructor
   virtual ~SoPlexBase();

   ///@}


   ///@name Access to the real LP
   ///@{

   /// returns number of rows
   int numRows() const;
   int numRowsReal() const;     /* For SCIP compatibility */
   int numRowsRational() const;

   /// Templated function that
   /// returns number of columns
   int numCols() const;
   int numColsReal() const;     /* For SCIP compatibility */
   int numColsRational() const;

   /// returns number of nonzeros
   int numNonzeros() const;

   int numNonzerosRational() const;

   /// returns smallest non-zero element in absolute value
   R minAbsNonzeroReal() const;

   /// returns biggest non-zero element in absolute value
   R maxAbsNonzeroReal() const;

   /// returns (unscaled) coefficient
   R coefReal(int row, int col) const;

   /// returns vector of row \p i, ignoring scaling
   const SVectorBase<R>& rowVectorRealInternal(int i) const;

   /// gets vector of row \p i
   void getRowVectorReal(int i, DSVectorBase<R>& row) const;

   /// returns right-hand side vector, ignoring scaling
   const VectorBase<R>& rhsRealInternal() const;

   /// gets right-hand side vector
   void getRhsReal(VectorBase<R>& rhs) const;

   /// returns right-hand side of row \p i
   R rhsReal(int i) const;

   /// returns left-hand side vector, ignoring scaling
   const VectorBase<R>& lhsRealInternal() const;

   /// gets left-hand side vector
   void getLhsReal(VectorBase<R>& lhs) const;

   /// returns left-hand side of row \p i
   R lhsReal(int i) const;

   /// returns inequality type of row \p i
   typename LPRowBase<R>::Type rowTypeReal(int i) const;

   /// returns vector of col \p i, ignoring scaling
   const SVectorBase<R>& colVectorRealInternal(int i) const;

   /// gets vector of col \p i
   void getColVectorReal(int i, DSVectorBase<R>& col) const;

   /// returns upper bound vector
   const VectorBase<R>& upperRealInternal() const;

   /// returns upper bound of column \p i
   R upperReal(int i) const;

   /// gets upper bound vector
   void getUpperReal(VectorBase<R>& upper) const;

   /// returns lower bound vector
   const VectorBase<R>& lowerRealInternal() const;

   /// returns lower bound of column \p i
   R lowerReal(int i) const;

   /// gets lower bound vector
   void getLowerReal(VectorBase<R>& lower) const;

   /// gets objective function vector
   void getObjReal(VectorBase<R>& obj) const;

   /// returns objective value of column \p i
   R objReal(int i) const;

   /// returns objective function vector after transformation to a maximization problem; since this is how it is stored
   /// internally, this is generally faster
   const VectorBase<R>& maxObjRealInternal() const;

   /// returns objective value of column \p i after transformation to a maximization problem; since this is how it is
   /// stored internally, this is generally faster
   R maxObjReal(int i) const;

   /// gets number of available dual norms
   void getNdualNorms(int& nnormsRow, int& nnormsCol) const;

   /// gets steepest edge norms and returns false if they are not available
   bool getDualNorms(int& nnormsRow, int& nnormsCol, R* norms) const;

   /// sets steepest edge norms and returns false if that's not possible
   bool setDualNorms(int nnormsRow, int nnormsCol, R* norms);

   /// pass integrality information about the variables to the solver
   void setIntegralityInformation(int ncols, int* intInfo);

   ///@}


   ///@name Access to the rational LP
   ///@{

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

   ///@}


   ///@name Modification of the real LP
   ///@{

   /// adds a single row
   void addRowReal(const LPRowBase<R>& lprow);

   /// adds multiple rows
   void addRowsReal(const LPRowSetBase<R>& lprowset);

   /// adds a single column
   void addColReal(const LPColBase<R>& lpcol);

   /// adds multiple columns
   void addColsReal(const LPColSetBase<R>& lpcolset);

   /// replaces row \p i with \p lprow
   void changeRowReal(int i, const LPRowBase<R>& lprow);

   /// changes left-hand side vector for constraints to \p lhs
   void changeLhsReal(const VectorBase<R>& lhs);

   /// changes left-hand side of row \p i to \p lhs
   void changeLhsReal(int i, const R& lhs);

   /// changes right-hand side vector to \p rhs
   void changeRhsReal(const VectorBase<R>& rhs);

   /// changes right-hand side of row \p i to \p rhs
   void changeRhsReal(int i, const R& rhs);

   /// changes left- and right-hand side vectors
   void changeRangeReal(const VectorBase<R>& lhs, const VectorBase<R>& rhs);

   /// changes left- and right-hand side of row \p i
   void changeRangeReal(int i, const R& lhs, const R& rhs);

   /// replaces column \p i with \p lpcol
   void changeColReal(int i, const LPColReal& lpcol);

   /// changes vector of lower bounds to \p lower
   void changeLowerReal(const VectorBase<R>& lower);

   /// changes lower bound of column i to \p lower
   void changeLowerReal(int i, const R& lower);

   /// changes vector of upper bounds to \p upper
   void changeUpperReal(const VectorBase<R>& upper);

   /// changes \p i 'th upper bound to \p upper
   void changeUpperReal(int i, const R& upper);

   /// changes vectors of column bounds to \p lower and \p upper
   void changeBoundsReal(const VectorBase<R>& lower, const VectorBase<R>& upper);

   /// changes bounds of column \p i to \p lower and \p upper
   void changeBoundsReal(int i, const R& lower, const R& upper);

   /// changes objective function vector to \p obj
   void changeObjReal(const VectorBase<R>& obj);

   /// changes objective coefficient of column i to \p obj
   void changeObjReal(int i, const R& obj);

   /// changes matrix entry in row \p i and column \p j to \p val
   void changeElementReal(int i, int j, const R& val);

   /// removes row \p i
   void removeRowReal(int i);

   /// removes all rows with an index \p i such that \p perm[i] < 0; upon completion, \p perm[i] >= 0 indicates the
   /// new index where row \p i has been moved to; note that \p perm must point to an array of size at least
   /// #numRows()
   void removeRowsReal(int perm[]);

   /// remove all rows with indices in array \p idx of size \p n; an array \p perm of size #numRows() may be passed
   /// as buffer memory
   void removeRowsReal(int idx[], int n, int perm[] = 0);

   /// removes rows \p start to \p end including both; an array \p perm of size #numRows() may be passed as buffer
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

   ///@}


   ///@name Modification of the rational LP
   ///@{

   /// adds a single row
   void addRowRational(const LPRowRational& lprow);

#ifdef SOPLEX_WITH_GMP
   /// adds a single row (GMP only method)
   void addRowRational(const mpq_t* lhs, const mpq_t* rowValues, const int* rowIndices,
                       const int rowSize, const mpq_t* rhs);

   /// adds a set of rows (GMP only method)
   void addRowsRational(const mpq_t* lhs, const mpq_t* rowValues, const int* rowIndices,
                        const int* rowStarts, const int* rowLengths, const int numRows, const int numValues,
                        const mpq_t* rhs);
#endif

   /// adds multiple rows
   void addRowsRational(const LPRowSetRational& lprowset);

   /// adds a single column
   void addColRational(const LPColRational& lpcol);

#ifdef SOPLEX_WITH_GMP
   /// adds a single column (GMP only method)
   void addColRational(const mpq_t* obj, const mpq_t* lower, const mpq_t* colValues,
                       const int* colIndices, const int colSize, const mpq_t* upper);

   /// adds a set of columns (GMP only method)
   void addColsRational(const mpq_t* obj, const mpq_t* lower, const mpq_t* colValues,
                        const int* colIndices, const int* colStarts, const int* colLengths, const int numCols,
                        const int numValues, const mpq_t* upper);
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
   /// changes left-hand side of row \p i to \p lhs (GMP only method)
   void changeLhsRational(int i, const mpq_t* lhs);
#endif

   /// changes right-hand side vector to \p rhs
   void changeRhsRational(const VectorRational& rhs);

#ifdef SOPLEX_WITH_GMP
   /// changes right-hand side vector to \p rhs (GMP only method)
   void changeRhsRational(const mpq_t* rhs, int rhsSize);
#endif

   /// changes right-hand side of row \p i to \p rhs
   void changeRhsRational(int i, const Rational& rhs);

   /// changes left- and right-hand side vectors
   void changeRangeRational(const VectorRational& lhs, const VectorRational& rhs);

   /// changes left- and right-hand side of row \p i
   void changeRangeRational(int i, const Rational& lhs, const Rational& rhs);

#ifdef SOPLEX_WITH_GMP
   /// changes left- and right-hand side of row \p i (GMP only method)
   void changeRangeRational(int i, const mpq_t* lhs, const mpq_t* rhs);
#endif

   /// replaces column \p i with \p lpcol
   void changeColRational(int i, const LPColRational& lpcol);

   /// changes vector of lower bounds to \p lower
   void changeLowerRational(const VectorRational& lower);

   /// changes lower bound of column i to \p lower
   void changeLowerRational(int i, const Rational& lower);

#ifdef SOPLEX_WITH_GMP
   /// changes lower bound of column i to \p lower (GMP only method)
   void changeLowerRational(int i, const mpq_t* lower);
#endif

   /// changes vector of upper bounds to \p upper
   void changeUpperRational(const VectorRational& upper);

   /// changes \p i 'th upper bound to \p upper
   void changeUpperRational(int i, const Rational& upper);

#ifdef SOPLEX_WITH_GMP
   /// changes upper bound of column i to \p upper (GMP only method)
   void changeUpperRational(int i, const mpq_t* upper);
#endif

   /// changes vectors of column bounds to \p lower and \p upper
   void changeBoundsRational(const VectorRational& lower, const VectorRational& upper);

   /// changes bounds of column \p i to \p lower and \p upper
   void changeBoundsRational(int i, const Rational& lower, const Rational& upper);

#ifdef SOPLEX_WITH_GMP
   /// changes bounds of column \p i to \p lower and \p upper (GMP only method)
   void changeBoundsRational(int i, const mpq_t* lower, const mpq_t* upper);
#endif

   /// changes objective function vector to \p obj
   void changeObjRational(const VectorRational& obj);

   /// changes objective coefficient of column i to \p obj
   void changeObjRational(int i, const Rational& obj);

#ifdef SOPLEX_WITH_GMP
   /// changes objective coefficient of column i to \p obj (GMP only method)
   void changeObjRational(int i, const mpq_t* obj);
#endif

   /// changes matrix entry in row \p i and column \p j to \p val
   void changeElementRational(int i, int j, const Rational& val);

#ifdef SOPLEX_WITH_GMP
   /// changes matrix entry in row \p i and column \p j to \p val (GMP only method)
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

   ///@}

   ///@name Solving and general solution query
   ///@{

   /// optimize the given LP
   typename SPxSolverBase<R>::Status optimize(volatile bool* interrupt = NULL);

   // old name for backwards compatibility
   typename SPxSolverBase<R>::Status solve(volatile bool* interrupt = NULL)
   {
      return optimize(interrupt);
   }

   /// returns the current solver status
   typename SPxSolverBase<R>::Status status() const;

   /// is stored primal solution feasible?
   bool isPrimalFeasible() const;

   /// is a solution available (not neccessarily feasible)?
   bool hasSol() const;

   /// deprecated: use #hasSol() instead
   bool hasPrimal() const
   {
      return hasSol();
   }

   /// deprecated: use #hasSol() instead
   bool hasDual() const
   {
      return hasSol();
   }

   /// is a primal unbounded ray available?
   bool hasPrimalRay() const;

   /// is stored dual solution feasible?
   bool isDualFeasible() const;

   /// is Farkas proof of infeasibility available?
   bool hasDualFarkas() const;

   /// sets the status to OPTIMAL in case the LP has been solved with unscaled violations
   bool ignoreUnscaledViolations()
   {
      if(_status == SPxSolverBase<R>::OPTIMAL_UNSCALED_VIOLATIONS)
      {
         _status = SPxSolverBase<R>::OPTIMAL;
         return true;
      }
      else
         return false;
   }
   ///@}


   ///@name Query for the real solution data
   ///@{

   /// returns the objective value if a primal solution is available
   R objValueReal();

   /// gets the primal solution vector if available; returns true on success
   bool getPrimal(VectorBase<R>& vector);
   bool getPrimalReal(R* p_vector, int size);      // For SCIP compatibility
   bool getPrimalRational(VectorRational& vector);

   /// Get the value for the \p i 'th row given the current primal solution
   bool getRowPrimalValue(R& value, int i);
   bool getRowPrimalValueRational(Rational& value, int i);

   /// Get the value for the row with indices in \p indices given the current primal solution
   bool getRowsPrimalValue(const std::set<int>& indices, VectorBase<R>& vector);
   bool getRowsPrimalValueReal(const std::set<int>& indices, R* p_vector, int i);
   bool getRowsPrimalValueRational(const std::set<int>& indices, VectorRational& vector);

   /// Get the value for all rows given the current primal solution
   bool getRowsPrimalValue(VectorBase<R>& vector);
   bool getRowsPrimalValueReal(R* p_vector, int i);
   bool getRowsPrimalValueRational(VectorRational& vector);

   /// gets the vector of slack values if available; returns true on success
   bool getSlacksReal(VectorBase<R>& vector);
   bool getSlacksReal(R* p_vector, int dim);

   /// gets the primal ray if available; returns true on success
   bool getPrimalRay(VectorBase<R>& vector);
   bool getPrimalRayReal(R* vector, int dim); /* For SCIP compatibility */
   bool getPrimalRayRational(VectorRational& vector);

   /// gets the dual solution vector if available; returns true on success
   bool getDual(VectorBase<R>& vector);
   bool getDualReal(R* p_vector, int dim); /* For SCIP compatibility */
   bool getDualRational(VectorRational& vector);

   /// gets the vector of reduced cost values if available; returns true on success
   bool getRedCost(VectorBase<R>& vector);
   bool getRedCostReal(R* vector, int dim); /* For SCIP compatibility */
   bool getRedCostRational(VectorRational& vector);

   /// gets the Farkas proof if available; returns true on success
   bool getDualFarkas(VectorBase<R>& vector);
   bool getDualFarkasReal(R* vector, int dim);
   bool getDualFarkasRational(VectorRational& vector);

   /// gets violation of bounds; returns true on success
   bool getBoundViolation(R& maxviol, R& sumviol);
   bool getBoundViolationRational(Rational& maxviol, Rational& sumviol);

   /// gets violation of constraints; returns true on success
   bool getRowViolation(R& maxviol, R& sumviol);
   bool getRowViolationRational(Rational& maxviol, Rational& sumviol);

   /// gets violation of reduced costs; returns true on success
   bool getRedCostViolation(R& maxviol, R& sumviol);
   bool getRedCostViolationRational(Rational& maxviol, Rational& sumviol);

   /// gets violation of dual multipliers; returns true on success
   bool getDualViolation(R& maxviol, R& sumviol);
   bool getDualViolationRational(Rational& maxviol, Rational& sumviol);

   ///@}


   ///@name Query for the rational solution data
   ///@{

   /// returns the objective value if a primal solution is available
   Rational objValueRational();

   /// gets the vector of slack values if available; returns true on success
   bool getSlacksRational(VectorRational& vector);

#ifdef SOPLEX_WITH_GMP
   /// gets the primal solution vector if available; returns true on success (GMP only method)
   bool getPrimalRational(mpq_t* vector, const int size);

   /// gets the vector of slack values if available; returns true on success (GMP only method)
   bool getSlacksRational(mpq_t* vector, const int size);

   /// gets the primal ray if LP is unbounded; returns true on success (GMP only method)
   bool getPrimalRayRational(mpq_t* vector, const int size);

   /// gets the dual solution vector if available; returns true on success (GMP only method)
   bool getDualRational(mpq_t* vector, const int size);

   /// gets the vector of reduced cost values if available; returns true on success (GMP only method)
   bool getRedCostRational(mpq_t* vector, const int size);

   /// gets the Farkas proof if LP is infeasible; returns true on success (GMP only method)
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

   ///@}


   ///@name Access and modification of basis information
   ///@{

   /// is an advanced starting basis available?
   bool hasBasis() const;

   /// returns the current basis status
   typename SPxBasisBase<R>::SPxStatus basisStatus() const;

   /// returns basis status for a single row
   typename  SPxSolverBase<R>::VarStatus basisRowStatus(int row) const;

   /// returns basis status for a single column
   typename SPxSolverBase<R>::VarStatus basisColStatus(int col) const;

   /// gets current basis via arrays of statuses
   void getBasis(typename SPxSolverBase<R>::VarStatus rows[],
                 typename SPxSolverBase<R>::VarStatus cols[]) const;

   /// gets the indices of the basic columns and rows; basic column n gives value n, basic row m gives value -1-m
   void getBasisInd(int* bind) const;

   /** compute one of several matrix metrics based on the diagonal of the LU factorization
    *  type = 0: max/min ratio
    *  type = 1: trace of U (sum of diagonal elements)
    *  type = 2: determinant (product of diagonal elements)
    */
   bool getBasisMetric(R& metric, int type = 0);

   /// computes an estimated condition number for the current basis matrix using the power method; returns true on success
   bool getEstimatedCondition(R& condition);

   /// computes the exact condition number for the current basis matrix using the power method; returns true on success
   bool getExactCondition(R& condition);

   /// computes row \p r of basis inverse; returns true on success
   /// @param r which row of the basis inverse is computed
   /// @param coef values of result vector (not packed but scattered)
   /// @param inds indices of result vector (NULL if not to be used)
   /// @param ninds number of nonzeros in result vector
   /// @param unscale determines whether the result should be unscaled according to the original LP data
   bool getBasisInverseRowReal(int r, R* coef, int* inds = NULL, int* ninds = NULL,
                               bool unscale = true);

   /// computes column \p c of basis inverse; returns true on success
   /// @param c which column of the basis inverse is computed
   /// @param coef values of result vector (not packed but scattered)
   /// @param inds indices of result vector (NULL if not to be used)
   /// @param ninds number of nonzeros in result vector
   /// @param unscale determines whether the result should be unscaled according to the original LP data
   bool getBasisInverseColReal(int c, R* coef, int* inds = NULL, int* ninds = NULL,
                               bool unscale = true);

   /// computes dense solution of basis matrix B * \p sol = \p rhs; returns true on success
   bool getBasisInverseTimesVecReal(R* rhs, R* sol, bool unscale = true);

   /// multiply with basis matrix; B * \p vec (inplace)
   /// @param vec (dense) vector to be multiplied with
   /// @param unscale determines whether the result should be unscaled according to the original LP data
   bool multBasis(R* vec, bool unscale = true);

   /// multiply with transpose of basis matrix; \p vec * B^T (inplace)
   /// @param vec (dense) vector to be multiplied with
   /// @param unscale determines whether the result should be unscaled according to the original LP data
   bool multBasisTranspose(R* vec, bool unscale = true);

   /// compute rational basis inverse; returns true on success
   bool computeBasisInverseRational();

   /// gets an array of indices for the columns of the rational basis matrix; bind[i] >= 0 means that the i-th column of
   /// the basis matrix contains variable bind[i]; bind[i] < 0 means that the i-th column of the basis matrix contains
   /// the slack variable for row -bind[i]-1; performs rational factorization if not available; returns true on success
   bool getBasisIndRational(DataArray<int>& bind);

   /// computes row r of basis inverse; performs rational factorization if not available; returns true on success
   bool getBasisInverseRowRational(const int r, SSVectorRational& vec);

   /// computes column c of basis inverse; performs rational factorization if not available; returns true on success
   bool getBasisInverseColRational(const int c, SSVectorRational& vec);

   /// computes solution of basis matrix B * sol = rhs; performs rational factorization if not available; returns true
   /// on success
   bool getBasisInverseTimesVecRational(const SVectorRational& rhs, SSVectorRational& sol);

   /// sets starting basis via arrays of statuses
   void setBasis(const typename SPxSolverBase<R>::VarStatus rows[],
                 const typename SPxSolverBase<R>::VarStatus cols[]);

   /// clears starting basis
   void clearBasis();

   ///@}


   ///@name Statistical information
   ///@{

   /// number of iterations since last call to solve
   int numIterations() const;

   /// number of precision boosts since last call to solve
   int numPrecisionBoosts() const;

   /// number of iterations in higher precision since last call to solve
   int numIterationsBoosted() const;

   /// time spen in higher precision since last call to solve
   Real precisionBoostTime() const;

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

   ///@}


   ///@name File I/O
   ///@{

   /// reads LP file in LP or MPS format according to READMODE parameter; gets row names, column names, and
   /// integer variables if desired; returns true on success
   bool readFile(const char* filename, NameSet* rowNames = 0, NameSet* colNames = 0,
                 DIdxSet* intVars = 0);

   /// Templated write function
   /// Real
   /// writes real LP to file; LP or MPS format is chosen from the extension in \p filename; if \p rowNames and \p
   /// colNames are \c NULL, default names are used; if \p intVars is not \c NULL, the variables contained in it are
   /// marked as integer; returns true on success
   /// Rational
   /// writes rational LP to file; LP or MPS format is chosen from the extension in \p filename; if \p rowNames and \p
   /// colNames are \c NULL, default names are used; if \p intVars is not \c NULL, the variables contained in it are
   /// marked as integer; returns true on success
   bool writeFile(const char* filename, const NameSet* rowNames = 0, const NameSet* colNames = 0,
                  const DIdxSet* intvars = 0, const bool unscale = true) const;

   bool writeFileRational(const char* filename, const NameSet* rowNames = 0,
                          const NameSet* colNames = 0, const DIdxSet* intvars = 0) const;

   /* For SCIP compatibility */
   bool writeFileReal(const char* filename, const NameSet* rowNames = 0, const NameSet* colNames = 0,
                      const DIdxSet* intvars = 0, const bool unscale = true) const;

   /// writes the dual of the real LP to file; LP or MPS format is chosen from the extension in \p filename;
   /// if \p rowNames and \p colNames are \c NULL, default names are used; if \p intVars is not \c NULL,
   /// the variables contained in it are marked as integer; returns true on success
   bool writeDualFileReal(const char* filename, const NameSet* rowNames = 0,
                          const NameSet* colNames = 0, const DIdxSet* intvars = 0) const;

   /// reads basis information from \p filename and returns true on success; if \p rowNames and \p colNames are \c NULL,
   /// default names are assumed; returns true on success
   bool readBasisFile(const char* filename, const NameSet* rowNames = 0, const NameSet* colNames = 0);

   /// writes basis information to \p filename; if \p rowNames and \p colNames are \c NULL, default names are used;
   /// returns true on success
   bool writeBasisFile(const char* filename, const NameSet* rowNames = 0, const NameSet* colNames = 0,
                       const bool cpxFormat = false) const;

   /// writes internal LP, basis information, and parameter settings; if \p rowNames and \p colNames are \c NULL,
   /// default names are used
   void writeStateReal(const char* filename, const NameSet* rowNames = 0, const NameSet* colNames = 0,
                       const bool cpxFormat = false) const;

   /// writes internal LP, basis information, and parameter settings; if \p rowNames and \p colNames are \c NULL,
   /// default names are used
   void writeStateRational(const char* filename, const NameSet* rowNames = 0,
                           const NameSet* colNames = 0, const bool cpxFormat = false) const;

   ///@}


   ///@name Parameters
   ///@{

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

      /// use bound flipping also for row representation?
      ROWBOUNDFLIPS = 8,

      /// use persistent scaling?
      PERSISTENTSCALING = 9,

      /// perturb the entire problem or only the relevant bounds of s single pivot?
      FULLPERTURBATION = 10,

      /// re-optimize the original problem to get a proof (ray) of infeasibility/unboundedness?
      ENSURERAY = 11,

      /// try to enforce that the optimal solution is a basic solution
      FORCEBASIC = 12,

      // enable presolver SingletonCols in PaPILO?
      SIMPLIFIER_SINGLETONCOLS = 13,

      // enable presolver ConstraintPropagation in PaPILO?
      SIMPLIFIER_CONSTRAINTPROPAGATION = 14,

      // enable presolver ParallelRowDetection in PaPILO?
      SIMPLIFIER_PARALLELROWDETECTION = 15,

      // enable presolver ParallelColDetection in PaPILO?
      SIMPLIFIER_PARALLELCOLDETECTION = 16,

      // enable presolver SingletonStuffing in PaPILO?
      SIMPLIFIER_SINGLETONSTUFFING = 17,

      // enable presolver DualFix in PaPILO?
      SIMPLIFIER_DUALFIX = 18,

      // enable presolver FixContinuous in PaPILO?
      SIMPLIFIER_FIXCONTINUOUS = 19,

      // enable presolver DominatedCols in PaPILO?
      SIMPLIFIER_DOMINATEDCOLS = 20,

      // enable iterative refinement ?
      ITERATIVE_REFINEMENT = 21,

      /// adapt tolerances to the multiprecision used
      ADAPT_TOLS_TO_MULTIPRECISION = 22,

      /// enable precision boosting ?
      PRECISION_BOOSTING = 23,

      /// boosted solver start from last basis
      BOOSTED_WARM_START = 24,

      /// try different settings when solve fails
      RECOVERY_MECHANISM = 25,

      /// store advanced and stable basis met before each simplex iteration, to better warm start
      STORE_BASIS_BEFORE_SIMPLEX_PIVOT = 26,

      /// number of boolean parameters
      BOOLPARAM_COUNT = 27
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

      /// maximum number of conjugate gradient iterations in least square scaling
      LEASTSQ_MAXROUNDS = 22,

      /// mode for solution polishing
      SOLUTION_POLISHING = 23,

      /// print condition number during the solve
      PRINTBASISMETRIC = 24,

      /// type of timer for statistics
      STATTIMER = 25,

      // maximum number of digits for the multiprecision type
      MULTIPRECISION_LIMIT = 26,

      ///@todo precision-boosting find better parameter name
      /// after how many simplex pivots do we store the advanced and stable basis, 1 = every iterations
      STORE_BASIS_SIMPLEX_FREQ = 27,

      /// number of integer parameters
      INTPARAM_COUNT = 28
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
      /// disabling presolving
      SIMPLIFIER_OFF = 0,

      /// using internal presolving methods
      SIMPLIFIER_INTERNAL = 3,

      /// using the presolve lib papilo
      SIMPLIFIER_PAPILO = 2,

      /// SoPlex chooses automatically (currently always "internal")
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

      /// least square scaling
      SCALER_LEASTSQ = 5,

      /// geometric mean scaling (max 8 rounds) followed by equilibrium scaling (rows and columns)
      SCALER_GEOEQUI = 6
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

   /// values for parameter SOLUTION_POLISHING
   enum
   {
      /// no solution polishing
      POLISHING_OFF = 0,

      /// maximize number of basic slack variables, i.e. more variables on bounds
      POLISHING_INTEGRALITY = 1,

      /// minimize number of basic slack variables, i.e. more variables between bounds
      POLISHING_FRACTIONALITY = 2
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

      /// refactor threshold for nonzeros in last factorized basis matrix compared to updated basis matrix
      REFAC_BASIS_NNZ = 19,

      /// refactor threshold for fill-in in current factor update compared to fill-in in last factorization
      REFAC_UPDATE_FILL = 20,

      /// refactor threshold for memory growth in factorization since last refactorization
      REFAC_MEM_FACTOR = 21,

      /// accuracy of conjugate gradient method in least squares scaling (higher value leads to more iterations)
      LEASTSQ_ACRCY = 22,

      /// objective offset
      OBJ_OFFSET = 23,

      /// minimal Markowitz threshold to control sparsity/stability in LU factorization
      MIN_MARKOWITZ = 24,

      /// minimal modification threshold to apply presolve reductions
      SIMPLIFIER_MODIFYROWFAC = 25,

      /// factor by which the precision of the floating-point solver is multiplied
      PRECISION_BOOSTING_FACTOR = 26,

      /// number of real parameters
      REALPARAM_COUNT = 27
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
   class Settings
   {
   public:
      static struct BoolParam
      {
         /// constructor
         BoolParam();
         /// array of names for boolean parameters
         std::string name[SoPlexBase<R>::BOOLPARAM_COUNT];
         /// array of descriptions for boolean parameters
         std::string description[SoPlexBase<R>::BOOLPARAM_COUNT];
         /// array of default values for boolean parameters
         bool defaultValue[SoPlexBase<R>::BOOLPARAM_COUNT];
      } boolParam;

      static struct IntParam
      {
         /// constructor
         IntParam();
         /// array of names for integer parameters
         std::string name[SoPlexBase<R>::INTPARAM_COUNT];
         /// array of descriptions for integer parameters
         std::string description[SoPlexBase<R>::INTPARAM_COUNT];
         /// array of default values for integer parameters
         int defaultValue[SoPlexBase<R>::INTPARAM_COUNT];
         /// array of lower bounds for int parameter values
         int lower[SoPlexBase<R>::INTPARAM_COUNT];
         /// array of upper bounds for int parameter values
         int upper[SoPlexBase<R>::INTPARAM_COUNT];
      } intParam;

      static struct RealParam
      {
         /// constructor
         RealParam();
         /// array of names for real parameters
         std::string name[SoPlexBase<R>::REALPARAM_COUNT];
         /// array of descriptions for real parameters
         std::string description[SoPlexBase<R>::REALPARAM_COUNT];
         /// array of default values for real parameters
         Real defaultValue[SoPlexBase<R>::REALPARAM_COUNT];
         /// array of lower bounds for real parameter values
         Real lower[SoPlexBase<R>::REALPARAM_COUNT];
         /// array of upper bounds for real parameter values
         Real upper[SoPlexBase<R>::REALPARAM_COUNT];
      } realParam;

#ifdef SOPLEX_WITH_RATIONALPARAM
      static struct RationalParam
      {
         /// constructor
         RationalParam();
         /// array of names for rational parameters
         std::string name[SoPlexBase<R>::RATIONALPARAM_COUNT];
         /// array of descriptions for rational parameters
         std::string description[SoPlexBase<R>::RATIONALPARAM_COUNT];
         /// array of default values for rational parameters
         Rational defaultValue[SoPlexBase<R>::RATIONALPARAM_COUNT];
         /// array of lower bounds for rational parameter values
         Rational lower[SoPlexBase<R>::RATIONALPARAM_COUNT];
         /// array of upper bounds for rational parameter values
         Rational upper[SoPlexBase<R>::RATIONALPARAM_COUNT];
      } rationalParam;
#endif

      /// array of current boolean parameter values
      bool _boolParamValues[SoPlexBase<R>::BOOLPARAM_COUNT];

      /// array of current integer parameter values
      int _intParamValues[SoPlexBase<R>::INTPARAM_COUNT];

      /// array of current real parameter values
      Real _realParamValues[SoPlexBase<R>::REALPARAM_COUNT];

#ifdef SOPLEX_WITH_RATIONALPARAM
      /// array of current rational parameter values
      Rational _rationalParamValues[SoPlexBase<R>::RATIONALPARAM_COUNT];
#endif

      /// default constructor initializing default settings
      Settings();

      /// copy constructor
      Settings(const Settings& settings);

      /// assignment operator
      Settings& operator=(const Settings& settings);
   };

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

   /// returns current tolerances
   const std::shared_ptr<Tolerances> tolerances() const;

   /// sets boolean parameter value; returns true on success
   bool setBoolParam(const BoolParam param, const bool value, const bool init = true);

   /// sets integer parameter value; returns true on success
   bool setIntParam(const IntParam param, const int value, const bool init = true);

   /// sets real parameter value; returns true on success
   bool setRealParam(const RealParam param, const Real value, const bool init = true);

#ifdef SOPLEX_WITH_RATIONALPARAM
   /// sets rational parameter value; returns true on success
   bool setRationalParam(const RationalParam param, const Rational value, const bool init = true);
#endif

   /// sets parameter settings; returns true on success
   bool setSettings(const Settings& newSettings, const bool init = true);

   /// resets default parameter settings
   void resetSettings(const bool quiet = false, const bool init = true);

   /// print non-default parameter values
   void printUserSettings();

   /// writes settings file; returns true on success
   bool saveSettingsFile(const char* filename, const bool onlyChanged = false,
                         int solvemode = 1) const;

   /// reads settings file; returns true on success
   bool loadSettingsFile(const char* filename);

   /// parses one setting string and returns true on success; note that string is modified
   bool parseSettingsString(char* str);

   ///@}


   ///@name Statistics
   ///@{

   /// set statistic timers to a certain type
   void setTimings(const Timer::TYPE ttype);

   /// prints solution statistics
   void printSolutionStatistics(std::ostream& os);

   /// prints statistics on solving process
   void printSolvingStatistics(std::ostream& os);

   /// prints short statistics
   void printShortStatistics(std::ostream& os);

   /// prints complete statistics
   void printStatistics(std::ostream& os);

   /// prints status

   void printStatus(std::ostream& os, typename SPxSolverBase<R>::Status status);

   ///@}


   ///@name Miscellaneous
   ///@{

   /// prints version and compilation options
   void printVersion() const;

   /// checks if real LP and rational LP are in sync; dimensions will always be compared,
   /// vector and matrix values only if the respective parameter is set to true.
   /// If quiet is set to true the function will only display which vectors are different.
   bool areLPsInSync(const bool checkVecVals = true, const bool checkMatVals = false,
                     const bool quiet = false) const;

   /// set the random seeds of the solver instance
   void setRandomSeed(unsigned int seed);

   /// returns the current random seed of the solver instance
   unsigned int randomSeed() const;

   ///@}

private:

   ///@name Statistics on solving process
   ///@{

   /// class of statistics
   class Statistics;

   /// statistics since last call to solveReal() or solveRational()
   Statistics* _statistics;

   ///@}


   ///@name Parameter settings
   ///@{

   Settings* _currentSettings;

   std::shared_ptr<Tolerances> _tolerances;

   Rational _rationalPosInfty;
   Rational _rationalNegInfty;
   Rational _rationalFeastol;
   Rational _rationalOpttol;
   Rational _rationalMaxscaleincr;

   ///@}


   ///@name Data for the real LP
   ///@{

   SPxSolverBase<R> _solver;
   SLUFactor<R> _slufactor;
   SPxMainSM<R> _simplifierMainSM;
   Presol<R> _simplifierPaPILO;
   SPxEquiliSC<R> _scalerUniequi;
   SPxEquiliSC<R> _scalerBiequi;
   SPxGeometSC<R> _scalerGeo1;
   SPxGeometSC<R> _scalerGeo8;
   SPxGeometSC<R> _scalerGeoequi;
   SPxLeastSqSC<R> _scalerLeastsq;
   SPxWeightST<R> _starterWeight;
   SPxSumST<R> _starterSum;
   SPxVectorST<R> _starterVector;
   SPxAutoPR<R> _pricerAuto;
   SPxDantzigPR<R> _pricerDantzig;
   SPxParMultPR<R> _pricerParMult;
   SPxDevexPR<R> _pricerDevex;
   SPxSteepPR<R> _pricerQuickSteep;
   SPxSteepExPR<R> _pricerSteep;
   SPxDefaultRT<R> _ratiotesterTextbook;
   SPxHarrisRT<R> _ratiotesterHarris;
   SPxFastRT<R> _ratiotesterFast;
   SPxBoundFlippingRT<R> _ratiotesterBoundFlipping;

   SPxLPBase<R>* _realLP;
   SPxSimplifier<R>* _simplifier;
   SPxScaler<R>* _scaler;
   SPxStarter<R>* _starter;

#ifdef SOPLEX_WITH_BOOST
#ifdef SOPLEX_WITH_MPFR
   //----------------------------- BOOSTED SOLVER -----------------------------
   // multiprecision type used for the boosted solver
   using BP = number<mpfr_float_backend<0>, et_off>;
#else
#ifdef SOPLEX_WITH_GMP
   using BP = number<gmp_float<50>, et_off>;
#else
   using BP = number<cpp_dec_float<50>, et_off>;
#endif
#endif
#else
   using BP = double;
#endif

   // boosted solver object
   SPxSolverBase<BP> _boostedSolver;

   // ------------- Main attributes for precision boosting

   int _initialPrecision   = 50; // initial number of digits for multiprecision
   bool _boostingLimitReached; // true if BP::default_precision() > max authorized number of digits
   bool _switchedToBoosted; // true if _boostedSolver is used instead of _solver to cope with the numerical failure of _solver
   // this attribute remembers wether we are testing feasibility (1), unboundedness (2) or neither (0)
   // it is used when storing/loading the right basis in precision boosting.
   // example: if _certificateMode == 1, it is the basis for the feasibility LP that should be stored/loaded.
   int _certificateMode;

   // ------------- Buffers for statistics of precision boosting

   // ideally these four attributes would be local variables, however the precision boosting loop
   // wraps the solve in a way that it is complicated to declare these variables locally.
   int _lastStallPrecBoosts; // number of previous stalling precision boosts
   bool _factorSolNewBasisPrecBoost; // false if the current basis has already been factorized (no new iterations have been done)
   int _nextRatrecPrecBoost; // the iteration during or after which rational reconstruction can be performed
   // buffer storing the number of iterations before a given precision boost
   // used to detect stalling (_prevIterations < _statistics->iterations)
   int _prevIterations;

   // ------------- Tolerances Ratios

   /// ratios for computing the tolerances for precision boosting
   /// ratio denotes the proportion of precision used by the tolerance
   /// e.g. ratio = 0.65, precision = 100 digits, new tol = 10^(0.65*100)
   Real _tolPrecisionRatio = 0.65;
   Real _epsZeroPrecisionRatio = 1.0;
   Real _epsFactorPrecisionRatio = 1.25;
   Real _epsUpdatePrecisionRatio = 1.0;
   Real _epsPivotPrecisionRatio = 0.625;

   // ------------- [Boosted] SLUFactor, Pricers, RatioTesters, Scalers, Simplifiers

   SLUFactor<BP> _boostedSlufactor;

   SPxAutoPR<BP> _boostedPricerAuto;
   SPxDantzigPR<BP> _boostedPricerDantzig;
   SPxParMultPR<BP> _boostedPricerParMult;
   SPxDevexPR<BP> _boostedPricerDevex;
   SPxSteepPR<BP> _boostedPricerQuickSteep;
   SPxSteepExPR<BP> _boostedPricerSteep;

   SPxDefaultRT<BP> _boostedRatiotesterTextbook;
   SPxHarrisRT<BP> _boostedRatiotesterHarris;
   SPxFastRT<BP> _boostedRatiotesterFast;
   SPxBoundFlippingRT<BP> _boostedRatiotesterBoundFlipping;

   SPxScaler<BP>* _boostedScaler;
   SPxSimplifier<BP>* _boostedSimplifier;

   SPxEquiliSC<BP> _boostedScalerUniequi;
   SPxEquiliSC<BP> _boostedScalerBiequi;
   SPxGeometSC<BP> _boostedScalerGeo1;
   SPxGeometSC<BP> _boostedScalerGeo8;
   SPxGeometSC<BP> _boostedScalerGeoequi;
   SPxLeastSqSC<BP> _boostedScalerLeastsq;

   SPxMainSM<BP> _boostedSimplifierMainSM;
   Presol<BP> _boostedSimplifierPaPILO;

   //--------------------------------------------------------------------------

   bool _isRealLPLoaded; // true indicates that the original LP is loaded in the _solver variable, hence all actions
   // are performed on the original LP.
   bool _isRealLPScaled;
   bool _applyPolishing;

   VectorBase<R> _manualLower;
   VectorBase<R> _manualUpper;
   VectorBase<R> _manualLhs;
   VectorBase<R> _manualRhs;
   VectorBase<R> _manualObj;
   SPxLPBase<R> _manualRealLP;

   ///@}


   ///@name Data for the rational LP
   ///@{

   SPxLPRational* _rationalLP;
   SLUFactorRational _rationalLUSolver;
   DataArray<int> _rationalLUSolverBind;

   LPColSetRational _slackCols;
   VectorRational _unboundedLower;
   VectorRational _unboundedUpper;
   VectorRational _unboundedLhs;
   VectorRational _unboundedRhs;
   DSVectorRational _tauColVector;
   VectorRational _feasObj;
   VectorRational _feasLhs;
   VectorRational _feasRhs;
   VectorRational _feasLower;
   VectorRational _feasUpper;
   VectorRational _modLower;
   VectorRational _modUpper;
   VectorRational _modLhs;
   VectorRational _modRhs;
   VectorRational _modObj;
   DSVectorRational _primalDualDiff;
   DataArray< typename SPxSolverBase<R>::VarStatus > _storedBasisStatusRows;
   DataArray< typename SPxSolverBase<R>::VarStatus > _storedBasisStatusCols;
   Array< UnitVectorRational* > _unitMatrixRational;
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

   ///@}


   ///@name Solution data
   ///@{

   typename SPxSolverBase<R>::Status _status;
   int _lastSolveMode;

   DataArray<typename SPxSolverBase<R>::VarStatus > _basisStatusRows;
   DataArray<typename  SPxSolverBase<R>::VarStatus > _basisStatusCols;

   // indicates wether an old basis is currently stored for warm start
   bool _hasOldBasis;
   bool _hasOldFeasBasis; // basis for testing feasibility
   bool _hasOldUnbdBasis; // basis for testing unboundedness

   // these vectors store the last basis met in precision boosting when not testing feasibility or unboundedness.
   DataArray<typename SPxSolverBase<R>::VarStatus > _oldBasisStatusRows;
   DataArray<typename  SPxSolverBase<R>::VarStatus > _oldBasisStatusCols;

   // these vectors store the last basis met when testing feasibility in precision boosting.
   DataArray<typename SPxSolverBase<R>::VarStatus > _oldFeasBasisStatusRows;
   DataArray<typename  SPxSolverBase<R>::VarStatus > _oldFeasBasisStatusCols;

   // these vectors store the last basis met when testing unboundedness in precision boosting.
   DataArray<typename SPxSolverBase<R>::VarStatus > _oldUnbdBasisStatusRows;
   DataArray<typename  SPxSolverBase<R>::VarStatus > _oldUnbdBasisStatusCols;

   // these vectors don't replace _basisStatusRows and _basisStatusCols
   // they aim to overcome the issue of having the enum VarStatus inside SPxSolverBase.
   // When calling setBasis or getBasis (from SPxSolverBase class), a specific conversion is needed.
   // Function: SPxSolverBase<BP>::setBasis(...)
   // Usage: copy _basisStatusRows(Cols) to _tmpBasisStatusRows(Cols) before calling
   // mysolver.setBasis(_tmpBasisStatusRows, _tmpBasisStatusCols)
   // Function: SPxSolverBase<BP>::getBasis(...)
   // Usage: copy _tmpBasisStatusRows(Cols) to _basisStatusRows(Cols) after calling
   // mysolver.getBasis(_tmpBasisStatusRows, _tmpBasisStatusCols, _basisStatusRows.size(), _basisStatusCols.size())
   DataArray<typename SPxSolverBase<BP>::VarStatus > _tmpBasisStatusRows;
   DataArray<typename  SPxSolverBase<BP>::VarStatus > _tmpBasisStatusCols;

   SolBase<R> _solReal;
   SolRational _solRational;
   SolRational _workSol;

   bool _hasBasis;
   bool _hasSolReal;
   bool _hasSolRational;

   ///@}

   ///@name Miscellaneous
   ///@{

   int  _optimizeCalls;
   int  _unscaleCalls;

   Rational _rationalPosone;
   Rational _rationalNegone;
   Rational _rationalZero;

   ///@}

   ///@name Constant helper methods
   ///@{

   /// extends sparse vector to hold newmax entries if and only if it holds no more free entries
   void _ensureDSVectorRationalMemory(DSVectorRational& vec, const int newmax) const;

   /// creates a permutation for removing rows/columns from an array of indices
   void _idxToPerm(int* idx, int idxSize, int* perm, int permSize) const;

   /// creates a permutation for removing rows/columns from a range of indices
   void _rangeToPerm(int start, int end, int* perm, int permSize) const;

   /// checks consistency for the boosted solver
   bool _isBoostedConsistent() const;

   /// checks consistency
   bool _isConsistent() const;

   /// should solving process be stopped?
   bool _isSolveStopped(bool& stoppedTime, bool& stoppedIter) const;

   /// determines RangeType from real bounds
   RangeType _rangeTypeReal(const R& lower, const R& upper) const;

   /// determines RangeType from rational bounds
   RangeType _rangeTypeRational(const Rational& lower, const Rational& upper) const;

   /// switches RANGETYPE_LOWER to RANGETYPE_UPPER and vice versa
   RangeType _switchRangeType(const RangeType& rangeType) const;

   /// checks whether RangeType corresponds to finite lower bound
   bool _lowerFinite(const RangeType& rangeType) const;

   /// checks whether RangeType corresponds to finite upper bound
   bool _upperFinite(const RangeType& rangeType) const;

   ///@}


   ///@name Non-constant helper methods
   ///@{

   /// adds a single row to the real LP and adjusts basis
   void _addRowReal(const LPRowBase<R>& lprow);

   /// adds a single row to the real LP and adjusts basis
   void _addRowReal(R lhs, const SVectorBase<R>& lprow, R rhs);

   /// adds multiple rows to the real LP and adjusts basis
   void _addRowsReal(const LPRowSetBase<R>& lprowset);

   /// adds a single column to the real LP and adjusts basis
   void _addColReal(const LPColReal& lpcol);

   /// adds a single column to the real LP and adjusts basis
   void _addColReal(R obj, R lower, const SVectorBase<R>& lpcol, R upper);

   /// adds multiple columns to the real LP and adjusts basis
   void _addColsReal(const LPColSetReal& lpcolset);

   /// replaces row \p i with \p lprow and adjusts basis
   void _changeRowReal(int i, const LPRowBase<R>& lprow);

   /// changes left-hand side vector for constraints to \p lhs and adjusts basis
   void _changeLhsReal(const VectorBase<R>& lhs);

   /// changes left-hand side of row \p i to \p lhs and adjusts basis
   void _changeLhsReal(int i, const R& lhs);

   /// changes right-hand side vector to \p rhs and adjusts basis
   void _changeRhsReal(const VectorBase<R>& rhs);

   /// changes right-hand side of row \p i to \p rhs and adjusts basis
   void _changeRhsReal(int i, const R& rhs);

   /// changes left- and right-hand side vectors and adjusts basis
   void _changeRangeReal(const VectorBase<R>& lhs, const VectorBase<R>& rhs);

   /// changes left- and right-hand side of row \p i and adjusts basis
   void _changeRangeReal(int i, const R& lhs, const R& rhs);

   /// replaces column \p i with \p lpcol and adjusts basis
   void _changeColReal(int i, const LPColReal& lpcol);

   /// changes vector of lower bounds to \p lower and adjusts basis
   void _changeLowerReal(const VectorBase<R>& lower);

   /// changes lower bound of column i to \p lower and adjusts basis
   void _changeLowerReal(int i, const R& lower);

   /// changes vector of upper bounds to \p upper and adjusts basis
   void _changeUpperReal(const VectorBase<R>& upper);

   /// changes \p i 'th upper bound to \p upper and adjusts basis
   void _changeUpperReal(int i, const R& upper);

   /// changes vectors of column bounds to \p lower and \p upper and adjusts basis
   void _changeBoundsReal(const VectorBase<R>& lower, const VectorBase<R>& upper);

   /// changes bounds of column \p i to \p lower and \p upper and adjusts basis
   void _changeBoundsReal(int i, const R& lower, const R& upper);

   /// changes matrix entry in row \p i and column \p j to \p val and adjusts basis
   void _changeElementReal(int i, int j, const R& val);

   /// removes row \p i and adjusts basis
   void _removeRowReal(int i);

   /// removes all rows with an index \p i such that \p perm[i] < 0; upon completion, \p perm[i] >= 0 indicates the
   /// new index where row \p i has been moved to; note that \p perm must point to an array of size at least
   /// #numRows()
   void _removeRowsReal(int perm[]);

   /// remove all rows with indices in array \p idx of size \p n; an array \p perm of size #numRows() may be passed
   /// as buffer memory
   void _removeRowsReal(int idx[], int n, int perm[]);

   /// removes rows \p start to \p end including both; an array \p perm of size #numRows() may be passed as buffer
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
   void _solveBoostedRealLPAndRecordStatistics(volatile bool* interrupt = NULL);

   /// call floating-point solver and update statistics on iterations etc.
   void _solveRealLPAndRecordStatistics(volatile bool* interrupt = NULL);

   /// reads real LP in LP or MPS format from file and returns true on success; gets row names, column names, and
   /// integer variables if desired
   bool _readFileReal(const char* filename, NameSet* rowNames = 0, NameSet* colNames = 0,
                      DIdxSet* intVars = 0);

   /// reads rational LP in LP or MPS format from file and returns true on success; gets row names, column names, and
   /// integer variables if desired
   bool _readFileRational(const char* filename, NameSet* rowNames = 0, NameSet* colNames = 0,
                          DIdxSet* intVars = 0);

   /// completes range type arrays after adding columns and/or rows
   void _completeRangeTypesRational();

   /// recomputes range types from scratch using real LP
   void _recomputeRangeTypesReal();

   /// recomputes range types from scratch using rational LP
   void _recomputeRangeTypesRational();

   /// synchronizes real LP with rational LP, i.e., copies (rounded) rational LP into real LP, without looking at the sync mode
   void _syncLPReal(bool time = true);

   /// synchronizes rational LP with real LP, i.e., copies real LP to rational LP, without looking at the sync mode
   void _syncLPRational(bool time = true);

   /// synchronizes rational solution with real solution, i.e., copies (rounded) rational solution to real solution
   void _syncRealSolution();

   /// synchronizes real solution with rational solution, i.e., copies real solution to rational solution
   void _syncRationalSolution();

   /// returns pointer to a constant unit vector available until destruction of the SoPlexBase class
   const UnitVectorRational* _unitVectorRational(const int i);

   /// parses one line in a settings file and returns true on success; note that the string is modified
   bool _parseSettingsLine(char* line, const int lineNumber);

   ///@}


   //**@name Private solving methods implemented in solverational.hpp */
   ///@{

   /// stores floating-point solution of original LP as current rational solution and ensure that solution vectors have right dimension; ensure that solution is aligned with basis
   template <typename T>
   void _storeRealSolutionAsRational(
      SolRational& sol,
      VectorBase<T>& primalReal,
      VectorBase<T>& dualReal,
      int& dualSize);

   /// computes violation of bounds during the refinement loop
   void _computeBoundsViolation(SolRational& sol, Rational& boundsViolation);

   /// computes violation of sides during the refinement loop
   void _computeSidesViolation(SolRational& sol, Rational& sideViolation);

   /// computes violation of reduced costs during the refinement loop
   void _computeReducedCostViolation(
      SolRational& sol,
      Rational& redCostViolation,
      const bool& maximizing);

   /// computes dual violation during the refinement loop
   void _computeDualViolation(
      SolRational& sol,
      Rational& dualViolation,
      const bool& maximizing);

   /// checks termination criteria for refinement loop
   bool _isRefinementOver(
      bool& primalFeasible,
      bool& dualFeasible,
      Rational& boundsViolation,
      Rational& sideViolation,
      Rational& redCostViolation,
      Rational& dualViolation,
      int minIRRoundsRemaining,
      bool& stoppedTime,
      bool& stoppedIter,
      int numFailedRefinements);

   /// checks refinement loop progress
   void _checkRefinementProgress(
      Rational& boundsViolation,
      Rational& sideViolation,
      Rational& redCostViolation,
      Rational& dualViolation,
      Rational& maxViolation,
      Rational& bestViolation,
      const Rational& violationImprovementFactor,
      int& numFailedRefinements);

   /// performs rational reconstruction and/or factorizationd
   void _ratrecAndOrRatfac(
      int& minIRRoundsRemaining,
      int& lastStallIterations,
      int& numberOfIterations,
      bool& factorSolNewBasis,
      int& nextRatrec,
      const Rational& errorCorrectionFactor,
      Rational& errorCorrection,
      Rational& maxViolation,
      SolRational& sol,
      bool& primalFeasible,
      bool& dualFeasible,
      bool& stoppedTime,
      bool& stoppedIter,
      bool& error,
      bool& breakAfter,
      bool& continueAfter);

   /// forces value of given nonbasic variable to bound
   void _forceNonbasicToBound(
      SolRational& sol,
      int& c,
      const int& maxDimRational,
      bool toLower);

   /// computes primal scaling factor; limit increase in scaling by tolerance used in floating point solve
   void _computePrimalScalingFactor(
      Rational& maxScale,
      Rational& primalScale,
      Rational& boundsViolation,
      Rational& sideViolation,
      Rational& redCostViolation);

   /// computes dual scaling factor; limit increase in scaling by tolerance used in floating point solve
   void _computeDualScalingFactor(
      Rational& maxScale,
      Rational& primalScale,
      Rational& dualScale,
      Rational& redCostViolation,
      Rational& dualViolation);

   /// applies scaled bounds
   template <typename T>
   void _applyScaledBounds(SPxSolverBase<T>& solver, Rational& primalScale);

   /// applies scaled sides
   template <typename T>
   void _applyScaledSides(SPxSolverBase<T>& solver, Rational& primalScale);

   /// applies scaled objective function
   template <typename T>
   void _applyScaledObj(SPxSolverBase<T>& solver, Rational& dualScale, SolRational& sol);

   /// evaluates result of solve. Return true if the algorithm must to stopped, false otherwise.
   template <typename T>
   bool _evaluateResult(
      SPxSolverBase<T>& solver,
      typename SPxSolverBase<T>::Status result,
      bool usingRefinedLP,
      SolRational& sol,
      VectorBase<T>& dualReal,
      bool& infeasible,
      bool& unbounded,
      bool& stoppedTime,
      bool& stoppedIter,
      bool& error);

   /// corrects primal solution and aligns with basis
   template <typename T>
   void _correctPrimalSolution(
      SolRational& sol,
      Rational& primalScale,
      int& primalSize,
      const int& maxDimRational,
      VectorBase<T>& primalReal);

   /// updates or recomputes slacks depending on which looks faster
   void _updateSlacks(SolRational& sol, int& primalSize);

   /// corrects dual solution and aligns with basis
   template <typename T>
   void _correctDualSolution(
      SPxSolverBase<T>& solver,
      SolRational& sol,
      const bool& maximizing,
      VectorBase<T>& dualReal,
      Rational& dualScale,
      int& dualSize,
      const int& maxDimRational);

   /// updates or recomputes reduced cost values depending on which looks faster; adding one to the length of the
   /// dual vector accounts for the objective function vector
   void _updateReducedCosts(SolRational& sol, int& dualSize, const int& numCorrectedPrimals);

   ///@todo precision-boosting move some place else
   /// converts the given DataArray of VarStatus to boostedPrecision
   void _convertDataArrayVarStatusToBoosted(
      DataArray< typename SPxSolverBase<R>::VarStatus >& base,
      DataArray< typename SPxSolverBase<BP>::VarStatus >& copy);

   ///@todo precision-boosting move some place else
   /// converts the given DataArray of VarStatus to R precision
   void _convertDataArrayVarStatusToRPrecision(
      DataArray< typename SPxSolverBase<BP>::VarStatus >& base,
      DataArray< typename SPxSolverBase<R>::VarStatus >& copy);

   /// disable initial precision solver and switch to boosted solver
   void _switchToBoosted();

   /// setup boosted solver before launching iteration
   void _setupBoostedSolver();

   /// increase the multiprecision, return false if maximum precision is reached, true otherwise
   bool _boostPrecision();

   /// reset the boosted precision to the default value
   void _resetBoostedPrecision();

   /// setup recovery mecanism using multiprecision, return false if maximum precision reached, true otherwise
   bool _setupBoostedSolverAfterRecovery();

   /// return true if slack basis has to be loaded for boosted solver
   bool _isBoostedStartingFromSlack(bool initialSolve = true);

   /// indicate if we are testing feasibility, unboundedness or neither
   void _switchToStandardMode();
   void _switchToFeasMode();
   void _switchToUnbdMode();

   /// check if we are testing feasibility, unboundedness or neither
   bool _inStandardMode();
   bool _inFeasMode();
   bool _inUnbdMode();

   // stores given basis in old basis attributes: _oldBasisStatusRows, _oldFeasBasisStatusRows, _oldUnbdBasisStatusRows (and ...Cols)
   void _storeBasisAsOldBasis(DataArray< typename SPxSolverBase<R>::VarStatus >& rows,
                              DataArray< typename SPxSolverBase<R>::VarStatus >& cols);

   // stores given basis in old basis attributes: _oldBasisStatusRows, _oldFeasBasisStatusRows, _oldUnbdBasisStatusRows (and ...Cols)
   void _storeBasisAsOldBasisBoosted(DataArray< typename SPxSolverBase<BP>::VarStatus >& rows,
                                     DataArray< typename SPxSolverBase<BP>::VarStatus >& cols);

   // get the last advanced and stable basis stored by the initial solver and store it as old basis, unsimplify basis if simplifier activated
   void _storeLastStableBasis(bool vanished);

   // get the last advanced and stable basis stored by the boosted solver and store it as old basis, unsimplify basis if simplifier activated
   void _storeLastStableBasisBoosted(bool vanished);

   // load old basis in solver. The old basis loaded depends on the certificate mode (feasibility, unboundedness, or neither)
   bool _loadBasisFromOldBasis(bool boosted);

   // update statistics for precision boosting
   void _updateBoostingStatistics();

   /// solves current problem using multiprecision floating-point solver
   /// return false if a new boosted iteration is necessary, true otherwise
   void _solveRealForRationalBoostedStable(
      SolRational& sol,
      bool& primalFeasible,
      bool& dualFeasible,
      bool& infeasible,
      bool& unbounded,
      bool& stoppedTime,
      bool& stoppedIter,
      bool& error,
      bool& needNewBoostedIt);

   /// solves current problem with iterative refinement and recovery mechanism using boosted solver
   void _performOptIRStableBoosted(
      SolRational& sol,
      bool acceptUnbounded,
      bool acceptInfeasible,
      int minIRRoundsRemaining,
      bool& primalFeasible,
      bool& dualFeasible,
      bool& infeasible,
      bool& unbounded,
      bool& stoppedTime,
      bool& stoppedIter,
      bool& error,
      bool& needNewBoostedIt);

   /// perform iterative refinement using the right precision
   void _performOptIRWrapper(
      SolRational& sol,
      bool acceptUnbounded,
      bool acceptInfeasible,
      int minIRRoundsRemaining,
      bool& primalFeasible,
      bool& dualFeasible,
      bool& infeasible,
      bool& unbounded,
      bool& stoppedTime,
      bool& stoppedIter,
      bool& error
   );

   /// solves current problem using double floating-point solver
   void _solveRealForRationalStable(
      SolRational& sol,
      bool& primalFeasible,
      bool& dualFeasible,
      bool& infeasible,
      bool& unbounded,
      bool& stoppedTime,
      bool& stoppedIter,
      bool& error);

   /// solves current problem with iterative refinement and recovery mechanism
   void _performOptIRStable(SolRational& sol,
                            bool acceptUnbounded,
                            bool acceptInfeasible,
                            int minIRRoundsRemaining,
                            bool& primalFeasible,
                            bool& dualFeasible,
                            bool& infeasible,
                            bool& unbounded,
                            bool& stoppedTime,
                            bool& stoppedIter,
                            bool& error);

   /// performs iterative refinement on the auxiliary problem for testing unboundedness
   void _performUnboundedIRStable(SolRational& sol, bool& hasUnboundedRay, bool& stoppedTime,
                                  bool& stoppedIter, bool& error);

   /// performs iterative refinement on the auxiliary problem for testing feasibility
   void _performFeasIRStable(SolRational& sol, bool& withDualFarkas, bool& stoppedTime,
                             bool& stoppedIter, bool& error);

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
   typename SPxSolverBase<R>::Status _solveRealForRational(bool fromscratch, VectorBase<R>& primal,
         VectorBase<R>& dual,
         DataArray< typename SPxSolverBase<R>::VarStatus >& basisStatusRows,
         DataArray< typename SPxSolverBase<R>::VarStatus >& basisStatusCols);

   /// solves real LP with recovery mechanism
   typename SPxSolverBase<R>::Status _solveRealStable(bool acceptUnbounded, bool acceptInfeasible,
         VectorBase<R>& primal, VectorBase<R>& dual,
         DataArray< typename SPxSolverBase<R>::VarStatus >& basisStatusRows,
         DataArray< typename SPxSolverBase<R>::VarStatus >& basisStatusCols,
         const bool forceNoSimplifier = false);

   /// solves real LP during iterative refinement
   void _solveRealForRationalBoosted(
      VectorBase<BP>& primal, VectorBase<BP>& dual,
      DataArray< typename SPxSolverBase<R>::VarStatus >& basisStatusRows,
      DataArray< typename SPxSolverBase<R>::VarStatus >& basisStatusCols,
      typename SPxSolverBase<BP>::Status& boostedResult, bool initialSolve);

   /// computes rational inverse of basis matrix as defined by _rationalLUSolverBind
   void _computeBasisInverseRational();

   /// factorizes rational basis matrix in column representation
   void _factorizeColumnRational(SolRational& sol,
                                 DataArray< typename SPxSolverBase<R>::VarStatus >& basisStatusRows,
                                 DataArray< typename SPxSolverBase<R>::VarStatus >& basisStatusCols, bool& stoppedTime,
                                 bool& stoppedIter, bool& error, bool& optimal);

   /// attempts rational reconstruction of primal-dual solution
   bool _reconstructSolutionRational(SolRational& sol,
                                     DataArray< typename SPxSolverBase<R>::VarStatus >& basisStatusRows,
                                     DataArray< typename SPxSolverBase<R>::VarStatus >& basisStatusCols,
                                     const Rational& denomBoundSquared);
   ///@}

   ///@name Private solving methods implemented in solvereal.cpp
   ///@{

   /// solves the templated LP
   void _optimize(volatile bool* interrupt = NULL);

   /// temporary fix for Rational
   void _optimizeRational(volatile bool* interrupt = NULL);

   /// checks result of the solving process and solves again without preprocessing if necessary
   void _evaluateSolutionReal(typename SPxSimplifier<R>::Result simplificationStatus);

   /// solves real LP with/without preprocessing
   void _preprocessAndSolveReal(bool applyPreprocessing, volatile bool* interrupt = NULL);

   /// loads original problem into solver and solves again after it has been solved to optimality with preprocessing
   void _resolveWithoutPreprocessing(typename SPxSimplifier<R>::Result simplificationStatus);

   /// verify computed solution and resolve if necessary
   void _verifySolutionReal();

   /// verify computed obj stop and resolve if necessary
   void _verifyObjLimitReal();

   /// stores solution of the real LP; before calling this, the real LP must be loaded in the solver and solved (again)
   void _storeSolutionReal(bool verify = true);

   /// stores solution from the simplifier because problem vanished in presolving step
   void _storeSolutionRealFromPresol();

   /// unscales stored solution to remove internal or external scaling of LP
   void _unscaleSolutionReal(SPxLPBase<R>& LP, bool persistent = true);

   /// load original LP and possibly setup a slack basis
   void _loadRealLP(bool initBasis);

   /// check scaling of LP
   void _checkScaling(SPxLPBase<R>* origLP) const;

   /// check correctness of (un)scaled basis matrix operations
   void _checkBasisScaling();

   /// check whether persistent scaling is supposed to be reapplied again after unscaling
   bool _reapplyPersistentScaling() const;

   /// checks the dual feasibility of the current basis
   bool checkBasisDualFeasibility(VectorBase<R> feasVec);
   ///@}
};

/* Backwards compatibility */
typedef SoPlexBase<Real> SoPlex;
// A header file containing all the general templated functions

} // namespace soplex

// General templated function
#include "soplex.hpp"
#include "soplex/solverational.hpp"
#include "soplex/testsoplex.hpp"
#include "soplex/solvereal.hpp"

#endif // _SOPLEX_H_
