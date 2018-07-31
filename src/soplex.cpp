/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the class library                   */
/*       SoPlex --- the Sequential object-oriented simPlex.                  */
/*                                                                           */
/*    Copyright (C) 1996-2017 Konrad-Zuse-Zentrum                            */
/*                            fuer Informationstechnik Berlin                */
/*                                                                           */
/*  SoPlex is distributed under the terms of the ZIB Academic Licence.       */
/*                                                                           */
/*  You should have received a copy of the ZIB Academic License              */
/*  along with SoPlex; see the file COPYING. If not email to soplex@zib.de.  */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/**@file  soplex.cpp
 * @brief Preconfigured SoPlex LP solver
 */

#ifndef SOPLEX_LEGACY
#include <assert.h>
#include "limits.h"
#include <iostream>

#ifndef _MSC_VER
#include <strings.h>
#endif

#include "soplex.h"
#include "soplex/spxfileio.h"
#include "soplex/statistics.h"
#include "soplex/mpsinput.h"

/// maximum length of lines in settings file
#define SET_MAX_LINE_LEN 500
/// default setting for LU refactorization interval
#define DEFAULT_REFACTOR_INTERVAL 200

#ifdef _MSC_VER
#define strncasecmp _strnicmp
#endif

namespace soplex
{
  template <>
  int SoPlex<Real>::intParam(const IntParam param) const;

  template <>
	bool SoPlex<Real>::areLPsInSync(const bool checkVecVals, const bool checkMatVals, const bool quiet) const;

  template <>
  typename SoPlex<Real>::Settings& SoPlex<Real>::Settings::operator=(const Settings& settings);

  template <>
  const SVectorReal& SoPlex<Real>::colVectorRealInternal(int i) const;

  template <>
  const Rational& SoPlex<Real>::maxObjRational(int i) const;

  template <>
  typename SPxSolver<Real>::VarStatus SoPlex<Real>::basisRowStatus(int row) const;

  template <>
  typename SPxSolver<Real>::VarStatus SoPlex<Real>::basisColStatus(int col) const;

  template <>
	bool SoPlex<Real>::setSettings(const Settings& newSettings, const bool init);

  template <>
	bool SoPlex<Real>::saveSettingsFile(const char* filename, const bool onlyChanged) const;

  template <>
  void SoPlex<Real>::printShortStatistics(std::ostream& os);

  template <>
  void SoPlex<Real>::printStatus(std::ostream& os, typename SPxSolver<Real>::Status stat);

  template <>
  void SoPlex<Real>::setRandomSeed(unsigned int seed);

  template <>
  void SoPlex<Real>::_idxToPerm(int* idx, int idxSize, int* perm, int permSize) const;

  template <>
  void SoPlex<Real>::_rangeToPerm(int start, int end, int* perm, int permSize) const;

  template <>
	bool SoPlex<Real>::_isConsistent() const;

  template <>
  typename SoPlex<Real>::RangeType SoPlex<Real>::_rangeTypeReal(const Real& lower, const Real& upper) const;

  template <>
  void SoPlex<Real>::_addRowReal(const LPRowReal& lprow);

  template <>
	bool SoPlex<Real>::setIntParam(const IntParam param, const int value, const bool init);

  template <>
  typename SoPlex<Real>::RangeType SoPlex<Real>::_rangeTypeRational(const Rational& lower, const Rational& upper) const;

  template <>
  void SoPlex<Real>::_addRowReal(Real lhs, const SVectorReal& lprow, Real rhs);

  template <>
  void SoPlex<Real>::_addRowsReal(const LPRowSetReal& lprowset);

  template <>
  void SoPlex<Real>::_addColReal(const LPColReal& lpcol);

  template <>
  void SoPlex<Real>::_addColReal(Real obj, Real lower, const SVectorReal& lpcol, Real upper);

  template <>
  void SoPlex<Real>::_addColsReal(const LPColSetReal& lpcolset);

  template <>
  void SoPlex<Real>::_changeRowReal(int i, const LPRowReal& lprow);

  template <>
  void SoPlex<Real>::_changeLhsReal(const VectorReal& lhs);

  template <>
  void SoPlex<Real>::_changeLhsReal(int i, const Real& lhs);

  template <>
  void SoPlex<Real>::_changeRhsReal(const VectorReal& rhs);

  template <>
  void SoPlex<Real>::_changeRhsReal(int i, const Real& rhs);

  template <>
  void SoPlex<Real>::_changeRangeReal(const VectorReal& lhs, const VectorReal& rhs);

  template <>
  void SoPlex<Real>::_changeRangeReal(int i, const Real& lhs, const Real& rhs);

  template <>
  void SoPlex<Real>::_changeColReal(int i, const LPColReal& lpcol);

  template <>
  void SoPlex<Real>::_changeLowerReal(const VectorReal& lower);

  template <>
  void SoPlex<Real>::_changeLowerReal(int i, const Real& lower);

  template <>
  void SoPlex<Real>::_changeUpperReal(const VectorReal& upper);

  template <>
  void SoPlex<Real>::_changeUpperReal(int i, const Real& upper);

  template <>
  void SoPlex<Real>::_changeBoundsReal(const VectorReal& lower, const VectorReal& upper);

  template <>
  void SoPlex<Real>::_changeBoundsReal(int i, const Real& lower, const Real& upper);

  template <>
  void SoPlex<Real>::_changeElementReal(int i, int j, const Real& val);

  template <>
  void SoPlex<Real>::_removeRowReal(int i);

  template <>
  void SoPlex<Real>::_removeRowsReal(int perm[]);

  template <>
  void SoPlex<Real>::_removeColReal(int i);

  template <>
  void SoPlex<Real>::_removeColsReal(int perm[]);

  template <>
  void SoPlex<Real>::_invalidateSolution();

  template <>
  void SoPlex<Real>::_ensureRationalLP();

  template <>
  void SoPlex<Real>::_ensureRealLPLoaded();

  template <>
	bool SoPlex<Real>::_readFileReal(const char* filename, NameSet* rowNames, NameSet* colNames, DIdxSet* intVars);

  template <>
  void SoPlex<Real>::_completeRangeTypesRational();

  template <>
  void SoPlex<Real>::_recomputeRangeTypesRational();

  template <>
  void SoPlex<Real>::_syncLPReal(bool time);

  template <>
  void SoPlex<Real>::_syncLPRational(bool time);

  template <>
  void SoPlex<Real>::_syncRealSolution();

  template <>
  void SoPlex<Real>::_syncRationalSolution();

  template <>
  const UnitVectorRational* SoPlex<Real>::_unitVectorRational(const int i);

  template <>
	bool SoPlex<Real>::_parseSettingsLine(char* line, const int lineNumber);

  template <>
	bool SoPlex<Real>::_readFileRational(const char* filename, NameSet* rowNames, NameSet* colNames, DIdxSet* intVars);

  template <>
  SoPlex<Real>::Settings::BoolParam::BoolParam() {
    // should lifting be used to reduce range of nonzero matrix coefficients?
    name[SoPlex<Real>::LIFTING] = "lifting";
    description[SoPlex<Real>::LIFTING] = "should lifting be used to reduce range of nonzero matrix coefficients?";
    defaultValue[SoPlex<Real>::LIFTING] = false;

    // should LP be transformed to equality form before a rational solve?
    name[SoPlex<Real>::EQTRANS] = "eqtrans";
    description[SoPlex<Real>::EQTRANS] = "should LP be transformed to equality form before a rational solve?";
    defaultValue[SoPlex<Real>::EQTRANS] = false;

    // should dual infeasibility be tested in order to try to return a dual solution even if primal infeasible?
    name[SoPlex<Real>::TESTDUALINF] = "testdualinf";
    description[SoPlex<Real>::TESTDUALINF] = "should dual infeasibility be tested in order to try to return a dual solution even if primal infeasible?";
    defaultValue[SoPlex<Real>::TESTDUALINF] = false;

    // should a rational factorization be performed after iterative refinement?
    name[SoPlex<Real>::RATFAC] = "ratfac";
    description[SoPlex<Real>::RATFAC] = "should a rational factorization be performed after iterative refinement?";
    defaultValue[SoPlex<Real>::RATFAC] = true;

    // should the decomposition based dual simplex be used to solve the LP? Setting this to true forces the solve mode to
    // SOLVEMODE_REAL and the basis representation to REPRESENTATION_ROW
    name[SoPlex<Real>::USEDECOMPDUALSIMPLEX] = "decompositiondualsimplex";
    description[SoPlex<Real>::USEDECOMPDUALSIMPLEX] = "should the decomposition based dual simplex be used to solve the LP?";
    defaultValue[SoPlex<Real>::USEDECOMPDUALSIMPLEX] = false;

    // should the degeneracy be computed for each basis?
    name[SoPlex<Real>::COMPUTEDEGEN] = "computedegen";
    description[SoPlex<Real>::COMPUTEDEGEN] = "should the degeneracy be computed for each basis?";
    defaultValue[SoPlex<Real>::COMPUTEDEGEN] = false;

    // should the dual of the complementary problem be used in the decomposition simplex?
    name[SoPlex<Real>::USECOMPDUAL] = "usecompdual";
    description[SoPlex<Real>::USECOMPDUAL] = "should the dual of the complementary problem be used in the decomposition simplex?";
    defaultValue[SoPlex<Real>::USECOMPDUAL] = false;

    /// should row and bound violations be computed explicitly in the update of reduced problem in the decomposition
    // simplex
    name[SoPlex<Real>::EXPLICITVIOL] = "explicitviol";
    description[SoPlex<Real>::EXPLICITVIOL] = "Should violations of the original problem be explicitly computed in the decomposition simplex?";
    defaultValue[SoPlex<Real>::EXPLICITVIOL] = false;

    // should cycling solutions be accepted during iterative refinement?
    name[SoPlex<Real>::ACCEPTCYCLING] = "acceptcycling";
    description[SoPlex<Real>::ACCEPTCYCLING] = "should cycling solutions be accepted during iterative refinement?";
    defaultValue[SoPlex<Real>::ACCEPTCYCLING] = false;

    // apply rational reconstruction after each iterative refinement?
    name[SoPlex<Real>::RATREC] = "ratrec";
    description[SoPlex<Real>::RATREC] = "apply rational reconstruction after each iterative refinement?";
    defaultValue[SoPlex<Real>::RATREC] = true;

    // round scaling factors for iterative refinement to powers of two?
    name[SoPlex<Real>::POWERSCALING] = "powerscaling";
    description[SoPlex<Real>::POWERSCALING] = "round scaling factors for iterative refinement to powers of two?";
    defaultValue[SoPlex<Real>::POWERSCALING] = true;

    // continue iterative refinement with exact basic solution if not optimal?
    name[SoPlex<Real>::RATFACJUMP] = "ratfacjump";
    description[SoPlex<Real>::RATFACJUMP] = "continue iterative refinement with exact basic solution if not optimal?";
    defaultValue[SoPlex<Real>::RATFACJUMP] = false;

    // use bound flipping also for row representation?
    name[SoPlex<Real>::ROWBOUNDFLIPS] = "rowboundflips";
    description[SoPlex<Real>::ROWBOUNDFLIPS] = "use bound flipping also for row representation?";
    defaultValue[SoPlex<Real>::ROWBOUNDFLIPS] = false;

    // use persistent scaling?
    name[SoPlex<Real>::PERSISTENTSCALING] = "persistentscaling";
    description[SoPlex<Real>::PERSISTENTSCALING] = "should persistent scaling be used?";
    defaultValue[SoPlex<Real>::PERSISTENTSCALING] = true;

    // perturb the entire problem or only the relevant bounds of s single pivot?
    name[SoPlex<Real>::FULLPERTURBATION] = "fullperturbation";
    description[SoPlex<Real>::FULLPERTURBATION] = "should perturbation be applied to the entire problem?";
    defaultValue[SoPlex<Real>::FULLPERTURBATION] = false;
  }

  template <class R>
  SoPlex<R>::Settings::IntParam::IntParam() {
    // objective sense
    name[SoPlex<R>::OBJSENSE] = "objsense";
    description[SoPlex<R>::OBJSENSE] = "objective sense (-1 - minimize, +1 - maximize)";
    lower[SoPlex<R>::OBJSENSE] = -1;
    upper[SoPlex<R>::OBJSENSE] = 1;
    defaultValue[SoPlex<R>::OBJSENSE] = SoPlex<R>::OBJSENSE_MAXIMIZE;

    // type of computational form, i.e., column or row representation
    name[SoPlex<R>::REPRESENTATION] = "representation";
    description[SoPlex<R>::REPRESENTATION] = "type of computational form (0 - auto, 1 - column representation, 2 - row representation)";
    lower[SoPlex<R>::REPRESENTATION] = 0;
    upper[SoPlex<R>::REPRESENTATION] = 2;
    defaultValue[SoPlex<R>::REPRESENTATION] = SoPlex<R>::REPRESENTATION_AUTO;

    // type of algorithm, i.e., primal or dual
    name[SoPlex<R>::ALGORITHM] = "algorithm";
    description[SoPlex<R>::ALGORITHM] = "type of algorithm (0 - primal, 1 - dual)";
    lower[SoPlex<R>::ALGORITHM] = 0;
    upper[SoPlex<R>::ALGORITHM] = 1;
    defaultValue[SoPlex<R>::ALGORITHM] = SoPlex<R>::ALGORITHM_DUAL;

    // type of LU update
    name[SoPlex<R>::FACTOR_UPDATE_TYPE] = "factor_update_type";
    description[SoPlex<R>::FACTOR_UPDATE_TYPE] = "type of LU update (0 - eta update, 1 - Forrest-Tomlin update)";
    lower[SoPlex<R>::FACTOR_UPDATE_TYPE] = 0;
    upper[SoPlex<R>::FACTOR_UPDATE_TYPE] = 1;
    defaultValue[SoPlex<R>::FACTOR_UPDATE_TYPE] = SoPlex<R>::FACTOR_UPDATE_TYPE_FT;

    // maximum number of updates without fresh factorization
    name[SoPlex<R>::FACTOR_UPDATE_MAX] = "factor_update_max";
    description[SoPlex<R>::FACTOR_UPDATE_MAX] = "maximum number of LU updates without fresh factorization (0 - auto)";
    lower[SoPlex<R>::FACTOR_UPDATE_MAX] = 0;
    upper[SoPlex<R>::FACTOR_UPDATE_MAX] = INT_MAX;
    defaultValue[SoPlex<R>::FACTOR_UPDATE_MAX] = 0;

    // iteration limit (-1 if unlimited)
    name[SoPlex<R>::ITERLIMIT] = "iterlimit";
    description[SoPlex<R>::ITERLIMIT] = "iteration limit (-1 - no limit)";
    lower[SoPlex<R>::ITERLIMIT] = -1;
    upper[SoPlex<R>::ITERLIMIT] = INT_MAX;
    defaultValue[SoPlex<R>::ITERLIMIT] = -1;

    // refinement limit (-1 if unlimited)
    name[SoPlex<R>::REFLIMIT] = "reflimit";
    description[SoPlex<R>::REFLIMIT] = "refinement limit (-1 - no limit)";
    lower[SoPlex<R>::REFLIMIT] = -1;
    upper[SoPlex<R>::REFLIMIT] = INT_MAX;
    defaultValue[SoPlex<R>::REFLIMIT] = -1;

    // stalling refinement limit (-1 if unlimited)
    name[SoPlex<R>::STALLREFLIMIT] = "stallreflimit";
    description[SoPlex<R>::STALLREFLIMIT] = "stalling refinement limit (-1 - no limit)";
    lower[SoPlex<R>::STALLREFLIMIT] = -1;
    upper[SoPlex<R>::STALLREFLIMIT] = INT_MAX;
    defaultValue[SoPlex<R>::STALLREFLIMIT] = -1;

    // display frequency
    name[SoPlex<R>::DISPLAYFREQ] = "displayfreq";
    description[SoPlex<R>::DISPLAYFREQ] = "display frequency";
    lower[SoPlex<R>::DISPLAYFREQ] = 1;
    upper[SoPlex<R>::DISPLAYFREQ] = INT_MAX;
    defaultValue[SoPlex<R>::DISPLAYFREQ] = 200;

    // verbosity level
    name[SoPlex<R>::VERBOSITY] = "verbosity";
    description[SoPlex<R>::VERBOSITY] = "verbosity level (0 - error, 1 - warning, 2 - debug, 3 - normal, 4 - high, 5 - full)";
    lower[SoPlex<R>::VERBOSITY] = 0;
    upper[SoPlex<R>::VERBOSITY] = 5;
    defaultValue[SoPlex<R>::VERBOSITY] = SoPlex<R>::VERBOSITY_NORMAL;
    //_intParamDefault[SoPlex<R>::VERBOSITY] = SoPlex<R>::VERBOSITY_FULL;

    // type of simplifier
    name[SoPlex<R>::SIMPLIFIER] = "simplifier";
    description[SoPlex<R>::SIMPLIFIER] = "simplifier (0 - off, 1 - auto)";
    lower[SoPlex<R>::SIMPLIFIER] = 0;
    upper[SoPlex<R>::SIMPLIFIER] = 1;
    defaultValue[SoPlex<R>::SIMPLIFIER] = SoPlex<R>::SIMPLIFIER_AUTO;

    // type of scaler
    name[SoPlex<R>::SCALER] = "scaler";
    description[SoPlex<R>::SCALER] = "scaling (0 - off, 1 - uni-equilibrium, 2 - bi-equilibrium, 3 - geometric, 4 - iterated geometric, 5 - least squares, 6 - geometric-equilibrium)";
    lower[SoPlex<R>::SCALER] = 0;
    upper[SoPlex<R>::SCALER] = 6;
    defaultValue[SoPlex<R>::SCALER] = SoPlex<R>::SCALER_BIEQUI;

    // type of starter used to create crash basis
    name[SoPlex<R>::STARTER] = "starter";
    description[SoPlex<R>::STARTER] = "crash basis generated when starting from scratch (0 - none, 1 - weight, 2 - sum, 3 - vector)";
    lower[SoPlex<R>::STARTER] = 0;
    upper[SoPlex<R>::STARTER] = 3;
    defaultValue[SoPlex<R>::STARTER] = SoPlex<R>::STARTER_OFF;

    // type of pricer
    name[SoPlex<R>::PRICER] = "pricer";
    description[SoPlex<R>::PRICER] = "pricing method (0 - auto, 1 - dantzig, 2 - parmult, 3 - devex, 4 - quicksteep, 5 - steep)";
    lower[SoPlex<R>::PRICER] = 0;
    upper[SoPlex<R>::PRICER] = 5;
    defaultValue[SoPlex<R>::PRICER] = SoPlex<R>::PRICER_AUTO;

    // type of ratio test
    name[SoPlex<R>::RATIOTESTER] = "ratiotester";
    description[SoPlex<R>::RATIOTESTER] = "method for ratio test (0 - textbook, 1 - harris, 2 - fast, 3 - boundflipping)";
    lower[SoPlex<R>::RATIOTESTER] = 0;
    upper[SoPlex<R>::RATIOTESTER] = 3;
    defaultValue[SoPlex<R>::RATIOTESTER] = SoPlex<R>::RATIOTESTER_BOUNDFLIPPING;

    // mode for synchronizing real and rational LP
    name[SoPlex<R>::SYNCMODE] = "syncmode";
    description[SoPlex<R>::SYNCMODE] = "mode for synchronizing real and rational LP (0 - store only real LP, 1 - auto, 2 - manual)";
    lower[SoPlex<R>::SYNCMODE] = 0;
    upper[SoPlex<R>::SYNCMODE] = 2;
    defaultValue[SoPlex<R>::SYNCMODE] = SoPlex<R>::SYNCMODE_ONLYREAL;

    // mode for reading LP files
    name[SoPlex<R>::READMODE] = "readmode";
    description[SoPlex<R>::READMODE] = "mode for reading LP files (0 - floating-point, 1 - rational)";
    lower[SoPlex<R>::READMODE] = 0;
    upper[SoPlex<R>::READMODE] = 1;
    defaultValue[SoPlex<R>::READMODE] = SoPlex<R>::READMODE_REAL;

    // mode for iterative refinement strategy
    name[SoPlex<R>::SOLVEMODE] = "solvemode";
    description[SoPlex<R>::SOLVEMODE] = "mode for iterative refinement strategy (0 - floating-point solve, 1 - auto, 2 - exact rational solve)";
    lower[SoPlex<R>::SOLVEMODE] = 0;
    upper[SoPlex<R>::SOLVEMODE] = 2;
    defaultValue[SoPlex<R>::SOLVEMODE] = SoPlex<R>::SOLVEMODE_AUTO;

    // mode for iterative refinement strategy
    name[SoPlex<R>::CHECKMODE] = "checkmode";
    description[SoPlex<R>::CHECKMODE] = "mode for a posteriori feasibility checks (0 - floating-point check, 1 - auto, 2 - exact rational check)";
    lower[SoPlex<R>::CHECKMODE] = 0;
    upper[SoPlex<R>::CHECKMODE] = 2;
    defaultValue[SoPlex<R>::CHECKMODE] = SoPlex<R>::CHECKMODE_AUTO;

    // type of timing
    name[SoPlex<R>::TIMER] = "timer";
    description[SoPlex<R>::TIMER] = "type of timer (1 - cputime, aka. usertime, 2 - wallclock time, 0 - no timing)";
    lower[SoPlex<R>::TIMER] = 0;
    upper[SoPlex<R>::TIMER] = 2;
    defaultValue[SoPlex<R>::TIMER] = SoPlex<R>::TIMER_CPU;

    // mode for hyper sparse pricing
    name[SoPlex<R>::HYPER_PRICING] = "hyperpricing";
    description[SoPlex<R>::HYPER_PRICING] = "mode for hyper sparse pricing (0 - off, 1 - auto, 2 - always)";
    lower[SoPlex<R>::HYPER_PRICING] = 0;
    upper[SoPlex<R>::HYPER_PRICING] = 2;
    defaultValue[SoPlex<R>::HYPER_PRICING] = SoPlex<R>::HYPER_PRICING_AUTO;

    // minimum number of stalling refinements since last pivot to trigger rational factorization
    name[SoPlex<R>::RATFAC_MINSTALLS] = "ratfac_minstalls";
    description[SoPlex<R>::RATFAC_MINSTALLS] = "minimum number of stalling refinements since last pivot to trigger rational factorization";
    lower[SoPlex<R>::RATFAC_MINSTALLS] = 0;
    upper[SoPlex<R>::RATFAC_MINSTALLS] = INT_MAX;
    defaultValue[SoPlex<R>::RATFAC_MINSTALLS] = 2;

    // maximum number of conjugate gradient iterations in least square scaling
    name[SoPlex<R>::LEASTSQ_MAXROUNDS] = "leastsq_maxrounds";
    description[SoPlex<R>::LEASTSQ_MAXROUNDS] = "maximum number of conjugate gradient iterations in least square scaling";
    lower[SoPlex<R>::LEASTSQ_MAXROUNDS] = 0;
    upper[SoPlex<R>::LEASTSQ_MAXROUNDS] = INT_MAX;
    defaultValue[SoPlex<R>::LEASTSQ_MAXROUNDS] = 50;

    // mode for solution polishing
    name[SoPlex<R>::SOLUTION_POLISHING] = "solution_polishing";
    description[SoPlex<R>::SOLUTION_POLISHING] = "mode for solution polishing (0 - off, 1 - max basic slack, 2 - min basic slack)";
    lower[SoPlex<R>::SOLUTION_POLISHING] = 0;
    upper[SoPlex<R>::SOLUTION_POLISHING] = 2;
    defaultValue[SoPlex<R>::SOLUTION_POLISHING] = SoPlex<R>::POLISHING_OFF;

    // the number of iterations before the decomposition simplex initialisation is terminated.
    name[SoPlex<R>::DECOMP_ITERLIMIT] = "decomp_iterlimit";
    description[SoPlex<R>::DECOMP_ITERLIMIT] = "the number of iterations before the decomposition simplex initialisation solve is terminated";
    lower[SoPlex<R>::DECOMP_ITERLIMIT] = 1;
    upper[SoPlex<R>::DECOMP_ITERLIMIT] = INT_MAX;
    defaultValue[SoPlex<R>::DECOMP_ITERLIMIT] = 100;

    // maximum number of violated rows added in each iteration of the decomposition simplex
    name[SoPlex<R>::DECOMP_MAXADDEDROWS] = "decomp_maxaddedrows";
    description[SoPlex<R>::DECOMP_MAXADDEDROWS] = "maximum number of rows that are added to the reduced problem when using the decomposition based simplex";
    lower[SoPlex<R>::DECOMP_MAXADDEDROWS] = 1;
    upper[SoPlex<R>::DECOMP_MAXADDEDROWS] = INT_MAX;
    defaultValue[SoPlex<R>::DECOMP_MAXADDEDROWS] = 500;

    // maximum number of violated rows added in each iteration of the decomposition simplex
    name[SoPlex<R>::DECOMP_DISPLAYFREQ] = "decomp_displayfreq";
    description[SoPlex<R>::DECOMP_DISPLAYFREQ] = "the frequency that the decomposition based simplex status output is displayed.";
    lower[SoPlex<R>::DECOMP_DISPLAYFREQ] = 1;
    upper[SoPlex<R>::DECOMP_DISPLAYFREQ] = INT_MAX;
    defaultValue[SoPlex<R>::DECOMP_DISPLAYFREQ] = 50;

    // the verbosity of the decomposition based simplex
    name[SoPlex<R>::DECOMP_VERBOSITY] = "decomp_verbosity";
    description[SoPlex<R>::DECOMP_VERBOSITY] = "the verbosity of decomposition based simplex (0 - error, 1 - warning, 2 - debug, 3 - normal, 4 - high, 5 - full).";
    lower[SoPlex<R>::DECOMP_VERBOSITY] = 1;
    upper[SoPlex<R>::DECOMP_VERBOSITY] = 5;
    defaultValue[SoPlex<R>::DECOMP_VERBOSITY] = VERBOSITY_ERROR;

    // printing condition number during the solve
    name[SoPlex<R>::PRINTCONDITION] = "printcondition";
    description[SoPlex<R>::PRINTCONDITION] = "print condition number during the solve (0 - off, 1 - ratio estimate , 2 - sum estimate, 3 - product estimate, 4 - exact)";
    lower[SoPlex<R>::PRINTCONDITION] = 0;
    upper[SoPlex<R>::PRINTCONDITION] = 4;
    defaultValue[SoPlex<R>::PRINTCONDITION] = 0;
  }

  template <class R>
  SoPlex<R>::Settings::RealParam::RealParam() {
    // primal feasibility tolerance
    name[SoPlex<R>::FEASTOL] = "feastol";
    description[SoPlex<R>::FEASTOL] = "primal feasibility tolerance";
    lower[SoPlex<R>::FEASTOL] = 0.0;
    upper[SoPlex<R>::FEASTOL] = 1.0;
    defaultValue[SoPlex<R>::FEASTOL] = 1e-6;

    // dual feasibility tolerance
    name[SoPlex<R>::OPTTOL] = "opttol";
    description[SoPlex<R>::OPTTOL] = "dual feasibility tolerance";
    lower[SoPlex<R>::OPTTOL] = 0.0;
    upper[SoPlex<R>::OPTTOL] = 1.0;
    defaultValue[SoPlex<R>::OPTTOL] = 1e-6;

    ///@todo define suitable values depending on Real type
    // general zero tolerance
    name[SoPlex<R>::EPSILON_ZERO] = "epsilon_zero";
    description[SoPlex<R>::EPSILON_ZERO] = "general zero tolerance";
    lower[SoPlex<R>::EPSILON_ZERO] = 0.0;
    upper[SoPlex<R>::EPSILON_ZERO] = 1.0;
    defaultValue[SoPlex<R>::EPSILON_ZERO] = DEFAULT_EPS_ZERO;

    ///@todo define suitable values depending on Real type
    // zero tolerance used in factorization
    name[SoPlex<R>::EPSILON_FACTORIZATION] = "epsilon_factorization";
    description[SoPlex<R>::EPSILON_FACTORIZATION] = "zero tolerance used in factorization";
    lower[SoPlex<R>::EPSILON_FACTORIZATION] = 0.0;
    upper[SoPlex<R>::EPSILON_FACTORIZATION] = 1.0;
    defaultValue[SoPlex<R>::EPSILON_FACTORIZATION] = DEFAULT_EPS_FACTOR;

    ///@todo define suitable values depending on Real type
    // zero tolerance used in update of the factorization
    name[SoPlex<R>::EPSILON_UPDATE] = "epsilon_update";
    description[SoPlex<R>::EPSILON_UPDATE] = "zero tolerance used in update of the factorization";
    lower[SoPlex<R>::EPSILON_UPDATE] = 0.0;
    upper[SoPlex<R>::EPSILON_UPDATE] = 1.0;
    defaultValue[SoPlex<R>::EPSILON_UPDATE] = DEFAULT_EPS_UPDATE;

    ///@todo define suitable values depending on Real type
    // pivot zero tolerance used in factorization
    name[SoPlex<R>::EPSILON_PIVOT] = "epsilon_pivot";
    description[SoPlex<R>::EPSILON_PIVOT] = "pivot zero tolerance used in factorization";
    lower[SoPlex<R>::EPSILON_PIVOT] = 0.0;
    upper[SoPlex<R>::EPSILON_PIVOT] = 1.0;
    defaultValue[SoPlex<R>::EPSILON_PIVOT] = DEFAULT_EPS_PIVOT;

    ///@todo define suitable values depending on Real type
    // infinity threshold
    name[SoPlex<R>::INFTY] = "infty";
    description[SoPlex<R>::INFTY] = "infinity threshold";
    lower[SoPlex<R>::INFTY] = 1e10;
    upper[SoPlex<R>::INFTY] = 1e100;
    defaultValue[SoPlex<R>::INFTY] = DEFAULT_INFINITY;

    // time limit in seconds (INFTY if unlimited)
    name[SoPlex<R>::TIMELIMIT] = "timelimit";
    description[SoPlex<R>::TIMELIMIT] = "time limit in seconds";
    lower[SoPlex<R>::TIMELIMIT] = 0.0;
    upper[SoPlex<R>::TIMELIMIT] = DEFAULT_INFINITY;
    defaultValue[SoPlex<R>::TIMELIMIT] = DEFAULT_INFINITY;

    // lower limit on objective value
    name[SoPlex<R>::OBJLIMIT_LOWER] = "objlimit_lower";
    description[SoPlex<R>::OBJLIMIT_LOWER] = "lower limit on objective value";
    lower[SoPlex<R>::OBJLIMIT_LOWER] = -DEFAULT_INFINITY;
    upper[SoPlex<R>::OBJLIMIT_LOWER] = DEFAULT_INFINITY;
    defaultValue[SoPlex<R>::OBJLIMIT_LOWER] = -DEFAULT_INFINITY;

    // upper limit on objective value
    name[SoPlex<R>::OBJLIMIT_UPPER] = "objlimit_upper";
    description[SoPlex<R>::OBJLIMIT_UPPER] = "upper limit on objective value";
    lower[SoPlex<R>::OBJLIMIT_UPPER] = -DEFAULT_INFINITY;
    upper[SoPlex<R>::OBJLIMIT_UPPER] = DEFAULT_INFINITY;
    defaultValue[SoPlex<R>::OBJLIMIT_UPPER] = DEFAULT_INFINITY;

    // working tolerance for feasibility in floating-point solver during iterative refinement
    name[SoPlex<R>::FPFEASTOL] = "fpfeastol";
    description[SoPlex<R>::FPFEASTOL] = "working tolerance for feasibility in floating-point solver during iterative refinement";
    lower[SoPlex<R>::FPFEASTOL] = 1e-12;
    upper[SoPlex<R>::FPFEASTOL] = 1.0;
    defaultValue[SoPlex<R>::FPFEASTOL] = 1e-9;

    // working tolerance for optimality in floating-point solver during iterative refinement
    name[SoPlex<R>::FPOPTTOL] = "fpopttol";
    description[SoPlex<R>::FPOPTTOL] = "working tolerance for optimality in floating-point solver during iterative refinement";
    lower[SoPlex<R>::FPOPTTOL] = 1e-12;
    upper[SoPlex<R>::FPOPTTOL] = 1.0;
    defaultValue[SoPlex<R>::FPOPTTOL] = 1e-9;

    // maximum increase of scaling factors between refinements
    name[SoPlex<R>::MAXSCALEINCR] = "maxscaleincr";
    description[SoPlex<R>::MAXSCALEINCR] = "maximum increase of scaling factors between refinements";
    lower[SoPlex<R>::MAXSCALEINCR] = 1.0;
    upper[SoPlex<R>::MAXSCALEINCR] = DEFAULT_INFINITY;
    defaultValue[SoPlex<R>::MAXSCALEINCR] = 1e25;

    // lower threshold in lifting (nonzero matrix coefficients with smaller absolute value will be reformulated)
    name[SoPlex<R>::LIFTMINVAL] = "liftminval";
    description[SoPlex<R>::LIFTMINVAL] = "lower threshold in lifting (nonzero matrix coefficients with smaller absolute value will be reformulated)";
    lower[SoPlex<R>::LIFTMINVAL] = 0.0;
    upper[SoPlex<R>::LIFTMINVAL] = 0.1;
    defaultValue[SoPlex<R>::LIFTMINVAL] = 0.000976562; // = 1/1024

    // upper threshold in lifting (nonzero matrix coefficients with larger absolute value will be reformulated)
    name[SoPlex<R>::LIFTMAXVAL] = "liftmaxval";
    description[SoPlex<R>::LIFTMAXVAL] = "lower threshold in lifting (nonzero matrix coefficients with smaller absolute value will be reformulated)";
    lower[SoPlex<R>::LIFTMAXVAL] = 10.0;
    upper[SoPlex<R>::LIFTMAXVAL] = DEFAULT_INFINITY;
    defaultValue[SoPlex<R>::LIFTMAXVAL] = 1024.0;

    // threshold for using sparse pricing (no. of violations need to be smaller than threshold * dimension of problem)
    name[SoPlex<R>::SPARSITY_THRESHOLD] = "sparsity_threshold";
    description[SoPlex<R>::SPARSITY_THRESHOLD] = "sparse pricing threshold (#violations < dimension * SPARSITY_THRESHOLD activates sparse pricing)";
    lower[SoPlex<R>::SPARSITY_THRESHOLD] = 0.0;
    upper[SoPlex<R>::SPARSITY_THRESHOLD] = 1.0;
    defaultValue[SoPlex<R>::SPARSITY_THRESHOLD] = 0.6;

    // threshold on number of rows vs. number of columns for switching from column to row representations in auto mode
    name[SoPlex<R>::REPRESENTATION_SWITCH] = "representation_switch";
    description[SoPlex<R>::REPRESENTATION_SWITCH] = "threshold on number of rows vs. number of columns for switching from column to row representations in auto mode";
    lower[SoPlex<R>::REPRESENTATION_SWITCH] = 0.0;
    upper[SoPlex<R>::REPRESENTATION_SWITCH] = DEFAULT_INFINITY;
    defaultValue[SoPlex<R>::REPRESENTATION_SWITCH] = 1.2;

    // geometric frequency at which to apply rational reconstruction
    name[SoPlex<R>::RATREC_FREQ] = "ratrec_freq";
    description[SoPlex<R>::RATREC_FREQ] = "geometric frequency at which to apply rational reconstruction";
    lower[SoPlex<R>::RATREC_FREQ] = 1.0;
    upper[SoPlex<R>::RATREC_FREQ] = DEFAULT_INFINITY;
    defaultValue[SoPlex<R>::RATREC_FREQ] = 1.2;

    // minimal reduction (sum of removed rows/cols) to continue simplification
    name[SoPlex<R>::MINRED] = "minred";
    description[SoPlex<R>::MINRED] = "minimal reduction (sum of removed rows/cols) to continue simplification";
    lower[SoPlex<R>::MINRED] = 0.0;
    upper[SoPlex<R>::MINRED] = 1.0;
    defaultValue[SoPlex<R>::MINRED] = 1e-4;

    // refactor threshold for nonzeros in last factorized basis matrix compared to updated basis matrix
    name[SoPlex<R>::REFAC_BASIS_NNZ] = "refac_basis_nnz";
    description[SoPlex<R>::REFAC_BASIS_NNZ] = "refactor threshold for nonzeros in last factorized basis matrix compared to updated basis matrix";
    lower[SoPlex<R>::REFAC_BASIS_NNZ] = 1.0;
    upper[SoPlex<R>::REFAC_BASIS_NNZ] = 100.0;
    defaultValue[SoPlex<R>::REFAC_BASIS_NNZ] = 10.0;

    // refactor threshold for fill-in in current factor update compared to fill-in in last factorization
    name[SoPlex<R>::REFAC_UPDATE_FILL] = "refac_update_fill";
    description[SoPlex<R>::REFAC_UPDATE_FILL] = "refactor threshold for fill-in in current factor update compared to fill-in in last factorization";
    lower[SoPlex<R>::REFAC_UPDATE_FILL] = 1.0;
    upper[SoPlex<R>::REFAC_UPDATE_FILL] = 100.0;
    defaultValue[SoPlex<R>::REFAC_UPDATE_FILL] = 5.0;

    // refactor threshold for memory growth in factorization since last refactorization
    name[SoPlex<R>::REFAC_MEM_FACTOR] = "refac_mem_factor";
    description[SoPlex<R>::REFAC_MEM_FACTOR] = "refactor threshold for memory growth in factorization since last refactorization";
    lower[SoPlex<R>::REFAC_MEM_FACTOR] = 1.0;
    upper[SoPlex<R>::REFAC_MEM_FACTOR] = 10.0;
    defaultValue[SoPlex<R>::REFAC_MEM_FACTOR] = 1.5;

    // accuracy of conjugate gradient method in least squares scaling (higher value leads to more iterations)
    name[SoPlex<R>::LEASTSQ_ACRCY] = "leastsq_acrcy";
    description[SoPlex<R>::LEASTSQ_ACRCY] = "accuracy of conjugate gradient method in least squares scaling (higher value leads to more iterations)";
    lower[SoPlex<R>::LEASTSQ_ACRCY] = 1.0;
    upper[SoPlex<R>::LEASTSQ_ACRCY] = DEFAULT_INFINITY;
    defaultValue[SoPlex<R>::LEASTSQ_ACRCY] = 1000.0;

    // objective offset
    name[SoPlex<R>::OBJ_OFFSET] = "obj_offset";
    description[SoPlex<R>::OBJ_OFFSET] = "objective offset to be used";
    lower[SoPlex<R>::OBJ_OFFSET] = -DEFAULT_INFINITY;
    upper[SoPlex<R>::OBJ_OFFSET] = DEFAULT_INFINITY;
    defaultValue[SoPlex<R>::OBJ_OFFSET] = 0.0;
  }

#ifdef SOPLEX_WITH_RATIONALPARAM
  SoPlex<R>::Settings::RationalParam::RationalParam() {}
#endif

  /// returns boolean parameter value
  /// #template #temp
  template <>
	bool SoPlex<Real>::boolParam(const BoolParam param) const
  {
    assert(param >= 0);
    assert(param < SoPlex<Real>::BOOLPARAM_COUNT);
    return _currentSettings->_boolParamValues[param];
  }

  template <>
  SoPlex<Real>::Settings::Settings()
  {
    for( int i = 0; i < SoPlex<Real>::BOOLPARAM_COUNT; i++ )
      _boolParamValues[i] = boolParam.defaultValue[i];

    for( int i = 0; i < SoPlex<Real>::INTPARAM_COUNT; i++ )
      _intParamValues[i] = intParam.defaultValue[i];

    for( int i = 0; i < SoPlex<Real>::REALPARAM_COUNT; i++ )
      _realParamValues[i] = realParam.defaultValue[i];

#ifdef SOPLEX_WITH_RATIONALPARAM
    for( int i = 0; i < SoPlex<Real>::RATIONALPARAM_COUNT; i++ )
      _rationalParamValues[i] = rationalParam.defaultValue[i];
#endif
  }

  template <>
  SoPlex<Real>::Settings::Settings(const Settings& settings)
  {
    *this = settings;
  }

  template <>
  typename SoPlex<Real>::Settings& SoPlex<Real>::Settings::operator=(const Settings& settings)
  {
    for( int i = 0; i < SoPlex<Real>::BOOLPARAM_COUNT; i++ )
      _boolParamValues[i] = settings._boolParamValues[i];

    for( int i = 0; i < SoPlex<Real>::INTPARAM_COUNT; i++ )
      _intParamValues[i] = settings._intParamValues[i];

    for( int i = 0; i < SoPlex<Real>::REALPARAM_COUNT; i++ )
      _realParamValues[i] = settings._realParamValues[i];

#ifdef SOPLEX_WITH_RATIONALPARAM
    for( int i = 0; i < SoPlex<Real>::RATIONALPARAM_COUNT; i++ )
      _rationalParamValues[i] = settings._rationalParamValues[i];
#endif

    return *this;
  }

  template <>
  typename SoPlex<Real>::Settings::BoolParam SoPlex<Real>::Settings::boolParam = BoolParam();
  template <>
  typename SoPlex<Real>::Settings::IntParam SoPlex<Real>::Settings::intParam = IntParam();
  template <>
  typename SoPlex<Real>::Settings::RealParam SoPlex<Real>::Settings::realParam = RealParam();
#ifdef SOPLEX_WITH_RATIONALPARAM
  template <>
  SoPlex<Real>::Settings::RationalParam SoPlex<Real>::Settings::rationalParam;
#endif

  /// default constructor
  template <>
  SoPlex<Real>::SoPlex()
    : _statistics(0)
    , _currentSettings(0)
    , _scalerUniequi(false)
    , _scalerBiequi(true)
    , _scalerGeo1(false, 1)
    , _scalerGeo8(false, 8)
    , _scalerGeoequi(true)
    , _scalerLeastsq()
    , _simplifier(0)
    , _scaler(0)
    , _starter(0)
    , _rationalLP(0)
    , _unitMatrixRational(0)
    , _status(SPxSolver<Real>::UNKNOWN)
    , _hasBasis(false)
    , _hasSolReal(false)
    , _hasSolRational(false)
    , _rationalPosone(1)
    , _rationalNegone(-1)
    , _rationalZero(0)
  {
    // transfer message handler
    _solver.setOutstream(spxout);
    _scalerUniequi.setOutstream(spxout);
    _scalerBiequi.setOutstream(spxout);
    _scalerGeo1.setOutstream(spxout);
    _scalerGeo8.setOutstream(spxout);
    _scalerGeoequi.setOutstream(spxout);
    _scalerLeastsq.setOutstream(spxout);

    // give lu factorization to solver
    _solver.setBasisSolver(&_slufactor);

    // the real LP is initially stored in the solver; the rational LP is constructed, when the parameter SYNCMODE is
    // initialized in setSettings() below
    _realLP = &_solver;
    _isRealLPLoaded = true;
    _isRealLPScaled = false;
    _applyPolishing = false;
    _optimizeCalls = 0;
    _unscaleCalls = 0;
    _realLP->setOutstream(spxout);
    _currentProb = DECOMP_ORIG;

    // initialize statistics
    spx_alloc(_statistics);
    _statistics = new (_statistics) Statistics();

    // initialize parameter settings to default
    spx_alloc(_currentSettings);
    _currentSettings = new (_currentSettings) Settings();
    setSettings(*_currentSettings, true);

    _lastSolveMode = intParam(SoPlex<Real>::SOLVEMODE);

    assert(_isConsistent());
  }

  /// returns the current solver status
  /// #template #temp 
  template <>
  typename SPxSolver<Real>::Status SoPlex<Real>::status() const
  {
    return _status;
  }

  /// assignment operator
  template <>
  SoPlex<Real>& SoPlex<Real>::operator=(const SoPlex<Real>& rhs)
  {
    if( this != &rhs )
      {
        // copy message handler
        spxout = rhs.spxout;

        // copy statistics
        *_statistics = *(rhs._statistics);

        // copy settings
        *_currentSettings = *(rhs._currentSettings);

        // copy solver components
        _solver = rhs._solver;
        _slufactor = rhs._slufactor;
        _simplifierMainSM = rhs._simplifierMainSM;
        _scalerUniequi = rhs._scalerUniequi;
        _scalerBiequi = rhs._scalerBiequi;
        _scalerGeo1 = rhs._scalerGeo1;
        _scalerGeo8 = rhs._scalerGeo8;
        _scalerGeoequi = rhs._scalerGeoequi;
        _scalerLeastsq = rhs._scalerLeastsq;
        _starterWeight = rhs._starterWeight;
        _starterSum = rhs._starterSum;
        _starterVector = rhs._starterVector;
        _pricerAuto = rhs._pricerAuto;
        _pricerDantzig = rhs._pricerDantzig;
        _pricerParMult = rhs._pricerParMult;
        _pricerDevex = rhs._pricerDevex;
        _pricerQuickSteep = rhs._pricerQuickSteep;
        _pricerSteep = rhs._pricerSteep;
        _ratiotesterTextbook = rhs._ratiotesterTextbook;
        _ratiotesterHarris = rhs._ratiotesterHarris;
        _ratiotesterFast = rhs._ratiotesterFast;
        _ratiotesterBoundFlipping = rhs._ratiotesterBoundFlipping;

        // copy solution data
        _status = rhs._status;
        _lastSolveMode = rhs._lastSolveMode;
        _basisStatusRows = rhs._basisStatusRows;
        _basisStatusCols = rhs._basisStatusCols;

        if( rhs._hasSolReal )
          _solReal = rhs._solReal;

        if( rhs._hasSolRational )
          _solRational = rhs._solRational;

        // set message handlers in members
        _solver.setOutstream(spxout);
        _scalerUniequi.setOutstream(spxout);
        _scalerBiequi.setOutstream(spxout);
        _scalerGeo1.setOutstream(spxout);
        _scalerGeo8.setOutstream(spxout);
        _scalerGeoequi.setOutstream(spxout);
        _scalerLeastsq.setOutstream(spxout);

        // transfer the lu solver
        _solver.setBasisSolver(&_slufactor);

        // initialize pointers for simplifier, scaler, and starter
        setIntParam(SoPlex<Real>::SIMPLIFIER, intParam(SoPlex<Real>::SIMPLIFIER), true);
        setIntParam(SoPlex<Real>::SCALER, intParam(SoPlex<Real>::SCALER), true);
        setIntParam(SoPlex<Real>::STARTER, intParam(SoPlex<Real>::STARTER), true);

        // copy real LP if different from the LP in the solver
        if( rhs._realLP != &(rhs._solver) )
          {
            _realLP = 0;
            spx_alloc(_realLP);
            _realLP = new (_realLP) SPxLPReal(*(rhs._realLP));
          }
        else
          _realLP = &_solver;

        // copy rational LP
        if( rhs._rationalLP == 0 )
          {
            assert(intParam(SoPlex<Real>::SYNCMODE) == SYNCMODE_ONLYREAL);
            _rationalLP = 0;
          }
        else
          {
            assert(intParam(SoPlex<Real>::SYNCMODE) != SYNCMODE_ONLYREAL);
            _rationalLP = 0;
            spx_alloc(_rationalLP);
            _rationalLP = new (_rationalLP) SPxLPRational(*rhs._rationalLP);
          }

        // copy rational factorization
        if( rhs._rationalLUSolver.status() == SLinSolverRational::OK )
          _rationalLUSolver = rhs._rationalLUSolver;

        // copy boolean flags
        _isRealLPLoaded = rhs._isRealLPLoaded;
        _isRealLPScaled = rhs._isRealLPScaled;
        _hasSolReal = rhs._hasSolReal;
        _hasSolRational = rhs._hasSolRational;
        _hasBasis = rhs._hasBasis;
        _applyPolishing = rhs._applyPolishing;

        // rational constants do not need to be assigned
        _rationalPosone = 1;
        _rationalNegone = -1;
        _rationalZero = 0;
      }

    assert(_isConsistent());

    return *this;
  }



  /// copy constructor
  ///@todo improve performance by implementing a separate copy constructor
  template <>
  SoPlex<Real>::SoPlex(const SoPlex<Real>& rhs)
  {
    // allocate memory as in default constructor
    _statistics = 0;
    spx_alloc(_statistics);
    _statistics = new (_statistics) Statistics();

    _currentSettings = 0;
    spx_alloc(_currentSettings);
    _currentSettings = new (_currentSettings) Settings();

    // call assignment operator
    *this = rhs;
  }



  /// destructor
  template <>
  SoPlex<Real>::~SoPlex()
  {
    assert(_isConsistent());

    // free settings
    _currentSettings->~Settings();
    spx_free(_currentSettings);

    // free statistics
    _statistics->~Statistics();
    spx_free(_statistics);

    // free real LP if different from the LP in the solver
    assert(_realLP != 0);
    if( _realLP != &_solver )
      {
        _realLP->~SPxLPReal();
        spx_free(_realLP);
      }

    // free rational LP
    if( _rationalLP != 0 )
      {
        _rationalLP->~SPxLPRational();
        spx_free(_rationalLP);
      }

    // free unit vectors
    for( int i = 0; i < _unitMatrixRational.size(); i++ )
      {
        if( _unitMatrixRational[i] != 0 )
          {
            _unitMatrixRational[i]->~UnitVectorRational();
            spx_free(_unitMatrixRational[i]);
          }
      }
  }



  /// returns number of rows #template #temp
  template <>
  int SoPlex<Real>::numRowsT() const
  {
    assert(_realLP != 0);
    return _realLP->nRows();
  }

  /// returns number of columns
  template <>
  int SoPlex<Real>::numColsT() const
  {
    assert(_realLP != 0);
    return _realLP->nCols();
  }



  /// returns number of nonzeros

  template <>
  int SoPlex<Real>::numNonzerosT() const
  {
    assert(_realLP != 0);
    return _realLP->nNzos();
  }



  /// returns smallest non-zero element in absolute value
  template <>
  Real SoPlex<Real>::minAbsNonzeroReal() const
  {
    assert(_realLP != 0);
    return _realLP->minAbsNzo();
  }



  /// returns biggest non-zero element in absolute value
  template <>
  Real SoPlex<Real>::maxAbsNonzeroReal() const
  {
    assert(_realLP != 0);
    return _realLP->maxAbsNzo();
  }



  /// returns (unscaled) coefficient
  template <>
  Real SoPlex<Real>::coefReal(int row, int col) const
  {
    if( _realLP->isScaled() )
      {
        assert(_scaler);
        return _scaler->getCoefUnscaled(*_realLP, row, col);
      }
    else
      return colVectorRealInternal(col)[row];
  }



  /// returns vector of row \p i, ignoring scaling
  template <>
  const SVectorReal& SoPlex<Real>::rowVectorRealInternal(int i) const
  {
    assert(_realLP != 0);
    return _realLP->rowVector(i);
  }



  /// gets vector of row \p i
  template <>
  void SoPlex<Real>::getRowVectorReal(int i, DSVectorReal& row) const
  {
    assert(_realLP);

    if( _realLP->isScaled() )
      {
        assert(_scaler);
        row.setMax(_realLP->rowVector(i).size());
        _scaler->getRowUnscaled(*_realLP, i, row);
      }
    else
      row = _realLP->rowVector(i);
  }



  /// returns right-hand side vector, ignoring scaling
  template <>
  const VectorReal& SoPlex<Real>::rhsRealInternal() const
  {
    assert(_realLP != 0);
    return _realLP->rhs();
  }



  /// gets right-hand side vector
  template <>
  void SoPlex<Real>::getRhsReal(DVectorReal& rhs) const
  {
    assert(_realLP);
    _realLP->getRhsUnscaled(rhs);
  }



  /// returns right-hand side of row \p i
  template <>
  Real SoPlex<Real>::rhsReal(int i) const
  {
    assert(_realLP != 0);
    return _realLP->rhsUnscaled(i);
  }



  /// returns left-hand side vector, ignoring scaling
  template <>
  const VectorReal& SoPlex<Real>::lhsRealInternal() const
  {
    assert(_realLP != 0);
    return _realLP->lhs();
  }



  /// gets left-hand side vector
  template <>
  void SoPlex<Real>::getLhsReal(DVectorReal& lhs) const
  {
    assert(_realLP);
    _realLP->getLhsUnscaled(lhs);
  }



  /// returns left-hand side of row \p i
  template <>
  Real SoPlex<Real>::lhsReal(int i) const
  {
    assert(_realLP != 0);
    return _realLP->lhsUnscaled(i);
  }



  /// returns inequality type of row \p i
  template <>
  LPRowReal::Type SoPlex<Real>::rowTypeReal(int i) const
  {
    assert(_realLP != 0);
    return _realLP->rowType(i);
  }



  /// returns vector of col \p i, ignoring scaling
  template <>
  const SVectorReal& SoPlex<Real>::colVectorRealInternal(int i) const
  {
    assert(_realLP != 0);
    return _realLP->colVector(i);
  }



  /// gets vector of col \p i
  template <>
  void SoPlex<Real>::getColVectorReal(int i, DSVectorReal& col) const
  {
    assert(_realLP);
    _realLP->getColVectorUnscaled(i, col);
  }



  /// returns upper bound vector
  template <>
  const VectorReal& SoPlex<Real>::upperRealInternal() const
  {
    assert(_realLP != 0);
    return _realLP->upper();
  }



  /// returns upper bound of column \p i
  template <>
  Real SoPlex<Real>::upperReal(int i) const
  {
    assert(_realLP != 0);
    return _realLP->upperUnscaled(i);
  }



  /// gets upper bound vector
  template <>
  void SoPlex<Real>::getUpperReal(DVectorReal& upper) const
  {
    assert(_realLP != 0);
    return _realLP->getUpperUnscaled(upper);
  }



  /// returns lower bound vector
  template <>
  const VectorReal& SoPlex<Real>::lowerRealInternal() const
  {
    assert(_realLP != 0);
    return _realLP->lower();
  }



  /// returns lower bound of column \p i
  template <>
  Real SoPlex<Real>::lowerReal(int i) const
  {
    assert(_realLP != 0);
    return _realLP->lowerUnscaled(i);
  }



  /// gets lower bound vector
  template <>
  void SoPlex<Real>::getLowerReal(DVectorReal& lower) const
  {
    assert(_realLP != 0);
    return _realLP->getLowerUnscaled(lower);
  }




  /// gets objective function vector
  template <>
  void SoPlex<Real>::getObjReal(VectorReal& obj) const
  {
    assert(_realLP != 0);
    _realLP->getObjUnscaled(obj);
  }



  /// returns objective value of column \p i
  template <>
  Real SoPlex<Real>::objReal(int i) const
  {
    assert(_realLP != 0);
    return _realLP->objUnscaled(i);
  }



  /// returns objective function vector after transformation to a maximization problem; since this is how it is stored
  /// internally, this is generally faster
  template <>
  const VectorReal& SoPlex<Real>::maxObjRealInternal() const
  {
    assert(_realLP != 0);
    return _realLP->maxObj();
  }



  /// returns objective value of column \p i after transformation to a maximization problem; since this is how it is
  /// stored internally, this is generally faster
  template <>
  Real SoPlex<Real>::maxObjReal(int i) const
  {
    assert(_realLP != 0);
    return _realLP->maxObjUnscaled(i);
  }



  /// gets number of available dual norms
  template <>
  void SoPlex<Real>::getNdualNorms(int& nnormsRow, int& nnormsCol) const
  {
    _solver.getNdualNorms(nnormsRow, nnormsCol);
  }



  /// gets steepest edge norms and returns false if they are not available
  template <>
  bool SoPlex<Real>::getDualNorms(int& nnormsRow, int& nnormsCol, Real* norms) const
  {
    return _solver.getDualNorms(nnormsRow, nnormsCol, norms);
  }



  /// sets steepest edge norms and returns false if that's not possible
  template <>
  bool SoPlex<Real>::setDualNorms(int nnormsRow, int nnormsCol, Real* norms)
  {
    return _solver.setDualNorms(nnormsRow, nnormsCol, norms);
  }



  /// pass integrality information about the variables to the solver
  template <>
  void SoPlex<Real>::setIntegralityInformation( int ncols, int* intInfo)
  {
    assert(ncols == _solver.nCols() || (ncols == 0 && intInfo == NULL));
    _solver.setIntegralityInformation(ncols, intInfo);
  }



  /// returns number of rows
  template <>
  int SoPlex<Rational>::numRowsT() const
  {
    assert(_rationalLP != 0);
    return _rationalLP->nRows();
  }

  template <>
  int SoPlex<Real>::numRowsRational() const
  {
    assert(_rationalLP != 0);
    return _rationalLP->nRows();
  }

  /// returns number of columns
  template <>
  int SoPlex<Rational>::numColsT() const
  {
    assert(_rationalLP != 0);
    return _rationalLP->nCols();
  }

  template <>
  int SoPlex<Real>::numColsRational() const
  {
    assert(_rationalLP != 0);
    return _rationalLP->nCols();
  }

  /// returns number of nonzeros
  template <>
  int SoPlex<Rational>::numNonzerosT() const
  {
    assert(_rationalLP != 0);
    return _rationalLP->nNzos();
  }



  /// returns smallest non-zero element in absolute value
  template <>
  Rational SoPlex<Real>::minAbsNonzeroRational() const
  {
    assert(_rationalLP != 0);
    return _rationalLP->minAbsNzo();
  }



  /// returns biggest non-zero element in absolute value
  template <>
  Rational SoPlex<Real>::maxAbsNonzeroRational() const
  {
    assert(_rationalLP != 0);
    return _rationalLP->maxAbsNzo();
  }



  /// gets row \p i
  template <>
  void SoPlex<Real>::getRowRational(int i, LPRowRational& lprow) const
  {
    assert(_rationalLP != 0);
    _rationalLP->getRow(i, lprow);
  }



  /// gets rows \p start, ..., \p end.
  template <>
  void SoPlex<Real>::getRowsRational(int start, int end, LPRowSetRational& lprowset) const
  {
    assert(_rationalLP != 0);
    _rationalLP->getRows(start, end, lprowset);
  }



  /// returns vector of row \p i
  template <>
  const SVectorRational& SoPlex<Real>::rowVectorRational(int i) const
  {
    assert(_rationalLP != 0);
    return _rationalLP->rowVector(i);
  }



  /// returns right-hand side vector
  template <>
  const VectorRational& SoPlex<Real>::rhsRational() const
  {
    assert(_rationalLP != 0);
    return _rationalLP->rhs();
  }



  /// returns right-hand side of row \p i
  template <>
  const Rational& SoPlex<Real>::rhsRational(int i) const
  {
    assert(_rationalLP != 0);
    return _rationalLP->rhs(i);
  }



  /// returns left-hand side vector
  template <>
  const VectorRational& SoPlex<Real>::lhsRational() const
  {
    assert(_rationalLP != 0);
    return _rationalLP->lhs();
  }



  /// returns left-hand side of row \p i
  template <>
  const Rational& SoPlex<Real>::lhsRational(int i) const
  {
    assert(_rationalLP != 0);
    return _rationalLP->lhs(i);
  }



  /// returns inequality type of row \p i
  template <>
  LPRowRational::Type SoPlex<Real>::rowTypeRational(int i) const
  {
    assert(_rationalLP != 0);
    return _rationalLP->rowType(i);
  }



  /// gets column \p i
  template <>
  void SoPlex<Real>::getColRational(int i, LPColRational& lpcol) const
  {
    assert(_rationalLP != 0);
    return _rationalLP->getCol(i, lpcol);
  }



  /// gets columns \p start, ..., \p end
  template <>
  void SoPlex<Real>::getColsRational(int start, int end, LPColSetRational& lpcolset) const
  {
    assert(_rationalLP != 0);
    return _rationalLP->getCols(start, end, lpcolset);
  }



  /// returns vector of column \p i
  template <>
  const SVectorRational& SoPlex<Real>::colVectorRational(int i) const
  {
    assert(_rationalLP != 0);
    return _rationalLP->colVector(i);
  }



  /// returns upper bound vector
  template <>
  const VectorRational& SoPlex<Real>::upperRational() const
  {
    assert(_rationalLP != 0);
    return _rationalLP->upper();
  }



  /// returns upper bound of column \p i
  template <>
  const Rational& SoPlex<Real>::upperRational(int i) const
  {
    assert(_rationalLP != 0);
    return _rationalLP->upper(i);
  }



  /// returns lower bound vector
  template <>
  const VectorRational& SoPlex<Real>::lowerRational() const
  {
    assert(_rationalLP != 0);
    return _rationalLP->lower();
  }

  /// returns integer parameter value
  // #template #baseclass
  template <>
	int SoPlex<Real>::intParam(const IntParam param) const
  {
    assert(param >= 0);
    assert(param < INTPARAM_COUNT);
    return _currentSettings->_intParamValues[param];
  }

  template <>
	int SoPlex<Rational>::intParam(const IntParam param) const
  {
    assert(param >= 0);
    assert(param < INTPARAM_COUNT);
    return _currentSettings->_intParamValues[param];
  }

  /// sets integer parameter value; returns true on success
  // #template #baseclass
  template <>
	bool SoPlex<Real>::setIntParam(const IntParam param, const int value, const bool init)
  {
    assert(param >= 0);
    assert(param < INTPARAM_COUNT);
    assert(init || _isConsistent());

    if( !init && value == intParam(param) )
      return true;

    // check for a valid parameter value wrt bounds
    if( value < _currentSettings->intParam.lower[param] || value > _currentSettings->intParam.upper[param] )
      return false;

    switch( param )
      {
        // objective sense
      case SoPlex<Real>::OBJSENSE:
        if( value != SoPlex<Real>::OBJSENSE_MAXIMIZE && value != SoPlex<Real>::OBJSENSE_MINIMIZE )
          return false;
        _realLP->changeSense(value == SoPlex<Real>::OBJSENSE_MAXIMIZE ? SPxLPReal::MAXIMIZE : SPxLPReal::MINIMIZE);
        if( _rationalLP != 0 )
          _rationalLP->changeSense(value == SoPlex<Real>::OBJSENSE_MAXIMIZE ? SPxLPRational::MAXIMIZE : SPxLPRational::MINIMIZE);
        _invalidateSolution();
        break;

        // type of computational form, i.e., column or row representation
      case SoPlex<Real>::REPRESENTATION:
        if( value != SoPlex<Real>::REPRESENTATION_COLUMN && value != SoPlex<Real>::REPRESENTATION_ROW && value != SoPlex<Real>::REPRESENTATION_AUTO )
          return false;
        break;

        // type of algorithm, i.e., primal or dual
      case SoPlex<Real>::ALGORITHM:
        // decide upon entering/leaving at solve time depending on representation
        break;

        // type of LU update
      case SoPlex<Real>::FACTOR_UPDATE_TYPE:
        if( value != SoPlex<Real>::FACTOR_UPDATE_TYPE_ETA && value != SoPlex<Real>::FACTOR_UPDATE_TYPE_FT )
          return false;
        _slufactor.setUtype(value == SoPlex<Real>::FACTOR_UPDATE_TYPE_ETA ? SLUFactor::ETA : SLUFactor::FOREST_TOMLIN);
        break;

        // maximum number of updates before fresh factorization
      case SoPlex<Real>::FACTOR_UPDATE_MAX:
        if( value == 0 )
          _solver.basis().setMaxUpdates(DEFAULT_REFACTOR_INTERVAL);
        else
          _solver.basis().setMaxUpdates(value);
        break;

        // iteration limit (-1 if unlimited)
      case SoPlex<Real>::ITERLIMIT:
        break;

        // refinement limit (-1 if unlimited)
      case SoPlex<Real>::REFLIMIT:
        break;

        // stalling refinement limit (-1 if unlimited)
      case SoPlex<Real>::STALLREFLIMIT:
        break;

        // display frequency
      case SoPlex<Real>::DISPLAYFREQ:
        _solver.setDisplayFreq(value);
        break;

        // verbosity level
      case SoPlex<Real>::VERBOSITY:
        switch(value)
          {
          case 0:
            spxout.setVerbosity(SPxOut::ERROR);
            break;
          case 1:
            spxout.setVerbosity(SPxOut::WARNING);
            break;
          case 2:
            spxout.setVerbosity(SPxOut::DEBUG);
            break;
          case 3:
            spxout.setVerbosity(SPxOut::INFO1);
            break;
          case 4:
            spxout.setVerbosity(SPxOut::INFO2);
            break;
          case 5:
            spxout.setVerbosity(SPxOut::INFO3);
            break;
          }
        break;

        // type of simplifier
      case SoPlex<Real>::SIMPLIFIER:
        switch( value )
          {
          case SIMPLIFIER_OFF:
            _simplifier = 0;
            break;
          case SIMPLIFIER_AUTO:
            _simplifier = &_simplifierMainSM;
            assert(_simplifier != 0);
            break;
          default:
            return false;
          }
        break;

        // type of scaler
      case SoPlex<Real>::SCALER:
        switch( value )
          {
          case SCALER_OFF:
            _scaler = 0;
            break;
          case SCALER_UNIEQUI:
            _scaler = &_scalerUniequi;
            break;
          case SCALER_BIEQUI:
            _scaler = &_scalerBiequi;
            break;
          case SCALER_GEO1:
            _scaler = &_scalerGeo1;
            break;
          case SCALER_GEO8:
            _scaler = &_scalerGeo8;
            break;
          case SCALER_LEASTSQ:
            _scaler = &_scalerLeastsq;
            break;
          case SCALER_GEOEQUI:
            _scaler = &_scalerGeoequi;
            break;
          default:
            return false;
          }
        break;

        // type of starter used to create crash basis
      case SoPlex<Real>::STARTER:
        switch( value )
          {
          case STARTER_OFF:
            _starter = 0;
            break;
          case STARTER_WEIGHT:
            _starter = &_starterWeight;
            break;
          case STARTER_SUM:
            _starter = &_starterSum;
            break;
          case STARTER_VECTOR:
            _starter = &_starterVector;
            break;
          default:
            return false;
          }
        break;

        // type of pricer
      case SoPlex<Real>::PRICER:
        switch( value )
          {
          case PRICER_AUTO:
            _solver.setPricer(&_pricerAuto);
            break;
          case PRICER_DANTZIG:
            _solver.setPricer(&_pricerDantzig);
            break;
          case PRICER_PARMULT:
            _solver.setPricer(&_pricerParMult);
            break;
          case PRICER_DEVEX:
            _solver.setPricer(&_pricerDevex);
            break;
          case PRICER_QUICKSTEEP:
            _solver.setPricer(&_pricerQuickSteep);
            break;
          case PRICER_STEEP:
            _solver.setPricer(&_pricerSteep);
            break;
          default:
            return false;
          }
        break;

        // mode for synchronizing real and rational LP
      case SoPlex<Real>::SYNCMODE:
        switch( value )
          {
          case SYNCMODE_ONLYREAL:
            if( _rationalLP != 0 )
              {
                _rationalLP->~SPxLPRational();
                spx_free(_rationalLP);
              }
            break;
          case SYNCMODE_AUTO:
            if( intParam(param) == SYNCMODE_ONLYREAL )
              _syncLPRational();
            break;
          case SYNCMODE_MANUAL:
            _ensureRationalLP();
            break;
          default:
            return false;
          }
        break;

        // mode for reading LP files; nothing to do but change the value if valid
      case SoPlex<Real>::READMODE:
        switch( value )
          {
          case READMODE_REAL:
          case READMODE_RATIONAL:
            break;
          default:
            return false;
          }
        break;

        // mode for iterative refinement strategy; nothing to do but change the value if valid
      case SoPlex<Real>::SOLVEMODE:
        switch( value )
          {
          case SOLVEMODE_REAL:
          case SOLVEMODE_AUTO:
          case SOLVEMODE_RATIONAL:
            break;
          default:
            return false;
          }
        break;

        // mode for a posteriori feasibility checks; nothing to do but change the value if valid
      case SoPlex<Real>::CHECKMODE:
        switch( value )
          {
          case CHECKMODE_REAL:
          case CHECKMODE_AUTO:
          case CHECKMODE_RATIONAL:
            break;
          default:
            return false;
          }
        break;

        // type of ratio test
      case SoPlex<Real>::RATIOTESTER:
        switch( value )
          {
          case RATIOTESTER_TEXTBOOK:
            _solver.setTester(&_ratiotesterTextbook);
            break;
          case RATIOTESTER_HARRIS:
            _solver.setTester(&_ratiotesterHarris);
            break;
          case RATIOTESTER_FAST:
            _solver.setTester(&_ratiotesterFast);
            break;
          case RATIOTESTER_BOUNDFLIPPING:
            _solver.setTester(&_ratiotesterBoundFlipping);
            break;
          default:
            return false;
          }
        break;

        // type of timer
      case SoPlex<Real>::TIMER:
        switch( value )
          {
          case TIMER_OFF:
            _solver.setTiming( Timer::OFF);
            break;
          case TIMER_CPU:
            _solver.setTiming( Timer::USER_TIME );
            break;
          case TIMER_WALLCLOCK:
            _solver.setTiming( Timer::WALLCLOCK_TIME);
            break;
          default:
            return false;
          }
        break;

        // mode of hyper pricing
      case SoPlex<Real>::HYPER_PRICING:
        switch( value )
          {
          case HYPER_PRICING_OFF:
          case HYPER_PRICING_AUTO:
          case HYPER_PRICING_ON:
            break;
          default:
            return false;
          }
        break;

        // minimum number of stalling refinements since last pivot to trigger rational factorization
      case SoPlex<Real>::RATFAC_MINSTALLS:
        break;

        // maximum number of conjugate gradient iterations in least square scaling
      case SoPlex<Real>::LEASTSQ_MAXROUNDS:
        if( _scaler )
          _scaler->setIntParam(value);
        break;

        // mode of solution polishing
      case SoPlex<Real>::SOLUTION_POLISHING:
        switch( value )
          {
          case POLISHING_OFF:
            _solver.setSolutionPolishing(SPxSolver<Real>::POLISH_OFF);
            break;
          case POLISHING_INTEGRALITY:
            _solver.setSolutionPolishing(SPxSolver<Real>::POLISH_INTEGRALITY);
            break;
          case POLISHING_FRACTIONALITY:
            _solver.setSolutionPolishing(SPxSolver<Real>::POLISH_FRACTIONALITY);
            break;
          default:
            return false;
          }
        break;

        // the decomposition based simplex parameter settings
      case DECOMP_ITERLIMIT:
        break;
      case DECOMP_MAXADDEDROWS:
        break;
      case DECOMP_DISPLAYFREQ:
        break;
      case DECOMP_VERBOSITY:
        break;

        // printing of condition n
      case PRINTCONDITION:
        _solver.setConditionInformation(value);
        break;

      default:
        return false;
      }

    _currentSettings->_intParamValues[param] = value;
    return true;
  }

  /// returns lower bound of column \p i
  template <>
  const Rational& SoPlex<Real>::lowerRational(int i) const
  {
    assert(_rationalLP != 0);
    return _rationalLP->lower(i);
  }



  /// gets objective function vector
  template <>
  void SoPlex<Real>::getObjRational(VectorRational& obj) const
  {
    assert(_rationalLP != 0);
    _rationalLP->getObj(obj);
  }



  /// gets objective value of column \p i
  template <>
  void SoPlex<Real>::getObjRational(int i, Rational& obj) const
  {
    obj = maxObjRational(i);
    if( intParam(SoPlex<Real>::OBJSENSE) == SoPlex<Real>::OBJSENSE_MINIMIZE )
      obj *= -1;
  }



  /// returns objective value of column \p i
  template <>
  Rational SoPlex<Real>::objRational(int i) const
  {
    assert(_rationalLP != 0);
    return _rationalLP->obj(i);
  }



  /// returns objective function vector after transformation to a maximization problem; since this is how it is stored
  /// internally, this is generally faster
  template <>
  const VectorRational& SoPlex<Real>::maxObjRational() const
  {
    assert(_rationalLP != 0);
    return _rationalLP->maxObj();
  }



  /// returns objective value of column \p i after transformation to a maximization problem; since this is how it is
  /// stored internally, this is generally faster
  template <>
  const Rational& SoPlex<Real>::maxObjRational(int i) const
  {
    assert(_rationalLP != 0);
    return _rationalLP->maxObj(i);
  }



  /// adds a single row
  template <>
  void SoPlex<Real>::addRowReal(const LPRowReal& lprow)
  {
    assert(_realLP != 0);

    _addRowReal(lprow);

    if( intParam(SoPlex<Real>::SYNCMODE) == SYNCMODE_AUTO )
      {
        _rationalLP->addRow(lprow);
        _completeRangeTypesRational();
      }

    _invalidateSolution();
  }



  /// adds multiple rows
  template <>
  void SoPlex<Real>::addRowsReal(const LPRowSetReal& lprowset)
  {
    assert(_realLP != 0);

    _addRowsReal(lprowset);

    if( intParam(SoPlex<Real>::SYNCMODE) == SYNCMODE_AUTO )
      {
        _rationalLP->addRows(lprowset);
        _completeRangeTypesRational();
      }

    _invalidateSolution();
  }



  /// adds a single column
  template <>
  void SoPlex<Real>::addColReal(const LPColReal& lpcol)
  {
    assert(_realLP != 0);

    _addColReal(lpcol);

    if( intParam(SoPlex<Real>::SYNCMODE) == SYNCMODE_AUTO )
      {
        _rationalLP->addCol(lpcol);
        _completeRangeTypesRational();
      }

    _invalidateSolution();
  }



  /// adds multiple columns
  template <>
  void SoPlex<Real>::addColsReal(const LPColSetReal& lpcolset)
  {
    assert(_realLP != 0);

    _addColsReal(lpcolset);

    if( intParam(SoPlex<Real>::SYNCMODE) == SYNCMODE_AUTO )
      {
        _rationalLP->addCols(lpcolset);
        _completeRangeTypesRational();
      }

    _invalidateSolution();
  }



  /// replaces row \p i with \p lprow
  template <>
  void SoPlex<Real>::changeRowReal(int i, const LPRowReal& lprow)
  {
    assert(_realLP != 0);

    _changeRowReal(i, lprow);

    if( intParam(SoPlex<Real>::SYNCMODE) == SYNCMODE_AUTO )
      {
        _rationalLP->changeRow(i, lprow);
        _rowTypes[i] = _rangeTypeReal(lprow.lhs(), lprow.rhs());
        _completeRangeTypesRational();
      }

    _invalidateSolution();
  }



  /// changes left-hand side vector for constraints to \p lhs
  template <>
  void SoPlex<Real>::changeLhsReal(const VectorReal& lhs)
  {
    assert(_realLP != 0);

    _changeLhsReal(lhs);

    if( intParam(SoPlex<Real>::SYNCMODE) == SYNCMODE_AUTO )
      {
        _rationalLP->changeLhs(DVectorRational(lhs));
        for( int i = 0; i < numRowsT(); i++ )
          _rowTypes[i] = _rangeTypeRational(_rationalLP->lhs(i), _rationalLP->rhs(i));
      }

    _invalidateSolution();
  }



  /// changes left-hand side of row \p i to \p lhs
  template <>
  void SoPlex<Real>::changeLhsReal(int i, const Real& lhs)
  {
    assert(_realLP != 0);

    _changeLhsReal(i, lhs);

    if( intParam(SoPlex<Real>::SYNCMODE) == SYNCMODE_AUTO )
      {
        _rationalLP->changeLhs(i, lhs);
        _rowTypes[i] = _rangeTypeRational(_rationalLP->lhs(i), _rationalLP->rhs(i));
      }

    _invalidateSolution();
  }



  /// changes right-hand side vector to \p rhs
  template <>
  void SoPlex<Real>::changeRhsReal(const VectorReal& rhs)
  {
    assert(_realLP != 0);

    _changeRhsReal(rhs);

    if( intParam(SoPlex<Real>::SYNCMODE) == SYNCMODE_AUTO )
      {
        _rationalLP->changeRhs(DVectorRational(rhs));
        for( int i = 0; i < numRowsT(); i++ )
          _rowTypes[i] = _rangeTypeRational(_rationalLP->lhs(i), _rationalLP->rhs(i));
      }

    _invalidateSolution();
  }



  /// changes right-hand side of row \p i to \p rhs
  template <>
  void SoPlex<Real>::changeRhsReal(int i, const Real& rhs)
  {
    assert(_realLP != 0);

    _changeRhsReal(i, rhs);

    if( intParam(SoPlex<Real>::SYNCMODE) == SYNCMODE_AUTO )
      {
        _rationalLP->changeRhs(i, rhs);
        _rowTypes[i] = _rangeTypeRational(_rationalLP->lhs(i), _rationalLP->rhs(i));
      }

    _invalidateSolution();
  }



  /// changes left- and right-hand side vectors
  template <>
  void SoPlex<Real>::changeRangeReal(const VectorReal& lhs, const VectorReal& rhs)
  {
    assert(_realLP != 0);

    _changeRangeReal(lhs, rhs);

    if( intParam(SoPlex<Real>::SYNCMODE) == SYNCMODE_AUTO )
      {
        _rationalLP->changeRange(DVectorRational(lhs), DVectorRational(rhs));
        for( int i = 0; i < numRowsT(); i++ )
          _rowTypes[i] = _rangeTypeReal(lhs[i], rhs[i]);
      }

    _invalidateSolution();
  }



  /// changes left- and right-hand side of row \p i
  template <>
  void SoPlex<Real>::changeRangeReal(int i, const Real& lhs, const Real& rhs)
  {
    assert(_realLP != 0);

    _changeRangeReal(i,lhs, rhs);

    if( intParam(SoPlex<Real>::SYNCMODE) == SYNCMODE_AUTO )
      {
        _rationalLP->changeRange(i, lhs, rhs);
        _rowTypes[i] = _rangeTypeReal(lhs, rhs);
      }

    _invalidateSolution();
  }



  /// replaces column \p i with \p lpcol
  template <>
  void SoPlex<Real>::changeColReal(int i, const LPColReal& lpcol)
  {
    assert(_realLP != 0);

    _changeColReal(i, lpcol);

    if( intParam(SoPlex<Real>::SYNCMODE) == SYNCMODE_AUTO )
      {
        _rationalLP->changeCol(i, lpcol);
        _colTypes[i] = _rangeTypeReal(lpcol.lower(), lpcol.upper());
        _completeRangeTypesRational();
      }

    _invalidateSolution();
  }



  /// changes vector of lower bounds to \p lower
  template <>
  void SoPlex<Real>::changeLowerReal(const VectorReal& lower)
  {
    assert(_realLP != 0);

    _changeLowerReal(lower);

    if( intParam(SoPlex<Real>::SYNCMODE) == SYNCMODE_AUTO )
      {
        _rationalLP->changeLower(DVectorRational(lower));
        for( int i = 0; i < numColsT(); i++ )
          _colTypes[i] = _rangeTypeRational(_rationalLP->lower(i), _rationalLP->upper(i));
      }


    _invalidateSolution();
  }



  /// changes lower bound of column i to \p lower
  template <>
  void SoPlex<Real>::changeLowerReal(int i, const Real& lower)
  {
    assert(_realLP != 0);

    _changeLowerReal(i, lower);

    if( intParam(SoPlex<Real>::SYNCMODE) == SYNCMODE_AUTO )
      {
        _rationalLP->changeLower(i, lower);
        _colTypes[i] = _rangeTypeRational(_rationalLP->lower(i), _rationalLP->upper(i));
      }

    _invalidateSolution();
  }



  /// changes vector of upper bounds to \p upper
  template <>
  void SoPlex<Real>::changeUpperReal(const VectorReal& upper)
  {
    assert(_realLP != 0);

    _changeUpperReal(upper);

    if( intParam(SoPlex<Real>::SYNCMODE) == SYNCMODE_AUTO )
      {
        _rationalLP->changeUpper(DVectorRational(upper));
        for( int i = 0; i < numColsT(); i++ )
          _colTypes[i] = _rangeTypeRational(_rationalLP->lower(i), _rationalLP->upper(i));
      }

    _invalidateSolution();
  }



  /// changes \p i 'th upper bound to \p upper
  template <>
  void SoPlex<Real>::changeUpperReal(int i, const Real& upper)
  {
    assert(_realLP != 0);

    _changeUpperReal(i, upper);

    if( intParam(SoPlex<Real>::SYNCMODE) == SYNCMODE_AUTO )
      {
        _rationalLP->changeUpper(i, upper);
        _colTypes[i] = _rangeTypeRational(_rationalLP->lower(i), _rationalLP->upper(i));
      }

    _invalidateSolution();
  }



  /// changes vectors of column bounds to \p lower and \p upper
  template <>
  void SoPlex<Real>::changeBoundsReal(const VectorReal& lower, const VectorReal& upper)
  {
    assert(_realLP != 0);

    _changeBoundsReal(lower, upper);

    if( intParam(SoPlex<Real>::SYNCMODE) == SYNCMODE_AUTO )
      {
        _rationalLP->changeBounds(DVectorRational(lower), DVectorRational(upper));
        for( int i = 0; i < numColsT(); i++ )
          _colTypes[i] = _rangeTypeReal(lower[i], upper[i]);
      }

    _invalidateSolution();
  }



  /// changes bounds of column \p i to \p lower and \p upper
  template <>
  void SoPlex<Real>::changeBoundsReal(int i, const Real& lower, const Real& upper)
  {
    assert(_realLP != 0);

    _changeBoundsReal(i, lower, upper);

    if( intParam(SoPlex<Real>::SYNCMODE) == SYNCMODE_AUTO )
      {
        _rationalLP->changeBounds(i, lower, upper);
        _colTypes[i] = _rangeTypeReal(lower, upper);
      }
    _invalidateSolution();
  }



  /// changes objective function vector to \p obj
  template <>
  void SoPlex<Real>::changeObjReal(const VectorReal& obj)
  {
    assert(_realLP != 0);

    bool scale = _realLP->isScaled();
    _realLP->changeObj(obj, scale);

    if( intParam(SoPlex<Real>::SYNCMODE) == SYNCMODE_AUTO )
      _rationalLP->changeObj(DVectorRational(obj));

    _invalidateSolution();
  }



  /// changes objective coefficient of column i to \p obj
  template <>
  void SoPlex<Real>::changeObjReal(int i, const Real& obj)
  {
    assert(_realLP != 0);

    bool scale = _realLP->isScaled();
    _realLP->changeObj(i, obj, scale);

    if( intParam(SoPlex<Real>::SYNCMODE) == SYNCMODE_AUTO )
      _rationalLP->changeObj(i, obj);

    _invalidateSolution();
  }



  /// changes matrix entry in row \p i and column \p j to \p val
  template <>
  void SoPlex<Real>::changeElementReal(int i, int j, const Real& val)
  {
    assert(_realLP != 0);

    _changeElementReal(i, j, val);

    if( intParam(SoPlex<Real>::SYNCMODE) == SYNCMODE_AUTO )
      _rationalLP->changeElement(i, j, val);

    _invalidateSolution();
  }



  /// removes row \p i
  template <>
  void SoPlex<Real>::removeRowReal(int i)
  {
    assert(_realLP != 0);

    _removeRowReal(i);

    if( intParam(SoPlex<Real>::SYNCMODE) == SYNCMODE_AUTO )
      {
        _rationalLP->removeRow(i);
        // only swap elements if not the last one was removed
        if( i < _rationalLP->nRows() )
          {
            _rowTypes[i] = _rowTypes[_rationalLP->nRows()];
            assert(_rowTypes[i] == _rangeTypeRational(lhsRational(i), rhsRational(i)));
          }
        _rowTypes.reSize(_rationalLP->nRows());
      }

    _invalidateSolution();
  }



  /// removes all rows with an index \p i such that \p perm[i] < 0; upon completion, \p perm[i] >= 0 indicates the
  /// new index where row \p i has been moved to; note that \p perm must point to an array of size at least
  /// #numRowsT()
  template <>
  void SoPlex<Real>::removeRowsReal(int perm[])
  {
    assert(_realLP != 0);

    const int oldsize = numRowsT();
    _removeRowsReal(perm);

    if( intParam(SoPlex<Real>::SYNCMODE) == SYNCMODE_AUTO )
      {
        _rationalLP->removeRows(perm);
        for( int i = 0; i < oldsize; i++ )
          {
            if( perm[i] >= 0 )
              _rowTypes[perm[i]] = _rowTypes[i];
          }
        _rowTypes.reSize(_rationalLP->nRows());
        for( int i = 0; i < numRowsT(); i++ )
          {
            assert(_rowTypes[i] == _rangeTypeRational(lhsRational(i), rhsRational(i)));
          }
      }

    _invalidateSolution();
  }



  /// remove all rows with indices in array \p idx of size \p n; an array \p perm of size #numRowsT() may be passed
  /// as buffer memory
  template <>
  void SoPlex<Real>::removeRowsReal(int idx[], int n, int perm[])
  {
    if( perm == 0 )
      {
        DataArray< int > p(numRowsT());
        _idxToPerm(idx, n, p.get_ptr(), numRowsT());
        SoPlex<Real>::removeRowsReal(p.get_ptr());
      }
    else
      {
        _idxToPerm(idx, n, perm, numRowsT());
        SoPlex<Real>::removeRowsReal(perm);
      }
  }



  /// removes rows \p start to \p end including both; an array \p perm of size #numRowsT() may be passed as buffer
  /// memory
  template <>
  void SoPlex<Real>::removeRowRangeReal(int start, int end, int perm[])
  {
    if( perm == 0 )
      {
        DataArray< int > p(numRowsT());
        _rangeToPerm(start, end, p.get_ptr(), numRowsT());
        SoPlex<Real>::removeRowsReal(p.get_ptr());
      }
    else
      {
        _rangeToPerm(start, end, perm, numRowsT());
        SoPlex<Real>::removeRowsReal(perm);
      }
  }



  /// removes column i
  template <>
  void SoPlex<Real>::removeColReal(int i)
  {
    assert(_realLP != 0);

    _removeColReal(i);

    if( intParam(SoPlex<Real>::SYNCMODE) == SYNCMODE_AUTO )
      {
        _rationalLP->removeCol(i);
        // only swap elements if not the last one was removed
        if( i < _rationalLP->nCols() )
          {
            _colTypes[i] = _colTypes[_rationalLP->nCols()];
            assert(_colTypes[i] == _rangeTypeRational(lowerRational(i), upperRational(i)));
          }
        _colTypes.reSize(_rationalLP->nCols());
      }

    _invalidateSolution();
  }



  /// removes all columns with an index \p i such that \p perm[i] < 0; upon completion, \p perm[i] >= 0 indicates the
  /// new index where column \p i has been moved to; note that \p perm must point to an array of size at least
  /// #numColsReal()
  template <>
  void SoPlex<Real>::removeColsReal(int perm[])
  {
    assert(_realLP != 0);

    const int oldsize = numColsT();
    _removeColsReal(perm);

    if( intParam(SoPlex<Real>::SYNCMODE) == SYNCMODE_AUTO )
      {
        _rationalLP->removeCols(perm);
        for( int i = 0; i < oldsize; i++ )
          {
            if( perm[i] >= 0 )
              _colTypes[perm[i]] = _colTypes[i];
          }
        _colTypes.reSize(_rationalLP->nCols());
        for( int i = 0; i < numColsT(); i++ )
          {
            assert(_colTypes[i] == _rangeTypeRational(lowerRational(i), upperRational(i)));
          }
      }

    _invalidateSolution();
  }



  /// remove all columns with indices in array \p idx of size \p n; an array \p perm of size #numColsReal() may be
  /// passed as buffer memory
  template <>
  void SoPlex<Real>::removeColsReal(int idx[], int n, int perm[])
  {
    if( perm == 0 )
      {
        DataArray< int > p(numColsT());
        _idxToPerm(idx, n, p.get_ptr(), numColsT());
        SoPlex<Real>::removeColsReal(p.get_ptr());
      }
    else
      {
        _idxToPerm(idx, n, perm, numColsT());
        SoPlex<Real>::removeColsReal(perm);
      }
  }



  /// removes columns \p start to \p end including both; an array \p perm of size #numColsT() may be passed as
  /// buffer memory
  template <>
  void SoPlex<Real>::removeColRangeReal(int start, int end, int perm[])
  {
    if( perm == 0 )
      {
        DataArray< int > p(numColsT());
        _rangeToPerm(start, end, p.get_ptr(), numColsT());
        SoPlex<Real>::removeColsReal(p.get_ptr());
      }
    else
      {
        _rangeToPerm(start, end, perm, numColsT());
        SoPlex<Real>::removeColsReal(perm);
      }
  }



  /// clears the LP
  template <>
  void SoPlex<Real>::clearLPReal()
  {
    assert(_realLP != 0);

    _realLP->clear();
    _hasBasis = false;
    _rationalLUSolver.clear();

    if( intParam(SoPlex<Real>::SYNCMODE) == SYNCMODE_AUTO )
      {
        _rationalLP->clear();
        _rowTypes.clear();
        _colTypes.clear();
      }

    _invalidateSolution();
  }



  /// synchronizes real LP with rational LP, i.e., copies (rounded) rational LP into real LP, if sync mode is manual
  template <>
  void SoPlex<Real>::syncLPReal()
  {
    assert(_isConsistent());

    if( intParam(SoPlex<Real>::SYNCMODE) == SYNCMODE_MANUAL )
      _syncLPReal();
  }



  /// adds a single row
  template <>
  void SoPlex<Real>::addRowRational(const LPRowRational& lprow)
  {
    assert(_rationalLP != 0);

    if( intParam(SoPlex<Real>::SYNCMODE) == SYNCMODE_ONLYREAL )
      return;

    _rationalLP->addRow(lprow);
    _completeRangeTypesRational();

    if( intParam(SoPlex<Real>::SYNCMODE) == SYNCMODE_AUTO )
      _addRowReal(lprow);

    _invalidateSolution();
  }



#ifdef SOPLEX_WITH_GMP
  /// adds a single row
  template <>
  void SoPlex<Real>::addRowRational(const mpq_t* lhs, const mpq_t* rowValues, const int* rowIndices, const int rowSize, const mpq_t* rhs)
  {
    assert(_rationalLP != 0);

    if( intParam(SoPlex<Real>::SYNCMODE) == SYNCMODE_ONLYREAL )
      return;

    _rationalLP->addRow(lhs, rowValues, rowIndices, rowSize, rhs);
    _completeRangeTypesRational();

    int i = numRowsT() - 1;
    if( intParam(SoPlex<Real>::SYNCMODE) == SYNCMODE_AUTO )
      _addRowReal(Real(lhsRational(i)), DSVectorReal(_rationalLP->rowVector(i)), Real(rhsRational(i)));

    _invalidateSolution();
  }



  /// adds a set of rows
  template <>
  void SoPlex<Real>::addRowsRational(const mpq_t* lhs, const mpq_t* rowValues, const int* rowIndices, const int* rowStarts, const int* rowLengths, const int numRows, const int numValues, const mpq_t* rhs)
  {
    assert(_rationalLP != 0);

    if( intParam(SoPlex<Real>::SYNCMODE) == SYNCMODE_ONLYREAL )
      return;

    _rationalLP->addRows(lhs, rowValues, rowIndices, rowStarts, rowLengths, numRows, numValues, rhs);
    _completeRangeTypesRational();

    if( intParam(SoPlex<Real>::SYNCMODE) == SYNCMODE_AUTO )
      {
        LPRowSetReal lprowset;
        for( int i = numRowsT() - numRows; i < numRowsT(); i++ )
          lprowset.add(Real(lhsRational(i)), DSVectorReal(_rationalLP->rowVector(i)), Real(rhsRational(i)));
        _addRowsReal(lprowset);
      }

    _invalidateSolution();
  }
#endif



  /// adds multiple rows
  template <>
  void SoPlex<Real>::addRowsRational(const LPRowSetRational& lprowset)
  {
    assert(_rationalLP != 0);

    if( intParam(SoPlex<Real>::SYNCMODE) == SYNCMODE_ONLYREAL )
      return;

    _rationalLP->addRows(lprowset);
    _completeRangeTypesRational();

    if( intParam(SoPlex<Real>::SYNCMODE) == SYNCMODE_AUTO )
      _addRowsReal(lprowset);

    _invalidateSolution();
  }



  /// adds a single column
  template <>
  void SoPlex<Real>::addColRational(const LPColRational& lpcol)
  {
    assert(_rationalLP != 0);

    if( intParam(SoPlex<Real>::SYNCMODE) == SYNCMODE_ONLYREAL )
      return;

    _rationalLP->addCol(lpcol);
    _completeRangeTypesRational();

    if( intParam(SoPlex<Real>::SYNCMODE) == SYNCMODE_AUTO )
      _addColReal(lpcol);

    _invalidateSolution();
  }



#ifdef SOPLEX_WITH_GMP
  /// adds a single column
  template <>
  void SoPlex<Real>::addColRational(const mpq_t* obj, const mpq_t* lower, const mpq_t* colValues, const int* colIndices, const int colSize, const mpq_t* upper)
  {
    assert(_rationalLP != 0);

    if( intParam(SoPlex<Real>::SYNCMODE) == SYNCMODE_ONLYREAL )
      return;

    _rationalLP->addCol(obj, lower, colValues, colIndices, colSize, upper);
    int i = numColsT() - 1;
    _completeRangeTypesRational();

    if( intParam(SoPlex<Real>::SYNCMODE) == SYNCMODE_AUTO )
      _addColReal(Real(maxObjRational(i)) * (intParam(SoPlex<Real>::OBJSENSE) == SoPlex<Real>::OBJSENSE_MAXIMIZE ? 1.0 : -1.0),
                  Real(lowerRational(i)), DSVectorReal(_rationalLP->colVector(i)), Real(upperRational(i)));

    _invalidateSolution();
  }



  /// adds a set of columns
  template <>
  void SoPlex<Real>::addColsRational(const mpq_t* obj, const mpq_t* lower, const mpq_t* colValues, const int* colIndices, const int* colStarts, const int* colLengths, const int numCols, const int numValues, const mpq_t* upper)
  {
    assert(_rationalLP != 0);

    if( intParam(SoPlex<Real>::SYNCMODE) == SYNCMODE_ONLYREAL )
      return;

    _rationalLP->addCols(obj, lower, colValues, colIndices, colStarts, colLengths, numCols, numValues, upper);
    _completeRangeTypesRational();

    if( intParam(SoPlex<Real>::SYNCMODE) == SYNCMODE_AUTO )
      {
        LPColSetReal lpcolset;
        for( int i = numColsT() - numCols; i < numColsT(); i++ )
          lpcolset.add(Real(maxObjRational(i)) * (intParam(SoPlex<Real>::OBJSENSE) == SoPlex<Real>::OBJSENSE_MAXIMIZE ? 1.0 : -1.0),
                       Real(lowerRational(i)), DSVectorReal(_rationalLP->colVector(i)), Real(upperRational(i)));
        _addColsReal(lpcolset);
      }

    _invalidateSolution();
  }
#endif



  /// adds multiple columns
  template <>
  void SoPlex<Real>::addColsRational(const LPColSetRational& lpcolset)
  {
    assert(_rationalLP != 0);

    if( intParam(SoPlex<Real>::SYNCMODE) == SYNCMODE_ONLYREAL )
      return;

    _rationalLP->addCols(lpcolset);
    _completeRangeTypesRational();

    if( intParam(SoPlex<Real>::SYNCMODE) == SYNCMODE_AUTO )
      _addColsReal(lpcolset);

    _invalidateSolution();
  }



  /// replaces row \p i with \p lprow
  template <>
  void SoPlex<Real>::changeRowRational(int i, const LPRowRational& lprow)
  {
    assert(_rationalLP != 0);

    if( intParam(SoPlex<Real>::SYNCMODE) == SYNCMODE_ONLYREAL )
      return;

    _rationalLP->changeRow(i, lprow);
    _rowTypes[i] = _rangeTypeRational(lprow.lhs(), lprow.rhs());
    _completeRangeTypesRational();

    if( intParam(SoPlex<Real>::SYNCMODE) == SYNCMODE_AUTO )
      _changeRowReal(i, lprow);

    _invalidateSolution();
  }



  /// changes left-hand side vector for constraints to \p lhs
  template <>
  void SoPlex<Real>::changeLhsRational(const VectorRational& lhs)
  {
    assert(_rationalLP != 0);

    if( intParam(SoPlex<Real>::SYNCMODE) == SYNCMODE_ONLYREAL )
      return;

    _rationalLP->changeLhs(lhs);
    for( int i = 0; i < numRowsT(); i++ )
      _rowTypes[i] = _rangeTypeRational(lhs[i], _rationalLP->rhs(i));

    if( intParam(SoPlex<Real>::SYNCMODE) == SYNCMODE_AUTO )
      _changeLhsReal(DVectorReal(lhs));

    _invalidateSolution();
  }



  /// changes left-hand side of row \p i to \p lhs
  template <>
  void SoPlex<Real>::changeLhsRational(int i, const Rational& lhs)
  {
    assert(_rationalLP != 0);

    if( intParam(SoPlex<Real>::SYNCMODE) == SYNCMODE_ONLYREAL )
      return;

    _rationalLP->changeLhs(i, lhs);
    _rowTypes[i] = _rangeTypeRational(lhs, _rationalLP->rhs(i));

    if( intParam(SoPlex<Real>::SYNCMODE) == SYNCMODE_AUTO )
      _changeLhsReal(i, Real(lhs));

    _invalidateSolution();
  }



#ifdef SOPLEX_WITH_GMP
  /// changes left-hand side of row \p i to \p lhs
  template <>
  void SoPlex<Real>::changeLhsRational(int i, const mpq_t* lhs)
  {
    assert(_rationalLP != 0);

    if( intParam(SoPlex<Real>::SYNCMODE) == SYNCMODE_ONLYREAL )
      return;

    _rationalLP->changeLhs(i, lhs);
    _rowTypes[i] = _rangeTypeRational(_rationalLP->lhs(i), _rationalLP->rhs(i));

    if( intParam(SoPlex<Real>::SYNCMODE) == SYNCMODE_AUTO )
      _changeLhsReal(i, Real(lhsRational(i)));

    _invalidateSolution();
  }
#endif



  /// changes right-hand side vector to \p rhs
  template <>
  void SoPlex<Real>::changeRhsRational(const VectorRational& rhs)
  {
    assert(_rationalLP != 0);

    if( intParam(SoPlex<Real>::SYNCMODE) == SYNCMODE_ONLYREAL )
      return;

    _rationalLP->changeRhs(rhs);
    for( int i = 0; i < numRowsT(); i++ )
      _rowTypes[i] = _rangeTypeRational(_rationalLP->lhs(i), rhs[i]);

    if( intParam(SoPlex<Real>::SYNCMODE) == SYNCMODE_AUTO )
      _changeRhsReal(DVectorReal(rhs));

    _invalidateSolution();
  }



#ifdef SOPLEX_WITH_GMP
  /// changes right-hand side vector to \p rhs
  template <>
  void SoPlex<Real>::changeRhsRational(const mpq_t* rhs, int rhsSize)
  {
    assert(_rationalLP != 0);

    if( intParam(SoPlex<Real>::SYNCMODE) == SYNCMODE_ONLYREAL )
      return;

    for( int i = 0; i < rhsSize; i++ )
      {
        _rationalLP->changeRhs(i, rhs[i]);
        _rowTypes[i] = _rangeTypeRational(_rationalLP->lhs(i), _rationalLP->rhs(i));
      }

    if( intParam(SoPlex<Real>::SYNCMODE) == SYNCMODE_AUTO )
      _changeRhsReal(DVectorReal(rhsRational()));

    _invalidateSolution();
  }
#endif



  /// changes right-hand side of row \p i to \p rhs
  template <>
  void SoPlex<Real>::changeRhsRational(int i, const Rational& rhs)
  {
    assert(_rationalLP != 0);

    if( intParam(SoPlex<Real>::SYNCMODE) == SYNCMODE_ONLYREAL )
      return;

    _rationalLP->changeRhs(i, rhs);
    _rowTypes[i] = _rangeTypeRational(_rationalLP->lhs(i), rhs);

    if( intParam(SoPlex<Real>::SYNCMODE) == SYNCMODE_AUTO )
      _changeRhsReal(i, Real(rhs));

    _invalidateSolution();
  }



  /// changes left- and right-hand side vectors
  template <>
  void SoPlex<Real>::changeRangeRational(const VectorRational& lhs, const VectorRational& rhs)
  {
    assert(_rationalLP != 0);

    if( intParam(SoPlex<Real>::SYNCMODE) == SYNCMODE_ONLYREAL )
      return;

    _rationalLP->changeRange(lhs, rhs);
    for( int i = 0; i < numRowsT(); i++ )
      _rowTypes[i] = _rangeTypeRational(lhs[i], rhs[i]);

    if( intParam(SoPlex<Real>::SYNCMODE) == SYNCMODE_AUTO )
      _changeRangeReal(DVectorReal(lhs), DVectorReal(rhs));

    _invalidateSolution();
  }



  /// changes left- and right-hand side of row \p i
  template <>
  void SoPlex<Real>::changeRangeRational(int i, const Rational& lhs, const Rational& rhs)
  {
    assert(_rationalLP != 0);

    if( intParam(SoPlex<Real>::SYNCMODE) == SYNCMODE_ONLYREAL )
      return;

    _rationalLP->changeRange(i, lhs, rhs);
    _rowTypes[i] = _rangeTypeRational(lhs, rhs);

    if( intParam(SoPlex<Real>::SYNCMODE) == SYNCMODE_AUTO )
      _changeRangeReal(i, Real(lhs), Real(rhs));

    _invalidateSolution();
  }



#ifdef SOPLEX_WITH_GMP
  /// changes left-hand side of row \p i to \p lhs
  template <>
  void SoPlex<Real>::changeRangeRational(int i, const mpq_t* lhs, const mpq_t* rhs)
  {
    assert(_rationalLP != 0);

    if( intParam(SoPlex<Real>::SYNCMODE) == SYNCMODE_ONLYREAL )
      return;

    _rationalLP->changeRange(i, lhs, rhs);
    _rowTypes[i] = _rangeTypeRational(_rationalLP->lhs(i), _rationalLP->rhs(i));

    if( intParam(SoPlex<Real>::SYNCMODE) == SYNCMODE_AUTO )
      _changeRangeReal(i, Real(lhsRational(i)), Real(rhsRational(i)));

    _invalidateSolution();
  }
#endif



  /// replaces column \p i with \p lpcol
  template <>
  void SoPlex<Real>::changeColRational(int i, const LPColRational& lpcol)
  {
    assert(_rationalLP != 0);

    if( intParam(SoPlex<Real>::SYNCMODE) == SYNCMODE_ONLYREAL )
      return;

    _rationalLP->changeCol(i, lpcol);
    _colTypes[i] = _rangeTypeRational(lpcol.lower(), lpcol.upper());
    _completeRangeTypesRational();

    if( intParam(SoPlex<Real>::SYNCMODE) == SYNCMODE_AUTO )
      _changeColReal(i, lpcol);

    _invalidateSolution();
  }



  /// changes vector of lower bounds to \p lower
  template <>
  void SoPlex<Real>::changeLowerRational(const VectorRational& lower)
  {
    assert(_rationalLP != 0);

    if( intParam(SoPlex<Real>::SYNCMODE) == SYNCMODE_ONLYREAL )
      return;

    _rationalLP->changeLower(lower);
    for( int i = 0; i < numColsT(); i++ )
      _colTypes[i] = _rangeTypeRational(lower[i], _rationalLP->upper(i));

    if( intParam(SoPlex<Real>::SYNCMODE) == SYNCMODE_AUTO )
      _changeLowerReal(DVectorReal(lower));

    _invalidateSolution();
  }



  /// changes lower bound of column i to \p lower
  template <>
  void SoPlex<Real>::changeLowerRational(int i, const Rational& lower)
  {
    assert(_rationalLP != 0);

    if( intParam(SoPlex<Real>::SYNCMODE) == SYNCMODE_ONLYREAL )
      return;

    _rationalLP->changeLower(i, lower);
    _colTypes[i] = _rangeTypeRational(lower, _rationalLP->upper(i));

    if( intParam(SoPlex<Real>::SYNCMODE) == SYNCMODE_AUTO )
      _changeLowerReal(i, Real(lower));

    _invalidateSolution();
  }



#ifdef SOPLEX_WITH_GMP
  /// changes lower bound of column i to \p lower
  template <>
  void SoPlex<Real>::changeLowerRational(int i, const mpq_t* lower)
  {
    assert(_rationalLP != 0);

    if( intParam(SoPlex<Real>::SYNCMODE) == SYNCMODE_ONLYREAL )
      return;

    _rationalLP->changeLower(i, lower);
    _colTypes[i] = _rangeTypeRational(_rationalLP->lower(i), _rationalLP->upper(i));

    if( intParam(SoPlex<Real>::SYNCMODE) == SYNCMODE_AUTO )
      _changeLowerReal(i, Real(lowerRational(i)));

    _invalidateSolution();
  }
#endif



  /// changes vector of upper bounds to \p upper
  template <>
  void SoPlex<Real>::changeUpperRational(const VectorRational& upper)
  {
    assert(_rationalLP != 0);

    if( intParam(SoPlex<Real>::SYNCMODE) == SYNCMODE_ONLYREAL )
      return;

    _rationalLP->changeUpper(upper);
    for( int i = 0; i < numColsT(); i++ )
      _colTypes[i] = _rangeTypeRational(_rationalLP->lower(i), upper[i]);

    if( intParam(SoPlex<Real>::SYNCMODE) == SYNCMODE_AUTO )
      _changeUpperReal(DVectorReal(upper));

    _invalidateSolution();
  }



  /// changes \p i 'th upper bound to \p upper
  template <>
  void SoPlex<Real>::changeUpperRational(int i, const Rational& upper)
  {
    assert(_rationalLP != 0);

    if( intParam(SoPlex<Real>::SYNCMODE) == SYNCMODE_ONLYREAL )
      return;

    _rationalLP->changeUpper(i, upper);
    _colTypes[i] = _rangeTypeRational(_rationalLP->lower(i), upper);

    if( intParam(SoPlex<Real>::SYNCMODE) == SYNCMODE_AUTO )
      _changeUpperReal(i, Real(upper));

    _invalidateSolution();
  }



#ifdef SOPLEX_WITH_GMP
  /// changes upper bound of column i to \p upper
  template <>
  void SoPlex<Real>::changeUpperRational(int i, const mpq_t* upper)
  {
    assert(_rationalLP != 0);

    if( intParam(SoPlex<Real>::SYNCMODE) == SYNCMODE_ONLYREAL )
      return;

    _rationalLP->changeUpper(i, upper);
    _colTypes[i] = _rangeTypeRational(_rationalLP->lower(i), _rationalLP->upper(i));

    if( intParam(SoPlex<Real>::SYNCMODE) == SYNCMODE_AUTO )
      _changeUpperReal(i, Real(upperRational(i)));

    _invalidateSolution();
  }
#endif



  /// changes vectors of column bounds to \p lower and \p upper
  template <>
  void SoPlex<Real>::changeBoundsRational(const VectorRational& lower, const VectorRational& upper)
  {
    assert(_rationalLP != 0);

    if( intParam(SoPlex<Real>::SYNCMODE) == SYNCMODE_ONLYREAL )
      return;

    _rationalLP->changeBounds(lower, upper);
    for( int i = 0; i < numColsT(); i++ )
      _colTypes[i] = _rangeTypeRational(lower[i], upper[i]);

    if( intParam(SoPlex<Real>::SYNCMODE) == SYNCMODE_AUTO )
      _changeBoundsReal(DVectorReal(lower), DVectorReal(upper));

    _invalidateSolution();
  }



  /// changes bounds of column \p i to \p lower and \p upper
  template <>
  void SoPlex<Real>::changeBoundsRational(int i, const Rational& lower, const Rational& upper)
  {
    assert(_rationalLP != 0);

    if( intParam(SoPlex<Real>::SYNCMODE) == SYNCMODE_ONLYREAL )
      return;

    _rationalLP->changeBounds(i, lower, upper);
    _colTypes[i] = _rangeTypeRational(lower, upper);

    if( intParam(SoPlex<Real>::SYNCMODE) == SYNCMODE_AUTO )
      _changeBoundsReal(i, Real(lower), Real(upper));

    _invalidateSolution();
  }



#ifdef SOPLEX_WITH_GMP
  /// changes bounds of column \p i to \p lower and \p upper
  template <>
  void SoPlex<Real>::changeBoundsRational(int i, const mpq_t* lower, const mpq_t* upper)
  {
    assert(_rationalLP != 0);

    if( intParam(SoPlex<Real>::SYNCMODE) == SYNCMODE_ONLYREAL )
      return;

    _rationalLP->changeBounds(i, lower, upper);
    _colTypes[i] = _rangeTypeRational(_rationalLP->lower(i), _rationalLP->upper(i));

    if( intParam(SoPlex<Real>::SYNCMODE) == SYNCMODE_AUTO )
      _changeBoundsReal(i, Real(lowerRational(i)), Real(upperRational(i)));

    _invalidateSolution();
  }
#endif



  /// changes objective function vector to \p obj
  template <>
  void SoPlex<Real>::changeObjRational(const VectorRational& obj)
  {
    assert(_rationalLP != 0);

    if( intParam(SoPlex<Real>::SYNCMODE) == SYNCMODE_ONLYREAL )
      return;

    _rationalLP->changeObj(obj);

    if( intParam(SoPlex<Real>::SYNCMODE) == SYNCMODE_AUTO )
      _realLP->changeObj(DVectorReal(obj));

    _invalidateSolution();
  }



  /// changes objective coefficient of column i to \p obj
  template <>
  void SoPlex<Real>::changeObjRational(int i, const Rational& obj)
  {
    assert(_rationalLP != 0);

    if( intParam(SoPlex<Real>::SYNCMODE) == SYNCMODE_ONLYREAL )
      return;

    _rationalLP->changeObj(i, obj);

    if( intParam(SoPlex<Real>::SYNCMODE) == SYNCMODE_AUTO )
      _realLP->changeObj(i, Real(obj));

    _invalidateSolution();
  }



#ifdef SOPLEX_WITH_GMP
  /// changes objective coefficient of column i to \p obj
  template <>
  void SoPlex<Real>::changeObjRational(int i, const mpq_t* obj)
  {
    assert(_rationalLP != 0);

    if( intParam(SoPlex<Real>::SYNCMODE) == SYNCMODE_ONLYREAL )
      return;

    _rationalLP->changeObj(i, obj);

    if( intParam(SoPlex<Real>::SYNCMODE) == SYNCMODE_AUTO )
      _realLP->changeObj(i, Real(objRational(i)));

    _invalidateSolution();
  }
#endif



  /// changes matrix entry in row \p i and column \p j to \p val
  template <>
  void SoPlex<Real>::changeElementRational(int i, int j, const Rational& val)
  {
    assert(_rationalLP != 0);

    if( intParam(SoPlex<Real>::SYNCMODE) == SYNCMODE_ONLYREAL )
      return;

    _rationalLP->changeElement(i, j, val);

    if( intParam(SoPlex<Real>::SYNCMODE) == SYNCMODE_AUTO )
      _changeElementReal(i, j, Real(val));

    _invalidateSolution();
  }


#ifdef SOPLEX_WITH_GMP
  /// changes matrix entry in row \p i and column \p j to \p val
  template <>
  void SoPlex<Real>::changeElementRational(int i, int j, const mpq_t* val)
  {
    assert(_rationalLP != 0);

    if( intParam(SoPlex<Real>::SYNCMODE) == SYNCMODE_ONLYREAL )
      return;

    _rationalLP->changeElement(i, j, val);

    if( intParam(SoPlex<Real>::SYNCMODE) == SYNCMODE_AUTO )
      _changeElementReal(i, j, mpq_get_d(*val));

    _invalidateSolution();
  }
#endif


  /// removes row \p i
  template <>
  void SoPlex<Real>::removeRowRational(int i)
  {
    assert(_rationalLP != 0);

    if( intParam(SoPlex<Real>::SYNCMODE) == SYNCMODE_ONLYREAL )
      return;

    _rationalLP->removeRow(i);
    // only swap elements if not the last one was removed
    if( i < _rationalLP->nRows() )
      {
        _rowTypes[i] = _rowTypes[_rationalLP->nRows()];
        assert(_rowTypes[i] == _rangeTypeRational(lhsRational(i), rhsRational(i)));
      }
    _rowTypes.reSize(_rationalLP->nRows());

    if( intParam(SoPlex<Real>::SYNCMODE) == SYNCMODE_AUTO )
      _removeRowReal(i);

    _invalidateSolution();
  }



  /// removes all rows with an index \p i such that \p perm[i] < 0; upon completion, \p perm[i] >= 0 indicates the new
  /// index where row \p i has been moved to; note that \p perm must point to an array of size at least
  /// #numRowsT()
  template <>
  void SoPlex<Real>::removeRowsRational(int perm[])
  {
    assert(_rationalLP != 0);

    if( intParam(SoPlex<Real>::SYNCMODE) == SYNCMODE_ONLYREAL )
      return;

    const int oldsize = numRowsT();
    _rationalLP->removeRows(perm);
    for( int i = 0; i < oldsize; i++ )
      {
        if( perm[i] >= 0 )
          _rowTypes[perm[i]] = _rowTypes[i];
      }
    _rowTypes.reSize(_rationalLP->nRows());
    for( int i = 0; i < numRowsT(); i++ )
      {
        assert(_rowTypes[i] == _rangeTypeRational(lhsRational(i), rhsRational(i)));
      }


    if( intParam(SoPlex<Real>::SYNCMODE) == SYNCMODE_AUTO )
      _removeRowsReal(perm);

    _invalidateSolution();
  }



  /// remove all rows with indices in array \p idx of size \p n; an array \p perm of size #numRowsT() may be
  /// passed as buffer memory
  template <>
  void SoPlex<Real>::removeRowsRational(int idx[], int n, int perm[])
  {
    if( perm == 0 )
      {
        DataArray< int > p(numRowsT());
        _idxToPerm(idx, n, p.get_ptr(), numRowsT());
        SoPlex<Real>::removeRowsRational(p.get_ptr());
      }
    else
      {
        _idxToPerm(idx, n, perm, numRowsT());
        SoPlex<Real>::removeRowsRational(perm);
      }
  }



  /// removes rows \p start to \p end including both; an array \p perm of size #numRowsT() may be passed as
  /// buffer memory
  template <>
  void SoPlex<Real>::removeRowRangeRational(int start, int end, int perm[])
  {
    if( perm == 0 )
      {
        DataArray< int > p(numRowsT());
        _rangeToPerm(start, end, p.get_ptr(), numRowsT());
        SoPlex<Real>::removeRowsRational(p.get_ptr());
      }
    else
      {
        _rangeToPerm(start, end, perm, numRowsT());
        SoPlex<Real>::removeRowsRational(perm);
      }
  }



  /// removes column i
  template <>
  void SoPlex<Real>::removeColRational(int i)
  {
    assert(_rationalLP != 0);

    if( intParam(SoPlex<Real>::SYNCMODE) == SYNCMODE_ONLYREAL )
      return;

    _rationalLP->removeCol(i);
    // only swap elements if not the last one was removed
    if( i < _rationalLP->nCols() )
      {
        _colTypes[i] = _colTypes[_rationalLP->nCols()];
        assert(_colTypes[i] == _rangeTypeRational(lowerRational(i), upperRational(i)));
      }
    _colTypes.reSize(_rationalLP->nCols());

    if( intParam(SoPlex<Real>::SYNCMODE) == SYNCMODE_AUTO )
      _removeColReal(i);

    _invalidateSolution();
  }



  /// removes all columns with an index \p i such that \p perm[i] < 0; upon completion, \p perm[i] >= 0 indicates the
  /// new index where column \p i has been moved to; note that \p perm must point to an array of size at least
  /// #numColsT()
  template <>
  void SoPlex<Real>::removeColsRational(int perm[])
  {
    assert(_rationalLP != 0);

    if( intParam(SoPlex<Real>::SYNCMODE) == SYNCMODE_ONLYREAL )
      return;

    const int oldsize = numColsT();
    _rationalLP->removeCols(perm);
    for( int i = 0; i < oldsize; i++ )
      {
        if( perm[i] >= 0 )
          _colTypes[perm[i]] = _colTypes[i];
      }
    _colTypes.reSize(_rationalLP->nCols());
    for( int i = 0; i < numColsT(); i++ )
      {
        assert(_colTypes[i] == _rangeTypeRational(lowerRational(i), upperRational(i)));
      }

    if( intParam(SoPlex<Real>::SYNCMODE) == SYNCMODE_AUTO )
      _removeColsReal(perm);

    _invalidateSolution();
  }



  /// remove all columns with indices in array \p idx of size \p n; an array \p perm of size #numColsT() may be
  /// passed as buffer memory
  template <>
  void SoPlex<Real>::removeColsRational(int idx[], int n, int perm[])
  {
    if( perm == 0 )
      {
        DataArray< int > p(numColsT());
        _idxToPerm(idx, n, p.get_ptr(), numColsT());
        SoPlex<Real>::removeColsRational(p.get_ptr());
      }
    else
      {
        _idxToPerm(idx, n, perm, numColsT());
        SoPlex<Real>::removeColsRational(perm);
      }
  }



  /// removes columns \p start to \p end including both; an array \p perm of size #numColsT() may be passed as
  /// buffer memory
  template <>
  void SoPlex<Real>::removeColRangeRational(int start, int end, int perm[])
  {
    if( perm == 0 )
      {
        DataArray< int > p(numColsT());
        _rangeToPerm(start, end, p.get_ptr(), numColsT());
        SoPlex<Real>::removeColsRational(p.get_ptr());
      }
    else
      {
        _rangeToPerm(start, end, perm, numColsT());
        SoPlex<Real>::removeColsRational(perm);
      }
  }



  /// clears the LP
  template <>
  void SoPlex<Real>::clearLPRational()
  {
    assert(_rationalLP != 0);

    if( intParam(SoPlex<Real>::SYNCMODE) == SYNCMODE_ONLYREAL )
      return;

    _rationalLP->clear();
    _rationalLUSolver.clear();
    _rowTypes.clear();
    _colTypes.clear();

    if( intParam(SoPlex<Real>::SYNCMODE) == SYNCMODE_AUTO )
      {
        _realLP->clear();
        _hasBasis = false;
      }

    _invalidateSolution();
  }



  /// synchronizes rational LP with real LP, i.e., copies real LP to rational LP, if sync mode is manual
  template <>
  void SoPlex<Real>::syncLPRational()
  {
    assert(_isConsistent());

    if( intParam(SoPlex<Real>::SYNCMODE) == SYNCMODE_MANUAL )
      _syncLPRational();
  }

  /// returns real parameter value
  /// #template #temp
  template <>
  Real SoPlex<Real>::realParam(const RealParam param) const
  {
    assert(param >= 0);
    assert(param < REALPARAM_COUNT);
    return _currentSettings->_realParamValues[param];
  }

  template <>
  Real SoPlex<Rational>::realParam(const RealParam param) const
  {
    std::exit(EXIT_FAILURE);
  }

  /// solves the LP
  /// #template needs to be rewritten #temp
  template <>
  typename SPxSolver<Real>::Status SoPlex<Real>::optimize()
  {
    assert(_isConsistent());

    // clear statistics
    _statistics->clearSolvingData();

    // the solution is no longer valid
    _invalidateSolution();

    // if the decomposition based dual simplex flag is set to true
    if ( boolParam(SoPlex<Real>::USEDECOMPDUALSIMPLEX) )
      {
        setIntParam(SoPlex<Real>::SOLVEMODE, SOLVEMODE_REAL);
        setIntParam(SoPlex<Real>::REPRESENTATION, REPRESENTATION_ROW);
        setIntParam(SoPlex<Real>::ALGORITHM, ALGORITHM_DUAL);
        //setBoolParam(SoPlex<Real>::PERSISTENTSCALING, false);

        _solver.setComputeDegenFlag(boolParam(COMPUTEDEGEN));

        _solveDecompositionDualSimplex();
      }
    // decide whether to solve the rational LP with iterative refinement or call the standard floating-point solver
    else if( intParam(SoPlex<Real>::SOLVEMODE) == SOLVEMODE_REAL || (intParam(SoPlex<Real>::SOLVEMODE) == SOLVEMODE_AUTO
                                                                  && GE(realParam(SoPlex<Real>::FEASTOL), 1e-9) && GE(realParam(SoPlex<Real>::OPTTOL), 1e-9)) )
      {
        // ensure that tolerances are reasonable for the floating-point solver
        if( realParam(SoPlex<Real>::FEASTOL) < _currentSettings->realParam.lower[SoPlex<Real>::FPFEASTOL] )
          {
            MSG_WARNING( spxout, spxout << "Cannot call floating-point solver with feasibility tolerance below "
                         << _currentSettings->realParam.lower[SoPlex<Real>::FPFEASTOL] << " - relaxing tolerance\n");
            _solver.setFeastol(_currentSettings->realParam.lower[SoPlex<Real>::FPFEASTOL]);
          }
        else
          _solver.setFeastol(realParam(SoPlex<Real>::FEASTOL));

        if( realParam(SoPlex<Real>::OPTTOL) < _currentSettings->realParam.lower[SoPlex<Real>::FPOPTTOL] )
          {
            MSG_WARNING( spxout, spxout << "Cannot call floating-point solver with optimality tolerance below "
                         << _currentSettings->realParam.lower[SoPlex<Real>::FPOPTTOL] << " - relaxing tolerance\n");
            _solver.setOpttol(_currentSettings->realParam.lower[SoPlex<Real>::FPOPTTOL]);
          }
        else
          _solver.setOpttol(realParam(SoPlex<Real>::OPTTOL));

        _solver.setComputeDegenFlag(boolParam(COMPUTEDEGEN));

        _optimizeT();
#ifdef SOPLEX_DEBUG // this check will remove scaling of the realLP
        _checkBasisScaling();
#endif
      }
    else if( intParam(SoPlex<Real>::SYNCMODE) == SYNCMODE_ONLYREAL )
      {
        _syncLPRational();
        _optimizeRational();
      }
    else if( intParam(SoPlex<Real>::SYNCMODE) == SYNCMODE_MANUAL )
      {
#ifdef ENABLE_ADDITIONAL_CHECKS
        assert(areLPsInSync(true, true, false));
#else
        assert(areLPsInSync(true, false, false));
#endif

        _optimizeRational();

#ifdef ENABLE_ADDITIONAL_CHECKS
        assert(areLPsInSync(true, true, false));
#else
        assert(areLPsInSync(true, false, false));
#endif
      }
    else
      {
#ifdef ENABLE_ADDITIONAL_CHECKS
        assert(areLPsInSync(true, true, false));
#else
        assert(areLPsInSync(true, false, false));
#endif

        _optimizeRational();
      }

    MSG_INFO1( spxout, spxout << "\n";
               printShortStatistics(spxout.getStream(SPxOut::INFO1));
               spxout << "\n" );


    return status();
  }

  /// is stored primal solution feasible?
  template <>
	bool SoPlex<Real>::isPrimalFeasible() const
  {
    return (_hasSolReal && _solReal.isPrimalFeasible());
  }

  template <>
	bool SoPlex<Rational>::isPrimalFeasible() const
  {
    return (_hasSolRational && _solRational.isPrimalFeasible());
  }



  /// is a primal feasible solution available?
  template <>
	bool SoPlex<Real>::hasPrimal() const
  {
    return _hasSolReal;
  }

  template <>
	bool SoPlex<Rational>::hasPrimal() const
  {
    return _hasSolRational;
  }

  /// is a primal unbounded ray available?
  template <>
	bool SoPlex<Real>::hasPrimalRay() const
  {
    return (_hasSolReal && _solReal.hasPrimalRay()) || (_hasSolRational && _solRational.hasPrimalRay());
  }



  /// Old code for isDualFeasible() #template remove later
  // template <>
	// bool SoPlex<Real>::isDualFeasible() const
  // {
  //   return (_hasSolReal && _solReal.isDualFeasible()) || (_hasSolRational && _solRational.isDualFeasible());

  // }

  /// is stored dual solution feasible?
  template<>
  bool SoPlex<Real>::isDualFeasible() const
  {
    return (_hasSolReal && _solReal.isDualFeasible());
  }

  template <>
  bool SoPlex<Rational>::isDualFeasible() const
  {
    return (_hasSolRational && _solRational.isDualFeasible());
  }

  /// is a dual feasible solution available?
  template <>
	bool SoPlex<Real>::hasDual() const
  {
    return _hasSolReal;
  }

  template <>
	bool SoPlex<Rational>::hasDual() const
  {
    return _hasSolRational;
  }

  /// is Farkas proof of infeasibility available?
  template <>
	bool SoPlex<Real>::hasDualFarkas() const
  {
    return (_hasSolReal && _solReal.hasDualFarkas()) || (_hasSolRational && _solRational.hasDualFarkas());
  }



  /// returns the objective value if a primal or dual solution is available
  template <>
  Real SoPlex<Real>::objValueReal()
  {
    assert(OBJSENSE_MAXIMIZE == 1);
    assert(OBJSENSE_MINIMIZE == -1);

    if( status() == SPxSolver<Real>::UNBOUNDED )
      return realParam(SoPlex<Real>::INFTY) * intParam(SoPlex<Real>::OBJSENSE);
    else if( status() == SPxSolver<Real>::INFEASIBLE )
      return -realParam(SoPlex<Real>::INFTY) * intParam(SoPlex<Real>::OBJSENSE);
    else if( hasPrimal() || hasDual() )
      {
        _syncRealSolution();
        return _solReal._objVal;
      }
    else
      return 0.0;
  }



  /// gets the primal solution vector if available; returns true on success
  template <>
	bool SoPlex<Real>::getPrimalT(VectorBase<Real>& vector)
  {
    if( hasPrimal() && vector.dim() >= numColsT() )
      {
        _syncRealSolution();
        _solReal.getPrimal(vector);
        return true;
      }
    else
      return false;
  }



  /// gets the vector of slack values if available; returns true on success
  template <>
	bool SoPlex<Real>::getSlacksReal(VectorReal& vector)
  {
    if( hasPrimal() && vector.dim() >= numRowsT() )
      {
        _syncRealSolution();
        _solReal.getSlacks(vector);
        return true;
      }
    else
      return false;
  }



  /// gets the primal ray if available; returns true on success
  template <>
	bool SoPlex<Real>::getPrimalRayT(VectorBase<Real>& vector)
  {
    if( hasPrimalRay() && vector.dim() >= numColsT() )
      {
        _syncRealSolution();
        _solReal.getPrimalRay(vector);
        return true;
      }
    else
      return false;
  }



  /// gets the dual solution vector if available; returns true on success
  template <>
	bool SoPlex<Real>::getDualT(VectorBase<Real>& vector)
  {
    if( hasDual() && vector.dim() >= numRowsT() )
      {
        _syncRealSolution();
        _solReal.getDual(vector);
        return true;
      }
    else
      return false;
  }



  /// gets the vector of reduced cost values if available; returns true on success
  template <>
	bool SoPlex<Real>::getRedCostT(VectorBase<Real>& vector)
  {
    if( hasDual() && vector.dim() >= numColsT() )
      {
        _syncRealSolution();
        _solReal.getRedCost(vector);
        return true;
      }
    else
      return false;
  }



  /// gets the Farkas proof if available; returns true on success
  template <>
	bool SoPlex<Real>::getDualFarkasT(VectorBase<Real>& vector)
  {
    if( hasDualFarkas() && vector.dim() >= numRowsT() )
      {
        _syncRealSolution();
        _solReal.getDualFarkas(vector);
        return true;
      }
    else
      return false;
  }



  /// gets violation of bounds; returns true on success
  template <>
	bool SoPlex<Real>::getBoundViolationT(Real& maxviol, Real& sumviol)
  {
    if( !isPrimalFeasible() )
      return false;

    _syncRealSolution();
    VectorReal& primal = _solReal._primal;
    assert(primal.dim() == numColsT());

    maxviol = 0.0;
    sumviol = 0.0;

    for( int i = numColsT() - 1; i >= 0; i-- )
      {
        Real lower = _realLP->lowerUnscaled(i);
        Real upper = _realLP->upperUnscaled(i);
        Real viol = lower - primal[i];
        if( viol > 0.0 )
          {
            sumviol += viol;
            if( viol > maxviol )
              maxviol = viol;
          }

        viol = primal[i] - upper;
        if( viol > 0.0 )
          {
            sumviol += viol;
            if( viol > maxviol )
              maxviol = viol;
          }
      }

    return true;
  }



  /// gets violation of constraints; returns true on success
  template <>
	bool SoPlex<Real>::getRowViolationT(Real& maxviol, Real& sumviol)
  {
    if( !isPrimalFeasible() )
      return false;

    _syncRealSolution();
    VectorReal& primal = _solReal._primal;
    assert(primal.dim() == numColsT());

    DVectorReal activity(numRowsT());
    _realLP->computePrimalActivity(primal, activity, true);
    maxviol = 0.0;
    sumviol = 0.0;

    for( int i = numRowsT() - 1; i >= 0; i-- )
      {
        Real lhs = _realLP->lhsUnscaled(i);
        Real rhs = _realLP->rhsUnscaled(i);

        Real viol = lhs - activity[i];
        if( viol > 0.0 )
          {
            sumviol += viol;
            if( viol > maxviol )
              maxviol = viol;
          }

        viol = activity[i] - rhs;
        if( viol > 0.0 )
          {
            sumviol += viol;
            if( viol > maxviol )
              maxviol = viol;
          }
      }

    return true;
  }

  /// is an advanced starting basis available?
  /// #template #baseclass
  template <>
	bool SoPlex<Real>::hasBasis() const
  {
    return _hasBasis;
  }

  template <>
	bool SoPlex<Rational>::hasBasis() const
  {
    return _hasBasis;
  }

  /// gets violation of reduced costs; returns true on success
  template <>
	bool SoPlex<Real>::getRedCostViolationT(Real& maxviol, Real& sumviol)
  {
    if( !isDualFeasible() || !hasBasis() )
      return false;

    _syncRealSolution();
    VectorReal& redcost = _solReal._redCost;
    assert(redcost.dim() == numColsT());

    maxviol = 0.0;
    sumviol = 0.0;

    for( int c = numColsT() - 1; c >= 0; c-- )
      {
        typename SPxSolver<Real>::VarStatus colStatus = basisColStatus(c);

        if( intParam(SoPlex<Real>::OBJSENSE) == OBJSENSE_MINIMIZE )
          {
            if( colStatus != SPxSolver<Real>::ON_UPPER && colStatus != SPxSolver<Real>::FIXED && redcost[c] < 0.0 )
              {
                sumviol += -redcost[c];
                if( redcost[c] < -maxviol )
                  maxviol = -redcost[c];
              }
            if( colStatus != SPxSolver<Real>::ON_LOWER && colStatus != SPxSolver<Real>::FIXED && redcost[c] > 0.0 )
              {
                sumviol += redcost[c];
                if( redcost[c] > maxviol )
                  maxviol = redcost[c];
              }
          }
        else
          {
            if( colStatus != SPxSolver<Real>::ON_UPPER && colStatus != SPxSolver<Real>::FIXED && redcost[c] > 0.0 )
              {
                sumviol += redcost[c];
                if( redcost[c] > maxviol )
                  maxviol = redcost[c];
              }
            if( colStatus != SPxSolver<Real>::ON_LOWER && colStatus != SPxSolver<Real>::FIXED && redcost[c] < 0.0 )
              {
                sumviol += -redcost[c];
                if( redcost[c] < -maxviol )
                  maxviol = -redcost[c];
              }
          }
      }

    return true;
  }



  /// gets violation of dual multipliers; returns true on success
  template <>
	bool SoPlex<Real>::getDualViolationT(Real& maxviol, Real& sumviol)
  {
    if( !isDualFeasible() || !hasBasis() )
      return false;

    _syncRealSolution();
    VectorReal& dual = _solReal._dual;
    assert(dual.dim() == numRowsT());

    maxviol = 0.0;
    sumviol = 0.0;

    for( int r = numRowsT() - 1; r >= 0; r-- )
      {
        typename SPxSolver<Real>::VarStatus rowStatus = basisRowStatus(r);

        if( intParam(SoPlex<Real>::OBJSENSE) == OBJSENSE_MINIMIZE )
          {
            if( rowStatus != SPxSolver<Real>::ON_UPPER && rowStatus != SPxSolver<Real>::FIXED && dual[r] < 0.0 )
              {
                sumviol += -dual[r];
                if( dual[r] < -maxviol )
                  maxviol = -dual[r];
              }
            if( rowStatus != SPxSolver<Real>::ON_LOWER && rowStatus != SPxSolver<Real>::FIXED && dual[r] > 0.0 )
              {
                sumviol += dual[r];
                if( dual[r] > maxviol )
                  maxviol = dual[r];
              }
          }
        else
          {
            if( rowStatus != SPxSolver<Real>::ON_UPPER && rowStatus != SPxSolver<Real>::FIXED && dual[r] > 0.0 )
              {
                sumviol += dual[r];
                if( dual[r] > maxviol )
                  maxviol = dual[r];
              }
            if( rowStatus != SPxSolver<Real>::ON_LOWER && rowStatus != SPxSolver<Real>::FIXED && dual[r] < 0.0 )
              {
                sumviol += -dual[r];
                if( dual[r] < -maxviol )
                  maxviol = -dual[r];
              }
          }
      }

    return true;
  }



  /// returns the objective value if a primal or dual solution is available
  template <>
  Rational SoPlex<Real>::objValueRational()
  {
    assert(OBJSENSE_MAXIMIZE == 1);
    assert(OBJSENSE_MINIMIZE == -1);

    if( status() == SPxSolver<Real>::UNBOUNDED )
      {
        if( intParam(SoPlex<Real>::OBJSENSE) == OBJSENSE_MAXIMIZE )
          return _rationalPosInfty;
        else
          return _rationalNegInfty;
      }
    else if( status() == SPxSolver<Real>::INFEASIBLE )
      {
        if( intParam(SoPlex<Real>::OBJSENSE) == OBJSENSE_MAXIMIZE )
          return _rationalNegInfty;
        else
          return _rationalPosInfty;
      }
    else if( hasPrimal() || hasDual() )
      {
        _syncRationalSolution();
        return _solRational._objVal;
      }
    else
      return _rationalZero;
  }


  /// gets the primal solution vector if available; returns true on success
  template <>
	bool SoPlex<Real>::getPrimalRational(VectorBase<Rational>& vector)
  {
    if( _rationalLP != 0 && hasPrimal() && vector.dim() >= numColsT() )
      {
        _syncRationalSolution();
        _solRational.getPrimal(vector);
        return true;
      }
    else
      return false;
  }

  /// gets the vector of slack values if available; returns true on success
  template <>
	bool SoPlex<Real>::getSlacksRational(VectorRational& vector)
  {
    if( _rationalLP != 0 && hasPrimal() && vector.dim() >= numRowsT() )
      {
        _syncRationalSolution();
        _solRational.getSlacks(vector);
        return true;
      }
    else
      return false;
  }



  /// gets the primal ray if LP is unbounded; returns true on success
  // template <>
	// bool SoPlex<Rational>::getPrimalRayT(VectorBase<Rational>& vector)
  // {
  //   if( _rationalLP != 0 && hasPrimalRay() && vector.dim() >= numColsT() )
  //     {
  //       _syncRationalSolution();
  //       _solRational.getPrimalRay(vector);
  //       return true;
  //     }
  //   else
  //     return false;
  // }

  template <>
	bool SoPlex<Real>::getPrimalRayRational(VectorBase<Rational>& vector)
  {
    if( _rationalLP != 0 && hasPrimalRay() && vector.dim() >= numColsT() )
      {
        _syncRationalSolution();
        _solRational.getPrimalRay(vector);
        return true;
      }
    else
      return false;
  }



  /// gets the dual solution vector if available; returns true on success
  template <>
	bool SoPlex<Real>::getDualRational(VectorBase<Rational>& vector)
  {
    if( _rationalLP != 0 && hasDual() && vector.dim() >= numRowsT() )
      {
        _syncRationalSolution();
        _solRational.getDual(vector);
        return true;
      }
    else
      return false;
  }



  /// gets the vector of reduced cost values if available; returns true on success
  template <>
	bool SoPlex<Real>::getRedCostRational(VectorRational& vector)
  {
    if( _rationalLP != 0 && hasDual() && vector.dim() >= numColsT() )
      {
        _syncRationalSolution();
        _solRational.getRedCost(vector);
        return true;
      }
    else
      return false;
  }

  /// gets the Farkas proof if LP is infeasible; returns true on success
  template <>
	bool SoPlex<Real>::getDualFarkasRational(VectorBase<Rational>& vector)
  {
    if( _rationalLP != 0 && hasDualFarkas() && vector.dim() >= numRowsT() )
      {
        _syncRationalSolution();
        _solRational.getDualFarkas(vector);
        return true;
      }
    else
      return false;
  }



  /// gets violation of bounds; returns true on success
  template <>
	bool SoPlex<Real>::getBoundViolationRational(Rational& maxviol, Rational& sumviol)
  {
    if( !isPrimalFeasible() )
      return false;

    // if we have to synchronize, we do not measure time, because this would affect the solving statistics
    if( intParam(SoPlex<Real>::SYNCMODE) == SYNCMODE_ONLYREAL )
      _syncLPRational(false);

    _syncRationalSolution();
    VectorRational& primal = _solRational._primal;
    assert(primal.dim() == numColsRational());

    maxviol = 0;
    sumviol = 0;

    for( int i = numColsT() - 1; i >= 0; i-- )
      {
        Rational viol = lowerRational(i) - primal[i];
        if( viol > 0 )
          {
            sumviol += viol;
            if( viol > maxviol )
              {
                maxviol = viol;
                MSG_DEBUG( std::cout << "increased bound violation for column " << i << ": " << rationalToString(viol)
                           << " lower: " << rationalToString(lowerRational(i))
                           << ", primal: " << rationalToString(primal[i]) << "\n" );
              }
          }

        viol = primal[i] - upperRational(i);
        if( viol > 0 )
          {
            sumviol += viol;
            if( viol > maxviol )
              {
                maxviol = viol;
                MSG_DEBUG( std::cout << "increased bound violation for column " << i << ": " << rationalToString(viol)
                           << " upper: " << rationalToString(upperRational(i))
                           << ", primal: " << rationalToString(primal[i]) << "\n" );
              }
          }
      }

    return true;
  }



  /// gets violation of constraints; returns true on success
  template <>
	bool SoPlex<Real>::getRowViolationRational(Rational& maxviol, Rational& sumviol)
  {
    if( !isPrimalFeasible() )
      return false;

    // if we have to synchronize, we do not measure time, because this would affect the solving statistics
    if( intParam(SoPlex<Real>::SYNCMODE) == SYNCMODE_ONLYREAL )
      _syncLPRational(false);

    _syncRationalSolution();
    VectorRational& primal = _solRational._primal;
    assert(primal.dim() == numColsRational());

    DVectorRational activity(numRowsRational());
    _rationalLP->computePrimalActivity(primal, activity);
    maxviol = 0;
    sumviol = 0;

    for( int i = numRowsT() - 1; i >= 0; i-- )
      {
        Rational viol = lhsRational(i) - activity[i];
        if( viol > 0 )
          {
            sumviol += viol;
            if( viol > maxviol )
              {
                maxviol = viol;
                MSG_DEBUG( std::cout << "increased constraint violation for row " << i << ": " << rationalToString(viol)
                           << " lhs: " << rationalToString(lhsRational(i))
                           << ", activity: " << rationalToString(activity[i]) << "\n" );
              }
          }

        viol = activity[i] - rhsRational(i);
        if( viol > 0 )
          {
            sumviol += viol;
            if( viol > maxviol )
              {
                maxviol = viol;
                MSG_DEBUG( std::cout << "increased constraint violation for row " << i << ": " << rationalToString(viol)
                           << " rhs: " << rationalToString(rhsRational(i))
                           << ", activity: " << rationalToString(activity[i]) << "\n" );
              }
          }
      }

    return true;
  }



  /// gets violation of reduced costs; returns true on success
  template <>
	bool SoPlex<Real>::getRedCostViolationRational(Rational& maxviol, Rational& sumviol)
  {
    if( !isPrimalFeasible() || !isDualFeasible() )
      return false;

    // if we have to synchronize, we do not measure time, because this would affect the solving statistics
    if( intParam(SoPlex<Real>::SYNCMODE) == SYNCMODE_ONLYREAL )
      _syncLPRational(false);

    _syncRationalSolution();
    VectorRational& redcost = _solRational._redCost;
    assert(redcost.dim() == numColsT());

    maxviol = 0;
    sumviol = 0;

    for( int c = numColsT() - 1; c >= 0; c-- )
      {
        assert(!_hasBasis || basisColStatus(c) != SPxSolver<Real>::UNDEFINED);

        if( _colTypes[c] == RANGETYPE_FIXED )
          {
            assert(lowerRational(c) == upperRational(c));
            continue;
          }

        assert(!_hasBasis || basisColStatus(c) != SPxSolver<Real>::ON_LOWER || _solRational._primal[c] == lowerRational(c));
        assert(!_hasBasis || basisColStatus(c) != SPxSolver<Real>::ON_UPPER || _solRational._primal[c] == upperRational(c));
        assert(!_hasBasis || basisColStatus(c) != SPxSolver<Real>::FIXED || _solRational._primal[c] == lowerRational(c));
        assert(!_hasBasis || basisColStatus(c) != SPxSolver<Real>::FIXED || _solRational._primal[c] == upperRational(c));

        if( intParam(SoPlex<Real>::OBJSENSE) == OBJSENSE_MINIMIZE )
          {
            if( _solRational._primal[c] != upperRational(c) && redcost[c] < 0 )
              {
                sumviol += -redcost[c];
                if( redcost[c] < -maxviol )
                  {
                    MSG_DEBUG( std::cout << "increased reduced cost violation for column " << c << " not on upper bound: " << rationalToString(-redcost[c]) << "\n" );
                    maxviol = -redcost[c];
                  }
              }
            if( _solRational._primal[c] != lowerRational(c) && redcost[c] > 0 )
              {
                sumviol += redcost[c];
                if( redcost[c] > maxviol )
                  {
                    MSG_DEBUG( std::cout << "increased reduced cost violation for column " << c << " not on lower bound: " << rationalToString(redcost[c]) << "\n" );
                    maxviol = redcost[c];
                  }
              }
          }
        else
          {
            if( _solRational._primal[c] != upperRational(c) && redcost[c] > 0 )
              {
                sumviol += redcost[c];
                if( redcost[c] > maxviol )
                  {
                    MSG_DEBUG( std::cout << "increased reduced cost violation for column " << c << " not on upper bound: " << rationalToString(redcost[c]) << "\n" );
                    maxviol = redcost[c];
                  }
              }
            if( _solRational._primal[c] != lowerRational(c) && redcost[c] < 0 )
              {
                sumviol += -redcost[c];
                if( redcost[c] < -maxviol )
                  {
                    MSG_DEBUG( std::cout << "increased reduced cost violation for column " << c << " not on lower bound: " << rationalToString(-redcost[c]) << "\n" );
                    maxviol = -redcost[c];
                  }
              }
          }
      }

    return true;
  }



  /// gets violation of dual multipliers; returns true on success
  template <>
	bool SoPlex<Real>::getDualViolationRational(Rational& maxviol, Rational& sumviol)
  {
    if( !isDualFeasible() || !isPrimalFeasible() )
      return false;

    // if we have to synchronize, we do not measure time, because this would affect the solving statistics
    if( intParam(SoPlex<Real>::SYNCMODE) == SYNCMODE_ONLYREAL )
      _syncLPRational(false);

    _syncRationalSolution();
    VectorRational& dual = _solRational._dual;
    assert(dual.dim() == numRowsT());

    maxviol = 0;
    sumviol = 0;

    for( int r = numRowsT() - 1; r >= 0; r-- )
      {
        assert(!_hasBasis || basisRowStatus(r) != SPxSolver<Real>::UNDEFINED);

        if( _rowTypes[r] == RANGETYPE_FIXED )
          {
            assert(lhsRational(r) == rhsRational(r));
            continue;
          }

        assert(!_hasBasis || basisRowStatus(r) != SPxSolver<Real>::ON_LOWER || _solRational._slacks[r] <= lhsRational(r) + _rationalFeastol);
        assert(!_hasBasis || basisRowStatus(r) != SPxSolver<Real>::ON_UPPER || _solRational._slacks[r] >= rhsRational(r) - _rationalFeastol);
        assert(!_hasBasis || basisRowStatus(r) != SPxSolver<Real>::FIXED || _solRational._slacks[r] <= lhsRational(r) + _rationalFeastol);
        assert(!_hasBasis || basisRowStatus(r) != SPxSolver<Real>::FIXED || _solRational._slacks[r] >= rhsRational(r) - _rationalFeastol);

        if( intParam(SoPlex<Real>::OBJSENSE) == OBJSENSE_MINIMIZE )
          {
            if( _solRational._slacks[r] < rhsRational(r) - _rationalFeastol && dual[r] < 0 )
              {
                sumviol += -dual[r];
                if( dual[r] < -maxviol )
                  {
                    MSG_DEBUG( std::cout << "increased dual violation for row " << r << " not on upper bound: " << rationalToString(-dual[r])
                               << " (slack = " << rationalToString(_solRational._slacks[r])
                               << ", status = " << basisRowStatus(r)
                               << ", lhs = " << rationalToString(lhsRational(r))
                               << ", rhs = " << rationalToString(rhsRational(r)) << ")\n" );
                    maxviol = -dual[r];
                  }
              }
            if( _solRational._slacks[r] > lhsRational(r) + _rationalFeastol && dual[r] > 0 )
              {
                sumviol += dual[r];
                if( dual[r] > maxviol )
                  {
                    MSG_DEBUG( std::cout << "increased dual violation for row " << r << " not on lower bound: " << rationalToString(dual[r])
                               << " (slack = " << rationalToString(_solRational._slacks[r])
                               << ", status = " << basisRowStatus(r)
                               << ", lhs = " << rationalToString(lhsRational(r))
                               << ", rhs = " << rationalToString(rhsRational(r)) << ")\n" );
                    maxviol = dual[r];
                  }
              }
          }
        else
          {
            if( _solRational._slacks[r] < rhsRational(r) - _rationalFeastol && dual[r] > 0 )
              {
                sumviol += dual[r];
                if( dual[r] > maxviol )
                  {
                    MSG_DEBUG( std::cout << "increased dual violation for row " << r << " not on upper bound: " << rationalToString(dual[r])
                               << " (slack = " << rationalToString(_solRational._slacks[r])
                               << ", status = " << basisRowStatus(r)
                               << ", lhs = " << rationalToString(lhsRational(r))
                               << ", rhs = " << rationalToString(rhsRational(r)) << ")\n" );
                    maxviol = dual[r];
                  }
              }
            if( _solRational._slacks[r] > lhsRational(r) + _rationalFeastol && dual[r] < 0 )
              {
                sumviol += -dual[r];
                if( dual[r] < -maxviol )
                  {
                    MSG_DEBUG( std::cout << "increased dual violation for row " << r << " not on lower bound: " << rationalToString(-dual[r])
                               << " (slack = " << rationalToString(_solRational._slacks[r])
                               << ", status = " << basisRowStatus(r)
                               << ", lhs = " << rationalToString(lhsRational(r))
                               << ", rhs = " << rationalToString(rhsRational(r)) << ")\n" );
                    maxviol = -dual[r];
                  }
              }
          }
      }

    return true;
  }



#ifdef SOPLEX_WITH_GMP
  /// gets the primal solution vector if available; returns true on success
  template <>
	bool SoPlex<Real>::getPrimalRational(mpq_t* vector, const int size)
  {
    assert(size >= numColsT());

    if( hasPrimal() )
      {
        _syncRationalSolution();
        for( int i = 0; i < numColsT(); i++ )
          mpq_set(vector[i], _solRational._primal[i].getMpqRef());
        return true;
      }
    else
      return false;
  }


  /// gets the vector of slack values if available; returns true on success
  template <>
	bool SoPlex<Real>::getSlacksRational(mpq_t* vector, const int size)
  {
    assert(size >= numRowsT());

    if( hasPrimal() )
      {
        _syncRationalSolution();
        for( int i = 0; i < numRowsT(); i++ )
          mpq_set(vector[i], _solRational._slacks[i].getMpqRef());
        return true;
      }
    else
      return false;
  }



  /// gets the primal ray if LP is unbounded; returns true on success
  template <>
	bool SoPlex<Real>::getPrimalRayRational(mpq_t* vector, const int size)
  {
    assert(size >= numColsT());

    if( hasPrimalRay() )
      {
        _syncRationalSolution();
        for( int i = 0; i < numColsT(); i++ )
          mpq_set(vector[i], _solRational._primalRay[i].getMpqRef());
        return true;
      }
    else
      return false;
  }



  /// gets the dual solution vector if available; returns true on success
  template <>
	bool SoPlex<Real>::getDualRational(mpq_t* vector, const int size)
  {
    assert(size >= numRowsT());

    if( hasDual() )
      {
        _syncRationalSolution();
        for( int i = 0; i < numRowsT(); i++ )
          mpq_set(vector[i], _solRational._dual[i].getMpqRef());
        return true;
      }
    else
      return false;
  }



  /// gets the vector of reduced cost values if available; returns true on success
  template <>
	bool SoPlex<Real>::getRedCostRational(mpq_t* vector, const int size)
  {
    assert(size >= numColsT());

    if( hasDual() )
      {
        _syncRationalSolution();
        for( int i = 0; i < numColsT(); i++ )
          mpq_set(vector[i], _solRational._redCost[i].getMpqRef());
        return true;
      }
    else
      return false;
  }



  /// gets the Farkas proof if LP is infeasible; returns true on success
  template <>
	bool SoPlex<Real>::getDualFarkasRational(mpq_t* vector, const int size)
  {
    assert(size >= numRowsT());

    if( hasDualFarkas() )
      {
        _syncRationalSolution();
        for( int i = 0; i < numRowsT(); i++ )
          mpq_set(vector[i], _solRational._dualFarkas[i].getMpqRef());
        return true;
      }
    else
      return false;
  }
#endif



  /// get size of primal solution
  template <>
  int SoPlex<Real>::totalSizePrimalRational(const int base)
  {
    if( hasPrimal() || hasPrimalRay() )
      {
        _syncRationalSolution();
        return _solRational.totalSizePrimal(base);
      }
    else
      return 0;
  }



  /// get size of dual solution
  template <>
  int SoPlex<Real>::totalSizeDualRational(const int base)
  {
    if( hasDual() || hasDualFarkas() )
      {
        _syncRationalSolution();
        return _solRational.totalSizeDual(base);
      }
    else
      return 0;
  }



  /// get size of least common multiple of denominators in primal solution
  template <>
  int SoPlex<Real>::dlcmSizePrimalRational(const int base)
  {
    if( hasPrimal() || hasPrimalRay() )
      {
        _syncRationalSolution();
        return _solRational.dlcmSizePrimal(base);
      }
    else
      return 0;
  }



  /// get size of least common multiple of denominators in dual solution
  template <>
  int SoPlex<Real>::dlcmSizeDualRational(const int base)
  {
    if( hasDual() || hasDualFarkas() )
      {
        _syncRationalSolution();
        return _solRational.dlcmSizeDual(base);
      }
    else
      return 0;
  }



  /// get size of largest denominator in primal solution
  template <>
  int SoPlex<Real>::dmaxSizePrimalRational(const int base)
  {
    if( hasPrimal() || hasPrimalRay() )
      {
        _syncRationalSolution();
        return _solRational.dmaxSizePrimal(base);
      }
    else
      return 0;
  }



  /// get size of largest denominator in dual solution
  template <>
	int SoPlex<Real>::dmaxSizeDualRational(const int base)
  {
    if( hasDual() || hasDualFarkas() )
      {
        _syncRationalSolution();
        return _solRational.dmaxSizeDual(base);
      }
    else
      return 0;
  }

  /// returns the current basis status
  template <>
  typename SPxBasis<Real>::SPxStatus SoPlex<Real>::basisStatus() const
  {
    if( !hasBasis() )
      return SPxBasis<Real>::NO_PROBLEM;
    else if( status() == SPxSolver<Real>::OPTIMAL )
      return SPxBasis<Real>::OPTIMAL;
    else if( status() == SPxSolver<Real>::UNBOUNDED )
      return SPxBasis<Real>::UNBOUNDED;
    else if( status() == SPxSolver<Real>::INFEASIBLE )
      return SPxBasis<Real>::INFEASIBLE;
    else if( hasPrimal() )
      return SPxBasis<Real>::PRIMAL;
    else if( hasDual() )
      return SPxBasis<Real>::DUAL;
    else
      return SPxBasis<Real>::REGULAR;
  }



  /// returns basis status for a single row
  template <>
  typename SPxSolver<Real>::VarStatus SoPlex<Real>::basisRowStatus(int row) const
  {
    assert(row >= 0);
    assert(row < numRowsT());

    // if no basis is available, return slack basis; if index is out of range, return basic status as for a newly
    // added row
    if( !hasBasis() || row < 0 || row >= numRowsT() )
      return SPxSolver<Real>::BASIC;
    // if the real LP is loaded, ask solver
    else if( _isRealLPLoaded )
      {
        return _solver.getBasisRowStatus(row);
      }
    // if the real LP is not loaded, the basis is stored in the basis arrays of this class
    else
      {
        assert(row < _basisStatusRows.size());
        return _basisStatusRows[row];
      }
  }



  /// returns basis status for a single column
  template <>
  typename SPxSolver<Real>::VarStatus SoPlex<Real>::basisColStatus(int col) const
  {
    assert(col >= 0);
    assert(col < numColsT());

    // if index is out of range, return nonbasic status as for a newly added unbounded column
    if( col < 0 || col >= numColsT() )
      {
        return SPxSolver<Real>::ZERO;
      }
    // if no basis is available, return slack basis
    else if( !hasBasis() )
      {
        if( lowerReal(col) > -realParam(SoPlex<Real>::INFTY) )
          return SPxSolver<Real>::ON_LOWER;
        else if( upperReal(col) < realParam(SoPlex<Real>::INFTY) )
          return SPxSolver<Real>::ON_UPPER;
        else
          return SPxSolver<Real>::ZERO;
      }
    // if the real LP is loaded, ask solver
    else if( _isRealLPLoaded )
      {
        return _solver.getBasisColStatus(col);
      }
    // if the real LP is not loaded, the basis is stored in the basis arrays of this class
    else
      {
        assert(col < _basisStatusCols.size());
        return _basisStatusCols[col];
      }
  }



  /// gets current basis
  /// #template #temp
  template <>
  void SoPlex<Real>::getBasis(typename SPxSolver<Real>::VarStatus rows[], typename SPxSolver<Real>::VarStatus cols[]) const
  {
    // if no basis is available, return slack basis
    if( !hasBasis() )
      {
        for( int i = numRowsT() - 1; i >= 0; i-- )
          rows[i] = SPxSolver<Real>::BASIC;

        for( int i = numColsT() - 1; i >= 0; i-- )
          {
            if( lowerReal(i) > -realParam(SoPlex<Real>::INFTY) )
              cols[i] = SPxSolver<Real>::ON_LOWER;
            else if( upperReal(i) < realParam(SoPlex<Real>::INFTY) )
              cols[i] = SPxSolver<Real>::ON_UPPER;
            else
              cols[i] = SPxSolver<Real>::ZERO;
          }
      }
    // if the real LP is loaded, ask solver
    else if( _isRealLPLoaded )
      {
        (void)_solver.getBasis(rows, cols);
      }
    // if the real LP is not loaded, the basis is stored in the basis arrays of this class
    else
      {
        assert(numRowsT() == _basisStatusRows.size());
        assert(numColsT() == _basisStatusCols.size());

        for( int i = numRowsT() - 1; i >= 0; i-- )
          rows[i] = _basisStatusRows[i];

        for( int i = numColsT() - 1; i >= 0; i-- )
          cols[i] = _basisStatusCols[i];
      }
  }



  /// returns the indices of the basic columns and rows; basic column n gives value n, basic row m gives value -1-m
  template <>
  void SoPlex<Real>::getBasisInd(int* bind) const
  {
    // if no basis is available, return slack basis
    if( !hasBasis() )
      {
        for( int i = 0; i < numRowsT(); ++i )
          bind[i] = -1 - i;
      }
    // if the real LP is not loaded, the basis is stored in the basis arrays of this class
    else if( !_isRealLPLoaded )
      {
        int k = 0;

        assert(numRowsT() == _basisStatusRows.size());
        assert(numColsT() == _basisStatusCols.size());

        for( int i = 0; i < numRowsT(); ++i )
          {
            if( _basisStatusRows[i] == SPxSolver<Real>::BASIC )
              {
                bind[k] = -1 - i;
                k++;
              }
          }

        for( int j = 0; j < numColsT(); ++j )
          {
            if( _basisStatusCols[j] == SPxSolver<Real>::BASIC )
              {
                bind[k] = j;
                k++;
              }
          }

        assert(k == numRowsT());
      }
    // if the real LP is loaded, the basis is stored in the solver and we need to distinguish between column and row
    // representation; ask the solver itself which representation it has, since the REPRESENTATION parameter of this
    // class might be set to automatic
    else if( _solver.rep() == SPxSolver<Real>::COLUMN )
      {
        for( int i = 0; i < numRowsT(); ++i )
          {
            SPxId id = _solver.basis().baseId(i);
            bind[i] = (id.isSPxColId() ? _solver.number(id) : - 1 - _solver.number(id));
          }
      }
    // for row representation, return the complement of the row basis; for this, we need to loop through all rows and columns
    else
      {
        assert(_solver.rep() == SPxSolver<Real>::ROW);

        int k = 0;

        for( int i = 0; i < numRowsT(); ++i )
          {
            if( !_solver.isRowBasic(i) )
              {
                bind[k] = -1 - i;
                k++;
              }
          }

        for( int j = 0; j < numColsT(); ++j )
          {
            if( !_solver.isColBasic(j) )
              {
                bind[k] = j;
                k++;
              }
          }

        assert(k == numRowsT());
      }
  }


  /// compute condition number estimate based on the diagonal of the LU factorization; returns true on success
  /// type = 0: max/min ratio
  /// type = 1: trace of U (sum of diagonal elements)
  /// type = 2: product of diagonal elements
  /// #template #temp 
  template <>
	bool SoPlex<Real>::getFastCondition(Real& condition, int type)
  {
    _ensureRealLPLoaded();
    if( !_isRealLPLoaded )
      return false;

    if( _solver.basis().status() == SPxBasis<Real>::NO_PROBLEM )
      return false;

    condition = _solver.basis().getFastCondition(type);

    return true;
  }

  /// computes an estimated condition number for the current basis matrix using the power method; returns true on success
  template <>
	bool SoPlex<Real>::getEstimatedCondition(Real& condition)
  {
    _ensureRealLPLoaded();
    if( !_isRealLPLoaded )
      return false;

    if( _solver.basis().status() == SPxBasis<Real>::NO_PROBLEM )
      return false;

    condition = _solver.basis().getEstimatedCondition();

    return true;
  }

  /// computes the exact condition number for the current basis matrix using the power method; returns true on success
  template <>
	bool SoPlex<Real>::getExactCondition(Real& condition)
  {
    _ensureRealLPLoaded();
    if( !_isRealLPLoaded )
      return false;

    if( _solver.basis().status() == SPxBasis<Real>::NO_PROBLEM )
      return false;

    condition = _solver.basis().getExactCondition();

    return true;
  }

  /// computes row r of basis inverse; returns true on success
  template <>
	bool SoPlex<Real>::getBasisInverseRowReal(int r, Real* coef, int* inds, int* ninds, bool unscale)
  {
    assert(r >= 0);
    assert(r < numRowsT());
    assert(coef != 0);

    if( !hasBasis() || r < 0 || r >= numRowsT() )
      return false;

    _ensureRealLPLoaded();

    if( !_isRealLPLoaded )
      return false;

    // we need to distinguish between column and row representation; ask the solver itself which representation it
    // has, since the REPRESENTATION parameter of this class might be set to automatic
    if( _solver.rep() == SPxSolver<Real>::COLUMN )
      {
        int idx;
        SSVectorReal x(numRowsT());
        try
          {
            /* unscaling required? */
            if( unscale && _solver.isScaled())
              {
                /* for information on the unscaling procedure see spxscaler.h */

                int scaleExp;
                DSVector rhs(_solver.unitVector(r));

                if( _solver.basis().baseId(r).isSPxColId() )
                  scaleExp = _scaler->getColScaleExp(_solver.number(_solver.basis().baseId(r)));
                else
                  scaleExp = - _scaler->getRowScaleExp(_solver.number(_solver.basis().baseId(r)));

                rhs *= spxLdexp(1.0, scaleExp);

                _solver.basis().coSolve(x, rhs);
                x.setup();
                int size = x.size();

                for( int i = 0; i < size; i++ )
                  {
                    scaleExp = _scaler->getRowScaleExp(x.index(i));
                    x.scaleValue(x.index(i), scaleExp);
                  }
              }
            else
              {
                _solver.basis().coSolve(x, _solver.unitVector(r));
              }
          }
        catch( const SPxException& E )
          {
            MSG_INFO1( spxout, spxout << "Caught exception <" << E.what() << "> while computing basis inverse row.\n" );
            return false;
          }
        // copy sparse data to dense result vector based on coef array
        if( ninds != 0 && inds != 0 )
          {
            // during solving SoPlex may have destroyed the sparsity structure so we need to restore it
            x.setup();
            *ninds = x.size();
            for( int i = 0; i < *ninds; ++i )
              {
                idx = x.index(i);
                coef[idx] = x[idx];
                // set sparsity pattern of coef array
                inds[i] = idx;
              }
          }
        else
          {
            VectorReal y(numRowsT(), coef);
            y = x;
            if( ninds != NULL )
              *ninds = -1;
          }
      }
    else
      {
        assert(_solver.rep() == SPxSolver<Real>::ROW);

        // @todo should rhs be a reference?
        DSVector rhs(numColsT());
        SSVector y(numColsT());
        int* bind = 0;
        int index;

        // get ordering of column basis matrix
        spx_alloc(bind, numRowsT());
        getBasisInd(bind);

        // get vector corresponding to requested index r
        index = bind[r];

        // r corresponds to a row vector
        if( index < 0 )
          {
            // transform index to actual row index
            index = -index - 1;

            // should be a valid row index and in the column basis matrix, i.e., not basic w.r.t. row representation
            assert(index >= 0);
            assert(index < numRowsT());
            assert(!_solver.isRowBasic(index));

            // get row vector
            rhs = _solver.rowVector(index);
            rhs *= -1.0;

            if( unscale && _solver.isScaled() )
              {
                for( int i = 0; i < rhs.size(); ++i)
                  rhs.value(i) = spxLdexp(rhs.value(i), -_scaler->getRowScaleExp(index));
              }
          }
        // r corresponds to a column vector
        else
          {
            // should be a valid column index and in the column basis matrix, i.e., not basic w.r.t. row representation
            assert(index < numColsT());
            assert(!_solver.isColBasic(index));

            // get unit vector
            rhs = UnitVectorReal(index);

            if( unscale && _solver.isScaled() )
              rhs *= spxLdexp(1.0, _scaler->getColScaleExp(index));
          }

        // solve system "y B = rhs", where B is the row basis matrix
        try
          {
            _solver.basis().solve(y, rhs);
          }
        catch( const SPxException& E )
          {
            MSG_INFO1( spxout, spxout << "Caught exception <" << E.what() << "> while computing basis inverse row.\n" );
            return false;
          }

        // initialize result vector x as zero
        memset(coef, 0, (unsigned int)numRowsT() * sizeof(Real));

        // add nonzero entries
        for( int i = 0; i < numColsT(); ++i )
          {
            SPxId id = _solver.basis().baseId(i);

            if( id.isSPxRowId() )
              {
                assert(_solver.number(id) >= 0);
                assert(_solver.number(id) < numRowsT());
                assert(bind[r] >= 0 || _solver.number(id) != index);

                int rowindex = _solver.number(id);
                coef[rowindex] = y[i];

                if( unscale && _solver.isScaled() )
                  coef[rowindex] = spxLdexp(y[i], _scaler->getRowScaleExp(rowindex));
              }
          }

        // if r corresponds to a row vector, we have to add a 1 at position r
        if( bind[r] < 0 )
          {
            assert(coef[index] == 0.0);
            coef[index] = 1.0;
          }

        // @todo implement returning of sparsity information like in column wise case
        if( ninds != NULL)
          *ninds = -1;

        // free memory
        spx_free(bind);
      }

    return true;
  }



  /// computes column c of basis inverse; returns true on success
  /// @todo does not work correctly for the row representation
  template <>
	bool SoPlex<Real>::getBasisInverseColReal(int c, Real* coef, int* inds, int* ninds, bool unscale)
  {
    assert(c >= 0);
    assert(c < numRowsT());
    assert(coef != 0);

    if( !hasBasis() || c < 0 || c >= numRowsT() )
      return false;

    _ensureRealLPLoaded();

    if( !_isRealLPLoaded )
      return false;

    // we need to distinguish between column and row representation; ask the solver itself which representation it
    // has, since the REPRESENTATION parameter of this class might be set to automatic
    if( _solver.rep() == SPxSolver<Real>::COLUMN )
      {
        int idx;
        SSVectorReal x(numRowsT());
        try
          {
            /* unscaling required? */
            if( unscale && _solver.isScaled())
              {
                /* for information on the unscaling procedure see spxscaler.h */

                int scaleExp =_scaler->getRowScaleExp(c);
                DSVector rhs(_solver.unitVector(c));
                rhs *= spxLdexp(1.0, scaleExp);

                _solver.basis().solve(x, rhs);

                x.setup();
                int size = x.size();

                for( int i = 0; i < size; i++ )
                  {
                    if( _solver.basis().baseId(x.index(i)).isSPxColId() )
                      {
                        idx = _solver.number(_solver.basis().baseId(x.index(i)));
                        scaleExp = _scaler->getColScaleExp(idx);
                        x.scaleValue(x.index(i), scaleExp);
                      }
                    else
                      {
                        idx = _solver.number(_solver.basis().baseId(x.index(i)));
                        scaleExp = - _scaler->getRowScaleExp(idx);
                        x.scaleValue(x.index(i), scaleExp);
                      }
                  }
              }
            else
              {
                _solver.basis().solve(x, _solver.unitVector(c));
              }
          }
        catch( const SPxException& E )
          {
            MSG_INFO1( spxout, spxout << "Caught exception <" << E.what() << "> while computing basis inverse row.\n" );
            return false;
          }
        // copy sparse data to dense result vector based on coef array
        if( ninds != 0 && inds != 0 )
          {
            // SoPlex may have destroyed the sparsity structure so we need to restore it
            x.setup();
            *ninds = x.size();
            for( int i = 0; i < *ninds; ++i )
              {
                idx = x.index(i);
                coef[idx] = x[idx];
                // set sparsity pattern of coef array
                inds[i] = idx;
              }
          }
        else
          {
            VectorReal y(numRowsT(), coef);
            y = x;
            if( ninds != 0 )
              *ninds = -1;
          }
      }
    else
      {
        assert(_solver.rep() == SPxSolver<Real>::ROW);

        // @todo should rhs be a reference?
        DSVectorReal rhs(numColsT());
        SSVectorReal y(numColsT());
        int* bind = 0;
        int index;

        // get ordering of column basis matrix
        spx_alloc(bind, numRowsT());
        getBasisInd(bind);

        // get vector corresponding to requested index c
        index = bind[c];

        // c corresponds to a row vector
        if( index < 0 )
          {
            // transform index to actual row index
            index = -index - 1;

            // should be a valid row index and in the column basis matrix, i.e., not basic w.r.t. row representation
            assert(index >= 0);
            assert(index < numRowsT());
            assert(!_solver.isRowBasic(index));

            // get row vector
            rhs = _solver.rowVector(index);
            rhs *= -1.0;
          }
        // c corresponds to a column vector
        else
          {
            // should be a valid column index and in the column basis matrix, i.e., not basic w.r.t. row representation
            assert(index < numColsT());
            assert(!_solver.isColBasic(index));

            // get unit vector
            rhs = UnitVectorReal(index);
          }

        // solve system "y B = rhs", where B is the row basis matrix
        try
          {
            /* unscaling required? */
            if( unscale && _solver.isScaled() )
              {
                int size = rhs.size();
                int scaleExp;

                for( int i = 0; i < size; i++ )
                  {
                    scaleExp = _scaler->getColScaleExp(i);
                    rhs.value(i) *= spxLdexp(1.0, scaleExp);
                  }

                _solver.basis().coSolve(y, rhs);

                int rowIdx;
                size = y.size();

                for( int i = 0; i < size; i++ )
                  {
                    assert(_solver.basis().baseId(y.index(i)).isSPxRowId());
                    rowIdx = _solver.basis().baseId(y.index(i)).getIdx();
                    scaleExp = _scaler->getRowScaleExp(rowIdx);
                    y.setValue(i, y.value(i) * spxLdexp(1.0, scaleExp));
                  }
              }
            else
              {
                _solver.basis().coSolve(y, rhs);
              }
          }
        catch( const SPxException& E )
          {
            MSG_INFO1( spxout, spxout << "Caught exception <" << E.what() << "> while computing basis inverse row.\n" );
            return false;
          }

        // initialize result vector x as zero
        memset(coef, 0, (unsigned int)numRowsT() * sizeof(Real));

        // add nonzero entries
        for( int i = 0; i < numColsT(); ++i )
          {
            SPxId id = _solver.basis().baseId(i);

            if( id.isSPxRowId() )
              {
                assert(_solver.number(id) >= 0);
                assert(_solver.number(id) < numRowsT());
                assert(bind[c] >= 0 || _solver.number(id) != index);

                coef[_solver.number(id)] = y[i];
              }
          }

        // if c corresponds to a row vector, we have to add a 1 at position c
        if( bind[c] < 0 )
          {
            assert(coef[index] == 0.0);
            coef[index] = 1.0;
          }

        // @todo implement returning of sparsity information like in column wise case
        if( ninds != NULL)
          *ninds = -1;

        // free memory
        spx_free(bind);
      }

    return true;
  }



  /// computes dense solution of basis matrix B * sol = rhs; returns true on success
  template <>
	bool SoPlex<Real>::getBasisInverseTimesVecReal(Real* rhs, Real* sol, bool unscale)
  {
    VectorReal v(numRowsT(), rhs);
    VectorReal x(numRowsT(), sol);

    if( !hasBasis() )
      return false;

    _ensureRealLPLoaded();

    if( !_isRealLPLoaded )
      return false;

    // we need to distinguish between column and row representation; ask the solver itself which representation it
    // has, since the REPRESENTATION parameter of this class might be set to automatic; in the column case we can use
    // the existing factorization
    if( _solver.rep() == SPxSolver<Real>::COLUMN )
      {
        // solve system "x = B^-1 * v"
        try
          {
            /* unscaling required? */
            if( unscale && _solver.isScaled())
              {
                /* for information on the unscaling procedure see spxscaler.h */
                int scaleExp;
                int idx;

                for( int i = 0; i < v.dim(); ++i)
                  {
                    if( isNotZero(v[i]) )
                      {
                        scaleExp =_scaler->getRowScaleExp(i);
                        v[i] = spxLdexp(v[i], scaleExp);
                      }
                  }

                _solver.basis().solve(x, v);

                for( int i = 0; i < x.dim(); i++ )
                  {
                    if( isNotZero(x[i]) )
                      {
                        idx = _solver.number(_solver.basis().baseId(i));
                        if( _solver.basis().baseId(i).isSPxColId() )
                          scaleExp = _scaler->getColScaleExp(idx);
                        else
                          scaleExp = - _scaler->getRowScaleExp(idx);
                        x[i] = spxLdexp(x[i], scaleExp);
                      }
                  }
              }
            else
              {
                _solver.basis().solve(x, v);
              }
          }
        catch( const SPxException& E )
          {
            MSG_INFO1( spxout, spxout << "Caught exception <" << E.what() << "> while solving with basis matrix.\n" );
            return false;
          }
      }
    else
      {
        assert(_solver.rep() == SPxSolver<Real>::ROW);

        DSVectorReal rowrhs(numColsT());
        SSVectorReal y(numColsT());
        int* bind = 0;

        bool adaptScaling = unscale && _realLP->isScaled();
        int scaleExp;
        int idx;

        // get ordering of column basis matrix
        spx_alloc(bind, numRowsT());
        getBasisInd(bind);

        // fill right-hand side for row-based system
        for( int i = 0; i < numColsT(); ++i )
          {
            SPxId id = _solver.basis().baseId(i);

            if( id.isSPxRowId() )
              {
                assert(_solver.number(id) >= 0);
                assert(_solver.number(id) < numRowsT());

                if( adaptScaling )
                  {
                    idx = _solver.number(id);
                    scaleExp = _scaler->getRowScaleExp(idx);
                    rowrhs.add(i, spxLdexp(v[idx], scaleExp));
                  }
                else
                  rowrhs.add(i, v[_solver.number(id)]);
              }
            else
              {
                assert(rowrhs[i] == 0.0);
              }
          }

        // solve system "B y = rowrhs", where B is the row basis matrix
        try
          {
            _solver.basis().coSolve(y, rowrhs);
          }
        catch( const SPxException& E )
          {
            MSG_INFO1( spxout, spxout << "Caught exception <" << E.what() << "> while solving with basis matrix.\n" );
            return false;
          }

        // fill result w.r.t. order given by bind
        for( int i = 0; i < numRowsT(); ++i )
          {
            int index;

            index = bind[i];

            if( index < 0 )
              {
                index = -index-1;

                // should be a valid row index and in the column basis matrix, i.e., not basic w.r.t. row representation
                assert(index >= 0);
                assert(index < numRowsT());
                assert(!_solver.isRowBasic(index));

                x[i] = v[index] - (rowVectorRealInternal(index) * Vector(numColsT(), y.get_ptr()));

                if( adaptScaling )
                  {
                    scaleExp = -_scaler->getRowScaleExp(index);
                    x[i] = spxLdexp(x[i], scaleExp);
                  }
              }
            else
              {
                // should be a valid column index and in the column basis matrix, i.e., not basic w.r.t. row representation
                assert(index >= 0);
                assert(index < numColsT());
                assert(!_solver.isColBasic(index));

                if( adaptScaling )
                  {
                    scaleExp = _scaler->getColScaleExp(index);
                    x[i] = spxLdexp(y[index], scaleExp);
                  }
                else
                  x[i] = y[index];
              }
          }

        // free memory
        spx_free(bind);
      }
    return true;
  }



  /// multiply with basis matrix; B * vec (inplace)
  template <>
	bool SoPlex<Real>::multBasis(Real* vec, bool unscale)
  {
    if( !hasBasis() )
      return false;

    _ensureRealLPLoaded();

    if( !_isRealLPLoaded )
      return false;

    if( _solver.rep() == SPxSolver<Real>::COLUMN )
      {
        int basisdim = numRowsT();

        // create Vector from input values
        Vector x(basisdim, vec);

        if( unscale && _solver.isScaled() )
          {
            /* for information on the unscaling procedure see spxscaler.h */

            int scaleExp;
            for( int i = 0; i < basisdim; ++i)
              {
                if( isNotZero(vec[i]) )
                  {
                    if( _solver.basis().baseId(i).isSPxColId() )
                      scaleExp = - _scaler->getColScaleExp(_solver.number(_solver.basis().baseId(i)));
                    else
                      scaleExp = _scaler->getRowScaleExp(_solver.number(_solver.basis().baseId(i)));

                    vec[i] = spxLdexp(vec[i], scaleExp);
                  }
              }
            _solver.basis().multBaseWith(x);
            for( int i = 0; i < basisdim; ++i)
              {
                scaleExp = _scaler->getRowScaleExp(i);
                vec[i] = spxLdexp(vec[i], -scaleExp);
              }
          }
        else
          _solver.basis().multBaseWith(x);
      }
    else
      {
        int colbasisdim = numRowsT();

        DSVector y(colbasisdim);

        y.clear();

        // create Vector from input values
        Vector x(colbasisdim, vec);

        int* bind = 0;
        int index;

        // get ordering of column basis matrix
        spx_alloc(bind, colbasisdim);
        getBasisInd(bind);

        // temporarily create the column basis and multiply every column with x
        for( int i = 0; i < colbasisdim; ++i)
          {
            if( isNotZero(x[i]) )
              {
                // get vector corresponding to requested index i
                index = bind[i];
                // r corresponds to a row vector
                if( index < 0 )
                  {
                    // transform index to actual row index
                    index = -index - 1;

                    // should be a valid row index and in the column basis matrix, i.e., not basic w.r.t. row representation
                    assert(index >= 0);
                    assert(index < numRowsT());
                    assert(!_solver.isRowBasic(index));

                    y.add(x[i] * UnitVectorReal(index));
                  }
                // r corresponds to a column vector
                else
                  {
                    // should be a valid column index and in the column basis matrix, i.e., not basic w.r.t. row representation
                    assert(index < numColsT());
                    assert(!_solver.isColBasic(index));

                    if( unscale && _solver.isScaled() )
                      {
                        DSVector col;
                        _solver.getColVectorUnscaled(index, col);
                        y.add(x[i] * col);
                      }

                    y.add(x[i] * _solver.colVector(index));
                  }
              }
          }
        spx_free(bind);
        x = y;
      }

    return true;
  }



  /// multiply with transpose of basis matrix; vec * B^T (inplace)
  template <>
	bool SoPlex<Real>::multBasisTranspose(Real* vec, bool unscale)
  {
    if( !hasBasis() )
      return false;

    _ensureRealLPLoaded();

    if( !_isRealLPLoaded )
      return false;

    if( _solver.rep() == SPxSolver<Real>::COLUMN )
      {
        int basisdim = numRowsT();

        // create Vector from input values
        Vector x(basisdim, vec);

        if( unscale && _solver.isScaled() )
          {
            /* for information on the unscaling procedure see spxscaler.h */

            int scaleExp;
            for( int i = 0; i < basisdim; ++i)
              {
                if( isNotZero(vec[i]) )
                  {
                    scaleExp = - _scaler->getRowScaleExp(i);
                    vec[i] = spxLdexp(vec[i], scaleExp);
                  }
              }
            _solver.basis().multWithBase(x);
            for( int i = 0; i < basisdim; ++i)
              {
                if( isNotZero(vec[i]) )
                  {
                    if( _solver.basis().baseId(i).isSPxColId() )
                      scaleExp = - _scaler->getColScaleExp(_solver.number(_solver.basis().baseId(i)));
                    else
                      scaleExp = _scaler->getRowScaleExp(_solver.number(_solver.basis().baseId(i)));

                    vec[i] = spxLdexp(vec[i], scaleExp);
                  }
              }
          }
        else
          _solver.basis().multWithBase(x);
      }
    else
      {
        int colbasisdim = numRowsT();

        DSVector y(colbasisdim);

        // create Vector from input values
        Vector x(colbasisdim, vec);

        int* bind = 0;
        int index;

        // get ordering of column basis matrix
        spx_alloc(bind, colbasisdim);
        getBasisInd(bind);

        // temporarily create the column basis and multiply every column with x
        for( int i = 0; i < colbasisdim; ++i)
          {
            // get vector corresponding to requested index i
            index = bind[i];
            // r corresponds to a row vector
            if( index < 0 )
              {
                // transform index to actual row index
                index = -index - 1;

                // should be a valid row index and in the column basis matrix, i.e., not basic w.r.t. row representation
                assert(index >= 0);
                assert(index < numRowsT());
                assert(!_solver.isRowBasic(index));

                y.add(i, x * UnitVectorReal(index));
              }
            // r corresponds to a column vector
            else
              {
                // should be a valid column index and in the column basis matrix, i.e., not basic w.r.t. row representation
                assert(index < numColsT());
                assert(!_solver.isColBasic(index));

                if( unscale && _solver.isScaled() )
                  {
                    DSVector col;
                    _solver.getColVectorUnscaled(index, col);
                    y.add(i, x * col);
                  }
                else
                  y.add(i, x * _solver.colVector(index));
              }
          }
        spx_free(bind);
        x = y;
      }

    return true;
  }



  /// compute rational basis inverse; returns true on success
  template <>
	bool SoPlex<Real>::computeBasisInverseRational()
  {
    if( !hasBasis() )
      {
        _rationalLUSolver.clear();
        assert(_rationalLUSolver.status() == SLinSolverRational::UNLOADED);
        return false;
      }

    if( _rationalLUSolver.status() == SLinSolverRational::UNLOADED
        || _rationalLUSolver.status() == SLinSolverRational::TIME )
      {
        _rationalLUSolverBind.reSize(numRowsT());
        getBasisInd(_rationalLUSolverBind.get_ptr());
        _computeBasisInverseRational();
      }

    if( _rationalLUSolver.status() == SLinSolverRational::OK )
      return true;

    return false;
  }



  /// gets an array of indices for the columns of the rational basis matrix; bind[i] >= 0 means that the i-th column of
  /// the basis matrix contains variable bind[i]; bind[i] < 0 means that the i-th column of the basis matrix contains
  /// the slack variable for row -bind[i]-1; performs rational factorization if not available; returns true on success
  template <>
	bool SoPlex<Real>::getBasisIndRational(DataArray<int>& bind)
  {
    if( _rationalLUSolver.status() != SLinSolverRational::OK )
      computeBasisInverseRational();

    if( _rationalLUSolver.status() != SLinSolverRational::OK )
      return false;

    bind = _rationalLUSolverBind;
    assert(bind.size() == numRowsT());
    return true;
  }



  /// computes row r of basis inverse; performs rational factorization if not available; returns true on success
  template <>
	bool SoPlex<Real>::getBasisInverseRowRational(const int r, SSVectorRational& vec)
  {
    if( _rationalLUSolver.status() != SLinSolverRational::OK )
      computeBasisInverseRational();

    if( _rationalLUSolver.status() != SLinSolverRational::OK )
      return false;

    try
      {
        vec.reDim(numRowsT());
        _rationalLUSolver.solveLeft(vec, *_unitVectorRational(r));
      }
    catch( const SPxException& E )
      {
        MSG_INFO1( spxout, spxout << "Caught exception <" << E.what() << "> while computing rational basis inverse row.\n" );
        return false;
      }
    return true;
  }



  /// computes column c of basis inverse; performs rational factorization if not available; returns true on success
  template <>
	bool SoPlex<Real>::getBasisInverseColRational(const int c, SSVectorRational& vec)
  {
    if( _rationalLUSolver.status() != SLinSolverRational::OK )
      computeBasisInverseRational();

    if( _rationalLUSolver.status() != SLinSolverRational::OK )
      return false;

    try
      {
        vec.reDim(numRowsT());
        _rationalLUSolver.solveRight(vec, *_unitVectorRational(c));
      }
    catch( const SPxException& E )
      {
        MSG_INFO1( spxout, spxout << "Caught exception <" << E.what() << "> while computing rational basis inverse column.\n" );
        return false;
      }
    return true;
  }



  /// computes solution of basis matrix B * sol = rhs; performs rational factorization if not available; returns true
  /// on success
  template <>
	bool SoPlex<Real>::getBasisInverseTimesVecRational(const SVectorRational& rhs, SSVectorRational& sol)
  {
    if( _rationalLUSolver.status() != SLinSolverRational::OK )
      computeBasisInverseRational();

    if( _rationalLUSolver.status() != SLinSolverRational::OK )
      return false;

    try
      {
        sol.reDim(numRowsT());
        _rationalLUSolver.solveRight(sol, rhs);
      }
    catch( const SPxException& E )
      {
        MSG_INFO1( spxout, spxout << "Caught exception <" << E.what() << "> during right solve with rational basis inverse.\n" );
        return false;
      }
    return true;
  }



  /// sets starting basis via arrays of statuses
  template <>
  void SoPlex<Real>::setBasis(const typename SPxSolver<Real>::VarStatus rows[], const typename SPxSolver<Real>::VarStatus cols[])
  {
    _rationalLUSolver.clear();

    if( _isRealLPLoaded )
      {
        assert(numRowsT() == _solver.nRows());
        assert(numColsT() == _solver.nCols());

        _solver.setBasis(rows, cols);
        _hasBasis = (_solver.basis().status() > SPxBasis<Real>::NO_PROBLEM);
      }
    else
      {
        _basisStatusRows.reSize(numRowsT());
        _basisStatusCols.reSize(numColsT());

        for( int i = numRowsT() - 1; i >= 0; i-- )
          _basisStatusRows[i] = rows[i];

        for( int j = numColsT() - 1; j >= 0; j-- )
          _basisStatusCols[j] = cols[j];

        _hasBasis = true;
      }
  }



  /// clears starting basis
  template <>
  void SoPlex<Real>::clearBasis()
  {
    _solver.reLoad();
    _status = _solver.status();
    _hasBasis = false;
    _rationalLUSolver.clear();
  }



  /// number of iterations since last call to solve
  template <>
	int SoPlex<Real>::numIterations() const
  {
    return _statistics->iterations;
  }



  /// time spent in last call to solve
  template <>
  Real SoPlex<Real>::solveTime() const
  {
    return _statistics->solvingTime->time();
  }



  /// statistical information in form of a string
  template <>
  std::string SoPlex<Real>::statisticString() const
  {
    std::stringstream s;
    s  << "Factorizations     : " << std::setw(10) << _statistics->luFactorizationsReal << std::endl
       << "  Time spent       : " << std::setw(10) << std::fixed << std::setprecision(2) << _statistics->luFactorizationTimeReal << std::endl
       << "Solves             : " << std::setw(10) << _statistics->luSolvesReal << std::endl
       << "  Time spent       : " << std::setw(10) << _statistics->luSolveTimeReal << std::endl
       << "Solution time      : " << std::setw(10) << std::fixed << std::setprecision(2) << solveTime() << std::endl
       << "Iterations         : " << std::setw(10) << numIterations() << std::endl;

    return s.str();
  }



  /// name of starter
  template <>
  const char* SoPlex<Real>::getStarterName()
  {
    if( _starter )
      return _starter->getName();
    else
      return "none";
  }



  /// name of simplifier
  template <>
  const char* SoPlex<Real>::getSimplifierName()
  {
    if( _simplifier )
      return _simplifier->getName();
    else
      return "none";
  }



  /// name of scaling method after simplifier
  template <>
  const char* SoPlex<Real>::getScalerName()
  {
    if( _scaler )
      return _scaler->getName();
    else
      return "none";
  }



  /// name of currently loaded pricer
  template <>
  const char* SoPlex<Real>::getPricerName()
  {
    return _solver.pricer()->getName();
  }



  /// name of currently loaded ratiotester
  template <>
  const char* SoPlex<Real>::getRatiotesterName()
  {
    return _solver.ratiotester()->getName();
  }



  /// reads LP file in LP or MPS format according to READMODE parameter; gets row names, column names, and
  /// integer variables if desired; returns true on success

  // template <>
	// bool SoPlex<Real>::readFile(const char* filename, NameSet* rowNames, NameSet* colNames, DIdxSet* intVars)
  // {
  //   bool success = false;
  //   if( intParam(SoPlex<Real>::READMODE) == READMODE_REAL )
  //     success = _readFileReal(filename, rowNames, colNames, intVars);
  //   else
  //     success = _readFileRational(filename, rowNames, colNames, intVars);

  //   // storing the row and column names for use in the DBDS print basis methods
  //   _rowNames = rowNames;
  //   _colNames = colNames;

  //   return success;
  // }

  template <>
  bool SoPlex<Real>::readFile(const char* filename, NameSet* rowNames, NameSet* colNames, DIdxSet* intVars)
  {
    bool success = false;
    if( intParam(SoPlex::READMODE) == READMODE_REAL )
      success = _readFileReal(filename, rowNames, colNames, intVars);
    else
      success = _readFileRational(filename, rowNames, colNames, intVars);

    // storing the row and column names for use in the DBDS print basis methods
    _rowNames = rowNames;
    _colNames = colNames;

    return success;
  }

  // // The proper templated version of the above function
  // template <>
	// bool SoPlex<Real>::readFile(const char* filename, NameSet* rowNames, NameSet* colNames, DIdxSet* intVars)
  // {
  //   bool success = false;

  //   success = _readFileReal(filename, rowNames, colNames, intVars);

  //   // storing the row and column names for use in the DBDS print basis methods
  //   _rowNames = rowNames;
  //   _colNames = colNames;

  //   return success;
  // }

  // template <>
	// bool SoPlex<Rational>::readFile(const char* filename, NameSet* rowNames, NameSet* colNames, DIdxSet* intVars)
  // {
  //   bool success = false;
  //   success = _readFileRational(filename, rowNames, colNames, intVars);

  //   // storing the row and column names for use in the DBDS print basis methods
  //   _rowNames = rowNames;
  //   _colNames = colNames;

  //   return success;
  // }


  /// writes real LP to file; LP or MPS format is chosen from the extension in \p filename; if \p rowNames and \p
  /// colNames are \c NULL, default names are used; if \p intVars is not \c NULL, the variables contained in it are
  /// marked as integer; returns true on success
  template <>
	bool SoPlex<Real>::writeFileT(const char* filename, const NameSet* rowNames, const NameSet* colNames, const DIdxSet* intVars, const bool unscale) const
  {
    ///@todo implement return value
    if( unscale && _realLP->isScaled() )
      {
        MSG_INFO3( spxout, spxout << "copy LP to write unscaled original problem" << std::endl; )
          SPxLPReal* origLP;
        origLP = 0;
        spx_alloc(origLP);
        origLP = new (origLP) SPxLPReal(*_realLP);
        origLP->unscaleLP();
        origLP->writeFile(filename, rowNames, colNames, intVars);
        origLP->~SPxLPReal();
        spx_free(origLP);
      }
    else
      _realLP->writeFile(filename, rowNames, colNames, intVars);

    return true;
  }



  /// writes rational LP to file; LP or MPS format is chosen from the extension in \p filename; if \p rowNames and \p
  /// colNames are \c NULL, default names are used; if \p intVars is not \c NULL, the variables contained in it are
  /// marked as integer; returns true on success
  /// Here unscale is just a junk variable that is used to match the type with the real write function
  template <>
	bool SoPlex<Rational>::writeFileT(const char* filename, const NameSet* rowNames, const NameSet* colNames, const DIdxSet* intVars, const bool unscale) const
  {
    if( intParam(SoPlex<Rational>::SYNCMODE) == SYNCMODE_ONLYREAL )
      return false;
    else
      {
        assert(_rationalLP != 0);
        _rationalLP->writeFile(filename, rowNames, colNames, intVars);

        ///@todo implement return value
        return true;
      }
  }



  /// writes the dual of the real LP to file; LP or MPS format is chosen from the extension in \p filename;
  /// if \p rowNames and \p colNames are \c NULL, default names are used; if \p intVars is not \c NULL,
  /// the variables contained in it are marked as integer; returns true on success
  template <>
	bool SoPlex<Real>::writeDualFileReal(const char* filename, const NameSet* rowNames, const NameSet* colNames, const DIdxSet* intVars) const
  {
    SPxLPReal dualLP;
    _realLP->buildDualProblem(dualLP);
    dualLP.setOutstream(spxout);

    // swap colnames and rownames
    dualLP.writeFile(filename, colNames, rowNames);
    return true;
  }



  /// reads basis information from \p filename and returns true on success; if \p rowNames and \p colNames are \c NULL,
  /// default names are assumed; returns true on success
  /// #template #baseclass
  template <>
  bool SoPlex<Real>::readBasisFile(const char* filename, const NameSet* rowNames, const NameSet* colNames)
  {
    clearBasis();
    /// @todo can't we just remove the else code?
#if 1
    assert(filename != 0);
    assert(_realLP != 0);

    // start timing
    _statistics->readingTime->start();

    // read
    if( !_isRealLPLoaded )
      {
        assert(_realLP != &_solver);

        _solver.loadLP(*_realLP);
        _realLP->~SPxLPReal();
        spx_free(_realLP);
        _realLP = &_solver;
        _isRealLPLoaded = true;
      }
    _hasBasis = _solver.readBasisFile(filename, rowNames, colNames);
    assert(_hasBasis == (_solver.basis().status() > SPxBasis<Real>::NO_PROBLEM));

    // stop timing
    _statistics->readingTime->stop();

    return _hasBasis;
#else
    // this is alternative code for reading bases without the SPxSolver class
    assert(filename != 0);

    // start timing
    _statistics->readingTime->start();

    // read
    spxifstream file(filename);

    if( !file )
      return false;

    // get problem size
    int numRows = numRowsT();
    int numCols = numColsT();

    // prepare column names
    const NameSet* colNamesPtr = colNames;
    NameSet* tmpColNames = 0;
    if( colNames == 0 )
      {
        std::stringstream name;

        spx_alloc(tmpColNames);
        tmpColNames = new (tmpColNames) NameSet();
        tmpColNames->reMax(numCols);

        for( int j = 0; j < numCols; ++j )
          {
            name << "x" << j;
            tmpColNames->add(name.str().c_str());
          }

        colNamesPtr = tmpColNames;
      }

    // prepare row names
    const NameSet* rowNamesPtr = rowNames;
    NameSet* tmpRowNames = 0;
    if( rowNamesPtr == 0 )
      {
        std::stringstream name;

        spx_alloc(tmpRowNames);
        tmpRowNames = new (tmpRowNames) NameSet();
        tmpRowNames->reMax(numRows);

        for( int i = 0; i < numRows; ++i )
          {
            name << "C" << i;
            tmpRowNames->add(name.str().c_str());
          }

        rowNamesPtr = tmpRowNames;
      }

    // initialize with default slack basis
    _basisStatusRows.reSize(numRows);
    _basisStatusCols.reSize(numCols);

    for( int i = 0; i < numRows; i++ )
      _basisStatusRows[i] = SPxSolver<Real>::BASIC;

    for( int i = 0; i < numCols; i++ )
      {
        if( lowerRealInternal(i) == upperRealInternal(i) )
          _basisStatusCols[i] = SPxSolver<Real>::FIXED;
        else if( lowerRealInternal(i) <= double(-realParam(SoPlex<Real>::INFTY)) && upperRealInternal(i) >= double(realParam(SoPlex<Real>::INFTY)) )
          _basisStatusCols[i] = SPxSolver<Real>::ZERO;
        else if( lowerRealInternal(i) <= double(-realParam(SoPlex<Real>::INFTY)) )
          _basisStatusCols[i] = SPxSolver<Real>::ON_UPPER;
        else
          _basisStatusCols[i] = SPxSolver<Real>::ON_LOWER;
      }

    // read basis
    MPSInput mps(file);
    if( mps.readLine() && (mps.field0() != 0) && !strcmp(mps.field0(), "NAME") )
      {
        while( mps.readLine() )
          {
            int c = -1;
            int r = -1;

            if( mps.field0() != 0 && !strcmp(mps.field0(), "ENDATA") )
              {
                mps.setSection(MPSInput::ENDATA);
                break;
              }

            if( mps.field1() == 0 || mps.field2() == 0 )
              break;

            if( (c = colNamesPtr->number(mps.field2())) < 0 )
              break;

            if( *mps.field1() == 'X' )
              {
                if( mps.field3() == 0 || (r = rowNamesPtr->number(mps.field3())) < 0 )
                  break;
              }

            if( !strcmp(mps.field1(), "XU") )
              {
                _basisStatusCols[c] = SPxSolver<Real>::BASIC;
                if( _rowTypes[r] == SoPlex<Real>::RANGETYPE_LOWER )
                  _basisStatusRows[r] = SPxSolver<Real>::ON_LOWER;
                else if( _rowTypes[r] == SoPlex<Real>::RANGETYPE_FIXED )
                  _basisStatusRows[r] = SPxSolver<Real>::FIXED;
                else
                  _basisStatusRows[r] = SPxSolver<Real>::ON_UPPER;
              }
            else if( !strcmp(mps.field1(), "XL") )
              {
                _basisStatusCols[c] = SPxSolver<Real>::BASIC;
                if( _rowTypes[r] == SoPlex<Real>::RANGETYPE_UPPER )
                  _basisStatusRows[r] = SPxSolver<Real>::ON_UPPER;
                else if( _rowTypes[r] == SoPlex<Real>::RANGETYPE_FIXED )
                  _basisStatusRows[r] = SPxSolver<Real>::FIXED;
                else
                  _basisStatusRows[r] = SPxSolver<Real>::ON_LOWER;
              }
            else if( !strcmp(mps.field1(), "UL") )
              {
                _basisStatusCols[c] = SPxSolver<Real>::ON_UPPER;
              }
            else if( !strcmp(mps.field1(), "LL") )
              {
                _basisStatusCols[c] = SPxSolver<Real>::ON_LOWER;
              }
            else
              {
                mps.syntaxError();
                break;
              }
          }
      }

    if( rowNames == 0 )
      {
        tmpRowNames->~NameSet();
        spx_free(tmpRowNames);
      }

    if( colNames == 0 )
      {
        tmpColNames->~NameSet();
        spx_free(tmpColNames);
      }

    _hasBasis = !mps.hasError();

    // stop timing
    _statistics->readingTime->stop();

    return _hasBasis;
#endif
  }



  /// writes basis information to \p filename; if \p rowNames and \p colNames are \c NULL, default names are used;
  /// returns true on success
  /// #baseclass
  template <>
	bool SoPlex<Real>::writeBasisFile(const char* filename, const NameSet* rowNames, const NameSet* colNames, const bool cpxFormat) const
  {
    assert(filename != 0);

    if( _isRealLPLoaded )
      return _solver.writeBasisFile(filename, rowNames, colNames, cpxFormat);
    else
      {
        std::ofstream file(filename);
        if( !file.good() )
          return false;

        file.setf(std::ios::left);
        file << "NAME  " << filename << "\n";

        // do not write basis if there is none
        if( !_hasBasis )
          {
            file << "ENDATA\n";
            return true;
          }

        // start writing
        int numRows = _basisStatusRows.size();
        int numCols = _basisStatusCols.size();
        int row = 0;

        for( int col = 0; col < numCols; col++ )
          {
            assert(_basisStatusCols[col] != SPxSolver<Real>::UNDEFINED);

            if( _basisStatusCols[col] == SPxSolver<Real>::BASIC )
              {
                // find nonbasic row
                for( ; row < numRows; row++ )
                  {
                    assert(_basisStatusRows[row] != SPxSolver<Real>::UNDEFINED);
                    if( _basisStatusRows[row] != SPxSolver<Real>::BASIC )
                      break;
                  }

                assert(row != numRows);

                if( _basisStatusRows[row] == SPxSolver<Real>::ON_UPPER && (!cpxFormat || _rowTypes[row] == SoPlex<Real>::RANGETYPE_BOXED) )
                  file << " XU ";
                else
                  file << " XL ";

                file << std::setw(8);
                if( colNames != 0 && colNames->has(col) )
                  file << (*colNames)[col];
                else
                  file << "x" << col;

                file << "       ";
                if( rowNames != 0 && rowNames->has(row) )
                  file << (*rowNames)[row];
                else
                  file << "C" << row;

                file << "\n";
                row++;
              }
            else
              {
                if( _basisStatusCols[col] == SPxSolver<Real>::ON_UPPER )
                  {
                    file << " UL ";

                    file << std::setw(8);
                    if( colNames != 0 && colNames->has(col) )
                      file << (*colNames)[col];
                    else
                      file << "x" << col;

                    file << "\n";
                  }
              }
          }

        file << "ENDATA\n";

#ifndef NDEBUG
        // check that the remaining rows are basic
        for( ; row < numRows; row++ )
          {
            assert(_basisStatusRows[row] == SPxSolver<Real>::BASIC);
          }
#endif

        return true;
      }
  }



  /// writes internal LP, basis information, and parameter settings; if \p rowNames and \p colNames are \c NULL,
  /// default names are used
  template <>
  void SoPlex<Real>::writeStateReal(const char* filename, const NameSet* rowNames, const NameSet* colNames, const bool cpxFormat) const
  {
    std::string ofname;

    // write parameter settings
    ofname = std::string(filename) + ".set";
    saveSettingsFile(ofname.c_str());

    // write problem in MPS/LP format
    ofname = std::string(filename) + ((cpxFormat) ? ".lp" : ".mps");
    writeFileT(ofname.c_str(), rowNames, colNames, 0);

    // write basis
    ofname = std::string(filename) + ".bas";
    writeBasisFile(ofname.c_str(), rowNames, colNames, cpxFormat);
  }



  /// writes internal LP, basis information, and parameter settings; if \p rowNames and \p colNames are \c NULL,
  /// default names are used
  template <>
  void SoPlex<Real>::writeStateRational(const char* filename, const NameSet* rowNames, const NameSet* colNames, const bool cpxFormat) const
  {
    std::string ofname;

    // write parameter settings
    ofname = std::string(filename) + ".set";
    saveSettingsFile(ofname.c_str());

    // write problem in MPS/LP format
    ofname = std::string(filename) + ((cpxFormat) ? ".lp" : ".mps");
    writeFileT(ofname.c_str(), rowNames, colNames, 0);

    // write basis
    ofname = std::string(filename) + ".bas";
    writeBasisFile(ofname.c_str(), rowNames, colNames, cpxFormat);
  }

  // /// returns integer parameter value
  // // #template #baseclass
  // template <>
	// int SoPlex<Real>::intParam(const IntParam param) const
  // {
  //   assert(param >= 0);
  //   assert(param < INTPARAM_COUNT);
  //   return _currentSettings->_intParamValues[param];
  // }

#ifdef SOPLEX_WITH_RATIONALPARAM
  /// returns rational parameter value
  Rational SoPlex<Real>::rationalParam(const RationalParam param) const
  {
    assert(param >= 0);
    assert(param < RATIONALPARAM_COUNT);
    return _currentSettings->_rationalParamValues[param];
  }
#endif



  /// returns current parameter settings
  template <>
  const typename SoPlex<Real>::Settings& SoPlex<Real>::settings() const
  {
    return *_currentSettings;
  }



  /// sets boolean parameter value; returns true on success
  template <>
	bool SoPlex<Real>::setBoolParam(const BoolParam param, const bool value, const bool init)
  {
    assert(param >= 0);
    assert(param < SoPlex<Real>::BOOLPARAM_COUNT);
    assert(init || _isConsistent());

    if( !init && value == boolParam(param) )
      return true;

    switch( param )
      {
      case LIFTING:
        break;
      case EQTRANS:
        break;
      case TESTDUALINF:
        break;
      case RATFAC:
        break;
      case USEDECOMPDUALSIMPLEX:
        break;
      case COMPUTEDEGEN:
        break;
      case USECOMPDUAL:
        break;
      case EXPLICITVIOL:
        break;
      case ACCEPTCYCLING:
        break;
      case RATREC:
        break;
      case POWERSCALING:
        break;
      case RATFACJUMP:
        break;
      case ROWBOUNDFLIPS:
        _ratiotesterBoundFlipping.useBoundFlipsRow(value);
        break;
      case PERSISTENTSCALING:
        break;
      case FULLPERTURBATION:
        _solver.useFullPerturbation(value);
        break;
      default:
        return false;
      }

    _currentSettings->_boolParamValues[param] = value;
    return true;
  }

  /// sets real parameter value; returns true on success
  template <>
	bool SoPlex<Real>::setRealParam(const RealParam param, const Real value, const bool init)
  {
    assert(param >= 0);
    assert(param < REALPARAM_COUNT);
    assert(init || _isConsistent());

    if( !init && value == realParam(param) )
      return true;

    if( value < _currentSettings->realParam.lower[param] || value > _currentSettings->realParam.upper[param] )
      return false;

    // required to set a different feastol or opttol
    Real tmp_value = value;

    switch( param )
      {
        // primal feasibility tolerance; passed to the floating point solver only when calling solve()
      case SoPlex<Real>::FEASTOL:
#ifndef SOPLEX_WITH_GMP
        if( value < DEFAULT_EPS_PIVOT )
          {
            MSG_WARNING( spxout, spxout << "Cannot set feasibility tolerance to small value " << value << " without GMP - using " << DEFAULT_EPS_PIVOT << ".\n");
            tmp_value = DEFAULT_EPS_PIVOT;
            _rationalFeastol = DEFAULT_EPS_PIVOT;
            break;
          }
#endif
        _rationalFeastol = value;
        break;

        // dual feasibility tolerance; passed to the floating point solver only when calling solve()
      case SoPlex<Real>::OPTTOL:
#ifndef SOPLEX_WITH_GMP
        if( value < DEFAULT_EPS_PIVOT )
          {
            MSG_WARNING( spxout, spxout << "Cannot set optimality tolerance to small value " << value << " without GMP - using " << DEFAULT_EPS_PIVOT << ".\n");
            tmp_value = DEFAULT_EPS_PIVOT;
            _rationalOpttol = DEFAULT_EPS_PIVOT;
            break;
          }
#endif
        _rationalOpttol = value;
        break;

        // general zero tolerance
      case SoPlex<Real>::EPSILON_ZERO:
        Param::setEpsilon(value);
        break;

        // zero tolerance used in factorization
      case SoPlex<Real>::EPSILON_FACTORIZATION:
        Param::setEpsilonFactorization(value);
        break;

        // zero tolerance used in update of the factorization
      case SoPlex<Real>::EPSILON_UPDATE:
        Param::setEpsilonUpdate(value);
        break;

        // pivot zero tolerance used in factorization (declare numerical singularity for small LU pivots)
      case SoPlex<Real>::EPSILON_PIVOT:
        Param::setEpsilonPivot(value);
        break;

        // infinity threshold
      case SoPlex<Real>::INFTY:
        _rationalPosInfty = value;
        _rationalNegInfty = -value;
        if( intParam(SoPlex<Real>::SYNCMODE) != SYNCMODE_ONLYREAL )
          _recomputeRangeTypesRational();
        break;

        // time limit in seconds (INFTY if unlimited)
      case SoPlex<Real>::TIMELIMIT:
        break;

        // lower limit on objective value is set in solveReal()
      case SoPlex<Real>::OBJLIMIT_LOWER:
        break;

        // upper limit on objective value is set in solveReal()
      case SoPlex<Real>::OBJLIMIT_UPPER:
        break;

        // working tolerance for feasibility in floating-point solver
      case SoPlex<Real>::FPFEASTOL:
        break;

        // working tolerance for optimality in floating-point solver
      case SoPlex<Real>::FPOPTTOL:
        break;

        // maximum increase of scaling factors between refinements
      case SoPlex<Real>::MAXSCALEINCR:
        _rationalMaxscaleincr = value;
        break;

        // lower threshold in lifting (nonzero matrix coefficients with smaller absolute value will be reformulated)
      case SoPlex<Real>::LIFTMINVAL:
        break;

        // upper threshold in lifting (nonzero matrix coefficients with larger absolute value will be reformulated)
      case SoPlex<Real>::LIFTMAXVAL:
        break;

        // threshold for sparse pricing
      case SoPlex<Real>::SPARSITY_THRESHOLD:
        break;

        // threshold on number of rows vs. number of columns for switching from column to row representations in auto mode
      case SoPlex<Real>::REPRESENTATION_SWITCH:
        break;

        // geometric frequency at which to apply rational reconstruction
      case SoPlex<Real>::RATREC_FREQ:
        break;

        // minimal reduction (sum of removed rows/cols) to continue simplification
      case SoPlex<Real>::MINRED:
        break;

      case SoPlex<Real>::REFAC_BASIS_NNZ:
        break;

      case SoPlex<Real>::REFAC_UPDATE_FILL:
        break;

      case SoPlex<Real>::REFAC_MEM_FACTOR:
        break;

        // accuracy of conjugate gradient method in least squares scaling (higher value leads to more iterations)
      case SoPlex<Real>::LEASTSQ_ACRCY:
        if( _scaler )
          _scaler->setRealParam(value);
        break;

        // objective offset
      case SoPlex<Real>::OBJ_OFFSET:
        if( _realLP )
          _realLP->changeObjOffset(value);
        if( _rationalLP )
          _rationalLP->changeObjOffset(value);
        break;

      default:
        return false;
      }

    _currentSettings->_realParamValues[param] = tmp_value;
    return true;
  }

  /// #template #temp

#ifdef SOPLEX_WITH_RATIONALPARAM
  /// sets rational parameter value; returns true on success
  template <>
	bool SoPlex<Real>::setRationalParam(const RationalParam param, const Rational value, const bool init)
  {
    assert(param >= 0);
    assert(param < RATIONALPARAM_COUNT);
    assert(init || _isConsistent());

    if( !init && value == rationalParam(param) )
      return true;

    if( value < _currentSettings->rationalParam.lower[param] || value > _currentSettings->rationalParam.upper[param] )
      return false;

    switch( param )
      {
      default:
        // currently, there are no rational-valued parameters
        return false;
      }

    _currentSettings->_rationalParamValues[param] = value;
    return true;
  }
#endif



  /// sets parameter settings; returns true on success
  template <>
	bool SoPlex<Real>::setSettings(const Settings& newSettings, const bool init)
  {
    assert(init || _isConsistent());

    bool success = true;

    *_currentSettings = newSettings;

    for( int i = 0; i < SoPlex<Real>::BOOLPARAM_COUNT; i++ )
      success &= setBoolParam((BoolParam)i, _currentSettings->_boolParamValues[i], init);

    for( int i = 0; i < SoPlex<Real>::INTPARAM_COUNT; i++ )
      success &= setIntParam((IntParam)i, _currentSettings->_intParamValues[i], init);

    for( int i = 0; i < SoPlex<Real>::REALPARAM_COUNT; i++ )
      success &= setRealParam((RealParam)i, _currentSettings->_realParamValues[i], init);

#ifdef SOPLEX_WITH_RATIONALPARAM
    for( int i = 0; i < SoPlex<Real>::RATIONALPARAM_COUNT; i++ )
      success &= setRationalParam((RationalParam)i, _currentSettings->_rationalParamValues[i], init);
#endif

    assert(_isConsistent());

    return success;
  }

  /// resets default parameter settings
  template <>
  void SoPlex<Real>::resetSettings(const bool quiet, const bool init)
  {
    for( int i = 0; i < SoPlex<Real>::BOOLPARAM_COUNT; i++ )
      setBoolParam((BoolParam)i, _currentSettings->boolParam.defaultValue[i], init);

    for( int i = 0; i < SoPlex<Real>::INTPARAM_COUNT; i++ )
      setIntParam((IntParam)i, _currentSettings->intParam.defaultValue[i], init);

    for( int i = 0; i < SoPlex<Real>::REALPARAM_COUNT; i++ )
      setRealParam((RealParam)i, _currentSettings->realParam.defaultValue[i], init);

#ifdef SOPLEX_WITH_RATIONALPARAM
    for( int i = 0; i < SoPlex<Real>::RATIONALPARAM_COUNT; i++ )
      success &= setRationalParam((RationalParam)i, _currentSettings->rationalParam.defaultValue[i], init);
#endif
  }


  /// print non-default parameter values
  /// #template #baseclass
  template <>
  void SoPlex<Real>::printUserSettings()
  {
    bool printedValue = false;

    SPxOut::setFixed(spxout.getCurrentStream());

    for( int i = 0; i < SoPlex<Real>::BOOLPARAM_COUNT; i++ )
      {
        if( _currentSettings->_boolParamValues[i] == _currentSettings->boolParam.defaultValue[i] )
          continue;

        spxout << "bool:" << _currentSettings->boolParam.name[i] << " = " << (_currentSettings->_boolParamValues[i] ? "true\n" : "false\n");
        printedValue = true;
      }

    for( int i = 0; i < SoPlex<Real>::INTPARAM_COUNT; i++ )
      {
        if( _currentSettings->_intParamValues[i] == _currentSettings->intParam.defaultValue[i] )
          continue;

        spxout << "int:" << _currentSettings->intParam.name[i] << " = " << _currentSettings->_intParamValues[i] << "\n";
        printedValue = true;
      }

    SPxOut::setScientific(spxout.getCurrentStream());

    for( int i = 0; i < SoPlex<Real>::REALPARAM_COUNT; i++ )
      {
        if( _currentSettings->_realParamValues[i] == _currentSettings->realParam.defaultValue[i] )
          continue;

        spxout << "real:" << _currentSettings->realParam.name[i] << " = " << _currentSettings->_realParamValues[i] << "\n";
        printedValue = true;
      }

#ifdef SOPLEX_WITH_RATIONALPARAM
    for( int i = 0; i < SoPlex<Real>::RATIONALPARAM_COUNT; i++ )
      {
        if( _currentSettings->_rationalParamValues[i] == _currentSettings->rationalParam.defaultValue[i] )
          continue;

        spxout << "rational:" << _currentSettings->rationalParam.name[i] << " = " << _currentSettings->_rationalParamValues[i] << "\n";
        printedValue = true;
      }
#endif

    if( _solver.random.getSeed() != DEFAULT_RANDOM_SEED )
      {
        spxout << "uint:random_seed = " << _solver.random.getSeed() << "\n";
        printedValue = true;
      }

    if( printedValue )
      spxout << std::endl;
  }



  /// writes settings file; returns true on success
  // #template #baseclass
  template <>
	bool SoPlex<Real>::saveSettingsFile(const char* filename, const bool onlyChanged) const
  {
    assert(filename != 0);

    std::ofstream file(filename);
    SPxOut::setScientific(file, 16);

    if( !file.good() )
      return false;

    file.setf(std::ios::left);

    SPxOut::setFixed(file);

    file << "# SoPlex version " << SOPLEX_VERSION / 100 << "." << (SOPLEX_VERSION / 10) % 10 << "." << SOPLEX_VERSION % 10;
#if SOPLEX_SUBVERSION > 0
    file << "." << SOPLEX_SUBVERSION;
#endif
    file << "\n";

    for( int i = 0; i < SoPlex<Real>::BOOLPARAM_COUNT; i++ )
      {
        if( onlyChanged && _currentSettings->_boolParamValues[i] == _currentSettings->boolParam.defaultValue[i] )
          continue;

        file << "\n";
        file << "# " << _currentSettings->boolParam.description[i] << "\n";
        file << "# range {true, false}, default " << (_currentSettings->boolParam.defaultValue[i] ? "true\n" : "false\n");
        file << "bool:" << _currentSettings->boolParam.name[i] << " = " << (_currentSettings->_boolParamValues[i] ? "true\n" : "false\n");
      }

    for( int i = 0; i < SoPlex<Real>::INTPARAM_COUNT; i++ )
      {
        if( onlyChanged && _currentSettings->_intParamValues[i] == _currentSettings->intParam.defaultValue[i] )
          continue;

        file << "\n";
        file << "# " << _currentSettings->intParam.description[i] << "\n";
        file << "# range [" << _currentSettings->intParam.lower[i] << "," << _currentSettings->intParam.upper[i]
             << "], default " << _currentSettings->intParam.defaultValue[i] << "\n";
        file << "int:" << _currentSettings->intParam.name[i] << " = " << _currentSettings->_intParamValues[i] << "\n";
      }

    SPxOut::setScientific(file);

    for( int i = 0; i < SoPlex<Real>::REALPARAM_COUNT; i++ )
      {
        if( onlyChanged && _currentSettings->_realParamValues[i] == _currentSettings->realParam.defaultValue[i] )
          continue;

        file << "\n";
        file << "# " << _currentSettings->realParam.description[i] << "\n";
        file << "# range [" << _currentSettings->realParam.lower[i] << "," << _currentSettings->realParam.upper[i]
             << "], default " << _currentSettings->realParam.defaultValue[i] << "\n";
        file << "real:" << _currentSettings->realParam.name[i] << " = " << _currentSettings->_realParamValues[i] << "\n";
      }

#ifdef SOPLEX_WITH_RATIONALPARAM
    for( int i = 0; i < SoPlex<Real>::RATIONALPARAM_COUNT; i++ )
      {
        if( onlyChanged && _currentSettings->_rationalParamValues[i] == _currentSettings->rationalParam.defaultValue[i] )
          continue;

        file << "\n";
        file << "# " << _currentSettings->rationalParam.description[i] << "\n";
        file << "# range [" << _currentSettings->rationalParam.lower[i] << "," << _currentSettings->rationalParam.upper[i]
             << "], default " << _currentSettings->rationalParam.defaultValue[i] << "\n";
        file << "rational:" << _currentSettings->rationalParam.name[i] << " = " << _currentSettings->_rationalParamValues[i] << "\n";
      }
#endif

    if( !onlyChanged || _solver.random.getSeed() != DEFAULT_RANDOM_SEED )
      {
        file << "\n";
        file << "# initial random seed used for perturbation\n";
        file << "# range [0, " << UINT_MAX << "], default "<< DEFAULT_RANDOM_SEED << "\n";
        file << "uint:random_seed = " << _solver.random.getSeed() << "\n";
      }

    return true;
  }



  /// reads settings file; returns true on success
  /// #template #baseclass
  template <>
	bool SoPlex<Real>::loadSettingsFile(const char* filename)
  {
    assert(filename != 0);

    // start timing
    _statistics->readingTime->start();

    MSG_INFO1( spxout, spxout << "Loading settings file <" << filename << "> . . .\n" );

    // open file
    spxifstream file(filename);

    if( !file )
      {
        MSG_INFO1( spxout, spxout << "Error opening settings file.\n" );
        return false;
      }

    // read file
    char line[SET_MAX_LINE_LEN];
    int lineNumber = 0;
    bool readError = false;
    bool parseError = false;

    while( !readError && !parseError)
      {
        lineNumber++;
        readError = !file.getline(line, sizeof(line));
        if( !readError )
          parseError = !_parseSettingsLine(line, lineNumber);
      }
    readError = readError && !file.eof();

    if( readError && strlen(line) == SET_MAX_LINE_LEN - 1 )
      {
        MSG_INFO1( spxout, spxout << "Error reading settings file: line " << lineNumber << " in settings file exceeds " << SET_MAX_LINE_LEN - 2 << " characters.\n" );
      }
    else if( readError )
      {
        MSG_INFO1( spxout, spxout << "Error reading settings file: line " << lineNumber << ".\n" );
      }

    // stop timing
    _statistics->readingTime->stop();

    return !readError;
  }

  /// parses one setting string and returns true on success
  /// #template #temp #baseclass
  template <>
	bool SoPlex<Real>::parseSettingsString(char* string)
  {
    assert(string != 0);
    if( string == 0 )
      return false;

    char parseString[SET_MAX_LINE_LEN];
    strncpy(parseString, string, SET_MAX_LINE_LEN-1);
    parseString[SET_MAX_LINE_LEN-1] = '\0';

    char* line = parseString;

    // find the start of the parameter type
    while( *line == ' ' || *line == '\t' || *line == '\r' )
      line++;
    if( *line == '\0' || *line == '\n' || *line == '#' )
      return true;
    char* paramTypeString = line;

    // find the end of the parameter type
    while( *line != ' ' && *line != '\t' && *line != '\r' && *line != '\n' && *line != '#' && *line != '\0' && *line != ':' )
      line++;
    if( *line == ':' )
      {
        *line = '\0';
        line++;
      }
    else
      {
        *line = '\0';
        line++;

        // search for the ':' char in the line
        while( *line == ' ' || *line == '\t' || *line == '\r' )
          line++;
        if( *line != ':' )
          {
            MSG_INFO1( spxout, spxout << "Error parsing setting string: no ':' separating parameter type and name.\n" );
            return false;
          }
        line++;
      }

    // find the start of the parameter name
    while( *line == ' ' || *line == '\t' || *line == '\r' )
      line++;
    if( *line == '\0' || *line == '\n' || *line == '#' )
      {
        MSG_INFO1( spxout, spxout << "Error parsing setting string: no parameter name.\n");
        return false;
      }
    char* paramName = line;

    // find the end of the parameter name
    while( *line != ' ' && *line != '\t' && *line != '\r' && *line != '\n' && *line != '#' && *line != '\0' && *line != '=' )
      line++;
    if( *line == '=' )
      {
        *line = '\0';
        line++;
      }
    else
      {
        *line = '\0';
        line++;

        // search for the '=' char in the line
        while( *line == ' ' || *line == '\t' || *line == '\r' )
          line++;
        if( *line != '=' )
          {
            MSG_INFO1( spxout, spxout << "Error parsing setting string: no '=' after parameter name.\n" );
            return false;
          }
        line++;
      }

    // find the start of the parameter value string
    while( *line == ' ' || *line == '\t' || *line == '\r' )
      line++;
    if( *line == '\0' || *line == '\n' || *line == '#' )
      {
        MSG_INFO1( spxout, spxout << "Error parsing setting string: no parameter value.\n");
        return false;
      }
    char* paramValueString = line;

    // find the end of the parameter value string
    while( *line != ' ' && *line != '\t' && *line != '\r' && *line != '\n' && *line != '#' && *line != '\0' )
      line++;
    if( *line != '\0' )
      {
        // check, if the rest of the line is clean
        *line = '\0';
        line++;
        while( *line == ' ' || *line == '\t' || *line == '\r' )
          line++;
        if( *line != '\0' && *line != '\n' && *line != '#' )
          {
            MSG_INFO1( spxout, spxout << "Error parsing setting string: additional character '" << *line << "' after parameter value.\n" );
            return false;
          }
      }

    // check whether we have a bool parameter
    if( strncmp(paramTypeString, "bool", 4) == 0 )
      {
        for( int param = 0; ; param++ )
          {
            if( param >= BOOLPARAM_COUNT )
              {
                MSG_INFO1( spxout, spxout << "Error parsing setting string: unknown parameter name <" << paramName << ">.\n" );
                return false;
              }
            else if( strncmp(paramName, _currentSettings->boolParam.name[param].c_str(), SET_MAX_LINE_LEN) == 0 )
              {
                if( strncasecmp(paramValueString, "true", 4) == 0
                    || strncasecmp(paramValueString, "TRUE", 4) == 0
                    || strncasecmp(paramValueString, "t", 4) == 0
                    || strncasecmp(paramValueString, "T", 4) == 0
                    || strtol(paramValueString, NULL, 4) == 1 )
                  {
                    setBoolParam((BoolParam)param, true);
                    break;
                  }
                else if( strncasecmp(paramValueString, "false", 5) == 0
                         || strncasecmp(paramValueString, "FALSE", 5) == 0
                         || strncasecmp(paramValueString, "f", 5) == 0
                         || strncasecmp(paramValueString, "F", 5) == 0
                         || strtol(paramValueString, NULL, 5) == 0 )
                  {
                    setBoolParam((BoolParam)param, false);
                    break;
                  }
                else
                  {
                    MSG_INFO1( spxout, spxout << "Error parsing setting string: invalid value <" << paramValueString << "> for bool parameter <" << paramName << ">.\n" );
                    return false;
                  }
              }
          }

        return true;
      }

    // check whether we have an integer parameter
    if( strncmp(paramTypeString, "int", 3) == 0 )
      {
        for( int param = 0; ; param++ )
          {
            if( param >= INTPARAM_COUNT )
              {
                MSG_INFO1( spxout, spxout << "Error parsing setting string: unknown parameter name <" << paramName << ">.\n" );
                return false;
              }
            else if( strncmp(paramName, _currentSettings->intParam.name[param].c_str(), SET_MAX_LINE_LEN) == 0 )
              {
                int value;

                if( sscanf(paramValueString, "%d", &value) == 1 && setIntParam((IntParam)param, value, false) )
                  break;
                else
                  {
                    MSG_INFO1( spxout, spxout << "Error parsing setting string: invalid value <" << paramValueString << "> for int parameter <" << paramName << ">.\n" );
                    return false;
                  }
              }
          }

        return true;
      }

    // check whether we have a real parameter
    if( strncmp(paramTypeString, "real", 4) == 0 )
      {
        for( int param = 0; ; param++ )
          {
            if( param >= REALPARAM_COUNT )
              {
                MSG_INFO1( spxout, spxout << "Error parsing setting string: unknown parameter name <" << paramName << ">.\n" );
                return false;
              }
            else if( strncmp(paramName, _currentSettings->realParam.name[param].c_str(), SET_MAX_LINE_LEN) == 0 )
              {
                Real value;

                if( sscanf(paramValueString, "%" REAL_FORMAT, &value) == 1 && setRealParam((RealParam)param, value) )
                  break;
                else
                  {
                    MSG_INFO1( spxout, spxout << "Error parsing setting string: invalid value <" << paramValueString << "> for real parameter <" << paramName << ">.\n" );
                    return false;
                  }
              }
          }

        return true;
      }

#ifdef SOPLEX_WITH_RATIONALPARAM
    // check whether we have a rational parameter
    if( strncmp(paramTypeString, "rational", 8) == 0 )
      {
        for( int param = 0; ; param++ )
          {
            if( param >= RATIONALPARAM_COUNT )
              {
                MSG_INFO1( spxout, spxout << "Error parsing setting string: unknown parameter name <" << paramName << ">.\n" );
                return false;
              }
            else if( strncmp(paramName, _currentSettings->rationalParam.name[param].c_str(), SET_MAX_LINE_LEN) == 0 )
              {
                Rational value;

                if( readStringRational(paramValueString, value) && setRationalParam((RationalParam)param, value) )
                  break;
                else
                  {
                    MSG_INFO1( spxout, spxout << "Error parsing setting string: invalid value <" << paramValueString << "> for rational parameter <" << paramName << ">.\n" );
                    return false;
                  }
              }
          }

        return true;
      }
#endif

    // check whether we have the random seed
    if( strncmp(paramTypeString, "uint", 4) == 0 )
      {
        if( strncmp(paramName, "random_seed", 11) == 0 )
          {
            unsigned int value;
            char format[SPX_MAXSTRLEN];
            spxSnprintf(format, sizeof(format), "%%%du", (int) sizeof(value)-1);

            if( sscanf(paramValueString, format, &value) == 1 )
              {

                setRandomSeed(value);
                return true;
              }
          }

        MSG_INFO1( spxout, spxout << "Error parsing setting string for uint parameter <random_seed>.\n" );
        return false;
      }

    MSG_INFO1( spxout, spxout << "Error parsing setting string: invalid parameter type <" << paramTypeString << "> for parameter <" << paramName << ">.\n" );

    return false;
  }




  /// prints solution statistics
  template <>
  void SoPlex<Real>::printSolutionStatistics(std::ostream& os)
  {
    SPxOut::setScientific(os);
    if( _lastSolveMode == SOLVEMODE_REAL )
      {
        os << "Solution (real)     : \n"
           << "  Objective value   : " << objValueReal() << "\n";
      }
    else if( _lastSolveMode == SOLVEMODE_RATIONAL )
      {
        os << "Solution (rational) : \n"
           << "  Objective value   : " << rationalToString(objValueRational()) << "\n";
        os << "Size (base 2/10)    : \n"
           << "  Total primal      : " << totalSizePrimalRational() << " / " << totalSizePrimalRational(10) << "\n"
           << "  Total dual        : " << totalSizeDualRational() << " / " << totalSizeDualRational(10) << "\n"
           << "  DLCM primal       : " << dlcmSizePrimalRational() << " / " << dlcmSizePrimalRational(10) << "\n"
           << "  DLCM dual         : " << dlcmSizeDualRational() << " / " << dlcmSizeDualRational(10) << "\n"
           << "  DMAX primal       : " << dmaxSizePrimalRational() << " / " << dmaxSizePrimalRational(10) << "\n"
           << "  DMAX dual         : " << dmaxSizeDualRational() << " / " << dmaxSizeDualRational(10) << "\n";
      }
    else
      {
        os << "Solution            : \n"
           << "  Objective value   : -\n";
      }

    if( intParam(CHECKMODE) == CHECKMODE_RATIONAL
        || (intParam(CHECKMODE) == CHECKMODE_AUTO && intParam(READMODE) == READMODE_RATIONAL) )
      {
        Rational maxviol;
        Rational sumviol;

        /// TODO: The following parts should be uncommented when the Rational parts needs to get implemented

        // os << "Violation (rational): \n";
        // if( getBoundViolationT(maxviol, sumviol) )
        //   os << "  Max/sum bound     : " << rationalToString(maxviol) << " / " << rationalToString(sumviol) << "\n";
        // else
        //   os << "  Max/sum bound     : - / -\n";
        // if( getRowViolationT(maxviol, sumviol) )
        //   os << "  Max/sum row       : " << rationalToString(maxviol) << " / " << rationalToString(sumviol) << "\n";
        // else
        //   os << "  Max/sum row       : - / -\n";
        // if( getRedCostViolationT(maxviol, sumviol) )
        //   os << "  Max/sum redcost   : " << rationalToString(maxviol) << " / " << rationalToString(sumviol) << "\n";
        // else
        //   os << "  Max/sum redcost   : - / -\n";
        // if( getDualViolationT(maxviol, sumviol) )
        //   os << "  Max/sum dual      : " << rationalToString(maxviol) << " / " << rationalToString(sumviol) << "\n";
        // else
        //   os << "  Max/sum dual      : - / -\n";
      }
    else
      {
        Real maxviol;
        Real sumviol;

        os << "Violations (real)   : \n";
        if( getBoundViolationT(maxviol, sumviol) )
          os << "  Max/sum bound     : " << maxviol << " / " << sumviol << "\n";
        else
          os << "  Max/sum bound     : - / -\n";
        if( getRowViolationT(maxviol, sumviol) )
          os << "  Max/sum row       : " << maxviol << " / " << sumviol << "\n";
        else
          os << "  Max/sum row       : - / -\n";
        if( getRedCostViolationT(maxviol, sumviol) )
          os << "  Max/sum redcost   : " << maxviol << " / " << sumviol << "\n";
        else
          os << "  Max/sum redcost   : - / -\n";
        if( getDualViolationT(maxviol, sumviol) )
          os << "  Max/sum dual      : " << maxviol << " / " << sumviol << "\n";
        else
          os << "  Max/sum dual      : - / -\n";
      }
  }


  /// prints statistics on solving process
  template <>
  void SoPlex<Real>::printSolvingStatistics(std::ostream& os)
  {
    assert(_statistics != 0);
    _statistics->print(os);
  }



  /// prints short statistics
  template <>
  void SoPlex<Real>::printShortStatistics(std::ostream& os)
  {
    printStatus(os, _status);
    SPxOut::setFixed(os, 2);
    os << "Solving time (sec)  : " << _statistics->solvingTime->time() << "\n"
       << "Iterations          : " << _statistics->iterations << "\n";
    SPxOut::setScientific(os);
    os << "Objective value     : " << objValueReal() << "\n";
  }



  /// prints complete statistics
  /// #template #baseclass
  template <>
  void SoPlex<Real>::printStatistics(std::ostream& os)
  {
    SPxOut::setFixed(os, 2);

    printStatus(os, _status);

    os << "Original problem    : \n";
    if ( boolParam(SoPlex<Real>::USEDECOMPDUALSIMPLEX) )
      printOriginalProblemStatistics(os);
    else
      {
        if( intParam(SoPlex<Real>::READMODE) == READMODE_REAL )
          _realLP->printProblemStatistics(os);
        else
          _rationalLP->printProblemStatistics(os);
      }

    os << "Objective sense     : " << (intParam(SoPlex<Real>::OBJSENSE) == SoPlex<Real>::OBJSENSE_MINIMIZE ? "minimize\n" : "maximize\n");
    printSolutionStatistics(os);
    printSolvingStatistics(os);
  }

  /// prints status
  template <>
  void SoPlex<Real>::printStatus(std::ostream& os, typename SPxSolver<Real>::Status stat)
  {
    os << "SoPlex status       : ";

    switch( stat )
      {
      case SPxSolver<Real>::ERROR:
        os << "error [unspecified]";
        break;
      case SPxSolver<Real>::NO_RATIOTESTER:
        os << "error [no ratiotester loaded]";
        break;
      case SPxSolver<Real>::NO_PRICER:
        os << "error [no pricer loaded]";
        break;
      case SPxSolver<Real>::NO_SOLVER:
        os << "error [no linear solver loaded]";
        break;
      case SPxSolver<Real>::NOT_INIT:
        os << "error [not initialized]";
        break;
      case SPxSolver<Real>::ABORT_CYCLING:
        os << "solving aborted [cycling]";
        break;
      case SPxSolver<Real>::ABORT_TIME:
        os << "solving aborted [time limit reached]";
        break;
      case SPxSolver<Real>::ABORT_ITER:
        os << "solving aborted [iteration limit reached]";
        break;
      case SPxSolver<Real>::ABORT_VALUE:
        os << "solving aborted [objective limit reached]";
        break;
      case SPxSolver<Real>::NO_PROBLEM:
        os << "no problem loaded";
        break;
      case SPxSolver<Real>::REGULAR:
        os << "basis is regular";
        break;
      case SPxSolver<Real>::SINGULAR:
        os << "basis is singular";
        break;
      case SPxSolver<Real>::OPTIMAL:
        os << "problem is solved [optimal]";
        break;
      case SPxSolver<Real>::UNBOUNDED:
        os << "problem is solved [unbounded]";
        break;
      case SPxSolver<Real>::INFEASIBLE:
        os << "problem is solved [infeasible]";
        break;
      case SPxSolver<Real>::INForUNBD:
        os << "problem is solved [infeasible or unbounded]";
        break;
      default:
      case SPxSolver<Real>::UNKNOWN:
        os << "unknown";
        break;
      }

    os << "\n";
  }



  /// prints version and compilation options
  template <>
  void SoPlex<Real>::printVersion() const
  {
    // do not use preprocessor directives within the MSG_INFO1 macro
#if (SOPLEX_SUBVERSION > 0)
    MSG_INFO1( spxout, spxout << "SoPlex version " << SOPLEX_VERSION/100
               << "." << (SOPLEX_VERSION % 100)/10
               << "." << SOPLEX_VERSION % 10
               << "." << SOPLEX_SUBVERSION );
#else
    MSG_INFO1( spxout, spxout << "SoPlex version " << SOPLEX_VERSION/100
               << "." << (SOPLEX_VERSION % 100)/10
               << "." << SOPLEX_VERSION % 10 );
#endif

#ifndef NDEBUG
    MSG_INFO1( spxout, spxout << " [mode: debug]" );
#else
    MSG_INFO1( spxout, spxout << " [mode: optimized]" );
#endif

    MSG_INFO1( spxout, spxout << " [precision: " << (int)sizeof(Real) << " byte]" );

#ifdef SOPLEX_WITH_GMP
#ifdef mpir_version
    MSG_INFO1( spxout, spxout << " [rational: MPIR " << mpir_version << "]" );
#else
    MSG_INFO1( spxout, spxout << " [rational: GMP " << gmp_version << "]" );
#endif
#else
    MSG_INFO1( spxout, spxout << " [rational: long double]" );
#endif

    MSG_INFO1( spxout, spxout << " [githash: " << getGitHash() << "]\n" );
  }


    template <>
  void SoPlex<Rational>::printVersion() const
  {
    // do not use preprocessor directives within the MSG_INFO1 macro
#if (SOPLEX_SUBVERSION > 0)
    MSG_INFO1( spxout, spxout << "SoPlex version " << SOPLEX_VERSION/100
               << "." << (SOPLEX_VERSION % 100)/10
               << "." << SOPLEX_VERSION % 10
               << "." << SOPLEX_SUBVERSION );
#else
    MSG_INFO1( spxout, spxout << "SoPlex version " << SOPLEX_VERSION/100
               << "." << (SOPLEX_VERSION % 100)/10
               << "." << SOPLEX_VERSION % 10 );
#endif

#ifndef NDEBUG
    MSG_INFO1( spxout, spxout << " [mode: debug]" );
#else
    MSG_INFO1( spxout, spxout << " [mode: optimized]" );
#endif

    MSG_INFO1( spxout, spxout << " [precision: " << (int)sizeof(Real) << " byte]" );

#ifdef SOPLEX_WITH_GMP
#ifdef mpir_version
    MSG_INFO1( spxout, spxout << " [rational: MPIR " << mpir_version << "]" );
#else
    MSG_INFO1( spxout, spxout << " [rational: GMP " << gmp_version << "]" );
#endif
#else
    MSG_INFO1( spxout, spxout << " [rational: long double]" );
#endif

    MSG_INFO1( spxout, spxout << " [githash: " << getGitHash() << "]\n" );
  }



  /// checks if real LP and rational LP are in sync; dimensions will always be compared,
  /// vector and matrix values only if the respective parameter is set to true.
  /// If quiet is set to true the function will only display which vectors are different.
  template <>
	bool SoPlex<Real>::areLPsInSync(const bool checkVecVals, const bool checkMatVals, const bool quiet) const
  {
    bool result = true;
    bool nRowsMatch = true;
    bool nColsMatch = true;
    bool rhsDimMatch = true;
    bool lhsDimMatch = true;
    bool maxObjDimMatch = true;
    bool upperDimMatch = true;
    bool lowerDimMatch = true;

    // compare number of Rows
    if( _realLP->nRows() != _rationalLP->nRows() )
      {
        MSG_INFO1( spxout, spxout << "The number of Rows in the Real LP does not match the one in the Rational LP."
                   << " Real LP: " << _realLP->nRows() << "  Rational LP: " << _rationalLP->nRows() << std::endl);
        result = false;
        nRowsMatch = false;
      }

    // compare number of Columns
    if( _realLP->nCols() != _rationalLP->nCols() )
      {
        MSG_INFO1( spxout, spxout << "The number of Columns in the Real LP does not match the one in the Rational LP."
                   << " Real LP: " << _realLP->nCols() << "  Rational LP: " << _rationalLP->nCols() << std::endl);
        result = false;
        nColsMatch = false;
      }

    // compare number of nonZeros
    if( _realLP->nNzos() != _rationalLP->nNzos() )
      {
        MSG_INFO1( spxout, spxout << "The number of nonZeros in the Real LP does not match the one in the Rational LP."
                   << " Real LP: " << _realLP->nNzos() << "  Rational LP: " << _rationalLP->nNzos() << std::endl);
        result = false;
      }

    // compare the dimensions of the right hand side vectors
    if( _realLP->rhs().dim() != _rationalLP->rhs().dim() )
      {
        MSG_INFO1( spxout, spxout << "The dimension of the right hand side vector of the Real LP does not match the one of the Rational LP."
                   << " Real LP: " << _realLP->rhs().dim() << "  Rational LP: " << _rationalLP->rhs().dim() << std::endl);
        result = false;
        rhsDimMatch = false;

      }

    // compare the dimensions of the left hand side vectors
    if( _realLP->lhs().dim() != _rationalLP->lhs().dim() )
      {
        MSG_INFO1( spxout, spxout << "The dimension of the left hand side vector of the Real LP does not match the one of the Rational LP."
                   << " Real LP: " << _realLP->lhs().dim() << "  Rational LP: " << _rationalLP->lhs().dim() << std::endl);
        result = false;
        lhsDimMatch = false;
      }

    // compare the dimensions of the objective function vectors
    if( _realLP->maxObj().dim() != _rationalLP->maxObj().dim() )
      {
        MSG_INFO1( spxout, spxout << "The dimension of the objective function vector of the Real LP does not match the one of the Rational LP."
                   << " Real LP: " << _realLP->maxObj().dim() << "  Rational LP: " << _rationalLP->maxObj().dim() << std::endl);
        result = false;
        maxObjDimMatch = false;
      }

    // compare the sense
    if( (int)_realLP->spxSense() != (int)_rationalLP->spxSense() )
      {
        MSG_INFO1( spxout, spxout << "The objective function sense of the Real LP does not match the one of the Rational LP."
                   << " Real LP: " << (_realLP->spxSense() == SPxLPReal::MINIMIZE ? "MIN" : "MAX")
                   << "  Rational LP: " << (_rationalLP->spxSense() == SPxLPRational::MINIMIZE ? "MIN" : "MAX") << std::endl);
        result = false;
      }

    // compare the dimensions of upper bound vectors
    if( _realLP->upper().dim() != _rationalLP->upper().dim() )
      {
        MSG_INFO1( spxout, spxout << "The dimension of the upper bound vector of the Real LP does not match the one of the Rational LP."
                   << " Real LP: " << _realLP->upper().dim() << "  Rational LP: " << _rationalLP->upper().dim() << std::endl);
        result = false;
        upperDimMatch = false;
      }

    // compare the dimensions of the objective function vectors
    if( _realLP->lower().dim() != _rationalLP->lower().dim() )
      {
        MSG_INFO1( spxout, spxout << "The dimension of the lower bound vector of the Real LP does not match the one of the Rational LP."
                   << " Real LP: " << _realLP->lower().dim() << "  Rational LP: " << _rationalLP->lower().dim() << std::endl);
        result = false;
        lowerDimMatch = false;
      }

    // compares the values of the rhs, lhs, maxObj, upper, lower vectors
    if( checkVecVals )
      {
        bool rhsValMatch = true;
        bool lhsValMatch = true;
        bool maxObjValMatch = true;
        bool upperValMatch = true;
        bool lowerValMatch = true;

        // compares the values of the right hand side vectors
        if( rhsDimMatch )
          {
            for( int i = 0; i < _realLP->rhs().dim(); i++ )
              {
                if( (GE(_realLP->rhs()[i], realParam(SoPlex<Real>::INFTY)) != (_rationalLP->rhs()[i] >= _rationalPosInfty))
                    || (LT(_realLP->rhs()[i], realParam(SoPlex<Real>::INFTY)) && _rationalLP->rhs()[i] < _rationalPosInfty
                        && !_rationalLP->rhs()[i].isAdjacentTo((double)_realLP->rhs()[i])) )
                  {
                    if( !quiet )
                      {
                        MSG_INFO1( spxout, spxout << "Entries number " << i << " of the right hand side vectors don't match."
                                   << " Real LP: " << _realLP->rhs()[i] << "  Rational LP: " << _rationalLP->rhs()[i] << std::endl);
                      }
                    rhsValMatch = false;
                    result = false;
                  }
              }

            if( !rhsValMatch && quiet )
              {
                MSG_INFO1( spxout, spxout << "The values of the right hand side vectors don't match." << std::endl );
              }
          }

        // compares the values of the left hand side vectors
        if( lhsDimMatch )
          {
            for( int i = 0; i < _realLP->lhs().dim(); i++ )
              {
                if( (LE(_realLP->lhs()[i], -realParam(SoPlex<Real>::INFTY)) != (_rationalLP->lhs()[i] <= _rationalNegInfty))
                    || (GT(_realLP->lhs()[i], -realParam(SoPlex<Real>::INFTY)) && _rationalLP->lhs()[i] > _rationalNegInfty
                        && !_rationalLP->lhs()[i].isAdjacentTo((double)_realLP->lhs()[i])) )
                  {
                    if( !quiet )
                      {
                        MSG_INFO1( spxout, spxout << "Entries number " << i << " of the left hand side vectors don't match."
                                   << " Real LP: " << _realLP->lhs()[i] << "  Rational LP: " << _rationalLP->lhs()[i] << std::endl);
                      }
                    lhsValMatch = false;
                    result = false;
                  }
              }

            if( !lhsValMatch && quiet )
              {
                MSG_INFO1( spxout, spxout << "The values of the left hand side vectors don't match." << std::endl );
              }
          }

        // compares the values of the objective function vectors
        if( maxObjDimMatch )
          {
            for( int i = 0; i < _realLP->maxObj().dim(); i++ )
              {
                if( !_rationalLP->maxObj()[i].isAdjacentTo((double)_realLP->maxObj()[i]) )
                  {
                    if( !quiet )
                      {
                        MSG_INFO1( spxout, spxout << "Entries number " << i << " of the objective function vectors don't match."
                                   << " Real LP: " << _realLP->maxObj()[i] << "  Rational LP: " << _rationalLP->maxObj()[i] << std::endl);
                      }
                    maxObjValMatch = false;
                    result = false;
                  }
              }

            if( !maxObjValMatch && quiet )
              {
                MSG_INFO1( spxout, spxout << "The values of the objective function vectors don't match." << std::endl );
              }
          }

        // compares the values of the upper bound vectors
        if( upperDimMatch )
          {
            for( int i = 0; i < _realLP->upper().dim(); i++ )
              {
                if( (GE(_realLP->upper()[i], realParam(SoPlex<Real>::INFTY)) != (_rationalLP->upper()[i] >= _rationalPosInfty))
                    || (LT(_realLP->upper()[i], realParam(SoPlex<Real>::INFTY)) && _rationalLP->upper()[i] < _rationalPosInfty
                        && !_rationalLP->upper()[i].isAdjacentTo((double)_realLP->upper()[i])) )
                  {
                    if( !quiet )
                      {
                        MSG_INFO1( spxout, spxout << "Entries number " << i << " of the upper bound vectors don't match."
                                   << " Real LP: " << _realLP->upper()[i] << "  Rational LP: " << _rationalLP->upper()[i] << std::endl);
                      }
                    upperValMatch = false;
                    result = false;
                  }
              }

            if( !upperValMatch && quiet )
              {
                MSG_INFO1( spxout, spxout << "The values of the upper bound vectors don't match." << std::endl );
              }
          }

        // compares the values of the lower bound vectors
        if( lowerDimMatch )
          {
            for( int i = 0; i < _realLP->lower().dim(); i++ )
              {
                if( (LE(_realLP->lower()[i], -realParam(SoPlex<Real>::INFTY)) != (_rationalLP->lower()[i] <= _rationalNegInfty))
                    || (GT(_realLP->lower()[i], -realParam(SoPlex<Real>::INFTY)) && _rationalLP->lower()[i] > _rationalNegInfty
                        && !_rationalLP->lower()[i].isAdjacentTo((double)_realLP->lower()[i])) )
                  {
                    if( !quiet )
                      {
                        MSG_INFO1( spxout, spxout << "Entries number " << i << " of the lower bound vectors don't match."
                                   << " Real LP: " << _realLP->lower()[i] << "  Rational LP: " << _rationalLP->lower()[i] << std::endl);
                      }
                    lowerValMatch = false;
                    result = false;
                  }
              }

            if( !lowerValMatch && quiet )
              {
                MSG_INFO1( spxout, spxout << "The values of the lower bound vectors don't match." << std::endl );
              }
          }
      }

    // compare the values of the matrix
    if( checkMatVals && nRowsMatch && nColsMatch )
      {
        bool matrixValMatch = true;

        for( int i = 0; i < _realLP->nCols() ; i++ )
          {
            for( int j = 0;j < _realLP->nRows() ; j++ )
              {
                if( !_rationalLP->colVector(i)[j].isAdjacentTo((double)_realLP->colVector(i)[j]) )
                  {
                    if( !quiet )
                      {
                        MSG_INFO1( spxout, spxout << "Entries number " << j << " of column number " << i << " don't match."
                                   << " Real LP: " << _realLP->colVector(i)[j] << "  Rational LP: " << _rationalLP->colVector(i)[j] << std::endl);
                      }
                    matrixValMatch = false;
                    result = false;
                  }
              }
          }

        if( !matrixValMatch && quiet )
          {
            MSG_INFO1( spxout, spxout << "The values of the matrices don't match." << std::endl );
          }
      }

    return result;
  }



  /// set the random seed of the solver instance
  template <>
  void SoPlex<Real>::setRandomSeed(unsigned int seed)
  {
    _solver.random.setSeed(seed);
  }



  /// returns the current random seed of the solver instance or the one stored in the settings
  template <>
  unsigned int  SoPlex<Real>::randomSeed() const
  {
    return _solver.random.getSeed();
  }



  /// extends sparse vector to hold newmax entries if and only if it holds no more free entries
  template <>
  void SoPlex<Real>::_ensureDSVectorRationalMemory(DSVectorRational& vec, const int newmax) const
  {
    assert(newmax > vec.size());
    if( vec.size() >= vec.max() )
      vec.setMax(newmax);
  }



  /// creates a permutation for removing rows/columns from an array of indices
  template <>
  void SoPlex<Real>::_idxToPerm(int* idx, int idxSize, int* perm, int permSize) const
  {
    assert(idx != 0);
    assert(idxSize >= 0);
    assert(perm != 0);
    assert(permSize >= 0);

    for( int i = 0; i < permSize; i++ )
      perm[i] = i;

    for( int i = 0; i < idxSize; i++ )
      {
        assert(idx[i] >= 0);
        assert(idx[i] < permSize);
        perm[idx[i]] = -1;
      }
  }



  /// creates a permutation for removing rows/columns from a range of indices
  template <>
  void SoPlex<Real>::_rangeToPerm(int start, int end, int* perm, int permSize) const
  {
    assert(perm != 0);
    assert(permSize >= 0);

    for( int i = 0; i < permSize; i++ )
      perm[i] = (i < start || i > end) ? i : -1;
  }



  /// checks consistency
  template <>
	bool SoPlex<Real>::_isConsistent() const
  {
    assert(_statistics != 0);
    assert(_currentSettings != 0);

    assert(_realLP != 0);
    assert(_rationalLP != 0 || intParam(SoPlex<Real>::SYNCMODE) == SYNCMODE_ONLYREAL);

    assert(_realLP != &_solver || _isRealLPLoaded);
    assert(_realLP == &_solver || !_isRealLPLoaded);

    assert(!_hasBasis || _isRealLPLoaded || _basisStatusRows.size() == numRowsT());
    assert(!_hasBasis || _isRealLPLoaded || _basisStatusCols.size() == numColsT());
    assert(_rationalLUSolver.status() == SLinSolverRational::UNLOADED || _hasBasis);
    assert(_rationalLUSolver.status() == SLinSolverRational::UNLOADED || _rationalLUSolver.dim() == _rationalLUSolverBind.size());
    assert(_rationalLUSolver.status() == SLinSolverRational::UNLOADED || _rationalLUSolver.dim() == numRowsT());

    assert(_rationalLP == 0 || _colTypes.size() == numColsRational());
    assert(_rationalLP == 0 || _rowTypes.size() == numRowsRational());

    return true;
  }



  /// should solving process be stopped?
  template <>
	bool SoPlex<Real>::_isSolveStopped(bool& stoppedTime, bool& stoppedIter) const
  {
    assert(_statistics != 0);

    stoppedTime = (realParam(TIMELIMIT) < realParam(INFTY) && _statistics->solvingTime->time() >= realParam(TIMELIMIT));
    stoppedIter = (intParam(ITERLIMIT) >= 0 && _statistics->iterations >= intParam(ITERLIMIT))
      || (intParam(REFLIMIT) >= 0 && _statistics->refinements >= intParam(REFLIMIT))
      || (intParam(STALLREFLIMIT) >= 0 && _statistics->stallRefinements >= intParam(STALLREFLIMIT));

    return stoppedTime || stoppedIter;
  }



  /// determines RangeType from real bounds
  template <>
  typename SoPlex<Real>::RangeType SoPlex<Real>::_rangeTypeReal(const Real& lower, const Real& upper) const
  {
    assert(lower <= upper);

    if( lower <= -infinity )
      {
        if( upper >= infinity )
          return RANGETYPE_FREE;
        else
          return RANGETYPE_UPPER;
      }
    else
      {
        if( upper >= infinity )
          return RANGETYPE_LOWER;
        else if( lower == upper )
          return RANGETYPE_FIXED;
        else
          return RANGETYPE_BOXED;
      }
  }



  /// determines RangeType from rational bounds
  template <>
  typename SoPlex<Real>::RangeType SoPlex<Real>::_rangeTypeRational(const Rational& lower, const Rational& upper) const
  {
    assert(lower <= upper);

    if( lower <= _rationalNegInfty )
      {
        if( upper >= _rationalPosInfty )
          return RANGETYPE_FREE;
        else
          return RANGETYPE_UPPER;
      }
    else
      {
        if( upper >= _rationalPosInfty )
          return RANGETYPE_LOWER;
        else if( lower == upper )
          return RANGETYPE_FIXED;
        else
          return RANGETYPE_BOXED;
      }
  }



  /// switches RANGETYPE_LOWER to RANGETYPE_UPPER and vice versa
  template <>
  typename SoPlex<Real>::RangeType SoPlex<Real>::_switchRangeType(const SoPlex<Real>::RangeType& rangeType) const
  {
    if( rangeType == RANGETYPE_LOWER )
      return RANGETYPE_UPPER;
    else if( rangeType == RANGETYPE_UPPER )
      return RANGETYPE_LOWER;
    else
      return rangeType;
  }



  /// checks whether RangeType corresponds to finite lower bound
  template <>
	bool SoPlex<Real>::_lowerFinite(const RangeType& rangeType) const
  {
    return (rangeType == RANGETYPE_LOWER || rangeType == RANGETYPE_BOXED || rangeType == RANGETYPE_FIXED);
  }



  /// checks whether RangeType corresponds to finite upper bound
  template <>
	bool SoPlex<Real>::_upperFinite(const RangeType& rangeType) const
  {
    return (rangeType == RANGETYPE_UPPER || rangeType == RANGETYPE_BOXED || rangeType == RANGETYPE_FIXED);
  }



  /// adds a single row to the real LP and adjusts basis
  template <>
  void SoPlex<Real>::_addRowReal(const LPRowReal& lprow)
  {
    assert(_realLP != 0);

    bool scale = _realLP->isScaled();
    _realLP->addRow(lprow, scale);

    if( _isRealLPLoaded )
      _hasBasis = (_solver.basis().status() > SPxBasis<Real>::NO_PROBLEM);
    else if( _hasBasis )
      _basisStatusRows.append(SPxSolver<Real>::BASIC);

    _rationalLUSolver.clear();
  }



  /// adds a single row to the real LP and adjusts basis
  template <>
  void SoPlex<Real>::_addRowReal(Real lhs, const SVectorReal& lprow, Real rhs)
  {
    assert(_realLP != 0);

    bool scale = _realLP->isScaled();
    _realLP->addRow(lhs, lprow, rhs, scale);

    if( _isRealLPLoaded )
      _hasBasis = (_solver.basis().status() > SPxBasis<Real>::NO_PROBLEM);
    else if( _hasBasis )
      _basisStatusRows.append(SPxSolver<Real>::BASIC);

    _rationalLUSolver.clear();
  }



  /// adds multiple rows to the real LP and adjusts basis
  template <>
  void SoPlex<Real>::_addRowsReal(const LPRowSetReal& lprowset)
  {
    assert(_realLP != 0);

    bool scale = _realLP->isScaled();
    _realLP->addRows(lprowset, scale);

    if( _isRealLPLoaded )
      _hasBasis = (_solver.basis().status() > SPxBasis<Real>::NO_PROBLEM);
    else if( _hasBasis )
      _basisStatusRows.append(lprowset.num(), SPxSolver<Real>::BASIC);

    _rationalLUSolver.clear();
  }


  /// adds a single column to the real LP and adjusts basis
  template <>
  void SoPlex<Real>::_addColReal(const LPColReal& lpcol)
  {
    assert(_realLP != 0);

    bool scale = _realLP->isScaled();
    _realLP->addCol(lpcol, scale);

    if( _isRealLPLoaded )
      _hasBasis = (_solver.basis().status() > SPxBasis<Real>::NO_PROBLEM);
    else if( _hasBasis )
      {
        if( lpcol.lower() > -realParam(SoPlex<Real>::INFTY) )
          _basisStatusCols.append(SPxSolver<Real>::ON_LOWER);
        else if( lpcol.upper() < realParam(SoPlex<Real>::INFTY) )
          _basisStatusCols.append(SPxSolver<Real>::ON_UPPER);
        else
          _basisStatusCols.append(SPxSolver<Real>::ZERO);
      }

    _rationalLUSolver.clear();
  }



  /// adds a single column to the real LP and adjusts basis
  template <>
  void SoPlex<Real>::_addColReal(Real obj, Real lower, const SVectorReal& lpcol, Real upper)
  {
    assert(_realLP != 0);

    bool scale = _realLP->isScaled();
    _realLP->addCol(obj, lower, lpcol, upper, scale);

    if( _isRealLPLoaded )
      _hasBasis = (_solver.basis().status() > SPxBasis<Real>::NO_PROBLEM);
    else if( _hasBasis )
      _basisStatusRows.append(SPxSolver<Real>::BASIC);

    _rationalLUSolver.clear();
  }



  /// adds multiple columns to the real LP and adjusts basis
  template <>
  void SoPlex<Real>::_addColsReal(const LPColSetReal& lpcolset)
  {
    assert(_realLP != 0);

    bool scale = _realLP->isScaled();
    _realLP->addCols(lpcolset, scale);

    if( _isRealLPLoaded )
      _hasBasis = (_solver.basis().status() > SPxBasis<Real>::NO_PROBLEM);
    else if( _hasBasis )
      {
        for( int i = 0; i < lpcolset.num(); i++ )
          {
            if( lpcolset.lower(i) > -realParam(SoPlex<Real>::INFTY) )
              _basisStatusCols.append(SPxSolver<Real>::ON_LOWER);
            else if( lpcolset.upper(i) < realParam(SoPlex<Real>::INFTY) )
              _basisStatusCols.append(SPxSolver<Real>::ON_UPPER);
            else
              _basisStatusCols.append(SPxSolver<Real>::ZERO);
          }
      }

    _rationalLUSolver.clear();
  }


  /// replaces row \p i with \p lprow and adjusts basis
  template <>
  void SoPlex<Real>::_changeRowReal(int i, const LPRowReal& lprow)
  {
    assert(_realLP != 0);

    bool scale = _realLP->isScaled();
    _realLP->changeRow(i, lprow, scale);

    if( _isRealLPLoaded )
      _hasBasis = (_solver.basis().status() > SPxBasis<Real>::NO_PROBLEM);
    else if( _hasBasis )
      {
        if( _basisStatusRows[i] != SPxSolver<Real>::BASIC )
          _hasBasis = false;
        else if( _basisStatusRows[i] == SPxSolver<Real>::ON_LOWER && lprow.lhs() <= -realParam(SoPlex<Real>::INFTY) )
          _basisStatusRows[i] = (lprow.rhs() < realParam(SoPlex<Real>::INFTY)) ? SPxSolver<Real>::ON_UPPER : SPxSolver<Real>::ZERO;
        else if( _basisStatusRows[i] == SPxSolver<Real>::ON_UPPER && lprow.rhs() >= realParam(SoPlex<Real>::INFTY) )
          _basisStatusRows[i] = (lprow.lhs() > -realParam(SoPlex<Real>::INFTY)) ? SPxSolver<Real>::ON_LOWER : SPxSolver<Real>::ZERO;
      }

    _rationalLUSolver.clear();
  }



  /// changes left-hand side vector for constraints to \p lhs and adjusts basis
  template <>
  void SoPlex<Real>::_changeLhsReal(const VectorReal& lhs)
  {
    assert(_realLP != 0);

    bool scale = _realLP->isScaled();
    _realLP->changeLhs(lhs, scale);

    if( _isRealLPLoaded )
      _hasBasis = (_solver.basis().status() > SPxBasis<Real>::NO_PROBLEM);
    else if( _hasBasis )
      {
        for( int i = numRowsT() - 1; i >= 0; i-- )
          {
            if( _basisStatusRows[i] == SPxSolver<Real>::ON_LOWER && lhs[i] <= -realParam(SoPlex<Real>::INFTY) )
              _basisStatusRows[i] = (rhsReal(i) < realParam(SoPlex<Real>::INFTY)) ? SPxSolver<Real>::ON_UPPER : SPxSolver<Real>::ZERO;
          }
      }
  }



  /// changes left-hand side of row \p i to \p lhs and adjusts basis
  template <>
  void SoPlex<Real>::_changeLhsReal(int i, const Real& lhs)
  {
    assert(_realLP != 0);

    bool scale = _realLP->isScaled();
    _realLP->changeLhs(i, lhs, scale);

    if( _isRealLPLoaded )
      {
        _hasBasis = (_solver.basis().status() > SPxBasis<Real>::NO_PROBLEM);
      }
    else if( _hasBasis && _basisStatusRows[i] == SPxSolver<Real>::ON_LOWER && lhs <= -realParam(SoPlex<Real>::INFTY) )
      _basisStatusRows[i] = (rhsReal(i) < realParam(SoPlex<Real>::INFTY)) ? SPxSolver<Real>::ON_UPPER : SPxSolver<Real>::ZERO;

  }



  /// changes right-hand side vector to \p rhs and adjusts basis
  template <>
  void SoPlex<Real>::_changeRhsReal(const VectorReal& rhs)
  {
    assert(_realLP != 0);

    bool scale = _realLP->isScaled();
    _realLP->changeRhs(rhs, scale);

    if( _isRealLPLoaded )
      {
        _hasBasis = (_solver.basis().status() > SPxBasis<Real>::NO_PROBLEM);
      }
    else if( _hasBasis )
      {
        for( int i = numRowsT() - 1; i >= 0; i-- )
          {
            if( _basisStatusRows[i] == SPxSolver<Real>::ON_UPPER && rhs[i] >= realParam(SoPlex<Real>::INFTY) )
              _basisStatusRows[i] = (lhsReal(i) > -realParam(SoPlex<Real>::INFTY)) ? SPxSolver<Real>::ON_LOWER : SPxSolver<Real>::ZERO;
          }
      }
  }



  /// changes right-hand side of row \p i to \p rhs and adjusts basis
  template <>
  void SoPlex<Real>::_changeRhsReal(int i, const Real& rhs)
  {
    assert(_realLP != 0);

    bool scale = _realLP->isScaled();
    _realLP->changeRhs(i, rhs, scale);

    if( _isRealLPLoaded )
      {
        _hasBasis = (_solver.basis().status() > SPxBasis<Real>::NO_PROBLEM);
      }
    else if( _hasBasis && _basisStatusRows[i] == SPxSolver<Real>::ON_UPPER && rhs >= realParam(SoPlex<Real>::INFTY) )
      _basisStatusRows[i] = (lhsReal(i) > -realParam(SoPlex<Real>::INFTY)) ? SPxSolver<Real>::ON_LOWER : SPxSolver<Real>::ZERO;
  }



  /// changes left- and right-hand side vectors and adjusts basis
  template <>
  void SoPlex<Real>::_changeRangeReal(const VectorReal& lhs, const VectorReal& rhs)
  {
    assert(_realLP != 0);

    bool scale = _realLP->isScaled();
    _realLP->changeRange(lhs, rhs, scale);

    if( _isRealLPLoaded )
      {
        _hasBasis = (_solver.basis().status() > SPxBasis<Real>::NO_PROBLEM);
      }
    else if( _hasBasis )
      {
        for( int i = numRowsT() - 1; i >= 0; i-- )
          {
            if( _basisStatusRows[i] == SPxSolver<Real>::ON_LOWER && lhs[i] <= -realParam(SoPlex<Real>::INFTY) )
              _basisStatusRows[i] = (rhs[i] < realParam(SoPlex<Real>::INFTY)) ? SPxSolver<Real>::ON_UPPER : SPxSolver<Real>::ZERO;
            else if( _basisStatusRows[i] == SPxSolver<Real>::ON_UPPER && rhs[i] >= realParam(SoPlex<Real>::INFTY) )
              _basisStatusRows[i] = (lhs[i] > -realParam(SoPlex<Real>::INFTY)) ? SPxSolver<Real>::ON_LOWER : SPxSolver<Real>::ZERO;
          }
      }
  }



  /// changes left- and right-hand side of row \p i and adjusts basis
  template <>
  void SoPlex<Real>::_changeRangeReal(int i, const Real& lhs, const Real& rhs)
  {
    assert(_realLP != 0);

    bool scale = _realLP->isScaled();
    _realLP->changeRange(i, lhs, rhs, scale);

    if( _isRealLPLoaded )
      {
        _hasBasis = (_solver.basis().status() > SPxBasis<Real>::NO_PROBLEM);
      }
    else if( _hasBasis )
      {
        if( _basisStatusRows[i] == SPxSolver<Real>::ON_LOWER && lhs <= -realParam(SoPlex<Real>::INFTY) )
          _basisStatusRows[i] = (rhs < realParam(SoPlex<Real>::INFTY)) ? SPxSolver<Real>::ON_UPPER : SPxSolver<Real>::ZERO;
        else if( _basisStatusRows[i] == SPxSolver<Real>::ON_UPPER && rhs >= realParam(SoPlex<Real>::INFTY) )
          _basisStatusRows[i] = (lhs > -realParam(SoPlex<Real>::INFTY)) ? SPxSolver<Real>::ON_LOWER : SPxSolver<Real>::ZERO;
      }
  }



  /// replaces column \p i with \p lpcol and adjusts basis
  template <>
  void SoPlex<Real>::_changeColReal(int i, const LPColReal& lpcol)
  {
    assert(_realLP != 0);

    bool scale = _realLP->isScaled();
    _realLP->changeCol(i, lpcol, scale);

    if( _isRealLPLoaded )
      {
        _hasBasis = (_solver.basis().status() > SPxBasis<Real>::NO_PROBLEM);
      }
    else if( _hasBasis )
      {
        if( _basisStatusCols[i] == SPxSolver<Real>::BASIC )
          _hasBasis = false;
        else if( _basisStatusCols[i] == SPxSolver<Real>::ON_LOWER && lpcol.lower() <= -realParam(SoPlex<Real>::INFTY) )
          _basisStatusCols[i] = (lpcol.upper() < realParam(SoPlex<Real>::INFTY)) ? SPxSolver<Real>::ON_UPPER : SPxSolver<Real>::ZERO;
        else if( _basisStatusCols[i] == SPxSolver<Real>::ON_UPPER && lpcol.upper() >= realParam(SoPlex<Real>::INFTY) )
          _basisStatusCols[i] = (lpcol.lower() > -realParam(SoPlex<Real>::INFTY)) ? SPxSolver<Real>::ON_LOWER : SPxSolver<Real>::ZERO;
      }

    _rationalLUSolver.clear();
  }



  /// changes vector of lower bounds to \p lower and adjusts basis
  template <>
  void SoPlex<Real>::_changeLowerReal(const VectorReal& lower)
  {
    assert(_realLP != 0);

    bool scale = _realLP->isScaled();
    _realLP->changeLower(lower, scale);

    if( _isRealLPLoaded )
      {
        _hasBasis = (_solver.basis().status() > SPxBasis<Real>::NO_PROBLEM);
      }
    else if( _hasBasis )
      {
        for( int i = numColsT() - 1; i >= 0; i-- )
          {
            if( _basisStatusCols[i] == SPxSolver<Real>::ON_LOWER && lower[i] <= -realParam(SoPlex<Real>::INFTY) )
              _basisStatusCols[i] = (upperReal(i) < realParam(SoPlex<Real>::INFTY)) ? SPxSolver<Real>::ON_UPPER : SPxSolver<Real>::ZERO;
          }
      }
  }



  /// changes lower bound of column i to \p lower and adjusts basis
  template <>
  void SoPlex<Real>::_changeLowerReal(int i, const Real& lower)
  {
    assert(_realLP != 0);

    bool scale = _realLP->isScaled();
    _realLP->changeLower(i, lower, scale);

    if( _isRealLPLoaded )
      {
        _hasBasis = (_solver.basis().status() > SPxBasis<Real>::NO_PROBLEM);
      }
    else if( _hasBasis && _basisStatusCols[i] == SPxSolver<Real>::ON_LOWER && lower <= -realParam(SoPlex<Real>::INFTY) )
      _basisStatusCols[i] = (upperReal(i) < realParam(SoPlex<Real>::INFTY)) ? SPxSolver<Real>::ON_UPPER : SPxSolver<Real>::ZERO;
  }



  /// changes vector of upper bounds to \p upper and adjusts basis
  template <>
  void SoPlex<Real>::_changeUpperReal(const VectorReal& upper)
  {
    assert(_realLP != 0);

    bool scale = _realLP->isScaled();
    _realLP->changeUpper(upper, scale);

    if( _isRealLPLoaded )
      {
        _hasBasis = (_solver.basis().status() > SPxBasis<Real>::NO_PROBLEM);
      }
    else if( _hasBasis )
      {
        for( int i = numColsT() - 1; i >= 0; i-- )
          {
            if( _basisStatusCols[i] == SPxSolver<Real>::ON_UPPER && upper[i] >= realParam(SoPlex<Real>::INFTY) )
              _basisStatusCols[i] = (lowerReal(i) > -realParam(SoPlex<Real>::INFTY)) ? SPxSolver<Real>::ON_LOWER : SPxSolver<Real>::ZERO;
          }
      }
  }



  /// changes \p i 'th upper bound to \p upper and adjusts basis
  template <>
  void SoPlex<Real>::_changeUpperReal(int i, const Real& upper)
  {
    assert(_realLP != 0);

    bool scale = _realLP->isScaled();
    _realLP->changeUpper(i, upper, scale);

    if( _isRealLPLoaded )
      {
        _hasBasis = (_solver.basis().status() > SPxBasis<Real>::NO_PROBLEM);
      }
    else if( _hasBasis &&  _basisStatusCols[i] == SPxSolver<Real>::ON_UPPER && upper >= realParam(SoPlex<Real>::INFTY) )
      _basisStatusCols[i] = (lowerReal(i) > -realParam(SoPlex<Real>::INFTY)) ? SPxSolver<Real>::ON_LOWER : SPxSolver<Real>::ZERO;
  }



  /// changes vectors of column bounds to \p lower and \p upper and adjusts basis
  template <>
  void SoPlex<Real>::_changeBoundsReal(const VectorReal& lower, const VectorReal& upper)
  {
    assert(_realLP != 0);

    bool scale = _realLP->isScaled();
    _realLP->changeBounds(lower, upper, scale);

    if( _isRealLPLoaded )
      {
        _hasBasis = (_solver.basis().status() > SPxBasis<Real>::NO_PROBLEM);
      }
    else if( _hasBasis )
      {
        for( int i = numColsT() - 1; i >= 0; i-- )
          {
            if( _basisStatusCols[i] == SPxSolver<Real>::ON_LOWER && lower[i] <= -realParam(SoPlex<Real>::INFTY) )
              _basisStatusCols[i] = (upper[i] < realParam(SoPlex<Real>::INFTY)) ? SPxSolver<Real>::ON_UPPER : SPxSolver<Real>::ZERO;
            else if( _basisStatusCols[i] == SPxSolver<Real>::ON_UPPER && upper[i] >= realParam(SoPlex<Real>::INFTY) )
              _basisStatusCols[i] = (lower[i] > -realParam(SoPlex<Real>::INFTY)) ? SPxSolver<Real>::ON_LOWER : SPxSolver<Real>::ZERO;
          }
      }
  }



  /// changes bounds of column \p i to \p lower and \p upper and adjusts basis
  template <>
  void SoPlex<Real>::_changeBoundsReal(int i, const Real& lower, const Real& upper)
  {
    assert(_realLP != 0);

    bool scale = _realLP->isScaled();
    _realLP->changeBounds(i, lower, upper, scale);

    if( _isRealLPLoaded )
      {
        _hasBasis = (_solver.basis().status() > SPxBasis<Real>::NO_PROBLEM);
      }
    else if( _hasBasis )
      {
        if( _basisStatusCols[i] == SPxSolver<Real>::ON_LOWER && lower <= -realParam(SoPlex<Real>::INFTY) )
          _basisStatusCols[i] = (upper < realParam(SoPlex<Real>::INFTY)) ? SPxSolver<Real>::ON_UPPER : SPxSolver<Real>::ZERO;
        else if( _basisStatusCols[i] == SPxSolver<Real>::ON_UPPER && upper >= realParam(SoPlex<Real>::INFTY) )
          _basisStatusCols[i] = (lower > -realParam(SoPlex<Real>::INFTY)) ? SPxSolver<Real>::ON_LOWER : SPxSolver<Real>::ZERO;
      }
  }



  /// changes matrix entry in row \p i and column \p j to \p val and adjusts basis
  template <>
  void SoPlex<Real>::_changeElementReal(int i, int j, const Real& val)
  {
    assert(_realLP != 0);

    bool scale = _realLP->isScaled();
    _realLP->changeElement(i, j, val, scale);

    if( _isRealLPLoaded )
      {
        _hasBasis = (_solver.basis().status() > SPxBasis<Real>::NO_PROBLEM);
      }
    else if( _hasBasis )
      {
        if( _basisStatusRows[i] != SPxSolver<Real>::BASIC && _basisStatusCols[i] == SPxSolver<Real>::BASIC )
          _hasBasis = false;
      }

    _rationalLUSolver.clear();
  }



  /// removes row \p i and adjusts basis
  template <>
  void SoPlex<Real>::_removeRowReal(int i)
  {
    assert(_realLP != 0);

    _realLP->removeRow(i);

    if( _isRealLPLoaded )
      {
        _hasBasis = (_solver.basis().status() > SPxBasis<Real>::NO_PROBLEM);
      }
    else if( _hasBasis )
      {
        if( _basisStatusRows[i] != SPxSolver<Real>::BASIC )
          _hasBasis = false;
        else
          {
            _basisStatusRows[i] = _basisStatusRows[_basisStatusRows.size() - 1];
            _basisStatusRows.removeLast();
          }
      }

    _rationalLUSolver.clear();
  }



  /// removes all rows with an index \p i such that \p perm[i] < 0; upon completion, \p perm[i] >= 0 indicates the
  /// new index where row \p i has been moved to; note that \p perm must point to an array of size at least
  /// #numRowsT()
  template <>
  void SoPlex<Real>::_removeRowsReal(int perm[])
  {
    assert(_realLP != 0);

    _realLP->removeRows(perm);

    if( _isRealLPLoaded )
      {
        _hasBasis = (_solver.basis().status() > SPxBasis<Real>::NO_PROBLEM);
      }
    else if( _hasBasis )
      {
        for( int i = numRowsT() - 1; i >= 0 && _hasBasis; i-- )
          {
            if( perm[i] < 0 && _basisStatusRows[i] != SPxSolver<Real>::BASIC )
              _hasBasis = false;
            else if( perm[i] >= 0 && perm[i] != i )
              {
                assert(perm[i] < numRowsT());
                assert(perm[perm[i]] < 0);

                _basisStatusRows[perm[i]] = _basisStatusRows[i];
              }
          }

        if( _hasBasis )
          _basisStatusRows.reSize(numRowsT());
      }

    _rationalLUSolver.clear();
  }



  /// removes column i
  template <>
  void SoPlex<Real>::_removeColReal(int i)
  {
    assert(_realLP != 0);

    _realLP->removeCol(i);

    if( _isRealLPLoaded )
      {
        _hasBasis = (_solver.basis().status() > SPxBasis<Real>::NO_PROBLEM);
      }
    else if( _hasBasis )
      {
        if( _basisStatusCols[i] == SPxSolver<Real>::BASIC )
          _hasBasis = false;
        else
          {
            _basisStatusCols[i] = _basisStatusCols[_basisStatusCols.size() - 1];
            _basisStatusCols.removeLast();
          }
      }

    _rationalLUSolver.clear();
  }



  /// removes all columns with an index \p i such that \p perm[i] < 0; upon completion, \p perm[i] >= 0 indicates the
  /// new index where column \p i has been moved to; note that \p perm must point to an array of size at least
  /// #numColsT()
  template <>
  void SoPlex<Real>::_removeColsReal(int perm[])
  {
    assert(_realLP != 0);

    _realLP->removeCols(perm);

    if( _isRealLPLoaded )
      {
        _hasBasis = (_solver.basis().status() > SPxBasis<Real>::NO_PROBLEM);
      }
    else if( _hasBasis )
      {
        for( int i = numColsT() - 1; i >= 0 && _hasBasis; i-- )
          {
            if( perm[i] < 0 && _basisStatusCols[i] == SPxSolver<Real>::BASIC )
              _hasBasis = false;
            else if( perm[i] >= 0 && perm[i] != i )
              {
                assert(perm[i] < numColsT());
                assert(perm[perm[i]] < 0);

                _basisStatusCols[perm[i]] = _basisStatusCols[i];
              }
          }

        if( _hasBasis )
          _basisStatusCols.reSize(numColsT());
      }

    _rationalLUSolver.clear();
  }



  /// invalidates solution
  template <>
  void SoPlex<Real>::_invalidateSolution()
  {
    ///@todo maybe this should be done individually at the places when this method is called
    _status = SPxSolver<Real>::UNKNOWN;

    _solReal.invalidate();
    _hasSolReal = false;

    _solRational.invalidate();
    _hasSolRational = false;
  }



  /// enables simplifier and scaler
  template <>
  void SoPlex<Real>::_enableSimplifierAndScaler()
  {
    // type of simplifier
    switch( intParam(SoPlex<Real>::SIMPLIFIER) )
      {
      case SIMPLIFIER_OFF:
        _simplifier = 0;
        break;
      case SIMPLIFIER_AUTO:
        _simplifier = &_simplifierMainSM;
        assert(_simplifier != 0);
        _simplifier->setMinReduction(realParam(MINRED));
        break;
      default:
        break;
      }

    // type of scaler
    switch( intParam(SoPlex<Real>::SCALER) )
      {
      case SCALER_OFF:
        _scaler = 0;
        break;
      case SCALER_UNIEQUI:
        _scaler = &_scalerUniequi;
        break;
      case SCALER_BIEQUI:
        _scaler = &_scalerBiequi;
        break;
      case SCALER_GEO1:
        _scaler = &_scalerGeo1;
        break;
      case SCALER_GEO8:
        _scaler = &_scalerGeo8;
        break;
      case SCALER_LEASTSQ:
        _scaler = &_scalerLeastsq;
        break;
      case SCALER_GEOEQUI:
        _scaler = &_scalerGeoequi;
        break;
      default:
        break;
      }
  }



  /// disables simplifier and scaler
  template <>
  void SoPlex<Real>::_disableSimplifierAndScaler()
  {
    _simplifier = 0;

    // preserve scaler when persistent scaling is used
    if( !_isRealLPScaled )
      _scaler = 0;
    else
      assert(boolParam(SoPlex<Real>::PERSISTENTSCALING));
  }



  /// ensures that the rational LP is available; performs no sync
  template <>
  void SoPlex<Real>::_ensureRationalLP()
  {
    if( _rationalLP == 0 )
      {
        spx_alloc(_rationalLP);
        _rationalLP = new (_rationalLP) SPxLPRational();
        _rationalLP->setOutstream(spxout);
      }
  }



  /// ensures that the real LP and the basis are loaded in the solver; performs no sync
  template <>
  void SoPlex<Real>::_ensureRealLPLoaded()
  {
    if( !_isRealLPLoaded )
      {
        assert(_realLP != &_solver);

        _solver.loadLP(*_realLP);
        _realLP->~SPxLPReal();
        spx_free(_realLP);
        _realLP = &_solver;
        _isRealLPLoaded = true;

        if( _hasBasis )
          {
            ///@todo this should not fail even if the basis is invalid (wrong dimension or wrong number of basic
            ///      entries); fix either in SPxSolver or in SPxBasis
            assert(_basisStatusRows.size() == numRowsT());
            assert(_basisStatusCols.size() == numColsT());
            _solver.setBasis(_basisStatusRows.get_const_ptr(), _basisStatusCols.get_const_ptr());
            _hasBasis = (_solver.basis().status() > SPxBasis<Real>::NO_PROBLEM);
          }
      }
  }



  /// call floating-point solver and update statistics on iterations etc.
  template <>
  void SoPlex<Real>::_solveRealLPAndRecordStatistics()
  {
    bool _hadBasis = _hasBasis;

    // set time and iteration limit
    if( intParam(SoPlex<Real>::ITERLIMIT) < realParam(SoPlex<Real>::INFTY) )
      _solver.setTerminationIter(intParam(SoPlex<Real>::ITERLIMIT) - _statistics->iterations);
    else
      _solver.setTerminationIter(-1);
    if( realParam(SoPlex<Real>::TIMELIMIT) < realParam(SoPlex<Real>::INFTY) )
      _solver.setTerminationTime(realParam(SoPlex<Real>::TIMELIMIT) - _statistics->solvingTime->time());
    else
      _solver.setTerminationTime(realParam(SoPlex<Real>::INFTY));

    // ensure that tolerances are not too small
    if( _solver.feastol() < 1e-12 )
      _solver.setFeastol(1e-12);
    if( _solver.opttol() < 1e-12 )
      _solver.setOpttol(1e-12);

    // set correct representation
    if( (intParam(SoPlex<Real>::REPRESENTATION) == SoPlex<Real>::REPRESENTATION_COLUMN
         || (intParam(SoPlex<Real>::REPRESENTATION) == SoPlex<Real>::REPRESENTATION_AUTO && (_solver.nCols() + 1) * realParam(SoPlex<Real>::REPRESENTATION_SWITCH) >= (_solver.nRows() + 1)))
        && _solver.rep() != SPxSolver<Real>::COLUMN )
      {
        _solver.setRep(SPxSolver<Real>::COLUMN);
      }
    else if( (intParam(SoPlex<Real>::REPRESENTATION) == SoPlex<Real>::REPRESENTATION_ROW
              || (intParam(SoPlex<Real>::REPRESENTATION) == SoPlex<Real>::REPRESENTATION_AUTO && (_solver.nCols() + 1) * realParam(SoPlex<Real>::REPRESENTATION_SWITCH) < (_solver.nRows() + 1)))
             &&_solver.rep() != SPxSolver<Real>::ROW )
      {
        _solver.setRep(SPxSolver<Real>::ROW);
      }

    // set correct type
    if( ((intParam(ALGORITHM) == SoPlex<Real>::ALGORITHM_PRIMAL && _solver.rep() == SPxSolver<Real>::COLUMN)
         || (intParam(ALGORITHM) == SoPlex<Real>::ALGORITHM_DUAL && _solver.rep() == SPxSolver<Real>::ROW))
        && _solver.type() != SPxSolver<Real>::ENTER )
      {
        _solver.setType(SPxSolver<Real>::ENTER);
      }
    else if( ((intParam(ALGORITHM) == SoPlex<Real>::ALGORITHM_DUAL && _solver.rep() == SPxSolver<Real>::COLUMN)
              || (intParam(ALGORITHM) == SoPlex<Real>::ALGORITHM_PRIMAL && _solver.rep() == SPxSolver<Real>::ROW))
             && _solver.type() != SPxSolver<Real>::LEAVE )
      {
        _solver.setType(SPxSolver<Real>::LEAVE);
      }

    // set pricing modes
    _solver.setSparsePricingFactor(realParam(SoPlex<Real>::SPARSITY_THRESHOLD));
    if( (intParam(SoPlex<Real>::HYPER_PRICING) == SoPlex<Real>::HYPER_PRICING_ON)
        || ((intParam(SoPlex<Real>::HYPER_PRICING) == SoPlex<Real>::HYPER_PRICING_AUTO)
            && (_solver.nRows() + _solver.nCols() > HYPERPRICINGTHRESHOLD )) )
      _solver.hyperPricing(true);
    else if( intParam(SoPlex<Real>::HYPER_PRICING) == SoPlex<Real>::HYPER_PRICING_OFF )
      _solver.hyperPricing(false);

    _solver.setNonzeroFactor(realParam(SoPlex<Real>::REFAC_BASIS_NNZ));
    _solver.setFillFactor(realParam(SoPlex<Real>::REFAC_UPDATE_FILL));
    _solver.setMemFactor(realParam(SoPlex<Real>::REFAC_MEM_FACTOR));

    // call floating-point solver and catch exceptions
    _statistics->simplexTime->start();
    try
      {
        _solver.solve();
      }
    catch( const SPxException& E )
      {
        MSG_INFO1( spxout, spxout << "Caught exception <" << E.what() << "> while solving real LP.\n" );
        _status = SPxSolver<Real>::ERROR;
      }
    catch( ... )
      {
        MSG_INFO1( spxout, spxout << "Caught unknown exception while solving real LP.\n" );
        _status = SPxSolver<Real>::ERROR;
      }
    _statistics->simplexTime->stop();

    // invalidate rational factorization of basis if pivots have been performed
    if( _solver.iterations() > 0 )
      _rationalLUSolver.clear();

    // record statistics
    _statistics->iterations += _solver.iterations();
    _statistics->iterationsPrimal += _solver.primalIterations();
    _statistics->iterationsFromBasis += _hadBasis ? _solver.iterations() : 0;
    _statistics->iterationsPolish += _solver.polishIterations();
    _statistics->boundflips += _solver.boundFlips();
    _statistics->luFactorizationTimeReal += _slufactor.getFactorTime();
    _statistics->luSolveTimeReal += _slufactor.getSolveTime();
    _statistics->luFactorizationsReal += _slufactor.getFactorCount();
    _statistics->luSolvesReal += _slufactor.getSolveCount();
    _slufactor.resetCounters();

    _statistics->degenPivotsPrimal += _solver.primalDegeneratePivots();
    _statistics->degenPivotsDual += _solver.dualDegeneratePivots();
    _statistics->sumDualDegen += _solver.sumDualDegeneracy();
    _statistics->sumPrimalDegen += _solver.sumPrimalDegeneracy();
  }



  /// reads real LP in LP or MPS format from file and returns true on success; gets row names, column names, and
  /// integer variables if desired
  template <>
	bool SoPlex<Real>::_readFileReal(const char* filename, NameSet* rowNames, NameSet* colNames, DIdxSet* intVars)
  {
    assert(_realLP != 0);

    // clear statistics
    _statistics->clearAllData();

    // update status
    clearBasis();
    _invalidateSolution();
    _status = SPxSolver<Real>::UNKNOWN;

    // start timing
    _statistics->readingTime->start();

    // read
    bool success = _realLP->readFile(filename, rowNames, colNames, intVars);

    // stop timing
    _statistics->readingTime->stop();

    if( success )
      {
        setIntParam(SoPlex<Real>::OBJSENSE, (_realLP->spxSense() == SPxLPReal::MAXIMIZE ? SoPlex<Real>::OBJSENSE_MAXIMIZE : SoPlex<Real>::OBJSENSE_MINIMIZE), true);
        _realLP->changeObjOffset(realParam(SoPlex<Real>::OBJ_OFFSET));

        // if sync mode is auto, we have to copy the (rounded) real LP to the rational LP; this is counted to sync time
        // and not to reading time
        if( intParam(SoPlex<Real>::SYNCMODE) == SYNCMODE_AUTO )
          _syncLPRational();
      }
    else
      clearLPReal();

    return success;
  }



  /// reads rational LP in LP or MPS format from file and returns true on success; gets row names, column names, and
  /// integer variables if desired
  template <>
	bool SoPlex<Real>::_readFileRational(const char* filename, NameSet* rowNames, NameSet* colNames, DIdxSet* intVars)
  {
    // clear statistics
    _statistics->clearAllData();

    // start timing
    _statistics->readingTime->start();

    // update status
    clearBasis();
    _invalidateSolution();
    _status = SPxSolver<Real>::UNKNOWN;

    // read
    _ensureRationalLP();
    bool success = _rationalLP->readFile(filename, rowNames, colNames, intVars);

    // stop timing
    _statistics->readingTime->stop();

    if( success )
      {
        setIntParam(SoPlex<Real>::OBJSENSE, (_rationalLP->spxSense() == SPxLPRational::MAXIMIZE ? SoPlex<Real>::OBJSENSE_MAXIMIZE : SoPlex<Real>::OBJSENSE_MINIMIZE), true);
        _rationalLP->changeObjOffset(realParam(SoPlex<Real>::OBJ_OFFSET));
        _recomputeRangeTypesRational();

        // if sync mode is auto, we have to copy the (rounded) real LP to the rational LP; this is counted to sync time
        // and not to reading time
        if( intParam(SoPlex<Real>::SYNCMODE) == SYNCMODE_AUTO )
          _syncLPReal();
        // if a rational LP file is read, but only the (rounded) real LP should be kept, we have to free the rational LP
        else if( intParam(SoPlex<Real>::SYNCMODE) == SYNCMODE_ONLYREAL )
          {
            _syncLPReal();
            _rationalLP->~SPxLPRational();
            spx_free(_rationalLP);
          }
      }
    else
      clearLPRational();

    return success;
  }



  /// completes range type arrays after adding columns and/or rows
  template <>
  void SoPlex<Real>::_completeRangeTypesRational()
  {
    // we use one method for bot columns and rows, because during column/row addition, rows/columns can be added
    // implicitly
    for( int i = _colTypes.size(); i < numColsT(); i++ )
      _colTypes.append(_rangeTypeRational(_rationalLP->lower(i), _rationalLP->upper(i)));
    for( int i = _rowTypes.size(); i < numRowsT(); i++ )
      _rowTypes.append(_rangeTypeRational(_rationalLP->lhs(i), _rationalLP->rhs(i)));
  }



  /// recomputes range types from scratch using real LP
  template <>
  void SoPlex<Real>::_recomputeRangeTypesReal()
  {
    _rowTypes.reSize(numRowsT());
    for( int i = 0; i < numRowsT(); i++ )
      _rowTypes[i] = _rangeTypeReal(_realLP->lhs(i), _realLP->rhs(i));
    _colTypes.reSize(numColsT());
    for( int i = 0; i < numColsT(); i++ )
      _colTypes[i] = _rangeTypeReal(_realLP->lower(i), _realLP->upper(i));
  }



  /// recomputes range types from scratch using rational LP
  template <>
  void SoPlex<Real>::_recomputeRangeTypesRational()
  {
    _rowTypes.reSize(numRowsT());
    for( int i = 0; i < numRowsT(); i++ )
      _rowTypes[i] = _rangeTypeRational(_rationalLP->lhs(i), _rationalLP->rhs(i));
    _colTypes.reSize(numColsT());
    for( int i = 0; i < numColsT(); i++ )
      _colTypes[i] = _rangeTypeRational(_rationalLP->lower(i), _rationalLP->upper(i));
  }



  /// synchronizes real LP with rational LP, i.e., copies (rounded) rational LP into real LP, without looking at the sync mode
  template <>
  void SoPlex<Real>::_syncLPReal(bool time)
  {
    // start timing
    if( time )
      _statistics->syncTime->start();

    // copy LP
    if( _isRealLPLoaded )
      _solver.loadLP((SPxLPReal)(*_rationalLP));
    else
      *_realLP = *_rationalLP;

    ///@todo try loading old basis
    _hasBasis = false;
    _rationalLUSolver.clear();

    // stop timing
    if( time )
      _statistics->syncTime->stop();
  }



  /// synchronizes rational LP with real LP, i.e., copies real LP to rational LP, without looking at the sync mode
  template <>
  void SoPlex<Real>::_syncLPRational(bool time)
  {
    // start timing
    if( time )
      _statistics->syncTime->start();

    // copy LP
    _ensureRationalLP();
    *_rationalLP = *_realLP;
    _recomputeRangeTypesRational();

    // stop timing
    if( time )
      _statistics->syncTime->stop();
  }



  /// synchronizes rational solution with real solution, i.e., copies (rounded) rational solution to real solution
  template <>
  void SoPlex<Real>::_syncRealSolution()
  {
    if( _hasSolRational && !_hasSolReal )
      {
        _solReal = _solRational;
        _hasSolReal = true;
      }
  }



  /// synchronizes real solution with rational solution, i.e., copies real solution to rational solution
  template <>
  void SoPlex<Real>::_syncRationalSolution()
  {
    if( _hasSolReal && !_hasSolRational )
      {
        _solRational = _solReal;
        _hasSolRational = true;
      }
  }



  /// returns pointer to a constant unit vector available until destruction of the SoPlex class
  template <>
  const UnitVectorRational* SoPlex<Real>::_unitVectorRational(const int i)
  {
    assert(i >= 0);

    if( i < 0 )
      return 0;
    else if( i >= _unitMatrixRational.size() )
      _unitMatrixRational.append(i + 1 - _unitMatrixRational.size(), (UnitVectorRational*)0);
    assert(i < _unitMatrixRational.size());

    if( _unitMatrixRational[i] == 0 )
      {
        spx_alloc(_unitMatrixRational[i]);
        new (_unitMatrixRational[i]) UnitVectorRational(i);
      }
    assert(_unitMatrixRational[i] != 0);

    return _unitMatrixRational[i];
  }



  /// parses one line in a settings file and returns true on success; note that the string is modified
  template <>
	bool SoPlex<Real>::_parseSettingsLine(char* line, const int lineNumber)
  {
    assert(line != 0);

    // find the start of the parameter type
    while( *line == ' ' || *line == '\t' || *line == '\r' )
      line++;
    if( *line == '\0' || *line == '\n' || *line == '#' )
      return true;
    char* paramTypeString = line;

    // find the end of the parameter type
    while( *line != ' ' && *line != '\t' && *line != '\r' && *line != '\n' && *line != '#' && *line != '\0' && *line != ':' )
      line++;
    if( *line == ':' )
      {
        *line = '\0';
        line++;
      }
    else
      {
        *line = '\0';
        line++;

        // search for the ':' char in the line
        while( *line == ' ' || *line == '\t' || *line == '\r' )
          line++;
        if( *line != ':' )
          {
            MSG_INFO1( spxout, spxout << "Error parsing settings file: no ':' separating parameter type and name in line " << lineNumber << ".\n" );
            return false;
          }
        line++;
      }

    // find the start of the parameter name
    while( *line == ' ' || *line == '\t' || *line == '\r' )
      line++;
    if( *line == '\0' || *line == '\n' || *line == '#' )
      {
        MSG_INFO1( spxout, spxout << "Error parsing settings file: no parameter name in line " << lineNumber << ".\n");
        return false;
      }
    char* paramName = line;

    // find the end of the parameter name
    while( *line != ' ' && *line != '\t' && *line != '\r' && *line != '\n' && *line != '#' && *line != '\0' && *line != '=' )
      line++;
    if( *line == '=' )
      {
        *line = '\0';
        line++;
      }
    else
      {
        *line = '\0';
        line++;

        // search for the '=' char in the line
        while( *line == ' ' || *line == '\t' || *line == '\r' )
          line++;
        if( *line != '=' )
          {
            MSG_INFO1( spxout, spxout << "Error parsing settings file: no '=' after parameter name in line " << lineNumber << ".\n" );
            return false;
          }
        line++;
      }

    // find the start of the parameter value string
    while( *line == ' ' || *line == '\t' || *line == '\r' )
      line++;
    if( *line == '\0' || *line == '\n' || *line == '#' )
      {
        MSG_INFO1( spxout, spxout << "Error parsing settings file: no parameter value in line " << lineNumber << ".\n");
        return false;
      }
    char* paramValueString = line;

    // find the end of the parameter value string
    while( *line != ' ' && *line != '\t' && *line != '\r' && *line != '\n' && *line != '#' && *line != '\0' )
      line++;
    if( *line != '\0' )
      {
        // check, if the rest of the line is clean
        *line = '\0';
        line++;
        while( *line == ' ' || *line == '\t' || *line == '\r' )
          line++;
        if( *line != '\0' && *line != '\n' && *line != '#' )
          {
            MSG_INFO1( spxout, spxout << "Error parsing settings file: additional character '" << *line << "' after parameter value in line " << lineNumber << ".\n" );
            return false;
          }
      }

    // check whether we have a bool parameter
    if( strncmp(paramTypeString, "bool", 4) == 0 )
      {
        for( int param = 0; ; param++ )
          {
            if( param >= SoPlex<Real>::BOOLPARAM_COUNT )
              {
                MSG_INFO1( spxout, spxout << "Error parsing settings file: unknown parameter name <" << paramName << "> in line " << lineNumber << ".\n" );
                return false;
              }
            else if( strncmp(paramName, _currentSettings->boolParam.name[param].c_str(), SET_MAX_LINE_LEN) == 0 )
              {
                if( strncasecmp(paramValueString, "true", 4) == 0
                    || strncasecmp(paramValueString, "TRUE", 4) == 0
                    || strncasecmp(paramValueString, "t", 4) == 0
                    || strncasecmp(paramValueString, "T", 4) == 0
                    || strtol(paramValueString, NULL, 4) == 1 )
                  {
                    setBoolParam((SoPlex<Real>::BoolParam)param, true);
                    break;
                  }
                else if( strncasecmp(paramValueString, "false", 5) == 0
                         || strncasecmp(paramValueString, "FALSE", 5) == 0
                         || strncasecmp(paramValueString, "f", 5) == 0
                         || strncasecmp(paramValueString, "F", 5) == 0
                         || strtol(paramValueString, NULL, 5) == 0 )
                  {
                    setBoolParam((SoPlex<Real>::BoolParam)param, false);
                    break;
                  }
                else
                  {
                    MSG_INFO1( spxout, spxout << "Error parsing settings file: invalid value <" << paramValueString << "> for bool parameter <" << paramName << "> in line " << lineNumber << ".\n" );
                    return false;
                  }
              }
          }

        return true;
      }

    // check whether we have an integer parameter
    if( strncmp(paramTypeString, "int", 3) == 0 )
      {
        for( int param = 0; ; param++ )
          {
            if( param >= SoPlex<Real>::INTPARAM_COUNT )
              {
                MSG_INFO1( spxout, spxout << "Error parsing settings file: unknown parameter name <" << paramName << "> in line " << lineNumber << ".\n" );
                return false;
              }
            else if( strncmp(paramName, _currentSettings->intParam.name[param].c_str(), SET_MAX_LINE_LEN) == 0 )
              {
                int value;

                if( sscanf(paramValueString, "%d", &value) == 1 && setIntParam((SoPlex<Real>::IntParam)param, value, false) )
                  break;
                else
                  {
                    MSG_INFO1( spxout, spxout << "Error parsing settings file: invalid value <" << paramValueString << "> for int parameter <" << paramName << "> in line " << lineNumber << ".\n" );
                    return false;
                  }
              }
          }

        return true;
      }

    // check whether we have a real parameter
    if( strncmp(paramTypeString, "real", 4) == 0 )
      {
        for( int param = 0; ; param++ )
          {
            if( param >= SoPlex<Real>::REALPARAM_COUNT )
              {
                MSG_INFO1( spxout, spxout << "Error parsing settings file: unknown parameter name <" << paramName << "> in line " << lineNumber << ".\n" );
                return false;
              }
            else if( strncmp(paramName, _currentSettings->realParam.name[param].c_str(), SET_MAX_LINE_LEN) == 0 )
              {
                Real value;

                if( sscanf(paramValueString, "%" REAL_FORMAT, &value) == 1 && setRealParam((SoPlex<Real>::RealParam)param, value) )
                  break;
                else
                  {
                    MSG_INFO1( spxout, spxout << "Error parsing settings file: invalid value <" << paramValueString << "> for real parameter <" << paramName << "> in line " << lineNumber << ".\n" );
                    return false;
                  }
              }
          }

        return true;
      }

#ifdef SOPLEX_WITH_RATIONALPARAM
    // check whether we have a rational parameter
    if( strncmp(paramTypeString, "rational", 8) == 0 )
      {
        for( int param = 0; ; param++ )
          {
            if( param >= SoPlex<Real>::RATIONALPARAM_COUNT )
              {
                MSG_INFO1( spxout, spxout << "Error parsing settings file: unknown parameter name <" << paramName << "> in line " << lineNumber << ".\n" );
                return false;
              }
            else if( strncmp(paramName, _currentSettings->rationalParam.name[param].c_str(), SET_MAX_LINE_LEN) == 0 )
              {
                Rational value;

                if( readStringRational(paramValueString, value) && setRationalParam((SoPlex<Real>::RationalParam)param, value) )
                  break;
                else
                  {
                    MSG_INFO1( spxout, spxout << "Error parsing settings file: invalid value <" << paramValueString << "> for rational parameter <" << paramName << "> in line " << lineNumber << ".\n" );
                    return false;
                  }
              }
          }

        return true;
      }
#endif

    // check whether we have the random seed
    if( strncmp(paramTypeString, "uint", 4) == 0 )
      {
        if( strncmp(paramName, "random_seed", 11) == 0 )
          {
            unsigned int value;
            char format[SPX_MAXSTRLEN];
            spxSnprintf(format, sizeof(format), "%%%du", (int) sizeof(value)-1);

            if( sscanf(paramValueString, format, &value) == 1 )
              {
                setRandomSeed(value);
                return true;
              }
          }

        MSG_INFO1( spxout, spxout << "Error parsing settings file for uint parameter <random_seed>.\n" );
        return false;
      }

    MSG_INFO1( spxout, spxout << "Error parsing settings file: invalid parameter type <" << paramTypeString << "> for parameter <" << paramName << "> in line " << lineNumber << ".\n" );

    return false;
  }
} // namespace soplex
#endif
