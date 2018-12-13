/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the class library                   */
/*       SoPlex --- the Sequential object-oriented simPlex.                  */
/*                                                                           */
/*    Copyright (C) 1996-2018 Konrad-Zuse-Zentrum                            */
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

#ifdef _MSC_VER
#define strncasecmp _strnicmp
#endif

namespace soplex
{
  // template <>
  // int SoPlexBase<Real>::intParam(const IntParam param) const;

  // // template <>
	// // bool SoPlexBase<Real>::areLPsInSync(const bool checkVecVals, const bool checkMatVals, const bool quiet) const;

  // template <>
  // typename SoPlexBase<Real>::Settings& SoPlexBase<Real>::Settings::operator=(const Settings& settings);

  // template <>
  // const SVectorReal& SoPlexBase<Real>::colVectorRealInternal(int i) const;

  // template <>
  // const Rational& SoPlexBase<Real>::maxObjRational(int i) const;

  // template <>
  // typename SPxSolverBase<Real>::VarStatus SoPlexBase<Real>::basisRowStatus(int row) const;

  // template <>
  // typename SPxSolverBase<Real>::VarStatus SoPlexBase<Real>::basisColStatus(int col) const;

  // // template <>
	// // bool SoPlexBase<Real>::setSettings(const Settings& newSettings, const bool init);

  template <>
	bool SoPlexBase<Real>::saveSettingsFile(const char* filename, const bool onlyChanged) const;

  template <>
  void SoPlexBase<Real>::printShortStatistics(std::ostream& os);

  // // template <>
  // // void SoPlexBase<Real>::printStatus(std::ostream& os, typename SPxSolverBase<Real>::Status stat);

  // // template <>
  // // void SoPlexBase<Real>::setRandomSeed(unsigned int seed);

  // template <>
  // void SoPlexBase<Real>::_idxToPerm(int* idx, int idxSize, int* perm, int permSize) const;

  // template <>
  // void SoPlexBase<Real>::_rangeToPerm(int start, int end, int* perm, int permSize) const;

  // template <>
	// bool SoPlexBase<Real>::_isConsistent() const;

  // // template <>
  // // typename SoPlexBase<Real>::RangeType SoPlexBase<Real>::_rangeTypeReal(const Real& lower, const Real& upper) const;

  // template <>
  // void SoPlexBase<Real>::_addRowReal(const LPRowReal& lprow);

  // template <>
	// bool SoPlexBase<Real>::setIntParam(const IntParam param, const int value, const bool init);

  // template <>
  // typename SoPlexBase<Real>::RangeType SoPlexBase<Real>::_rangeTypeRational(const Rational& lower, const Rational& upper) const;

  // template <>
  // void SoPlexBase<Real>::_addRowReal(Real lhs, const SVectorReal& lprow, Real rhs);

  // template <>
  // void SoPlexBase<Real>::_addRowsReal(const LPRowSetReal& lprowset);

  // template <>
  // void SoPlexBase<Real>::_addColReal(const LPColReal& lpcol);

  // template <>
  // void SoPlexBase<Real>::_addColReal(Real obj, Real lower, const SVectorReal& lpcol, Real upper);

  // template <>
  // void SoPlexBase<Real>::_addColsReal(const LPColSetReal& lpcolset);

  // template <>
  // void SoPlexBase<Real>::_changeRowReal(int i, const LPRowReal& lprow);

  // template <>
  // void SoPlexBase<Real>::_changeLhsReal(const VectorReal& lhs);

  // template <>
  // void SoPlexBase<Real>::_changeLhsReal(int i, const Real& lhs);

  // template <>
  // void SoPlexBase<Real>::_changeRhsReal(const VectorReal& rhs);

  // template <>
  // void SoPlexBase<Real>::_changeRhsReal(int i, const Real& rhs);

  // template <>
  // void SoPlexBase<Real>::_changeRangeReal(const VectorReal& lhs, const VectorReal& rhs);

  // template <>
  // void SoPlexBase<Real>::_changeRangeReal(int i, const Real& lhs, const Real& rhs);

  // template <>
  // void SoPlexBase<Real>::_changeColReal(int i, const LPColReal& lpcol);

  // template <>
  // void SoPlexBase<Real>::_changeLowerReal(const VectorReal& lower);

  // template <>
  // void SoPlexBase<Real>::_changeLowerReal(int i, const Real& lower);

  // template <>
  // void SoPlexBase<Real>::_changeUpperReal(const VectorReal& upper);

  // template <>
  // void SoPlexBase<Real>::_changeUpperReal(int i, const Real& upper);

  // template <>
  // void SoPlexBase<Real>::_changeBoundsReal(const VectorReal& lower, const VectorReal& upper);

  // template <>
  // void SoPlexBase<Real>::_changeBoundsReal(int i, const Real& lower, const Real& upper);

  // template <>
  // void SoPlexBase<Real>::_changeElementReal(int i, int j, const Real& val);

  // template <>
  // void SoPlexBase<Real>::_removeRowReal(int i);

  // template <>
  // void SoPlexBase<Real>::_removeRowsReal(int perm[]);

  // template <>
  // void SoPlexBase<Real>::_removeColReal(int i);

  // template <>
  // void SoPlexBase<Real>::_removeColsReal(int perm[]);

  // template <>
  // void SoPlexBase<Real>::_invalidateSolution();

  // template <>
  // void SoPlexBase<Real>::_ensureRationalLP();

  // template <>
  // void SoPlexBase<Real>::_ensureRealLPLoaded();

  // template <>
	// bool SoPlexBase<Real>::_readFileReal(const char* filename, NameSet* rowNames, NameSet* colNames, DIdxSet* intVars);

  // template <>
  // void SoPlexBase<Real>::_completeRangeTypesRational();

  // template <>
  // void SoPlexBase<Real>::_recomputeRangeTypesRational();

  // template <>
  // void SoPlexBase<Real>::_syncLPReal(bool time);

  // // template <>
  // // void SoPlexBase<Real>::_syncLPRational(bool time);

  // template <>
  // void SoPlexBase<Real>::_syncRealSolution();

  // template <>
  // void SoPlexBase<Real>::_syncRationalSolution();

  // template <>
  // const UnitVectorRational* SoPlexBase<Real>::_unitVectorRational(const int i);

  // template <>
	// bool SoPlexBase<Real>::_parseSettingsLine(char* line, const int lineNumber);

  // template <>
	// bool SoPlexBase<Real>::_readFileRational(const char* filename, NameSet* rowNames, NameSet* colNames, DIdxSet* intVars);

  // template <>
	// bool SoPlexBase<Real>::boolParam(const BoolParam param) const;

  template <>
  SoPlexBase<Real>::Settings::BoolParam::BoolParam();
  template <>
  typename SoPlexBase<Real>::Settings::BoolParam SoPlexBase<Real>::Settings::boolParam = BoolParam();

  // template <>
  // typename SPxSolverBase<Real>::Status SoPlexBase<Real>::status() const;

  // template <>
	// bool SoPlexBase<Real>::hasBasis() const;

  template <>
  typename SoPlexBase<Real>::Settings::IntParam SoPlexBase<Real>::Settings::intParam = IntParam();

  template <>
  typename SoPlexBase<Real>::Settings::RealParam SoPlexBase<Real>::Settings::realParam = RealParam();

  template <>
  Real SoPlexBase<Real>::realParam(const RealParam param) const;

  template <>
  SoPlexBase<Real>::Settings::Settings();


  template <>
  SoPlexBase<Real>::Settings::BoolParam::BoolParam() {
    // should lifting be used to reduce range of nonzero matrix coefficients?
    name[SoPlexBase<Real>::LIFTING] = "lifting";
    description[SoPlexBase<Real>::LIFTING] = "should lifting be used to reduce range of nonzero matrix coefficients?";
    defaultValue[SoPlexBase<Real>::LIFTING] = false;

    // should LP be transformed to equality form before a rational solve?
    name[SoPlexBase<Real>::EQTRANS] = "eqtrans";
    description[SoPlexBase<Real>::EQTRANS] = "should LP be transformed to equality form before a rational solve?";
    defaultValue[SoPlexBase<Real>::EQTRANS] = false;

    // should dual infeasibility be tested in order to try to return a dual solution even if primal infeasible?
    name[SoPlexBase<Real>::TESTDUALINF] = "testdualinf";
    description[SoPlexBase<Real>::TESTDUALINF] = "should dual infeasibility be tested in order to try to return a dual solution even if primal infeasible?";
    defaultValue[SoPlexBase<Real>::TESTDUALINF] = false;

    // should a rational factorization be performed after iterative refinement?
    name[SoPlexBase<Real>::RATFAC] = "ratfac";
    description[SoPlexBase<Real>::RATFAC] = "should a rational factorization be performed after iterative refinement?";
    defaultValue[SoPlexBase<Real>::RATFAC] = true;

    // should the decomposition based dual simplex be used to solve the LP? Setting this to true forces the solve mode to
    // SOLVEMODE_REAL and the basis representation to REPRESENTATION_ROW
    name[SoPlexBase<Real>::USEDECOMPDUALSIMPLEX] = "decompositiondualsimplex";
    description[SoPlexBase<Real>::USEDECOMPDUALSIMPLEX] = "should the decomposition based dual simplex be used to solve the LP?";
    defaultValue[SoPlexBase<Real>::USEDECOMPDUALSIMPLEX] = false;

    // should the degeneracy be computed for each basis?
    name[SoPlexBase<Real>::COMPUTEDEGEN] = "computedegen";
    description[SoPlexBase<Real>::COMPUTEDEGEN] = "should the degeneracy be computed for each basis?";
    defaultValue[SoPlexBase<Real>::COMPUTEDEGEN] = false;

    // should the dual of the complementary problem be used in the decomposition simplex?
    name[SoPlexBase<Real>::USECOMPDUAL] = "usecompdual";
    description[SoPlexBase<Real>::USECOMPDUAL] = "should the dual of the complementary problem be used in the decomposition simplex?";
    defaultValue[SoPlexBase<Real>::USECOMPDUAL] = false;

    /// should row and bound violations be computed explicitly in the update of reduced problem in the decomposition
    // simplex
    name[SoPlexBase<Real>::EXPLICITVIOL] = "explicitviol";
    description[SoPlexBase<Real>::EXPLICITVIOL] = "Should violations of the original problem be explicitly computed in the decomposition simplex?";
    defaultValue[SoPlexBase<Real>::EXPLICITVIOL] = false;

    // should cycling solutions be accepted during iterative refinement?
    name[SoPlexBase<Real>::ACCEPTCYCLING] = "acceptcycling";
    description[SoPlexBase<Real>::ACCEPTCYCLING] = "should cycling solutions be accepted during iterative refinement?";
    defaultValue[SoPlexBase<Real>::ACCEPTCYCLING] = false;

    // apply rational reconstruction after each iterative refinement?
    name[SoPlexBase<Real>::RATREC] = "ratrec";
    description[SoPlexBase<Real>::RATREC] = "apply rational reconstruction after each iterative refinement?";
    defaultValue[SoPlexBase<Real>::RATREC] = true;

    // round scaling factors for iterative refinement to powers of two?
    name[SoPlexBase<Real>::POWERSCALING] = "powerscaling";
    description[SoPlexBase<Real>::POWERSCALING] = "round scaling factors for iterative refinement to powers of two?";
    defaultValue[SoPlexBase<Real>::POWERSCALING] = true;

    // continue iterative refinement with exact basic solution if not optimal?
    name[SoPlexBase<Real>::RATFACJUMP] = "ratfacjump";
    description[SoPlexBase<Real>::RATFACJUMP] = "continue iterative refinement with exact basic solution if not optimal?";
    defaultValue[SoPlexBase<Real>::RATFACJUMP] = false;

    // use bound flipping also for row representation?
    name[SoPlexBase<Real>::ROWBOUNDFLIPS] = "rowboundflips";
    description[SoPlexBase<Real>::ROWBOUNDFLIPS] = "use bound flipping also for row representation?";
    defaultValue[SoPlexBase<Real>::ROWBOUNDFLIPS] = false;

    // use persistent scaling?
    name[SoPlexBase<Real>::PERSISTENTSCALING] = "persistentscaling";
    description[SoPlexBase<Real>::PERSISTENTSCALING] = "should persistent scaling be used?";
    defaultValue[SoPlexBase<Real>::PERSISTENTSCALING] = true;

    // perturb the entire problem or only the relevant bounds of s single pivot?
    name[SoPlexBase<Real>::FULLPERTURBATION] = "fullperturbation";
    description[SoPlexBase<Real>::FULLPERTURBATION] = "should perturbation be applied to the entire problem?";
    defaultValue[SoPlexBase<Real>::FULLPERTURBATION] = false;

    /// re-optimize the original problem to get a proof of infeasibility/unboundedness?
    name[SoPlexBase<Real>::ENSURERAY] = "ensureray";
    description[SoPlexBase<Real>::ENSURERAY] = "re-optimize the original problem to get a proof (ray) of infeasibility/unboundedness?";
    defaultValue[SoPlexBase<Real>::ENSURERAY] = false;
  }

  template <class R>
  SoPlexBase<R>::Settings::IntParam::IntParam() {
    // objective sense
    name[SoPlexBase<R>::OBJSENSE] = "objsense";
    description[SoPlexBase<R>::OBJSENSE] = "objective sense (-1 - minimize, +1 - maximize)";
    lower[SoPlexBase<R>::OBJSENSE] = -1;
    upper[SoPlexBase<R>::OBJSENSE] = 1;
    defaultValue[SoPlexBase<R>::OBJSENSE] = SoPlexBase<R>::OBJSENSE_MAXIMIZE;

    // type of computational form, i.e., column or row representation
    name[SoPlexBase<R>::REPRESENTATION] = "representation";
    description[SoPlexBase<R>::REPRESENTATION] = "type of computational form (0 - auto, 1 - column representation, 2 - row representation)";
    lower[SoPlexBase<R>::REPRESENTATION] = 0;
    upper[SoPlexBase<R>::REPRESENTATION] = 2;
    defaultValue[SoPlexBase<R>::REPRESENTATION] = SoPlexBase<R>::REPRESENTATION_AUTO;

    // type of algorithm, i.e., primal or dual
    name[SoPlexBase<R>::ALGORITHM] = "algorithm";
    description[SoPlexBase<R>::ALGORITHM] = "type of algorithm (0 - primal, 1 - dual)";
    lower[SoPlexBase<R>::ALGORITHM] = 0;
    upper[SoPlexBase<R>::ALGORITHM] = 1;
    defaultValue[SoPlexBase<R>::ALGORITHM] = SoPlexBase<R>::ALGORITHM_DUAL;

    // type of LU update
    name[SoPlexBase<R>::FACTOR_UPDATE_TYPE] = "factor_update_type";
    description[SoPlexBase<R>::FACTOR_UPDATE_TYPE] = "type of LU update (0 - eta update, 1 - Forrest-Tomlin update)";
    lower[SoPlexBase<R>::FACTOR_UPDATE_TYPE] = 0;
    upper[SoPlexBase<R>::FACTOR_UPDATE_TYPE] = 1;
    defaultValue[SoPlexBase<R>::FACTOR_UPDATE_TYPE] = SoPlexBase<R>::FACTOR_UPDATE_TYPE_FT;

    // maximum number of updates without fresh factorization
    name[SoPlexBase<R>::FACTOR_UPDATE_MAX] = "factor_update_max";
    description[SoPlexBase<R>::FACTOR_UPDATE_MAX] = "maximum number of LU updates without fresh factorization (0 - auto)";
    lower[SoPlexBase<R>::FACTOR_UPDATE_MAX] = 0;
    upper[SoPlexBase<R>::FACTOR_UPDATE_MAX] = INT_MAX;
    defaultValue[SoPlexBase<R>::FACTOR_UPDATE_MAX] = 0;

    // iteration limit (-1 if unlimited)
    name[SoPlexBase<R>::ITERLIMIT] = "iterlimit";
    description[SoPlexBase<R>::ITERLIMIT] = "iteration limit (-1 - no limit)";
    lower[SoPlexBase<R>::ITERLIMIT] = -1;
    upper[SoPlexBase<R>::ITERLIMIT] = INT_MAX;
    defaultValue[SoPlexBase<R>::ITERLIMIT] = -1;

    // refinement limit (-1 if unlimited)
    name[SoPlexBase<R>::REFLIMIT] = "reflimit";
    description[SoPlexBase<R>::REFLIMIT] = "refinement limit (-1 - no limit)";
    lower[SoPlexBase<R>::REFLIMIT] = -1;
    upper[SoPlexBase<R>::REFLIMIT] = INT_MAX;
    defaultValue[SoPlexBase<R>::REFLIMIT] = -1;

    // stalling refinement limit (-1 if unlimited)
    name[SoPlexBase<R>::STALLREFLIMIT] = "stallreflimit";
    description[SoPlexBase<R>::STALLREFLIMIT] = "stalling refinement limit (-1 - no limit)";
    lower[SoPlexBase<R>::STALLREFLIMIT] = -1;
    upper[SoPlexBase<R>::STALLREFLIMIT] = INT_MAX;
    defaultValue[SoPlexBase<R>::STALLREFLIMIT] = -1;

    // display frequency
    name[SoPlexBase<R>::DISPLAYFREQ] = "displayfreq";
    description[SoPlexBase<R>::DISPLAYFREQ] = "display frequency";
    lower[SoPlexBase<R>::DISPLAYFREQ] = 1;
    upper[SoPlexBase<R>::DISPLAYFREQ] = INT_MAX;
    defaultValue[SoPlexBase<R>::DISPLAYFREQ] = 200;

    // verbosity level
    name[SoPlexBase<R>::VERBOSITY] = "verbosity";
    description[SoPlexBase<R>::VERBOSITY] = "verbosity level (0 - error, 1 - warning, 2 - debug, 3 - normal, 4 - high, 5 - full)";
    lower[SoPlexBase<R>::VERBOSITY] = 0;
    upper[SoPlexBase<R>::VERBOSITY] = 5;
    defaultValue[SoPlexBase<R>::VERBOSITY] = SoPlexBase<R>::VERBOSITY_NORMAL;
    //_intParamDefault[SoPlexBase<R>::VERBOSITY] = SoPlexBase<R>::VERBOSITY_FULL;

    // type of simplifier
    name[SoPlexBase<R>::SIMPLIFIER] = "simplifier";
    description[SoPlexBase<R>::SIMPLIFIER] = "simplifier (0 - off, 1 - auto)";
    lower[SoPlexBase<R>::SIMPLIFIER] = 0;
    upper[SoPlexBase<R>::SIMPLIFIER] = 1;
    defaultValue[SoPlexBase<R>::SIMPLIFIER] = SoPlexBase<R>::SIMPLIFIER_AUTO;

    // type of scaler
    name[SoPlexBase<R>::SCALER] = "scaler";
    description[SoPlexBase<R>::SCALER] = "scaling (0 - off, 1 - uni-equilibrium, 2 - bi-equilibrium, 3 - geometric, 4 - iterated geometric, 5 - least squares, 6 - geometric-equilibrium)";
    lower[SoPlexBase<R>::SCALER] = 0;
    upper[SoPlexBase<R>::SCALER] = 6;
    defaultValue[SoPlexBase<R>::SCALER] = SoPlexBase<R>::SCALER_BIEQUI;

    // type of starter used to create crash basis
    name[SoPlexBase<R>::STARTER] = "starter";
    description[SoPlexBase<R>::STARTER] = "crash basis generated when starting from scratch (0 - none, 1 - weight, 2 - sum, 3 - vector)";
    lower[SoPlexBase<R>::STARTER] = 0;
    upper[SoPlexBase<R>::STARTER] = 3;
    defaultValue[SoPlexBase<R>::STARTER] = SoPlexBase<R>::STARTER_OFF;

    // type of pricer
    name[SoPlexBase<R>::PRICER] = "pricer";
    description[SoPlexBase<R>::PRICER] = "pricing method (0 - auto, 1 - dantzig, 2 - parmult, 3 - devex, 4 - quicksteep, 5 - steep)";
    lower[SoPlexBase<R>::PRICER] = 0;
    upper[SoPlexBase<R>::PRICER] = 5;
    defaultValue[SoPlexBase<R>::PRICER] = SoPlexBase<R>::PRICER_AUTO;

    // type of ratio test
    name[SoPlexBase<R>::RATIOTESTER] = "ratiotester";
    description[SoPlexBase<R>::RATIOTESTER] = "method for ratio test (0 - textbook, 1 - harris, 2 - fast, 3 - boundflipping)";
    lower[SoPlexBase<R>::RATIOTESTER] = 0;
    upper[SoPlexBase<R>::RATIOTESTER] = 3;
    defaultValue[SoPlexBase<R>::RATIOTESTER] = SoPlexBase<R>::RATIOTESTER_BOUNDFLIPPING;

    // mode for synchronizing real and rational LP
    name[SoPlexBase<R>::SYNCMODE] = "syncmode";
    description[SoPlexBase<R>::SYNCMODE] = "mode for synchronizing real and rational LP (0 - store only real LP, 1 - auto, 2 - manual)";
    lower[SoPlexBase<R>::SYNCMODE] = 0;
    upper[SoPlexBase<R>::SYNCMODE] = 2;
    defaultValue[SoPlexBase<R>::SYNCMODE] = SoPlexBase<R>::SYNCMODE_ONLYREAL;

    // mode for reading LP files
    name[SoPlexBase<R>::READMODE] = "readmode";
    description[SoPlexBase<R>::READMODE] = "mode for reading LP files (0 - floating-point, 1 - rational)";
    lower[SoPlexBase<R>::READMODE] = 0;
    upper[SoPlexBase<R>::READMODE] = 1;
    defaultValue[SoPlexBase<R>::READMODE] = SoPlexBase<R>::READMODE_REAL;

    // mode for iterative refinement strategy
    name[SoPlexBase<R>::SOLVEMODE] = "solvemode";
    description[SoPlexBase<R>::SOLVEMODE] = "mode for iterative refinement strategy (0 - floating-point solve, 1 - auto, 2 - exact rational solve)";
    lower[SoPlexBase<R>::SOLVEMODE] = 0;
    upper[SoPlexBase<R>::SOLVEMODE] = 2;
    defaultValue[SoPlexBase<R>::SOLVEMODE] = SoPlexBase<R>::SOLVEMODE_AUTO;

    // mode for iterative refinement strategy
    name[SoPlexBase<R>::CHECKMODE] = "checkmode";
    description[SoPlexBase<R>::CHECKMODE] = "mode for a posteriori feasibility checks (0 - floating-point check, 1 - auto, 2 - exact rational check)";
    lower[SoPlexBase<R>::CHECKMODE] = 0;
    upper[SoPlexBase<R>::CHECKMODE] = 2;
    defaultValue[SoPlexBase<R>::CHECKMODE] = SoPlexBase<R>::CHECKMODE_AUTO;

    // type of timing
    name[SoPlexBase<R>::TIMER] = "timer";
    description[SoPlexBase<R>::TIMER] = "type of timer (1 - cputime, aka. usertime, 2 - wallclock time, 0 - no timing)";
    lower[SoPlexBase<R>::TIMER] = 0;
    upper[SoPlexBase<R>::TIMER] = 2;
    defaultValue[SoPlexBase<R>::TIMER] = SoPlexBase<R>::TIMER_CPU;

    // mode for hyper sparse pricing
    name[SoPlexBase<R>::HYPER_PRICING] = "hyperpricing";
    description[SoPlexBase<R>::HYPER_PRICING] = "mode for hyper sparse pricing (0 - off, 1 - auto, 2 - always)";
    lower[SoPlexBase<R>::HYPER_PRICING] = 0;
    upper[SoPlexBase<R>::HYPER_PRICING] = 2;
    defaultValue[SoPlexBase<R>::HYPER_PRICING] = SoPlexBase<R>::HYPER_PRICING_AUTO;

    // minimum number of stalling refinements since last pivot to trigger rational factorization
    name[SoPlexBase<R>::RATFAC_MINSTALLS] = "ratfac_minstalls";
    description[SoPlexBase<R>::RATFAC_MINSTALLS] = "minimum number of stalling refinements since last pivot to trigger rational factorization";
    lower[SoPlexBase<R>::RATFAC_MINSTALLS] = 0;
    upper[SoPlexBase<R>::RATFAC_MINSTALLS] = INT_MAX;
    defaultValue[SoPlexBase<R>::RATFAC_MINSTALLS] = 2;

    // maximum number of conjugate gradient iterations in least square scaling
    name[SoPlexBase<R>::LEASTSQ_MAXROUNDS] = "leastsq_maxrounds";
    description[SoPlexBase<R>::LEASTSQ_MAXROUNDS] = "maximum number of conjugate gradient iterations in least square scaling";
    lower[SoPlexBase<R>::LEASTSQ_MAXROUNDS] = 0;
    upper[SoPlexBase<R>::LEASTSQ_MAXROUNDS] = INT_MAX;
    defaultValue[SoPlexBase<R>::LEASTSQ_MAXROUNDS] = 50;

    // mode for solution polishing
    name[SoPlexBase<R>::SOLUTION_POLISHING] = "solution_polishing";
    description[SoPlexBase<R>::SOLUTION_POLISHING] = "mode for solution polishing (0 - off, 1 - max basic slack, 2 - min basic slack)";
    lower[SoPlexBase<R>::SOLUTION_POLISHING] = 0;
    upper[SoPlexBase<R>::SOLUTION_POLISHING] = 2;
    defaultValue[SoPlexBase<R>::SOLUTION_POLISHING] = SoPlexBase<R>::POLISHING_OFF;

    // the number of iterations before the decomposition simplex initialisation is terminated.
    name[SoPlexBase<R>::DECOMP_ITERLIMIT] = "decomp_iterlimit";
    description[SoPlexBase<R>::DECOMP_ITERLIMIT] = "the number of iterations before the decomposition simplex initialisation solve is terminated";
    lower[SoPlexBase<R>::DECOMP_ITERLIMIT] = 1;
    upper[SoPlexBase<R>::DECOMP_ITERLIMIT] = INT_MAX;
    defaultValue[SoPlexBase<R>::DECOMP_ITERLIMIT] = 100;

    // maximum number of violated rows added in each iteration of the decomposition simplex
    name[SoPlexBase<R>::DECOMP_MAXADDEDROWS] = "decomp_maxaddedrows";
    description[SoPlexBase<R>::DECOMP_MAXADDEDROWS] = "maximum number of rows that are added to the reduced problem when using the decomposition based simplex";
    lower[SoPlexBase<R>::DECOMP_MAXADDEDROWS] = 1;
    upper[SoPlexBase<R>::DECOMP_MAXADDEDROWS] = INT_MAX;
    defaultValue[SoPlexBase<R>::DECOMP_MAXADDEDROWS] = 500;

    // maximum number of violated rows added in each iteration of the decomposition simplex
    name[SoPlexBase<R>::DECOMP_DISPLAYFREQ] = "decomp_displayfreq";
    description[SoPlexBase<R>::DECOMP_DISPLAYFREQ] = "the frequency that the decomposition based simplex status output is displayed.";
    lower[SoPlexBase<R>::DECOMP_DISPLAYFREQ] = 1;
    upper[SoPlexBase<R>::DECOMP_DISPLAYFREQ] = INT_MAX;
    defaultValue[SoPlexBase<R>::DECOMP_DISPLAYFREQ] = 50;

    // the verbosity of the decomposition based simplex
    name[SoPlexBase<R>::DECOMP_VERBOSITY] = "decomp_verbosity";
    description[SoPlexBase<R>::DECOMP_VERBOSITY] = "the verbosity of decomposition based simplex (0 - error, 1 - warning, 2 - debug, 3 - normal, 4 - high, 5 - full).";
    lower[SoPlexBase<R>::DECOMP_VERBOSITY] = 1;
    upper[SoPlexBase<R>::DECOMP_VERBOSITY] = 5;
    defaultValue[SoPlexBase<R>::DECOMP_VERBOSITY] = VERBOSITY_ERROR;

    // printing condition number during the solve
    name[SoPlexBase<Real>::PRINTBASISMETRIC] = "printbasismetric";
    description[SoPlexBase<Real>::PRINTBASISMETRIC] = "print basis metric during the solve (-1 - off, 0 - condition estimate , 1 - trace, 2 - determinant, 3 - condition)";
    lower[SoPlexBase<Real>::PRINTBASISMETRIC] = -1;
    upper[SoPlexBase<Real>::PRINTBASISMETRIC] = 3;
    defaultValue[SoPlexBase<Real>::PRINTBASISMETRIC] = -1;
  }

  template <class R>
  SoPlexBase<R>::Settings::RealParam::RealParam() {
    // primal feasibility tolerance
    name[SoPlexBase<R>::FEASTOL] = "feastol";
    description[SoPlexBase<R>::FEASTOL] = "primal feasibility tolerance";
    lower[SoPlexBase<R>::FEASTOL] = 0.0;
    upper[SoPlexBase<R>::FEASTOL] = 1.0;
    defaultValue[SoPlexBase<R>::FEASTOL] = 1e-6;

    // dual feasibility tolerance
    name[SoPlexBase<R>::OPTTOL] = "opttol";
    description[SoPlexBase<R>::OPTTOL] = "dual feasibility tolerance";
    lower[SoPlexBase<R>::OPTTOL] = 0.0;
    upper[SoPlexBase<R>::OPTTOL] = 1.0;
    defaultValue[SoPlexBase<R>::OPTTOL] = 1e-6;

    ///@todo define suitable values depending on Real type
    // general zero tolerance
    name[SoPlexBase<R>::EPSILON_ZERO] = "epsilon_zero";
    description[SoPlexBase<R>::EPSILON_ZERO] = "general zero tolerance";
    lower[SoPlexBase<R>::EPSILON_ZERO] = 0.0;
    upper[SoPlexBase<R>::EPSILON_ZERO] = 1.0;
    defaultValue[SoPlexBase<R>::EPSILON_ZERO] = DEFAULT_EPS_ZERO;

    ///@todo define suitable values depending on Real type
    // zero tolerance used in factorization
    name[SoPlexBase<R>::EPSILON_FACTORIZATION] = "epsilon_factorization";
    description[SoPlexBase<R>::EPSILON_FACTORIZATION] = "zero tolerance used in factorization";
    lower[SoPlexBase<R>::EPSILON_FACTORIZATION] = 0.0;
    upper[SoPlexBase<R>::EPSILON_FACTORIZATION] = 1.0;
    defaultValue[SoPlexBase<R>::EPSILON_FACTORIZATION] = DEFAULT_EPS_FACTOR;

    ///@todo define suitable values depending on Real type
    // zero tolerance used in update of the factorization
    name[SoPlexBase<R>::EPSILON_UPDATE] = "epsilon_update";
    description[SoPlexBase<R>::EPSILON_UPDATE] = "zero tolerance used in update of the factorization";
    lower[SoPlexBase<R>::EPSILON_UPDATE] = 0.0;
    upper[SoPlexBase<R>::EPSILON_UPDATE] = 1.0;
    defaultValue[SoPlexBase<R>::EPSILON_UPDATE] = DEFAULT_EPS_UPDATE;

    ///@todo define suitable values depending on Real type
    // pivot zero tolerance used in factorization
    name[SoPlexBase<R>::EPSILON_PIVOT] = "epsilon_pivot";
    description[SoPlexBase<R>::EPSILON_PIVOT] = "pivot zero tolerance used in factorization";
    lower[SoPlexBase<R>::EPSILON_PIVOT] = 0.0;
    upper[SoPlexBase<R>::EPSILON_PIVOT] = 1.0;
    defaultValue[SoPlexBase<R>::EPSILON_PIVOT] = DEFAULT_EPS_PIVOT;

    ///@todo define suitable values depending on Real type
    // infinity threshold
    name[SoPlexBase<R>::INFTY] = "infty";
    description[SoPlexBase<R>::INFTY] = "infinity threshold";
    lower[SoPlexBase<R>::INFTY] = 1e10;
    upper[SoPlexBase<R>::INFTY] = 1e100;
    defaultValue[SoPlexBase<R>::INFTY] = DEFAULT_INFINITY;

    // time limit in seconds (INFTY if unlimited)
    name[SoPlexBase<R>::TIMELIMIT] = "timelimit";
    description[SoPlexBase<R>::TIMELIMIT] = "time limit in seconds";
    lower[SoPlexBase<R>::TIMELIMIT] = 0.0;
    upper[SoPlexBase<R>::TIMELIMIT] = DEFAULT_INFINITY;
    defaultValue[SoPlexBase<R>::TIMELIMIT] = DEFAULT_INFINITY;

    // lower limit on objective value
    name[SoPlexBase<R>::OBJLIMIT_LOWER] = "objlimit_lower";
    description[SoPlexBase<R>::OBJLIMIT_LOWER] = "lower limit on objective value";
    lower[SoPlexBase<R>::OBJLIMIT_LOWER] = -DEFAULT_INFINITY;
    upper[SoPlexBase<R>::OBJLIMIT_LOWER] = DEFAULT_INFINITY;
    defaultValue[SoPlexBase<R>::OBJLIMIT_LOWER] = -DEFAULT_INFINITY;

    // upper limit on objective value
    name[SoPlexBase<R>::OBJLIMIT_UPPER] = "objlimit_upper";
    description[SoPlexBase<R>::OBJLIMIT_UPPER] = "upper limit on objective value";
    lower[SoPlexBase<R>::OBJLIMIT_UPPER] = -DEFAULT_INFINITY;
    upper[SoPlexBase<R>::OBJLIMIT_UPPER] = DEFAULT_INFINITY;
    defaultValue[SoPlexBase<R>::OBJLIMIT_UPPER] = DEFAULT_INFINITY;

    // working tolerance for feasibility in floating-point solver during iterative refinement
    name[SoPlexBase<R>::FPFEASTOL] = "fpfeastol";
    description[SoPlexBase<R>::FPFEASTOL] = "working tolerance for feasibility in floating-point solver during iterative refinement";
    lower[SoPlexBase<R>::FPFEASTOL] = 1e-12;
    upper[SoPlexBase<R>::FPFEASTOL] = 1.0;
    defaultValue[SoPlexBase<R>::FPFEASTOL] = 1e-9;

    // working tolerance for optimality in floating-point solver during iterative refinement
    name[SoPlexBase<R>::FPOPTTOL] = "fpopttol";
    description[SoPlexBase<R>::FPOPTTOL] = "working tolerance for optimality in floating-point solver during iterative refinement";
    lower[SoPlexBase<R>::FPOPTTOL] = 1e-12;
    upper[SoPlexBase<R>::FPOPTTOL] = 1.0;
    defaultValue[SoPlexBase<R>::FPOPTTOL] = 1e-9;

    // maximum increase of scaling factors between refinements
    name[SoPlexBase<R>::MAXSCALEINCR] = "maxscaleincr";
    description[SoPlexBase<R>::MAXSCALEINCR] = "maximum increase of scaling factors between refinements";
    lower[SoPlexBase<R>::MAXSCALEINCR] = 1.0;
    upper[SoPlexBase<R>::MAXSCALEINCR] = DEFAULT_INFINITY;
    defaultValue[SoPlexBase<R>::MAXSCALEINCR] = 1e25;

    // lower threshold in lifting (nonzero matrix coefficients with smaller absolute value will be reformulated)
    name[SoPlexBase<R>::LIFTMINVAL] = "liftminval";
    description[SoPlexBase<R>::LIFTMINVAL] = "lower threshold in lifting (nonzero matrix coefficients with smaller absolute value will be reformulated)";
    lower[SoPlexBase<R>::LIFTMINVAL] = 0.0;
    upper[SoPlexBase<R>::LIFTMINVAL] = 0.1;
    defaultValue[SoPlexBase<R>::LIFTMINVAL] = 0.000976562; // = 1/1024

    // upper threshold in lifting (nonzero matrix coefficients with larger absolute value will be reformulated)
    name[SoPlexBase<R>::LIFTMAXVAL] = "liftmaxval";
    description[SoPlexBase<R>::LIFTMAXVAL] = "lower threshold in lifting (nonzero matrix coefficients with smaller absolute value will be reformulated)";
    lower[SoPlexBase<R>::LIFTMAXVAL] = 10.0;
    upper[SoPlexBase<R>::LIFTMAXVAL] = DEFAULT_INFINITY;
    defaultValue[SoPlexBase<R>::LIFTMAXVAL] = 1024.0;

    // threshold for using sparse pricing (no. of violations need to be smaller than threshold * dimension of problem)
    name[SoPlexBase<R>::SPARSITY_THRESHOLD] = "sparsity_threshold";
    description[SoPlexBase<R>::SPARSITY_THRESHOLD] = "sparse pricing threshold (#violations < dimension * SPARSITY_THRESHOLD activates sparse pricing)";
    lower[SoPlexBase<R>::SPARSITY_THRESHOLD] = 0.0;
    upper[SoPlexBase<R>::SPARSITY_THRESHOLD] = 1.0;
    defaultValue[SoPlexBase<R>::SPARSITY_THRESHOLD] = 0.6;

    // threshold on number of rows vs. number of columns for switching from column to row representations in auto mode
    name[SoPlexBase<R>::REPRESENTATION_SWITCH] = "representation_switch";
    description[SoPlexBase<R>::REPRESENTATION_SWITCH] = "threshold on number of rows vs. number of columns for switching from column to row representations in auto mode";
    lower[SoPlexBase<R>::REPRESENTATION_SWITCH] = 0.0;
    upper[SoPlexBase<R>::REPRESENTATION_SWITCH] = DEFAULT_INFINITY;
    defaultValue[SoPlexBase<R>::REPRESENTATION_SWITCH] = 1.2;

    // geometric frequency at which to apply rational reconstruction
    name[SoPlexBase<R>::RATREC_FREQ] = "ratrec_freq";
    description[SoPlexBase<R>::RATREC_FREQ] = "geometric frequency at which to apply rational reconstruction";
    lower[SoPlexBase<R>::RATREC_FREQ] = 1.0;
    upper[SoPlexBase<R>::RATREC_FREQ] = DEFAULT_INFINITY;
    defaultValue[SoPlexBase<R>::RATREC_FREQ] = 1.2;

    // minimal reduction (sum of removed rows/cols) to continue simplification
    name[SoPlexBase<R>::MINRED] = "minred";
    description[SoPlexBase<R>::MINRED] = "minimal reduction (sum of removed rows/cols) to continue simplification";
    lower[SoPlexBase<R>::MINRED] = 0.0;
    upper[SoPlexBase<R>::MINRED] = 1.0;
    defaultValue[SoPlexBase<R>::MINRED] = 1e-4;

    // refactor threshold for nonzeros in last factorized basis matrix compared to updated basis matrix
    name[SoPlexBase<R>::REFAC_BASIS_NNZ] = "refac_basis_nnz";
    description[SoPlexBase<R>::REFAC_BASIS_NNZ] = "refactor threshold for nonzeros in last factorized basis matrix compared to updated basis matrix";
    lower[SoPlexBase<R>::REFAC_BASIS_NNZ] = 1.0;
    upper[SoPlexBase<R>::REFAC_BASIS_NNZ] = 100.0;
    defaultValue[SoPlexBase<R>::REFAC_BASIS_NNZ] = 10.0;

    // refactor threshold for fill-in in current factor update compared to fill-in in last factorization
    name[SoPlexBase<R>::REFAC_UPDATE_FILL] = "refac_update_fill";
    description[SoPlexBase<R>::REFAC_UPDATE_FILL] = "refactor threshold for fill-in in current factor update compared to fill-in in last factorization";
    lower[SoPlexBase<R>::REFAC_UPDATE_FILL] = 1.0;
    upper[SoPlexBase<R>::REFAC_UPDATE_FILL] = 100.0;
    defaultValue[SoPlexBase<R>::REFAC_UPDATE_FILL] = 5.0;

    // refactor threshold for memory growth in factorization since last refactorization
    name[SoPlexBase<R>::REFAC_MEM_FACTOR] = "refac_mem_factor";
    description[SoPlexBase<R>::REFAC_MEM_FACTOR] = "refactor threshold for memory growth in factorization since last refactorization";
    lower[SoPlexBase<R>::REFAC_MEM_FACTOR] = 1.0;
    upper[SoPlexBase<R>::REFAC_MEM_FACTOR] = 10.0;
    defaultValue[SoPlexBase<R>::REFAC_MEM_FACTOR] = 1.5;

    // accuracy of conjugate gradient method in least squares scaling (higher value leads to more iterations)
    name[SoPlexBase<R>::LEASTSQ_ACRCY] = "leastsq_acrcy";
    description[SoPlexBase<R>::LEASTSQ_ACRCY] = "accuracy of conjugate gradient method in least squares scaling (higher value leads to more iterations)";
    lower[SoPlexBase<R>::LEASTSQ_ACRCY] = 1.0;
    upper[SoPlexBase<R>::LEASTSQ_ACRCY] = DEFAULT_INFINITY;
    defaultValue[SoPlexBase<R>::LEASTSQ_ACRCY] = 1000.0;

    // objective offset
    name[SoPlexBase<R>::OBJ_OFFSET] = "obj_offset";
    description[SoPlexBase<R>::OBJ_OFFSET] = "objective offset to be used";
    lower[SoPlexBase<R>::OBJ_OFFSET] = -DEFAULT_INFINITY;
    upper[SoPlexBase<R>::OBJ_OFFSET] = DEFAULT_INFINITY;
    defaultValue[SoPlexBase<R>::OBJ_OFFSET] = 0.0;
  }

  template <>
  typename SoPlexBase<Real>::Settings& SoPlexBase<Real>::Settings::operator=(const Settings& settings)
  {
    for( int i = 0; i < SoPlexBase<Real>::BOOLPARAM_COUNT; i++ )
      _boolParamValues[i] = settings._boolParamValues[i];

    for( int i = 0; i < SoPlexBase<Real>::INTPARAM_COUNT; i++ )
      _intParamValues[i] = settings._intParamValues[i];

    for( int i = 0; i < SoPlexBase<Real>::REALPARAM_COUNT; i++ )
      _realParamValues[i] = settings._realParamValues[i];

#ifdef SOPLEX_WITH_RATIONALPARAM
    for( int i = 0; i < SoPlexBase<Real>::RATIONALPARAM_COUNT; i++ )
      _rationalParamValues[i] = settings._rationalParamValues[i];
#endif

    return *this;
  }


#ifdef SOPLEX_WITH_RATIONALPARAM
  SoPlexBase<R>::Settings::RationalParam::RationalParam() {}
#endif


  template <>
  SoPlexBase<Real>::Settings::Settings()
  {
    for( int i = 0; i < SoPlexBase<Real>::BOOLPARAM_COUNT; i++ )
      _boolParamValues[i] = boolParam.defaultValue[i];

    for( int i = 0; i < SoPlexBase<Real>::INTPARAM_COUNT; i++ )
      _intParamValues[i] = intParam.defaultValue[i];

    for( int i = 0; i < SoPlexBase<Real>::REALPARAM_COUNT; i++ )
      _realParamValues[i] = realParam.defaultValue[i];

#ifdef SOPLEX_WITH_RATIONALPARAM
    for( int i = 0; i < SoPlexBase<Real>::RATIONALPARAM_COUNT; i++ )
      _rationalParamValues[i] = rationalParam.defaultValue[i];
#endif
  }

  template <>
  SoPlexBase<Real>::Settings::Settings(const Settings& settings)
  {
    *this = settings;
  }


#ifdef SOPLEX_WITH_RATIONALPARAM
  template <>
  SoPlexBase<Real>::Settings::RationalParam SoPlexBase<Real>::Settings::rationalParam;
#endif


  /// copy constructor
  ///@todo improve performance by implementing a separate copy constructor
  template <>
  SoPlexBase<Real>::SoPlexBase(const SoPlexBase<Real>& rhs)
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
  SoPlexBase<Real>::~SoPlexBase()
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

  // For SCIP compatibility
  template <>
  int SoPlexBase<Real>::numRowsReal() const
  {
    return numRows();
  }


  // Wrapper for reverse compatibility
  template <>
  int SoPlexBase<Real>::numColsReal() const
  {
    return numCols();
  }



  /// solves the LP
  /// Real specialization of the optimize function
  template <>
  typename SPxSolverBase<Real>::Status SoPlexBase<Real>::optimize()
  {
    assert(_isConsistent());

    // clear statistics
    _statistics->clearSolvingData();

    // the solution is no longer valid
    _invalidateSolution();

    // if the decomposition based dual simplex flag is set to true
    if ( boolParam(SoPlexBase<Real>::USEDECOMPDUALSIMPLEX) )
      {
        setIntParam(SoPlexBase<Real>::SOLVEMODE, SOLVEMODE_REAL);
        setIntParam(SoPlexBase<Real>::REPRESENTATION, REPRESENTATION_ROW);
        setIntParam(SoPlexBase<Real>::ALGORITHM, ALGORITHM_DUAL);
        //setBoolParam(SoPlexBase<Real>::PERSISTENTSCALING, false);

        _solver.setComputeDegenFlag(boolParam(COMPUTEDEGEN));

        _solveDecompositionDualSimplex();
      }
    // decide whether to solve the rational LP with iterative refinement or call the standard floating-point solver
    else if( intParam(SoPlexBase<Real>::SOLVEMODE) == SOLVEMODE_REAL || (intParam(SoPlexBase<Real>::SOLVEMODE) == SOLVEMODE_AUTO
                                                                  && GE(realParam(SoPlexBase<Real>::FEASTOL), 1e-9) && GE(realParam(SoPlexBase<Real>::OPTTOL), 1e-9)) )
      {
        // ensure that tolerances are reasonable for the floating-point solver
        if( realParam(SoPlexBase<Real>::FEASTOL) < _currentSettings->realParam.lower[SoPlexBase<Real>::FPFEASTOL] )
          {
            MSG_WARNING( spxout, spxout << "Cannot call floating-point solver with feasibility tolerance below "
                         << _currentSettings->realParam.lower[SoPlexBase<Real>::FPFEASTOL] << " - relaxing tolerance\n");
            _solver.setFeastol(_currentSettings->realParam.lower[SoPlexBase<Real>::FPFEASTOL]);
          }
        else
          _solver.setFeastol(realParam(SoPlexBase<Real>::FEASTOL));

        if( realParam(SoPlexBase<Real>::OPTTOL) < _currentSettings->realParam.lower[SoPlexBase<Real>::FPOPTTOL] )
          {
            MSG_WARNING( spxout, spxout << "Cannot call floating-point solver with optimality tolerance below "
                         << _currentSettings->realParam.lower[SoPlexBase<Real>::FPOPTTOL] << " - relaxing tolerance\n");
            _solver.setOpttol(_currentSettings->realParam.lower[SoPlexBase<Real>::FPOPTTOL]);
          }
        else
          _solver.setOpttol(realParam(SoPlexBase<Real>::OPTTOL));

        _solver.setComputeDegenFlag(boolParam(COMPUTEDEGEN));

        _optimize();
#ifdef SOPLEX_DEBUG // this check will remove scaling of the realLP
        _checkBasisScaling();
#endif
      }
    else if( intParam(SoPlexBase<Real>::SYNCMODE) == SYNCMODE_ONLYREAL )
      {
        _syncLPRational();
        _optimizeRational();
      }
    else if( intParam(SoPlexBase<Real>::SYNCMODE) == SYNCMODE_MANUAL )
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


  // Wrapper function for reverse compatibility
  template <>
	bool SoPlexBase<Real>::getPrimalReal(VectorBase<Real>& vector)
  {
    return getPrimal(vector);
  }


  // Wrapper function for reverse compatibility
  template <>
  bool SoPlexBase<Real>::getPrimalRayReal(VectorBase<Real>& vector)
  {
    return getPrimalRay(vector);
  }


  // Wrapper function for reverse compatibility
  template <>
	bool SoPlexBase<Real>::getDualReal(VectorBase<Real>& vector) // For SCIP
  {
    return getDual(vector);
  }



  /// wrapper for backwards compatibility
  template <>
	bool SoPlexBase<Real>::getRedCostReal(VectorBase<Real>& vector) // For SCIP compatibility
  {
    return getRedCost(vector);
  }




  // Wrapping the function for reverse Compatibility
  template <>
	bool SoPlexBase<Real>::getDualFarkasReal(VectorBase<Real>& vector)
  {
    return getDualFarkas(vector);
  }



#ifdef SOPLEX_WITH_GMP
  /// gets the primal solution vector if available; returns true on success
  template <>
	bool SoPlexBase<Real>::getPrimalRational(mpq_t* vector, const int size)
  {
    assert(size >= numColsRational());

    if( hasPrimal() )
      {
        _syncRationalSolution();
        for( int i = 0; i < numColsRational(); i++ )
          mpq_set(vector[i], _solRational._primal[i].getMpqRef());
        return true;
      }
    else
      return false;
  }


  /// gets the vector of slack values if available; returns true on success
  template <>
	bool SoPlexBase<Real>::getSlacksRational(mpq_t* vector, const int size)
  {
    assert(size >= numRowsRational());

    if( hasPrimal() )
      {
        _syncRationalSolution();
        for( int i = 0; i < numRowsRational(); i++ )
          mpq_set(vector[i], _solRational._slacks[i].getMpqRef());
        return true;
      }
    else
      return false;
  }



  /// gets the primal ray if LP is unbounded; returns true on success
  template <>
	bool SoPlexBase<Real>::getPrimalRayRational(mpq_t* vector, const int size)
  {
    assert(size >= numColsRational());

    if( hasPrimalRay() )
      {
        _syncRationalSolution();
        for( int i = 0; i < numColsRational(); i++ )
          mpq_set(vector[i], _solRational._primalRay[i].getMpqRef());
        return true;
      }
    else
      return false;
  }



  /// gets the dual solution vector if available; returns true on success
  template <>
	bool SoPlexBase<Real>::getDualRational(mpq_t* vector, const int size)
  {
    assert(size >= numRowsRational());

    if( hasDual() )
      {
        _syncRationalSolution();
        for( int i = 0; i < numRowsRational(); i++ )
          mpq_set(vector[i], _solRational._dual[i].getMpqRef());
        return true;
      }
    else
      return false;
  }



  /// gets the vector of reduced cost values if available; returns true on success
  template <>
	bool SoPlexBase<Real>::getRedCostRational(mpq_t* vector, const int size)
  {
    assert(size >= numColsRational());

    if( hasDual() )
      {
        _syncRationalSolution();
        for( int i = 0; i < numColsRational(); i++ )
          mpq_set(vector[i], _solRational._redCost[i].getMpqRef());
        return true;
      }
    else
      return false;
  }



  /// gets the Farkas proof if LP is infeasible; returns true on success
  template <>
	bool SoPlexBase<Real>::getDualFarkasRational(mpq_t* vector, const int size)
  {
    assert(size >= numRowsRational());

    if( hasDualFarkas() )
      {
        _syncRationalSolution();
        for( int i = 0; i < numRowsRational(); i++ )
          mpq_set(vector[i], _solRational._dualFarkas[i].getMpqRef());
        return true;
      }
    else
      return false;
  }
#endif






  /// reads LP file in LP or MPS format according to READMODE parameter; gets row names, column names, and
  /// integer variables if desired; returns true on success

  template <>
  bool SoPlexBase<Real>::readFile(const char* filename, NameSet* rowNames, NameSet* colNames, DIdxSet* intVars)
  {
    bool success = false;
    if( intParam(SoPlexBase<Real>::READMODE) == READMODE_REAL )
      success = _readFileReal(filename, rowNames, colNames, intVars);
    else
      success = _readFileRational(filename, rowNames, colNames, intVars);

    // storing the row and column names for use in the DBDS print basis methods
    _rowNames = rowNames;
    _colNames = colNames;

    return success;
  }

  /// writes real LP to file; LP or MPS format is chosen from the extension in \p filename; if \p rowNames and \p
  /// colNames are \c NULL, default names are used; if \p intVars is not \c NULL, the variables contained in it are
  /// marked as integer; returns true on success
  template <>
	bool SoPlexBase<Real>::writeFile(const char* filename, const NameSet* rowNames, const NameSet* colNames, const DIdxSet* intVars, const bool unscale) const
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
        origLP->writeFileLPBase(filename, rowNames, colNames, intVars);
        origLP->~SPxLPReal();
        spx_free(origLP);
      }
    else
      _realLP->writeFileLPBase(filename, rowNames, colNames, intVars);

    return true;
  }

  // Alias for writeFile; SCIP
  template <>
	bool SoPlexBase<Real>::writeFileReal(const char* filename, const NameSet* rowNames, const NameSet* colNames, const DIdxSet* intVars, const bool unscale) const
  {
    return writeFile(filename, rowNames, colNames, intVars, unscale);
  }

  /// writes rational LP to file; LP or MPS format is chosen from the extension in \p filename; if \p rowNames and \p
  /// colNames are \c NULL, default names are used; if \p intVars is not \c NULL, the variables contained in it are
  /// marked as integer; returns true on success
  /// Here unscale is just a junk variable that is used to match the type with the real write function
  template <>
	bool SoPlexBase<Real>::writeFileRational(const char* filename, const NameSet* rowNames, const NameSet* colNames, const DIdxSet* intVars) const
  {
    if( intParam(SoPlexBase<Real>::SYNCMODE) == SYNCMODE_ONLYREAL )
      return false;
    else
      {
        assert(_rationalLP != 0);
        _rationalLP->writeFileLPBase(filename, rowNames, colNames, intVars);

        ///@todo implement return value
        return true;
      }
  }



  /// writes the dual of the real LP to file; LP or MPS format is chosen from the extension in \p filename;
  /// if \p rowNames and \p colNames are \c NULL, default names are used; if \p intVars is not \c NULL,
  /// the variables contained in it are marked as integer; returns true on success
  template <>
	bool SoPlexBase<Real>::writeDualFileReal(const char* filename, const NameSet* rowNames, const NameSet* colNames, const DIdxSet* intVars) const
  {
    SPxLPReal dualLP;
    _realLP->buildDualProblem(dualLP);
    dualLP.setOutstream(spxout);

    // swap colnames and rownames
    dualLP.writeFileLPBase(filename, colNames, rowNames);
    return true;
  }



  /// reads basis information from \p filename and returns true on success; if \p rowNames and \p colNames are \c NULL,
  /// default names are assumed; returns true on success
  template <>
  bool SoPlexBase<Real>::readBasisFile(const char* filename, const NameSet* rowNames, const NameSet* colNames)
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
    assert(_hasBasis == (_solver.basis().status() > SPxBasisBase<Real>::NO_PROBLEM));

    // stop timing
    _statistics->readingTime->stop();

    return _hasBasis;
#else
    // this is alternative code for reading bases without the SPxSolverBase class
    assert(filename != 0);

    // start timing
    _statistics->readingTime->start();

    // read
    spxifstream file(filename);

    if( !file )
      return false;

    // get problem size
    int numRows = numRows();
    int numCols = numCols();

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
      _basisStatusRows[i] = SPxSolverBase<Real>::BASIC;

    for( int i = 0; i < numCols; i++ )
      {
        if( lowerRealInternal(i) == upperRealInternal(i) )
          _basisStatusCols[i] = SPxSolverBase<Real>::FIXED;
        else if( lowerRealInternal(i) <= double(-realParam(SoPlexBase<Real>::INFTY)) && upperRealInternal(i) >= double(realParam(SoPlexBase<Real>::INFTY)) )
          _basisStatusCols[i] = SPxSolverBase<Real>::ZERO;
        else if( lowerRealInternal(i) <= double(-realParam(SoPlexBase<Real>::INFTY)) )
          _basisStatusCols[i] = SPxSolverBase<Real>::ON_UPPER;
        else
          _basisStatusCols[i] = SPxSolverBase<Real>::ON_LOWER;
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
                _basisStatusCols[c] = SPxSolverBase<Real>::BASIC;
                if( _rowTypes[r] == SoPlexBase<Real>::RANGETYPE_LOWER )
                  _basisStatusRows[r] = SPxSolverBase<Real>::ON_LOWER;
                else if( _rowTypes[r] == SoPlexBase<Real>::RANGETYPE_FIXED )
                  _basisStatusRows[r] = SPxSolverBase<Real>::FIXED;
                else
                  _basisStatusRows[r] = SPxSolverBase<Real>::ON_UPPER;
              }
            else if( !strcmp(mps.field1(), "XL") )
              {
                _basisStatusCols[c] = SPxSolverBase<Real>::BASIC;
                if( _rowTypes[r] == SoPlexBase<Real>::RANGETYPE_UPPER )
                  _basisStatusRows[r] = SPxSolverBase<Real>::ON_UPPER;
                else if( _rowTypes[r] == SoPlexBase<Real>::RANGETYPE_FIXED )
                  _basisStatusRows[r] = SPxSolverBase<Real>::FIXED;
                else
                  _basisStatusRows[r] = SPxSolverBase<Real>::ON_LOWER;
              }
            else if( !strcmp(mps.field1(), "UL") )
              {
                _basisStatusCols[c] = SPxSolverBase<Real>::ON_UPPER;
              }
            else if( !strcmp(mps.field1(), "LL") )
              {
                _basisStatusCols[c] = SPxSolverBase<Real>::ON_LOWER;
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
  template <>
	bool SoPlexBase<Real>::writeBasisFile(const char* filename, const NameSet* rowNames, const NameSet* colNames, const bool cpxFormat) const
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
            assert(_basisStatusCols[col] != SPxSolverBase<Real>::UNDEFINED);

            if( _basisStatusCols[col] == SPxSolverBase<Real>::BASIC )
              {
                // find nonbasic row
                for( ; row < numRows; row++ )
                  {
                    assert(_basisStatusRows[row] != SPxSolverBase<Real>::UNDEFINED);
                    if( _basisStatusRows[row] != SPxSolverBase<Real>::BASIC )
                      break;
                  }

                assert(row != numRows);

                if( _basisStatusRows[row] == SPxSolverBase<Real>::ON_UPPER && (!cpxFormat || _rowTypes[row] == SoPlexBase<Real>::RANGETYPE_BOXED) )
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
                if( _basisStatusCols[col] == SPxSolverBase<Real>::ON_UPPER )
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
            assert(_basisStatusRows[row] == SPxSolverBase<Real>::BASIC);
          }
#endif

        return true;
      }
  }



  /// writes internal LP, basis information, and parameter settings; if \p rowNames and \p colNames are \c NULL,
  /// default names are used
  template <>
  void SoPlexBase<Real>::writeStateReal(const char* filename, const NameSet* rowNames, const NameSet* colNames, const bool cpxFormat) const
  {
    std::string ofname;

    // write parameter settings
    ofname = std::string(filename) + ".set";
    saveSettingsFile(ofname.c_str());

    // write problem in MPS/LP format
    ofname = std::string(filename) + ((cpxFormat) ? ".lp" : ".mps");
    writeFile(ofname.c_str(), rowNames, colNames, 0);

    // write basis
    ofname = std::string(filename) + ".bas";
    writeBasisFile(ofname.c_str(), rowNames, colNames, cpxFormat);
  }



  /// writes internal LP, basis information, and parameter settings; if \p rowNames and \p colNames are \c NULL,
  /// default names are used
  template <>
  void SoPlexBase<Real>::writeStateRational(const char* filename, const NameSet* rowNames, const NameSet* colNames, const bool cpxFormat) const
      {
    std::string ofname;

    // write parameter settings
    ofname = std::string(filename) + ".set";
    saveSettingsFile(ofname.c_str());

    // write problem in MPS/LP format
    ofname = std::string(filename) + ((cpxFormat) ? ".lp" : ".mps");
    writeFileRational(ofname.c_str(), rowNames, colNames, 0);

    // write basis
    ofname = std::string(filename) + ".bas";
    writeBasisFile(ofname.c_str(), rowNames, colNames, cpxFormat);
  }




  /// prints solution statistics
  template <>
  void SoPlexBase<Real>::printSolutionStatistics(std::ostream& os)
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

        os << "Violation (rational): \n";
        if( getBoundViolationRational(maxviol, sumviol) )
          os << "  Max/sum bound     : " << rationalToString(maxviol) << " / " << rationalToString(sumviol) << "\n";
        else
          os << "  Max/sum bound     : - / -\n";
        if( getRowViolationRational(maxviol, sumviol) )
          os << "  Max/sum row       : " << rationalToString(maxviol) << " / " << rationalToString(sumviol) << "\n";
        else
          os << "  Max/sum row       : - / -\n";
        if( getRedCostViolationRational(maxviol, sumviol) )
          os << "  Max/sum redcost   : " << rationalToString(maxviol) << " / " << rationalToString(sumviol) << "\n";
        else
          os << "  Max/sum redcost   : - / -\n";
        if( getDualViolationRational(maxviol, sumviol) )
          os << "  Max/sum dual      : " << rationalToString(maxviol) << " / " << rationalToString(sumviol) << "\n";
        else
          os << "  Max/sum dual      : - / -\n";
      }
    else
      {
        Real maxviol;
        Real sumviol;

        os << "Violations (real)   : \n";
        if( getBoundViolation(maxviol, sumviol) )
          os << "  Max/sum bound     : " << maxviol << " / " << sumviol << "\n";
        else
          os << "  Max/sum bound     : - / -\n";
        if( getRowViolation(maxviol, sumviol) )
          os << "  Max/sum row       : " << maxviol << " / " << sumviol << "\n";
        else
          os << "  Max/sum row       : - / -\n";
        if( getRedCostViolation(maxviol, sumviol) )
          os << "  Max/sum redcost   : " << maxviol << " / " << sumviol << "\n";
        else
          os << "  Max/sum redcost   : - / -\n";
        if( getDualViolation(maxviol, sumviol) )
          os << "  Max/sum dual      : " << maxviol << " / " << sumviol << "\n";
        else
          os << "  Max/sum dual      : - / -\n";
      }
  }


  /// prints statistics on solving process
  template <>
  void SoPlexBase<Real>::printSolvingStatistics(std::ostream& os)
  {
    assert(_statistics != 0);
    _statistics->print(os);
  }



  /// prints short statistics
  template <>
  void SoPlexBase<Real>::printShortStatistics(std::ostream& os)
  {
    printStatus(os, _status);
    SPxOut::setFixed(os, 2);
    os << "Solving time (sec)  : " << this->_statistics->solvingTime->time() << "\n"
       << "Iterations          : " << this->_statistics->iterations << "\n";
    SPxOut::setScientific(os);
    os << "Objective value     : " << objValueReal() << "\n";
  }



  /// prints complete statistics
  template <>
  void SoPlexBase<Real>::printStatistics(std::ostream& os)
  {
    SPxOut::setFixed(os, 2);

    printStatus(os, _status);

    os << "Original problem    : \n";
    if ( boolParam(SoPlexBase<Real>::USEDECOMPDUALSIMPLEX) )
      printOriginalProblemStatistics(os);
    else
      {
        if( intParam(SoPlexBase<Real>::READMODE) == READMODE_REAL )
          _realLP->printProblemStatistics(os);
        else
          _rationalLP->printProblemStatistics(os);
      }

    os << "Objective sense     : " << (intParam(SoPlexBase<Real>::OBJSENSE) == SoPlexBase<Real>::OBJSENSE_MINIMIZE ? "minimize\n" : "maximize\n");
    printSolutionStatistics(os);
    printSolvingStatistics(os);
  }

} // namespace soplex
