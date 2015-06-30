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
#include "spxfileio.h"
#include "spxgithash.h"
#include "statistics.h"
#include "mpsinput.h"

/// maximum length of lines in settings file
#define SET_MAX_LINE_LEN 500

#ifdef _MSC_VER
#define strncasecmp _strnicmp
#endif

namespace soplex
{
   /// class of parameter settings
   class SoPlex::Settings
   {
   public:
      /// array of names for boolean parameters
      static std::string _boolParamName[SoPlex::BOOLPARAM_COUNT];

      /// array of names for integer parameters
      static std::string _intParamName[SoPlex::INTPARAM_COUNT];

      /// array of names for real parameters
      static std::string _realParamName[SoPlex::REALPARAM_COUNT];

#ifdef SOPLEX_WITH_RATIONALPARAM
      /// array of names for rational parameters
      static std::string _rationalParamName[SoPlex::RATIONALPARAM_COUNT];
#endif

      /// array of descriptions for boolean parameters
      static std::string _boolParamDescription[SoPlex::BOOLPARAM_COUNT];

      /// array of descriptions for integer parameters
      static std::string _intParamDescription[SoPlex::INTPARAM_COUNT];

      /// array of descriptions for real parameters
      static std::string _realParamDescription[SoPlex::REALPARAM_COUNT];

#ifdef SOPLEX_WITH_RATIONALPARAM
      /// array of descriptions for rational parameters
      static std::string _rationalParamDescription[SoPlex::RATIONALPARAM_COUNT];
#endif

      /// array of default values for boolean parameters
      static bool _boolParamDefault[SoPlex::BOOLPARAM_COUNT];

      /// array of default values for integer parameters
      static int _intParamDefault[SoPlex::INTPARAM_COUNT];

      /// array of default values for real parameters
      static Real _realParamDefault[SoPlex::REALPARAM_COUNT];

#ifdef SOPLEX_WITH_RATIONALPARAM
      /// array of default values for rational parameters
      static Rational _rationalParamDefault[SoPlex::RATIONALPARAM_COUNT];
#endif

      /// array of lower bounds for int parameter values
      static int _intParamLower[SoPlex::INTPARAM_COUNT];

      /// array of upper bounds for int parameter values
      static int _intParamUpper[SoPlex::INTPARAM_COUNT];

      /// array of lower bounds for real parameter values
      static Real _realParamLower[SoPlex::REALPARAM_COUNT];

      /// array of upper bounds for real parameter values
      static Real _realParamUpper[SoPlex::REALPARAM_COUNT];

#ifdef SOPLEX_WITH_RATIONALPARAM
      /// array of lower bounds for rational parameter values
      static Rational _rationalParamLower[SoPlex::RATIONALPARAM_COUNT];

      /// array of upper bounds for rational parameter values
      static Rational _rationalParamUpper[SoPlex::RATIONALPARAM_COUNT];
#endif

      /// have static arrays been initialized?
      static bool _defaultsAndBoundsInitialized;

      /// array of current boolean parameter values
      bool _boolParamValues[SoPlex::BOOLPARAM_COUNT];

      /// array of current integer parameter values
      int _intParamValues[SoPlex::INTPARAM_COUNT];

      /// array of current real parameter values
      Real _realParamValues[SoPlex::REALPARAM_COUNT];

#ifdef SOPLEX_WITH_RATIONALPARAM
      /// array of current rational parameter values
      Rational _rationalParamValues[SoPlex::RATIONALPARAM_COUNT];
#endif

      /// default constructor initializing default settings
      Settings()
      {
         if( !_defaultsAndBoundsInitialized )
         {
            // should lifting be used to reduce range of nonzero matrix coefficients?
            _boolParamName[SoPlex::LIFTING] = "lifting";
            _boolParamDescription[SoPlex::LIFTING] = "should lifting be used to reduce range of nonzero matrix coefficients?";
            _boolParamDefault[SoPlex::LIFTING] = false;

            // should LP be transformed to equality form before a rational solve?
            _boolParamName[SoPlex::EQTRANS] = "eqtrans";
            _boolParamDescription[SoPlex::EQTRANS] = "should LP be transformed to equality form before a rational solve?";
            _boolParamDefault[SoPlex::EQTRANS] = false;

            // should dual infeasibility be tested in order to try to return a dual solution even if primal infeasible?
            _boolParamName[SoPlex::TESTDUALINF] = "testdualinf";
            _boolParamDescription[SoPlex::TESTDUALINF] = "should dual infeasibility be tested in order to try to return a dual solution even if primal infeasible?";
            _boolParamDefault[SoPlex::TESTDUALINF] = false;

            // should a rational factorization be performed after iterative refinement?
            _boolParamName[SoPlex::RATFAC] = "ratfac";
            _boolParamDescription[SoPlex::RATFAC] = "should a rational factorization be performed after iterative refinement?";
            _boolParamDefault[SoPlex::RATFAC] = true;

            // should cycling solutions be accepted during iterative refinement?
            _boolParamName[SoPlex::ACCEPTCYCLING] = "acceptcycling";
            _boolParamDescription[SoPlex::ACCEPTCYCLING] = "should cycling solutions be accepted during iterative refinement?";
            _boolParamDefault[SoPlex::ACCEPTCYCLING] = false;

            // apply rational reconstruction after each iterative refinement?
            _boolParamName[SoPlex::RATREC] = "ratrec";
            _boolParamDescription[SoPlex::RATREC] = "apply rational reconstruction after each iterative refinement?";
            _boolParamDefault[SoPlex::RATREC] = true;

            // round scaling factors for iterative refinement to powers of two?
            _boolParamName[SoPlex::POWERSCALING] = "powerscaling";
            _boolParamDescription[SoPlex::POWERSCALING] = "round scaling factors for iterative refinement to powers of two?";
            _boolParamDefault[SoPlex::POWERSCALING] = true;

            // continue iterative refinement with exact basic solution if not optimal?
            _boolParamName[SoPlex::RATFACJUMP] = "ratfacjump";
            _boolParamDescription[SoPlex::RATFACJUMP] = "continue iterative refinement with exact basic solution if not optimal?";
            _boolParamDefault[SoPlex::RATFACJUMP] = false;

            // should feasibility be tested with relaxed bounds and sides?
            _boolParamName[SoPlex::FEASRELAX] = "feasrelax";
            _boolParamDescription[SoPlex::FEASRELAX] = "should feasibility be tested with relaxed bounds and sides?";
            _boolParamDefault[SoPlex::FEASRELAX] = false;

            // use bound flipping also for row representation?
            _boolParamName[SoPlex::ROWBOUNDFLIPS] = "rowboundflips";
            _boolParamDescription[SoPlex::ROWBOUNDFLIPS] = "use bound flipping also for row representation?";
            _boolParamDefault[SoPlex::ROWBOUNDFLIPS] = false;

            // objective sense
            _intParamName[SoPlex::OBJSENSE] = "objsense";
            _intParamDescription[SoPlex::OBJSENSE] = "objective sense (-1 - minimize, +1 - maximize)";
            _intParamLower[SoPlex::OBJSENSE] = -1;
            _intParamUpper[SoPlex::OBJSENSE] = 1;
            _intParamDefault[SoPlex::OBJSENSE] = SoPlex::OBJSENSE_MAXIMIZE;

            // type of computational form, i.e., column or row representation
            _intParamName[SoPlex::REPRESENTATION] = "representation";
            _intParamDescription[SoPlex::REPRESENTATION] = "type of computational form (0 - auto, 1 - column representation, 2 - row representation)";
            _intParamLower[SoPlex::REPRESENTATION] = 0;
            _intParamUpper[SoPlex::REPRESENTATION] = 2;
            _intParamDefault[SoPlex::REPRESENTATION] = SoPlex::REPRESENTATION_AUTO;

            // type of algorithm, i.e., primal or dual
            _intParamName[SoPlex::ALGORITHM] = "algorithm";
            _intParamDescription[SoPlex::ALGORITHM] = "type of algorithm (0 - primal, 1 - dual)";
            _intParamLower[SoPlex::ALGORITHM] = 0;
            _intParamUpper[SoPlex::ALGORITHM] = 1;
            _intParamDefault[SoPlex::ALGORITHM] = SoPlex::ALGORITHM_DUAL;

            // type of LU update
            _intParamName[SoPlex::FACTOR_UPDATE_TYPE] = "factor_update_type";
            _intParamDescription[SoPlex::FACTOR_UPDATE_TYPE] = "type of LU update (0 - eta update, 1 - Forrest-Tomlin update)";
            _intParamLower[SoPlex::FACTOR_UPDATE_TYPE] = 0;
            _intParamUpper[SoPlex::FACTOR_UPDATE_TYPE] = 1;
            _intParamDefault[SoPlex::FACTOR_UPDATE_TYPE] = SoPlex::FACTOR_UPDATE_TYPE_FT;

            ///@todo which value?
            // maximum number of updates without fresh factorization
            _intParamName[SoPlex::FACTOR_UPDATE_MAX] = "factor_update_max";
            _intParamDescription[SoPlex::FACTOR_UPDATE_MAX] = "maximum number of LU updates without fresh factorization";
            _intParamLower[SoPlex::FACTOR_UPDATE_MAX] = 0;
            _intParamUpper[SoPlex::FACTOR_UPDATE_MAX] = INT_MAX;
            _intParamDefault[SoPlex::FACTOR_UPDATE_MAX] = 200;

            // iteration limit (-1 if unlimited)
            _intParamName[SoPlex::ITERLIMIT] = "iterlimit";
            _intParamDescription[SoPlex::ITERLIMIT] = "iteration limit (-1 - no limit)";
            _intParamLower[SoPlex::ITERLIMIT] = -1;
            _intParamUpper[SoPlex::ITERLIMIT] = INT_MAX;
            _intParamDefault[SoPlex::ITERLIMIT] = -1;

            // refinement limit (-1 if unlimited)
            _intParamName[SoPlex::REFLIMIT] = "reflimit";
            _intParamDescription[SoPlex::REFLIMIT] = "refinement limit (-1 - no limit)";
            _intParamLower[SoPlex::REFLIMIT] = -1;
            _intParamUpper[SoPlex::REFLIMIT] = INT_MAX;
            _intParamDefault[SoPlex::REFLIMIT] = -1;

            // stalling refinement limit (-1 if unlimited)
            _intParamName[SoPlex::STALLREFLIMIT] = "stallreflimit";
            _intParamDescription[SoPlex::STALLREFLIMIT] = "stalling refinement limit (-1 - no limit)";
            _intParamLower[SoPlex::STALLREFLIMIT] = -1;
            _intParamUpper[SoPlex::STALLREFLIMIT] = INT_MAX;
            _intParamDefault[SoPlex::STALLREFLIMIT] = -1;

            // display frequency
            _intParamName[SoPlex::DISPLAYFREQ] = "displayfreq";
            _intParamDescription[SoPlex::DISPLAYFREQ] = "display frequency";
            _intParamLower[SoPlex::DISPLAYFREQ] = 1;
            _intParamUpper[SoPlex::DISPLAYFREQ] = INT_MAX;
            _intParamDefault[SoPlex::DISPLAYFREQ] = 200;

            // verbosity level
            _intParamName[SoPlex::VERBOSITY] = "verbosity";
            _intParamDescription[SoPlex::VERBOSITY] = "verbosity level (0 - error, 1 - warning, 2 - debug, 3 - normal, 4 - high, 5 - full)";
            _intParamLower[SoPlex::VERBOSITY] = 0;
            _intParamUpper[SoPlex::VERBOSITY] = 5;
            _intParamDefault[SoPlex::VERBOSITY] = SoPlex::VERBOSITY_NORMAL;

            // type of simplifier
            _intParamName[SoPlex::SIMPLIFIER] = "simplifier";
            _intParamDescription[SoPlex::SIMPLIFIER] = "simplifier (0 - off, 1 - auto)";
            _intParamLower[SoPlex::SIMPLIFIER] = 0;
            _intParamUpper[SoPlex::SIMPLIFIER] = 1;
            _intParamDefault[SoPlex::SIMPLIFIER] = SoPlex::SIMPLIFIER_AUTO;

            // type of scaler
            _intParamName[SoPlex::SCALER] = "scaler";
            _intParamDescription[SoPlex::SCALER] = "scaling (0 - off, 1 - uni-equilibrium, 2 - bi-equilibrium, 3 - geometric, 4 - iterated geometric)";
            _intParamLower[SoPlex::SCALER] = 0;
            _intParamUpper[SoPlex::SCALER] = 4;
            _intParamDefault[SoPlex::SCALER] = SoPlex::SCALER_BIEQUI;

            // type of starter used to create crash basis
            _intParamName[SoPlex::STARTER] = "starter";
            _intParamDescription[SoPlex::STARTER] = "crash basis generated when starting from scratch (0 - none, 1 - weight, 2 - sum, 3 - vector)";
            _intParamLower[SoPlex::STARTER] = 0;
            _intParamUpper[SoPlex::STARTER] = 3;
            _intParamDefault[SoPlex::STARTER] = SoPlex::STARTER_OFF;

            // type of pricer
            _intParamName[SoPlex::PRICER] = "pricer";
            _intParamDescription[SoPlex::PRICER] = "pricing method (0 - auto, 1 - dantzig, 2 - parmult, 3 - devex, 4 - quicksteep, 5 - steep)";
            _intParamLower[SoPlex::PRICER] = 0;
            _intParamUpper[SoPlex::PRICER] = 5;
            _intParamDefault[SoPlex::PRICER] = SoPlex::PRICER_AUTO;

            // type of ratio test
            _intParamName[SoPlex::RATIOTESTER] = "ratiotester";
            _intParamDescription[SoPlex::RATIOTESTER] = "method for ratio test (0 - textbook, 1 - harris, 2 - fast, 3 - boundflipping)";
            _intParamLower[SoPlex::RATIOTESTER] = 0;
            _intParamUpper[SoPlex::RATIOTESTER] = 3;
            _intParamDefault[SoPlex::RATIOTESTER] = SoPlex::RATIOTESTER_BOUNDFLIPPING;

            // mode for synchronizing real and rational LP
            _intParamName[SoPlex::SYNCMODE] = "syncmode";
            _intParamDescription[SoPlex::SYNCMODE] = "mode for synchronizing real and rational LP (0 - store only real LP, 1 - auto, 2 - manual)";
            _intParamLower[SoPlex::SYNCMODE] = 0;
            _intParamUpper[SoPlex::SYNCMODE] = 2;
            _intParamDefault[SoPlex::SYNCMODE] = SoPlex::SYNCMODE_ONLYREAL;

            // mode for reading LP files
            _intParamName[SoPlex::READMODE] = "readmode";
            _intParamDescription[SoPlex::READMODE] = "mode for reading LP files (0 - floating-point, 1 - rational)";
            _intParamLower[SoPlex::READMODE] = 0;
            _intParamUpper[SoPlex::READMODE] = 1;
            _intParamDefault[SoPlex::READMODE] = SoPlex::READMODE_REAL;

            // mode for iterative refinement strategy
            _intParamName[SoPlex::SOLVEMODE] = "solvemode";
            _intParamDescription[SoPlex::SOLVEMODE] = "mode for iterative refinement strategy (0 - floating-point solve, 1 - auto, 2 - exact rational solve)";
            _intParamLower[SoPlex::SOLVEMODE] = 0;
            _intParamUpper[SoPlex::SOLVEMODE] = 2;
            _intParamDefault[SoPlex::SOLVEMODE] = SoPlex::SOLVEMODE_AUTO;

            // mode for iterative refinement strategy
            _intParamName[SoPlex::CHECKMODE] = "checkmode";
            _intParamDescription[SoPlex::CHECKMODE] = "mode for a posteriori feasibility checks (0 - floating-point check, 1 - auto, 2 - exact rational check)";
            _intParamLower[SoPlex::CHECKMODE] = 0;
            _intParamUpper[SoPlex::CHECKMODE] = 2;
            _intParamDefault[SoPlex::CHECKMODE] = SoPlex::CHECKMODE_AUTO;

            // type of timing
            _intParamName[SoPlex::TIMER] = "timer";
            _intParamDescription[SoPlex::TIMER] = "type of timer (1 - cputime, aka. usertime, 2 - wallclock time, 0 - no timing)";
            _intParamLower[SoPlex::TIMER] = 0;
            _intParamUpper[SoPlex::TIMER] = 2;
            _intParamDefault[SoPlex::TIMER] = SoPlex::TIMER_CPU;

            // mode for hyper sparse pricing
            _intParamName[SoPlex::HYPER_PRICING] = "hyperpricing";
            _intParamDescription[SoPlex::HYPER_PRICING] = "mode for hyper sparse pricing (0 - off, 1 - auto, 2 - always)";
            _intParamLower[SoPlex::HYPER_PRICING] = 0;
            _intParamUpper[SoPlex::HYPER_PRICING] = 2;
            _intParamDefault[SoPlex::HYPER_PRICING] = SoPlex::HYPER_PRICING_AUTO;

            // minimum number of stalling refinements since last pivot to trigger rational factorization
            _intParamName[SoPlex::RATFAC_MINSTALLS] = "ratfac_minstalls";
            _intParamDescription[SoPlex::RATFAC_MINSTALLS] = "minimum number of stalling refinements since last pivot to trigger rational factorization";
            _intParamLower[SoPlex::RATFAC_MINSTALLS] = 0;
            _intParamUpper[SoPlex::RATFAC_MINSTALLS] = INT_MAX;
            _intParamDefault[SoPlex::RATFAC_MINSTALLS] = 2;

            // primal feasibility tolerance
            _realParamName[SoPlex::FEASTOL] = "feastol";
            _realParamDescription[SoPlex::FEASTOL] = "primal feasibility tolerance";
            _realParamLower[SoPlex::FEASTOL] = 0.0;
            _realParamUpper[SoPlex::FEASTOL] = 1.0;
            _realParamDefault[SoPlex::FEASTOL] = 1e-6;

            // dual feasibility tolerance
            _realParamName[SoPlex::OPTTOL] = "opttol";
            _realParamDescription[SoPlex::OPTTOL] = "dual feasibility tolerance";
            _realParamLower[SoPlex::OPTTOL] = 0.0;
            _realParamUpper[SoPlex::OPTTOL] = 1.0;
            _realParamDefault[SoPlex::OPTTOL] = 1e-6;

            ///@todo define suitable values depending on Real type
            // general zero tolerance
            _realParamName[SoPlex::EPSILON_ZERO] = "epsilon_zero";
            _realParamDescription[SoPlex::EPSILON_ZERO] = "general zero tolerance";
            _realParamLower[SoPlex::EPSILON_ZERO] = 0.0;
            _realParamUpper[SoPlex::EPSILON_ZERO] = 1.0;
            _realParamDefault[SoPlex::EPSILON_ZERO] = DEFAULT_EPS_ZERO;

            ///@todo define suitable values depending on Real type
            // zero tolerance used in factorization
            _realParamName[SoPlex::EPSILON_FACTORIZATION] = "epsilon_factorization";
            _realParamDescription[SoPlex::EPSILON_FACTORIZATION] = "zero tolerance used in factorization";
            _realParamLower[SoPlex::EPSILON_FACTORIZATION] = 0.0;
            _realParamUpper[SoPlex::EPSILON_FACTORIZATION] = 1.0;
            _realParamDefault[SoPlex::EPSILON_FACTORIZATION] = DEFAULT_EPS_FACTOR;

            ///@todo define suitable values depending on Real type
            // zero tolerance used in update of the factorization
            _realParamName[SoPlex::EPSILON_UPDATE] = "epsilon_update";
            _realParamDescription[SoPlex::EPSILON_UPDATE] = "zero tolerance used in update of the factorization";
            _realParamLower[SoPlex::EPSILON_UPDATE] = 0.0;
            _realParamUpper[SoPlex::EPSILON_UPDATE] = 1.0;
            _realParamDefault[SoPlex::EPSILON_UPDATE] = DEFAULT_EPS_UPDATE;

            ///@todo define suitable values depending on Real type
            // pivot zero tolerance used in factorization
            _realParamName[SoPlex::EPSILON_PIVOT] = "epsilon_pivot";
            _realParamDescription[SoPlex::EPSILON_PIVOT] = "pivot zero tolerance used in factorization";
            _realParamLower[SoPlex::EPSILON_PIVOT] = 0.0;
            _realParamUpper[SoPlex::EPSILON_PIVOT] = 1.0;
            _realParamDefault[SoPlex::EPSILON_PIVOT] = DEFAULT_EPS_PIVOT;

            ///@todo define suitable values depending on Real type
            // infinity threshold
            _realParamName[SoPlex::INFTY] = "infty";
            _realParamDescription[SoPlex::INFTY] = "infinity threshold";
            _realParamLower[SoPlex::INFTY] = 1e10;
            _realParamUpper[SoPlex::INFTY] = 1e100;
            _realParamDefault[SoPlex::INFTY] = DEFAULT_INFINITY;

            // time limit in seconds (INFTY if unlimited)
            _realParamName[SoPlex::TIMELIMIT] = "timelimit";
            _realParamDescription[SoPlex::TIMELIMIT] = "time limit in seconds";
            _realParamLower[SoPlex::TIMELIMIT] = 0.0;
            _realParamUpper[SoPlex::TIMELIMIT] = DEFAULT_INFINITY;
            _realParamDefault[SoPlex::TIMELIMIT] = DEFAULT_INFINITY;

            // lower limit on objective value
            _realParamName[SoPlex::OBJLIMIT_LOWER] = "objlimit_lower";
            _realParamDescription[SoPlex::OBJLIMIT_LOWER] = "lower limit on objective value";
            _realParamLower[SoPlex::OBJLIMIT_LOWER] = -DEFAULT_INFINITY;
            _realParamUpper[SoPlex::OBJLIMIT_LOWER] = DEFAULT_INFINITY;
            _realParamDefault[SoPlex::OBJLIMIT_LOWER] = -DEFAULT_INFINITY;

            // upper limit on objective value
            _realParamName[SoPlex::OBJLIMIT_UPPER] = "objlimit_upper";
            _realParamDescription[SoPlex::OBJLIMIT_UPPER] = "upper limit on objective value";
            _realParamLower[SoPlex::OBJLIMIT_UPPER] = -DEFAULT_INFINITY;
            _realParamUpper[SoPlex::OBJLIMIT_UPPER] = DEFAULT_INFINITY;
            _realParamDefault[SoPlex::OBJLIMIT_UPPER] = DEFAULT_INFINITY;

            // working tolerance for feasibility in floating-point solver during iterative refinement
            _realParamName[SoPlex::FPFEASTOL] = "fpfeastol";
            _realParamDescription[SoPlex::FPFEASTOL] = "working tolerance for feasibility in floating-point solver during iterative refinement";
            _realParamLower[SoPlex::FPFEASTOL] = 1e-12;
            _realParamUpper[SoPlex::FPFEASTOL] = 1.0;
            _realParamDefault[SoPlex::FPFEASTOL] = 1e-9;

            // working tolerance for optimality in floating-point solver during iterative refinement
            _realParamName[SoPlex::FPOPTTOL] = "fpopttol";
            _realParamDescription[SoPlex::FPOPTTOL] = "working tolerance for optimality in floating-point solver during iterative refinement";
            _realParamLower[SoPlex::FPOPTTOL] = 1e-12;
            _realParamUpper[SoPlex::FPOPTTOL] = 1.0;
            _realParamDefault[SoPlex::FPOPTTOL] = 1e-9;

            // maximum increase of scaling factors between refinements
            _realParamName[SoPlex::MAXSCALEINCR] = "maxscaleincr";
            _realParamDescription[SoPlex::MAXSCALEINCR] = "maximum increase of scaling factors between refinements";
            _realParamLower[SoPlex::MAXSCALEINCR] = 1.0;
            _realParamUpper[SoPlex::MAXSCALEINCR] = DEFAULT_INFINITY;
            _realParamDefault[SoPlex::MAXSCALEINCR] = 1e25;

            // lower threshold in lifting (nonzero matrix coefficients with smaller absolute value will be reformulated)
            _realParamName[SoPlex::LIFTMINVAL] = "liftminval";
            _realParamDescription[SoPlex::LIFTMINVAL] = "lower threshold in lifting (nonzero matrix coefficients with smaller absolute value will be reformulated)";
            _realParamLower[SoPlex::LIFTMINVAL] = 0.0;
            _realParamUpper[SoPlex::LIFTMINVAL] = 0.1;
            _realParamDefault[SoPlex::LIFTMINVAL] = 0.000976562; // = 1/1024

            // upper threshold in lifting (nonzero matrix coefficients with larger absolute value will be reformulated)
            _realParamName[SoPlex::LIFTMAXVAL] = "liftmaxval";
            _realParamDescription[SoPlex::LIFTMAXVAL] = "lower threshold in lifting (nonzero matrix coefficients with smaller absolute value will be reformulated)";
            _realParamLower[SoPlex::LIFTMAXVAL] = 10.0;
            _realParamUpper[SoPlex::LIFTMAXVAL] = DEFAULT_INFINITY;
            _realParamDefault[SoPlex::LIFTMAXVAL] = 1024.0;

            // threshold for using sparse pricing (no. of violations need to be smaller than threshold * dimension of problem)
            _realParamName[SoPlex::SPARSITY_THRESHOLD] = "sparsity_threshold";
            _realParamDescription[SoPlex::SPARSITY_THRESHOLD] = "sparse pricing threshold (#violations < dimension * SPARSITY_THRESHOLD activates sparse pricing)";
            _realParamLower[SoPlex::SPARSITY_THRESHOLD] = 0.0;
            _realParamUpper[SoPlex::SPARSITY_THRESHOLD] = 1.0;
            _realParamDefault[SoPlex::SPARSITY_THRESHOLD] = 0.6;

            // threshold on number of rows vs. number of columns for switching from column to row representations in auto mode
            _realParamName[SoPlex::REPRESENTATION_SWITCH] = "representation_switch";
            _realParamDescription[SoPlex::REPRESENTATION_SWITCH] = "threshold on number of rows vs. number of columns for switching from column to row representations in auto mode";
            _realParamLower[SoPlex::REPRESENTATION_SWITCH] = 0.0;
            _realParamUpper[SoPlex::REPRESENTATION_SWITCH] = DEFAULT_INFINITY;
            _realParamDefault[SoPlex::REPRESENTATION_SWITCH] = DEFAULT_INFINITY;

            // geometric frequency at which to apply rational reconstruction
            _realParamName[SoPlex::RATREC_FREQ] = "ratrec_freq";
            _realParamDescription[SoPlex::RATREC_FREQ] = "geometric frequency at which to apply rational reconstruction";
            _realParamLower[SoPlex::RATREC_FREQ] = 1.0;
            _realParamUpper[SoPlex::RATREC_FREQ] = DEFAULT_INFINITY;
            _realParamDefault[SoPlex::RATREC_FREQ] = 1.2;

            // minimal reduction (sum of removed rows/cols) to continue simplification
            _realParamName[SoPlex::MINRED] = "minred";
            _realParamDescription[SoPlex::MINRED] = "minimal reduction (sum of removed rows/cols) to continue simplification";
            _realParamLower[SoPlex::MINRED] = 0.0;
            _realParamUpper[SoPlex::MINRED] = 1.0;
            _realParamDefault[SoPlex::MINRED] = 1e-4;

            _defaultsAndBoundsInitialized = true;
         }

         for( int i = 0; i < SoPlex::BOOLPARAM_COUNT; i++ )
            _boolParamValues[i] = _boolParamDefault[i];

         for( int i = 0; i < SoPlex::INTPARAM_COUNT; i++ )
            _intParamValues[i] = _intParamDefault[i];

         for( int i = 0; i < SoPlex::REALPARAM_COUNT; i++ )
            _realParamValues[i] = _realParamDefault[i];

#ifdef SOPLEX_WITH_RATIONALPARAM
         for( int i = 0; i < SoPlex::RATIONALPARAM_COUNT; i++ )
            _rationalParamValues[i] = _rationalParamDefault[i];
#endif
      }

      /// copy constructor
      Settings(const Settings& settings)
      {
         *this = settings;
      }

      /// assignment operator
      Settings& operator=(const Settings& settings)
      {
         for( int i = 0; i < SoPlex::BOOLPARAM_COUNT; i++ )
            _boolParamValues[i] = settings._boolParamValues[i];

         for( int i = 0; i < SoPlex::INTPARAM_COUNT; i++ )
            _intParamValues[i] = settings._intParamValues[i];

         for( int i = 0; i < SoPlex::REALPARAM_COUNT; i++ )
            _realParamValues[i] = settings._realParamValues[i];

#ifdef SOPLEX_WITH_RATIONALPARAM
         for( int i = 0; i < SoPlex::RATIONALPARAM_COUNT; i++ )
            _rationalParamValues[i] = settings._rationalParamValues[i];
#endif

         return *this;
      }
   };



   bool SoPlex::Settings::_defaultsAndBoundsInitialized = false;



   std::string SoPlex::Settings::_boolParamName[SoPlex::BOOLPARAM_COUNT];
   std::string SoPlex::Settings::_boolParamDescription[SoPlex::BOOLPARAM_COUNT];
   bool SoPlex::Settings::_boolParamDefault[SoPlex::BOOLPARAM_COUNT];



   std::string SoPlex::Settings::_intParamName[SoPlex::INTPARAM_COUNT];
   std::string SoPlex::Settings::_intParamDescription[SoPlex::INTPARAM_COUNT];
   int SoPlex::Settings::_intParamLower[SoPlex::INTPARAM_COUNT];
   int SoPlex::Settings::_intParamUpper[SoPlex::INTPARAM_COUNT];
   int SoPlex::Settings::_intParamDefault[SoPlex::INTPARAM_COUNT];



   std::string SoPlex::Settings::_realParamName[SoPlex::REALPARAM_COUNT];
   std::string SoPlex::Settings::_realParamDescription[SoPlex::REALPARAM_COUNT];
   Real SoPlex::Settings::_realParamLower[SoPlex::REALPARAM_COUNT];
   Real SoPlex::Settings::_realParamUpper[SoPlex::REALPARAM_COUNT];
   Real SoPlex::Settings::_realParamDefault[SoPlex::REALPARAM_COUNT];



#ifdef SOPLEX_WITH_RATIONALPARAM
   std::string SoPlex::Settings::_rationalParamName[SoPlex::RATIONALPARAM_COUNT];
   std::string SoPlex::Settings::_rationalParamDescription[SoPlex::RATIONALPARAM_COUNT];
   Rational SoPlex::Settings::_rationalParamLower[SoPlex::RATIONALPARAM_COUNT];
   Rational SoPlex::Settings::_rationalParamUpper[SoPlex::RATIONALPARAM_COUNT];
   Rational SoPlex::Settings::_rationalParamDefault[SoPlex::RATIONALPARAM_COUNT];
#endif



   /// default constructor
   SoPlex::SoPlex()
      : _statistics(0)
      , _currentSettings(0)
      , _scalerUniequi(false)
      , _scalerBiequi(true)
      , _scalerGeo1(1)
      , _scalerGeo8(8)
      , _simplifier(0)
      , _scaler(0)
      , _starter(0)
      , _rationalLP(0)
      , _unitMatrixRational(0)
      , _status(SPxSolver::UNKNOWN)
      , _hasBasis(false)
      , _hasSolReal(false)
      , _hasSolRational(false)
   {
      // transfer message handler
      _solver.setOutstream(spxout);
      _scalerUniequi.setOutstream(spxout);
      _scalerBiequi.setOutstream(spxout);
      _scalerGeo1.setOutstream(spxout);
      _scalerGeo8.setOutstream(spxout);

      // give lu factorization to solver
      _solver.setSolver(&_slufactor);

      // the real LP is initially stored in the solver; the rational LP is constructed, when the parameter SYNCMODE is
      // initialized in setSettings() below
      _realLP = &_solver;
      _isRealLPLoaded = true;
      _realLP->setOutstream(spxout);

      // initialize statistics
      spx_alloc(_statistics);
      _statistics = new (_statistics) Statistics();

      // initialize parameter settings to default
      spx_alloc(_currentSettings);
      _currentSettings = new (_currentSettings) Settings();
      setSettings(*_currentSettings, true, true);

      _lastSolveMode = intParam(SoPlex::SOLVEMODE);

      assert(_isConsistent());
   }



   /// assignment operator
   SoPlex& SoPlex::operator=(const SoPlex& rhs)
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

         // transfer the lu solver
         _solver.setSolver(&_slufactor);

         // initialize pointers for simplifier, scaler, and starter
         setIntParam(SoPlex::SIMPLIFIER, intParam(SoPlex::SIMPLIFIER), true, true);
         setIntParam(SoPlex::SCALER, intParam(SoPlex::SCALER), true, true);
         setIntParam(SoPlex::STARTER, intParam(SoPlex::STARTER), true, true);

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
            assert(intParam(SoPlex::SYNCMODE) == SYNCMODE_ONLYREAL);
            _rationalLP = 0;
         }
         else
         {
            assert(intParam(SoPlex::SYNCMODE) != SYNCMODE_ONLYREAL);
            _rationalLP = 0;
            spx_alloc(_rationalLP);
            _rationalLP = new (_rationalLP) SPxLPRational(*rhs._rationalLP);
         }

         // copy boolean flags
         _isRealLPLoaded = rhs._isRealLPLoaded;
         _hasSolReal = rhs._hasSolReal;
         _hasSolRational = rhs._hasSolRational;
         _hasBasis = rhs._hasBasis;
      }

      assert(_isConsistent());

      return *this;
   }



   /// copy constructor
   ///@todo improve performance by implementing a separate copy constructor
   SoPlex::SoPlex(const SoPlex& rhs)
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
   SoPlex::~SoPlex()
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



   /// returns number of rows
   int SoPlex::numRowsReal() const
   {
      assert(_realLP != 0);
      return _realLP->nRows();
   }



   /// returns number of columns
   int SoPlex::numColsReal() const
   {
      assert(_realLP != 0);
      return _realLP->nCols();
   }



   /// returns number of nonzeros
   int SoPlex::numNonzerosReal() const
   {
      assert(_realLP != 0);
      return _realLP->nNzos();
   }



   /// returns smallest non-zero element in absolute value
   Real SoPlex::minAbsNonzeroReal() const
   {
      assert(_realLP != 0);
      return _realLP->minAbsNzo();
   }



   /// returns biggest non-zero element in absolute value
   Real SoPlex::maxAbsNonzeroReal() const
   {
      assert(_realLP != 0);
      return _realLP->maxAbsNzo();
   }



   /// gets row \p i
   void SoPlex::getRowReal(int i, LPRowReal& lprow) const
   {
      assert(_realLP != 0);
      _realLP->getRow(i, lprow);
   }



   /// gets rows \p start, ..., \p end.
   void SoPlex::getRowsReal(int start, int end, LPRowSetReal& lprowset) const
   {
      assert(_realLP != 0);
      _realLP->getRows(start, end, lprowset);
   }



   /// returns vector of row \p i
   const SVectorReal& SoPlex::rowVectorReal(int i) const
   {
      assert(_realLP != 0);
      return _realLP->rowVector(i);
   }



   /// returns right-hand side vector
   const VectorReal& SoPlex::rhsReal() const
   {
      assert(_realLP != 0);
      return _realLP->rhs();
   }



   /// returns right-hand side of row \p i
   Real SoPlex::rhsReal(int i) const
   {
      assert(_realLP != 0);
      return _realLP->rhs(i);
   }



   /// returns left-hand side vector
   const VectorReal& SoPlex::lhsReal() const
   {
      assert(_realLP != 0);
      return _realLP->lhs();
   }



   /// returns left-hand side of row \p i
   Real SoPlex::lhsReal(int i) const
   {
      assert(_realLP != 0);
      return _realLP->lhs(i);
   }



   /// returns inequality type of row \p i
   LPRowReal::Type SoPlex::rowTypeReal(int i) const
   {
      assert(_realLP != 0);
      return _realLP->rowType(i);
   }



   /// gets column \p i
   void SoPlex::getColReal(int i, LPColReal& lpcol) const
   {
      assert(_realLP != 0);
      return _realLP->getCol(i, lpcol);
   }



   /// gets columns \p start, ..., \p end
   void SoPlex::getColsReal(int start, int end, LPColSetReal& lpcolset) const
   {
      assert(_realLP != 0);
      return _realLP->getCols(start, end, lpcolset);
   }



   /// returns vector of column \p i
   const SVectorReal& SoPlex::colVectorReal(int i) const
   {
      assert(_realLP != 0);
      return _realLP->colVector(i);
   }



   /// returns upper bound vector
   const VectorReal& SoPlex::upperReal() const
   {
      assert(_realLP != 0);
      return _realLP->upper();
   }



   /// returns upper bound of column \p i
   Real SoPlex::upperReal(int i) const
   {
      assert(_realLP != 0);
      return _realLP->upper(i);
   }



   /// returns lower bound vector
   const VectorReal& SoPlex::lowerReal() const
   {
      assert(_realLP != 0);
      return _realLP->lower();
   }



   /// returns lower bound of column \p i
   Real SoPlex::lowerReal(int i) const
   {
      assert(_realLP != 0);
      return _realLP->lower(i);
   }



   /// gets objective function vector
   void SoPlex::getObjReal(VectorReal& obj) const
   {
      assert(_realLP != 0);
      _realLP->getObj(obj);
   }



   /// returns objective value of column \p i
   Real SoPlex::objReal(int i) const
   {
      assert(_realLP != 0);
      return _realLP->obj(i);
   }



   /// returns objective function vector after transformation to a maximization problem; since this is how it is stored
   /// internally, this is generally faster
   const VectorReal& SoPlex::maxObjReal() const
   {
      assert(_realLP != 0);
      return _realLP->maxObj();
   }



   /// returns objective value of column \p i after transformation to a maximization problem; since this is how it is
   /// stored internally, this is generally faster
   Real SoPlex::maxObjReal(int i) const
   {
      assert(_realLP != 0);
      return _realLP->maxObj(i);
   }



   /// gets number of available dual norms
   void SoPlex::getNdualNorms(int& nnormsRow, int& nnormsCol) const
   {
      _solver.getNdualNorms(nnormsRow, nnormsCol);
   }



   /// gets steepest edge norms and returns false if they are not available
   bool SoPlex::getDualNorms(int& nnormsRow, int& nnormsCol, Real* norms) const
   {
      return _solver.getDualNorms(nnormsRow, nnormsCol, norms);
   }



   /// sets steepest edge norms and returns false if that's not possible
   bool SoPlex::setDualNorms(int nnormsRow, int nnormsCol, Real* norms)
   {
      return _solver.setDualNorms(nnormsRow, nnormsCol, norms);
   }



   /// returns number of rows
   int SoPlex::numRowsRational() const
   {
      assert(_rationalLP != 0);
      return _rationalLP->nRows();
   }



   /// returns number of columns
   int SoPlex::numColsRational() const
   {
      assert(_rationalLP != 0);
      return _rationalLP->nCols();
   }



   /// returns number of nonzeros
   int SoPlex::numNonzerosRational() const
   {
      assert(_rationalLP != 0);
      return _rationalLP->nNzos();
   }



   /// returns smallest non-zero element in absolute value
   Rational SoPlex::minAbsNonzeroRational() const
   {
      assert(_rationalLP != 0);
      return _rationalLP->minAbsNzo();
   }



   /// returns biggest non-zero element in absolute value
   Rational SoPlex::maxAbsNonzeroRational() const
   {
      assert(_rationalLP != 0);
      return _rationalLP->maxAbsNzo();
   }



   /// gets row \p i
   void SoPlex::getRowRational(int i, LPRowRational& lprow) const
   {
      assert(_rationalLP != 0);
      _rationalLP->getRow(i, lprow);
   }



   /// gets rows \p start, ..., \p end.
   void SoPlex::getRowsRational(int start, int end, LPRowSetRational& lprowset) const
   {
      assert(_rationalLP != 0);
      _rationalLP->getRows(start, end, lprowset);
   }



   /// returns vector of row \p i
   const SVectorRational& SoPlex::rowVectorRational(int i) const
   {
      assert(_rationalLP != 0);
      return _rationalLP->rowVector(i);
   }



   /// returns right-hand side vector
   const VectorRational& SoPlex::rhsRational() const
   {
      assert(_rationalLP != 0);
      return _rationalLP->rhs();
   }



   /// returns right-hand side of row \p i
   const Rational& SoPlex::rhsRational(int i) const
   {
      assert(_rationalLP != 0);
      return _rationalLP->rhs(i);
   }



   /// returns left-hand side vector
   const VectorRational& SoPlex::lhsRational() const
   {
      assert(_rationalLP != 0);
      return _rationalLP->lhs();
   }



   /// returns left-hand side of row \p i
   const Rational& SoPlex::lhsRational(int i) const
   {
      assert(_rationalLP != 0);
      return _rationalLP->lhs(i);
   }



   /// returns inequality type of row \p i
   LPRowRational::Type SoPlex::rowTypeRational(int i) const
   {
      assert(_rationalLP != 0);
      return _rationalLP->rowType(i);
   }



   /// gets column \p i
   void SoPlex::getColRational(int i, LPColRational& lpcol) const
   {
      assert(_rationalLP != 0);
      return _rationalLP->getCol(i, lpcol);
   }



   /// gets columns \p start, ..., \p end
   void SoPlex::getColsRational(int start, int end, LPColSetRational& lpcolset) const
   {
      assert(_rationalLP != 0);
      return _rationalLP->getCols(start, end, lpcolset);
   }



   /// returns vector of column \p i
   const SVectorRational& SoPlex::colVectorRational(int i) const
   {
      assert(_rationalLP != 0);
      return _rationalLP->colVector(i);
   }



   /// returns upper bound vector
   const VectorRational& SoPlex::upperRational() const
   {
      assert(_rationalLP != 0);
      return _rationalLP->upper();
   }



   /// returns upper bound of column \p i
   const Rational& SoPlex::upperRational(int i) const
   {
      assert(_rationalLP != 0);
      return _rationalLP->upper(i);
   }



   /// returns lower bound vector
   const VectorRational& SoPlex::lowerRational() const
   {
      assert(_rationalLP != 0);
      return _rationalLP->lower();
   }



   /// returns lower bound of column \p i
   const Rational& SoPlex::lowerRational(int i) const
   {
      assert(_rationalLP != 0);
      return _rationalLP->lower(i);
   }



   /// gets objective function vector
   void SoPlex::getObjRational(VectorRational& obj) const
   {
      assert(_rationalLP != 0);
      _rationalLP->getObj(obj);
   }



   /// gets objective value of column \p i
   void SoPlex::getObjRational(int i, Rational& obj) const
   {
      obj = maxObjRational(i);
      if( intParam(SoPlex::OBJSENSE) == SoPlex::OBJSENSE_MINIMIZE )
         obj *= -1;
   }



   /// returns objective value of column \p i
   Rational SoPlex::objRational(int i) const
   {
      assert(_rationalLP != 0);
      return _rationalLP->obj(i);
   }



   /// returns objective function vector after transformation to a maximization problem; since this is how it is stored
   /// internally, this is generally faster
   const VectorRational& SoPlex::maxObjRational() const
   {
      assert(_rationalLP != 0);
      return _rationalLP->maxObj();
   }



   /// returns objective value of column \p i after transformation to a maximization problem; since this is how it is
   /// stored internally, this is generally faster
   const Rational& SoPlex::maxObjRational(int i) const
   {
      assert(_rationalLP != 0);
      return _rationalLP->maxObj(i);
   }



   /// adds a single row
   void SoPlex::addRowReal(const LPRowReal& lprow)
   {
      assert(_realLP != 0);

      _addRowReal(lprow);

      if( intParam(SoPlex::SYNCMODE) == SYNCMODE_AUTO )
      {
         _rationalLP->addRow(lprow);
         _rowTypes.append(_rangeTypeReal(lprow.lhs(), lprow.rhs()));
      }

      _invalidateSolution();
   }



   /// adds multiple rows
   void SoPlex::addRowsReal(const LPRowSetReal& lprowset)
   {
      assert(_realLP != 0);

      _addRowsReal(lprowset);

      if( intParam(SoPlex::SYNCMODE) == SYNCMODE_AUTO )
      {
         _rationalLP->addRows(lprowset);
         for( int i = 0; i < lprowset.num(); i++ )
            _rowTypes.append(_rangeTypeReal(lprowset.lhs(i), lprowset.rhs(i)));
      }

      _invalidateSolution();
   }



   /// adds a single column
   void SoPlex::addColReal(const LPColReal& lpcol)
   {
      assert(_realLP != 0);

      _addColReal(lpcol);

      if( intParam(SoPlex::SYNCMODE) == SYNCMODE_AUTO )
      {
         _rationalLP->addCol(lpcol);
         _colTypes.append(_rangeTypeReal(lpcol.lower(), lpcol.upper()));
      }

      _invalidateSolution();
   }



   /// adds multiple columns
   void SoPlex::addColsReal(const LPColSetReal& lpcolset)
   {
      assert(_realLP != 0);

      _addColsReal(lpcolset);

      if( intParam(SoPlex::SYNCMODE) == SYNCMODE_AUTO )
      {
         _rationalLP->addCols(lpcolset);
         for( int i = 0; i < lpcolset.num(); i++ )
            _colTypes.append(_rangeTypeReal(lpcolset.lower(i), lpcolset.upper(i)));
      }

      _invalidateSolution();
   }



   /// replaces row \p i with \p lprow
   void SoPlex::changeRowReal(int i, const LPRowReal& lprow)
   {
      assert(_realLP != 0);

      _changeRowReal(i, lprow);

      if( intParam(SoPlex::SYNCMODE) == SYNCMODE_AUTO )
      {
         _rationalLP->changeRow(i, lprow);
         _rowTypes[i] = _rangeTypeReal(lprow.lhs(), lprow.rhs());
      }

      _invalidateSolution();
   }



   /// changes left-hand side vector for constraints to \p lhs
   void SoPlex::changeLhsReal(const VectorReal& lhs)
   {
      assert(_realLP != 0);

      _changeLhsReal(lhs);

      if( intParam(SoPlex::SYNCMODE) == SYNCMODE_AUTO )
      {
         _rationalLP->changeLhs(DVectorRational(lhs));
         for( int i = 0; i < numRowsRational(); i++ )
            _rowTypes[i] = _rangeTypeRational(_rationalLP->lhs(i), _rationalLP->rhs(i));
      }

      _invalidateSolution();
   }



   /// changes left-hand side of row \p i to \p lhs
   void SoPlex::changeLhsReal(int i, const Real& lhs)
   {
      assert(_realLP != 0);

      _changeLhsReal(i, lhs);

      if( intParam(SoPlex::SYNCMODE) == SYNCMODE_AUTO )
      {
         _rationalLP->changeLhs(i, lhs);
         _rowTypes[i] = _rangeTypeRational(_rationalLP->lhs(i), _rationalLP->rhs(i));
      }

      _invalidateSolution();
   }



   /// changes right-hand side vector to \p rhs
   void SoPlex::changeRhsReal(const VectorReal& rhs)
   {
      assert(_realLP != 0);

      _changeRhsReal(rhs);

      if( intParam(SoPlex::SYNCMODE) == SYNCMODE_AUTO )
      {
         _rationalLP->changeRhs(DVectorRational(rhs));
         for( int i = 0; i < numRowsRational(); i++ )
            _rowTypes[i] = _rangeTypeRational(_rationalLP->lhs(i), _rationalLP->rhs(i));
      }

      _invalidateSolution();
   }



   /// changes right-hand side of row \p i to \p rhs
   void SoPlex::changeRhsReal(int i, const Real& rhs)
   {
      assert(_realLP != 0);

      _changeRhsReal(i, rhs);

      if( intParam(SoPlex::SYNCMODE) == SYNCMODE_AUTO )
      {
         _rationalLP->changeRhs(i, rhs);
         _rowTypes[i] = _rangeTypeRational(_rationalLP->lhs(i), _rationalLP->rhs(i));
      }

      _invalidateSolution();
   }



   /// changes left- and right-hand side vectors
   void SoPlex::changeRangeReal(const VectorReal& lhs, const VectorReal& rhs)
   {
      assert(_realLP != 0);

      _changeRangeReal(lhs, rhs);

      if( intParam(SoPlex::SYNCMODE) == SYNCMODE_AUTO )
      {
         _rationalLP->changeRange(DVectorRational(lhs), DVectorRational(rhs));
         for( int i = 0; i < numRowsRational(); i++ )
            _rowTypes[i] = _rangeTypeReal(lhs[i], rhs[i]);
      }

      _invalidateSolution();
   }



   /// changes left- and right-hand side of row \p i
   void SoPlex::changeRangeReal(int i, const Real& lhs, const Real& rhs)
   {
      assert(_realLP != 0);

      _changeRangeReal(i,lhs, rhs);

      if( intParam(SoPlex::SYNCMODE) == SYNCMODE_AUTO )
      {
         _rationalLP->changeRange(i, lhs, rhs);
         _rowTypes[i] = _rangeTypeReal(lhs, rhs);
      }

      _invalidateSolution();
   }



   /// replaces column \p i with \p lpcol
   void SoPlex::changeColReal(int i, const LPColReal& lpcol)
   {
      assert(_realLP != 0);

      _changeColReal(i, lpcol);

      if( intParam(SoPlex::SYNCMODE) == SYNCMODE_AUTO )
      {
         _rationalLP->changeCol(i, lpcol);
         _colTypes[i] = _rangeTypeReal(lpcol.lower(), lpcol.upper());
      }

      _invalidateSolution();
   }



   /// changes vector of lower bounds to \p lower
   void SoPlex::changeLowerReal(const VectorReal& lower)
   {
      assert(_realLP != 0);

      _changeLowerReal(lower);

      if( intParam(SoPlex::SYNCMODE) == SYNCMODE_AUTO )
      {
         _rationalLP->changeLower(DVectorRational(lower));
         for( int i = 0; i < numColsRational(); i++ )
            _colTypes[i] = _rangeTypeRational(_rationalLP->lower(i), _rationalLP->upper(i));
      }


      _invalidateSolution();
   }



   /// changes lower bound of column i to \p lower
   void SoPlex::changeLowerReal(int i, const Real& lower)
   {
      assert(_realLP != 0);

      _changeLowerReal(i, lower);

      if( intParam(SoPlex::SYNCMODE) == SYNCMODE_AUTO )
      {
         _rationalLP->changeLower(i, lower);
         _colTypes[i] = _rangeTypeRational(_rationalLP->lower(i), _rationalLP->upper(i));
      }

      _invalidateSolution();
   }



   /// changes vector of upper bounds to \p upper
   void SoPlex::changeUpperReal(const VectorReal& upper)
   {
      assert(_realLP != 0);

      _changeUpperReal(upper);

      if( intParam(SoPlex::SYNCMODE) == SYNCMODE_AUTO )
      {
         _rationalLP->changeUpper(DVectorRational(upper));
         for( int i = 0; i < numColsRational(); i++ )
            _colTypes[i] = _rangeTypeRational(_rationalLP->lower(i), _rationalLP->upper(i));
      }

      _invalidateSolution();
   }



   /// changes \p i 'th upper bound to \p upper
   void SoPlex::changeUpperReal(int i, const Real& upper)
   {
      assert(_realLP != 0);

      _changeUpperReal(i, upper);

      if( intParam(SoPlex::SYNCMODE) == SYNCMODE_AUTO )
      {
         _rationalLP->changeUpper(i, upper);
         _colTypes[i] = _rangeTypeRational(_rationalLP->lower(i), _rationalLP->upper(i));
      }

      _invalidateSolution();
   }



   /// changes vectors of column bounds to \p lower and \p upper
   void SoPlex::changeBoundsReal(const VectorReal& lower, const VectorReal& upper)
   {
      assert(_realLP != 0);

      _changeBoundsReal(lower, upper);

      if( intParam(SoPlex::SYNCMODE) == SYNCMODE_AUTO )
      {
         _rationalLP->changeBounds(DVectorRational(lower), DVectorRational(upper));
         for( int i = 0; i < numColsRational(); i++ )
            _colTypes[i] = _rangeTypeReal(lower[i], upper[i]);
      }

      _invalidateSolution();
   }



   /// changes bounds of column \p i to \p lower and \p upper
   void SoPlex::changeBoundsReal(int i, const Real& lower, const Real& upper)
   {
      assert(_realLP != 0);

      _changeBoundsReal(i, lower, upper);

      if( intParam(SoPlex::SYNCMODE) == SYNCMODE_AUTO )
      {
         _rationalLP->changeBounds(i, lower, upper);
         _colTypes[i] = _rangeTypeReal(lower, upper);
      }
      _invalidateSolution();
   }



   /// changes objective function vector to \p obj
   void SoPlex::changeObjReal(const VectorReal& obj)
   {
      assert(_realLP != 0);

      _realLP->changeObj(obj);

      if( intParam(SoPlex::SYNCMODE) == SYNCMODE_AUTO )
         _rationalLP->changeObj(DVectorRational(obj));

      _invalidateSolution();
   }



   /// changes objective coefficient of column i to \p obj
   void SoPlex::changeObjReal(int i, const Real& obj)
   {
      assert(_realLP != 0);

      _realLP->changeObj(i, obj);

      if( intParam(SoPlex::SYNCMODE) == SYNCMODE_AUTO )
         _rationalLP->changeObj(i, obj);

      _invalidateSolution();
   }



   /// changes matrix entry in row \p i and column \p j to \p val
   void SoPlex::changeElementReal(int i, int j, const Real& val)
   {
      assert(_realLP != 0);

      _changeElementReal(i, j, val);

      if( intParam(SoPlex::SYNCMODE) == SYNCMODE_AUTO )
         _rationalLP->changeElement(i, j, val);

      _invalidateSolution();
   }



   /// removes row \p i
   void SoPlex::removeRowReal(int i)
   {
      assert(_realLP != 0);

      _removeRowReal(i);

      if( intParam(SoPlex::SYNCMODE) == SYNCMODE_AUTO )
      {
         _rationalLP->removeRow(i);
         _rowTypes[i] = _rowTypes[_rationalLP->nRows()];
         _rowTypes.reSize(_rationalLP->nRows());
         assert(_rowTypes[i] == _rangeTypeRational(lhsRational(i), rhsRational(i)));
      }

      _invalidateSolution();
   }



   /// removes all rows with an index \p i such that \p perm[i] < 0; upon completion, \p perm[i] >= 0 indicates the
   /// new index where row \p i has been moved to; note that \p perm must point to an array of size at least
   /// #numRowsReal()
   void SoPlex::removeRowsReal(int perm[])
   {
      assert(_realLP != 0);

      const int oldsize = numRowsReal();
      _removeRowsReal(perm);

      if( intParam(SoPlex::SYNCMODE) == SYNCMODE_AUTO )
      {
         _rationalLP->removeRows(perm);
         for( int i = 0; i < oldsize; i++ )
         {
            if( perm[i] >= 0 )
               _rowTypes[perm[i]] = _rowTypes[i];
         }
         _rowTypes.reSize(_rationalLP->nRows());
         for( int i = 0; i < numRowsRational(); i++ )
         {
            assert(_rowTypes[i] == _rangeTypeRational(lhsRational(i), rhsRational(i)));
         }
      }

      _invalidateSolution();
   }



   /// remove all rows with indices in array \p idx of size \p n; an array \p perm of size #numRowsReal() may be passed
   /// as buffer memory
   void SoPlex::removeRowsReal(int idx[], int n, int perm[])
   {
      if( perm == 0 )
      {
         DataArray< int > p(numRowsReal());
         _idxToPerm(idx, n, p.get_ptr(), numRowsReal());
         SoPlex::removeRowsReal(p.get_ptr());
      }
      else
      {
         _idxToPerm(idx, n, perm, numRowsReal());
         SoPlex::removeRowsReal(perm);
      }
   }



   /// removes rows \p start to \p end including both; an array \p perm of size #numRowsReal() may be passed as buffer
   /// memory
   void SoPlex::removeRowRangeReal(int start, int end, int perm[])
   {
      if( perm == 0 )
      {
         DataArray< int > p(numRowsReal());
         _rangeToPerm(start, end, p.get_ptr(), numRowsReal());
         SoPlex::removeRowsReal(p.get_ptr());
      }
      else
      {
         _rangeToPerm(start, end, perm, numRowsReal());
         SoPlex::removeRowsReal(perm);
      }
   }



   /// removes column i
   void SoPlex::removeColReal(int i)
   {
      assert(_realLP != 0);

      _removeColReal(i);

      if( intParam(SoPlex::SYNCMODE) == SYNCMODE_AUTO )
      {
         _rationalLP->removeCol(i);
         _colTypes[i] = _colTypes[_rationalLP->nCols()];
         _colTypes.reSize(_rationalLP->nCols());
         assert(_colTypes[i] == _rangeTypeRational(lowerRational(i), upperRational(i)));
      }

      _invalidateSolution();
   }



   /// removes all columns with an index \p i such that \p perm[i] < 0; upon completion, \p perm[i] >= 0 indicates the
   /// new index where column \p i has been moved to; note that \p perm must point to an array of size at least
   /// #numColsReal()
   void SoPlex::removeColsReal(int perm[])
   {
      assert(_realLP != 0);

      const int oldsize = numColsReal();
      _removeColsReal(perm);

      if( intParam(SoPlex::SYNCMODE) == SYNCMODE_AUTO )
      {
         _rationalLP->removeCols(perm);
         for( int i = 0; i < oldsize; i++ )
         {
            if( perm[i] >= 0 )
               _colTypes[perm[i]] = _colTypes[i];
         }
         _colTypes.reSize(_rationalLP->nCols());
         for( int i = 0; i < numColsRational(); i++ )
         {
            assert(_colTypes[i] == _rangeTypeRational(lowerRational(i), upperRational(i)));
         }
      }

      _invalidateSolution();
   }



   /// remove all columns with indices in array \p idx of size \p n; an array \p perm of size #numColsReal() may be
   /// passed as buffer memory
   void SoPlex::removeColsReal(int idx[], int n, int perm[])
   {
      if( perm == 0 )
      {
         DataArray< int > p(numColsReal());
         _idxToPerm(idx, n, p.get_ptr(), numColsReal());
         SoPlex::removeColsReal(p.get_ptr());
      }
      else
      {
         _idxToPerm(idx, n, perm, numColsReal());
         SoPlex::removeColsReal(perm);
      }
   }



   /// removes columns \p start to \p end including both; an array \p perm of size #numColsReal() may be passed as
   /// buffer memory
   void SoPlex::removeColRangeReal(int start, int end, int perm[])
   {
      if( perm == 0 )
      {
         DataArray< int > p(numColsReal());
         _rangeToPerm(start, end, p.get_ptr(), numColsReal());
         SoPlex::removeColsReal(p.get_ptr());
      }
      else
      {
         _rangeToPerm(start, end, perm, numColsReal());
         SoPlex::removeColsReal(perm);
      }
   }



   /// clears the LP
   void SoPlex::clearLPReal()
   {
      assert(_realLP != 0);

      _realLP->clear();
      _hasBasis = false;

      if( intParam(SoPlex::SYNCMODE) == SYNCMODE_AUTO )
      {
         _rationalLP->clear();
         _rowTypes.clear();
         _colTypes.clear();
      }

      _invalidateSolution();
   }



   /// synchronizes real LP with rational LP, i.e., copies (rounded) rational LP into real LP, if sync mode is manual
   void SoPlex::syncLPReal()
   {
      assert(_isConsistent());

      if( intParam(SoPlex::SYNCMODE) == SYNCMODE_MANUAL )
         _syncLPReal();
   }



   /// adds a single row
   void SoPlex::addRowRational(const LPRowRational& lprow)
   {
      assert(_rationalLP != 0);

      if( intParam(SoPlex::SYNCMODE) == SYNCMODE_ONLYREAL )
         return;

      _rationalLP->addRow(lprow);
      _rowTypes.append(_rangeTypeRational(lprow.lhs(), lprow.rhs()));

      if( intParam(SoPlex::SYNCMODE) == SYNCMODE_AUTO )
         _addRowReal(lprow);

      _invalidateSolution();
   }



#ifdef SOPLEX_WITH_GMP
   /// adds a single row
   void SoPlex::addRowRational(const mpq_t* lhs, const mpq_t* rowValues, const int* rowIndices, const int rowSize, const mpq_t* rhs)
   {
      assert(_rationalLP != 0);

      if( intParam(SoPlex::SYNCMODE) == SYNCMODE_ONLYREAL )
         return;

      _rationalLP->addRow(lhs, rowValues, rowIndices, rowSize, rhs);
      int i = numRowsRational() - 1;
      _rowTypes.append(_rangeTypeRational(_rationalLP->lhs(i), _rationalLP->rhs(i)));

      if( intParam(SoPlex::SYNCMODE) == SYNCMODE_AUTO )
         _addRowReal(Real(lhsRational(i)), DSVectorReal(_rationalLP->rowVector(i)), Real(rhsRational(i)));

      _invalidateSolution();
   }



   /// adds a set of rows
   void SoPlex::addRowsRational(const mpq_t* lhs, const mpq_t* rowValues, const int* rowIndices, const int* rowStarts, const int* rowLengths, const int numRows, const int numValues, const mpq_t* rhs)
   {
      assert(_rationalLP != 0);

      if( intParam(SoPlex::SYNCMODE) == SYNCMODE_ONLYREAL )
         return;

      _rationalLP->addRows(lhs, rowValues, rowIndices, rowStarts, rowLengths, numRows, numValues, rhs);
      for( int i = numRowsRational() - numRows; i < numRowsRational(); i++ )
         _rowTypes.append(_rangeTypeRational(_rationalLP->lhs(i), _rationalLP->rhs(i)));

      if( intParam(SoPlex::SYNCMODE) == SYNCMODE_AUTO )
      {
         LPRowSetReal lprowset;
         for( int i = numRowsRational() - numRows; i < numRowsRational(); i++ )
            lprowset.add(Real(lhsRational(i)), DSVectorReal(_rationalLP->rowVector(i)), Real(rhsRational(i)));
         _addRowsReal(lprowset);
      }

      _invalidateSolution();
   }
#endif



   /// adds multiple rows
   void SoPlex::addRowsRational(const LPRowSetRational& lprowset)
   {
      assert(_rationalLP != 0);

      if( intParam(SoPlex::SYNCMODE) == SYNCMODE_ONLYREAL )
         return;

      _rationalLP->addRows(lprowset);
      for( int i = 0; i < lprowset.num(); i++ )
         _rowTypes.append(_rangeTypeRational(lprowset.lhs(i), lprowset.rhs(i)));

      if( intParam(SoPlex::SYNCMODE) == SYNCMODE_AUTO )
         _addRowsReal(lprowset);

      _invalidateSolution();
   }



   /// adds a single column
   void SoPlex::addColRational(const LPColRational& lpcol)
   {
      assert(_rationalLP != 0);

      if( intParam(SoPlex::SYNCMODE) == SYNCMODE_ONLYREAL )
         return;

      _rationalLP->addCol(lpcol);
      _colTypes.append(_rangeTypeRational(lpcol.lower(), lpcol.upper()));

      if( intParam(SoPlex::SYNCMODE) == SYNCMODE_AUTO )
         _addColReal(lpcol);

      _invalidateSolution();
   }



#ifdef SOPLEX_WITH_GMP
   /// adds a single column
   void SoPlex::addColRational(const mpq_t* obj, const mpq_t* lower, const mpq_t* colValues, const int* colIndices, const int colSize, const mpq_t* upper)
   {
      assert(_rationalLP != 0);

      if( intParam(SoPlex::SYNCMODE) == SYNCMODE_ONLYREAL )
         return;

      _rationalLP->addCol(obj, lower, colValues, colIndices, colSize, upper);
      int i = numColsRational() - 1;
      _colTypes.append(_rangeTypeRational(lowerRational(i), upperRational(i)));

      if( intParam(SoPlex::SYNCMODE) == SYNCMODE_AUTO )
         _addColReal(Real(maxObjRational(i)) * (intParam(SoPlex::OBJSENSE) == SoPlex::OBJSENSE_MAXIMIZE ? 1.0 : -1.0),
            Real(lowerRational(i)), DSVectorReal(_rationalLP->colVector(i)), Real(upperRational(i)));

      _invalidateSolution();
   }



   /// adds a set of columns
   void SoPlex::addColsRational(const mpq_t* obj, const mpq_t* lower, const mpq_t* colValues, const int* colIndices, const int* colStarts, const int* colLengths, const int numCols, const int numValues, const mpq_t* upper)
   {
      assert(_rationalLP != 0);

      if( intParam(SoPlex::SYNCMODE) == SYNCMODE_ONLYREAL )
         return;

      _rationalLP->addCols(obj, lower, colValues, colIndices, colStarts, colLengths, numCols, numValues, upper);
      for( int i = numColsRational() - numCols; i < numColsRational(); i++ )
         _colTypes.append(_rangeTypeRational(lowerRational(i), upperRational(i)));

      if( intParam(SoPlex::SYNCMODE) == SYNCMODE_AUTO )
      {
         LPColSetReal lpcolset;
         for( int i = numColsRational() - numCols; i < numColsRational(); i++ )
            lpcolset.add(Real(maxObjRational(i)) * (intParam(SoPlex::OBJSENSE) == SoPlex::OBJSENSE_MAXIMIZE ? 1.0 : -1.0),
               Real(lowerRational(i)), DSVectorReal(_rationalLP->colVector(i)), Real(upperRational(i)));
         _addColsReal(lpcolset);
      }

      _invalidateSolution();
   }
#endif



   /// adds multiple columns
   void SoPlex::addColsRational(const LPColSetRational& lpcolset)
   {
      assert(_rationalLP != 0);

      if( intParam(SoPlex::SYNCMODE) == SYNCMODE_ONLYREAL )
         return;

      _rationalLP->addCols(lpcolset);
      for( int i = 0; i < lpcolset.num(); i++ )
         _colTypes.append(_rangeTypeRational(lpcolset.lower(i), lpcolset.upper(i)));

      if( intParam(SoPlex::SYNCMODE) == SYNCMODE_AUTO )
         _addColsReal(lpcolset);

      _invalidateSolution();
   }



   /// replaces row \p i with \p lprow
   void SoPlex::changeRowRational(int i, const LPRowRational& lprow)
   {
      assert(_rationalLP != 0);

      if( intParam(SoPlex::SYNCMODE) == SYNCMODE_ONLYREAL )
         return;

      _rationalLP->changeRow(i, lprow);
      _rowTypes[i] = _rangeTypeRational(lprow.lhs(), lprow.rhs());

      if( intParam(SoPlex::SYNCMODE) == SYNCMODE_AUTO )
         _changeRowReal(i, lprow);

      _invalidateSolution();
   }



   /// changes left-hand side vector for constraints to \p lhs
   void SoPlex::changeLhsRational(const VectorRational& lhs)
   {
      assert(_rationalLP != 0);

      if( intParam(SoPlex::SYNCMODE) == SYNCMODE_ONLYREAL )
         return;

      _rationalLP->changeLhs(lhs);
      for( int i = 0; i < numRowsRational(); i++ )
         _rowTypes[i] = _rangeTypeRational(lhs[i], _rationalLP->rhs(i));

      if( intParam(SoPlex::SYNCMODE) == SYNCMODE_AUTO )
         _changeLhsReal(DVectorReal(lhs));

      _invalidateSolution();
   }



   /// changes left-hand side of row \p i to \p lhs
   void SoPlex::changeLhsRational(int i, const Rational& lhs)
   {
      assert(_rationalLP != 0);

      if( intParam(SoPlex::SYNCMODE) == SYNCMODE_ONLYREAL )
         return;

      _rationalLP->changeLhs(i, lhs);
      _rowTypes[i] = _rangeTypeRational(lhs, _rationalLP->rhs(i));

      if( intParam(SoPlex::SYNCMODE) == SYNCMODE_AUTO )
         _changeLhsReal(i, Real(lhs));

      _invalidateSolution();
   }



#ifdef SOPLEX_WITH_GMP
   /// changes left-hand side of row \p i to \p lhs
   void SoPlex::changeLhsRational(int i, const mpq_t* lhs)
   {
      assert(_rationalLP != 0);

      if( intParam(SoPlex::SYNCMODE) == SYNCMODE_ONLYREAL )
         return;

      _rationalLP->changeLhs(i, lhs);
      _rowTypes[i] = _rangeTypeRational(_rationalLP->lhs(i), _rationalLP->rhs(i));

      if( intParam(SoPlex::SYNCMODE) == SYNCMODE_AUTO )
         _changeLhsReal(i, Real(lhsRational(i)));

      _invalidateSolution();
   }
#endif



   /// changes right-hand side vector to \p rhs
   void SoPlex::changeRhsRational(const VectorRational& rhs)
   {
      assert(_rationalLP != 0);

      if( intParam(SoPlex::SYNCMODE) == SYNCMODE_ONLYREAL )
         return;

      _rationalLP->changeRhs(rhs);
      for( int i = 0; i < numRowsRational(); i++ )
         _rowTypes[i] = _rangeTypeRational(_rationalLP->lhs(i), rhs[i]);

      if( intParam(SoPlex::SYNCMODE) == SYNCMODE_AUTO )
         _changeRhsReal(DVectorReal(rhs));

      _invalidateSolution();
   }



#ifdef SOPLEX_WITH_GMP
   /// changes right-hand side vector to \p rhs
   void SoPlex::changeRhsRational(const mpq_t* rhs, int rhsSize)
   {
      assert(_rationalLP != 0);

      if( intParam(SoPlex::SYNCMODE) == SYNCMODE_ONLYREAL )
         return;

      for( int i = 0; i < rhsSize; i++ )
      {
         _rationalLP->changeRhs(i, rhs[i]);
         _rowTypes[i] = _rangeTypeRational(_rationalLP->lhs(i), _rationalLP->rhs(i));
      }

      if( intParam(SoPlex::SYNCMODE) == SYNCMODE_AUTO )
         _changeRhsReal(DVectorReal(rhsRational()));

      _invalidateSolution();
   }
#endif



   /// changes right-hand side of row \p i to \p rhs
   void SoPlex::changeRhsRational(int i, const Rational& rhs)
   {
      assert(_rationalLP != 0);

      if( intParam(SoPlex::SYNCMODE) == SYNCMODE_ONLYREAL )
         return;

      _rationalLP->changeRhs(i, rhs);
      _rowTypes[i] = _rangeTypeRational(_rationalLP->lhs(i), rhs);

      if( intParam(SoPlex::SYNCMODE) == SYNCMODE_AUTO )
         _changeRhsReal(i, Real(rhs));

      _invalidateSolution();
   }



   /// changes left- and right-hand side vectors
   void SoPlex::changeRangeRational(const VectorRational& lhs, const VectorRational& rhs)
   {
      assert(_rationalLP != 0);

      if( intParam(SoPlex::SYNCMODE) == SYNCMODE_ONLYREAL )
         return;

      _rationalLP->changeRange(lhs, rhs);
      for( int i = 0; i < numRowsRational(); i++ )
         _rowTypes[i] = _rangeTypeRational(lhs[i], rhs[i]);

      if( intParam(SoPlex::SYNCMODE) == SYNCMODE_AUTO )
         _changeRangeReal(DVectorReal(lhs), DVectorReal(rhs));

      _invalidateSolution();
   }



   /// changes left- and right-hand side of row \p i
   void SoPlex::changeRangeRational(int i, const Rational& lhs, const Rational& rhs)
   {
      assert(_rationalLP != 0);

      if( intParam(SoPlex::SYNCMODE) == SYNCMODE_ONLYREAL )
         return;

      _rationalLP->changeRange(i, lhs, rhs);
      _rowTypes[i] = _rangeTypeRational(lhs, rhs);

      if( intParam(SoPlex::SYNCMODE) == SYNCMODE_AUTO )
         _changeRangeReal(i, Real(lhs), Real(rhs));

      _invalidateSolution();
   }



#ifdef SOPLEX_WITH_GMP
   /// changes left-hand side of row \p i to \p lhs
   void SoPlex::changeRangeRational(int i, const mpq_t* lhs, const mpq_t* rhs)
   {
      assert(_rationalLP != 0);

      if( intParam(SoPlex::SYNCMODE) == SYNCMODE_ONLYREAL )
         return;

      _rationalLP->changeRange(i, lhs, rhs);
      _rowTypes[i] = _rangeTypeRational(_rationalLP->lhs(i), _rationalLP->rhs(i));

      if( intParam(SoPlex::SYNCMODE) == SYNCMODE_AUTO )
         _changeRangeReal(i, Real(lhsRational(i)), Real(rhsRational(i)));

      _invalidateSolution();
   }
#endif



   /// replaces column \p i with \p lpcol
   void SoPlex::changeColRational(int i, const LPColRational& lpcol)
   {
      assert(_rationalLP != 0);

      if( intParam(SoPlex::SYNCMODE) == SYNCMODE_ONLYREAL )
         return;

      _rationalLP->changeCol(i, lpcol);
      _rowTypes[i] = _rangeTypeRational(lpcol.lower(), lpcol.upper());

      if( intParam(SoPlex::SYNCMODE) == SYNCMODE_AUTO )
         _changeColReal(i, lpcol);

      _invalidateSolution();
   }



   /// changes vector of lower bounds to \p lower
   void SoPlex::changeLowerRational(const VectorRational& lower)
   {
      assert(_rationalLP != 0);

      if( intParam(SoPlex::SYNCMODE) == SYNCMODE_ONLYREAL )
         return;

      _rationalLP->changeLower(lower);
      for( int i = 0; i < numColsRational(); i++ )
         _colTypes[i] = _rangeTypeRational(lower[i], _rationalLP->upper(i));

      if( intParam(SoPlex::SYNCMODE) == SYNCMODE_AUTO )
         _changeLowerReal(DVectorReal(lower));

      _invalidateSolution();
   }



   /// changes lower bound of column i to \p lower
   void SoPlex::changeLowerRational(int i, const Rational& lower)
   {
      assert(_rationalLP != 0);

      if( intParam(SoPlex::SYNCMODE) == SYNCMODE_ONLYREAL )
         return;

      _rationalLP->changeLower(i, lower);
      _colTypes[i] = _rangeTypeRational(lower, _rationalLP->upper(i));

      if( intParam(SoPlex::SYNCMODE) == SYNCMODE_AUTO )
         _changeLowerReal(i, Real(lower));

      _invalidateSolution();
   }



#ifdef SOPLEX_WITH_GMP
   /// changes lower bound of column i to \p lower
   void SoPlex::changeLowerRational(int i, const mpq_t* lower)
   {
      assert(_rationalLP != 0);

      if( intParam(SoPlex::SYNCMODE) == SYNCMODE_ONLYREAL )
         return;

      _rationalLP->changeLower(i, lower);
      _colTypes[i] = _rangeTypeRational(_rationalLP->lower(i), _rationalLP->upper(i));

      if( intParam(SoPlex::SYNCMODE) == SYNCMODE_AUTO )
         _changeLowerReal(i, Real(lowerRational(i)));

      _invalidateSolution();
   }
#endif



   /// changes vector of upper bounds to \p upper
   void SoPlex::changeUpperRational(const VectorRational& upper)
   {
      assert(_rationalLP != 0);

      if( intParam(SoPlex::SYNCMODE) == SYNCMODE_ONLYREAL )
         return;

      _rationalLP->changeUpper(upper);
      for( int i = 0; i < numColsRational(); i++ )
         _colTypes[i] = _rangeTypeRational(_rationalLP->lower(i), upper[i]);

      if( intParam(SoPlex::SYNCMODE) == SYNCMODE_AUTO )
         _changeUpperReal(DVectorReal(upper));

      _invalidateSolution();
   }



   /// changes \p i 'th upper bound to \p upper
   void SoPlex::changeUpperRational(int i, const Rational& upper)
   {
      assert(_rationalLP != 0);

      if( intParam(SoPlex::SYNCMODE) == SYNCMODE_ONLYREAL )
         return;

      _rationalLP->changeUpper(i, upper);
      _colTypes[i] = _rangeTypeRational(_rationalLP->lower(i), upper);

      if( intParam(SoPlex::SYNCMODE) == SYNCMODE_AUTO )
         _changeUpperReal(i, Real(upper));

      _invalidateSolution();
   }



#ifdef SOPLEX_WITH_GMP
   /// changes upper bound of column i to \p upper
   void SoPlex::changeUpperRational(int i, const mpq_t* upper)
   {
      assert(_rationalLP != 0);

      if( intParam(SoPlex::SYNCMODE) == SYNCMODE_ONLYREAL )
         return;

      _rationalLP->changeUpper(i, upper);
      _colTypes[i] = _rangeTypeRational(_rationalLP->lower(i), _rationalLP->upper(i));

      if( intParam(SoPlex::SYNCMODE) == SYNCMODE_AUTO )
         _changeUpperReal(i, Real(upperRational(i)));

      _invalidateSolution();
   }
#endif



   /// changes vectors of column bounds to \p lower and \p upper
   void SoPlex::changeBoundsRational(const VectorRational& lower, const VectorRational& upper)
   {
      assert(_rationalLP != 0);

      if( intParam(SoPlex::SYNCMODE) == SYNCMODE_ONLYREAL )
         return;

      _rationalLP->changeBounds(lower, upper);
      for( int i = 0; i < numColsRational(); i++ )
         _colTypes[i] = _rangeTypeRational(lower[i], upper[i]);

      if( intParam(SoPlex::SYNCMODE) == SYNCMODE_AUTO )
         _changeBoundsReal(DVectorReal(lower), DVectorReal(upper));

      _invalidateSolution();
   }



   /// changes bounds of column \p i to \p lower and \p upper
   void SoPlex::changeBoundsRational(int i, const Rational& lower, const Rational& upper)
   {
      assert(_rationalLP != 0);

      if( intParam(SoPlex::SYNCMODE) == SYNCMODE_ONLYREAL )
         return;

      _rationalLP->changeBounds(i, lower, upper);
      _colTypes[i] = _rangeTypeRational(lower, upper);

      if( intParam(SoPlex::SYNCMODE) == SYNCMODE_AUTO )
         _changeBoundsReal(i, Real(lower), Real(upper));

      _invalidateSolution();
   }



#ifdef SOPLEX_WITH_GMP
   /// changes bounds of column \p i to \p lower and \p upper
   void SoPlex::changeBoundsRational(int i, const mpq_t* lower, const mpq_t* upper)
   {
      assert(_rationalLP != 0);

      if( intParam(SoPlex::SYNCMODE) == SYNCMODE_ONLYREAL )
         return;

      _rationalLP->changeBounds(i, lower, upper);
      _colTypes[i] = _rangeTypeRational(_rationalLP->lower(i), _rationalLP->upper(i));

      if( intParam(SoPlex::SYNCMODE) == SYNCMODE_AUTO )
         _changeBoundsReal(i, Real(lowerRational(i)), Real(upperRational(i)));

      _invalidateSolution();
   }
#endif



   /// changes objective function vector to \p obj
   void SoPlex::changeObjRational(const VectorRational& obj)
   {
      assert(_rationalLP != 0);

      if( intParam(SoPlex::SYNCMODE) == SYNCMODE_ONLYREAL )
         return;

      _rationalLP->changeObj(obj);

      if( intParam(SoPlex::SYNCMODE) == SYNCMODE_AUTO )
         _realLP->changeObj(DVectorReal(obj));

      _invalidateSolution();
   }



   /// changes objective coefficient of column i to \p obj
   void SoPlex::changeObjRational(int i, const Rational& obj)
   {
      assert(_rationalLP != 0);

      if( intParam(SoPlex::SYNCMODE) == SYNCMODE_ONLYREAL )
         return;

      _rationalLP->changeObj(i, obj);

      if( intParam(SoPlex::SYNCMODE) == SYNCMODE_AUTO )
         _realLP->changeObj(i, Real(obj));

      _invalidateSolution();
   }



#ifdef SOPLEX_WITH_GMP
   /// changes objective coefficient of column i to \p obj
   void SoPlex::changeObjRational(int i, const mpq_t* obj)
   {
      assert(_rationalLP != 0);

      if( intParam(SoPlex::SYNCMODE) == SYNCMODE_ONLYREAL )
         return;

      _rationalLP->changeObj(i, obj);

      if( intParam(SoPlex::SYNCMODE) == SYNCMODE_AUTO )
         _realLP->changeObj(i, Real(objRational(i)));

      _invalidateSolution();
   }
#endif



   /// changes matrix entry in row \p i and column \p j to \p val
   void SoPlex::changeElementRational(int i, int j, const Rational& val)
   {
      assert(_rationalLP != 0);

      if( intParam(SoPlex::SYNCMODE) == SYNCMODE_ONLYREAL )
         return;

      _rationalLP->changeElement(i, j, val);

      if( intParam(SoPlex::SYNCMODE) == SYNCMODE_AUTO )
         _changeElementReal(i, j, Real(val));

      _invalidateSolution();
   }


#ifdef SOPLEX_WITH_GMP
   /// changes matrix entry in row \p i and column \p j to \p val
   void SoPlex::changeElementRational(int i, int j, const mpq_t* val)
   {
      assert(_rationalLP != 0);

      if( intParam(SoPlex::SYNCMODE) == SYNCMODE_ONLYREAL )
         return;

      _rationalLP->changeElement(i, j, val);

      if( intParam(SoPlex::SYNCMODE) == SYNCMODE_AUTO )
         _changeElementReal(i, j, mpq_get_d(*val));

      _invalidateSolution();
   }
#endif


   /// removes row \p i
   void SoPlex::removeRowRational(int i)
   {
      assert(_rationalLP != 0);

      if( intParam(SoPlex::SYNCMODE) == SYNCMODE_ONLYREAL )
         return;

      _rationalLP->removeRow(i);
      _rowTypes[i] = _rowTypes[_rationalLP->nRows()];
      _rowTypes.reSize(_rationalLP->nRows());
      assert(_rowTypes[i] == _rangeTypeRational(lhsRational(i), rhsRational(i)));

      if( intParam(SoPlex::SYNCMODE) == SYNCMODE_AUTO )
         _removeRowReal(i);

      _invalidateSolution();
   }



   /// removes all rows with an index \p i such that \p perm[i] < 0; upon completion, \p perm[i] >= 0 indicates the new
   /// index where row \p i has been moved to; note that \p perm must point to an array of size at least
   /// #numRowsRational()
   void SoPlex::removeRowsRational(int perm[])
   {
      assert(_rationalLP != 0);

      if( intParam(SoPlex::SYNCMODE) == SYNCMODE_ONLYREAL )
         return;

      const int oldsize = numRowsRational();
      _rationalLP->removeRows(perm);
      for( int i = 0; i < oldsize; i++ )
      {
         if( perm[i] >= 0 )
            _rowTypes[perm[i]] = _rowTypes[i];
      }
      _rowTypes.reSize(_rationalLP->nRows());
      for( int i = 0; i < numRowsRational(); i++ )
      {
         assert(_rowTypes[i] == _rangeTypeRational(lhsRational(i), rhsRational(i)));
      }


      if( intParam(SoPlex::SYNCMODE) == SYNCMODE_AUTO )
         _removeRowsReal(perm);

      _invalidateSolution();
   }



   /// remove all rows with indices in array \p idx of size \p n; an array \p perm of size #numRowsRational() may be
   /// passed as buffer memory
   void SoPlex::removeRowsRational(int idx[], int n, int perm[])
   {
      if( perm == 0 )
      {
         DataArray< int > p(numRowsRational());
         _idxToPerm(idx, n, p.get_ptr(), numRowsRational());
         SoPlex::removeRowsRational(p.get_ptr());
      }
      else
      {
         _idxToPerm(idx, n, perm, numRowsRational());
         SoPlex::removeRowsRational(perm);
      }
   }



   /// removes rows \p start to \p end including both; an array \p perm of size #numRowsRational() may be passed as
   /// buffer memory
   void SoPlex::removeRowRangeRational(int start, int end, int perm[])
   {
      if( perm == 0 )
      {
         DataArray< int > p(numRowsRational());
         _rangeToPerm(start, end, p.get_ptr(), numRowsRational());
         SoPlex::removeRowsRational(p.get_ptr());
      }
      else
      {
         _rangeToPerm(start, end, perm, numRowsRational());
         SoPlex::removeRowsRational(perm);
      }
   }



   /// removes column i
   void SoPlex::removeColRational(int i)
   {
      assert(_rationalLP != 0);

      if( intParam(SoPlex::SYNCMODE) == SYNCMODE_ONLYREAL )
         return;

      _rationalLP->removeCol(i);
      _colTypes[i] = _colTypes[_rationalLP->nCols()];
      _colTypes.reSize(_rationalLP->nCols());
      assert(_colTypes[i] == _rangeTypeRational(lowerRational(i), upperRational(i)));

      if( intParam(SoPlex::SYNCMODE) == SYNCMODE_AUTO )
         _removeColReal(i);

      _invalidateSolution();
   }



   /// removes all columns with an index \p i such that \p perm[i] < 0; upon completion, \p perm[i] >= 0 indicates the
   /// new index where column \p i has been moved to; note that \p perm must point to an array of size at least
   /// #numColsRational()
   void SoPlex::removeColsRational(int perm[])
   {
      assert(_rationalLP != 0);

      if( intParam(SoPlex::SYNCMODE) == SYNCMODE_ONLYREAL )
         return;

      const int oldsize = numColsRational();
      _rationalLP->removeCols(perm);
      for( int i = 0; i < oldsize; i++ )
      {
         if( perm[i] >= 0 )
            _colTypes[perm[i]] = _colTypes[i];
      }
      _colTypes.reSize(_rationalLP->nCols());
      for( int i = 0; i < numColsRational(); i++ )
      {
         assert(_colTypes[i] == _rangeTypeRational(lowerRational(i), upperRational(i)));
      }

      if( intParam(SoPlex::SYNCMODE) == SYNCMODE_AUTO )
         _removeColsReal(perm);

      _invalidateSolution();
   }



   /// remove all columns with indices in array \p idx of size \p n; an array \p perm of size #numColsRational() may be
   /// passed as buffer memory
   void SoPlex::removeColsRational(int idx[], int n, int perm[])
   {
      if( perm == 0 )
      {
         DataArray< int > p(numColsRational());
         _idxToPerm(idx, n, p.get_ptr(), numColsRational());
         SoPlex::removeColsRational(p.get_ptr());
      }
      else
      {
         _idxToPerm(idx, n, perm, numColsRational());
         SoPlex::removeColsRational(perm);
      }
   }



   /// removes columns \p start to \p end including both; an array \p perm of size #numColsRational() may be passed as
   /// buffer memory
   void SoPlex::removeColRangeRational(int start, int end, int perm[])
   {
      if( perm == 0 )
      {
         DataArray< int > p(numColsRational());
         _rangeToPerm(start, end, p.get_ptr(), numColsRational());
         SoPlex::removeColsRational(p.get_ptr());
      }
      else
      {
         _rangeToPerm(start, end, perm, numColsRational());
         SoPlex::removeColsRational(perm);
      }
   }



   /// clears the LP
   void SoPlex::clearLPRational()
   {
      assert(_rationalLP != 0);

      if( intParam(SoPlex::SYNCMODE) == SYNCMODE_ONLYREAL )
         return;

      _rationalLP->clear();
      _rowTypes.clear();
      _colTypes.clear();

      if( intParam(SoPlex::SYNCMODE) == SYNCMODE_AUTO )
      {
         _realLP->clear();
         _hasBasis = false;
      }

      _invalidateSolution();
   }



   /// synchronizes rational LP with real LP, i.e., copies real LP to rational LP, if sync mode is manual
   void SoPlex::syncLPRational()
   {
      assert(_isConsistent());

      if( intParam(SoPlex::SYNCMODE) == SYNCMODE_MANUAL )
         _syncLPRational();
   }



   /// solves the LP
   SPxSolver::Status SoPlex::solve()
   {
      assert(_isConsistent());

      // clear statistics
      _statistics->clearSolvingData();

      // the solution is no longer valid
      _invalidateSolution();

      // decide whether to solve the rational LP with iterative refinement or call the standard floating-point solver
      if( intParam(SoPlex::SOLVEMODE) == SOLVEMODE_REAL || (intParam(SoPlex::SOLVEMODE) == SOLVEMODE_AUTO
             && GE(realParam(SoPlex::FEASTOL), 1e-9) && GE(realParam(SoPlex::OPTTOL), 1e-9)) )
      {
         // ensure that tolerances are reasonable for the floating-point solver
         if( realParam(SoPlex::FEASTOL) < _currentSettings->_realParamLower[SoPlex::FPFEASTOL] )
         {
            MSG_WARNING( spxout, spxout << "Cannot call floating-point solver with feasibility tolerance below "
               << _currentSettings->_realParamLower[SoPlex::FPFEASTOL] << " - relaxing tolerance\n");
            _solver.setFeastol(_currentSettings->_realParamLower[SoPlex::FPFEASTOL]);
         }
         else
            _solver.setFeastol(realParam(SoPlex::FEASTOL));

         if( realParam(SoPlex::OPTTOL) < _currentSettings->_realParamLower[SoPlex::FPOPTTOL] )
         {
            MSG_WARNING( spxout, spxout << "Cannot call floating-point solver with optimality tolerance below "
               << _currentSettings->_realParamLower[SoPlex::FPOPTTOL] << " - relaxing tolerance\n");
            _solver.setOpttol(_currentSettings->_realParamLower[SoPlex::FPOPTTOL]);
         }
         else
            _solver.setOpttol(realParam(SoPlex::OPTTOL));

         _solveReal();
      }
      else if( intParam(SoPlex::SYNCMODE) == SYNCMODE_ONLYREAL )
      {
         _syncLPRational();
         _solveRational();
      }
      else if( intParam(SoPlex::SYNCMODE) == SYNCMODE_MANUAL )
      {
#ifdef ENABLE_ADDITIONAL_CHECKS
         assert(areLPsInSync(true, true, false));
#else
         assert(areLPsInSync(true, false, false));
#endif

         _solveRational();

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

         _solveRational();
      }

      MSG_INFO1( spxout, spxout << "\n";
         printShortStatistics(spxout.getStream(SPxOut::INFO1));
         spxout << "\n" );

      return status();
   }



   /// returns the current solver status
   SPxSolver::Status SoPlex::status() const
   {
      return _status;
   }



   /// is a primal feasible solution available?
   bool SoPlex::hasPrimal() const
   {
      return (_hasSolReal && _solReal.hasPrimal()) || (_hasSolRational && _solRational.hasPrimal());
   }



   /// is a primal unbounded ray available?
   bool SoPlex::hasPrimalRay() const
   {
      return (_hasSolReal && _solReal.hasPrimalRay()) || (_hasSolRational && _solRational.hasPrimalRay());
   }



   /// is a dual feasible solution available?
   bool SoPlex::hasDual() const
   {
      return (_hasSolReal && _solReal.hasDual()) || (_hasSolRational && _solRational.hasDual());
   }



   /// is Farkas proof of infeasibility available?
   bool SoPlex::hasDualFarkas() const
   {
      return (_hasSolReal && _solReal.hasDualFarkas()) || (_hasSolRational && _solRational.hasDualFarkas());
   }



   /// returns the objective value if a primal or dual solution is available
   Real SoPlex::objValueReal()
   {
      assert(OBJSENSE_MAXIMIZE == 1);
      assert(OBJSENSE_MINIMIZE == -1);

      if( status() == SPxSolver::UNBOUNDED )
         return realParam(SoPlex::INFTY) * intParam(SoPlex::OBJSENSE);
      else if( status() == SPxSolver::INFEASIBLE )
         return -realParam(SoPlex::INFTY) * intParam(SoPlex::OBJSENSE);
      else if( hasPrimal() )
      {
         _syncRealSolution();
         return _solReal._primalObjVal;
      }
      else if( hasDual() )
      {
         _syncRealSolution();
         return _solReal._dualObjVal;
      }
      else
         return 0.0;
   }



   /// gets the primal solution vector if available; returns true on success
   bool SoPlex::getPrimalReal(VectorReal& vector)
   {
      if( hasPrimal() && vector.dim() >= numColsReal() )
      {
         _syncRealSolution();
         _solReal.getPrimal(vector);
         return true;
      }
      else
         return false;
   }



   /// gets the vector of slack values if available; returns true on success
   bool SoPlex::getSlacksReal(VectorReal& vector)
   {
      if( hasPrimal() && vector.dim() >= numRowsReal() )
      {
         _syncRealSolution();
         _solReal.getSlacks(vector);
         return true;
      }
      else
         return false;
   }



   /// gets the primal ray if available; returns true on success
   bool SoPlex::getPrimalRayReal(VectorReal& vector)
   {
      if( hasPrimalRay() && vector.dim() >= numColsReal() )
      {
         _syncRealSolution();
         _solReal.getPrimalRay(vector);
         return true;
      }
      else
         return false;
   }



   /// gets the dual solution vector if available; returns true on success
   bool SoPlex::getDualReal(VectorReal& vector)
   {
      if( hasDual() && vector.dim() >= numRowsReal() )
      {
         _syncRealSolution();
         _solReal.getDual(vector);
         return true;
      }
      else
         return false;
   }



   /// gets the vector of reduced cost values if available; returns true on success
   bool SoPlex::getRedCostReal(VectorReal& vector)
   {
      if( hasDual() && vector.dim() >= numColsReal() )
      {
         _syncRealSolution();
         _solReal.getRedCost(vector);
         return true;
      }
      else
         return false;
   }



   /// gets the Farkas proof if available; returns true on success
   bool SoPlex::getDualFarkasReal(VectorReal& vector)
   {
      if( hasDualFarkas() && vector.dim() >= numRowsReal() )
      {
         _syncRealSolution();
         _solReal.getDualFarkas(vector);
         return true;
      }
      else
         return false;
   }



   /// gets violation of bounds; returns true on success
   bool SoPlex::getBoundViolationReal(Real& maxviol, Real& sumviol)
   {
      if( !hasPrimal() )
         return false;

      _syncRealSolution();
      VectorReal& primal = _solReal._primal;
      assert(primal.dim() == numColsReal());

      maxviol = 0.0;
      sumviol = 0.0;

      for( int i = numColsReal() - 1; i >= 0; i-- )
      {
         Real viol = lowerReal(i) - primal[i];
         if( viol > 0.0 )
         {
            sumviol += viol;
            if( viol > maxviol )
               maxviol = viol;
         }

         viol = primal[i] - upperReal(i);
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
   bool SoPlex::getRowViolationReal(Real& maxviol, Real& sumviol)
   {
      if( !hasPrimal() )
         return false;

      _syncRealSolution();
      VectorReal& primal = _solReal._primal;
      assert(primal.dim() == numColsReal());

      DVectorReal activity(numRowsReal());
      _realLP->computePrimalActivity(primal, activity);
      maxviol = 0.0;
      sumviol = 0.0;

      for( int i = numRowsReal() - 1; i >= 0; i-- )
      {
         Real viol = lhsReal(i) - activity[i];
         if( viol > 0.0 )
         {
            sumviol += viol;
            if( viol > maxviol )
               maxviol = viol;
         }

         viol = activity[i] - rhsReal(i);
         if( viol > 0.0 )
         {
            sumviol += viol;
            if( viol > maxviol )
               maxviol = viol;
         }
      }

      return true;
   }



   /// gets violation of reduced costs; returns true on success
   bool SoPlex::getRedCostViolationReal(Real& maxviol, Real& sumviol)
   {
      if( !hasDual() || !hasBasis() )
         return false;

      _syncRealSolution();
      VectorReal& redcost = _solReal._redCost;
      assert(redcost.dim() == numColsReal());

      maxviol = 0.0;
      sumviol = 0.0;

      for( int c = numColsReal() - 1; c >= 0; c-- )
      {
         SPxSolver::VarStatus colStatus = basisColStatus(c);

         if( intParam(SoPlex::OBJSENSE) == OBJSENSE_MINIMIZE )
         {
            if( colStatus != SPxSolver::ON_UPPER && colStatus != SPxSolver::FIXED && redcost[c] < 0.0 )
            {
               sumviol += -redcost[c];
               if( redcost[c] < -maxviol )
                  maxviol = -redcost[c];
            }
            if( colStatus != SPxSolver::ON_LOWER && colStatus != SPxSolver::FIXED && redcost[c] > 0.0 )
            {
               sumviol += redcost[c];
               if( redcost[c] > maxviol )
                  maxviol = redcost[c];
            }
         }
         else
         {
            if( colStatus != SPxSolver::ON_UPPER && colStatus != SPxSolver::FIXED && redcost[c] > 0.0 )
            {
               sumviol += redcost[c];
               if( redcost[c] > maxviol )
                  maxviol = redcost[c];
            }
            if( colStatus != SPxSolver::ON_LOWER && colStatus != SPxSolver::FIXED && redcost[c] < 0.0 )
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
   bool SoPlex::getDualViolationReal(Real& maxviol, Real& sumviol)
   {
      if( !hasDual() || !hasBasis() )
         return false;

      _syncRealSolution();
      VectorReal& dual = _solReal._dual;
      assert(dual.dim() == numRowsReal());

      maxviol = 0.0;
      sumviol = 0.0;

      for( int r = numRowsReal() - 1; r >= 0; r-- )
      {
         SPxSolver::VarStatus rowStatus = basisRowStatus(r);

         if( intParam(SoPlex::OBJSENSE) == OBJSENSE_MINIMIZE )
         {
            if( rowStatus != SPxSolver::ON_UPPER && rowStatus != SPxSolver::FIXED && dual[r] < 0.0 )
            {
               sumviol += -dual[r];
               if( dual[r] < -maxviol )
                  maxviol = -dual[r];
            }
            if( rowStatus != SPxSolver::ON_LOWER && rowStatus != SPxSolver::FIXED && dual[r] > 0.0 )
            {
               sumviol += dual[r];
               if( dual[r] > maxviol )
                  maxviol = dual[r];
            }
         }
         else
         {
            if( rowStatus != SPxSolver::ON_UPPER && rowStatus != SPxSolver::FIXED && dual[r] > 0.0 )
            {
               sumviol += dual[r];
               if( dual[r] > maxviol )
                  maxviol = dual[r];
            }
            if( rowStatus != SPxSolver::ON_LOWER && rowStatus != SPxSolver::FIXED && dual[r] < 0.0 )
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
   Rational SoPlex::objValueRational()
   {
      assert(OBJSENSE_MAXIMIZE == 1);
      assert(OBJSENSE_MINIMIZE == -1);

      if( status() == SPxSolver::UNBOUNDED )
      {
         if( intParam(SoPlex::OBJSENSE) == OBJSENSE_MAXIMIZE )
            return _rationalPosInfty;
         else
            return _rationalNegInfty;
      }
      else if( status() == SPxSolver::INFEASIBLE )
      {
         if( intParam(SoPlex::OBJSENSE) == OBJSENSE_MAXIMIZE )
            return _rationalNegInfty;
         else
            return _rationalPosInfty;
      }
      else if( hasPrimal() )
      {
         _syncRationalSolution();
         return _solRational._primalObjVal;
      }
      else if( hasDual() )
      {
         _syncRationalSolution();
         return _solRational._dualObjVal;
      }
      else
         return Rational::ZERO;
   }



   /// gets the primal solution vector if available; returns true on success
   bool SoPlex::getPrimalRational(VectorRational& vector)
   {
      if( hasPrimal() && vector.dim() >= numColsRational() )
      {
         _syncRationalSolution();
         _solRational.getPrimal(vector);
         return true;
      }
      else
         return false;
   }



   /// gets the vector of slack values if available; returns true on success
   bool SoPlex::getSlacksRational(VectorRational& vector)
   {
      if( hasPrimal() && vector.dim() >= numRowsRational() )
      {
         _syncRationalSolution();
         _solRational.getSlacks(vector);
         return true;
      }
      else
         return false;
   }



   /// gets the primal ray if LP is unbounded; returns true on success
   bool SoPlex::getPrimalRayRational(VectorRational& vector)
   {
      if( hasPrimalRay() && vector.dim() >= numColsRational() )
      {
         _syncRationalSolution();
         _solRational.getPrimalRay(vector);
         return true;
      }
      else
         return false;
   }



   /// gets the dual solution vector if available; returns true on success
   bool SoPlex::getDualRational(VectorRational& vector)
   {
      if( hasDual() && vector.dim() >= numRowsRational() )
      {
         _syncRationalSolution();
         _solRational.getDual(vector);
         return true;
      }
      else
         return false;
   }



   /// gets the vector of reduced cost values if available; returns true on success
   bool SoPlex::getRedCostRational(VectorRational& vector)
   {
      if( hasDual() && vector.dim() >= numColsRational() )
      {
         _syncRationalSolution();
         _solRational.getRedCost(vector);
         return true;
      }
      else
         return false;
   }



   /// gets the Farkas proof if LP is infeasible; returns true on success
   bool SoPlex::getDualFarkasRational(VectorRational& vector)
   {
      if( hasDualFarkas() && vector.dim() >= numRowsRational() )
      {
         _syncRationalSolution();
         _solRational.getDualFarkas(vector);
         return true;
      }
      else
         return false;
   }



   /// gets violation of bounds; returns true on success
   bool SoPlex::getBoundViolationRational(Rational& maxviol, Rational& sumviol)
   {
      if( !hasPrimal() )
         return false;

      // if we have to synchronize, we do not measure time, because this would affect the solving statistics
      if( intParam(SoPlex::SYNCMODE) == SYNCMODE_ONLYREAL )
         _syncLPRational(false);

      _syncRationalSolution();
      VectorRational& primal = _solRational._primal;
      assert(primal.dim() == numColsRational());

      maxviol = 0;
      sumviol = 0;

      for( int i = numColsRational() - 1; i >= 0; i-- )
      {
         Rational viol = lowerRational(i) - primal[i];
         if( viol > 0 )
         {
            sumviol += viol;
            if( viol > maxviol )
               maxviol = viol;
         }

         viol = primal[i] - upperRational(i);
         if( viol > 0 )
         {
            sumviol += viol;
            if( viol > maxviol )
               maxviol = viol;
         }
      }

      return true;
   }



   /// gets violation of constraints; returns true on success
   bool SoPlex::getRowViolationRational(Rational& maxviol, Rational& sumviol)
   {
      if( !hasPrimal() )
         return false;

      // if we have to synchronize, we do not measure time, because this would affect the solving statistics
      if( intParam(SoPlex::SYNCMODE) == SYNCMODE_ONLYREAL )
         _syncLPRational(false);

      _syncRationalSolution();
      VectorRational& primal = _solRational._primal;
      assert(primal.dim() == numColsRational());

      DVectorRational activity(numRowsRational());
      _rationalLP->computePrimalActivity(primal, activity);
      maxviol = 0;
      sumviol = 0;

      for( int i = numRowsRational() - 1; i >= 0; i-- )
      {
         Rational viol = lhsRational(i) - activity[i];
         if( viol > 0 )
         {
            sumviol += viol;
            if( viol > maxviol )
               maxviol = viol;
         }

         viol = activity[i] - rhsRational(i);
         if( viol > 0 )
         {
            sumviol += viol;
            if( viol > maxviol )
               maxviol = viol;
         }
      }

      return true;
   }



   /// gets violation of reduced costs; returns true on success
   bool SoPlex::getRedCostViolationRational(Rational& maxviol, Rational& sumviol)
   {
      if( !hasDual() || !hasPrimal() )
         return false;

      // if we have to synchronize, we do not measure time, because this would affect the solving statistics
      if( intParam(SoPlex::SYNCMODE) == SYNCMODE_ONLYREAL )
         _syncLPRational(false);

      _syncRationalSolution();
      VectorRational& redcost = _solRational._redCost;
      assert(redcost.dim() == numColsRational());

      maxviol = 0;
      sumviol = 0;

      for( int c = numColsReal() - 1; c >= 0; c-- )
      {
         assert(!_hasBasis || basisColStatus(c) != SPxSolver::UNDEFINED);

         if( _colTypes[c] == RANGETYPE_FIXED )
         {
            assert(lowerRational(c) == upperRational(c));
            continue;
         }

         assert(!_hasBasis || basisColStatus(c) != SPxSolver::ON_LOWER || _solRational._primal[c] == lowerRational(c));
         assert(!_hasBasis || basisColStatus(c) != SPxSolver::ON_UPPER || _solRational._primal[c] == upperRational(c));
         assert(!_hasBasis || basisColStatus(c) != SPxSolver::FIXED || _solRational._primal[c] == lowerRational(c));
         assert(!_hasBasis || basisColStatus(c) != SPxSolver::FIXED || _solRational._primal[c] == upperRational(c));

         if( intParam(SoPlex::OBJSENSE) == OBJSENSE_MINIMIZE )
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
   bool SoPlex::getDualViolationRational(Rational& maxviol, Rational& sumviol)
   {
      if( !hasDual() || !hasPrimal() )
         return false;

      // if we have to synchronize, we do not measure time, because this would affect the solving statistics
      if( intParam(SoPlex::SYNCMODE) == SYNCMODE_ONLYREAL )
         _syncLPRational(false);

      _syncRationalSolution();
      VectorRational& dual = _solRational._dual;
      assert(dual.dim() == numRowsRational());

      maxviol = 0;
      sumviol = 0;

      for( int r = numRowsReal() - 1; r >= 0; r-- )
      {
         assert(!_hasBasis || basisRowStatus(r) != SPxSolver::UNDEFINED);

         if( _rowTypes[r] == RANGETYPE_FIXED )
         {
            assert(lhsRational(r) == rhsRational(r));
            continue;
         }

         assert(!_hasBasis || basisRowStatus(r) != SPxSolver::ON_LOWER || _solRational._slacks[r] <= lhsRational(r) + _rationalFeastol);
         assert(!_hasBasis || basisRowStatus(r) != SPxSolver::ON_UPPER || _solRational._slacks[r] >= rhsRational(r) - _rationalFeastol);
         assert(!_hasBasis || basisRowStatus(r) != SPxSolver::FIXED || _solRational._slacks[r] <= lhsRational(r) + _rationalFeastol);
         assert(!_hasBasis || basisRowStatus(r) != SPxSolver::FIXED || _solRational._slacks[r] >= rhsRational(r) - _rationalFeastol);

         if( intParam(SoPlex::OBJSENSE) == OBJSENSE_MINIMIZE )
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
   bool SoPlex::getPrimalRational(mpq_t* vector, const int size)
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
   bool SoPlex::getSlacksRational(mpq_t* vector, const int size)
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
   bool SoPlex::getPrimalRayRational(mpq_t* vector, const int size)
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
   bool SoPlex::getDualRational(mpq_t* vector, const int size)
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
   bool SoPlex::getRedCostRational(mpq_t* vector, const int size)
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
   bool SoPlex::getDualFarkasRational(mpq_t* vector, const int size)
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



   /// get size of primal solution
   int SoPlex::totalSizePrimalRational(const int base)
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
   int SoPlex::totalSizeDualRational(const int base)
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
   int SoPlex::dlcmSizePrimalRational(const int base)
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
   int SoPlex::dlcmSizeDualRational(const int base)
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
   int SoPlex::dmaxSizePrimalRational(const int base)
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
   int SoPlex::dmaxSizeDualRational(const int base)
   {
      if( hasDual() || hasDualFarkas() )
      {
         _syncRationalSolution();
         return _solRational.dmaxSizeDual(base);
      }
      else
         return 0;
   }



   /// is an advanced starting basis available?
   bool SoPlex::hasBasis() const
   {
      return _hasBasis;
   }



   /// returns the current basis status
   SPxBasis::SPxStatus SoPlex::basisStatus() const
   {
      if( !hasBasis() )
         return SPxBasis::NO_PROBLEM;
      else if( status() == SPxSolver::OPTIMAL )
         return SPxBasis::OPTIMAL;
      else if( status() == SPxSolver::UNBOUNDED )
         return SPxBasis::UNBOUNDED;
      else if( status() == SPxSolver::INFEASIBLE )
         return SPxBasis::INFEASIBLE;
      else if( hasPrimal() )
         return SPxBasis::PRIMAL;
      else if( hasDual() )
         return SPxBasis::DUAL;
      else
         return SPxBasis::REGULAR;
   }



   /// returns basis status for a single row
   SPxSolver::VarStatus SoPlex::basisRowStatus(int row) const
   {
      assert(row >= 0);
      assert(row < numRowsReal());

      // if no basis is available, return slack basis; if index is out of range, return basic status as for a newly
      // added row
      if( !hasBasis() || row < 0 || row >= numRowsReal() )
         return SPxSolver::BASIC;
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
   SPxSolver::VarStatus SoPlex::basisColStatus(int col) const
   {
      assert(col >= 0);
      assert(col < numColsReal());

      // if index is out of range, return nonbasic status as for a newly added unbounded column
      if( col < 0 || col >= numColsReal() )
      {
         return SPxSolver::ZERO;
      }
      // if no basis is available, return slack basis
      else if( !hasBasis() )
      {
         if( lowerReal(col) > -realParam(SoPlex::INFTY) )
            return SPxSolver::ON_LOWER;
         else if( upperReal(col) < realParam(SoPlex::INFTY) )
            return SPxSolver::ON_UPPER;
         else
            return SPxSolver::ZERO;
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
   void SoPlex::getBasis(SPxSolver::VarStatus rows[], SPxSolver::VarStatus cols[]) const
   {
      // if no basis is available, return slack basis
      if( !hasBasis() )
      {
         for( int i = numRowsReal() - 1; i >= 0; i-- )
            rows[i] = SPxSolver::BASIC;

         for( int i = numColsReal() - 1; i >= 0; i-- )
         {
            if( lowerReal(i) > -realParam(SoPlex::INFTY) )
               cols[i] = SPxSolver::ON_LOWER;
            else if( upperReal(i) < realParam(SoPlex::INFTY) )
               cols[i] = SPxSolver::ON_UPPER;
            else
               cols[i] = SPxSolver::ZERO;
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
         assert(numRowsReal() == _basisStatusRows.size());
         assert(numColsReal() == _basisStatusCols.size());

         for( int i = numRowsReal() - 1; i >= 0; i-- )
            rows[i] = _basisStatusRows[i];

         for( int i = numColsReal() - 1; i >= 0; i-- )
            cols[i] = _basisStatusCols[i];
      }
   }



   /// returns the indices of the basic columns and rows; basic column n gives value n, basic row m gives value -1-m
   void SoPlex::getBasisInd(int* bind) const
   {
      // if no basis is available, return slack basis
      if( !hasBasis() )
      {
         for( int i = 0; i < numRowsReal(); ++i )
            bind[i] = -1 - i;
      }
      // if the real LP is not loaded, the basis is stored in the basis arrays of this class
      else if( !_isRealLPLoaded )
      {
         int k = 0;

         assert(numRowsReal() == _basisStatusRows.size());
         assert(numColsReal() == _basisStatusCols.size());

         for( int i = 0; i < numRowsReal(); ++i )
         {
            if( _basisStatusRows[i] == SPxSolver::BASIC )
            {
               bind[k] = -1 - i;
               k++;
            }
         }

         for( int j = 0; j < numColsReal(); ++j )
         {
            if( _basisStatusCols[j] == SPxSolver::BASIC )
            {
               bind[k] = j;
               k++;
            }
         }

         assert(k == numRowsReal());
      }
      // if the real LP is loaded, the basis is stored in the solver and we need to distinguish between column and row
      // representation; ask the solver itself which representation it has, since the REPRESENTATION parameter of this
      // class might be set to automatic
      else if( _solver.rep() == SPxSolver::COLUMN )
      {
         for( int i = 0; i < numRowsReal(); ++i )
         {
            SPxId id = _solver.basis().baseId(i);
            bind[i] = (id.isSPxColId() ? _solver.number(id) : - 1 - _solver.number(id));
         }
      }
      // for row representation, return the complement of the basis; for this, we need to loop through all rows and columns
      else
      {
         assert(_solver.rep() == SPxSolver::ROW);

         int k = 0;

         for( int i = 0; i < numRowsReal(); ++i )
         {
            if( !_solver.isRowBasic(i) )
            {
               bind[k] = -1 - i;
               k++;
            }
         }

         for( int j = 0; j < numColsReal(); ++j )
         {
            if( !_solver.isColBasic(j) )
            {
               bind[k] = j;
               k++;
            }
         }

         assert(k == numRowsReal());
      }
   }



   /// computes an estimated condition number for the current basis matrix using the power method; returns true on success
   bool SoPlex::getEstimatedCondition(Real& condition)
   {
      _ensureRealLPLoaded();
      if( !_isRealLPLoaded )
         return false;

      if( _solver.basis().status() == SPxBasis::NO_PROBLEM )
         return false;

      condition = _solver.basis().getEstimatedCondition();

      return true;
   }

   /// computes the exact condition number for the current basis matrix using the power method; returns true on success
   bool SoPlex::getExactCondition(Real& condition)
   {
      _ensureRealLPLoaded();
      if( !_isRealLPLoaded )
         return false;

      if( _solver.basis().status() == SPxBasis::NO_PROBLEM )
         return false;

      condition = _solver.basis().getExactCondition();

      return true;
   }

   /// computes row r of basis inverse; returns true on success
   bool SoPlex::getBasisInverseRowReal(int r, Real* coef, int* inds, int* ninds)
   {
      assert(r >= 0);
      assert(r < numRowsReal());
      assert(coef != 0);

      if( !hasBasis() || r < 0 || r >= numRowsReal() )
         return false;

      _ensureRealLPLoaded();

      if( !_isRealLPLoaded )
         return false;

      // we need to distinguish between column and row representation; ask the solver itself which representation it
      // has, since the REPRESENTATION parameter of this class might be set to automatic
      if( _solver.rep() == SPxSolver::COLUMN )
      {
         int idx;
         SSVectorReal x(numRowsReal());
         try
         {
            _solver.basis().coSolve(x, _solver.unitVector(r));
         }
         catch( const SPxException& E )
         {
            MSG_ERROR( std::cerr << "Caught exception <" << E.what() << "> while computing basis inverse row.\n" );
            return false;
         }
         // copy sparse data to dense result vector based on coef array
         if( ninds != NULL && inds != NULL )
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
            VectorReal y(numRowsReal(), coef);
            y = x;
            if( ninds != NULL )
               *ninds = -1;
         }
      }
      else
      {
         assert(_solver.rep() == SPxSolver::ROW);

         // @todo should rhs be a reference?
         DSVector rhs(numColsReal());
         SSVector y(numColsReal());
         int* bind = 0;
         int index;

         // get ordering of column basis matrix
         spx_alloc(bind, numRowsReal());
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
            assert(index < numRowsReal());
            assert(!_solver.isRowBasic(index));

            // get row vector
            rhs = _solver.rowVector(index);
            rhs *= -1.0;
         }
         // r corresponds to a column vector
         else
         {
            // should be a valid column index and in the column basis matrix, i.e., not basic w.r.t. row representation
            assert(index < numColsReal());
            assert(!_solver.isColBasic(index));

            // get unit vector
            rhs = UnitVectorReal(index);
         }

         // solve system "y B = rhs", where B is the row basis matrix
         try
         {
            _solver.basis().solve(y, rhs);
         }
         catch( const SPxException& E )
         {
            MSG_ERROR( std::cerr << "Caught exception <" << E.what() << "> while computing basis inverse row.\n" );
            return false;
         }

         // initialize result vector x as zero
         memset(coef, 0, (unsigned int)numRowsReal() * sizeof(Real));

         // add nonzero entries
         for( int i = 0; i < numColsReal(); ++i )
         {
            SPxId id = _solver.basis().baseId(i);

            if( id.isSPxRowId() )
            {
               assert(_solver.number(id) >= 0);
               assert(_solver.number(id) < numRowsReal());
               assert(bind[r] >= 0 || _solver.number(id) != index);

               coef[_solver.number(id)] = y[i];
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
   bool SoPlex::getBasisInverseColReal(int c, Real* coef, int* inds, int* ninds)
   {
      assert(c >= 0);
      assert(c < numRowsReal());
      assert(coef != 0);

      if( !hasBasis() || c < 0 || c >= numRowsReal() )
         return false;

      _ensureRealLPLoaded();

      if( !_isRealLPLoaded )
         return false;

      // we need to distinguish between column and row representation; ask the solver itself which representation it
      // has, since the REPRESENTATION parameter of this class might be set to automatic
      if( _solver.rep() == SPxSolver::COLUMN )
      {
         int idx;
         SSVectorReal x(numColsReal());
         try
         {
            _solver.basis().solve(x, _solver.unitVector(c));
         }
         catch( const SPxException& E )
         {
            MSG_ERROR( std::cerr << "Caught exception <" << E.what() << "> while computing basis inverse row.\n" );
            return false;
         }
         // copy sparse data to dense result vector based on coef array
         if( ninds != NULL && inds != NULL )
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
            VectorReal y(numColsReal(), coef);
            y = x;
            if( ninds != NULL )
               *ninds = -1;
         }
      }
      else
      {
         assert(_solver.rep() == SPxSolver::ROW);

         // @todo should rhs be a reference?
         DSVectorReal rhs(numColsReal());
         SSVectorReal y(numColsReal());
         int* bind = 0;
         int index;

         // get ordering of column basis matrix
         spx_alloc(bind, numRowsReal());
         getBasisInd(bind);

         // get vector corresponding to requested index r
         index = bind[c];

         // c corresponds to a row vector
         if( index < 0 )
         {
            // transform index to actual row index
            index = -index - 1;

            // should be a valid row index and in the column basis matrix, i.e., not basic w.r.t. row representation
            assert(index >= 0);
            assert(index < numRowsReal());
            assert(!_solver.isRowBasic(index));

            // get row vector
            rhs = _solver.rowVector(index);
            rhs *= -1.0;
         }
         // c corresponds to a column vector
         else
         {
            // should be a valid column index and in the column basis matrix, i.e., not basic w.r.t. row representation
            assert(index < numColsReal());
            assert(!_solver.isColBasic(index));

            // get unit vector
            rhs = UnitVectorReal(index);
         }

         // solve system "y B = rhs", where B is the row basis matrix
         try
         {
            _solver.basis().coSolve(y, rhs);
         }
         catch( const SPxException& E )
         {
            MSG_ERROR( std::cerr << "Caught exception <" << E.what() << "> while computing basis inverse row.\n" );
            return false;
         }

         // initialize result vector x as zero
         memset(coef, 0, (unsigned int)numRowsReal() * sizeof(Real));

         // add nonzero entries
         for( int i = 0; i < numColsReal(); ++i )
         {
            SPxId id = _solver.basis().baseId(i);

            if( id.isSPxRowId() )
            {
               assert(_solver.number(id) >= 0);
               assert(_solver.number(id) < numRowsReal());
               assert(bind[c] >= 0 || _solver.number(id) != index);

               coef[_solver.number(id)] = y[i];
            }
         }

         // if r corresponds to a row vector, we have to add a 1 at position r
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
   bool SoPlex::getBasisInverseTimesVecReal(Real* rhs, Real* sol)
   {
      VectorReal v(numRowsReal(), rhs);
      VectorReal x(numRowsReal(), sol);

      if( !hasBasis() )
         return false;

      _ensureRealLPLoaded();

      if( !_isRealLPLoaded )
         return false;

      // we need to distinguish between column and row representation; ask the solver itself which representation it
      // has, since the REPRESENTATION parameter of this class might be set to automatic; in the column case we can use
      // the existing factorization
      if( _solver.rep() == SPxSolver::COLUMN )
      {
         // solve system "x = B^-1 * A_c" to get c'th column of B^-1 * A
         try
         {
            _solver.basis().solve(x, v);
         }
         catch( const SPxException& E )
         {
            MSG_ERROR( std::cerr << "Caught exception <" << E.what() << "> while solving with basis matrix.\n" );
            return false;
         }
      }
      else
      {
         assert(_solver.rep() == SPxSolver::ROW);

         DSVectorReal rowrhs(numColsReal());
         SSVectorReal y(numColsReal());
         int* bind = 0;

         // get ordering of column basis matrix
         spx_alloc(bind, numRowsReal());
         getBasisInd(bind);

         // fill right-hand side for row-based system
         for( int i = 0; i < numColsReal(); ++i )
         {
            SPxId id = _solver.basis().baseId(i);

            if( id.isSPxRowId() )
            {
               assert(_solver.number(id) >= 0);
               assert(_solver.number(id) < numRowsReal());

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
            MSG_ERROR( std::cerr << "Caught exception <" << E.what() << "> while solving with basis matrix.\n" );
            return false;
         }

         // fill result w.r.t. order given by bind
         for( int i = 0; i < numRowsReal(); ++i )
         {
            int index;

            index = bind[i];

            if( index < 0 )
            {
               index = -index-1;

               // should be a valid row index and in the column basis matrix, i.e., not basic w.r.t. row representation
               assert(index >= 0);
               assert(index < numRowsReal());
               assert(!_solver.isRowBasic(index));

               x[i] = v[index] - (rowVectorReal(index) * Vector(numColsReal(), y.get_ptr()));
            }
            else
            {
               // should be a valid column index and in the column basis matrix, i.e., not basic w.r.t. row representation
               assert(index >= 0);
               assert(index < numColsReal());
               assert(!_solver.isColBasic(index));

               x[i] = y[index];
            }
         }

         // free memory
         spx_free(bind);
      }
      return true;
   }



   /// sets starting basis via arrays of statuses
   void SoPlex::setBasis(SPxSolver::VarStatus rows[], SPxSolver::VarStatus cols[])
   {
      if( _isRealLPLoaded )
      {
         assert(numRowsReal() == _solver.nRows());
         assert(numColsReal() == _solver.nCols());

         _solver.setBasis(rows, cols);
         _hasBasis = (_solver.basis().status() > SPxBasis::NO_PROBLEM);
      }
      else
      {
         _basisStatusRows.reSize(numRowsReal());
         _basisStatusCols.reSize(numColsReal());

         for( int i = numRowsReal() - 1; i >= 0; i-- )
            _basisStatusRows[i] = rows[i];

         for( int j = numColsReal() - 1; j >= 0; j-- )
            _basisStatusCols[j] = cols[j];

         _hasBasis = true;
      }
   }



   /// clears starting basis
   void SoPlex::clearBasis()
   {
      _solver.reLoad();
      _status = _solver.status();
      _hasBasis = false;
   }



   /// number of iterations since last call to solve
   int SoPlex::numIterations() const
   {
      return _statistics->iterations;
   }



   /// time spent in last call to solve
   Real SoPlex::solveTime() const
   {
       return _statistics->solvingTime->time();
   }



   /// statistical information in form of a string
   std::string SoPlex::statisticString() const
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
   const char* SoPlex::getStarterName()
   {
      if( _starter )
         return _starter->getName();
      else
         return "none";
   }



   /// name of simplifier
   const char* SoPlex::getSimplifierName()
   {
      if( _simplifier )
         return _simplifier->getName();
      else
         return "none";
   }



   /// name of scaling method after simplifier
   const char* SoPlex::getScalerName()
   {
      if( _scaler )
         return _scaler->getName();
      else
         return "none";
   }



   /// name of currently loaded pricer
   const char* SoPlex::getPricerName()
   {
      return _solver.pricer()->getName();
   }



   /// name of currently loaded ratiotester
   const char* SoPlex::getRatiotesterName()
   {
      return _solver.ratiotester()->getName();
   }



   /// reads LP file in LP or MPS format according to READMODE parameter; gets row names, column names, and
   /// integer variables if desired; returns true on success
   bool SoPlex::readFile(const char* filename, NameSet* rowNames, NameSet* colNames, DIdxSet* intVars)
   {
      if( intParam(SoPlex::READMODE) == READMODE_REAL )
         return _readFileReal(filename, rowNames, colNames, intVars);
      else
         return _readFileRational(filename, rowNames, colNames, intVars);
   }

   /// writes real LP to file; LP or MPS format is chosen from the extension in \p filename; if \p rowNames and \p
   /// colNames are \c NULL, default names are used; if \p intVars is not \c NULL, the variables contained in it are
   /// marked as integer; returns true on success
   bool SoPlex::writeFileReal(const char* filename, const NameSet* rowNames, const NameSet* colNames, const DIdxSet* intVars) const
   {
      ///@todo implement return value
      _realLP->writeFile(filename, rowNames, colNames, intVars);
      return true;
   }



   /// writes rational LP to file; LP or MPS format is chosen from the extension in \p filename; if \p rowNames and \p
   /// colNames are \c NULL, default names are used; if \p intVars is not \c NULL, the variables contained in it are
   /// marked as integer; returns true on success
   bool SoPlex::writeFileRational(const char* filename, const NameSet* rowNames, const NameSet* colNames, const DIdxSet* intVars) const
   {
      if( intParam(SoPlex::SYNCMODE) == SYNCMODE_ONLYREAL )
         return false;
      else
      {
         assert(_rationalLP != 0);
         _rationalLP->writeFile(filename, rowNames, colNames, intVars);

         ///@todo implement return value
         return true;
      }
   }



   /// reads basis information from \p filename and returns true on success; if \p rowNames and \p colNames are \c NULL,
   /// default names are assumed; returns true on success
   bool SoPlex::readBasisFile(const char* filename, const NameSet* rowNames, const NameSet* colNames)
   {
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
      assert(_hasBasis == (_solver.basis().status() > SPxBasis::NO_PROBLEM));

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
      int numRows = numRowsReal();
      int numCols = numColsReal();

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
         _basisStatusRows[i] = SPxSolver::BASIC;

      for( int i = 0; i < numCols; i++ )
      {
         if( lowerReal(i) == upperReal(i) )
            _basisStatusCols[i] = SPxSolver::FIXED;
         else if( lowerReal(i) <= double(-realParam(SoPlex::INFTY)) && upperReal(i) >= double(realParam(SoPlex::INFTY)) )
            _basisStatusCols[i] = SPxSolver::ZERO;
         else if( lowerReal(i) <= double(-realParam(SoPlex::INFTY)) )
            _basisStatusCols[i] = SPxSolver::ON_UPPER;
         else
            _basisStatusCols[i] = SPxSolver::ON_LOWER;
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
               _basisStatusCols[c] = SPxSolver::BASIC;
               if( _rowTypes[r] == SoPlex::RANGETYPE_LOWER )
                  _basisStatusRows[r] = SPxSolver::ON_LOWER;
               else if( _rowTypes[r] == SoPlex::RANGETYPE_FIXED )
                  _basisStatusRows[r] = SPxSolver::FIXED;
               else
                  _basisStatusRows[r] = SPxSolver::ON_UPPER;
            }
            else if( !strcmp(mps.field1(), "XL") )
            {
               _basisStatusCols[c] = SPxSolver::BASIC;
               if( _rowTypes[r] == SoPlex::RANGETYPE_UPPER )
                  _basisStatusRows[r] = SPxSolver::ON_UPPER;
               else if( _rowTypes[r] == SoPlex::RANGETYPE_FIXED )
                  _basisStatusRows[r] = SPxSolver::FIXED;
               else
                  _basisStatusRows[r] = SPxSolver::ON_LOWER;
            }
            else if( !strcmp(mps.field1(), "UL") )
            {
               _basisStatusCols[c] = SPxSolver::ON_UPPER;
            }
            else if( !strcmp(mps.field1(), "LL") )
            {
               _basisStatusCols[c] = SPxSolver::ON_LOWER;
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
   bool SoPlex::writeBasisFile(const char* filename, const NameSet* rowNames, const NameSet* colNames, const bool cpxFormat) const
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
            assert(_basisStatusCols[col] != SPxSolver::UNDEFINED);

            if( _basisStatusCols[col] == SPxSolver::BASIC )
            {
               // find nonbasic row
               for( ; row < numRows; row++ )
               {
                  assert(_basisStatusRows[row] != SPxSolver::UNDEFINED);
                  if( _basisStatusRows[row] != SPxSolver::BASIC )
                     break;
               }

               assert(row != numRows);

               if( _basisStatusRows[row] == SPxSolver::ON_UPPER && (!cpxFormat || _rowTypes[row] == SoPlex::RANGETYPE_BOXED) )
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
               if( _basisStatusCols[col] == SPxSolver::ON_UPPER )
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
            assert(_basisStatusRows[row] == SPxSolver::BASIC);
         }
#endif

         return true;
      }
   }



   /// writes internal LP, basis information, and parameter settings; if \p rowNames and \p colNames are \c NULL,
   /// default names are used
   void SoPlex::writeStateReal(const char* filename, const NameSet* rowNames, const NameSet* colNames, const bool cpxFormat) const
   {
      std::string ofname;

      // write parameter settings
      ofname = std::string(filename) + ".set";
      saveSettingsFile(ofname.c_str());

      // write problem in MPS format
      ofname = std::string(filename) + ".mps";
      writeFileReal(ofname.c_str(), rowNames, colNames, 0);

      // write basis
      ofname = std::string(filename) + ".bas";
      writeBasisFile(ofname.c_str(), rowNames, colNames, cpxFormat);
   }



   /// writes internal LP, basis information, and parameter settings; if \p rowNames and \p colNames are \c NULL,
   /// default names are used
   void SoPlex::writeStateRational(const char* filename, const NameSet* rowNames, const NameSet* colNames, const bool cpxFormat) const
   {
      std::string ofname;

      // write parameter settings
      ofname = std::string(filename) + ".set";
      saveSettingsFile(ofname.c_str());

      // write problem in MPS format
      ofname = std::string(filename) + ".mps";
      writeFileRational(ofname.c_str(), rowNames, colNames, 0);

      // write basis
      ofname = std::string(filename) + ".bas";
      writeBasisFile(ofname.c_str(), rowNames, colNames, cpxFormat);
   }



   /// returns boolean parameter value
   bool SoPlex::boolParam(const BoolParam param) const
   {
      assert(param >= 0);
      assert(param < SoPlex::BOOLPARAM_COUNT);
      return _currentSettings->_boolParamValues[param];
   }



   /// returns integer parameter value
   int SoPlex::intParam(const IntParam param) const
   {
      assert(param >= 0);
      assert(param < INTPARAM_COUNT);
      return _currentSettings->_intParamValues[param];
   }



   /// returns real parameter value
   Real SoPlex::realParam(const RealParam param) const
   {
      assert(param >= 0);
      assert(param < REALPARAM_COUNT);
      return _currentSettings->_realParamValues[param];
   }



#ifdef SOPLEX_WITH_RATIONALPARAM
   /// returns rational parameter value
   Rational SoPlex::rationalParam(const RationalParam param) const
   {
      assert(param >= 0);
      assert(param < RATIONALPARAM_COUNT);
      return _currentSettings->_rationalParamValues[param];
   }
#endif



   /// returns current parameter settings
   const SoPlex::Settings& SoPlex::settings() const
   {
      return *_currentSettings;
   }



   /// sets boolean parameter value; returns true on success
   bool SoPlex::setBoolParam(const BoolParam param, const bool value, const bool quiet, const bool init)
   {
      assert(param >= 0);
      assert(param < SoPlex::BOOLPARAM_COUNT);
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
      case ACCEPTCYCLING:
         break;
      case RATREC:
         break;
      case POWERSCALING:
         break;
      case RATFACJUMP:
         break;
      case FEASRELAX:
         break;
      case ROWBOUNDFLIPS:
         _ratiotesterBoundFlipping.useBoundFlipsRow(value);
         break;
      default:
         return false;
      }

      _currentSettings->_boolParamValues[param] = value;
      return true;
   }



   /// sets integer parameter value; returns true on success
   bool SoPlex::setIntParam(const IntParam param, const int value, const bool quiet, const bool init)
   {
      assert(param >= 0);
      assert(param < INTPARAM_COUNT);
      assert(init || _isConsistent());

      if( !init && value == intParam(param) )
         return true;

      // check for a valid parameter value wrt bounds
      if( value < _currentSettings->_intParamLower[param] || value > _currentSettings->_intParamUpper[param] )
         return false;

      switch( param )
      {
      // objective sense
      case SoPlex::OBJSENSE:
         if( value != SoPlex::OBJSENSE_MAXIMIZE && value != SoPlex::OBJSENSE_MINIMIZE )
            return false;
         _realLP->changeSense(value == SoPlex::OBJSENSE_MAXIMIZE ? SPxLPReal::MAXIMIZE : SPxLPReal::MINIMIZE);
         if( _rationalLP != 0 )
            _rationalLP->changeSense(value == SoPlex::OBJSENSE_MAXIMIZE ? SPxLPRational::MAXIMIZE : SPxLPRational::MINIMIZE);
         _invalidateSolution();
         break;

      // type of computational form, i.e., column or row representation
      case SoPlex::REPRESENTATION:
         if( value != SoPlex::REPRESENTATION_COLUMN && value != SoPlex::REPRESENTATION_ROW && value != SoPlex::REPRESENTATION_AUTO )
            return false;
         break;

      // type of algorithm, i.e., primal or dual
      case SoPlex::ALGORITHM:
         // decide upon entering/leaving at solve time depending on representation
         break;

      // type of LU update
      case SoPlex::FACTOR_UPDATE_TYPE:
         if( value != SoPlex::FACTOR_UPDATE_TYPE_ETA && value != SoPlex::FACTOR_UPDATE_TYPE_FT )
            return false;
         _slufactor.setUtype(value == SoPlex::FACTOR_UPDATE_TYPE_ETA ? SLUFactor::ETA : SLUFactor::FOREST_TOMLIN);
         break;

      // maximum number of updates before fresh factorization
      case SoPlex::FACTOR_UPDATE_MAX:
         _solver.basis().setMaxUpdates(value);
         break;

      // iteration limit (-1 if unlimited)
      case SoPlex::ITERLIMIT:
         break;

      // refinement limit (-1 if unlimited)
      case SoPlex::REFLIMIT:
         break;

      // stalling refinement limit (-1 if unlimited)
      case SoPlex::STALLREFLIMIT:
         break;

      // display frequency
      case SoPlex::DISPLAYFREQ:
         _solver.setDisplayFreq(value);
         break;

      // verbosity level
      case SoPlex::VERBOSITY:
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
      case SoPlex::SIMPLIFIER:
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
      case SoPlex::SCALER:
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
         default:
            return false;
         }
         break;

      // type of starter used to create crash basis
      case SoPlex::STARTER:
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
      case SoPlex::PRICER:
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
      case SoPlex::SYNCMODE:
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
      case SoPlex::READMODE:
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
      case SoPlex::SOLVEMODE:
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
      case SoPlex::CHECKMODE:
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
      case SoPlex::RATIOTESTER:
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
      case SoPlex::TIMER:
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
      case SoPlex::HYPER_PRICING:
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
      case SoPlex::RATFAC_MINSTALLS:
         break;

      default:
         return false;
      }

      _currentSettings->_intParamValues[param] = value;
      return true;
   }



   /// sets real parameter value; returns true on success
   bool SoPlex::setRealParam(const RealParam param, const Real value, const bool quiet, const bool init)
   {
      assert(param >= 0);
      assert(param < REALPARAM_COUNT);
      assert(init || _isConsistent());

      if( !init && value == realParam(param) )
         return true;

      if( value < _currentSettings->_realParamLower[param] || value > _currentSettings->_realParamUpper[param] )
         return false;

      switch( param )
      {
      // primal feasibility tolerance; passed to the floating point solver only when calling solve()
      case SoPlex::FEASTOL:
         _rationalFeastol = value;
         break;

      // dual feasibility tolerance; passed to the floating point solver only when calling solve()
      case SoPlex::OPTTOL:
         _rationalOpttol = value;
         break;

      // general zero tolerance
      case SoPlex::EPSILON_ZERO:
         Param::setEpsilon(value);
         break;

      // zero tolerance used in factorization
      case SoPlex::EPSILON_FACTORIZATION:
         Param::setEpsilonFactorization(value);
         break;

      // zero tolerance used in update of the factorization
      case SoPlex::EPSILON_UPDATE:
         Param::setEpsilonUpdate(value);
         break;

      // pivot zero tolerance used in factorization (declare numerical singularity for small LU pivots)
      case SoPlex::EPSILON_PIVOT:
         Param::setEpsilonPivot(value);
         break;

      // infinity threshold
      case SoPlex::INFTY:
         _rationalPosInfty = value;
         _rationalNegInfty = -value;
         if( intParam(SoPlex::SYNCMODE) != SYNCMODE_ONLYREAL )
            _recomputeRangeTypesRational();
         break;

      // time limit in seconds (INFTY if unlimited)
      case SoPlex::TIMELIMIT:
         break;

      // lower limit on objective value is set in solveReal()
      case SoPlex::OBJLIMIT_LOWER:
         break;

      // upper limit on objective value is set in solveReal()
      case SoPlex::OBJLIMIT_UPPER:
         break;

      // working tolerance for feasibility in floating-point solver
      case SoPlex::FPFEASTOL:
         break;

      // working tolerance for optimality in floating-point solver
      case SoPlex::FPOPTTOL:
         break;

      // maximum increase of scaling factors between refinements
      case SoPlex::MAXSCALEINCR:
         _rationalMaxscaleincr = value;
         break;

      // lower threshold in lifting (nonzero matrix coefficients with smaller absolute value will be reformulated)
      case SoPlex::LIFTMINVAL:
         break;

      // upper threshold in lifting (nonzero matrix coefficients with larger absolute value will be reformulated)
      case SoPlex::LIFTMAXVAL:
         break;

      // threshold for sparse pricing
      case SoPlex::SPARSITY_THRESHOLD:
         break;

      // threshold on number of rows vs. number of columns for switching from column to row representations in auto mode
      case SoPlex::REPRESENTATION_SWITCH:
         break;

      // geometric frequency at which to apply rational reconstruction
      case SoPlex::RATREC_FREQ:
         break;

      // minimal reduction (sum of removed rows/cols) to continue simplification
      case SoPlex::MINRED:
         break;

      default:
         return false;
      }

      _currentSettings->_realParamValues[param] = value;
      return true;
   }



#ifdef SOPLEX_WITH_RATIONALPARAM
   /// sets rational parameter value; returns true on success
   bool SoPlex::setRationalParam(const RationalParam param, const Rational value, const bool quiet, const bool init)
   {
      assert(param >= 0);
      assert(param < RATIONALPARAM_COUNT);
      assert(init || _isConsistent());

      if( !init && value == rationalParam(param) )
         return true;

      if( value < _currentSettings->_rationalParamLower[param] || value > _currentSettings->_rationalParamUpper[param] )
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
   bool SoPlex::setSettings(const Settings& newSettings, const bool quiet, const bool init)
   {
      assert(init || _isConsistent());

      bool success = true;

      *_currentSettings = newSettings;

      for( int i = 0; i < SoPlex::BOOLPARAM_COUNT; i++ )
         success &= setBoolParam((BoolParam)i, _currentSettings->_boolParamValues[i], quiet, init);

      for( int i = 0; i < SoPlex::INTPARAM_COUNT; i++ )
         success &= setIntParam((IntParam)i, _currentSettings->_intParamValues[i], quiet, init);

      for( int i = 0; i < SoPlex::REALPARAM_COUNT; i++ )
         success &= setRealParam((RealParam)i, _currentSettings->_realParamValues[i], quiet, init);

#ifdef SOPLEX_WITH_RATIONALPARAM
      for( int i = 0; i < SoPlex::RATIONALPARAM_COUNT; i++ )
         success &= setRationalParam((RationalParam)i, _currentSettings->_rationalParamValues[i], quiet, init);
#endif

      assert(_isConsistent());

      return success;
   }



   /// print non-default parameter values
   void SoPlex::printUserSettings()
   {
      bool printedValue = false;

      for( int i = 0; i < SoPlex::BOOLPARAM_COUNT; i++ )
      {
         if( _currentSettings->_boolParamValues[i] == _currentSettings->_boolParamDefault[i] )
            continue;

         spxout << "bool:" << _currentSettings->_boolParamName[i] << " = " << (_currentSettings->_boolParamValues[i] ? "true\n" : "false\n");
         printedValue = true;
      }

      for( int i = 0; i < SoPlex::INTPARAM_COUNT; i++ )
      {
         if( _currentSettings->_intParamValues[i] == _currentSettings->_intParamDefault[i] )
            continue;

         spxout << "int:" << _currentSettings->_intParamName[i] << " = " << _currentSettings->_intParamValues[i] << "\n";
         printedValue = true;
      }

      for( int i = 0; i < SoPlex::REALPARAM_COUNT; i++ )
      {
         if( _currentSettings->_realParamValues[i] == _currentSettings->_realParamDefault[i] )
            continue;

         spxout << "real:" << _currentSettings->_realParamName[i] << " = " << _currentSettings->_realParamValues[i] << "\n";
         printedValue = true;
      }

#ifdef SOPLEX_WITH_RATIONALPARAM
      for( int i = 0; i < SoPlex::RATIONALPARAM_COUNT; i++ )
      {
         if( _currentSettings->_rationalParamValues[i] == _currentSettings->_rationalParamDefault[i] )
            continue;

         spxout << "rational:" << _currentSettings->_rationalParamName[i] << " = " << _currentSettings->_rationalParamValues[i] << "\n";
         printedValue = true;
      }
#endif
      if( printedValue )
         spxout << std::endl;
   }



   /// writes settings file; returns true on success
   bool SoPlex::saveSettingsFile(const char* filename, const bool onlyChanged) const
   {
      assert(filename != 0);

      std::ofstream file(filename);
      if( !file.good() )
         return false;

      file.setf(std::ios::left);
      file << "# SoPlex version " << SOPLEX_VERSION / 100 << "." << (SOPLEX_VERSION / 10) % 10 << "." << SOPLEX_VERSION % 10 << "." << SOPLEX_SUBVERSION << "\n";

      for( int i = 0; i < SoPlex::BOOLPARAM_COUNT; i++ )
      {
         if( onlyChanged && _currentSettings->_boolParamValues[i] == _currentSettings->_boolParamDefault[i] )
            continue;

         file << "\n";
         file << "# " << _currentSettings->_boolParamDescription[i] << "\n";
         file << "# range {true, false}, default " << (_currentSettings->_boolParamDefault[i] ? "true\n" : "false\n");
         file << "bool:" << _currentSettings->_boolParamName[i] << " = " << (_currentSettings->_boolParamValues[i] ? "true\n" : "false\n");
      }

      for( int i = 0; i < SoPlex::INTPARAM_COUNT; i++ )
      {
         if( onlyChanged && _currentSettings->_intParamValues[i] == _currentSettings->_intParamDefault[i] )
            continue;

         file << "\n";
         file << "# " << _currentSettings->_intParamDescription[i] << "\n";
         file << "# range [-2147483648,2147483647], default " << _currentSettings->_intParamDefault[i] << "\n";
         file << "int:" << _currentSettings->_intParamName[i] << " = " << _currentSettings->_intParamValues[i] << "\n";
      }

      for( int i = 0; i < SoPlex::REALPARAM_COUNT; i++ )
      {
         if( onlyChanged && _currentSettings->_realParamValues[i] == _currentSettings->_realParamDefault[i] )
            continue;

         file << "\n";
         file << "# " << _currentSettings->_realParamDescription[i] << "\n";
         file << "# range [" << _currentSettings->_realParamLower[i] << "," << _currentSettings->_realParamUpper[i]
            << "], default " << _currentSettings->_realParamDefault[i] << "\n";
         file << "real:" << _currentSettings->_realParamName[i] << " = " << _currentSettings->_realParamValues[i] << "\n";
      }

#ifdef SOPLEX_WITH_RATIONALPARAM
      for( int i = 0; i < SoPlex::RATIONALPARAM_COUNT; i++ )
      {
         if( onlyChanged && _currentSettings->_rationalParamValues[i] == _currentSettings->_rationalParamDefault[i] )
            continue;

         file << "\n";
         file << "# " << _currentSettings->_rationalParamDescription[i] << "\n";
         file << "# range [" << _currentSettings->_rationalParamLower[i] << "," << _currentSettings->_rationalParamUpper[i]
            << "], default " << _currentSettings->_rationalParamDefault[i] << "\n";
         file << "rational:" << _currentSettings->_rationalParamName[i] << " = " << _currentSettings->_rationalParamValues[i] << "\n";
      }
#endif

      return true;
   }



   /// reads settings file; returns true on success
   bool SoPlex::loadSettingsFile(const char* filename)
   {
      assert(filename != 0);

      // start timing
      _statistics->readingTime->start();

      MSG_INFO1( spxout, spxout << "Loading settings file <" << filename << "> . . .\n" );

      // open file
      spxifstream file(filename);

      if( !file )
      {
         MSG_ERROR( std::cerr << "Error opening settings file.\n" );
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
         MSG_ERROR( std::cerr << "Error reading settings file: line " << lineNumber << " in settings file exceeds " << SET_MAX_LINE_LEN - 2 << " characters.\n" );
      }
      else if( readError )
      {
         MSG_ERROR( std::cerr << "Error reading settings file: line " << lineNumber << ".\n" );
      }

      // stop timing
      _statistics->readingTime->stop();

      return !readError && !parseError;
   }

   /// parses one setting string and returns true on success
   bool SoPlex::parseSettingsString(char* string)
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
            MSG_ERROR( std::cerr << "Error parsing setting string: no ':' separating parameter type and name.\n" );
            return false;
         }
         line++;
      }

      // find the start of the parameter name
      while( *line == ' ' || *line == '\t' || *line == '\r' )
         line++;
      if( *line == '\0' || *line == '\n' || *line == '#' )
      {
         MSG_ERROR( std::cerr << "Error parsing setting string: no parameter name.\n");
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
            MSG_ERROR( std::cerr << "Error parsing setting string: no '=' after parameter name.\n" );
            return false;
         }
         line++;
      }

      // find the start of the parameter value string
      while( *line == ' ' || *line == '\t' || *line == '\r' )
         line++;
      if( *line == '\0' || *line == '\n' || *line == '#' )
      {
         MSG_ERROR( std::cerr << "Error parsing setting string: no parameter value.\n");
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
            MSG_ERROR( std::cerr << "Error parsing setting string: additional character '" << *line << "' after parameter value.\n" );
            return false;
         }
      }

      // check whether we have a bool parameter
      if( strncmp(paramTypeString, "bool", 4) == 0 )
      {
         for( int param = 0; ; param++ )
         {
            if( param >= SoPlex::BOOLPARAM_COUNT )
            {
               MSG_ERROR( std::cerr << "Error parsing setting string: unknown parameter name <" << paramName << ">.\n" );
               return false;
            }
            else if( strncmp(paramName, _currentSettings->_boolParamName[param].c_str(), SET_MAX_LINE_LEN) == 0 )
            {
               if( strncasecmp(paramValueString, "true", 4) == 0
                  || strncasecmp(paramValueString, "TRUE", 4) == 0
                  || strncasecmp(paramValueString, "t", 4) == 0
                  || strncasecmp(paramValueString, "T", 4) == 0
                  || strtol(paramValueString, NULL, 4) == 1 )
               {
                  setBoolParam((SoPlex::BoolParam)param, true);
                  break;
               }
               else if( strncasecmp(paramValueString, "false", 5) == 0
                  || strncasecmp(paramValueString, "FALSE", 5) == 0
                  || strncasecmp(paramValueString, "f", 5) == 0
                  || strncasecmp(paramValueString, "F", 5) == 0
                  || strtol(paramValueString, NULL, 5) == 0 )
               {
                  setBoolParam((SoPlex::BoolParam)param, false);
                  break;
               }
               else
               {
                  MSG_ERROR( std::cerr << "Error parsing setting string: invalid value <" << paramValueString << "> for bool parameter <" << paramName << ">.\n" );
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
            if( param >= SoPlex::INTPARAM_COUNT )
            {
               MSG_ERROR( std::cerr << "Error parsing setting string: unknown parameter name <" << paramName << ">.\n" );
               return false;
            }
            else if( strncmp(paramName, _currentSettings->_intParamName[param].c_str(), SET_MAX_LINE_LEN) == 0 )
            {
               int value;

               if( sscanf(paramValueString, "%d", &value) == 1 && setIntParam((SoPlex::IntParam)param, value) )
                  break;
               else
               {
                  MSG_ERROR( std::cerr << "Error parsing setting string: invalid value <" << paramValueString << "> for int parameter <" << paramName << ">.\n" );
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
            if( param >= SoPlex::REALPARAM_COUNT )
            {
               MSG_ERROR( std::cerr << "Error parsing setting string: unknown parameter name <" << paramName << ">.\n" );
               return false;
            }
            else if( strncmp(paramName, _currentSettings->_realParamName[param].c_str(), SET_MAX_LINE_LEN) == 0 )
            {
               Real value;

               if( sscanf(paramValueString, "%" REAL_FORMAT, &value) == 1 && setRealParam((SoPlex::RealParam)param, value) )
                  break;
               else
               {
                  MSG_ERROR( std::cerr << "Error parsing setting string: invalid value <" << paramValueString << "> for real parameter <" << paramName << ">.\n" );
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
            if( param >= SoPlex::RATIONALPARAM_COUNT )
            {
               MSG_ERROR( std::cerr << "Error parsing setting string: unknown parameter name <" << paramName << ">.\n" );
               return false;
            }
            else if( strncmp(paramName, _currentSettings->_rationalParamName[param].c_str(), SET_MAX_LINE_LEN) == 0 )
            {
               Rational value;

               if( readStringRational(paramValueString, value) && setRationalParam((SoPlex::RationalParam)param, value) )
                  break;
               else
               {
                  MSG_ERROR( std::cerr << "Error parsing setting string: invalid value <" << paramValueString << "> for rational parameter <" << paramName << ">.\n" );
                  return false;
               }
            }
         }

         return true;
      }
#endif

      MSG_ERROR( std::cerr << "Error parsing setting string: invalid parameter type <" << paramTypeString << "> for parameter <" << paramName << ">.\n" );

      return false;
   }




   /// prints solution statistics
   void SoPlex::printSolutionStatistics(std::ostream& os)
   {
      if( _lastSolveMode == SOLVEMODE_REAL )
      {
         os << std::scientific << std::setprecision(8)
            << "Solution (real)     : \n"
            << "  Objective Value   : " << objValueReal() << "\n";
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

      if( intParam(SoPlex::CHECKMODE) == CHECKMODE_RATIONAL
         || (intParam(SoPlex::CHECKMODE) == CHECKMODE_AUTO && intParam(SoPlex::READMODE) == READMODE_RATIONAL) )
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
         if( getBoundViolationReal(maxviol, sumviol) )
            os << "  Max/sum bound     : " << rationalToString(maxviol) << " / " << rationalToString(sumviol) << "\n";
         else
            os << "  Max/sum bound     : - / -\n";
         if( getRowViolationReal(maxviol, sumviol) )
            os << "  Max/sum row       : " << rationalToString(maxviol) << " / " << rationalToString(sumviol) << "\n";
         else
            os << "  Max/sum row       : - / -\n";
         if( getRedCostViolationReal(maxviol, sumviol) )
            os << "  Max/sum redcost   : " << rationalToString(maxviol) << " / " << rationalToString(sumviol) << "\n";
         else
            os << "  Max/sum redcost   : - / -\n";
         if( getDualViolationReal(maxviol, sumviol) )
            os << "  Max/sum dual      : " << rationalToString(maxviol) << " / " << rationalToString(sumviol) << "\n";
         else
            os << "  Max/sum dual      : - / -\n";
      }
   }



   /// prints statistics on solving process
   void SoPlex::printSolvingStatistics(std::ostream& os)
   {
      assert(_statistics != 0);
      _statistics->print(os);
   }



   /// prints short statistics
   void SoPlex::printShortStatistics(std::ostream& os)
   {
      printStatus(os, _status);
      os << "Solving time (sec)  : " << std::fixed << std::setprecision(2) << _statistics->solvingTime->time() << "\n"
         << "Iterations          : " << _statistics->iterations << "\n"
         << "Objective value     : " << std::scientific << std::setprecision(8) << objValueReal() << std::fixed << "\n";
   }



   /// prints complete statistics
   void SoPlex::printStatistics(std::ostream& os)
   {
      os << std::setprecision(2);

      printStatus(os, _status);

      os << "Original problem    : \n";
      if( intParam(SoPlex::READMODE) == READMODE_REAL )
         _realLP->printProblemStatistics(os);
      else
         _rationalLP->printProblemStatistics(os);

      os << "Objective sense     : " << (intParam(SoPlex::OBJSENSE) == SoPlex::OBJSENSE_MINIMIZE ? "minimize\n" : "maximize\n");
      printSolutionStatistics(os);
      printSolvingStatistics(os);
   }



   /// prints status
   void SoPlex::printStatus(std::ostream& os, SPxSolver::Status stat)
   {
      os << "SoPlex status       : ";

      switch( stat )
      {
      case SPxSolver::ERROR:
         os << "error [unspecified]";
         break;
      case SPxSolver::NO_RATIOTESTER:
         os << "error [no ratiotester loaded]";
         break;
      case SPxSolver::NO_PRICER:
         os << "error [no pricer loaded]";
         break;
      case SPxSolver::NO_SOLVER:
         os << "error [no linear solver loaded]";
         break;
      case SPxSolver::NOT_INIT:
         os << "error [not initialized]";
         break;
      case SPxSolver::ABORT_CYCLING:
         os << "solving aborted [cycling]";
         break;
      case SPxSolver::ABORT_TIME:
         os << "solving aborted [time limit reached]";
         break;
      case SPxSolver::ABORT_ITER:
         os << "solving aborted [iteration limit reached]";
         break;
      case SPxSolver::ABORT_VALUE:
         os << "solving aborted [objective limit reached]";
         break;
      case SPxSolver::NO_PROBLEM:
         os << "no problem loaded";
         break;
      case SPxSolver::REGULAR:
         os << "basis is regular";
         break;
      case SPxSolver::SINGULAR:
         os << "basis is singular";
         break;
      case SPxSolver::OPTIMAL:
         os << "problem is solved [optimal]";
         break;
      case SPxSolver::UNBOUNDED:
         os << "problem is solved [unbounded]";
         break;
      case SPxSolver::INFEASIBLE:
         os << "problem is solved [infeasible]";
         break;
      case SPxSolver::INForUNBD:
         os << "problem is solved [infeasible or unbounded]";
         break;
      default:
      case SPxSolver::UNKNOWN:
         os << "unknown";
         break;
      }

      os << "\n";
   }



   /// prints version and compilation options
   void SoPlex::printVersion() const
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
   bool SoPlex::areLPsInSync(const bool checkVecVals, const bool checkMatVals, const bool quiet) const
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
         MSG_ERROR( std::cerr << "The number of Rows in the Real LP does not match the one in the Rational LP."
               << " Real LP: " << _realLP->nRows() << "  Rational LP: " << _rationalLP->nRows() << std::endl);
         result = false;
         nRowsMatch = false;
      }

      // compare number of Columns
      if( _realLP->nCols() != _rationalLP->nCols() )
      {
         MSG_ERROR( std::cerr << "The number of Columns in the Real LP does not match the one in the Rational LP."
               << " Real LP: " << _realLP->nCols() << "  Rational LP: " << _rationalLP->nCols() << std::endl);
         result = false;
         nColsMatch = false;
      }

      // compare number of nonZeros
      if( _realLP->nNzos() != _rationalLP->nNzos() )
      {
         MSG_ERROR( std::cerr << "The number of nonZeros in the Real LP does not match the one in the Rational LP."
               << " Real LP: " << _realLP->nNzos() << "  Rational LP: " << _rationalLP->nNzos() << std::endl);
         result = false;
      }

      // compare the dimensions of the right hand side vectors
      if( _realLP->rhs().dim() != _rationalLP->rhs().dim() )
      {
         MSG_ERROR( std::cerr << "The dimension of the right hand side vector of the Real LP does not match the one of the Rational LP."
               << " Real LP: " << _realLP->rhs().dim() << "  Rational LP: " << _rationalLP->rhs().dim() << std::endl);
         result = false;
         rhsDimMatch = false;

      }

      // compare the dimensions of the left hand side vectors
      if( _realLP->lhs().dim() != _rationalLP->lhs().dim() )
      {
         MSG_ERROR( std::cerr << "The dimension of the left hand side vector of the Real LP does not match the one of the Rational LP."
               << " Real LP: " << _realLP->lhs().dim() << "  Rational LP: " << _rationalLP->lhs().dim() << std::endl);
         result = false;
         lhsDimMatch = false;
      }

      // compare the dimensions of the objective function vectors
      if( _realLP->maxObj().dim() != _rationalLP->maxObj().dim() )
      {
         MSG_ERROR( std::cerr << "The dimension of the objective function vector of the Real LP does not match the one of the Rational LP."
               << " Real LP: " << _realLP->maxObj().dim() << "  Rational LP: " << _rationalLP->maxObj().dim() << std::endl);
         result = false;
         maxObjDimMatch = false;
      }

      // compare the sense
      if( (int)_realLP->spxSense() != (int)_rationalLP->spxSense() )
         {
            MSG_ERROR( std::cerr << "The objective function sense of the Real LP does not match the one of the Rational LP."
               << " Real LP: " << (_realLP->spxSense() == SPxLPReal::MINIMIZE ? "MIN" : "MAX")
               << "  Rational LP: " << (_rationalLP->spxSense() == SPxLPRational::MINIMIZE ? "MIN" : "MAX") << std::endl);
            result = false;
         }

      // compare the dimensions of upper bound vectors
      if( _realLP->upper().dim() != _rationalLP->upper().dim() )
      {
         MSG_ERROR( std::cerr << "The dimension of the upper bound vector of the Real LP does not match the one of the Rational LP."
               << " Real LP: " << _realLP->upper().dim() << "  Rational LP: " << _rationalLP->upper().dim() << std::endl);
         result = false;
         upperDimMatch = false;
      }

      // compare the dimensions of the objective function vectors
      if( _realLP->lower().dim() != _rationalLP->lower().dim() )
      {
         MSG_ERROR( std::cerr << "The dimension of the lower bound vector of the Real LP does not match the one of the Rational LP."
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
               if( (GE(_realLP->rhs()[i], realParam(SoPlex::INFTY)) != (_rationalLP->rhs()[i] >= _rationalPosInfty))
                  || (LT(_realLP->rhs()[i], realParam(SoPlex::INFTY)) && _rationalLP->rhs()[i] < _rationalPosInfty
                     && !_rationalLP->rhs()[i].isAdjacentTo((double)_realLP->rhs()[i])) )
               {
                  if( !quiet )
                  {
                     MSG_ERROR( std::cerr << "Entries number " << i << " of the right hand side vectors don't match."
                           << " Real LP: " << _realLP->rhs()[i] << "  Rational LP: " << _rationalLP->rhs()[i] << std::endl);
                  }
                  rhsValMatch = false;
                  result = false;
               }
            }

            if( !rhsValMatch && quiet )
            {
               MSG_ERROR( std::cerr << "The values of the right hand side vectors don't match." << std::endl );
            }
         }

         // compares the values of the left hand side vectors
         if( lhsDimMatch )
         {
            for( int i = 0; i < _realLP->lhs().dim(); i++ )
            {
               if( (LE(_realLP->lhs()[i], -realParam(SoPlex::INFTY)) != (_rationalLP->lhs()[i] <= _rationalNegInfty))
                  || (GT(_realLP->lhs()[i], -realParam(SoPlex::INFTY)) && _rationalLP->lhs()[i] > _rationalNegInfty
                     && !_rationalLP->lhs()[i].isAdjacentTo((double)_realLP->lhs()[i])) )
               {
                  if( !quiet )
                  {
                     MSG_ERROR( std::cerr << "Entries number " << i << " of the left hand side vectors don't match."
                           << " Real LP: " << _realLP->lhs()[i] << "  Rational LP: " << _rationalLP->lhs()[i] << std::endl);
                  }
                  lhsValMatch = false;
                  result = false;
               }
            }

            if( !lhsValMatch && quiet )
            {
               MSG_ERROR( std::cerr << "The values of the left hand side vectors don't match." << std::endl );
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
                     MSG_ERROR( std::cerr << "Entries number " << i << " of the objective function vectors don't match."
                           << " Real LP: " << _realLP->maxObj()[i] << "  Rational LP: " << _rationalLP->maxObj()[i] << std::endl);
                  }
                  maxObjValMatch = false;
                  result = false;
               }
            }

            if( !maxObjValMatch && quiet )
            {
               MSG_ERROR( std::cerr << "The values of the objective function vectors don't match." << std::endl );
            }
         }

         // compares the values of the upper bound vectors
         if( upperDimMatch )
         {
            for( int i = 0; i < _realLP->upper().dim(); i++ )
            {
               if( (GE(_realLP->upper()[i], realParam(SoPlex::INFTY)) != (_rationalLP->upper()[i] >= _rationalPosInfty))
                  || (LT(_realLP->upper()[i], realParam(SoPlex::INFTY)) && _rationalLP->upper()[i] < _rationalPosInfty
                     && !_rationalLP->upper()[i].isAdjacentTo((double)_realLP->upper()[i])) )
               {
                  if( !quiet )
                  {
                     MSG_ERROR( std::cerr << "Entries number " << i << " of the upper bound vectors don't match."
                           << " Real LP: " << _realLP->upper()[i] << "  Rational LP: " << _rationalLP->upper()[i] << std::endl);
                  }
                  upperValMatch = false;
                  result = false;
               }
            }

            if( !upperValMatch && quiet )
            {
               MSG_ERROR( std::cerr << "The values of the upper bound vectors don't match." << std::endl );
            }
         }

         // compares the values of the lower bound vectors
         if( lowerDimMatch )
         {
            for( int i = 0; i < _realLP->lower().dim(); i++ )
            {
               if( (LE(_realLP->lower()[i], -realParam(SoPlex::INFTY)) != (_rationalLP->lower()[i] <= _rationalNegInfty))
                  || (GT(_realLP->lower()[i], -realParam(SoPlex::INFTY)) && _rationalLP->lower()[i] > _rationalNegInfty
                     && !_rationalLP->lower()[i].isAdjacentTo((double)_realLP->lower()[i])) )
               {
                  if( !quiet )
                  {
                     MSG_ERROR( std::cerr << "Entries number " << i << " of the lower bound vectors don't match."
                           << " Real LP: " << _realLP->lower()[i] << "  Rational LP: " << _rationalLP->lower()[i] << std::endl);
                  }
                  lowerValMatch = false;
                  result = false;
               }
            }

            if( !lowerValMatch && quiet )
            {
               MSG_ERROR( std::cerr << "The values of the lower bound vectors don't match." << std::endl );
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
                     MSG_ERROR( std::cerr << "Entries number " << j << " of column number " << i << " don't match."
                           << " Real LP: " << _realLP->colVector(i)[j] << "  Rational LP: " << _rationalLP->colVector(i)[j] << std::endl);
                  }
                  matrixValMatch = false;
                  result = false;
               }
            }
         }

         if( !matrixValMatch && quiet )
         {
            MSG_ERROR( std::cerr << "The values of the matrices don't match." << std::endl );
         }
      }

      return result;
   }



   /// extends sparse vector to hold newmax entries if and only if it holds no more free entries
   void SoPlex::_ensureDSVectorRationalMemory(DSVectorRational& vec, const int newmax) const
   {
      assert(newmax > vec.size());
      if( vec.size() >= vec.max() )
         vec.setMax(newmax);
   }



   /// creates a permutation for removing rows/columns from an array of indices
   void SoPlex::_idxToPerm(int* idx, int idxSize, int* perm, int permSize) const
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
   void SoPlex::_rangeToPerm(int start, int end, int* perm, int permSize) const
   {
      assert(perm != 0);
      assert(permSize >= 0);

      for( int i = 0; i < permSize; i++ )
         perm[i] = (i < start || i > end) ? i : -1;
   }



   /// checks consistency
   bool SoPlex::_isConsistent() const
   {
      assert(_statistics != 0);
      assert(_currentSettings != 0);

      assert(_realLP != 0);
      assert(_rationalLP != 0 || intParam(SoPlex::SYNCMODE) == SYNCMODE_ONLYREAL);

      assert(_realLP != &_solver || _isRealLPLoaded);
      assert(_realLP == &_solver || !_isRealLPLoaded);

      assert(!_hasBasis || _isRealLPLoaded || _basisStatusRows.size() == numRowsReal());
      assert(!_hasBasis || _isRealLPLoaded || _basisStatusCols.size() == numColsReal());

      return true;
   }



   /// should solving process be stopped?
   bool SoPlex::_isSolveStopped() const
   {
      assert(_statistics != 0);

      return (realParam(TIMELIMIT) < realParam(INFTY) && _statistics->solvingTime->time() >= realParam(TIMELIMIT))
         || (intParam(ITERLIMIT) >= 0 && _statistics->iterations >= intParam(ITERLIMIT))
         || (intParam(REFLIMIT) >= 0 && _statistics->refinements >= intParam(REFLIMIT))
         || (intParam(STALLREFLIMIT) >= 0 && _statistics->stallRefinements >= intParam(STALLREFLIMIT));
   }



   /// determines RangeType from real bounds
   SoPlex::RangeType SoPlex::_rangeTypeReal(const Real& lower, const Real& upper) const
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
   SoPlex::RangeType SoPlex::_rangeTypeRational(const Rational& lower, const Rational& upper) const
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
   SoPlex::RangeType SoPlex::_switchRangeType(const SoPlex::RangeType& rangeType) const
   {
      if( rangeType == RANGETYPE_LOWER )
         return RANGETYPE_UPPER;
      else if( rangeType == RANGETYPE_UPPER )
         return RANGETYPE_LOWER;
      else
         return rangeType;
   }



   /// checks whether RangeType corresponds to finite lower bound
   bool SoPlex::_lowerFinite(const RangeType& rangeType) const
   {
      return (rangeType == RANGETYPE_LOWER || rangeType == RANGETYPE_BOXED || rangeType == RANGETYPE_FIXED);
   }



   /// checks whether RangeType corresponds to finite upper bound
   bool SoPlex::_upperFinite(const RangeType& rangeType) const
   {
      return (rangeType == RANGETYPE_UPPER || rangeType == RANGETYPE_BOXED || rangeType == RANGETYPE_FIXED);
   }



   /// adds a single row to the real LP and adjusts basis
   void SoPlex::_addRowReal(const LPRowReal& lprow)
   {
      assert(_realLP != 0);

      _realLP->addRow(lprow);

      if( _isRealLPLoaded )
         _hasBasis = (_solver.basis().status() > SPxBasis::NO_PROBLEM);
      else if( _hasBasis )
         _basisStatusRows.append(SPxSolver::BASIC);
   }



   /// adds a single row to the real LP and adjusts basis
   void SoPlex::_addRowReal(Real lhs, const SVectorReal& lprow, Real rhs)
   {
      assert(_realLP != 0);

      _realLP->addRow(lhs, lprow, rhs);

      if( _isRealLPLoaded )
         _hasBasis = (_solver.basis().status() > SPxBasis::NO_PROBLEM);
      else if( _hasBasis )
         _basisStatusRows.append(SPxSolver::BASIC);
   }



   /// adds multiple rows to the real LP and adjusts basis
   void SoPlex::_addRowsReal(const LPRowSetReal& lprowset)
   {
      assert(_realLP != 0);

      _realLP->addRows(lprowset);

      if( _isRealLPLoaded )
         _hasBasis = (_solver.basis().status() > SPxBasis::NO_PROBLEM);
      else if( _hasBasis )
         _basisStatusRows.append(lprowset.num(), SPxSolver::BASIC);
   }


   /// adds a single column to the real LP and adjusts basis
   void SoPlex::_addColReal(const LPColReal& lpcol)
   {
      assert(_realLP != 0);

      _realLP->addCol(lpcol);

      if( _isRealLPLoaded )
         _hasBasis = (_solver.basis().status() > SPxBasis::NO_PROBLEM);
      else if( _hasBasis )
      {
         if( lpcol.lower() > -realParam(SoPlex::INFTY) )
            _basisStatusCols.append(SPxSolver::ON_LOWER);
         else if( lpcol.upper() < realParam(SoPlex::INFTY) )
            _basisStatusCols.append(SPxSolver::ON_UPPER);
         else
            _basisStatusCols.append(SPxSolver::ZERO);
      }
   }



   /// adds a single column to the real LP and adjusts basis
   void SoPlex::_addColReal(Real obj, Real lower, const SVectorReal& lpcol, Real upper)
   {
      assert(_realLP != 0);

      _realLP->addCol(obj, lower, lpcol, upper);

      if( _isRealLPLoaded )
         _hasBasis = (_solver.basis().status() > SPxBasis::NO_PROBLEM);
      else if( _hasBasis )
         _basisStatusRows.append(SPxSolver::BASIC);
   }



   /// adds multiple columns to the real LP and adjusts basis
   void SoPlex::_addColsReal(const LPColSetReal& lpcolset)
   {
      assert(_realLP != 0);

      _realLP->addCols(lpcolset);

      if( _isRealLPLoaded )
         _hasBasis = (_solver.basis().status() > SPxBasis::NO_PROBLEM);
      else if( _hasBasis )
      {
         for( int i = 0; i < lpcolset.num(); i++ )
         {
            if( lpcolset.lower(i) > -realParam(SoPlex::INFTY) )
               _basisStatusCols.append(SPxSolver::ON_LOWER);
            else if( lpcolset.upper(i) < realParam(SoPlex::INFTY) )
               _basisStatusCols.append(SPxSolver::ON_UPPER);
            else
               _basisStatusCols.append(SPxSolver::ZERO);
         }
      }
   }


   /// replaces row \p i with \p lprow and adjusts basis
   void SoPlex::_changeRowReal(int i, const LPRowReal& lprow)
   {
      assert(_realLP != 0);

      _realLP->changeRow(i, lprow);

      if( _isRealLPLoaded )
         _hasBasis = (_solver.basis().status() > SPxBasis::NO_PROBLEM);
      else if( _hasBasis )
      {
         if( _basisStatusRows[i] != SPxSolver::BASIC )
            _hasBasis = false;
         else if( _basisStatusRows[i] == SPxSolver::ON_LOWER && lprow.lhs() <= -realParam(SoPlex::INFTY) )
            _basisStatusRows[i] = (lprow.rhs() < realParam(SoPlex::INFTY)) ? SPxSolver::ON_UPPER : SPxSolver::ZERO;
         else if( _basisStatusRows[i] == SPxSolver::ON_UPPER && lprow.rhs() >= realParam(SoPlex::INFTY) )
            _basisStatusRows[i] = (lprow.lhs() > -realParam(SoPlex::INFTY)) ? SPxSolver::ON_LOWER : SPxSolver::ZERO;
      }
   }



   /// changes left-hand side vector for constraints to \p lhs and adjusts basis
   void SoPlex::_changeLhsReal(const VectorReal& lhs)
   {
      assert(_realLP != 0);

      _realLP->changeLhs(lhs);

      if( _isRealLPLoaded )
         _hasBasis = (_solver.basis().status() > SPxBasis::NO_PROBLEM);
      else if( _hasBasis )
      {
         for( int i = numRowsReal() - 1; i >= 0; i-- )
         {
            if( _basisStatusRows[i] == SPxSolver::ON_LOWER && lhs[i] <= -realParam(SoPlex::INFTY) )
               _basisStatusRows[i] = (rhsReal(i) < realParam(SoPlex::INFTY)) ? SPxSolver::ON_UPPER : SPxSolver::ZERO;
         }
      }
   }



   /// changes left-hand side of row \p i to \p lhs and adjusts basis
   void SoPlex::_changeLhsReal(int i, const Real& lhs)
   {
      assert(_realLP != 0);

      _realLP->changeLhs(i, lhs);

      if( _isRealLPLoaded )
      {
         _hasBasis = (_solver.basis().status() > SPxBasis::NO_PROBLEM);
      }
      else if( _hasBasis && _basisStatusRows[i] == SPxSolver::ON_LOWER && lhs <= -realParam(SoPlex::INFTY) )
         _basisStatusRows[i] = (rhsReal(i) < realParam(SoPlex::INFTY)) ? SPxSolver::ON_UPPER : SPxSolver::ZERO;

   }



   /// changes right-hand side vector to \p rhs and adjusts basis
   void SoPlex::_changeRhsReal(const VectorReal& rhs)
   {
      assert(_realLP != 0);

      _realLP->changeRhs(rhs);

      if( _isRealLPLoaded )
      {
         _hasBasis = (_solver.basis().status() > SPxBasis::NO_PROBLEM);
      }
      else if( _hasBasis )
      {
         for( int i = numRowsReal() - 1; i >= 0; i-- )
         {
            if( _basisStatusRows[i] == SPxSolver::ON_UPPER && rhs[i] >= realParam(SoPlex::INFTY) )
               _basisStatusRows[i] = (lhsReal(i) > -realParam(SoPlex::INFTY)) ? SPxSolver::ON_LOWER : SPxSolver::ZERO;
         }
      }
   }



   /// changes right-hand side of row \p i to \p rhs and adjusts basis
   void SoPlex::_changeRhsReal(int i, const Real& rhs)
   {
      assert(_realLP != 0);

      _realLP->changeRhs(i, rhs);

      if( _isRealLPLoaded )
      {
         _hasBasis = (_solver.basis().status() > SPxBasis::NO_PROBLEM);
      }
      else if( _hasBasis && _basisStatusRows[i] == SPxSolver::ON_UPPER && rhs >= realParam(SoPlex::INFTY) )
         _basisStatusRows[i] = (lhsReal(i) > -realParam(SoPlex::INFTY)) ? SPxSolver::ON_LOWER : SPxSolver::ZERO;
   }



   /// changes left- and right-hand side vectors and adjusts basis
   void SoPlex::_changeRangeReal(const VectorReal& lhs, const VectorReal& rhs)
   {
      assert(_realLP != 0);

      _realLP->changeRange(lhs, rhs);

      if( _isRealLPLoaded )
      {
         _hasBasis = (_solver.basis().status() > SPxBasis::NO_PROBLEM);
      }
      else if( _hasBasis )
      {
         for( int i = numRowsReal() - 1; i >= 0; i-- )
         {
            if( _basisStatusRows[i] == SPxSolver::ON_LOWER && lhs[i] <= -realParam(SoPlex::INFTY) )
               _basisStatusRows[i] = (rhs[i] < realParam(SoPlex::INFTY)) ? SPxSolver::ON_UPPER : SPxSolver::ZERO;
            else if( _basisStatusRows[i] == SPxSolver::ON_UPPER && rhs[i] >= realParam(SoPlex::INFTY) )
               _basisStatusRows[i] = (lhs[i] > -realParam(SoPlex::INFTY)) ? SPxSolver::ON_LOWER : SPxSolver::ZERO;
         }
      }
   }



   /// changes left- and right-hand side of row \p i and adjusts basis
   void SoPlex::_changeRangeReal(int i, const Real& lhs, const Real& rhs)
   {
      assert(_realLP != 0);

      _realLP->changeRange(i, lhs, rhs);

      if( _isRealLPLoaded )
      {
         _hasBasis = (_solver.basis().status() > SPxBasis::NO_PROBLEM);
      }
      else if( _hasBasis )
      {
         if( _basisStatusRows[i] == SPxSolver::ON_LOWER && lhs <= -realParam(SoPlex::INFTY) )
            _basisStatusRows[i] = (rhs < realParam(SoPlex::INFTY)) ? SPxSolver::ON_UPPER : SPxSolver::ZERO;
         else if( _basisStatusRows[i] == SPxSolver::ON_UPPER && rhs >= realParam(SoPlex::INFTY) )
            _basisStatusRows[i] = (lhs > -realParam(SoPlex::INFTY)) ? SPxSolver::ON_LOWER : SPxSolver::ZERO;
      }
   }



   /// replaces column \p i with \p lpcol and adjusts basis
   void SoPlex::_changeColReal(int i, const LPColReal& lpcol)
   {
      assert(_realLP != 0);

      _realLP->changeCol(i, lpcol);

      if( _isRealLPLoaded )
      {
         _hasBasis = (_solver.basis().status() > SPxBasis::NO_PROBLEM);
      }
      else if( _hasBasis )
      {
         if( _basisStatusCols[i] == SPxSolver::BASIC )
            _hasBasis = false;
         else if( _basisStatusCols[i] == SPxSolver::ON_LOWER && lpcol.lower() <= -realParam(SoPlex::INFTY) )
            _basisStatusCols[i] = (lpcol.upper() < realParam(SoPlex::INFTY)) ? SPxSolver::ON_UPPER : SPxSolver::ZERO;
         else if( _basisStatusCols[i] == SPxSolver::ON_UPPER && lpcol.upper() >= realParam(SoPlex::INFTY) )
            _basisStatusCols[i] = (lpcol.lower() > -realParam(SoPlex::INFTY)) ? SPxSolver::ON_LOWER : SPxSolver::ZERO;
      }
   }



   /// changes vector of lower bounds to \p lower and adjusts basis
   void SoPlex::_changeLowerReal(const VectorReal& lower)
   {
      assert(_realLP != 0);

      _realLP->changeLower(lower);

      if( _isRealLPLoaded )
      {
         _hasBasis = (_solver.basis().status() > SPxBasis::NO_PROBLEM);
      }
      else if( _hasBasis )
      {
         for( int i = numColsReal() - 1; i >= 0; i-- )
         {
            if( _basisStatusCols[i] == SPxSolver::ON_LOWER && lower[i] <= -realParam(SoPlex::INFTY) )
               _basisStatusCols[i] = (upperReal(i) < realParam(SoPlex::INFTY)) ? SPxSolver::ON_UPPER : SPxSolver::ZERO;
         }
      }
   }



   /// changes lower bound of column i to \p lower and adjusts basis
   void SoPlex::_changeLowerReal(int i, const Real& lower)
   {
      assert(_realLP != 0);

      _realLP->changeLower(i, lower);

      if( _isRealLPLoaded )
      {
         _hasBasis = (_solver.basis().status() > SPxBasis::NO_PROBLEM);
      }
      else if( _hasBasis && _basisStatusCols[i] == SPxSolver::ON_LOWER && lower <= -realParam(SoPlex::INFTY) )
         _basisStatusCols[i] = (upperReal(i) < realParam(SoPlex::INFTY)) ? SPxSolver::ON_UPPER : SPxSolver::ZERO;
   }



   /// changes vector of upper bounds to \p upper and adjusts basis
   void SoPlex::_changeUpperReal(const VectorReal& upper)
   {
      assert(_realLP != 0);

      _realLP->changeUpper(upper);

      if( _isRealLPLoaded )
      {
         _hasBasis = (_solver.basis().status() > SPxBasis::NO_PROBLEM);
      }
      else if( _hasBasis )
      {
         for( int i = numColsReal() - 1; i >= 0; i-- )
         {
            if( _basisStatusCols[i] == SPxSolver::ON_UPPER && upper[i] >= realParam(SoPlex::INFTY) )
               _basisStatusCols[i] = (lowerReal(i) > -realParam(SoPlex::INFTY)) ? SPxSolver::ON_LOWER : SPxSolver::ZERO;
         }
      }
   }



   /// changes \p i 'th upper bound to \p upper and adjusts basis
   void SoPlex::_changeUpperReal(int i, const Real& upper)
   {
      assert(_realLP != 0);

      _realLP->changeUpper(i, upper);

      if( _isRealLPLoaded )
      {
         _hasBasis = (_solver.basis().status() > SPxBasis::NO_PROBLEM);
      }
      else if( _hasBasis &&  _basisStatusCols[i] == SPxSolver::ON_UPPER && upper >= realParam(SoPlex::INFTY) )
         _basisStatusCols[i] = (lowerReal(i) > -realParam(SoPlex::INFTY)) ? SPxSolver::ON_LOWER : SPxSolver::ZERO;
   }



   /// changes vectors of column bounds to \p lower and \p upper and adjusts basis
   void SoPlex::_changeBoundsReal(const VectorReal& lower, const VectorReal& upper)
   {
      assert(_realLP != 0);

      _realLP->changeBounds(lower, upper);

      if( _isRealLPLoaded )
      {
         _hasBasis = (_solver.basis().status() > SPxBasis::NO_PROBLEM);
      }
      else if( _hasBasis )
      {
         for( int i = numColsReal() - 1; i >= 0; i-- )
         {
            if( _basisStatusCols[i] == SPxSolver::ON_LOWER && lower[i] <= -realParam(SoPlex::INFTY) )
               _basisStatusCols[i] = (upper[i] < realParam(SoPlex::INFTY)) ? SPxSolver::ON_UPPER : SPxSolver::ZERO;
            else if( _basisStatusCols[i] == SPxSolver::ON_UPPER && upper[i] >= realParam(SoPlex::INFTY) )
               _basisStatusCols[i] = (lower[i] > -realParam(SoPlex::INFTY)) ? SPxSolver::ON_LOWER : SPxSolver::ZERO;
         }
      }
   }



   /// changes bounds of column \p i to \p lower and \p upper and adjusts basis
   void SoPlex::_changeBoundsReal(int i, const Real& lower, const Real& upper)
   {
      assert(_realLP != 0);

      _realLP->changeBounds(i, lower, upper);

      if( _isRealLPLoaded )
      {
         _hasBasis = (_solver.basis().status() > SPxBasis::NO_PROBLEM);
      }
      else if( _hasBasis )
      {
         if( _basisStatusCols[i] == SPxSolver::ON_LOWER && lower <= -realParam(SoPlex::INFTY) )
            _basisStatusCols[i] = (upper < realParam(SoPlex::INFTY)) ? SPxSolver::ON_UPPER : SPxSolver::ZERO;
         else if( _basisStatusCols[i] == SPxSolver::ON_UPPER && upper >= realParam(SoPlex::INFTY) )
            _basisStatusCols[i] = (lower > -realParam(SoPlex::INFTY)) ? SPxSolver::ON_LOWER : SPxSolver::ZERO;
      }
   }



   /// changes matrix entry in row \p i and column \p j to \p val and adjusts basis
   void SoPlex::_changeElementReal(int i, int j, const Real& val)
   {
      assert(_realLP != 0);

      _realLP->changeElement(i, j, val);

      if( _isRealLPLoaded )
      {
         _hasBasis = (_solver.basis().status() > SPxBasis::NO_PROBLEM);
      }
      else if( _hasBasis )
      {
         if( _basisStatusRows[i] != SPxSolver::BASIC && _basisStatusCols[i] == SPxSolver::BASIC )
            _hasBasis = false;
      }
   }



   /// removes row \p i and adjusts basis
   void SoPlex::_removeRowReal(int i)
   {
      assert(_realLP != 0);

      _realLP->removeRow(i);

      if( _isRealLPLoaded )
      {
         _hasBasis = (_solver.basis().status() > SPxBasis::NO_PROBLEM);
      }
      else if( _hasBasis )
      {
         if( _basisStatusRows[i] != SPxSolver::BASIC )
            _hasBasis = false;
         else
         {
            _basisStatusRows[i] = _basisStatusRows[_basisStatusRows.size() - 1];
            _basisStatusRows.removeLast();
         }
      }
   }



   /// removes all rows with an index \p i such that \p perm[i] < 0; upon completion, \p perm[i] >= 0 indicates the
   /// new index where row \p i has been moved to; note that \p perm must point to an array of size at least
   /// #numRowsReal()
   void SoPlex::_removeRowsReal(int perm[])
   {
      assert(_realLP != 0);

      _realLP->removeRows(perm);

      if( _isRealLPLoaded )
      {
         _hasBasis = (_solver.basis().status() > SPxBasis::NO_PROBLEM);
      }
      else if( _hasBasis )
      {
         for( int i = numRowsReal() - 1; i >= 0 && _hasBasis; i-- )
         {
            if( perm[i] < 0 && _basisStatusRows[i] != SPxSolver::BASIC )
               _hasBasis = false;
            else if( perm[i] >= 0 && perm[i] != i )
            {
               assert(perm[i] < numRowsReal());
               assert(perm[perm[i]] < 0);

               _basisStatusRows[perm[i]] = _basisStatusRows[i];
            }
         }

         if( _hasBasis )
            _basisStatusRows.reSize(numRowsReal());
      }
   }



   /// removes column i
   void SoPlex::_removeColReal(int i)
   {
      assert(_realLP != 0);

      _realLP->removeCol(i);

      if( _isRealLPLoaded )
      {
         _hasBasis = (_solver.basis().status() > SPxBasis::NO_PROBLEM);
      }
      else if( _hasBasis )
      {
         if( _basisStatusCols[i] == SPxSolver::BASIC )
            _hasBasis = false;
         else
         {
            _basisStatusCols[i] = _basisStatusCols[_basisStatusCols.size() - 1];
            _basisStatusCols.removeLast();
         }
      }
   }



   /// removes all columns with an index \p i such that \p perm[i] < 0; upon completion, \p perm[i] >= 0 indicates the
   /// new index where column \p i has been moved to; note that \p perm must point to an array of size at least
   /// #numColsReal()
   void SoPlex::_removeColsReal(int perm[])
   {
      assert(_realLP != 0);

      _realLP->removeCols(perm);

      if( _isRealLPLoaded )
      {
         _hasBasis = (_solver.basis().status() > SPxBasis::NO_PROBLEM);
      }
      else if( _hasBasis )
      {
         for( int i = numColsReal() - 1; i >= 0 && _hasBasis; i-- )
         {
            if( perm[i] < 0 && _basisStatusCols[i] == SPxSolver::BASIC )
               _hasBasis = false;
            else if( perm[i] >= 0 && perm[i] != i )
            {
               assert(perm[i] < numColsReal());
               assert(perm[perm[i]] < 0);

               _basisStatusCols[perm[i]] = _basisStatusCols[i];
            }
         }

         if( _hasBasis )
            _basisStatusCols.reSize(numColsReal());
      }
   }



   /// invalidates solution
   void SoPlex::_invalidateSolution()
   {
      ///@todo maybe this should be done individually at the places when this method is called
      _status = SPxSolver::UNKNOWN;

      _solReal.invalidate();
      _hasSolReal = false;

      _solRational.invalidate();
      _hasSolRational = false;
   }



   /// enables simplifier and scaler
   void SoPlex::_enableSimplifierAndScaler()
   {
      // type of simplifier
      switch( intParam(SoPlex::SIMPLIFIER) )
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
      switch( intParam(SoPlex::SCALER) )
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
      default:
         break;
      }
   }



   /// disables simplifier and scaler
   void SoPlex::_disableSimplifierAndScaler()
   {
      _simplifier = 0;
      _scaler = 0;
   }



   /// ensures that the rational LP is available; performs no sync
   void SoPlex::_ensureRationalLP()
   {
      if( _rationalLP == 0 )
      {
         spx_alloc(_rationalLP);
         _rationalLP = new (_rationalLP) SPxLPRational();
         _rationalLP->setOutstream(spxout);
      }
   }



   /// ensures that the real LP and the basis are loaded in the solver; performs no sync
   void SoPlex::_ensureRealLPLoaded()
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
            assert(_basisStatusRows.size() == numRowsReal());
            assert(_basisStatusCols.size() == numColsReal());
            _solver.setBasis(_basisStatusRows.get_const_ptr(), _basisStatusCols.get_const_ptr());
            _hasBasis = (_solver.basis().status() > SPxBasis::NO_PROBLEM);
         }
      }
   }



   /// call floating-point solver and update statistics on iterations etc.
   void SoPlex::_solveRealLPAndRecordStatistics()
   {
      bool _hadBasis = _hasBasis;

      // set time and iteration limit
      if( intParam(SoPlex::ITERLIMIT) >= 0 )
         _solver.setTerminationIter(intParam(SoPlex::ITERLIMIT) - _statistics->iterations);
      if( realParam(SoPlex::TIMELIMIT) < realParam(SoPlex::INFTY) )
         _solver.setTerminationTime(realParam(SoPlex::TIMELIMIT) - _statistics->solvingTime->time());

      // ensure that tolerances are not too small
      if( _solver.feastol() < 1e-12 )
         _solver.setFeastol(1e-12);
      if( _solver.opttol() < 1e-12 )
         _solver.setOpttol(1e-12);

      // set correct representation
      if( (intParam(SoPlex::REPRESENTATION) == SoPlex::REPRESENTATION_COLUMN
            || (intParam(SoPlex::REPRESENTATION) == SoPlex::REPRESENTATION_AUTO && (_solver.nCols() + 1) * realParam(SoPlex::REPRESENTATION_SWITCH) >= (_solver.nRows() + 1)))
         && _solver.rep() != SPxSolver::COLUMN )
      {
         _solver.setRep(SPxSolver::COLUMN);
      }
      else if( (intParam(SoPlex::REPRESENTATION) == SoPlex::REPRESENTATION_ROW
            || (intParam(SoPlex::REPRESENTATION) == SoPlex::REPRESENTATION_AUTO && (_solver.nCols() + 1) * realParam(SoPlex::REPRESENTATION_SWITCH) < (_solver.nRows() + 1)))
         &&_solver.rep() != SPxSolver::ROW )
      {
         _solver.setRep(SPxSolver::ROW);
      }

      // set correct type
      if( ((intParam(ALGORITHM) == SoPlex::ALGORITHM_PRIMAL && _solver.rep() == SPxSolver::COLUMN)
            || (intParam(ALGORITHM) == SoPlex::ALGORITHM_DUAL && _solver.rep() == SPxSolver::ROW))
         && _solver.type() != SPxSolver::ENTER )
      {
         _solver.setType(SPxSolver::ENTER);
      }
      else if( ((intParam(ALGORITHM) == SoPlex::ALGORITHM_DUAL && _solver.rep() == SPxSolver::COLUMN)
            || (intParam(ALGORITHM) == SoPlex::ALGORITHM_PRIMAL && _solver.rep() == SPxSolver::ROW))
         && _solver.type() != SPxSolver::LEAVE )
      {
         _solver.setType(SPxSolver::LEAVE);
      }

      // set pricing modes
      _solver.setSparsePricingFactor(realParam(SoPlex::SPARSITY_THRESHOLD));
      if( (intParam(SoPlex::HYPER_PRICING) == SoPlex::HYPER_PRICING_ON)
            || ((intParam(SoPlex::HYPER_PRICING) == SoPlex::HYPER_PRICING_AUTO)
            && (_solver.nRows() + _solver.nCols() > HYPERPRICINGTHRESHOLD )) )
         _solver.hyperPricing(true);
      else if( intParam(SoPlex::HYPER_PRICING) == SoPlex::HYPER_PRICING_OFF )
         _solver.hyperPricing(false);

      // call floating-point solver and catch exceptions
      _statistics->simplexTime->start();
      try
      {
         _solver.solve();
      }
      catch( const SPxException& E )
      {
         MSG_ERROR( std::cerr << "Caught exception <" << E.what() << "> while solving real LP.\n" );
         _status = SPxSolver::ERROR;
      }
      catch( ... )
      {
         MSG_ERROR( std::cerr << "Caught unknown exception while solving real LP.\n" );
         _status = SPxSolver::ERROR;
      }
      _statistics->simplexTime->stop();

      // record statistics
      _statistics->iterations += _solver.iterations();
      _statistics->iterationsPrimal += _solver.primalIterations();
      _statistics->iterationsFromBasis += _hadBasis ? _solver.iterations() : 0;
      _statistics->boundflips += _solver.boundFlips();
      _statistics->luFactorizationTimeReal += _slufactor.getFactorTime();
      _statistics->luSolveTimeReal += _slufactor.getSolveTime();
      _statistics->luFactorizationsReal += _slufactor.getFactorCount();
      _statistics->luSolvesReal += _slufactor.getSolveCount();
      _slufactor.resetCounters();
   }



   /// reads real LP in LP or MPS format from file and returns true on success; gets row names, column names, and
   /// integer variables if desired
   bool SoPlex::_readFileReal(const char* filename, NameSet* rowNames, NameSet* colNames, DIdxSet* intVars)
   {
      assert(_realLP != 0);

      // clear statistics
      _statistics->clearAllData();

      // update status
      _status = SPxSolver::UNKNOWN;
      _invalidateSolution();
      _hasBasis = false;

      // start timing
      _statistics->readingTime->start();

      // read
      bool success = _realLP->readFile(filename, rowNames, colNames, intVars);

      // stop timing
      _statistics->readingTime->stop();

      if( success )
      {
         setIntParam(SoPlex::OBJSENSE, (_realLP->spxSense() == SPxLPReal::MAXIMIZE ? SoPlex::OBJSENSE_MAXIMIZE : SoPlex::OBJSENSE_MINIMIZE), true, true);
         _realLP->changeObjOffset(0.0);

         // if sync mode is auto, we have to copy the (rounded) real LP to the rational LP; this is counted to sync time
         // and not to reading time
         if( intParam(SoPlex::SYNCMODE) == SYNCMODE_AUTO )
            _syncLPRational();
      }
      else
         clearLPReal();

      return success;
   }



   /// reads rational LP in LP or MPS format from file and returns true on success; gets row names, column names, and
   /// integer variables if desired
   bool SoPlex::_readFileRational(const char* filename, NameSet* rowNames, NameSet* colNames, DIdxSet* intVars)
   {
      // clear statistics
      _statistics->clearAllData();

      // start timing
      _statistics->readingTime->start();

      // update status
      _status = SPxSolver::UNKNOWN;
      _invalidateSolution();
      _hasBasis = false;

      // read
      _ensureRationalLP();
      bool success = _rationalLP->readFile(filename, rowNames, colNames, intVars);

      // stop timing
      _statistics->readingTime->stop();

      if( success )
      {
         setIntParam(SoPlex::OBJSENSE, (_rationalLP->spxSense() == SPxLPRational::MAXIMIZE ? SoPlex::OBJSENSE_MAXIMIZE : SoPlex::OBJSENSE_MINIMIZE), true, true);
         _rationalLP->changeObjOffset(0);
         _recomputeRangeTypesRational();

         // if sync mode is auto, we have to copy the (rounded) real LP to the rational LP; this is counted to sync time
         // and not to reading time
         if( intParam(SoPlex::SYNCMODE) == SYNCMODE_AUTO )
            _syncLPReal();
         // if a rational LP file is read, but only the (rounded) real LP should be kept, we have to free the rational LP
         else if( intParam(SoPlex::SYNCMODE) == SYNCMODE_ONLYREAL )
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



   /// recomputes range types from scratch using real LP
   void SoPlex::_recomputeRangeTypesReal()
   {
      _rowTypes.reSize(numRowsReal());
      for( int i = 0; i < numRowsReal(); i++ )
         _rowTypes[i] = _rangeTypeReal(_realLP->lhs(i), _realLP->rhs(i));
      _colTypes.reSize(numColsReal());
      for( int i = 0; i < numColsReal(); i++ )
         _colTypes[i] = _rangeTypeReal(_realLP->lower(i), _realLP->upper(i));
   }



   /// recomputes range types from scratch using rational LP
   void SoPlex::_recomputeRangeTypesRational()
   {
      _rowTypes.reSize(numRowsRational());
      for( int i = 0; i < numRowsRational(); i++ )
         _rowTypes[i] = _rangeTypeRational(_rationalLP->lhs(i), _rationalLP->rhs(i));
      _colTypes.reSize(numColsRational());
      for( int i = 0; i < numColsRational(); i++ )
         _colTypes[i] = _rangeTypeRational(_rationalLP->lower(i), _rationalLP->upper(i));
   }



   /// synchronizes real LP with rational LP, i.e., copies (rounded) rational LP into real LP, without looking at the sync mode
   void SoPlex::_syncLPReal(bool time)
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

      // stop timing
      if( time )
         _statistics->syncTime->stop();
   }



   /// synchronizes rational LP with real LP, i.e., copies real LP to rational LP, without looking at the sync mode
   void SoPlex::_syncLPRational(bool time)
   {
      // start timing
      if( time )
         _statistics->syncTime->start();

      // copy LP
      _ensureRationalLP();
      *_rationalLP = *_realLP;
      _recomputeRangeTypesReal();

      // stop timing
      if( time )
         _statistics->syncTime->stop();
   }



   /// synchronizes real solution with rational solution, i.e., copies real solution to rational solution
   void SoPlex::_syncRealSolution()
   {
      if( _hasSolRational && !_hasSolReal )
      {
         _solReal = _solRational;
         _hasSolReal = true;
      }
   }



   /// synchronizes rational solution with real solution, i.e., copies (rounded) rational solution to real solution
   void SoPlex::_syncRationalSolution()
   {
      if( _hasSolReal && !_hasSolRational )
      {
         _solRational = _solReal;
         _hasSolRational = true;
      }
   }



   /// returns pointer to a constant unit vector available until destruction of the SoPlex class
   const UnitVectorRational* SoPlex::_unitVectorRational(const int i)
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
   bool SoPlex::_parseSettingsLine(char* line, const int lineNumber)
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
            MSG_ERROR( std::cerr << "Error parsing settings file: no ':' separating parameter type and name in line " << lineNumber << ".\n" );
            return false;
         }
         line++;
      }

      // find the start of the parameter name
      while( *line == ' ' || *line == '\t' || *line == '\r' )
         line++;
      if( *line == '\0' || *line == '\n' || *line == '#' )
      {
         MSG_ERROR( std::cerr << "Error parsing settings file: no parameter name in line " << lineNumber << ".\n");
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
            MSG_ERROR( std::cerr << "Error parsing settings file: no '=' after parameter name in line " << lineNumber << ".\n" );
            return false;
         }
         line++;
      }

      // find the start of the parameter value string
      while( *line == ' ' || *line == '\t' || *line == '\r' )
         line++;
      if( *line == '\0' || *line == '\n' || *line == '#' )
      {
         MSG_ERROR( std::cerr << "Error parsing settings file: no parameter value in line " << lineNumber << ".\n");
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
            MSG_ERROR( std::cerr << "Error parsing settings file: additional character '" << *line << "' after parameter value in line " << lineNumber << ".\n" );
            return false;
         }
      }

      // check whether we have a bool parameter
      if( strncmp(paramTypeString, "bool", 4) == 0 )
      {
         for( int param = 0; ; param++ )
         {
            if( param >= SoPlex::BOOLPARAM_COUNT )
            {
               MSG_ERROR( std::cerr << "Error parsing settings file: unknown parameter name <" << paramName << "> in line " << lineNumber << ".\n" );
               return false;
            }
            else if( strncmp(paramName, _currentSettings->_boolParamName[param].c_str(), SET_MAX_LINE_LEN) == 0 )
            {
               if( strncasecmp(paramValueString, "true", 4) == 0
                  || strncasecmp(paramValueString, "TRUE", 4) == 0
                  || strncasecmp(paramValueString, "t", 4) == 0
                  || strncasecmp(paramValueString, "T", 4) == 0
                  || strtol(paramValueString, NULL, 4) == 1 )
               {
                  setBoolParam((SoPlex::BoolParam)param, true);
                  break;
               }
               else if( strncasecmp(paramValueString, "false", 5) == 0
                  || strncasecmp(paramValueString, "FALSE", 5) == 0
                  || strncasecmp(paramValueString, "f", 5) == 0
                  || strncasecmp(paramValueString, "F", 5) == 0
                  || strtol(paramValueString, NULL, 5) == 0 )
               {
                  setBoolParam((SoPlex::BoolParam)param, false);
                  break;
               }
               else
               {
                  MSG_ERROR( std::cerr << "Error parsing settings file: invalid value <" << paramValueString << "> for bool parameter <" << paramName << "> in line " << lineNumber << ".\n" );
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
            if( param >= SoPlex::INTPARAM_COUNT )
            {
               MSG_ERROR( std::cerr << "Error parsing settings file: unknown parameter name <" << paramName << "> in line " << lineNumber << ".\n" );
               return false;
            }
            else if( strncmp(paramName, _currentSettings->_intParamName[param].c_str(), SET_MAX_LINE_LEN) == 0 )
            {
               int value;

               if( sscanf(paramValueString, "%d", &value) == 1 && setIntParam((SoPlex::IntParam)param, value) )
                  break;
               else
               {
                  MSG_ERROR( std::cerr << "Error parsing settings file: invalid value <" << paramValueString << "> for int parameter <" << paramName << "> in line " << lineNumber << ".\n" );
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
            if( param >= SoPlex::REALPARAM_COUNT )
            {
               MSG_ERROR( std::cerr << "Error parsing settings file: unknown parameter name <" << paramName << "> in line " << lineNumber << ".\n" );
               return false;
            }
            else if( strncmp(paramName, _currentSettings->_realParamName[param].c_str(), SET_MAX_LINE_LEN) == 0 )
            {
               Real value;

               if( sscanf(paramValueString, "%" REAL_FORMAT, &value) == 1 && setRealParam((SoPlex::RealParam)param, value) )
                  break;
               else
               {
                  MSG_ERROR( std::cerr << "Error parsing settings file: invalid value <" << paramValueString << "> for real parameter <" << paramName << "> in line " << lineNumber << ".\n" );
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
            if( param >= SoPlex::RATIONALPARAM_COUNT )
            {
               MSG_ERROR( std::cerr << "Error parsing settings file: unknown parameter name <" << paramName << "> in line " << lineNumber << ".\n" );
               return false;
            }
            else if( strncmp(paramName, _currentSettings->_rationalParamName[param].c_str(), SET_MAX_LINE_LEN) == 0 )
            {
               Rational value;

               if( readStringRational(paramValueString, value) && setRationalParam((SoPlex::RationalParam)param, value) )
                  break;
               else
               {
                  MSG_ERROR( std::cerr << "Error parsing settings file: invalid value <" << paramValueString << "> for rational parameter <" << paramName << "> in line " << lineNumber << ".\n" );
                  return false;
               }
            }
         }

         return true;
      }
#endif

      MSG_ERROR( std::cerr << "Error parsing settings file: invalid parameter type <" << paramTypeString << "> for parameter <" << paramName << "> in line " << lineNumber << ".\n" );

      return false;
   }
} // namespace soplex
#endif
