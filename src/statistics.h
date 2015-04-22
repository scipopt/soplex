/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the class library                   */
/*       SoPlex --- the Sequential object-oriented simPlex.                  */
/*                                                                           */
/*    Copyright (C) 1996-2014 Konrad-Zuse-Zentrum                            */
/*                            fuer Informationstechnik Berlin                */
/*                                                                           */
/*  SoPlex is distributed under the terms of the ZIB Academic Licence.       */
/*                                                                           */
/*  You should have received a copy of the ZIB Academic License              */
/*  along with SoPlex; see the file COPYING. If not email to soplex@zib.de.  */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/**@file  statistics.h
 * @brief Class for collecting statistical information
 */
#ifndef _STATISTICS_H_
#define _STATISTICS_H_

#ifndef SOPLEX_LEGACY
#include <iostream>

#include "soplex.h"
#include "timer.h"

namespace soplex
{
   /**@class   Statistics
    * @brief   Class for collecting statistical information
    * @ingroup Algo
    */
   class SoPlex::Statistics
   {
   public:

      //**@name Construction, resetting, printing */
      //@{

      /// default constructor
      Statistics();

      /// clears all statistics
      void clearAllData();

      /// clears statistics on solving process
      void clearSolvingData();

      /// prints statistics
      void print(std::ostream& os);

      //@}


      //**@name Data */
      //@{

      Timer readingTime; ///< reading time not included in solving time
      Timer solvingTime; ///< solving time
      Timer preprocessingTime; ///< preprocessing time
      Timer simplexTime; ///< simplex time
      Timer syncTime; ///< time for synchronization between real and rational LP (included in solving time)
      Timer transformTime; ///< time for transforming LPs (included in solving time)
      Real luFactorizationTime; ///< time for factorizing bases matrices
      Real luSolveTime; ///< time for solving linear systems
      int iterations; ///< number of iterations/pivots
      int iterationsPrimal; ///< number of iterations with Primal
      int iterationsFromBasis; ///< number of iterations from Basis
      int boundflips; ///< number of dual bound flips
      int luFactorizations; ///< number of basis matrix factorizations
      int luSolves; ///< number of (forward and backward) solves with basis matrix
      int refinements; ///< number of refinement steps
      int stallRefinements; ///< number of refinement steps without subsequent pivots

      // Improved dual simplex statistics
      int callsReducedProb;      ///< number of times the reduced problem is solved. This includes the initial solve.
      int iterationsInit;        ///< number of iterations in the initial LP
      int iterationsRedProb;     ///< number of iterations of the reduced problem
      int iterationsCompProb;    ///< number of iterations of the complementary problem
      int degenPivotsPrimal;     ///< number of primal degenerate pivots
      int degenPivotsDual;       ///< number of dual degenerate pivots
      int degenPivotCandPrimal;  ///< number of pivoting candidates that will produce a degenerate step in the primal
      int degenPivotCandDual;    ///< number of pivoting candidates that will produce a degenerate step in the dual

      //@}
   };
} // namespace soplex
#endif
#endif // _STATISTICS_H_
