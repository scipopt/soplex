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

#include <iostream>
#include <assert.h>

#include "statistics.h"

namespace soplex
{
   /// default constructor
   SoPlex2::Statistics::Statistics()
   {
      clearAllData();
   }

   /// clears all statistics
   void SoPlex2::Statistics::clearAllData()
   {
      readingTime.reset();
      clearSolvingData();
   }

   /// clears statistics on solving process
   void SoPlex2::Statistics::clearSolvingData()
   {
      solvingTime.reset();
      preprocessingTime.reset();
      simplexTime.reset();
      syncTime.reset();
      transformTime.reset();
      luFactorizationTime = 0.0;
      luSolveTime = 0.0;
      iterations = 0;
      luFactorizations = 0;
      luSolves = 0;
      refinements = 0;
      stallRefinements = 0;
   }

   /// prints statistics
   void SoPlex2::Statistics::print(std::ostream& os)
   {
      Real solTime = solvingTime.userTime();
      Real totTime = readingTime.userTime() + solTime;
      Real otherTime = solTime - syncTime.userTime() - transformTime.userTime() - preprocessingTime.userTime() - simplexTime.userTime();

      os << std::fixed << std::setprecision(2);

      os << "Total time         : " << totTime << " seconds\n"
         << "  Reading          : " << readingTime.userTime() << "\n"
         << "  Solving          : " << solTime << "\n"
         << "  Preprocessing    : " << preprocessingTime.userTime() << " (" << 100 * (preprocessingTime.userTime() / solTime) << "% of solving time)\n"
         << "  Simplex          : " << simplexTime.userTime() << " (" << 100 * (simplexTime.userTime() / solTime) << "% of solving time)\n"
         << "  Synchronization  : " << syncTime.userTime() << " (" << 100 * (syncTime.userTime() / solTime) << "% of solving time)\n"
         << "  Transformation   : " << transformTime.userTime() << " (" << 100*transformTime.userTime() / solTime << "% of solving time)\n"
         << "  Other            : " << otherTime << " (" << 100*otherTime / solTime << "% of solving time)\n";

      os << "Refinements        : " << refinements << "\n"
         << "  Stalling         : " << stallRefinements << "\n";

      os << "Iterations         : " << iterations << "\n"
         << "  From scratch     : " << iterations - iterationsFromBasis << " (" << 100*double((iterations - iterationsFromBasis))/double(iterations) << "%)\n"
         << "  From basis       : " << iterationsFromBasis << " (" << 100*double(iterationsFromBasis)/double(iterations) << "%)\n"
         << "  Primal           : " << iterationsPrimal << " (" << 100*double(iterationsPrimal)/double(iterations) << "%)\n"
         << "  Dual             : " << iterations - iterationsPrimal << " (" << 100*double((iterations - iterationsPrimal))/double(iterations) << "%)\n";

      os << "LU factorizations  : " << luFactorizations << "\n"
         << "  Factor. frequency: ";
      if( luFactorizations > 0 )
         os << double(iterations) / double(luFactorizations) << " iterations per factorization\n";
      else
         os << "-\n";
      os << "  Factor. time     : " << luFactorizationTime << " seconds\n";

      os << "LU solves          : " << luSolves << "\n"
         << "  Solve frequency  : ";
      if( luSolves > 0 )
         os << double(luSolves) / double(iterations) << " solves per iteration\n";
      else
         os << "-\n";
      os << "  Solve time       : " << luSolveTime << " seconds\n";
   }
} // namespace soplex

//-----------------------------------------------------------------------------
//Emacs Local Variables:
//Emacs mode:c++
//Emacs c-basic-offset:3
//Emacs tab-width:8
//Emacs indent-tabs-mode:nil
//Emacs End:
//-----------------------------------------------------------------------------
