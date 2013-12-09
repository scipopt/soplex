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
      Real otherTime = solTime - syncTime.userTime() - transformTime.userTime();

      os << std::fixed << std::setprecision(2);

      os << "Total time         : " << totTime << " seconds\n"
         << "  Reading          : " << readingTime.userTime() << "\n"
         << "  Solving          : " << solTime << "\n"
         << "  Preprocessing    : " << "?" << " (?% of solving time)\n"
         << "  Simplex          : " << "?" << " (?% of solving time)\n"
         << "  Synchronization  : " << syncTime.userTime() << " (" << syncTime.userTime() / solTime << "% of solving time)\n"
         << "  Transformation   : " << transformTime.userTime() << " (" << transformTime.userTime() / solTime << "% of solving time)\n"
         << "  Other            : " << otherTime << " (" << otherTime / solTime << "% of solving time)\n";

      os << "Refinements        : " << refinements << "\n"
         << "  Stalling         : " << stallRefinements << "\n";

      os << "Iterations         : " << iterations << "\n"
         << "  From scratch     : " << "?" << " (?%)\n"
         << "  From basis       : " << "?" << " (?%)\n"
         << "  Primal           : " << "?" << " (?%)\n"
         << "  Dual             : " << "?" << " (?%)\n"
         << "  Column rep.      : " << "?" << " (?%)\n"
         << "  Row rep.         : " << "?" << " (?%)\n";

      os << "LU                 : " << "\n"
         << "  Factorizations   : " << luFactorizations << "\n"
         << "  Factor. frequency: ";
      if( luFactorizations > 0 )
         os << double(iterations) / double(luFactorizations) << " iterations per factorization\n";
      else
         os << "-\n";
      os << "  Factor. time     : " << luFactorizationTime << " seconds\n"
         << "  Solves           : " << luSolves << "\n"
         << "  Solve time       : " << luSolveTime << " seconds\n";
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
