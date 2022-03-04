/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the class library                   */
/*       SoPlex --- the Sequential object-oriented simPlex.                  */
/*                                                                           */
/*    Copyright (C) 1996-2022 Konrad-Zuse-Zentrum                            */
/*                            fuer Informationstechnik Berlin                */
/*                                                                           */
/*  SoPlex is distributed under the terms of the ZIB Academic Licence.       */
/*                                                                           */
/*  You should have received a copy of the ZIB Academic License              */
/*  along with SoPlex; see the file COPYING. If not email to soplex@zib.de.  */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/**@file  validation.hpp
 * @brief Validation object for soplex solutions
 */

#include "soplex/validation.h"

namespace soplex
{

template <class R>
bool Validation<R>::updateExternalSolution(const std::string& solution)
{
   validate = true;
   validatesolution = solution;

   if(solution == "+infinity")
   {
      return true;
   }
   else if(solution == "-infinity")
   {
      return true;
   }
   else
   {
      char* tailptr;
      strtod(solution.c_str(), &tailptr);

      if(*tailptr)
      {
         //conversion failed because the input wasn't a number
         return false;
      }
   }

   return true;
}


/// updates the tolerance used for validation
template <class R>
bool Validation<R>::updateValidationTolerance(const std::string& tolerance)
{
   char* tailptr;
   validatetolerance = strtod(tolerance.c_str(), &tailptr);

   if(*tailptr)
   {
      //conversion failed because the input wasn't a number
      return false;
   }

   return true;
}

template <class R>
void Validation<R>::validateSolveReal(SoPlexBase<R>& soplex)
{
   bool passedValidation = true;
   std::string reason = "";
   R objViolation = 0.0;
   R maxBoundViolation = 0.0;
   R maxRowViolation = 0.0;
   R maxRedCostViolation = 0.0;
   R maxDualViolation = 0.0;
   R sumBoundViolation = 0.0;
   R sumRowViolation = 0.0;
   R sumRedCostViolation = 0.0;
   R sumDualViolation = 0.0;
   R sol;

   std::ostream& os = soplex.spxout.getStream(SPxOut::INFO1);

   if(validatesolution == "+infinity")
   {
      sol = soplex.realParam(SoPlexBase<R>::INFTY);
   }
   else if(validatesolution == "-infinity")
   {
      sol = -soplex.realParam(SoPlexBase<R>::INFTY);
   }
   else
   {
      sol = atof(validatesolution.c_str());
   }

   objViolation = spxAbs(sol - soplex.objValueReal());

   // skip check in case presolving detected infeasibility/unboundedness
   if(SPxSolverBase<R>::INForUNBD == soplex.status() &&
         (sol == soplex.realParam(SoPlexBase<R>::INFTY)
          || sol == -soplex.realParam(SoPlexBase<R>::INFTY)))
      objViolation = 0.0;

   if(! EQ(objViolation, R(0.0), validatetolerance))
   {
      passedValidation = false;
      reason += "Objective Violation; ";
   }

   if(SPxSolverBase<R>::OPTIMAL == soplex.status())
   {
      soplex.getBoundViolation(maxBoundViolation, sumBoundViolation);
      soplex.getRowViolation(maxRowViolation, sumRowViolation);
      soplex.getRedCostViolation(maxRedCostViolation, sumRedCostViolation);
      soplex.getDualViolation(maxDualViolation, sumDualViolation);

      if(! LE(maxBoundViolation, validatetolerance))
      {
         passedValidation = false;
         reason += "Bound Violation; ";
      }

      if(! LE(maxRowViolation, validatetolerance))
      {
         passedValidation = false;
         reason += "Row Violation; ";
      }

      if(! LE(maxRedCostViolation, validatetolerance))
      {
         passedValidation = false;
         reason += "Reduced Cost Violation; ";
      }

      if(! LE(maxDualViolation, validatetolerance))
      {
         passedValidation = false;
         reason += "Dual Violation; ";
      }
   }

   os << "\n";
   os << "Validation          :";

   if(passedValidation)
      os << " Success\n";
   else
   {
      reason[reason.length() - 2] = ']';
      os << " Fail [" + reason + "\n";
   }

   os << "   Objective        : " << std::scientific << std::setprecision(
         8) << objViolation << std::fixed << "\n";
   os << "   Bound            : " << std::scientific << std::setprecision(
         8) << maxBoundViolation << std::fixed << "\n";
   os << "   Row              : " << std::scientific << std::setprecision(
         8) << maxRowViolation << std::fixed << "\n";
   os << "   Reduced Cost     : " << std::scientific << std::setprecision(
         8) << maxRedCostViolation << std::fixed << "\n";
   os << "   Dual             : " << std::scientific << std::setprecision(
         8) << maxDualViolation << std::fixed << "\n";
}

} // namespace soplex
