/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the class library                   */
/*       SoPlex --- the Sequential object-oriented simPlex.                  */
/*                                                                           */
/*    Copyright (C) 1996-2019 Konrad-Zuse-Zentrum                            */
/*                            fuer Informationstechnik Berlin                */
/*                                                                           */
/*  SoPlex is distributed under the terms of the ZIB Academic Licence.       */
/*                                                                           */
/*  You should have received a copy of the ZIB Academic License              */
/*  along with SoPlex; see the file COPYING. If not email to soplex@zib.de.  */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/**@file  validation.cpp
 * @brief Validation object for soplex solutions
 */

#include "soplex/validation.h"

#ifdef SOPLEX_WITH_BOOST
#include "boost/lexical_cast.hpp"
#endif

namespace soplex
{

#ifdef SOPLEX_WITH_BOOST
/// updates the external solution used for validation
template <>
bool Validation<Real>::updateExternalSolution(const std::string& solution)
{
   validate = true;
   validatesolution = solution;

   if(solution == "+infinity")
   {
      return true;
   }
   else if(solution == "-infinity")
   {
      return false;
   }
   else
   {
      // This will throw boost::bad_lexical cast if bad cast. Will be caught
      // by the catch in soplexmain.cpp
      boost::lexical_cast<double>(solution);
   }

   return true;
}



/// updates the tolerance used for validation
template <>
bool Validation<Real>::updateValidationTolerance(const std::string& tolerance)
{
   // Will throw boost::bad_lexical_cast if conversion fails
   validatetolerance = boost::lexical_cast<double>(tolerance);
   return true;
}



/// validates the soplex solution using the external solution
template <>
void Validation<Real>::validateSolveReal(SoPlexBase<Real>& soplex)
{
   bool passedValidation = true;
   std::string reason = "";
   Real objViolation = 0.0;
   Real maxBoundViolation = 0.0;
   Real maxRowViolation = 0.0;
   Real maxRedCostViolation = 0.0;
   Real maxDualViolation = 0.0;
   Real sumBoundViolation = 0.0;
   Real sumRowViolation = 0.0;
   Real sumRedCostViolation = 0.0;
   Real sumDualViolation = 0.0;
   Real sol;

   std::ostream& os = soplex.spxout.getStream(SPxOut::INFO1);

   if(validatesolution == "+infinity")
   {
      sol = soplex.realParam(SoPlexBase<Real>::INFTY);
   }
   else if(validatesolution == "-infinity")
   {
      sol = -soplex.realParam(SoPlexBase<Real>::INFTY);
   }
   else
   {
      // This will not throw here because it was checked in updateExternalSolution()
      sol = boost::lexical_cast<double>(validatesolution);
   }

   objViolation = spxAbs(sol - soplex.objValueReal());

   // skip check in case presolving detected infeasibility/unboundedness
   if(SPxSolverBase<Real>::INForUNBD == soplex.status() &&
         (sol == soplex.realParam(SoPlexBase<Real>::INFTY)
          || sol == -soplex.realParam(SoPlexBase<Real>::INFTY)))
      objViolation = 0.0;

   if(! EQ(objViolation, 0.0, validatetolerance))
   {
      passedValidation = false;
      reason += "Objective Violation; ";
   }

   if(SPxSolverBase<Real>::OPTIMAL == soplex.status())
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
#endif

} /* namespace soplex */
