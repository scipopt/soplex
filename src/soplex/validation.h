/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the class library                   */
/*       SoPlex --- the Sequential object-oriented simPlex.                  */
/*                                                                           */
/*    Copyright (C) 1996-2020 Konrad-Zuse-Zentrum                            */
/*                            fuer Informationstechnik Berlin                */
/*                                                                           */
/*  SoPlex is distributed under the terms of the ZIB Academic Licence.       */
/*                                                                           */
/*  You should have received a copy of the ZIB Academic License              */
/*  along with SoPlex; see the file COPYING. If not email to soplex@zib.de.  */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/**@file  validation.h
 * @brief Validation object for soplex solutions
 */

#ifndef SRC_VALIDATION_H_
#define SRC_VALIDATION_H_

#include "soplex.h"

namespace soplex
{

template <class R>
class Validation
{
public:

   /// should the soplex solution be validated?
   bool           validate;

   /// external solution used for validation
   std::string          validatesolution;

   /// tolerance used for validation
   R         validatetolerance;

   /// default constructor
   Validation()
   {
      validate = false;
      validatetolerance = 1e-5;
   }

   /// default destructor
   ~Validation()
   {
      ;
   }

   /// updates the external solution used for validation
   bool updateExternalSolution(const std::string& solution);

   /// updates the tolerance used for validation
   bool updateValidationTolerance(const std::string& tolerance);

   /// validates the soplex solution using the external solution
   void validateSolveReal(SoPlexBase<R>& soplex);
};

} /* namespace soplex */

// For general templated functions
#include "validation.hpp"

#endif /* SRC_VALIDATION_H_ */
