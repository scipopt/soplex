/*
 * validation.cpp
 *
 *  Created on: 15.08.2017
 *      Author: bzfviern
 */

#include "validation.h"

namespace soplex {



bool Validation::updateExternalSolution(char* solution)
{
   validatesolution = solution;
   validate = true;
   return true;
}

bool Validation::updateValidationTolerance(Real tolerance)
{
   validatetolerance = tolerance;
   validate = true;
   return true;
}

} /* namespace soplex */
