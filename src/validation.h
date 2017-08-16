/*
 * validation.h
 *
 *  Created on: 15.08.2017
 *      Author: bzfviern
 */

#ifndef SRC_VALIDATION_H_
#define SRC_VALIDATION_H_

#include "soplex.h"

namespace soplex {

class Validation
{
public:
   bool           validate;
   double         validatetolerance;
   char*          validatesolution;
   Validation()
   {
      validate = false;
      validatetolerance = 1e-5;
      validatesolution = 0;
   }
   ~Validation()
   {
      ;
   }
   bool updateExternalSolution(char* solution);
   bool updateValidationTolerance(Real tolerance);
   bool validateSolveReal(SoPlex* soplex);
};

} /* namespace soplex */

#endif /* SRC_VALIDATION_H_ */
