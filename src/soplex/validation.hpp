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

/**@file  validation.hpp
 * @brief General templated functions for SoPlex
 */

/// updates the external solution used for validation
template <class R>
bool Validation<R>::updateExternalSolution(char* solution)
{
  validate = true;
  validatesolution = solution;

  if( strncmp(solution, "+infinity", 9 ) == 0 )
    return true;
  else if ( strncmp(solution, "-infinity", 9) == 0 )
    return true;
  else
    {
      char* tailptr;
      strtod(solution, &tailptr);
      if (*tailptr) {
        //conversion failed because the input wasn't a number
        return false;
      }
    }
  return true;
}



/// updates the tolerance used for validation
template <class R>
bool Validation<R>::updateValidationTolerance(char* tolerance)
{
  char* tailptr;
  validatetolerance = strtod(tolerance, &tailptr);
  if (*tailptr) {
    //conversion failed because the input wasn't a number
    return false;
  }
  return true;
}

