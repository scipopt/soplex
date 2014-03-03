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

/**@file  ratrecon.h
 * @brief Rational reconstruction of solution vector
 */


#ifndef _RATRECON_H_
#define _RATRECON_H_

#include <math.h>
#include <assert.h>
#include <iostream>
#include <string.h>

#include "rational.h"
#include "sol.h"
#include "spxdefines.h"
#include "basevectors.h"
#include "gmp.h"
#include "gmpxx.h"

namespace soplex
{
   bool reconstructSol(SolRational& solution);
   bool reconstructVector(VectorBase<Rational>& input);

} // namespace soplex
#endif // _RATRECON_H_
