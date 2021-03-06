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

/**@file  ratrecon.h
 * @brief Rational reconstruction of solution vector
 */


#ifndef _RATRECON_H_
#define _RATRECON_H_

#include "soplex/spxdefines.h"
#include "soplex/rational.h"
#include "soplex/sol.h"
#include "soplex/basevectors.h"
#include "soplex/didxset.h"

namespace soplex
{
/** reconstruct a rational vector */
bool reconstructVector(VectorRational& input, const Rational& denomBoundSquared,
                       const DIdxSet* indexSet = 0);

/** reconstruct a rational solution */
bool reconstructSol(SolRational& solution);
} // namespace soplex

#include "soplex/ratrecon.hpp"

#endif // _RATRECON_H_
