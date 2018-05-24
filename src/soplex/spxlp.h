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

/**@file  spxlp.h
 * @brief Saving LPs in a form suitable for SoPlex.
 */

#ifndef _SPXLP_H_
#define _SPXLP_H_

#include "soplex/spxdefines.h"
#include "soplex/spxlpbase.h"
#include "soplex/vector.h" // for compatibility
#include "soplex/svector.h" // for compatibility
#include "soplex/svset.h" // for compatibility
#include "soplex/lprowset.h" // for compatibility
#include "soplex/lpcolset.h" // for compatibility
#include "soplex/lprow.h" // for compatibility
#include "soplex/lpcol.h" // for compatibility

namespace soplex
{
typedef SPxLPBase< Real > SPxLP;
typedef SPxLPBase< Real > SPxLPReal;
typedef SPxLPBase< Rational > SPxLPRational;
} // namespace soplex
#endif // _SPXLP_H_
