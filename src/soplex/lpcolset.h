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

/**@file  lpcolset.h
 * @brief Set of LP columns.
 */
#ifndef _LPCOLSET_H_
#define _LPCOLSET_H_

#include "soplex/spxdefines.h"
#include "soplex/lpcolsetbase.h"

namespace soplex
{
typedef LPColSetBase< Real > LPColSet;
typedef LPColSetBase< Real > LPColSetReal;
typedef LPColSetBase< Rational > LPColSetRational;
} // namespace soplex
#endif // _LPCOLSET_H_
