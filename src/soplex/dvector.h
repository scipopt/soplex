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

/**@file  dvector.h
 * @brief Dynamic vectors.
 */
#ifndef _DVECTOR_H_
#define _DVECTOR_H_

#include "soplex/spxdefines.h"
#include "soplex/basevectors.h"
#include "soplex/vector.h" // for compatibility

// This file exists for reverse compatibility with SCIP. This isn't currently
// needed in SoPlex. DVector used to be typedefs from DVectorBase<T>, but
// DVectorBase has been replaced by VectorBase.
namespace soplex
{
typedef VectorBase< Real > DVector;
typedef VectorBase< Real > DVectorReal;
typedef VectorBase< Rational > DVectorRational;
} // namespace soplex
#endif // _DVECTOR_H_
