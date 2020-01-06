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

/**@file  spxscaler.hpp
 * @brief LP scaling base class.
 */

#include <cmath>

#include <iostream>
#include <assert.h>
#include "soplex/spxscaler.h"
#include "soplex/spxlp.h"
#include "soplex/dsvector.h"
#include <limits>

namespace soplex
{


template <>
int SPxScaler<Rational>::computeScaleExp(const SVectorBase<Rational>& vec,
      const DataArray<int>& oldScaleExp) const
{
   // This should never be called.
   // Refer to issue 164 in soplex gitlab

   assert(false);

   return 0;
}

} // namespace soplex
