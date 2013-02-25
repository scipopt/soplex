/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the class library                   */
/*       SoPlex --- the Sequential object-oriented simPlex.                  */
/*                                                                           */
/*    Copyright (C) 1996-2013 Konrad-Zuse-Zentrum                            */
/*                            fuer Informationstechnik Berlin                */
/*                                                                           */
/*  SoPlex is distributed under the terms of the ZIB Academic Licence.       */
/*                                                                           */
/*  You should have received a copy of the ZIB Academic License              */
/*  along with SoPlex; see the file COPYING. If not email to soplex@zib.de.  */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/**@file  rational.cpp
 * @brief Wrapper for GMP types.
 */

#include <math.h>

#include "rational.h"

namespace soplex
{
#ifdef SOPLEX_WITH_GMP

/// Negation.
Rational operator-(const Rational& q)
{
   Rational res = q;
   res *= -1;
   return res;
}

/// Division.
Rational operator/(const Rational& p, const Rational& q)
{
   Rational res = p;
   res /= q;
   return res;
}

/// print Rational with limited floating point precision
std::ostream& operator<<(std::ostream& os, const Rational& q)
{
   os << mpf_class(q);
   return os;
}

#endif
} // namespace soplex

//-----------------------------------------------------------------------------
//Emacs Local Variables:
//Emacs mode:c++
//Emacs c-basic-offset:3
//Emacs tab-width:8
//Emacs indent-tabs-mode:nil
//Emacs End:
//-----------------------------------------------------------------------------
