/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the class library                   */
/*       SoPlex --- the Sequential object-oriented simPlex.                  */
/*                                                                           */
/*    Copyright (C) 1997-1999 Roland Wunderling                              */
/*                  1997-2001 Konrad-Zuse-Zentrum                            */
/*                            fuer Informationstechnik Berlin                */
/*                                                                           */
/*  SoPlex is distributed under the terms of the ZIB Academic Licence.       */
/*                                                                           */
/*  You should have received a copy of the ZIB Academic License              */
/*  along with SoPlex; see the file COPYING. If not email to soplex@zib.de.  */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
#pragma ident "@(#) $Id: lprow.cpp,v 1.3 2001/11/11 20:27:31 bzfkocht Exp $"

/*      \Section{Complex Methods}
 */

/*  Import system include files
 */
#include <math.h>
#include <assert.h>
#include <iostream>


/*  and class header files
 */
#include "lprow.h"

namespace soplex
{


double LPRow::infinity = 1e+100;

//@ ----------------------------------------------------------------------------
LPRow::Type LPRow::type() const
{
   if (rhs() >= infinity)
      return GREATER_EQUAL;
   if (lhs() <= -infinity)
      return LESS_EQUAL;
   if (lhs() == rhs())
      return EQUAL;
   return RANGE;
}

// TK13OCT1998
void LPRow::setType(
   LPRow::Type type)
{
   switch (type)
   {
   case LESS_EQUAL:
      lhs() = -infinity;
      break;
   case EQUAL:
      if (lhs() > -infinity)
         rhs() = lhs();
      else
         lhs() = rhs();
      break;
   case GREATER_EQUAL:
      rhs() = infinity;
      break;
   case RANGE :
      std::cerr << __FILE__ << __LINE__
      << "RANGE not supported in LPRow::setType()";
      abort();
   default:
      abort();
   }
}

double LPRow::value() const
{
   assert(type() != RANGE);

   return (rhs() < infinity) ? rhs() : lhs();
}

LPRow::LPRow(
   const SVector& rowVector,
   LPRow::Type type,
   double value)
   : vec(rowVector)
{
   switch (type)
   {
   case LESS_EQUAL:
      left = -infinity;
      right = value;
      break;
   case EQUAL:
      left = value;
      right = value;
      break;
   case GREATER_EQUAL:
      left = value;
      right = infinity;
      break;
   default:
      abort();
   }
}
} // namespace soplex

//-----------------------------------------------------------------------------
//Emacs Local Variables:
//Emacs mode:c++
//Emacs c-basic-offset:3
//Emacs tab-width:8
//Emacs indent-tabs-mode:nil
//Emacs End:
//-----------------------------------------------------------------------------
