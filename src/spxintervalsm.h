/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the class library                   */
/*       SoPlex --- the Sequential object-oriented simPlex.                  */
/*                                                                           */
/*    Copyright (C) 1997-1999 Roland Wunderling                              */
/*                  1997-2002 Konrad-Zuse-Zentrum                            */
/*                            fuer Informationstechnik Berlin                */
/*                                                                           */
/*  SoPlex is distributed under the terms of the ZIB Academic Licence.       */
/*                                                                           */
/*  You should have received a copy of the ZIB Academic License              */
/*  along with SoPlex; see the file COPYING. If not email to soplex@zib.de.  */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
#pragma ident "@(#) $Id: spxintervalsm.h,v 1.2 2003/01/10 12:46:14 bzfkocht Exp $"

/**@file  spxintervalsm.h
 * @brief Remove all values outside the interval [#epsilon()..#infinity] 
 *        from the LP.
 */
#ifndef _SPXINTERVALSM_H_
#define _SPXINTERVALSM_H_

#include "spxdefines.h"
#include "spxsimplifier.h"

namespace soplex
{
/**@brief   Simplifier for forcing the LP values into a [#epsilon()..#infinity] interval.
 * @ingroup Algo
 *
 * This simplifier sets all values smaller than #epsilon() to zero and all values
 * bigger then infinity/5 to infinity.
 * 
 */
class SPxIntervalSM : public SPxSimplifier
{
public:
   /// default constructor
   SPxIntervalSM() 
      : SPxSimplifier("Interval")
   {}   
   /// just do it!
   Result simplify(SPxLP& lp, Real eps, Real delta);
};
} // namespace soplex
#endif // _SPXINTERVALSM_H_

//-----------------------------------------------------------------------------
//Emacs Local Variables:
//Emacs mode:c++
//Emacs c-basic-offset:3
//Emacs tab-width:8
//Emacs indent-tabs-mode:nil
//Emacs End:
//-----------------------------------------------------------------------------
