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
#pragma ident "@(#) $Id: spxequilisc.h,v 1.4 2003/01/12 13:09:40 bzfkocht Exp $"

/**@file  spxequilisc.h
 * @brief LP euilibrium scaling.
 */
#ifndef _SPXEQUILISC_H_
#define _SPXEQUILISC_H_

#include <assert.h>

#include "spxdefines.h"
#include "spxscaler.h"

namespace soplex
{
/**@brief Equilibrium row/column scaling.
   @ingroup Algo

   This #SPxScaler implementation performs equilibrium scaling of the 
   LPs rows and columns.
*/
class SPxEquiliSC : public SPxScaler
{
protected:
   ///@return maxi
   virtual Real computeScale(Real /*mini*/, Real maxi) const;

public:
   /// Scale the loaded #SPxLP.
   virtual void scale(SPxLP& lp);

   /// default constructor.
   explicit SPxEquiliSC(bool colFirst = true, bool doBoth = true);
};
} // namespace soplex
#endif // _SPXEQUILISC_H_

//-----------------------------------------------------------------------------
//Emacs Local Variables:
//Emacs mode:c++
//Emacs c-basic-offset:3
//Emacs tab-width:8
//Emacs indent-tabs-mode:nil
//Emacs End:
//-----------------------------------------------------------------------------
