/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the class library                   */
/*       SoPlex --- the Sequential object-oriented simPlex.                  */
/*                                                                           */
/*    Copyright (C) 1997-1999 Roland Wunderling                              */
/*                  1997-2007 Konrad-Zuse-Zentrum                            */
/*                            fuer Informationstechnik Berlin                */
/*                                                                           */
/*  SoPlex is distributed under the terms of the ZIB Academic Licence.       */
/*                                                                           */
/*  You should have received a copy of the ZIB Academic License              */
/*  along with SoPlex; see the file COPYING. If not email to soplex@zib.de.  */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
#pragma ident "@(#) $Id: spxequilisc.h,v 1.8 2008/09/22 10:02:45 bzftuchs Exp $"

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

   This SPxScaler implementation performs equilibrium scaling of the 
   LPs rows and columns.
*/
class SPxEquiliSC : public SPxScaler
{
public:

   //-------------------------------------
   /**@name Construction / destruction */
   //@{
   /// default constructor (this scaler makes no use of inherited member m_colFirst)
   explicit SPxEquiliSC(bool doBoth = true);
   /// destructor
   virtual ~SPxEquiliSC()
   {}
   //@}

   //-------------------------------------
   /**@name Scaling */
   //@{
   /// Scale the loaded SPxLP.
   virtual void scale(SPxLP& lp);
   //@}

protected:

   //-------------------------------------
   /**@name Protected helpers */
   //@{
   /// Does nothing but returning \p maxi.
   virtual Real computeScale(Real /*mini*/, Real maxi) const;
   //@}

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
