/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the class library                   */
/*       SoPlex --- the Sequential object-oriented simPlex.                  */
/*                                                                           */
/*    Copyright (C) 1997-1999 Roland Wunderling                              */
/*                  1997-2008 Konrad-Zuse-Zentrum                            */
/*                            fuer Informationstechnik Berlin                */
/*                                                                           */
/*  SoPlex is distributed under the terms of the ZIB Academic Licence.       */
/*                                                                           */
/*  You should have received a copy of the ZIB Academic License              */
/*  along with SoPlex; see the file COPYING. If not email to soplex@zib.de.  */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
#pragma ident "@(#) $Id: spxdantzigpr.h,v 1.4 2008/09/22 20:43:18 bzfpfets Exp $"

/**@file  spxdantzigpr.h
 * @brief Dantzig pricer.
 */
#ifndef _SPXDEFAULTPR_H_
#define _SPXDEFAULTPR_H_

#include <assert.h>

#include "spxpricer.h"

namespace soplex
{

/**@brief   Dantzig pricer.
   @ingroup Algo

   Class SPxDantzigPR is an implementation class of an SPxPricer implementing
   Dantzig's default pricing strategy, i.e., maximal/minimal reduced cost or
   maximally violated constraint.

   See SPxPricer for a class documentation.
*/
class SPxDantzigPR : public SPxPricer
{
public:

   //-------------------------------------
   /**@name Constructors / destructors */
   //@{
   /// default constructor
   SPxDantzigPR() 
      : SPxPricer("Dantzig")
   {}   
   /// destructor
   virtual ~SPxDantzigPR()
   {}
   //@}

   //-------------------------------------
   /**@name Select enter/leave */
   //@{
   ///
   virtual int selectLeave();
   ///
   virtual SPxId selectEnter();
   //@}
};
} // namespace soplex
#endif // _SPXDEFAULTPRR_H_

//-----------------------------------------------------------------------------
//Emacs Local Variables:
//Emacs mode:c++
//Emacs c-basic-offset:3
//Emacs tab-width:8
//Emacs indent-tabs-mode:nil
//Emacs End:
//-----------------------------------------------------------------------------
