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
#pragma ident "@(#) $Id: spxgeneralsm.h,v 1.2 2001/11/06 23:31:04 bzfkocht Exp $"


#ifndef _SPXGENERALSM_H_
#define _SPXGENERALSM_H_

//@ ----------------------------------------------------------------------------
/*  \Section{Imports}
    Import required system include files ...
 */
#include <assert.h>


/*  ... and class header files
 */

#include "spxredundantsm.h"
#include "spxaggregatesm.h"
#include "spxrem1sm.h"
#include "spxscale.h"

namespace soplex
{



//@ ----------------------------------------------------------------------------
/* \Section{Class Declaration}
 */

/** General LP preprocessing.
    This \Ref{SPxSimplifier} iterativly applies a number of preprocessors to its
    loaded \Ref{SPxLP}.
 */
class SPxGeneralSM : public SPxSimplifier
{
private:
   SPxRem1SM rem1;
   SPxRedundantSM redu;
   SPxAggregateSM aggr;
   SPxScale scale;

   SPxLP* lp;

public:
   ///
   virtual void load(SPxLP*);
   ///
   virtual void unload();
   ///
   virtual SPxLP* loadedLP() const;
   ///
   virtual int simplify();
   ///
   virtual void unsimplify();
   ///
   virtual double value(double x);
   ///
   SPxGeneralSM()
   {}
  ;

   virtual ~SPxGeneralSM()
   {}
  ;
};

} // namespace soplex
#endif // _SPXGENERALSM_H_

//-----------------------------------------------------------------------------
//Emacs Local Variables:
//Emacs mode:c++
//Emacs c-basic-offset:3
//Emacs tab-width:8
//Emacs indent-tabs-mode:nil
//Emacs End:
//-----------------------------------------------------------------------------
