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
#pragma ident "@(#) $Id: spxgeneralsm.h,v 1.10 2003/01/05 19:03:16 bzfkocht Exp $"

/**@file  spxgeneralsm.h
 * @brief General LP preprocessing.
 */
#ifndef _SPXGENERALSM_H_
#define _SPXGENERALSM_H_

#include <assert.h>

#include "spxdefines.h"
#include "spxintervalsm.h"
#include "spxaggregatesm.h"
#include "spxrem1sm.h"

namespace soplex
{
/**@brief General LP preprocessing.
   @ingroup Algo

   This #SPxSimplifier iterativly applies a number of preprocessors to its
   loaded #SPxLP.
*/
class SPxGeneralSM : public SPxSimplifier
{
private:
   SPxIntervalSM  m_inter;    ///< remove too small values.
   SPxRem1SM      m_rem1;     ///< remove row/column singletons.
   SPxAggregateSM m_aggr;     ///< do variable aggregation.
   Real           m_repth;    ///< repetition threshold.

public:
   /// default constructor
   explicit SPxGeneralSM(Real repth = 0.95) 
      : SPxSimplifier("General")
      , m_repth(repth)
   {}   
   /// destructor.
   virtual ~SPxGeneralSM()
   {}  
   /// Simplify #SPxLP.
   virtual Result simplify(SPxLP& lp);
   /// returns a reference to the unsimplified primal solution.
   virtual const Vector& unsimplifiedPrimal(const Vector& x);
   /// returns a reference to the unsimplified dual solution. 
   virtual const Vector& unsimplifiedDual(const Vector& pi);
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
