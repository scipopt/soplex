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
#pragma ident "@(#) $Id: spxgeneralsm.h,v 1.9 2002/04/14 12:41:54 bzfkocht Exp $"

/**@file  spxgeneralsm.h
 * @brief General LP preprocessing.
 */
#ifndef _SPXGENERALSM_H_
#define _SPXGENERALSM_H_

#include <assert.h>

#include "spxdefines.h"
#include "spxredundantsm.h"
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
   SPxRem1SM      m_rem1;   ///< remove row/column singletons .
   SPxRedundantSM m_redu;   ///< remove redundant rows/columns.
   SPxAggregateSM m_aggr;   ///< do variable aggregation.
   Real           m_repth;  ///< repetition threashold.

public:
   /// Load the #SPxLP to be simplified.
   virtual void load(SPxLP* p_lp);

   /// Unload the #SPxLP.
   virtual void unload();

   /// Simplify #SPxLP.
   virtual int simplify();

   /// Reverse what #simplify() had done.
   virtual void unsimplify();

   /// objective value for unsimplified LP.
   virtual Real value(Real x);

   /// default constructor.
   SPxGeneralSM(Real repth = 0.95)
      : m_repth(repth)
   {}
  
   /// destructor.
   virtual ~SPxGeneralSM()
   {}  
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
