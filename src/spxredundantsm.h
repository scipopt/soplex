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
#pragma ident "@(#) $Id: spxredundantsm.h,v 1.12 2003/01/05 19:03:17 bzfkocht Exp $"

/**@file  spxredundantsm.h
 * @brief Remove redundant row and columns.
 */
#ifndef _SPXREDUNDANTSM_H_
#define _SPXREDUNDANTSM_H_

#include <assert.h>

#include "spxdefines.h"
#include "spxsimplifier.h"

namespace soplex
{
/**@brief   Remove redundant row and columns.
   @ingroup Algo

   This #SPxSimplifier tries to eliminate redundant rows and columns from
   its loaded #SPxLP.
 */
class SPxRedundantSM : public SPxSimplifier
{
private:
   DVector        prim;   ///< unsimplified primal solution vector.
   DVector        dual;   ///< unsimplified dual solution vector.
   DataArray<int> cperm;  ///< column permutation vector.
   DataArray<int> rperm;  ///< row permutation vector.

private:
   Result treat_cols(SPxLP& lp);
   Result treat_rows(SPxLP& lp);

public:
   // default constructor
   SPxRedundantSM() 
      : SPxSimplifier("Redundant")
   {}   
   /// destructor.
   virtual ~SPxRedundantSM()
   {}  
   /// Remove redundant rows and columns.
   virtual Result simplify(SPxLP& lp);
   /// returns a reference to the unsimplified primal solution.
   virtual const Vector& unsimplifiedPrimal(const Vector& x);
   /// returns a reference to the unsimplified dual solution. 
   virtual const Vector& unsimplifiedDual(const Vector& pi);
};
} // namespace soplex
#endif // _SPXREDUNDANTSM_H_

//-----------------------------------------------------------------------------
//Emacs Local Variables:
//Emacs mode:c++
//Emacs c-basic-offset:3
//Emacs tab-width:8
//Emacs indent-tabs-mode:nil
//Emacs End:
//-----------------------------------------------------------------------------
