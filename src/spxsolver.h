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
#pragma ident "@(#) $Id: spxsolver.h,v 1.2 2001/11/06 23:31:05 bzfkocht Exp $"


#ifndef _SPXSOLVER_H_
#define _SPXSOLVER_H_

//@ ----------------------------------------------------------------------------
/*  \Section{Imports}
    Import required system include files ...
 */
#include <assert.h>


/*  ... and class header files
 */

#include "spxsteeppr.h"
#include "spxfastrt.h"
#include "spxweightst.h"
#include "slufactor.h"

namespace soplex
{

//@ ----------------------------------------------------------------------------
/* \Section{Class Declaration}
 */

/// preconfigured \Ref{SoPlex} LP-solver.
class SPxSolver : public SoPlex
{
private:
   SPxFastRT rt;
   SPxSteepPR pr;
   SPxWeightST st;
   SLUFactor slu;

public:
   SPxSolver (Type type = LEAVE, SoPlex::Representation rep = SoPlex::ROW);
};

} // namespace soplex
#endif // _SPXSOLVER_H_

//-----------------------------------------------------------------------------
//Emacs Local Variables:
//Emacs mode:c++
//Emacs c-basic-offset:3
//Emacs tab-width:8
//Emacs indent-tabs-mode:nil
//Emacs End:
//-----------------------------------------------------------------------------
