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
#pragma ident "@(#) $Id: spxsumst.h,v 1.2 2001/11/06 23:31:05 bzfkocht Exp $"


#ifndef _SPXSUMST_H_
#define _SPXSUMST_H_

//@ ----------------------------------------------------------------------------
/*  \Section{Imports}
    Import required system include files ...
 */
#include <assert.h>


/*  ... and class header files
 */

#include "spxvectorst.h"

namespace soplex
{





//@ ----------------------------------------------------------------------------
/* \Section{Class Declaration}
 */

/** Simple heuristic #SPxStarter#.
    Testing version of an #SPxVectorST# using a very simplistic heuristic to
    build up an approximated solution vector.
 */
class SPxSumST : public SPxVectorST
{
private:
protected:
   void setupWeights(SoPlex& base);
public:
   ///
   SPxSumST()
   {}

};

} // namespace soplex
#endif // _SPXSUMST_H_

//-----------------------------------------------------------------------------
//Emacs Local Variables:
//Emacs mode:c++
//Emacs c-basic-offset:3
//Emacs tab-width:8
//Emacs indent-tabs-mode:nil
//Emacs End:
//-----------------------------------------------------------------------------
