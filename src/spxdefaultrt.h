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
#pragma ident "@(#) $Id: spxdefaultrt.h,v 1.5 2001/11/29 14:43:46 bzfpfend Exp $"

/**@file  spxdefaultrt.h
 * @brief Textbook ratio test for #SoPlex.
 */
#ifndef _SPXDEFAULTRT_H_
#define _SPXDEFAULTRT_H_


#include <assert.h>

#include "spxratiotester.h"

namespace soplex
{

/**@brief   Textbook ratio test for #SoPlex.
   @ingroup Algo
   
   Class #SPxDefaultRT provides an implementation of the textbook ratio test
   as a derived class of #SPxRatioTester. This class is not intended for
   reliably solving LPs (even though it does the job for ``numerically simple''
   LPs). Instead, it should serve as a demonstration of how to write ratio
   tester classes.

   See #SPxRatioTester for a class documentation.
*/
class SPxDefaultRT : public SPxRatioTester
{
protected:
   SoPlex* thesolver;

   ///
   int selectLeave(double& val, int start, int incr);

   ///
   SoPlex::Id selectEnter(double& val, int start1, int incr1, int start2,
                          int incr2);

public:
   ///
   SoPlex* solver() const
   {
      return thesolver;
   }

   ///
   void load(SoPlex* p_solver)
   {
      thesolver = p_solver;
   }

   ///
   void clear()
   {
      thesolver = 0;
   }

   ///
   int selectLeave(double& val);

   ///
   SoPlex::Id selectEnter(double& val);

   ///
   void setType(SoPlex::Type)
   {}

};

} // namespace soplex
#endif // _SPXDEFAULTRT_H_

//-----------------------------------------------------------------------------
//Emacs Local Variables:
//Emacs mode:c++
//Emacs c-basic-offset:3
//Emacs tab-width:8
//Emacs indent-tabs-mode:nil
//Emacs End:
//-----------------------------------------------------------------------------
