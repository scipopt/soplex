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
#pragma ident "@(#) $Id: spxharrisrt.h,v 1.6 2001/12/25 14:25:56 bzfkocht Exp $"

/**@file  spxharrisrt.h
 * @brief Harris pricing with shifting.
 */
#ifndef _SPXHARRISRT_H_
#define _SPXHARRISRT_H_


#include <assert.h>

#include "spxratiotester.h"

namespace soplex
{

/**@brief   Harris pricing with shifting.
   @ingroup Algo
   
   Class #SPxHarrisRT is a stable implementation of a #SPxRatioTester class
   along the lines of Harris' two phase algorithm. Additionally it uses
   shifting of bounds in order to avoid cycling.

   See #SPxRatioTester for a class documentation.
*/
class SPxHarrisRT : public SPxRatioTester
{
protected:
   SoPlex* thesolver;

   int maxDelta(
      double* /*max*/,       ///< max abs value in upd
      double* val,           ///< initial and chosen value
      int num,               ///< # of indices in idx
      const int* idx,        ///< nonzero indices in upd
      const double* upd,     ///< update vector for vec
      const double* vec,     ///< current vector
      const double* low,     ///< lower bounds for vec
      const double* up,      ///< upper bounds for vec
      double delta,          ///< allowed bound violation
      double epsilon,        ///< what is 0?
      double infinity);      ///< what is $\infty$?

   int minDelta(
      double* /*max*/,       ///< max abs value in upd
      double* val,           ///< initial and chosen value
      int num,               ///< of indices in idx
      const int* idx,        ///< nonzero indices in upd
      const double* upd,     ///< update vector for vec
      const double* vec,     ///< current vector
      const double* low,     ///< lower bounds for vec
      const double* up,      ///< upper bounds for vec
      double delta,          ///< allowed bound violation
      double epsilon,        ///< what is 0?
      double infinity);      ///< what is $\infty$?

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
#endif // _SPXHARRISRT_H_

//-----------------------------------------------------------------------------
//Emacs Local Variables:
//Emacs mode:c++
//Emacs c-basic-offset:3
//Emacs tab-width:8
//Emacs indent-tabs-mode:nil
//Emacs End:
//-----------------------------------------------------------------------------
