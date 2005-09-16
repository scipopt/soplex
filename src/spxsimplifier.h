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
#pragma ident "@(#) $Id: spxsimplifier.h,v 1.14 2005/09/16 12:42:36 bzfhille Exp $"

/**@file  spxsimplifier.h
 * @brief LP simplification base class.
 */
#ifndef _SPXSIMPLIFIER_H_
#define _SPXSIMPLIFIER_H_

#include <assert.h>

#include "spxdefines.h"
#include "timer.h"
#include "spxlp.h"

namespace soplex
{
/**@brief   LP simplification abstract base class.
   @ingroup Algo

   Instances of classes derived from #SPxSimplifier may be loaded to #SoPlex in
   order to simplify LPs before solving them. #SoPlex# will call #simplify()
   on its self. Generally any #SPxLP can be given to 
   a #SPxSimplifier for #simplify()%ing it. The simplification can not be undone,
   but given an primal/dual solution for the simplified #SPxLP, the simplifier
   can reconstruct the primal/dual solution of the unsimplified LP.
*/
class SPxSimplifier
{
protected:
   const char* m_name;
   Timer       m_timeUsed;
   int         m_remRows;
   int         m_remCols;
   int         m_remNzos;
   int         m_chgBnds;
   int         m_chgLRhs;

public:
   enum Result
   {
      OKAY       =  0,
      INFEASIBLE =  1,
      UNBOUNDED  =  2,
      VANISHED   =  3
   };
   /// constructor
   explicit SPxSimplifier(const char* p_name)
      : m_name(p_name)
      , m_remRows(0)
      , m_remCols(0)
      , m_remNzos(0)
      , m_chgBnds(0)
      , m_chgLRhs(0)
   {}
   /// destructor.
   virtual ~SPxSimplifier()
   {
      m_name = 0;
   }
   /// get name of simplifier.
   virtual const char* getName() const
   {
      return m_name;
   }
   virtual Real timeUsed() const
   {
      return m_timeUsed.userTime();
   }
   /// simplify #SPxLP \p lp. 
   /**
    * @return 
    *  <TABLE>
    *  <TR><TD>#OKAY      </TD><TD>if this could be done,</TD></TR>
    *  <TR><TD>#UNBOUNDED </TD><TD>if primal unboundedness was detected or</TD></TR>
    *  <TR><TD>#INFEASIBLE</TD><TD>if primal infeasibility was detected.</TD></TR>
    *  </TABLE>
    */
   virtual Result simplify(SPxLP& lp, Real eps, Real delta) = 0;

   /// returns a reference to the unsimplified primal solution.
   virtual const Vector& unsimplifiedPrimal(const Vector& x)
   {
      return x;
   }
   /// returns a reference to the unsimplified dual solution. 
   virtual const Vector& unsimplifiedDual(const Vector& pi)
   {
      return pi;
   }
#ifndef NO_CONSISTENCY_CHECKS
   /// consistency check
   virtual bool isConsistent() const
   {
      return true;
   }
#endif
};
} // namespace soplex
#endif // _SPXSIMPLIFIER_H_

//-----------------------------------------------------------------------------
//Emacs Local Variables:
//Emacs mode:c++
//Emacs c-basic-offset:3
//Emacs tab-width:8
//Emacs indent-tabs-mode:nil
//Emacs End:
//-----------------------------------------------------------------------------

