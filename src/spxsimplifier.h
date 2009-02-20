/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the class library                   */
/*       SoPlex --- the Sequential object-oriented simPlex.                  */
/*                                                                           */
/*    Copyright (C) 1997-1999 Roland Wunderling                              */
/*                  1997-2009 Konrad-Zuse-Zentrum                            */
/*                            fuer Informationstechnik Berlin                */
/*                                                                           */
/*  SoPlex is distributed under the terms of the ZIB Academic Licence.       */
/*                                                                           */
/*  You should have received a copy of the ZIB Academic License              */
/*  along with SoPlex; see the file COPYING. If not email to soplex@zib.de.  */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
#pragma ident "@(#) $Id: spxsimplifier.h,v 1.21 2009/02/20 01:06:37 bzfgleix Exp $"

/**@file  spxsimplifier.h
 * @brief LP simplification base class.
 */
#ifndef _SPXSIMPLIFIER_H_
#define _SPXSIMPLIFIER_H_

#include <assert.h>

#include "spxdefines.h"
#include "timer.h"
#include "spxlp.h"
#include "spxsolver.h"

namespace soplex
{
/**@brief   LP simplification abstract base class.
   @ingroup Algo

   Instances of classes derived from SPxSimplifier may be loaded to SoPlex in
   order to simplify LPs before solving them. SoPlex will call #simplify()
   on itself. Generally any SPxLP can be given to 
   a SPxSimplifier for #simplify()%ing it. The simplification cannot be undone,
   but given an primal/dual solution for the simplified SPxLP, the simplifier
   can reconstruct the primal/dual solution of the unsimplified LP.
*/
class SPxSimplifier
{
protected:

   //-------------------------------------
   /**@name Protected Data */
   //@{
   /// name of the simplifier
   const char* m_name;
   /// user time used for simplification
   Timer       m_timeUsed;
   /// number of removed rows
   int         m_remRows;
   /// number of removed columns
   int         m_remCols;
   /// number of removed nonzero coefficients
   int         m_remNzos;
   /// number of changed bounds
   int         m_chgBnds;
   /// number of change right-hand sides
   int         m_chgLRhs;
   //@}

public:

   //-------------------------------------
   /**@name Types */
   //@{
   /// Result of the simplification.
   enum Result
   {
      OKAY            =  0,  ///< simplification could be done
      INFEASIBLE      =  1,  ///< primal infeasibility was detected
      DUAL_INFEASIBLE =  2,  ///< dual infeasibility was detected
      UNBOUNDED       =  3,  ///< primal unboundedness was detected
      VANISHED        =  4   ///< the problem was so much simplified that it vanished
   };
   //@}

   //-------------------------------------
   /**@name Types */
   //@{
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
   //@}

   //-------------------------------------
   /**@name Access / modfication */
   //@{
   /// get name of simplifier.
   virtual const char* getName() const
   {
      return m_name;
   }
   virtual Real timeUsed() const
   {
      return m_timeUsed.userTime();
   }
   //@}

   //-------------------------------------
   /**@name Simplifying / unsimplifying */
   //@{
   /// simplify SPxLP \p lp. 
   virtual Result simplify(SPxLP& lp, Real eps, Real delta) = 0;
   /// reconstructs an optimal solution for the unsimplified LP.
   virtual void unsimplify(const Vector&, const Vector&, const Vector&, const Vector&,
                           const SPxSolver::VarStatus[], const SPxSolver::VarStatus[]) {}
   /// specifies whether an optimal solution has already been unsimplified.
   virtual bool isUnsimplified() const
   {
      return false;
   }
   /// returns a reference to the unsimplified primal solution.
   virtual const Vector& unsimplifiedPrimal() = 0;

   /// returns a reference to the unsimplified dual solution.
   virtual const Vector& unsimplifiedDual() = 0;

   /// returns a reference to the unsimplified slack values.
   virtual const Vector& unsimplifiedSlacks() = 0;

   /// returns a reference to the unsimplified reduced costs.
   virtual const Vector& unsimplifiedRedCost() = 0;

   /// gets basis status for a single row.
   virtual SPxSolver::VarStatus getBasisRowStatus(int) const = 0;

   /// gets basis status for a single column.
   virtual SPxSolver::VarStatus getBasisColStatus(int) const = 0;

   /// get optimal basis.
   virtual void getBasis(SPxSolver::VarStatus[], SPxSolver::VarStatus[]) const = 0;
   //@}

#ifndef NO_CONSISTENCY_CHECKS
   //-------------------------------------
   /**@name Consistency check */
   //@{
   /// consistency check
   virtual bool isConsistent() const
   {
      return true;
   }
   //@}
#endif

private:

   //-------------------------------------
   /**@name Blocked */
   //@{
   /// copy constructor
   SPxSimplifier( const SPxSimplifier& );
   /// assignment operator
   SPxSimplifier& operator=( const SPxSimplifier& );
   //@}
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

