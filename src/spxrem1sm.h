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
#pragma ident "@(#) $Id: spxrem1sm.h,v 1.10 2003/01/10 12:46:14 bzfkocht Exp $"

/**@file  spxrem1sm.h
 * @brief Remove singletons from LP.
 */
#ifndef _SPXREM1SM_H_
#define _SPXREM1SM_H_

#include <assert.h>

#include "spxdefines.h"
#include "spxsimplifier.h"

namespace soplex
{
/**@brief   LP simplifier for removing uneccessary row/columns
   @ingroup Algo

   This #SPxSimplifier only removes rows and columns.
   For example those who containing one nonzero value only. 
   Also empty rows and columns are removed.   
*/
class SPxRem1SM : public SPxSimplifier
{
private:
   DVector        m_prim;       ///< unsimplified primal solution vector.
   DVector        m_dual;       ///< unsimplified dual solution vector.
   DataArray<int> m_cperm;      ///< column permutation vector.
   DataArray<int> m_rperm;      ///< row permutation vector.
   DSVector       m_pval;       ///< fixed variable values.
   Real           m_epsilon;    ///< epsilon zero
   Real           m_delta;      ///< maximum bound violation
private:
   ///
   void fixColumn(SPxLP& lp, int i);
   ///
   bool removeRows(SPxLP& lp, DataArray<int>& rem, int num, const char* msg);
   ///
   bool removeCols(SPxLP& lp, DataArray<int>& rem, int num, const char* msg);
   ///
   Result redundantRows(SPxLP& lp, bool& again);
   ///
   Result redundantCols(SPxLP& lp, bool& again);
   ///
   Result simpleRows(SPxLP& lp, bool& again);
   ///
   Result simpleCols(SPxLP& lp, bool& again);
   ///
   Real epsZero() const
   {
      return m_epsilon;
   }
   ///
   Real deltaBnd() const
   {
      return m_delta;
   }

public:
   /// default constructor
   SPxRem1SM() 
      : SPxSimplifier("Rem1")
   {}   
   /// destructor.
   virtual ~SPxRem1SM()
   {}  
   /// Remove singletons from the LP.
   virtual Result simplify(SPxLP& lp, Real eps, Real delta);
   /// returns a reference to the unsimplified primal solution.
   virtual const Vector& unsimplifiedPrimal(const Vector& x);
   /// returns a reference to the unsimplified dual solution. 
   virtual const Vector& unsimplifiedDual(const Vector& pi);
};
} // namespace soplex
#endif // _SPXREM1SM_H_

//-----------------------------------------------------------------------------
//Emacs Local Variables:
//Emacs mode:c++
//Emacs c-basic-offset:3
//Emacs tab-width:8
//Emacs indent-tabs-mode:nil
//Emacs End:
//-----------------------------------------------------------------------------
