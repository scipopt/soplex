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
#pragma ident "@(#) $Id: spxredundantsm.h,v 1.16 2005/11/08 19:56:52 bzforlow Exp $"

/**@file  spxredundantsm.h
 * @brief Remove singletons from LP.
 */
#ifndef _SPXREDUNDANTSM_H_
#define _SPXREDUNDANTSM_H_

#include <assert.h>

#include "spxdefines.h"
#include "spxsimplifier.h"

namespace soplex
{
/**@brief   LP simplifier for removing unnecessary rows or columns
   @ingroup Algo

   This SPxSimplifier removes redundant rows and columns and bounds.
   Also infeasibility or unboundedness may be detected.
   Removed are:
   - empty rows or columns
   - unconstrained constraints
   - row singletons
   - fixed columns
   - some columns singletons
   - rows with all fixed variables due to implied bounds
   - redundant rhs/lhs
   - redundant column bounds
   - duplicate rows
   - columns which are implicitly fixed to a bound
   - columns with redundant bounds
*/
class SPxRedundantSM : public SPxSimplifier
{
private:

   //------------------------------------
   //**@name Types */
   //@{
   /// comparison operator for rows
   class RowHash
   {
   public: 
      /// constructor
      RowHash() {}
      /// destructor
      ~RowHash() {}
      /// comparison operator
      int operator() (const RowHash& rh1, const RowHash& rh2) const
      {
         // rh1.hid - rh2.hid is a bad idea, because they are unsigned

         if (rh1.hid < rh2.hid)
            return -1;
         if (rh1.hid > rh2.hid)
            return  1;

         assert(rh1.hid == rh2.hid);

         return rh1.row - rh2.row;
      }

      int          row;  ///< row no.
      unsigned int hid;  ///< hash id
   };
   //@}

private:

   //------------------------------------
   //**@name Data */
   //@{
   /// 
   DVector        m_prim;       ///< unsimplified primal solution vector.
   DVector        m_dual;       ///< unsimplified dual solution vector.
   DataArray<int> m_cperm;      ///< column permutation vector.
   DataArray<int> m_rperm;      ///< row permutation vector.
   DSVector       m_pval;       ///< fixed variable values.
   Real           m_epsilon;    ///< epsilon zero
   Real           m_delta;      ///< maximum bound violation
   //@}

private:

   //------------------------------------
   //**@name Private helpers */
   //@{
   ///
   void fixColumn(SPxLP& lp, int i);
   ///
   void removeRows(SPxLP& lp, DataArray<int>& rem, int num);
   ///
   void removeCols(SPxLP& lp, DataArray<int>& rem, int num);
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
   //@}

public:

   //------------------------------------
   //**@name Constructors / destructors */
   //@{
   /// default constructor
   SPxRedundantSM() 
      : SPxSimplifier("Redundant")
   {}   
   /// destructor.
   virtual ~SPxRedundantSM()
   {}  
   //@}

   //------------------------------------
   //**@name LP simplification */
   //@{
   /// Remove singletons from the LP.
   virtual Result simplify(SPxLP& lp, Real eps, Real delta);
   /// returns a reference to the unsimplified primal solution.
   virtual const Vector& unsimplifiedPrimal(const Vector& x);
   /// returns a reference to the unsimplified dual solution. 
   virtual const Vector& unsimplifiedDual(const Vector& pi);
   //@}

   //------------------------------------
   //**@name Blocked */
   //@{
   /// copy constructor
   SPxRedundantSM( const SPxRedundantSM& );
   /// assignment operator
   SPxRedundantSM& operator=( const SPxRedundantSM& );
   //@}
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
