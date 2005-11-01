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
#pragma ident "@(#) $Id: spxaggregatesm.h,v 1.14 2005/11/01 15:35:40 bzforlow Exp $"

/**@file  spxaggregatesm.h
 * @brief LP variable aggregation.
 */
#ifndef _SPXAGGREGATESM_H_
#define _SPXAGGREGATESM_H_

#include <assert.h>

#include "spxdefines.h"
#include "spxsimplifier.h"

namespace soplex
{
/** @class SPxAggregateSM
    @ingroup Algo

    @brief This #SPxSimplifier does variable aggregation.
 */
class SPxAggregateSM : public SPxSimplifier
{
public:

   //------------------------------------
   //**@name Types */
   //@{
   /// 
   struct RowCnt
   {
      int row;
      int size;
   };
   /// Compares two RowCnts.
   struct Compare
   {
   public:
      /// constructor
      Compare() {}
      /// comparison
      inline
      int operator()( const RowCnt& i1, const RowCnt& i2) const
      {
         return i1.size - i2.size;
      }
   };
   //@}

private:

   //------------------------------------
   //**@name Data */
   //@{
   Real stability;   ///< stability factor, e.g. 0.01.   
   Real maxFill;     ///< ???  
   //@}

   //------------------------------------
   //**@name Private helpers */
   //@{
   /// ???
   int eliminate(SPxLP& lp, const SVector& row, Real b);
   //@}

public:

   //------------------------------------
   //**@name Construction / destruction */
   //@{
   /// default constructor
   SPxAggregateSM() 
      : SPxSimplifier("Aggregate")
   {}   
   /// destructor.
   virtual ~SPxAggregateSM()
   {}  
   //@}

   //------------------------------------
   //**@name Aggregation */
   //@{
   /// Aggregate variables.
   virtual Result simplify(SPxLP& lp, Real /*eps*/, Real /*delta*/);
   /// returns a reference to the unsimplified primal solution.
   virtual const Vector& unsimplifiedPrimal(const Vector& x);
   /// returns a reference to the unsimplified dual solution. 
   virtual const Vector& unsimplifiedDual(const Vector& pi);
   //@}
};
} // namespace soplex
#endif // _SPXAGGREGATESM_H_

//-----------------------------------------------------------------------------
//Emacs Local Variables:
//Emacs mode:c++
//Emacs c-basic-offset:3
//Emacs tab-width:8
//Emacs indent-tabs-mode:nil
//Emacs End:
//-----------------------------------------------------------------------------
