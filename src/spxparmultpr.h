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
#pragma ident "@(#) $Id: spxparmultpr.h,v 1.18 2009/02/20 01:06:37 bzfgleix Exp $"

/**@file  spxparmultpr.h
 * @brief Partial multiple pricing.
 */
#ifndef _SPXPARMULTPR_H_
#define _SPXPARMULTPR_H_

#include <assert.h>

#include "spxdefines.h"
#include "spxpricer.h"
#include "dataarray.h"
#include "array.h"
#include "ssvector.h"

namespace soplex
{

/**@brief   Partial multiple pricing.
   @ingroup Algo

   Class SPxParMultPr is an implementation class for SPxPricer implementing
   Dantzig's default pricing strategy with partial multiple pricing.
   Partial multiple pricing applies to the entering Simplex only. A set of
   #partialSize eligible pivot indices is selected (partial pricing). In the
   following Simplex iterations pricing is restricted to these indices
   (multiple pricing) until no more eliiable pivots are available. Partial
   multiple pricing significantly reduces the computation time for computing
   the matrix-vector-product in the Simplex algorithm.

   See SPxPricer for a class documentation.
*/
class SPxParMultPR : public SPxPricer
{
private:

   //-------------------------------------
   /**@name Private types */
   //@{
   /// Helper structure.
   struct SPxParMultPr_Tmp
   {
      ///
      SPxId id;
      ///
      Real test;
   };
   //@}

   //-------------------------------------
   /**@name Helper data */
   //@{
   ///
   DataArray < SPxParMultPr_Tmp > pricSet;
   ///
   int multiParts;
   ///
   int used;
   ///
   int min;
   ///
   int last;
   /// Set size for partial pricing.
   static int partialSize;
   //@}

public:

   //-------------------------------------
   /**@name Construction / destruction */
   //@{
   /// default constructor
   SPxParMultPR() 
      : SPxPricer("ParMult")
   {}   
   /// destructor
   virtual ~SPxParMultPR()
   {}
   //@}

   //-------------------------------------
   /**@name Interface */
   //@{
   /// set the solver
   virtual void load(SPxSolver* solver);
   /// set entering or leaving algorithm
   virtual void setType(SPxSolver::Type tp);
   /// 
   virtual int selectLeave();
   ///
   virtual SPxId selectEnter();
   //@}

   //-------------------------------------
   /**@name Blocked */
   //@{
   /// copy constructor
   SPxParMultPR( const SPxParMultPR& );
   /// assignment operator
   SPxParMultPR& operator=( const SPxParMultPR& );
   //@}
};


} // namespace soplex
#endif // _SPXPARMULTPRR_H_

//-----------------------------------------------------------------------------
//Emacs Local Variables:
//Emacs mode:c++
//Emacs c-basic-offset:3
//Emacs tab-width:8
//Emacs indent-tabs-mode:nil
//Emacs End:
//-----------------------------------------------------------------------------
