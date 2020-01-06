/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the class library                   */
/*       SoPlex --- the Sequential object-oriented simPlex.                  */
/*                                                                           */
/*    Copyright (C) 1996-2020 Konrad-Zuse-Zentrum                            */
/*                            fuer Informationstechnik Berlin                */
/*                                                                           */
/*  SoPlex is distributed under the terms of the ZIB Academic Licence.       */
/*                                                                           */
/*  You should have received a copy of the ZIB Academic License              */
/*  along with SoPlex; see the file COPYING. If not email to soplex@zib.de.  */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */


/**@file  spxsteepexpr.h
 * @brief Steepest edge pricer with exact initialization of weights.
 */
#ifndef _SPXSTEEPEXPR_H_
#define _SPXSTEEPEXPR_H_


#include <assert.h>

#include "soplex/spxdefines.h"
#include "soplex/spxsteeppr.h"

namespace soplex
{

/**@brief   Steepest edge pricer.
   @ingroup Algo

   Class SPxSteepExPR implements a steepest edge pricer to be used with
   SoPlex. Exact initialization of weights is used.

   See SPxPricer for a class documentation.
*/
template <class R>
class SPxSteepExPR : public SPxSteepPR<R>
{

public:

   //-------------------------------------
   /**@name Construction / destruction */
   ///@{
   ///
   SPxSteepExPR()
      : SPxSteepPR<R>("SteepEx", SPxSteepPR<R>::EXACT)
   {
      assert(this->isConsistent());
   }
   /// copy constructor
   SPxSteepExPR(const SPxSteepExPR& old)
      : SPxSteepPR<R>(old)
   {
      assert(this->isConsistent());
   }
   /// assignment operator
   SPxSteepExPR& operator=(const SPxSteepExPR& rhs)
   {
      if(this != &rhs)
      {
         SPxSteepPR<R>::operator=(rhs);

         assert(this->isConsistent());
      }

      return *this;
   }
   /// destructor
   virtual ~SPxSteepExPR()
   {}
   /// clone function for polymorphism
   inline virtual SPxSteepPR<R>* clone()  const
   {
      return new SPxSteepExPR(*this);
   }
   ///@}
};

} // namespace soplex

#include "spxsteeppr.hpp"

#endif // _SPXSTEEPPR_H_
