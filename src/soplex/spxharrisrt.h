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

/**@file  spxharrisrt.h
 * @brief Harris pricing with shifting.
 */
#ifndef _SPXHARRISRT_H_
#define _SPXHARRISRT_H_

#include <assert.h>

#include "soplex/spxdefines.h"
#include "soplex/spxratiotester.h"

namespace soplex
{

/**@brief   Harris pricing with shifting.
   @ingroup Algo

   Class SPxHarrisRT is a stable implementation of a SPxRatioTester class
   along the lines of Harris' two phase algorithm. Additionally it uses
   shifting of bounds in order to avoid cycling.

   See SPxRatioTester for a class documentation.
*/
/**@todo HarrisRT leads to cycling in dcmulti.sub.lp */
template <class R>
class SPxHarrisRT : public SPxRatioTester<R>
{
private:

   //-------------------------------------
   /**@name Private helpers */
   ///@{
   ///
   R degenerateEps() const;

   ///
   int maxDelta(
      R* /*max*/,        ///< max abs value in \p upd
      R* val,            ///< initial and chosen value
      int num,              ///< number of indices in \p idx
      const int* idx,       ///< nonzero indices in \p upd
      const R* upd,      ///< update VectorBase<R> for \p vec
      const R* vec,      ///< current vector
      const R* low,      ///< lower bounds for \p vec
      const R* up,       ///< upper bounds for \p vec
      R epsilon          ///< what is 0?
   ) const;

   ///
   int minDelta(
      R* /*max*/,        ///< max abs value in \p upd
      R* val,            ///< initial and chosen value
      int num,              ///< of indices in \p idx
      const int* idx,       ///< nonzero indices in \p upd
      const R* upd,      ///< update VectorBase<R> for \p vec
      const R* vec,      ///< current vector
      const R* low,      ///< lower bounds for \p vec
      const R* up,       ///< upper bounds for \p vec
      R epsilon          ///< what is 0?
   ) const;
   ///@}

public:

   //-------------------------------------
   /**@name Construction / destruction */
   ///@{
   /// default constructor
   SPxHarrisRT()
      : SPxRatioTester<R>("Harris")
   {}
   /// copy constructor
   SPxHarrisRT(const SPxHarrisRT& old)
      : SPxRatioTester<R>(old)
   {}
   /// assignment operator
   SPxHarrisRT& operator=(const SPxHarrisRT& rhs)
   {
      if(this != &rhs)
      {
         SPxRatioTester<R>::operator=(rhs);
      }

      return *this;
   }
   /// destructor
   virtual ~SPxHarrisRT()
   {}
   /// clone function for polymorphism
   inline virtual SPxRatioTester<R>* clone() const
   {
      return new SPxHarrisRT(*this);
   }
   ///@}

   //-------------------------------------
   /**@name Leave / enter */
   ///@{
   ///
   virtual int selectLeave(R& val, R, bool);
   ///
   virtual SPxId selectEnter(R& val, int, bool);
   ///@}

};

} // namespace soplex
// For the general template
#include "spxharrisrt.hpp"


#endif // _SPXHARRISRT_H_
