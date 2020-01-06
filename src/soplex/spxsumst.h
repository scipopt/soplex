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


/**@file  spxsumst.h
 * @brief Simple heuristic SPxStarter.
 */
#ifndef _SPXSUMST_H_
#define _SPXSUMST_H_


#include <assert.h>

#include "soplex/spxvectorst.h"

namespace soplex
{

/**@brief   Simple heuristic SPxStarter.
   @ingroup Algo

   Testing version of an SPxVectorST using a very simplistic heuristic to
   build up an approximated solution vector.
*/
template <class R>
class SPxSumST : public SPxVectorST<R>
{
protected:

   //-------------------------------------
   /**@name Protected helpers */
   ///@{
   /// sets up variable weights.
   void setupWeights(SPxSolverBase<R>& base);
   ///@}

public:

   //-------------------------------------
   /**@name Construction / destruction */
   ///@{
   /// default constructor.
   SPxSumST()
   {
      this->m_name = "Sum";
   }
   /// copy constructor
   SPxSumST(const SPxSumST& old)
      : SPxVectorST<R>(old)
   {
      assert(this->isConsistent());
   }
   /// assignment operator
   SPxSumST& operator=(const SPxSumST& rhs)
   {
      if(this != &rhs)
      {
         SPxVectorST<R>::operator=(rhs);

         assert(this->isConsistent());
      }

      return *this;
   }
   /// destructor.
   virtual ~SPxSumST()
   {}
   /// clone function for polymorphism
   inline virtual SPxStarter<R>* clone() const
   {
      return new SPxSumST(*this);
   }
   ///@}

};

} // namespace soplex

#include "spxsumst.hpp"
#endif // _SPXSUMST_H_
