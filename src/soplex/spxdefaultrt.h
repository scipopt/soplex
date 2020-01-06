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

/**@file  spxdefaultrt.h
 * @brief Textbook ratio test for SoPlex.
 */
#ifndef _SPXDEFAULTRT_H_
#define _SPXDEFAULTRT_H_


#include <assert.h>

#include "soplex/spxdefines.h"
#include "soplex/spxratiotester.h"

namespace soplex
{

/**@brief   Textbook ratio test for SoPlex.
   @ingroup Algo

   Class SPxDefaultRT provides an implementation of the textbook ratio test
   as a derived class of SPxRatioTester. This class is not intended for
   reliably solving LPs (even though it does the job for ``numerically simple''
   LPs). Instead, it should serve as a demonstration of how to write ratio
   tester classes.

   See SPxRatioTester for a class documentation.
*/
template <class R>
class SPxDefaultRT : public SPxRatioTester<R>
{
public:

   //-------------------------------------
   /**@name Construction / destruction */
   ///@{
   /// default constructor
   SPxDefaultRT()
      : SPxRatioTester<R>("Default")
   {}
   /// copy constructor
   SPxDefaultRT(const SPxDefaultRT& old)
      : SPxRatioTester<R>(old)
   {}
   /// assignment operator
   SPxDefaultRT& operator=(const SPxDefaultRT& rhs)
   {
      if(this != &rhs)
      {
         SPxRatioTester<R>::operator=(rhs);
      }

      return *this;
   }
   /// destructor
   virtual ~SPxDefaultRT()
   {}
   /// clone function for polymorphism
   inline virtual SPxRatioTester<R>* clone() const
   {
      return new SPxDefaultRT(*this);
   }
   ///@}

   //-------------------------------------
   /**@name Select enter/leave */
   ///@{
   ///
   virtual int selectLeave(R& val, R, bool);
   ///
   virtual SPxId selectEnter(R& val, int, bool);
};

} // namespace soplex

#include "spxdefaultrt.hpp"

#endif // _SPXDEFAULTRT_H_
