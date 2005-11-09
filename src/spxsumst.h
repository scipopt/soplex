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
#pragma ident "@(#) $Id: spxsumst.h,v 1.8 2005/11/09 13:53:50 bzforlow Exp $"


/**@file  spxsumst.h
 * @brief Simple heuristic SPxStarter.
 */
#ifndef _SPXSUMST_H_
#define _SPXSUMST_H_


#include <assert.h>

#include "spxvectorst.h"

namespace soplex
{

/**@brief   Simple heuristic SPxStarter.
   @ingroup Algo

   Testing version of an SPxVectorST using a very simplistic heuristic to
   build up an approximated solution vector.
*/
class SPxSumST : public SPxVectorST
{
protected:

   //-------------------------------------
   /**@name Protected helpers */
   //@{
   /// sets up variable weights.
   void setupWeights(SPxSolver& base);
   //@}

public:

   //-------------------------------------
   /**@name Construction / destruction */
   //@{
   /// default constructor.
   SPxSumST()
   {
      m_name = "Sum";
   }
   /// destructor.
   virtual ~SPxSumST()
   {}  
   //@}

private:

   //-------------------------------------
   /**@name Blocked */
   //@{
   /// copy constructor
   SPxSumST( const SPxSumST& );
   /// assignment operator
   SPxSumST& operator=( const SPxSumST& );
   //@}
};

} // namespace soplex
#endif // _SPXSUMST_H_

//-----------------------------------------------------------------------------
//Emacs Local Variables:
//Emacs mode:c++
//Emacs c-basic-offset:3
//Emacs tab-width:8
//Emacs indent-tabs-mode:nil
//Emacs End:
//-----------------------------------------------------------------------------


