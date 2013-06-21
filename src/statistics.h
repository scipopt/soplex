/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the class library                   */
/*       SoPlex --- the Sequential object-oriented simPlex.                  */
/*                                                                           */
/*    Copyright (C) 1996-2013 Konrad-Zuse-Zentrum                            */
/*                            fuer Informationstechnik Berlin                */
/*                                                                           */
/*  SoPlex is distributed under the terms of the ZIB Academic Licence.       */
/*                                                                           */
/*  You should have received a copy of the ZIB Academic License              */
/*  along with SoPlex; see the file COPYING. If not email to soplex@zib.de.  */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/**@file  statistics.h
 * @brief Class for collecting statistical information
 */
#ifndef _STATISTICS_H_
#define _STATISTICS_H_

#include <iostream>

#include "soplex2.h"

namespace soplex
{
   /**@class   Statistics
    * @brief   Class for collecting statistical information
    * @ingroup Algo
    */
   class SoPlex2::Statistics
   {
   public:

      //**@name Construction, resetting, printing */
      //@{

      /// default constructor
      Statistics();

      /// clears statistics
      void clear();

      /// prints statistics
      void print(std::ostream& os);
      //@}


      //**@name Data */
      //@{

      Real solvingTime;
      int iterations;
      int refinements;
      int stallRefinements;

      //@}
   };
} // namespace soplex
#endif // _STATISTICS_H_

//-----------------------------------------------------------------------------
//Emacs Local Variables:
//Emacs mode:c++
//Emacs c-basic-offset:3
//Emacs tab-width:8
//Emacs indent-tabs-mode:nil
//Emacs End:
//-----------------------------------------------------------------------------
