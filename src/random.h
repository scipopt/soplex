/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the class library                   */
/*       SoPlex --- the Sequential object-oriented simPlex.                  */
/*                                                                           */
/*    Copyright (C) 1996-2016 Konrad-Zuse-Zentrum                            */
/*                            fuer Informationstechnik Berlin                */
/*                                                                           */
/*  SoPlex is distributed under the terms of the ZIB Academic Licence.       */
/*                                                                           */
/*  You should have received a copy of the ZIB Academic License              */
/*  along with SoPlex; see the file COPYING. If not email to soplex@zib.de.  */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/**@file  random.h
 * @brief Random numbers.
 */
#ifndef _RANDOM_H_
#define _RANDOM_H_

#include <assert.h>

namespace soplex
{
#define  RSTEP    1103515245
#define  RDIVIDE  65536
#define  RADD     12345
#define  RMULT    32768

/**@brief   Random numbers.
   @ingroup Elementary

   Class Random provides random Real variables, i.e. a value variable that
   gives another value each time it is accessed. It may be used just like an
   ordinary Real by means of an overloaded cast operator Real()%.
   This is largely the same implementation as rand() from std lib.
*/
class Random
{
private:

   //--------------------------------------
   /**@name Data */
   //@{
   unsigned int seed;    ///< random seed.
   //@}

   //--------------------------------------
   /**@name Helpers */
   //@{
   /// increases rand seed and returns a pseudo random Real value in [0,1).
   Real next_random()
   {
      seed = seed * RSTEP + RADD;
      Real i = int ((seed / RDIVIDE) % RMULT);
      return ( i / RMULT );
   }

   //@}

public:

   //--------------------------------------
   /**@name Access */
   //@{
   /// returns next random number.
   Real next(Real minimum = 0.0, Real maximum = 1.0)
   {
      return (minimum + (maximum - minimum) * next_random());
   }

   /// returns current seed
   unsigned int getSeed() const
   {
      return seed;
   }

   //@}

   //--------------------------------------
   /**@name Modification */
   //@{
   /// resets seed for next random number.
   void setSeed(unsigned int randomseed)
   {
      seed = randomseed;
   }
   //@}


   //--------------------------------------
   /**@name Constructors / destructors */
   //@{
   /// default constructor.
   /** Constructs a new (pseudo) Random variable using \p randomseed as seed for the random
       variable's sequence.
   */
   explicit
   Random(unsigned int randomseed = 0)
   {
      setSeed(randomseed);
   }

   /// destructor
   ~Random()
   {}
   //@}
};

} // namespace soplex
#endif // _RANDOM_H_
