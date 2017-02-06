/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the class library                   */
/*       SoPlex --- the Sequential object-oriented simPlex.                  */
/*                                                                           */
/*    Copyright (C) 1996-2017 Konrad-Zuse-Zentrum                            */
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

/* initial seeds for KISS random number generator */
#define  DEFAULT_LIN  123456789
#define  DEFAULT_XOR  362436000
#define  DEFAULT_MWC  521288629
#define  DEFAULT_CST  7654321

/* defines for linear congruential generator */
#define  RSTEP    1103515245
#define  RADD     12345

/* defines for shifting an unsigned integer into [0,1) */
#define  RDIVIDE  65536
#define  RMULT    32768

/**@brief   Random numbers.
   @ingroup Elementary

   Class Random provides random Real variables, i.e. a value variable that
   gives another value each time it is accessed. It may be used just like an
   ordinary Real by means of an overloaded cast operator Real()%.

   This is an implementation of KISS random number generator developed by George Marsaglia.
   KISS is combination of three different random number generators:
    - Linear congruential generator
    - Xorshift
    - Lag-1 Multiply-with-carry

   KISS has a period of 2^123 and passes all statistical test part of BigCrush-Test of TestU01 [1].

   [1] http://dl.acm.org/citation.cfm?doid=1268776.1268777
*/
class Random
{
private:

   //--------------------------------------
   /**@name Data */
   //@{
   unsigned long seedshift;  ///< initial shift for random seeds.
   unsigned long lin_seed;   ///< random seed for linear congruential RNS.
   unsigned long xor_seed;   ///< random seed for XOR-shift RNS.
   unsigned long mwc_seed;   ///< random seed Multiple-with-carry RNS.
   unsigned long cst_seed;   ///< random seed shifted for mwc_seed.
   //@}

   //--------------------------------------
   /**@name Helpers */
   //@{
   /// executes KISS random number generator and returns a pseudo random Real value in [0,1).
   Real next_random()
   {
      unsigned long long t;

      /* linear congruential */
      lin_seed = lin_seed * RSTEP + RADD;

      /* Xorshift */
      xor_seed ^= (xor_seed << 13);
      xor_seed ^= (xor_seed >> 17);
      xor_seed ^= (xor_seed << 5);

      /* Multiply-with-carry */
      t = 698769069ULL * mwc_seed + cst_seed;
      cst_seed = (unsigned long) (t >> 32);
      mwc_seed = (unsigned long) t;

      Real i = (Real) (((lin_seed + xor_seed + mwc_seed) / RDIVIDE) % RMULT);
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

   /// returns the initial seed shift
   unsigned long getSeed() const
   {
      return seedshift;
   }

   //@}

   //--------------------------------------
   /**@name Modification */
   //@{
   /// initialize all seeds of the random number generator.
   void setSeed(unsigned long initshift)
   {
      seedshift = initshift;
      lin_seed = (unsigned long)(DEFAULT_LIN + initshift);
      xor_seed = (unsigned long)(DEFAULT_XOR + initshift);
      mwc_seed = (unsigned long)(DEFAULT_MWC + initshift);
      cst_seed = (unsigned long)(DEFAULT_CST + initshift);
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
   Random(unsigned long randomseed = 0)
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
