/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the class library                   */
/*       SoPlex --- the Sequential object-oriented simPlex.                  */
/*                                                                           */
/*  Copyright (c) 1996-2025 Zuse Institute Berlin (ZIB)                      */
/*                                                                           */
/*  Licensed under the Apache License, Version 2.0 (the "License");          */
/*  you may not use this file except in compliance with the License.         */
/*  You may obtain a copy of the License at                                  */
/*                                                                           */
/*      http://www.apache.org/licenses/LICENSE-2.0                           */
/*                                                                           */
/*  Unless required by applicable law or agreed to in writing, software      */
/*  distributed under the License is distributed on an "AS IS" BASIS,        */
/*  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied. */
/*  See the License for the specific language governing permissions and      */
/*  limitations under the License.                                           */
/*                                                                           */
/*  You should have received a copy of the Apache-2.0 license                */
/*  along with SoPlex; see the file LICENSE. If not email to soplex@zib.de.  */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/**@file  random.h
 * @brief Random numbers.
 */
#ifndef _RANDOM_H_
#define _RANDOM_H_

#include <assert.h>
#include <stdint.h>
#include <algorithm>

namespace soplex
{

/* initial seeds for KISS random number generator */
#define SOPLEX_DEFAULT_LIN  UINT32_C(123456789)
#define SOPLEX_DEFAULT_XOR  UINT32_C(362436000)
#define SOPLEX_DEFAULT_MWC  UINT32_C(521288629)
#define SOPLEX_DEFAULT_CST  UINT32_C(7654321)

/* defines for linear congruential generator */
#define SOPLEX_RSTEP    UINT64_C(1103515245)
#define SOPLEX_RADD     UINT64_C(12345)

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
   ///@{
   uint32_t seedshift;  ///< initial shift for random seeds.
   uint32_t lin_seed;   ///< random seed for linear congruential RNS.
   uint32_t xor_seed;   ///< random seed for XOR-shift RNS.
   uint32_t mwc_seed;   ///< random seed Multiple-with-carry RNS.
   uint32_t cst_seed;   ///< random seed shifted for mwc_seed.
   ///@}

   //--------------------------------------
   /**@name Helpers */
   ///@{
   /// executes KISS random number generator and returns a pseudo random Real value in [0,1].
   Real next_random()
   {
      uint64_t t;

      /* linear congruential */
      lin_seed = (uint32_t)(lin_seed * SOPLEX_RSTEP + SOPLEX_RADD);

      /* Xorshift */
      xor_seed ^= (xor_seed << 13);
      xor_seed ^= (xor_seed >> 17);
      xor_seed ^= (xor_seed << 5);

      /* Multiply-with-carry */
      t = UINT64_C(698769069) * mwc_seed + cst_seed;
      cst_seed = (uint32_t)(t >> 32);
      mwc_seed = (uint32_t) t;

      return (lin_seed + xor_seed + mwc_seed) / (Real)UINT32_MAX;
   }

   ///@}

public:

   //--------------------------------------
   /**@name Access */
   ///@{
   /// returns next random number.
   Real next(Real minimum = 0.0, Real maximum = 1.0)
   {
      Real randnumber = next_random();

      /* we multiply minimum and maximum separately by randnumber in order to avoid overflow if they are more than
       * std::numeric_limits<Real>::max() apart
       */
      return minimum * (1.0 - randnumber) + maximum * randnumber;
   }

   /// returns the initial seed shift
   uint32_t getSeed() const
   {
      return seedshift;
   }

   ///@}

   //--------------------------------------
   /**@name Modification */
   ///@{
   /// initialize all seeds of the random number generator.
   void setSeed(uint32_t initshift)
   {
      seedshift = initshift;

      /* use SOPLEX_MAX to avoid zero after over flowing */
      lin_seed = SOPLEX_MAX(SOPLEX_DEFAULT_LIN + initshift, 1u);
      xor_seed = SOPLEX_MAX(SOPLEX_DEFAULT_XOR + initshift, 1u);
      mwc_seed = SOPLEX_MAX(SOPLEX_DEFAULT_MWC + initshift, 1u);
      cst_seed = SOPLEX_DEFAULT_CST + initshift;

      assert(lin_seed > 0);
      assert(xor_seed > 0);
      assert(mwc_seed > 0);

      /* advance state once to have more random values */
      (void) next_random();
   }

   ///@}


   //--------------------------------------
   /**@name Constructors / destructors */
   ///@{
   /// default constructor.
   /** Constructs a new (pseudo) Random variable using \p randomseed as seed for the random
       variable's sequence.
   */
   explicit
   Random(uint32_t randomseed = 0)
   {
      setSeed(randomseed);
   }

   /// destructor
   ~Random()
   {}
   ///@}
};

} // namespace soplex
#endif // _RANDOM_H_
