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


/**@file  spxsteeppr.h
 * @brief Steepest edge pricer.
 */
#ifndef _SPXSTEEPPR_H_
#define _SPXSTEEPPR_H_


#include <assert.h>

#include "soplex/spxdefines.h"
#include "soplex/spxpricer.h"
#include "soplex/random.h"

namespace soplex
{

/**@brief   Steepest edge pricer.
   @ingroup Algo

   Class SPxSteepPR implements a steepest edge pricer to be used with
   SoPlex.

   See SPxPricer for a class documentation.
*/
template <class R>
class SPxSteepPR : public SPxPricer<R>
{
public:

   //-------------------------------------
   /**@name Types */
   ///@{
   /// How to setup the direction multipliers.
   /** Possible settings are #EXACT for starting with exactly computed
       values, or #DEFAULT for starting with multipliers set to 1. The
       latter is the default.
   */
   enum Setup
   {
      EXACT,   ///< starting with exactly computed values
      DEFAULT  ///< starting with multipliers set to 1
   };
   ///@}
   /// setup steepest edge weights
   void setupWeights(typename SPxSolverBase<R>::Type type);

private:

   //-------------------------------------
   /**@name Data */
   ///@{
   /// working vector
   SSVectorBase<R>  workVec;
   /// working vector
   SSVectorBase<R>  workRhs;
   /// temporary array of precomputed pricing values
   Array<typename SPxPricer<R>::IdxElement> prices;
   /// temporary array of precomputed pricing values
   Array<typename SPxPricer<R>::IdxElement> pricesCo;
   /// array of best pricing candidates
   DIdxSet bestPrices;
   /// array of best pricing candidates
   DIdxSet bestPricesCo;
   ///
   R pi_p;
   /// setup type.
   Setup setup;
   /// has a refinement step already been tried?
   bool refined;
   ///@}

   //-------------------------------------
   /// prepare data structures for hyper sparse pricing
   int buildBestPriceVectorLeave(R feastol);
   /// implementation of full pricing
   int selectLeaveX(R tol);
   /// implementation of sparse pricing in the leaving Simplex
   int selectLeaveSparse(R tol);
   /// implementation of hyper sparse pricing in the leaving Simplex
   int selectLeaveHyper(R tol);
   /// build up vector of pricing values for later use
   SPxId buildBestPriceVectorEnterDim(R& best, R feastol);
   SPxId buildBestPriceVectorEnterCoDim(R& best, R feastol);
   /// choose the best entering index among columns and rows but prefer sparsity
   SPxId selectEnterX(R tol);
   /// implementation of sparse pricing for the entering Simplex (slack variables)
   SPxId selectEnterSparseDim(R& best, R tol);
   /// implementation of sparse pricing for the entering Simplex
   SPxId selectEnterSparseCoDim(R& best, R tol);
   /// implementation of selectEnter() in dense case (slack variables)
   SPxId selectEnterDenseDim(R& best, R tol);
   /// implementation of selectEnter() in dense case
   SPxId selectEnterDenseCoDim(R& best, R tol);
   /// implementation of hyper sparse pricing in the entering Simplex
   SPxId selectEnterHyperDim(R& best, R feastol);
   /// implementation of hyper sparse pricing in the entering Simplex
   SPxId selectEnterHyperCoDim(R& best, R feastol);

public:

   //-------------------------------------
   /**@name Construction / destruction */
   ///@{
   ///
   SPxSteepPR(const char* name = "Steep", Setup mode = DEFAULT)
      : SPxPricer<R>(name)
      , workVec(0, nullptr)
      , workRhs(0, nullptr)
      , pi_p(1.0)
      , setup(mode)
      , refined(false)
   {
      assert(isConsistent());
   }
   /// copy constructor
   SPxSteepPR(const SPxSteepPR& old)
      : SPxPricer<R>(old)
      , workVec(old.workVec)
      , workRhs(old.workRhs)
      , pi_p(old.pi_p)
      , setup(old.setup)
      , refined(old.refined)
   {
      assert(isConsistent());
   }
   /// assignment operator
   SPxSteepPR& operator=(const SPxSteepPR& rhs)
   {
      if(this != &rhs)
      {
         SPxPricer<R>::operator=(rhs);
         workVec = rhs.workVec;
         workRhs = rhs.workRhs;
         pi_p = rhs.pi_p;
         setup = rhs.setup;
         refined = rhs.refined;

         assert(isConsistent());
      }

      return *this;
   }
   /// destructor
   virtual ~SPxSteepPR()
   {}
   /// clone function for polymorphism
   inline virtual SPxPricer<R>* clone()  const
   {
      return new SPxSteepPR(*this);
   }
   ///@}

   //-------------------------------------
   /**@name Access / modification */
   ///@{
   /// sets the solver
   virtual void load(SPxSolverBase<R>* base);
   /// clear solver and preferences
   virtual void clear();
   /// set entering/leaving algorithm
   virtual void setType(typename SPxSolverBase<R>::Type);
   /// set row/column representation
   virtual void setRep(typename SPxSolverBase<R>::Representation rep);
   ///
   virtual int selectLeave();
   ///
   virtual void left4(int n, SPxId id);
   ///
   virtual SPxId selectEnter();
   ///
   virtual void entered4(SPxId id, int n);
   /// \p n vectors have been added to loaded LP.
   virtual void addedVecs(int n);
   /// \p n covectors have been added to loaded LP.
   virtual void addedCoVecs(int n);
   /// \p the i'th vector has been removed from the loaded LP.
   virtual void removedVec(int i);
   /// \p the i'th covector has been removed from the loaded LP.
   virtual void removedCoVec(int i);
   /// \p n vectors have been removed from loaded LP.
   virtual void removedVecs(const int perm[]);
   /// \p n covectors have been removed from loaded LP.
   virtual void removedCoVecs(const int perm[]);
   ///@}

   //-------------------------------------
   /**@name Consistency check */
   ///@{
   ///
   virtual bool isConsistent() const;
   ///@}
};

} // namespace soplex
#endif // _SPXSTEEPPR_H_
