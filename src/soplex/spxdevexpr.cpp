/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the class library                   */
/*       SoPlex --- the Sequential object-oriented simPlex.                  */
/*                                                                           */
/*    Copyright (C) 1996-2018 Konrad-Zuse-Zentrum                            */
/*                            fuer Informationstechnik Berlin                */
/*                                                                           */
/*  SoPlex is distributed under the terms of the ZIB Academic Licence.       */
/*                                                                           */
/*  You should have received a copy of the ZIB Academic License              */
/*  along with SoPlex; see the file COPYING. If not email to soplex@zib.de.  */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

#include "soplex/spxdefines.h"
#include "soplex/spxdevexpr.h"

#define DEVEX_REFINETOL 2.0

namespace soplex
{
  // Definition of signature to avoid the specialization after instantiation error

  template <>
  void SPxDevexPR<Real>::setRep(typename SPxSolverBase<Real>::Representation);

  template <>
  int SPxDevexPR<Real>::selectLeaveX(Real feastol, int start, int incr);

  template <>
  bool SPxDevexPR<Real>::isConsistent() const;

  template <>
  int SPxDevexPR<Real>::selectLeaveSparse(Real feastol);

  template <>
  int SPxDevexPR<Real>::selectLeaveHyper(Real feastol);

  template <>
  SPxId SPxDevexPR<Real>::selectEnterX(Real tol);

  template <>
  void SPxDevexPR<Real>::addedVecs (int n);

  template <>
  void SPxDevexPR<Real>::addedCoVecs(int n);

  template <>
  SPxId SPxDevexPR<Real>::selectEnterHyperDim(Real& best, Real feastol);

  template <>
  SPxId SPxDevexPR<Real>::selectEnterHyperCoDim(Real& best, Real feastol);

  template <>
  SPxId SPxDevexPR<Real>::selectEnterSparseDim(Real& best, Real feastol);

  template <>
  SPxId SPxDevexPR<Real>::selectEnterSparseCoDim(Real& best, Real feastol);

  template <>
  SPxId SPxDevexPR<Real>::selectEnterDenseDim(Real& best, Real feastol, int start, int incr);

  template <>
  SPxId SPxDevexPR<Real>::selectEnterDenseCoDim(Real& best, Real feastol, int start, int incr);

  template <>
  void SPxDevexPR<Real>::load(SPxSolverBase<Real>* base)
  {
    this->thesolver = base;
    setRep(base->rep());
    assert(isConsistent());
  }

  template <>
  bool SPxDevexPR<Real>::isConsistent() const
  {
#ifdef ENABLE_CONSISTENCY_CHECKS
    if (this->thesolver != 0)
      if (weights.dim() != this->thesolver->coDim()
          || coWeights.dim() != this->thesolver->dim())
        return MSGinconsistent("SPxDevexPR");
#endif

    return true;
  }

  template <>
  void SPxDevexPR<Real>::setupWeights(typename SPxSolverBase<Real>::Type tp)
  {
    int i;
    int coWeightSize = 0;
    int weightSize = 0;

    DVector& weights = this->thesolver->weights;
    DVector& coWeights = this->thesolver->coWeights;

    if( tp == SPxSolverBase<Real>::ENTER )
      {
        coWeights.reDim(this->thesolver->dim(), false);
        for( i = this->thesolver->dim() - 1; i >= coWeightSize; --i )
          coWeights[i] = 2.0;

        weights.reDim(this->thesolver->coDim(), false);
        for( i = this->thesolver->coDim() - 1; i >= weightSize; --i )
          weights[i] = 2.0;
      }
    else
      {
        coWeights.reDim(this->thesolver->dim(), false);
        for( i = this->thesolver->dim() - 1; i >= coWeightSize; --i )
          coWeights[i] = 1.0;
      }
    this->thesolver->weightsAreSetup = true;
  }

  template <>
  void SPxDevexPR<Real>::setType(typename SPxSolverBase<Real>::Type tp)
  {
    setupWeights(tp);
    refined = false;

    bestPrices.clear();
    bestPrices.setMax(this->thesolver->dim());
    prices.reMax(this->thesolver->dim());

    if( tp == SPxSolverBase<Real>::ENTER )
      {
        bestPricesCo.clear();
        bestPricesCo.setMax(this->thesolver->coDim());
        pricesCo.reMax(this->thesolver->coDim());
      }

    assert(isConsistent());
  }

  /**@todo suspicious: Shouldn't the relation between dim, coDim, Vecs,
   *       and CoVecs be influenced by the representation ?
   */
  template <>
  void SPxDevexPR<Real>::setRep(typename SPxSolverBase<Real>::Representation)
  {
    if (this->thesolver != 0)
      {
        // resize weights and initialize new entries
        addedVecs(this->thesolver->coDim());
        addedCoVecs(this->thesolver->dim());
        assert(isConsistent());
      }
  }

  Real inline computePrice(Real viol, Real weight, Real tol)
  {
    if( weight < tol )
      return viol * viol / tol;
    else
      return viol * viol / weight;
  }

  template <>
  int SPxDevexPR<Real>::buildBestPriceVectorLeave( Real feastol )
  {
    int idx;
    int nsorted;
    Real fTesti;
    const Real* fTest = this->thesolver->fTest().get_const_ptr();
    const Real* cpen = this->thesolver->coWeights.get_const_ptr();
    typename SPxPricer<Real>::IdxElement price;
    prices.clear();
    bestPrices.clear();

    // TODO we should check infeasiblities for duplicates or loop over dimension
    //      bestPrices may then also contain duplicates!
    // construct vector of all prices
    for (int i = this->thesolver->infeasibilities.size() - 1; i >= 0; --i)
      {
        idx = this->thesolver->infeasibilities.index(i);
        fTesti = fTest[idx];
        if (fTesti < -feastol)
          {
            this->thesolver->isInfeasible[idx] = this->VIOLATED;
            price.idx = idx;
            price.val = computePrice(fTesti, cpen[idx], feastol);
            prices.append(price);
          }
      }
    // set up structures for the quicksort implementation
    this->compare.elements = prices.get_const_ptr();
    // do a partial sort to move the best ones to the front
    // TODO this can be done more efficiently, since we only need the indices
    nsorted = SPxQuicksortPart(prices.get_ptr(), this->compare, 0, prices.size(), HYPERPRICINGSIZE);
    // copy indices of best values to bestPrices
    for( int i = 0; i < nsorted; ++i )
      {
        bestPrices.addIdx(prices[i].idx);
        this->thesolver->isInfeasible[prices[i].idx] = this->VIOLATED_AND_CHECKED;
      }

    if( nsorted > 0 )
      return prices[0].idx;
    else
      return -1;
  }

  template <>
  int SPxDevexPR<Real>::selectLeave()
  {
    int retid;

    if (this->thesolver->hyperPricingLeave && this->thesolver->sparsePricingLeave)
      {
        if ( bestPrices.size() < 2 || this->thesolver->basis().lastUpdate() == 0 )
          {
            // call init method to build up price-vector and return index of largest price
            retid = buildBestPriceVectorLeave(this->theeps);
          }
        else
          retid = selectLeaveHyper(this->theeps);
      }
    else if (this->thesolver->sparsePricingLeave)
      retid = selectLeaveSparse(this->theeps);
    else
      retid = selectLeaveX(this->theeps);

    if ( retid < 0 && !refined )
      {
        refined = true;
        MSG_INFO3( (*this->thesolver->spxout), (*this->thesolver->spxout) << "WDEVEX02 trying refinement step..\n"; )
          retid = selectLeaveX(this->theeps/DEVEX_REFINETOL);
      }

    assert(retid < this->thesolver->dim());

    return retid;
  }

  template <>
  int SPxDevexPR<Real>::selectLeaveX(Real feastol, int start, int incr)
  {
    Real x;

    const Real* fTest = this->thesolver->fTest().get_const_ptr();
    const Real* cpen = this->thesolver->coWeights.get_const_ptr();
    Real best = 0;
    int bstI = -1;
    int end = this->thesolver->coWeights.dim();

    for (; start < end; start += incr)
      {
        if (fTest[start] < -feastol)
          {
            x = computePrice(fTest[start], cpen[start], feastol);
            if (x > best)
              {
                best = x;
                bstI = start;
                last = cpen[start];
              }
          }
      }
    return bstI;
  }

  template <>
  int SPxDevexPR<Real>::selectLeaveSparse(Real feastol)
  {
    Real x;

    const Real* fTest = this->thesolver->fTest().get_const_ptr();
    const Real* cpen = this->thesolver->coWeights.get_const_ptr();
    Real best = 0;
    int bstI = -1;
    int idx = -1;

    for (int i = this->thesolver->infeasibilities.size() - 1; i >= 0; --i)
      {
        idx = this->thesolver->infeasibilities.index(i);
        x = fTest[idx];
        if (x < -feastol)
          {
            x = computePrice(x, cpen[idx], feastol);
            if (x > best)
              {
                best = x;
                bstI = idx;
                last = cpen[idx];
              }
          }
        else
          {
            this->thesolver->infeasibilities.remove(i);
            assert(this->thesolver->isInfeasible[idx] == VIOLATED || this->thesolver->isInfeasible[idx] == VIOLATED_AND_CHECKED);
            this->thesolver->isInfeasible[idx] = this->NOT_VIOLATED;
          }
      }
    return bstI;
  }

  template <>
  int SPxDevexPR<Real>::selectLeaveHyper(Real feastol)
  {
    Real x;

    const Real* fTest = this->thesolver->fTest().get_const_ptr();
    const Real* cpen = this->thesolver->coWeights.get_const_ptr();
    Real best = 0;
    Real leastBest = infinity;
    int bstI = -1;
    int idx = -1;

    // find the best price from the short candidate list
    for( int i = bestPrices.size() - 1; i >= 0; --i )
      {
        idx = bestPrices.index(i);
        x = fTest[idx];
        if( x < -feastol )
          {
            x = computePrice(x, cpen[idx], feastol);
            if( x > best )
              {
                best = x;
                bstI = idx;
                last = cpen[idx];
              }
            // get the smallest price of candidate list
            if( x < leastBest )
              leastBest = x;
          }
        else
          {
            bestPrices.remove(i);
            this->thesolver->isInfeasible[idx] = this->NOT_VIOLATED;
          }
      }

    // make sure we do not skip potential candidates due to a high leastBest value
    if( leastBest == R(infinity) )
      {
        assert(bestPrices.size() == 0);
        leastBest = 0;
      }

    // scan the updated indices for a better price
    for( int i = this->thesolver->updateViols.size() - 1; i >= 0; --i )
      {
        idx = this->thesolver->updateViols.index(i);
        // only look at indeces that were not checked already
        if( this->thesolver->isInfeasible[idx] == this->VIOLATED )
          {
            x = fTest[idx];
            assert(x < -feastol);
            x = computePrice(x, cpen[idx], feastol);
            if( x > leastBest )
              {
                if( x > best )
                  {
                    best = x;
                    bstI = idx;
                    last = cpen[idx];
                  }
                // put index into candidate list
                this->thesolver->isInfeasible[idx] = this->VIOLATED_AND_CHECKED;
                bestPrices.addIdx(idx);
              }
          }
      }

    return bstI;
  }

  template <>
  void SPxDevexPR<Real>::left4(int n, SPxId id)
  {
    DVector& coWeights = this->thesolver->coWeights;
    if (id.isValid())
      {
        int i, j;
        Real x;
        const Real* rhoVec = this->thesolver->fVec().delta().values();
        Real rhov_1 = 1 / rhoVec[n];
        Real beta_q = this->thesolver->coPvec().delta().length2() * rhov_1 * rhov_1;

#ifndef NDEBUG
        if (spxAbs(rhoVec[n]) < theeps)
          {
            MSG_INFO3( (*this->thesolver->spxout), (*this->thesolver->spxout) << "WDEVEX01: rhoVec = "
                       << rhoVec[n] << " with smaller absolute value than theeps = " << theeps << std::endl; )
              }
#endif  // NDEBUG

        //  Update #coPenalty# vector
        const IdxSet& rhoIdx = this->thesolver->fVec().idx();
        int len = this->thesolver->fVec().idx().size();
        for (i = len - 1; i >= 0; --i)
          {
            j = rhoIdx.index(i);
            x = rhoVec[j] * rhoVec[j] * beta_q;
            // if(x > coPenalty[j])
            coWeights[j] += x;
          }

        coWeights[n] = beta_q;
      }
  }

  template <>
  SPxId SPxDevexPR<Real>::buildBestPriceVectorEnterDim( Real& best, Real feastol )
  {
    int idx;
    int nsorted;
    Real x;
    const Real* coTest = this->thesolver->coTest().get_const_ptr();
    const Real* cpen = this->thesolver->coWeights.get_const_ptr();
    typename SPxPricer<Real>::IdxElement price;
    prices.clear();
    bestPrices.clear();

    // construct vector of all prices
    for (int i = this->thesolver->infeasibilities.size() - 1; i >= 0; --i)
      {
        idx = this->thesolver->infeasibilities.index(i);
        x = coTest[idx];
        if ( x < -feastol)
          {
            this->thesolver->isInfeasible[idx] = this->VIOLATED;
            price.idx = idx;
            price.val = computePrice(x, cpen[idx], feastol);
            prices.append(price);
          }
        else
          {
            this->thesolver->infeasibilities.remove(i);
            this->thesolver->isInfeasible[idx] = this->NOT_VIOLATED;
          }
      }
    // set up structures for the quicksort implementation
    this->compare.elements = prices.get_const_ptr();
    // do a partial sort to move the best ones to the front
    // TODO this can be done more efficiently, since we only need the indices
    nsorted = SPxQuicksortPart(prices.get_ptr(), this->compare, 0, prices.size(), HYPERPRICINGSIZE);
    // copy indices of best values to bestPrices
    for( int i = 0; i < nsorted; ++i )
      {
        bestPrices.addIdx(prices[i].idx);
        this->thesolver->isInfeasible[prices[i].idx] = this->VIOLATED_AND_CHECKED;
      }

    if( nsorted > 0 )
      {
        best = prices[0].val;
        return this->thesolver->coId(prices[0].idx);
      }
    else
      return SPxId();
  }

  template <>
  SPxId SPxDevexPR<Real>::buildBestPriceVectorEnterCoDim( Real& best, Real feastol )
  {
    int idx;
    int nsorted;
    Real x;
    const Real* test = this->thesolver->test().get_const_ptr();
    const Real* pen = this->thesolver->weights.get_const_ptr();
    typename SPxPricer<Real>::IdxElement price;
    pricesCo.clear();
    bestPricesCo.clear();

    // construct vector of all prices
    for (int i = this->thesolver->infeasibilitiesCo.size() - 1; i >= 0; --i)
      {
        idx = this->thesolver->infeasibilitiesCo.index(i);
        x = test[idx];
        if ( x < -feastol)
          {
            this->thesolver->isInfeasibleCo[idx] = this->VIOLATED;
            price.idx = idx;
            price.val = computePrice(x, pen[idx], feastol);
            pricesCo.append(price);
          }
        else
          {
            this->thesolver->infeasibilitiesCo.remove(i);
            this->thesolver->isInfeasibleCo[idx] = this->NOT_VIOLATED;
          }
      }
    // set up structures for the quicksort implementation
    this->compare.elements = pricesCo.get_const_ptr();
    // do a partial sort to move the best ones to the front
    // TODO this can be done more efficiently, since we only need the indices
    nsorted = SPxQuicksortPart(pricesCo.get_ptr(), this->compare, 0, pricesCo.size(), HYPERPRICINGSIZE);
    // copy indices of best values to bestPrices
    for( int i = 0; i < nsorted; ++i )
      {
        bestPricesCo.addIdx(pricesCo[i].idx);
        this->thesolver->isInfeasibleCo[pricesCo[i].idx] = this->VIOLATED_AND_CHECKED;
      }

    if( nsorted > 0 )
      {
        best = pricesCo[0].val;
        return this->thesolver->id(pricesCo[0].idx);
      }
    else
      return SPxId();
  }

  template <>
  SPxId SPxDevexPR<Real>::selectEnter()
  {
    assert(this->thesolver != 0);

    SPxId enterId;

    enterId = selectEnterX(this->theeps);

    if( !enterId.isValid() && !refined )
      {
        refined = true;
        MSG_INFO3( (*this->thesolver->spxout), (*this->thesolver->spxout) << "WDEVEX02 trying refinement step..\n"; )
          enterId = selectEnterX(this->theeps/DEVEX_REFINETOL);
      }

    return enterId;
  }

  // choose the best entering index among columns and rows but prefer sparsity
  template <>
  SPxId SPxDevexPR<Real>::selectEnterX(Real tol)
  {
    SPxId enterId;
    SPxId enterCoId;
    Real best;
    Real bestCo;

    best = 0;
    bestCo = 0;
    last = 1.0;

    // avoid uninitialized value later on in entered4X()
    last = 1.0;

    if( this->thesolver->hyperPricingEnter && !refined )
      {
        if( bestPrices.size() < 2 || this->thesolver->basis().lastUpdate() == 0 )
          enterCoId = (this->thesolver->sparsePricingEnter) ? buildBestPriceVectorEnterDim(best, tol) : selectEnterDenseDim(best, tol);
        else
          enterCoId = (this->thesolver->sparsePricingEnter) ? selectEnterHyperDim(best, tol) : selectEnterDenseDim(best, tol);

        if( bestPricesCo.size() < 2 || this->thesolver->basis().lastUpdate() == 0 )
          enterId = (this->thesolver->sparsePricingEnterCo) ? buildBestPriceVectorEnterCoDim(bestCo, tol) : selectEnterDenseCoDim(bestCo, tol);
        else
          enterId = (this->thesolver->sparsePricingEnterCo) ? selectEnterHyperCoDim(bestCo, tol) : selectEnterDenseCoDim(bestCo, tol);
      }
    else
      {
        enterCoId = (this->thesolver->sparsePricingEnter && !refined) ? selectEnterSparseDim(best, tol) : selectEnterDenseDim(best, tol);
        enterId = (this->thesolver->sparsePricingEnterCo && !refined) ? selectEnterSparseCoDim(bestCo, tol) : selectEnterDenseCoDim(bestCo, tol);
      }

    // prefer coIds to increase the number of unit vectors in the basis matrix, i.e., rows in colrep and cols in rowrep
    if( enterCoId.isValid() && (best > SPARSITY_TRADEOFF * bestCo || !enterId.isValid()) )
      return enterCoId;
    else
      return enterId;
  }

  template <>
  SPxId SPxDevexPR<Real>::selectEnterHyperDim(Real& best, Real feastol)
  {
    const Real* cTest = this->thesolver->coTest().get_const_ptr();
    const Real* cpen = this->thesolver->coWeights.get_const_ptr();
    Real leastBest = R(infinity);
    Real x;
    int enterIdx = -1;
    int idx;

    // find the best price from short candidate list
    for( int i = bestPrices.size() - 1; i >= 0; --i )
      {
        idx = bestPrices.index(i);
        x = cTest[idx];
        if( x < -feastol )
          {
            x = computePrice(x, cpen[idx], feastol);
            if( x > best )
              {
                best = x;
                enterIdx = idx;
                last = cpen[idx];
              }
            if( x < leastBest )
              leastBest = x;
          }
        else
          {
            bestPrices.remove(i);
            this->thesolver->isInfeasible[idx] = this->NOT_VIOLATED;
          }
      }

    // make sure we do not skip potential candidates due to a high leastBest value
    if( leastBest == R(infinity) )
      {
        assert(bestPrices.size() == 0);
        leastBest = 0;
      }

    // scan the updated indeces for a better price
    for( int i = this->thesolver->updateViols.size() -1; i >= 0; --i )
      {
        idx = this->thesolver->updateViols.index(i);
        // only look at indeces that were not checked already
        if( this->thesolver->isInfeasible[idx] == this->VIOLATED )
          {
            x = cTest[idx];
            if( x < -feastol )
              {
                x = computePrice(x, cpen[idx], feastol);
                if(x > leastBest)
                  {
                    if( x > best )
                      {
                        best = x;
                        enterIdx = idx;
                        last = cpen[idx];
                      }
                    // put index into candidate list
                    this->thesolver->isInfeasible[idx] = this->VIOLATED_AND_CHECKED;
                    bestPrices.addIdx(idx);
                  }
              }
            else
              {
                this->thesolver->isInfeasible[idx] = this->NOT_VIOLATED;
              }
          }
      }

    if (enterIdx >= 0)
      return this->thesolver->coId(enterIdx);
    else
      return SPxId();
  }

  template <>
  SPxId SPxDevexPR<Real>::selectEnterHyperCoDim(Real& best, Real feastol)
  {
    const Real* test = this->thesolver->test().get_const_ptr();
    const Real* pen = this->thesolver->weights.get_const_ptr();
    Real leastBest = R(infinity);
    Real x;
    int enterIdx = -1;
    int idx;

    // find the best price from short candidate list
    for( int i = bestPricesCo.size() - 1; i >= 0; --i )
      {
        idx = bestPricesCo.index(i);
        x = test[idx];
        if( x < -feastol )
          {
            x = computePrice(x, pen[idx], feastol);
            if( x > best )
              {
                best = x;
                enterIdx = idx;
                last = pen[idx];
              }
            if( x < leastBest )
              leastBest = x;
          }
        else
          {
            bestPricesCo.remove(i);
            this->thesolver->isInfeasibleCo[idx] = this->NOT_VIOLATED;
          }
      }
    // make sure we do not skip potential candidates due to a high leastBest value
    if( leastBest == R(infinity) )
      {
        assert(bestPricesCo.size() == 0);
        leastBest = 0;
      }

    //scan the updated indeces for a better price
    for( int i = this->thesolver->updateViolsCo.size() -1; i >= 0; --i )
      {
        idx = this->thesolver->updateViolsCo.index(i);
        // only look at indeces that were not checked already
        if( this->thesolver->isInfeasibleCo[idx] == this->VIOLATED )
          {
            x = test[idx];
            if( x < -feastol )
              {
                x = computePrice(x, pen[idx], feastol);
                if( x > leastBest )
                  {
                    if( x > best )
                      {
                        best = x;
                        enterIdx = idx;
                        last = pen[idx];
                      }
                    // put index into candidate list
                    this->thesolver->isInfeasibleCo[idx] = this->VIOLATED_AND_CHECKED;
                    bestPricesCo.addIdx(idx);
                  }
              }
            else
              {
                this->thesolver->isInfeasibleCo[idx] = this->NOT_VIOLATED;
              }
          }
      }

    if( enterIdx >= 0 )
      return this->thesolver->id(enterIdx);
    else
      return SPxId();
  }


  template <>
  SPxId SPxDevexPR<Real>::selectEnterSparseDim(Real& best, Real feastol)
  {
    const Real* cTest = this->thesolver->coTest().get_const_ptr();
    const Real* cpen = this->thesolver->coWeights.get_const_ptr();
    int enterIdx = -1;
    int idx;
    Real x;

    assert(this->thesolver->coWeights.dim() == this->thesolver->coTest().dim());
    for(int i = this->thesolver->infeasibilities.size() -1; i >= 0; --i)
      {
        idx = this->thesolver->infeasibilities.index(i);
        x = cTest[idx];
        if( x < -feastol )
          {
            x = computePrice(x, cpen[idx], feastol);
            if (x > best)
              {
                best = x;
                enterIdx = idx;
                last = cpen[idx];
              }
          }
        else
          {
            this->thesolver->infeasibilities.remove(i);
            this->thesolver->isInfeasible[idx] = this->NOT_VIOLATED;
          }
      }
    if (enterIdx >= 0)
      return this->thesolver->coId(enterIdx);

    return SPxId();
  }


  template <>
  SPxId SPxDevexPR<Real>::selectEnterSparseCoDim(Real& best, Real feastol)
  {
    const Real* test = this->thesolver->test().get_const_ptr();
    const Real* pen = this->thesolver->weights.get_const_ptr();
    int enterIdx = -1;
    int idx;
    Real x;

    assert(this->thesolver->weights.dim() == this->thesolver->test().dim());
    for (int i = this->thesolver->infeasibilitiesCo.size() -1; i >= 0; --i)
      {
        idx = this->thesolver->infeasibilitiesCo.index(i);
        x = test[idx];
        if (x < -feastol)
          {
            x = computePrice(x, pen[idx], feastol);
            if (x > best)
              {
                best = x;
                enterIdx = idx;
                last = pen[idx];
              }
          }
        else
          {
            this->thesolver->infeasibilitiesCo.remove(i);
            this->thesolver->isInfeasibleCo[idx] = this->NOT_VIOLATED;
          }
      }

    if (enterIdx >= 0)
      return this->thesolver->id(enterIdx);

    return SPxId();
  }


  template <>
  SPxId SPxDevexPR<Real>::selectEnterDenseDim(Real& best, Real feastol, int start, int incr)
  {
    const Real* cTest = this->thesolver->coTest().get_const_ptr();
    const Real* cpen = this->thesolver->coWeights.get_const_ptr();
    int end = this->thesolver->coWeights.dim();
    int enterIdx = -1;
    Real x;

    assert(end == this->thesolver->coTest().dim());
    for (; start < end; start += incr)
      {
        x = cTest[start];
        if( x < -feastol )
          {
            x = computePrice(x, cpen[start], feastol);
            if (x > best)
              {
                best = x;
                enterIdx = start;
                last = cpen[start];
              }
          }
      }

    if (enterIdx >= 0)
      return this->thesolver->coId(enterIdx);

    return SPxId();
  }

  template <>
  SPxId SPxDevexPR<Real>::selectEnterDenseCoDim(Real& best, Real feastol, int start, int incr)
  {
    const Real* test = this->thesolver->test().get_const_ptr();
    const Real* pen = this->thesolver->weights.get_const_ptr();
    int end = this->thesolver->weights.dim();
    int enterIdx = -1;
    Real x;

    assert(end == this->thesolver->test().dim());
    for (; start < end; start += incr)
      {
        x = test[start];
        if (test[start] < -feastol)
          {
            x = computePrice(x, pen[start], feastol);
            if (x > best)
              {
                best = x;
                enterIdx = start;
                last = pen[start];
              }
          }
      }

    if (enterIdx >= 0)
      return this->thesolver->id(enterIdx);

    return SPxId();
  }


  /**@todo suspicious: the pricer should be informed, that variable id
     has entered the basis at position n, but the id is not used here
     (this is true for all pricers)
  */
  template <>
  void SPxDevexPR<Real>::entered4(SPxId /*id*/, int n)
  {
    DVector& weights = this->thesolver->weights;
    DVector& coWeights = this->thesolver->coWeights;

    if (n >= 0 && n < this->thesolver->dim())
      {
        const Real* pVec = this->thesolver->pVec().delta().values();
        const IdxSet& pIdx = this->thesolver->pVec().idx();
        const Real* coPvec = this->thesolver->coPvec().delta().values();
        const IdxSet& coPidx = this->thesolver->coPvec().idx();
        Real xi_p = 1 / this->thesolver->fVec().delta()[n];
        int i, j;

        assert(this->thesolver->fVec().delta()[n] > this->thesolver->epsilon()
               || this->thesolver->fVec().delta()[n] < -this->thesolver->epsilon());

        xi_p = xi_p * xi_p * last;

        for (j = coPidx.size() - 1; j >= 0; --j)
          {
            i = coPidx.index(j);
            coWeights[i] += xi_p * coPvec[i] * coPvec[i];
            if (coWeights[i] <= 1 || coWeights[i] > 1e+6)
              {
                setupWeights(SPxSolverBase<Real>::ENTER);
                return;
              }
          }

        for (j = pIdx.size() - 1; j >= 0; --j)
          {
            i = pIdx.index(j);
            weights[i] += xi_p * pVec[i] * pVec[i];
            if (weights[i] <= 1 || weights[i] > 1e+6)
              {
                setupWeights(SPxSolverBase<Real>::ENTER);
                return;
              }
          }
      }
  }

  template <>
  void SPxDevexPR<Real>::addedVecs (int n)
  {
    int initval = (this->thesolver->type() == SPxSolverBase<Real>::ENTER) ? 2 : 1;
    DVector& weights = this->thesolver->weights;
    n = weights.dim();
    weights.reDim (this->thesolver->coDim());
    for( int i = weights.dim()-1; i >= n; --i )
      weights[i] = initval;
  }

  template <>
  void SPxDevexPR<Real>::addedCoVecs(int n)
  {
    int initval = (this->thesolver->type() == SPxSolverBase<Real>::ENTER) ? 2 : 1;
    DVector& coWeights = this->thesolver->coWeights;
    n = coWeights.dim();
    coWeights.reDim(this->thesolver->dim());
    for( int i = coWeights.dim()-1; i >= n; --i )
      coWeights[i] = initval;
  }

} // namespace soplex
