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

#include <assert.h>
#include <iostream>

// #define EQ_PREF 1000 

#include "spxdefines.h"
#include "spxdantzigpr.h"

namespace soplex
{

  template <class R>
  int SPxDantzigPR<R>::selectLeave()
  {
    assert(thesolver != 0);

    if( this->thesolver->sparsePricingLeave )
      return selectLeaveSparse();

    //    const Real* up  = this->thesolver->ubBound();
    //    const Real* low = this->thesolver->lbBound();

    Real best = -this->theeps;
    int  n    = -1;

    for(int i = this->thesolver->dim() - 1; i >= 0; --i)
      {
        Real x = this->thesolver->fTest()[i];

        if (x < -this->theeps)
          {
            // x *= EQ_PREF * (1 + (up[i] == low[i]));
            if (x < best)
              {
                n    = i;
                best = x;
              }
          }
      }
    return n;
  }

  template <class R>
  int SPxDantzigPR<R>::selectLeaveSparse()
  {
    assert(thesolver != 0);

    Real best   = -this->theeps;
    Real x;
    int  n      = -1;
    int  index;

    for(int i = this->thesolver->infeasibilities.size() - 1; i >= 0; --i)
      {
        index = this->thesolver->infeasibilities.index(i);
        x = this->thesolver->fTest()[index];
        if (x < -this->theeps)
          {
            if (x < best)
              {
                n    = index;
                best = x;
              }
          }
        else
          {
            this->thesolver->infeasibilities.remove(i);
            assert(this->thesolver->isInfeasible[index] > 0);
            this->thesolver->isInfeasible[index] = 0;
          }
      }
    return n;
  }

  template <class R>
  SPxId SPxDantzigPR<R>::selectEnter()
  {
    assert(thesolver != 0);

    // const SPxBasis<R>::Desc&    ds   = this->thesolver->basis().desc();

    SPxId enterId;
    enterId = selectEnterX();

    return enterId;
  }

  template <class R>
  SPxId SPxDantzigPR<R>::selectEnterX()
  {
    SPxId enterId;
    SPxId enterIdCo;
    Real best;
    Real bestCo;

    best = -this->theeps;
    bestCo = -this->theeps;
    enterId = (this->thesolver->sparsePricingEnter) ? selectEnterSparseDim(best,enterId) : selectEnterDenseDim(best,enterId);
    enterIdCo = (this->thesolver->sparsePricingEnterCo) ? selectEnterSparseCoDim(bestCo,enterId) : selectEnterDenseCoDim(bestCo,enterId);

    // prefer slack indices to reduce nonzeros in basis matrix
    if( enterId.isValid() && (best > SPARSITY_TRADEOFF * bestCo || !enterIdCo.isValid()) )
      return enterId;
    else
      return enterIdCo;
  }


  template <class R>
  SPxId SPxDantzigPR<R>::selectEnterSparseDim(Real& best,SPxId& enterId)
  {
    assert(thesolver != 0);

    int idx;
    Real x;

    for (int i = this->thesolver->infeasibilities.size() - 1; i >= 0; --i)
      {
        idx = this->thesolver->infeasibilities.index(i);
        x = this->thesolver->coTest()[idx];

        if (x < -this->theeps)
          {
            // x *= EQ_PREF * (1 + (ds.coStatus(i) == SPxBasis<R>::Desc::P_FREE
            //                || ds.coStatus(i) == SPxBasis<R>::Desc::D_FREE));
            if (x < best)
              {
                enterId = this->thesolver->coId(idx);
                best = x;
              }
          }
        else
          {
            this->thesolver->infeasibilities.remove(i);

            assert(this->thesolver->isInfeasible[idx]);
            this->thesolver->isInfeasible[idx] = 0;
          }
      }
    return enterId;
  }

  template <class R>
  SPxId SPxDantzigPR<R>::selectEnterSparseCoDim(Real& best,SPxId& enterId)
  {
    assert(thesolver != 0);

    int idx;
    Real x;

    for (int i = this->thesolver->infeasibilitiesCo.size() - 1; i >= 0; --i)
      {
        idx = this->thesolver->infeasibilitiesCo.index(i);
        x = this->thesolver->test()[idx];

        if (x < -this->theeps)
          {
            // x *= EQ_PREF * (1 + (ds.coStatus(i) == SPxBasis<R>::Desc::P_FREE
            //                || ds.coStatus(i) == SPxBasis<R>::Desc::D_FREE));
            if (x < best)
              {
                enterId = this->thesolver->id(idx);
                best = x;
              }
          }
        else
          {
            this->thesolver->infeasibilitiesCo.remove(i);
            assert(this->thesolver->isInfeasibleCo[idx] > 0);
            this->thesolver->isInfeasibleCo[idx] = 0;
          }
      }
    return enterId;
  }

  template <class R>
  SPxId SPxDantzigPR<R>::selectEnterDenseDim(Real& best,SPxId& enterId)
  {
    assert(thesolver != 0);

    Real x;

    for (int i = this->thesolver->dim() - 1; i >= 0; --i)
      {
        x = this->thesolver->coTest()[i];

        if (x < -this->theeps)
          {
            // x *= EQ_PREF * (1 + (ds.coStatus(i) == SPxBasis<R>::Desc::P_FREE
            //                || ds.coStatus(i) == SPxBasis<R>::Desc::D_FREE));
            if (x < best)
              {
                enterId   = this->thesolver->coId(i);
                best = x;
              }
          }
      }
    return enterId;
  }

  template <class R>
  SPxId SPxDantzigPR<R>::selectEnterDenseCoDim(Real& best,SPxId& enterId)
  {
    assert(thesolver != 0);

    Real x;

    for (int i = this->thesolver->coDim() - 1; i >= 0; --i)
      {
        x = this->thesolver->test()[i];

        if (x < -this->theeps)
          {
            // x *= EQ_PREF * (1 + (ds.status(i) == SPxBasis<R>::Desc::P_FREE
            //                || ds.status(i) == SPxBasis<R>::Desc::D_FREE));
            if (x < best)
              {
                enterId   = this->thesolver->id(i);
                best = x;
              }
          }
      }
    return enterId;
  }

} // namespace soplex
