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

#include "soplex/spxdefines.h"
#include "soplex/spxweightpr.h"
#include "soplex/exceptions.h"

namespace soplex
{
  template <>
  void SPxWeightPR<Real>::computeLeavePenalty(int start, int end);

  template <>
  bool SPxWeightPR<Real>::isConsistent() const;

  template <>
  void SPxWeightPR<Real>::setRep(typename SPxSolver<Real>::Representation rep)
  {
    if (rep == SPxSolver<Real>::ROW)
      {
        penalty = rPenalty.get_const_ptr();
        coPenalty = cPenalty.get_const_ptr();
      }
    else
      {
        penalty = cPenalty.get_const_ptr();
        coPenalty = rPenalty.get_const_ptr();
      }
  }

  template <>
  void SPxWeightPR<Real>::setType(typename SPxSolver<Real>::Type tp)
  {
    if (this->thesolver && tp == SPxSolver<Real>::LEAVE)
      {
        leavePenalty.reDim( this->thesolver->dim() );
        computeLeavePenalty( 0, this->thesolver->dim() );
      }
  }

  template <>
  void SPxWeightPR<Real>::computeLeavePenalty(int start, int end)
  {
    const SPxBasis<Real>& basis = this->solver()->basis();

    for (int i = start; i < end; ++i)
      {
        SPxId id = basis.baseId(i);
        if (id.type() == SPxId::ROW_ID)
          leavePenalty[i] = rPenalty[ this->thesolver->number(id) ];
        else
          leavePenalty[i] = cPenalty[ this->thesolver->number(id) ];
      }
  }

  template <>
  void SPxWeightPR<Real>::computeRP(int start, int end)
  {
    for (int i = start; i < end; ++i)
      {
        /**@todo TK04NOV98 here is a bug.
         *       this->solver()->rowVector(i).length() could be zero, so
         *       this->solver()->rowVector(i).length2() is also zero and we
         *       get an arithmetic exception.
         */
        assert(this->solver()->rowVector(i).length() > 0);

        rPenalty[i] = (this->solver()->rowVector(i) * this->solver()->maxObj()) * objlength
          / this->solver()->rowVector(i).length2();
        ASSERT_WARN( "WWGTPR01", rPenalty[i] > -1 - this->solver()->epsilon() );
      }
  }

  template <>
  void SPxWeightPR<Real>::computeCP(int start, int end)
  {
    for (int i = start; i < end; ++i)
      {
        cPenalty[i] = this->solver()->maxObj(i) * objlength;
        ASSERT_WARN( "WWGTPR02", cPenalty[i] > -1 - this->solver()->epsilon() );
      }
  }

  template <>
  void SPxWeightPR<Real>::load(SPxSolver<Real>* base)
  {
    this->thesolver = base;

    rPenalty.reDim(base->nRows());
    cPenalty.reDim(base->nCols());

    objlength = 1 / this->solver()->maxObj().length();
    computeCP(0, base->nCols());
    computeRP(0, base->nRows());
  }

  template <>
  int SPxWeightPR<Real>::selectLeave()
  {
    const Real* test = this->thesolver->fTest().get_const_ptr();
    Real type = 1 - 2 * (this->thesolver->rep() == SPxSolver<Real>::COLUMN ? 1 : 0);
    Real best = type * infinity;
    int lastIdx = -1;
    Real x;
    int i;

    for (i = this->solver()->dim() - 1; i >= 0; --i)
      {
        x = test[i];
        if (x < -this->theeps)
          {
            x *= leavePenalty[i];
            if (type * (x-best) < 0.0)
              {
                best = x;
                lastIdx = i;
              }
          }
      }
    assert(isConsistent());
    return lastIdx;
  }

#if 0
  /**@todo remove this code */
  // ??? This is the old (buggy) version
  template <>
  int SPxWeightPR<Real>::selectLeave()
  {
    const Real* test = this->thesolver->fTest().get_const_ptr();
    Real type = 1 - 2 * (this->thesolver->rep() == SPxSolver<Real>::COLUMN);
    Real best = type * infinity;
    int lastIdx = -1;
    Real x;
    int i;

    for (i = this->solver()->dim() - 1; i >= 0; --i)
      {
        x = test[i];
        if (x < -theeps)
          {
            x *= leavePenalty[i];
            if (x < best)
              {
                best = x;
                lastIdx = i;
              }
          }
      }
    assert(isConsistent());
    return lastIdx;
  }
#endif

  template <>
  SPxId SPxWeightPR<Real>::selectEnter()
  {
    const Vector& rTest = (this->solver()->rep() == SPxSolver<Real>::ROW)
      ? this->solver()->test() : this->solver()->coTest();
    const Vector& cTest = (this->solver()->rep() == SPxSolver<Real>::ROW)
      ? this->solver()->coTest() : this->solver()->test();
    const typename SPxBasis<Real>::Desc& ds = this->solver()->basis().desc();
    Real best = infinity;
    SPxId lastId;
    Real x;
    int i;

    for (i = this->solver()->nRows() - 1; i >= 0; --i)
      {
        x = rTest[i];
        if (x < -this->theeps)
          {
            x *= -x;
            switch (ds.rowStatus(i))
              {
              case SPxBasis<Real>::Desc::P_ON_LOWER :
              case SPxBasis<Real>::Desc::D_ON_LOWER :
                x *= 1 + rPenalty[i];
                break;
              case SPxBasis<Real>::Desc::P_ON_UPPER :
              case SPxBasis<Real>::Desc::D_ON_UPPER :
                x *= 1 - rPenalty[i];
                break;
              case SPxBasis<Real>::Desc::P_FREE :
              case SPxBasis<Real>::Desc::D_FREE :
                return SPxId(this->solver()->rId(i));
              case SPxBasis<Real>::Desc::D_ON_BOTH :
                if (this->solver()->pVec()[i] > this->solver()->upBound()[i])
                  x *= 1 + rPenalty[i];
                else
                  x *= 1 - rPenalty[i];
                break;
              case SPxBasis<Real>::Desc::D_UNDEFINED :
              case SPxBasis<Real>::Desc::P_FIXED :
              default:
                throw SPxInternalCodeException("XWGTPR01 This should never happen.");
              }
            if (x < best)
              {
                best = x;
                lastId = this->solver()->rId(i);
              }
          }
      }

    for (i = this->solver()->nCols() - 1; i >= 0; --i)
      {
        x = cTest[i];
        if (x < -this->theeps)
          {
            x *= -x;
            switch (ds.colStatus(i))
              {
              case SPxBasis<Real>::Desc::P_ON_LOWER :
              case SPxBasis<Real>::Desc::D_ON_LOWER :
                x *= 1 + cPenalty[i];
                break;
              case SPxBasis<Real>::Desc::P_ON_UPPER :
              case SPxBasis<Real>::Desc::D_ON_UPPER :
                x *= 1 - cPenalty[i];
                break;
              case SPxBasis<Real>::Desc::P_FREE :
              case SPxBasis<Real>::Desc::D_FREE :
                return SPxId(this->solver()->cId(i));
              case SPxBasis<Real>::Desc::D_ON_BOTH :
                if (this->solver()->coPvec()[i] > this->solver()->ucBound()[i])
                  x *= 1 + cPenalty[i];
                else
                  x *= 1 - cPenalty[i];
                break;
              case SPxBasis<Real>::Desc::P_FIXED :
              case SPxBasis<Real>::Desc::D_UNDEFINED :
              default:
                throw SPxInternalCodeException("XWGTPR02 This should never happen.");
              }
            if (x < best)
              {
                best = x;
                lastId = this->solver()->cId(i);
              }
          }
      }
    assert(isConsistent());
    return lastId;
  }

  template <>
  void SPxWeightPR<Real>::addedVecs(int)
  {
    if (this->solver()->rep() == SPxSolver<Real>::ROW)
      {
        int start = rPenalty.dim();
        rPenalty.reDim(this->solver()->nRows());
        computeRP(start, this->solver()->nRows());
      }
    else
      {
        int start = cPenalty.dim();
        cPenalty.reDim(this->solver()->nCols());
        computeCP(start, this->solver()->nCols());
      }
    if (this->solver()->type() == SPxSolver<Real>::LEAVE)
      {
        int start = leavePenalty.dim();
        leavePenalty.reDim( this->solver()->dim() );
        computeLeavePenalty( start, this->solver()->dim() );
      }
  }

  template <>
  void SPxWeightPR<Real>::addedCoVecs(int)
  {
    if (this->solver()->rep() == SPxSolver<Real>::COLUMN)
      {
        int start = rPenalty.dim();
        rPenalty.reDim(this->solver()->nRows());
        computeRP(start, this->solver()->nRows());
      }
    else
      {
        int start = cPenalty.dim();
        cPenalty.reDim(this->solver()->nCols());
        computeCP(start, this->solver()->nCols());
      }
    if (this->solver()->type() == SPxSolver<Real>::LEAVE)
      {
        int start = leavePenalty.dim();
        leavePenalty.reDim( this->solver()->dim() );
        computeLeavePenalty( start, this->solver()->dim() );
      }
  }

  template <>
  void SPxWeightPR<Real>::removedVec(int i)
  {
    assert(this->solver() != 0);

    if (this->solver()->rep() == SPxSolver<Real>::ROW)
      {
        rPenalty[i] = rPenalty[rPenalty.dim()];
        rPenalty.reDim(this->solver()->nRows());
      }
    else
      {
        cPenalty[i] = cPenalty[cPenalty.dim()];
        cPenalty.reDim(this->solver()->nCols());
      }
  }

  template <>
  void SPxWeightPR<Real>::removedVecs(const int perm[])
  {
    assert(this->solver() != 0);

    if (this->solver()->rep() == SPxSolver<Real>::ROW)
      {
        int j = rPenalty.dim();
        for (int i = 0; i < j; ++i)
          {
            if (perm[i] >= 0)
              rPenalty[perm[i]] = rPenalty[i];
          }
        rPenalty.reDim(this->solver()->nRows());
      }
    else
      {
        int j = cPenalty.dim();
        for (int i = 0; i < j; ++i)
          {
            if (perm[i] >= 0)
              cPenalty[perm[i]] = cPenalty[i];
          }
        cPenalty.reDim(this->solver()->nCols());
      }
  }

  template <>
  void SPxWeightPR<Real>::removedCoVec(int i)
  {
    assert(this->solver() != 0);

    if (this->solver()->rep() == SPxSolver<Real>::COLUMN)
      {
        rPenalty[i] = rPenalty[rPenalty.dim()];
        rPenalty.reDim(this->solver()->nRows());
      }
    else
      {
        cPenalty[i] = cPenalty[cPenalty.dim()];
        cPenalty.reDim(this->solver()->nCols());
      }
  }

  template <>
  void SPxWeightPR<Real>::removedCoVecs(const int perm[])
  {
    assert(this->solver() != 0);

    if (this->solver()->rep() == SPxSolver<Real>::COLUMN)
      {
        int j = rPenalty.dim();
        for (int i = 0; i < j; ++i)
          {
            if (perm[i] >= 0)
              rPenalty[perm[i]] = rPenalty[i];
          }
        rPenalty.reDim(this->solver()->nRows());
      }
    else
      {
        int j = cPenalty.dim();
        for (int i = 0; i < j; ++i)
          {
            if (perm[i] >= 0)
              cPenalty[perm[i]] = cPenalty[i];
          }
        cPenalty.reDim(this->solver()->nCols());
      }
  }

  template <>
  bool SPxWeightPR<Real>::isConsistent() const
  {
#ifdef ENABLE_CONSISTENCY_CHECKS
    if (this->solver() != 0)
      {
        if (rPenalty.dim() != this->solver()->nRows())
          return MSGinconsistent("SPxWeightPR");
        if (cPenalty.dim() != this->solver()->nCols())
          return MSGinconsistent("SPxWeightPR");
      }
#endif

    return true;
  }
} // namespace soplex
