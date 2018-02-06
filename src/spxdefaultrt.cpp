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

#include "spxdefines.h"
#include "spxdefaultrt.h"

namespace soplex
{
  /**
   * Here comes the ratio test for selecting a variable to leave the basis. 
   * It is assumed that Vec.delta() and fVec.idx() have been setup
   * correctly!
   *
   * The leaving variable is selected such that the update of fVec() (using
   * fVec.value() * fVec.delta()) keeps the basis feasible within
   * solver()->entertol(). Hence, fVec.value() must be chosen such that one
   * updated value of theFvec just reaches its bound and no other one exceeds
   * them by more than solver()->entertol(). Further, fVec.value() must have the
   * same sign as argument \p val.
   *
   * The return value of selectLeave() is the number of a variable in the
   * basis selected to leave the basis. -1 indicates that no variable could be
   * selected. Otherwise, parameter \p val contains the chosen fVec.value().
   */
  template <class R>
  int SPxDefaultRT<R>::selectLeave(Real& val, Real, bool)
  {
    this->solver()->fVec().delta().setup();

    const Real*   vec = this->solver()->fVec().get_const_ptr();
    const Real*   upd = this->solver()->fVec().delta().values();
    const IdxSet& idx = this->solver()->fVec().idx();
    const Real*   ub  = this->solver()->ubBound().get_const_ptr();
    const Real*   lb  = this->solver()->lbBound().get_const_ptr();

    Real epsilon = this->solver()->epsilon();
    int  leave   = -1;

    Real x;
    int  i;
    int  j;

    // PARALLEL the j loop could be parallelized
    if (val > 0)
      {
        // Loop over NZEs of delta vector.
        for( j = 0; j < idx.size(); ++j)
          {
            i = idx.index(j);
            x = upd[i];
 
            if (x > epsilon)
              {
                if (ub[i] < infinity)
                  {
                    Real y = (ub[i] - vec[i] + this->delta) / x;

                    if (y < val)
                      {
                        leave = i;
                        val   = y;
                      }
                  }
              }
            else if (x < -epsilon)
              {
                if (lb[i] > -infinity)
                  {
                    Real y = (lb[i] - vec[i] - this->delta) / x;

                    if (y < val)
                      {
                        leave = i;
                        val   = y;
                      }
                  }
              }
          }
        if (leave >= 0)
          {
            x   = upd[leave];

            // BH 2005-11-30: It may well happen that the basis is degenerate and the
            // selected leaving variable is (at most this->delta) beyond its bound. (This
            // happens for instance on LP/netlib/adlittle.mps with setting -r -t0.) 
            // In that case we do a pivot step with length zero to avoid difficulties.
            if ( ( x > epsilon  && vec[leave] >= ub[leave] ) ||
                 ( x < -epsilon && vec[leave] <= lb[leave] ) )
              {
                val = 0.0;
              }
            else 
              {
                val = (x > epsilon) ? ub[leave] : lb[leave];
                val = (val - vec[leave]) / x;
              }
          }

        ASSERT_WARN( "WDEFRT01", val > -epsilon );
      }
    else
      {
        for( j = 0; j < idx.size(); ++j)
          {
            i = idx.index(j);
            x = upd[i];

            if (x < -epsilon)
              {
                if (ub[i] < infinity)
                  {
                    Real y = (ub[i] - vec[i] + this->delta) / x;

                    if (y > val)
                      {
                        leave = i;
                        val   = y;
                      }
                  }
              }
            else if (x > epsilon)
              {
                if (lb[i] > -infinity)
                  {
                    Real y = (lb[i] - vec[i] - this->delta) / x;

                    if (y > val)
                      {
                        leave = i;
                        val   = y;
                      }
                  }
              }
          }
        if (leave >= 0)
          {
            x   = upd[leave];

            // See comment above.
            if ( ( x < -epsilon && vec[leave] >= ub[leave] ) || 
                 ( x > epsilon  && vec[leave] <= lb[leave] ) )
              {
                val = 0.0;
              }
            else 
              {
                val = (x < epsilon) ? ub[leave] : lb[leave];
                val = (val - vec[leave]) / x;
              }
          }

        ASSERT_WARN( "WDEFRT02", val < epsilon );
      }
    return leave;
  }

  /**
     Here comes the ratio test. It is assumed that theCoPvec.this->delta() and
     theCoPvec.idx() have been setup correctly!
  */
  template <class R>
  SPxId SPxDefaultRT<R>::selectEnter(Real& max, int, bool)
  {
    this->solver()->coPvec().delta().setup();
    this->solver()->pVec().delta().setup();

    const Real*   pvec = this->solver()->pVec().get_const_ptr();
    const Real*   pupd = this->solver()->pVec().delta().values();
    const IdxSet& pidx = this->solver()->pVec().idx();
    const Real*   lpb  = this->solver()->lpBound().get_const_ptr();
    const Real*   upb  = this->solver()->upBound().get_const_ptr();

    const Real*   cvec = this->solver()->coPvec().get_const_ptr();
    const Real*   cupd = this->solver()->coPvec().delta().values();
    const IdxSet& cidx = this->solver()->coPvec().idx();
    const Real*   lcb  = this->solver()->lcBound().get_const_ptr();
    const Real*   ucb  = this->solver()->ucBound().get_const_ptr();

    Real epsilon = this->solver()->epsilon();
    Real val     = max;
    int  pnum    = -1;
    int  cnum    = -1;

    SPxId enterId;
    int   i;
    int   j;
    Real  x;

    // PARALLEL the j loops could be parallelized
    if (val > 0)
      {
        for( j = 0; j < pidx.size(); ++j )
          {
            i = pidx.index(j);
            x = pupd[i];

            if (x > epsilon)
              {
                if (upb[i] < infinity)
                  {
                    Real y = (upb[i] - pvec[i] + this->delta) / x;
                        
                    if (y < val)
                      {
                        enterId = this->solver()->id(i);
                        val     = y;
                        pnum    = j;
                      }
                  }
              }
            else if (x < -epsilon)
              {
                if (lpb[i] > -infinity)
                  {
                    Real y = (lpb[i] - pvec[i] - this->delta) / x;

                    if (y < val)
                      {
                        enterId = this->solver()->id(i);
                        val     = y;
                        pnum    = j;
                      }
                  }
              }
          }
        for( j = 0; j < cidx.size(); ++j )
          {
            i = cidx.index(j);
            x = cupd[i];

            if (x > epsilon)
              {
                if (ucb[i] < infinity)
                  {
                    Real y = (ucb[i] - cvec[i] + this->delta) / x;

                    if (y < val)
                      {
                        enterId = this->solver()->coId(i);
                        val     = y;
                        cnum    = j;
                      }
                  }
              }
            else if (x < -epsilon)
              {
                if (lcb[i] > -infinity)
                  {
                    Real y = (lcb[i] - cvec[i] - this->delta) / x;

                    if (y < val)
                      {
                        enterId = this->solver()->coId(i);
                        val     = y;
                        cnum    = j;
                      }
                  }
              }
          }
        if (cnum >= 0)
          {
            i   = cidx.index(cnum);
            x   = cupd[i];
            val = (x > epsilon) ? ucb[i] : lcb[i];
            val = (val - cvec[i]) / x;
          }
        else if (pnum >= 0)
          {
            i   = pidx.index(pnum);
            x   = pupd[i];
            val = (x > epsilon) ? upb[i] : lpb[i];
            val = (val - pvec[i]) / x;
          }
      }
    else
      {
        for( j = 0; j < pidx.size(); ++j )
          {
            i = pidx.index(j);
            x = pupd[i];

            if (x > epsilon)
              {
                if (lpb[i] > -infinity)
                  {
                    Real y = (lpb[i] - pvec[i] - this->delta) / x;

                    if (y > val)
                      {
                        enterId = this->solver()->id(i);
                        val     = y;
                        pnum    = j;
                      }
                  }
              }
            else if (x < -epsilon)
              {
                if (upb[i] < infinity)
                  {
                    Real y = (upb[i] - pvec[i] + this->delta) / x;

                    if (y > val)
                      {
                        enterId = this->solver()->id(i);
                        val     = y;
                        pnum    = j;
                      }
                  }
              }
          }
        for( j = 0; j < cidx.size(); ++j )
          {
            i = cidx.index(j);
            x = cupd[i];

            if (x > epsilon)
              {
                if (lcb[i] > -infinity)
                  {
                    Real y = (lcb[i] - cvec[i] - this->delta) / x;

                    if (y > val)
                      {
                        enterId = this->solver()->coId(i);
                        val     = y;
                        cnum    = j;
                      }
                  }
              }
            else if (x < -epsilon)
              {
                if (ucb[i] < infinity)
                  {
                    Real y = (ucb[i] - cvec[i] + this->delta) / x;

                    if (y > val)
                      {
                        enterId = this->solver()->coId(i);
                        val     = y;
                        cnum    = j;
                      }
                  }
              }
          }
        if (cnum >= 0)
          {
            i   = cidx.index(cnum);
            x   = cupd[i];
            val = (x < epsilon) ? ucb[i] : lcb[i];
            val = (val - cvec[i]) / x;
          }
        else if (pnum >= 0)
          {
            i   = pidx.index(pnum);
            x   = pupd[i];
            val = (x < epsilon) ? upb[i] : lpb[i];
            val = (val - pvec[i]) / x;
          }
      }

    if (enterId.isValid() && this->solver()->isBasic(enterId))
      {
        MSG_DEBUG( std::cout << "DDEFRT01 isValid() && isBasic(): max=" << max
                   << std::endl; )
          if (cnum >= 0)
            this->solver()->coPvec().delta().clearNum(cnum);
          else if( pnum >= 0 )
            this->solver()->pVec().delta().clearNum(pnum);
        return SPxDefaultRT<R>::selectEnter(max, 0, false);
      }

    MSG_DEBUG(
              if( !enterId.isValid() )
                std::cout << "DDEFRT02 !isValid(): max=" << max << ", x=" << x << std::endl;
              )
      max = val;

    return enterId;
  }

} // namespace soplex
