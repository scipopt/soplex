/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the class library                   */
/*       SoPlex --- the Sequential object-oriented simPlex.                  */
/*                                                                           */
/*    Copyright (C) 1997-1999 Roland Wunderling                              */
/*                  1997-2001 Konrad-Zuse-Zentrum                            */
/*                            fuer Informationstechnik Berlin                */
/*                                                                           */
/*  SoPlex is distributed under the terms of the ZIB Academic Licence.       */
/*                                                                           */
/*  You should have received a copy of the ZIB Academic License              */
/*  along with SoPlex; see the file COPYING. If not email to soplex@zib.de.  */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
#pragma ident "@(#) $Id: spxweightst.cpp,v 1.4 2001/12/28 14:55:13 bzfkocht Exp $"

/*      \Section{Complex Methods}
 */

/*  Import system include files
 */
#include <assert.h>
#include <iostream>

/*  and class header files
 */
#include "spxweightst.h"
#include "svset.h"
#include "sorter.h"

namespace soplex
{
#define EPS     1e-6
#define STABLE  1e-3

//@ ----------------------------------------------------------------------------
int SPxWeightST::isConsistent() const
{
   return rowWeight.isConsistent()
          && colWeight.isConsistent()
          && rowRight.isConsistent()
          && colUp.isConsistent()
          && SPxStarter::isConsistent();
}

//@ ----------------------------------------------------------------------------
/* \SubSection{Generating the Starting Basis}
    The generation of a starting basis works as follows: After setting up the
    preference arrays #weight# and #coWeight#, #Id#s are selected to be dual in
    a greedy manner.  Initially, the first #Id# is taken. Then the next #Id# is
    checked wheter its vector is linearly dependend of the vectors of the #Id#s
    allready selected.  If not, it is added to them. This is iterated until a
    full matrix has been constructed.
 
    Testing for linear independence is done very much along the lines of LU
    factorization. A vector is taken, and updated with all previous L-vectors.
    Then the maximal absolut element is selected as pivot element for computing
    the L-vector. This is stored for further processing.
 */

/*
    The following two functions set the status of |id| to primal or dual,
    respectively.
 */
void SPxWeightST::setPrimalStatus(SPxBasis::Desc& desc, const SoPlex& base, const SoPlex::Id& id)
{
   if (id.isSPxRowId())
   {
      int n = base.number(SoPlex::SPxRowId(id));
      if (base.rhs(n) >= SPxLP::infinity)
      {
         if (base.lhs(n) <= -SPxLP::infinity)
            desc.rowStatus(n) = SPxBasis::Desc::P_FREE;
         else
            desc.rowStatus(n) = SPxBasis::Desc::P_ON_LOWER;
      }
      else
      {
         if (base.lhs(n) <= -SPxLP::infinity)
            desc.rowStatus(n) = SPxBasis::Desc::P_ON_UPPER;
         else if (base.lhs(n) >= base.rhs(n) - base.epsilon())
            desc.rowStatus(n) = SPxBasis::Desc::P_FIXED;
         else if (rowRight[n])
            desc.rowStatus(n) = SPxBasis::Desc::P_ON_UPPER;
         else
            desc.rowStatus(n) = SPxBasis::Desc::P_ON_LOWER;
      }
   }
   else
   {
      int n = base.number(SoPlex::SPxColId(id));
      if (base.SPxLP::upper(n) >= SPxLP::infinity)
      {
         if (base.SPxLP::lower(n) <= -SPxLP::infinity)
            desc.colStatus(n) = SPxBasis::Desc::P_FREE;
         else
            desc.colStatus(n) = SPxBasis::Desc::P_ON_LOWER;
      }
      else
      {
         if (base.SPxLP::lower(n) <= -SPxLP::infinity)
            desc.colStatus(n) = SPxBasis::Desc::P_ON_UPPER;
         else if (base.SPxLP::lower(n) >= base.SPxLP::upper(n) - base.epsilon())
            desc.colStatus(n) = SPxBasis::Desc::P_FIXED;
         else if (colUp[n])
            desc.colStatus(n) = SPxBasis::Desc::P_ON_UPPER;
         else
            desc.colStatus(n) = SPxBasis::Desc::P_ON_LOWER;
      }
   }
}

static void setDualStatus(SPxBasis::Desc& desc, const SoPlex& base, const SoPlex::Id& id)
{
   if (id.isSPxRowId())
   {
      int n = base.number(SoPlex::SPxRowId(id));
      desc.rowStatus(n) = base.basis().dualRowStatus(n);
   }
   else
   {
      int n = base.number(SoPlex::SPxColId(id));
      desc.colStatus(n) = base.basis().dualColStatus(n);
   }
}

/*
    The following method initializes |pref| such that it contains the set of
    |SoPlex::Id|s ordered following |rowWeight| and |colWeight|. For the sorting
    we take the following approach: first we sort the rows, then the columns.
    Finally we perform a mergesort of both.
 */
struct Compare
{
   const SoPlex* base;
   const double* weight;
   double operator()(int i1, int i2)
   {
      return weight[i1] - weight[i2];
   }
};

static void initPrefs(DataArray < SoPlex::Id > & pref,
                       const SoPlex& base,
                       const DataArray < double > & rowWeight,
                       const DataArray < double > & colWeight
                    )
{
   int i, j, k;
   DataArray < int > row(base.nRows());
   DataArray < int > col(base.nCols());

   Compare compare;
   compare.base = &base;

   for (i = base.nRows(); i--;)
      row[i] = i;
   compare.weight = rowWeight.get_const_ptr();

   sorter_qsort(row.get_ptr(), row.size(), compare); // Sort rows

   for (i = base.nCols(); i--;)
      col[i] = i;
   compare.weight = colWeight.get_const_ptr();

   sorter_qsort(col.get_ptr(), col.size(), compare); // Sort column

   for (i = j = k = 0; k < pref.size();)            // merge sort
   {
      if (rowWeight[row[i]] < colWeight[col[j]])
      {
         pref[k++] = base.rId(row[i++]);
         if (i >= base.nRows())
            while (k < pref.size())
               pref[k++] = base.cId(col[j++]);
      }
      else
      {
         pref[k++] = base.cId(col[j++]);
         if (j >= base.nCols())
            while (k < pref.size())
               pref[k++] = base.rId(row[i++]);
      }
   }
   assert(i == base.nRows());
   assert(j == base.nCols());
}

void SPxWeightST::generate(SoPlex& base)
{
   SoPlex::Id tmpId;

   forbidden.reSize(base.dim());
   rowWeight.reSize(base.nRows());
   colWeight.reSize(base.nCols());
   rowRight.reSize (base.nRows());
   colUp.reSize (base.nCols());
   if (base.rep() == SoPlex::COLUMN)
   {
      weight = &colWeight;
      coWeight = &rowWeight;
   }
   else
   {
      weight = &rowWeight;
      coWeight = &colWeight;
   }
   assert(weight->size() == base.coDim());
   assert(coWeight->size() == base.dim());

   setupWeights(base);

   SPxBasis::Desc desc;
   desc.reSize(base.nRows(), base.nCols());

   DataArray < SoPlex::Id > pref(base.nRows() + base.nCols());
   initPrefs(pref, base, rowWeight, colWeight);

   int i, deltai, j, sel;
   for (i = forbidden.size(); --i >= 0;)
      forbidden[i] = 0;

   if (base.rep() == SoPlex::COLUMN)
   {
      i = 0;
      deltai = 1;
   }
   else
   {
      i = pref.size() - 1;
      deltai = -1;
   }
   {
      int dim = base.dim();
      double max = 0;

      for (; i >= 0 && i < pref.size(); i += deltai)
      {
         tmpId = pref[i];
         const SVector& bVec = base.vector(tmpId);
         sel = -1;
         if (bVec.size() == 1)
         {
            int idx = bVec.index(0);
            if (forbidden[idx] < 2)
            {
               sel = idx;
               dim += (forbidden[idx] > 0);
            }
         }
         else
         {
            max = bVec.maxAbs();

            int best = base.nRows();
            for (j = bVec.size(); --j >= 0;)
            {
               double x = bVec.value(j);
               int k = bVec.index(j);
               int l = base.coVector(k).size();
               if (!forbidden[k] && (x > STABLE*max || -x > STABLE*max)
                    && l < best)
               {
                  best = l;
                  sel = k;
               }
            }
         }

         if (sel >= 0)
         {
#ifndef NDEBUG
            if (pref[i].type() == SPxLP::Id::ROWID)
               std::cerr << " r" << base.number(pref[i]);
            else
               std::cerr << " c" << base.number(pref[i]);
#endif  // NDEBUG
            forbidden[sel] = 2;
            if (base.rep() == SoPlex::COLUMN)
               setDualStatus(desc, base, pref[i]);
            else
               setPrimalStatus(desc, base, pref[i]);

            for (j = bVec.size(); --j >= 0;)
            {
               double x = bVec.value(j);
               int k = bVec.index(j);
               if (!forbidden[k] && (x > EPS*max || -x > EPS*max))
               {
                  forbidden[k] = 1;
                  --dim;
               }
            }

            if (--dim == 0)
            {
               //@ for(++i; i < pref.size(); ++i)
               if (base.rep() == SoPlex::COLUMN)
               {
                  for (i += deltai; i >= 0 && i < pref.size(); i += deltai)
                     setPrimalStatus(desc, base, pref[i]);
                  for (i = forbidden.size(); --i >= 0;)
                  {
                     if (forbidden[i] < 2)
                        setDualStatus(desc, base, base.coId(i));
                  }
               }
               else
               {
                  for (i += deltai; i >= 0 && i < pref.size(); i += deltai)
                     setDualStatus(desc, base, pref[i]);
                  for (i = forbidden.size(); --i >= 0;)
                  {
                     if (forbidden[i] < 2)
                        setPrimalStatus(desc, base, base.coId(i));
                  }
               }
               break;
            }
         }
         else if (base.rep() == SoPlex::COLUMN)
            setPrimalStatus(desc, base, pref[i]);
         else
            setDualStatus(desc, base, pref[i]);
#ifndef NDEBUG
         {
            int n, m;
            for (n = 0, m = forbidden.size(); n < forbidden.size(); ++n)
               m -= (forbidden[n] != 0);
            assert(m == dim);
         }
#endif  // NDEBUG
      }
      assert(dim == 0);
   }

   base.load(desc);
#ifdef  TEST
   base.init();

   int changed = 0;
   const Vector& pvec = base.pVec();
   for (i = pvec.dim() - 1; i >= 0; --i)
   {
      if (desc.colStatus(i) == SPxBasis::Desc::P_ON_UPPER
           && base.lower(i) > -SPxLP::infinity
           && pvec[i] > base.maxObj(i))
      {
         changed = 1;
         desc.colStatus(i) = SPxBasis::Desc::P_ON_LOWER;
      }
      else if (desc.colStatus(i) == SPxBasis::Desc::P_ON_LOWER
                && base.upper(i) < SPxLP::infinity
                && pvec[i] < base.maxObj(i))
      {
         changed = 1;
         desc.colStatus(i) = SPxBasis::Desc::P_ON_UPPER;
      }
   }

   if (changed)
   {
      std::cerr << "changed basis\n";
      base.load(desc);
   }
   else
      std::cerr << "nothing changed\n";
#endif  // TEST
}


//@ ----------------------------------------------------------------------------
/* \SubSection{Computation of Weights}
 */
void SPxWeightST::setupWeights(SoPlex& bse)
{
   const SoPlex& base = bse;
   const Vector& obj = bse.maxObj();
   const Vector& low = bse.SPxLP::lower();
   const Vector& up = bse.SPxLP::upper();
   const Vector& rhs = bse.rhs();
   const Vector& lhs = bse.lhs();
   double len1;
   double nne, ax, bx;
   double r_fixed, r_dbl_bounded, r_bounded, r_free;
   double c_fixed, c_dbl_bounded, c_bounded, c_free;
   double x, l, u, n;
   int i;

   x = 1;
   for (i = bse.nCols(); i-- > 0;)
   {
      if (up[i] < SPxLP::infinity)
      {
         if (up[i] > x)
            x = up[i];
         else if (-up[i] > x)
            x = -up[i];
      }
      if (low[i] > -SPxLP::infinity)
      {
         if (low[i] > x)
            x = low[i];
         else if (-low[i] > x)
            x = -low[i];
      }
   }
   for (i = bse.nRows(); i-- > 0;)
   {
      if (rhs[i] < SPxLP::infinity)
      {
         if (rhs[i] > x)
            x = rhs[i];
         else if (-rhs[i] > x)
            x = -rhs[i];
      }
      if (lhs[i] > -SPxLP::infinity)
      {
         if (lhs[i] > x)
            x = lhs[i];
         else if (-lhs[i] > x)
            x = -lhs[i];
      }
   }

   if (bse.rep() * bse.type() > 0)            // primal simplex
   {
      bx = 1e+0 / x;
      ax = 1e-3 / obj.maxAbs();
      nne = ax / lhs.dim();  // 1e-4 * ax;
      c_fixed = 1e+5;
      r_fixed = 1e+4;
      c_dbl_bounded = 1e+1;
      r_dbl_bounded = 0;
      c_bounded = 1e+1;
      r_bounded = 0;
      c_free = -1e+4;
      r_free = -1e+5;

      for (i = bse.nCols(); i-- > 0;)
      {
         n = nne * (bse.colVector(i).size() - 1);
         x = ax * obj[i];
         u = bx * up [i];
         l = bx * low[i];

         if (up[i] < SPxLP::infinity)
         {
            if (low[i] > up[i] - base.epsilon())
               colWeight[i] = c_fixed + n + (x > 0) ? x : -x;
            else if (low[i] > -SPxLP::infinity)
            {
               colWeight[i] = c_dbl_bounded + l - u + n;
               l = (l < 0) ? -l : l;
               u = (u < 0) ? -u : u;
               if (u < l)
               {
                  colUp[i] = 1;
                  colWeight[i] += x;
               }
               else
               {
                  colUp[i] = 0;
                  colWeight[i] -= x;
               }
            }
            else
            {
               colWeight[i] = c_bounded - u + x + n;
               colUp[i] = 1;
            }
         }
         else
         {
            if (low[i] > -SPxLP::infinity)
            {
               colWeight[i] = c_bounded + l + n - x;
               colUp[i] = 0;
            }
            else
               colWeight[i] = c_free + n - (x > 0) ? x : -x;
         }
      }

      for (i = bse.nRows(); i-- > 0;)
      {
         if (rhs[i] < SPxLP::infinity)
         {
            if (lhs[i] > rhs[i] - base.epsilon())
               rowWeight[i] = r_fixed;
            else if (lhs[i] > -SPxLP::infinity)
            {
               u = bx * rhs[i];
               l = bx * lhs[i];
               rowWeight[i] = r_dbl_bounded + l - u;
               l = (l < 0) ? -l : l;
               u = (u < 0) ? -u : u;
               rowRight[i] = (u < l);
            }
            else
            {
               rowWeight[i] = r_bounded - bx * rhs[i];
               rowRight[i] = 1;
            }
         }
         else
         {
            if (lhs[i] > -SPxLP::infinity)
            {
               rowWeight[i] = r_bounded + bx * lhs[i];
               rowRight[i] = 0;
            }
            else
               rowWeight[i] = r_free;
         }
      }
   }

   else
   {
      assert(bse.rep() * bse.type() < 0);           // dual simplex
      ax = 1e+0 / obj.maxAbs();
      bx = 1e-2 / x;
      nne = 1e-4 * bx;
      c_fixed = 1e+5;
      r_fixed = 1e+4;
      c_dbl_bounded = 1;
      r_dbl_bounded = 0;
      c_bounded = 0;
      r_bounded = 0;
      c_free = -1e+4;
      r_free = -1e+5;

      for (i = bse.nCols(); i-- > 0;)
      {
         n = nne * (bse.colVector(i).size() - 1);
         x = ax * obj[i];
         u = bx * up [i];
         l = bx * low[i];

         if (up[i] < SPxLP::infinity)
         {
            if (low[i] > up[i] - base.epsilon())
               colWeight[i] = c_fixed + n + (x > 0) ? x : -x;
            else if (low[i] > -SPxLP::infinity)
            {
               if (x > 0)
               {
                  colWeight[i] = c_dbl_bounded + x - u + n;
                  colUp[i] = 1;
               }
               else
               {
                  colWeight[i] = c_dbl_bounded - x + l + n;
                  colUp[i] = 0;
               }
            }
            else
            {
               colWeight[i] = c_bounded - u + x + n;
               colUp[i] = 1;
            }
         }
         else
         {
            if (low[i] > -SPxLP::infinity)
            {
               colWeight[i] = c_bounded - x + l + n;
               colUp[i] = 0;
            }
            else
               colWeight[i] = c_free + n - (x > 0) ? x : -x;
         }
      }

      for (i = bse.nRows(); i-- > 0;)
      {
         len1 = 1;  // (bse.rowVector(i).length() + bse.epsilon());
         n = 0;  // nne * (bse.rowVector(i).size() - 1);
         u = bx * len1 * rhs[i];
         l = bx * len1 * lhs[i];
         x = ax * len1 * (obj * bse.rowVector(i));

         if (rhs[i] < SPxLP::infinity)
         {
            if (lhs[i] > rhs[i] - base.epsilon())
               rowWeight[i] = r_fixed + n + (x > 0) ? x : -x;
            else if (lhs[i] > -SPxLP::infinity)
            {
               if (x > 0)
               {
                  rowWeight[i] = r_dbl_bounded + x - u + n;
                  rowRight[i] = 1;
               }
               else
               {
                  rowWeight[i] = r_dbl_bounded - x + l + n;
                  rowRight[i] = 0;
               }
            }
            else
            {
               rowWeight[i] = r_bounded - u + n + x;
               rowRight[i] = 1;
            }
         }
         else
         {
            if (lhs[i] > -SPxLP::infinity)
            {
               rowWeight[i] = r_bounded + l + n - x;
               rowRight[i] = 0;
            }
            else
               rowWeight[i] = r_free + n - (x > 0) ? x : -x;
         }
      }
   }
}
} // namespace soplex

//-----------------------------------------------------------------------------
//Emacs Local Variables:
//Emacs mode:c++
//Emacs c-basic-offset:3
//Emacs tab-width:8
//Emacs indent-tabs-mode:nil
//Emacs End:
//-----------------------------------------------------------------------------
