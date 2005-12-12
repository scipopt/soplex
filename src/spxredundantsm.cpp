/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the class library                   */
/*       SoPlex --- the Sequential object-oriented simPlex.                  */
/*                                                                           */
/*    Copyright (C) 1997-1999 Roland Wunderling                              */
/*                  1997-2002 Konrad-Zuse-Zentrum                            */
/*                            fuer Informationstechnik Berlin                */
/*                                                                           */
/*  SoPlex is distributed under the terms of the ZIB Academic Licence.       */
/*                                                                           */
/*  You should have received a copy of the ZIB Academic License              */
/*  along with SoPlex; see the file COPYING. If not email to soplex@zib.de.  */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
#pragma ident "@(#) $Id: spxredundantsm.cpp,v 1.33 2005/12/12 20:22:50 bzforlow Exp $"

//#define DEBUGGING 1

#include <iostream>

#include "spxdefines.h"
#include "spxredundantsm.h"
#include "dataarray.h"
#include "sorter.h"
#include "spxout.h"

#define DISPERSE(x) (1664525U * (x) + 1013904223U)

namespace soplex
{
void SPxRedundantSM::fixColumn(SPxLP& lp, int i)
{
   METHOD( "SPxRedundantSM::fixColumn" );

   assert(EQrel(lp.lower(i), lp.upper(i), deltaBnd()));

   Real x = lp.lower(i);

   MSG_INFO3( spxout << "IREDSM01 fixed col " << i 
                        << " lower= " << std::setprecision(16) << lp.lower(i)
                        << " upper= " << std::setprecision(16) << lp.upper(i)
                        << std::endl; )

   m_pval.add(m_cperm[i], x);
   
   if (isNotZero(x, epsZero()))
   {
      const SVector& col = lp.colVector(i);

      for(int j = 0; j < col.size(); ++j)
      {
         int k = col.index(j);
         
         if (lp.rhs(k) < infinity)
         {
            Real y     = x * col.value(j);
            Real scale = maxAbs(lp.rhs(k), y);
            Real rhs   = (lp.rhs(k) / scale) - (y / scale);
                  
            if (isZero(rhs, epsZero()))
               rhs = 0.0;
            else
               rhs *= scale;

            MSG_INFO3( spxout << "IREDSM02 \trhs " << k 
                                 << " r= "   << std::setprecision(16) << rhs 
                                 << " rhs= " << std::setprecision(16) << lp.rhs(k) 
                                 << " x= "   << std::setprecision(16) << col.value(j) 
                                 << std::endl; )

            lp.changeRhs(k, rhs);
         }
         if (lp.lhs(k) > -infinity)
         {
            Real y     = x * col.value(j);
            Real scale = maxAbs(lp.lhs(k), y);
            Real lhs   = (lp.lhs(k) / scale) - (y / scale);

            if (isZero(lhs, epsZero()))
               lhs = 0.0;
            else
               lhs *= scale;
                  
            MSG_INFO3( spxout << "IREDSM03 \tlhs " << k 
                                 << " l= "   << std::setprecision(16) << lhs 
                                 << " lhs= " << std::setprecision(16) << lp.lhs(k) 
                                 << " x= "   << std::setprecision(16) << col.value(j) 
                                 << std::endl; )

            lp.changeLhs(k, lhs);
         }
      }
   }
}

SPxSimplifier::Result SPxRedundantSM::redundantRows(SPxLP& lp, bool& again)
{
   DataArray<int>     rem(lp.nRows());
   DataArray<RowHash> rowhash(lp.nRows());
   int                remRows = 0;
   int                remNzos = 0;
   int                chgLRhs = 0;
   int                chgBnds = 0;
   int                i;
   int                j;

   for(i = 0; i < lp.nRows(); ++i )
   {
      //      std::cout << "row " << i << std::endl;

      const SVector& row    = lp.rowVector(i);       
      Real           upbnd  = 0.0;
      Real           lobnd  = 0.0;
      int            upcnt  = 0;
      int            locnt  = 0;

      rowhash[i].row = i;
      rowhash[i].hid = DISPERSE(row.size());
      rem[i]         = 0;

      for(j = 0; j < row.size(); ++j )
      {
         Real x = row.value(j);
         int  k = row.index(j);

         rowhash[i].hid = rowhash[i].hid * 31 + k; 
         
         // std::cout << "\t\tx= " << x << " lower= " << lp.lower(k) 
         // << " upper= " << lp.upper(k) << std::endl;

         // warn since this unhandled case may slip through unnoticed otherwise
         ASSERT_WARN( "WREDSM28", isNotZero(x, epsZero()) );

         if (x > 0.0)
         {
            if (lp.upper(k) >= infinity)
               upcnt++;
            else
               upbnd += lp.upper(k) * x;

            if (lp.lower(k) <= -infinity)
               locnt++;
            else
               lobnd += lp.lower(k) * x;
         }
         else if (x < 0.0)
         {
            if (lp.upper(k) >= infinity)
               locnt++;
            else
               lobnd += lp.upper(k) * x;

            if (lp.lower(k) <= -infinity)
               upcnt++;
            else
               upbnd += lp.lower(k) * x;
         }
      }
      // infeasible ?
      if (  (LTrel(lp.rhs(i), lobnd, deltaBnd()) && locnt == 0) 
         || (GTrel(lp.lhs(i), upbnd, deltaBnd()) && upcnt == 0))
      {
         MSG_INFO3( spxout << "IREDSM04 infeasible row " << i 
                              << " lo= " << lobnd
                              << " up= " << upbnd 
                              << " lhs= " << lp.lhs(i) 
                              << " rhs= " << lp.rhs(i)
                              << std::endl; )
         return INFEASIBLE;
      }
      // forcing equality constraint ?
      if (EQ(lp.lhs(i), lp.rhs(i), deltaBnd()))
      {
         // all fixed on upper bound ?
         if (upcnt == 0 && EQrel(lp.rhs(i), upbnd, deltaBnd()))
         {
            MSG_INFO3( spxout << "IREDSM05 rhs fixed on upbnd row " << i
                                 << " rhs= " << lp.rhs(i)
                                 << " up= " << upbnd 
                                 << std::endl; )

            for(j = 0; j < row.size(); ++j )
            {
               Real x = row.value(j);
               int  k = row.index(j);

               ASSERT_WARN( "WREDSM29", isNotZero(x, epsZero()) );
               ASSERT_WARN( "WREDSM30", GE(lp.lower(k), 0.0) || LE(lp.upper(k), 0.0) );

               if (x > 0.0 && GE(lp.lower(k), 0.0))
                  lp.changeLower(k, lp.upper(k));
               else
                  lp.changeUpper(k, lp.lower(k));               
            }
            rem[i]         = -1;
            rowhash[i].row = -1;
            rowhash[i].hid = 0;
            remRows++;
            remNzos += row.size();
            continue;
         }
         // all fixed on lower bound ?
         if (locnt == 0 && EQrel(lp.lhs(i), lobnd, deltaBnd()))
         {
            MSG_INFO3( spxout << "IREDSM06 rhs fixed on lowbnd row " << i
                                 << " lhs= " << lp.lhs(i)
                                 << " lo= " << lobnd 
                                 << std::endl; )

            for(j = 0; j < row.size(); ++j )
            {
               Real x = row.value(j);
               int  k = row.index(j);

               ASSERT_WARN( "WREDSM31", isNotZero(x, epsZero()) );
               ASSERT_WARN( "WREDSM32", GE(lp.lower(k), 0.0) || LE(lp.upper(k), 0.0) );

               if (x > 0.0 && GE(lp.lower(k), 0.0))
                  lp.changeUpper(k, lp.lower(k));
               else
                  lp.changeLower(k, lp.upper(k));
            }
            rem[i]         = -1;
            rowhash[i].row = -1;
            rowhash[i].hid = 0;
            remRows++;
            remNzos += row.size();
            continue;
         }
      }

      // redundant rhs ?
      if (lp.rhs(i) <  infinity && upcnt == 0 && GErel(lp.rhs(i), upbnd, deltaBnd()))
      {
         MSG_INFO3( spxout << "IREDSM07 redundant rhs row " << i
                              << " rhs= " << lp.rhs(i)
                              << " up= "  << upbnd 
                              << std::endl; )

         lp.changeRhs(i, infinity);
         chgLRhs++;
      }
      // redundant lhs ?
      if (lp.lhs(i) > -infinity && locnt == 0 && LErel(lp.lhs(i), lobnd, deltaBnd()))
      {
         MSG_INFO3( spxout << "IREDSM08 redundant lhs row " << i
                              << " lhs= " << lp.lhs(i)
                              << " lo= "  << lobnd 
                              << std::endl; )

         lp.changeLhs(i, -infinity);
         chgLRhs++;
      }
      // unconstraint constraints are also handled by simpleRows
      // but since they might come up here we do it again.
      if (lp.rhs(i) >= infinity && lp.lhs(i) <= -infinity)
      {
         MSG_INFO3( spxout << "IREDSM09 unconstraint row " << i 
                              << " removed" << std::endl; )

         rem[i]         = -1;
         rowhash[i].row = -1;
         rowhash[i].hid = 0;
         remRows++;
         remNzos += row.size();
         continue;
      }
      if (upcnt <= 1 || locnt <= 1)
      {
         for(j = 0; j < row.size(); ++j )
         {
            Real x = row.value(j);
            int  k = row.index(j);

            // warn since this unhandled case may slip through unnoticed otherwise
            ASSERT_WARN( "WREDSM33", isNotZero(x) );

            if (x > 0.0)
            {
               if (lp.lhs(i) > -infinity && lp.lower(k) > -infinity && upcnt <= 1 && NErel(lp.lhs(i), upbnd))
               {
                  Real y     = -infinity;
                  Real scale = maxAbs(lp.lhs(i), upbnd);
                  Real z     = (lp.lhs(i) / scale) - (upbnd / scale);

                  assert(upcnt > 0 || lp.upper(k) < infinity);
                  
                  if (upcnt == 0)
                     y = lp.upper(k) + z * scale / x;
                  else if (lp.upper(k) >= infinity)
                     y = z * scale / x;

                  if (isZero(y, deltaBnd()))
                     y = 0.0;

                  if (GE(y, lp.lower(k)))
                  {
                     MSG_INFO3( spxout << "IREDSM10 dominated bound row " << i
                                          << " col " << k
                                          << " removed y= " << y
                                          << " lower= " << lp.lower(k)
                                          << std::endl; )

                     locnt++;
                     lobnd -= lp.lower(k) * x;
                     lp.changeLower(k, -infinity);
                     chgBnds++;
                  }
               }
               if (lp.rhs(i) < infinity && lp.upper(k) < infinity && locnt <= 1 && NErel(lp.rhs(i), lobnd))
               {
                  Real y = infinity;
                  Real scale = maxAbs(lp.rhs(i), lobnd);
                  Real z     = (lp.rhs(i) / scale) - (lobnd / scale);
                  
                  assert(locnt > 0 || lp.lower(k) > -infinity);

                  if (locnt == 0)
                     y = lp.lower(k) + z * scale / x;
                  else if (lp.lower(k) <= -infinity)
                     y = z * scale / x;
                  
                  if (isZero(y, deltaBnd()))
                     y = 0.0;

                  if (LE(y, lp.upper(k)))
                  {
                     MSG_INFO3( spxout << "IREDSM11 dominated bound row " << i
                                          << " col " << k
                                          << " removed y= " << y
                                          << " lower= " << lp.lower(k)
                                          << std::endl; )

                     upcnt++;
                     upbnd -= lp.upper(k) * x;
                     lp.changeUpper(k, infinity);
                     chgBnds++;
                  }
               }
            }
            else if (x < 0.0)
            {
               if (lp.lhs(i) >= -infinity && lp.upper(k) < infinity && upcnt <= 1 && NErel(lp.lhs(i), upbnd))
               {
                  Real y = infinity;
                  Real scale = maxAbs(lp.lhs(i), upbnd);
                  Real z     = (lp.lhs(i) / scale) - (upbnd / scale);
                  
                  assert(upcnt > 0 || lp.lower(k) > -infinity);

                  if (upcnt == 0)
                     y = lp.lower(k) + z * scale / x;
                  else if (lp.lower(k) <= -infinity)
                     y = z * scale / x;
                  
                  if (isZero(y, deltaBnd()))
                     y = 0.0;

                  if (LE(y, lp.upper(k)))
                  {
                     MSG_INFO3( spxout << "IREDSM12 dominated bound row " << i
                                          << " col " << k
                                          << " removed y= " << y
                                          << " lower= " << lp.lower(k)
                                          << std::endl; )

                     locnt++;
                     lobnd -= lp.upper(k) * x;
                     lp.changeUpper(k, infinity);
                     chgBnds++;
                  }
               }
               if (lp.rhs(i) <= infinity && lp.lower(k) > -infinity && locnt <= 1 && NErel(lp.rhs(i), lobnd))
               {
                  Real y = -infinity;
                  Real scale = maxAbs(lp.rhs(i), lobnd);
                  Real z     = (lp.rhs(i) / scale) - (lobnd / scale);

                  assert(locnt > 0 || lp.upper(k) < infinity);

                  if (locnt == 0)
                     y = lp.upper(k) + z * scale / x;
                  else if (lp.upper(k) >= infinity)
                     y = z * scale / x;
                  
                  if (isZero(y, deltaBnd()))
                     y = 0.0;

                  if (GE(y, lp.lower(k)))
                  {
                     MSG_INFO3( spxout << "IREDSM13 dominated bound row " << i
                                          << " col " << k
                                          << " removed y= " << y
                                          << " lower= " << lp.lower(k)
                                          << std::endl; )

                     upcnt++;
                     upbnd -= lp.lower(k) * x;
                     lp.changeLower(k, -infinity);
                     chgBnds++;
                  }
               }
            }
         }
      }
   }
   // --- Check for dublicate rows ------------------------------------------------------
   RowHash compare;

   sorter_qsort(rowhash.get_ptr(), rowhash.size(), compare);

   for(i = 0; i < lp.nRows() - 1; ++i)
   {
      if (rowhash[i].row < 0 || rowhash[i].hid != rowhash[i + 1].hid)
         continue;

      int ri1 = rowhash[i    ].row;
      int ri2 = rowhash[i + 1].row;

      assert(ri1 >= 0);
      assert(ri1 < lp.nRows());
      // This has to be true, because if case of the sorting.
      assert(ri2 >= 0);
      assert(ri2 < lp.nRows());
      assert(rowhash[i].hid == rowhash[i + 1].hid);

      if (lp.rowVector(ri1).size() != lp.rowVector(ri2).size())
         continue;

      SVector row1(lp.rowVector(ri1));       
      SVector row2(lp.rowVector(ri2));       
         
      row1.sort();
      row2.sort();
      
      if (row1.index(0) != row2.index(0))
         break;

      Real x1 = row1.value(0);
      Real x2 = row2.value(0);
      Real alpha;

      if (isZero(x1, epsZero()))
         alpha = 1.0;
      else if (isZero(x2, epsZero()))
         alpha = 1.0;
      else
         alpha = x1 / x2;

      ASSERT_WARN( "WREDSM34", isNotZero(alpha) );

      for(j = 0; j < row1.size(); ++j)
      {
         int  k1 = row1.index(j);
         int  k2 = row2.index(j);

         if (k1 != k2)
            break;

         x1 = row1.value(j);
         x2 = row2.value(j);

         if (NErel(x1, x2 * alpha))
            break;
      }
      if (j < row1.size())
         continue;
      
      // ok, now row[i] = alpha * row[i+1]

      Real maxlhs;
      Real minrhs;

      if (alpha > 0.0)
      {
         // compute maxlhs    
         if (lp.lhs(ri1) <= -infinity)
            maxlhs = lp.lhs(ri2);
         else if (lp.lhs(ri2) <= -infinity)
            maxlhs = lp.lhs(ri1);
         else
            maxlhs = (lp.lhs(ri1) / alpha >= lp.lhs(ri2)) ? lp.lhs(ri1) / alpha : lp.lhs(ri2);
         
         // compute minrhs
         if (lp.rhs(ri1) >= infinity)
            minrhs = lp.rhs(ri2);
         else if (lp.rhs(ri2) >= infinity)
            minrhs = lp.rhs(ri1);
         else
            minrhs = (lp.rhs(ri1) / alpha <= lp.rhs(ri2)) ? lp.rhs(ri1) / alpha : lp.rhs(ri2);
      }
      else
      {
         assert(alpha < 0.0);

         // compute maxlhs    
         if (lp.rhs(ri1) >= infinity)
            maxlhs = lp.lhs(ri2);
         else if (lp.lhs(ri2) <= -infinity)
            maxlhs = -lp.rhs(ri1);
         else
            maxlhs = (lp.rhs(ri1) / alpha >= lp.lhs(ri2)) ? lp.rhs(ri1) / alpha : lp.lhs(ri2);

         // compute minrhs
         if (lp.lhs(ri1) <= -infinity)
            minrhs = lp.rhs(ri2);
         else if (lp.rhs(ri2) >= infinity)
            minrhs = -lp.lhs(ri1);
         else
            minrhs = (lp.lhs(ri1) / alpha <= lp.rhs(ri2)) ? lp.lhs(ri1) / alpha : lp.rhs(ri2);
      }
      MSG_INFO3( spxout << "IREDSM14 duplicate rows " << ri1 << "/" << ri2
                           << " alpha= " << alpha
                           << " lhs= " << lp.lhs(ri1)
                           << " rhs= " << lp.rhs(ri1)
                           << " lhs= " << lp.lhs(ri2)
                           << " rhs= " << lp.rhs(ri2)
                           << " maxlhs= " << maxlhs
                           << " minrhs= " << minrhs
                           << std::endl; )

      lp.changeLhs(ri2, maxlhs);
      lp.changeRhs(ri2, minrhs);

      rem[ri1] = -1;
      remRows++;
      remNzos += row1.size();
   }
   if (remRows > 0)
      removeRows(lp, rem, remRows);

   if (remRows + remNzos + chgLRhs + chgBnds > 0)
   {
      again      = true;
      m_remRows += remRows;
      m_remNzos += remNzos;
      m_chgLRhs += chgLRhs;
      m_chgBnds += chgBnds;

      MSG_INFO2( spxout << "IREDSM15 redundant row simplifier removed "
                           << remRows << " rows, "
                           << remNzos << " nzos, changed "
                           << chgBnds << " col bounds, " 
                           << chgLRhs << " row bounds"
                           << std::endl; )
   }
   return OKAY;
}

SPxSimplifier::Result SPxRedundantSM::redundantCols(SPxLP& lp, bool& again)
{
   DataArray<int> rem(lp.nCols());
   DataArray<int> tmp(lp.nCols());
   int            remCols = 0;
   int            remNzos = 0;
   int            chgBnds = 0;

   for( int i = 0; i < lp.nCols(); ++i )
   {
      rem[i] = 0;

      const SVector& col = lp.colVector(i);

      if (NErel(lp.upper(i), lp.lower(i), deltaBnd()))
      {
         // test if all coefficents are going in one direction
         int upcnt = 0;
         int locnt = 0;
         int j;

         for(j = 0; j < col.size(); ++j )
         {            
            if (upcnt > 0 && locnt > 0)
               break;

            Real x = col.value(j);
            int  k = col.index(j);

            // warn since this unhandled case may slip through unnoticed otherwise
            ASSERT_WARN( "WREDSM35", isNotZero(x) );

            if (x > 0.0)
            {
               upcnt += (lp.rhs(k) <  infinity) ? 1 : 0;
               locnt += (lp.lhs(k) > -infinity) ? 1 : 0;
            }
            else if (x < 0.0)
            {
               locnt += (lp.rhs(k) <  infinity) ? 1 : 0;
               upcnt += (lp.lhs(k) > -infinity) ? 1 : 0;
            }
         }
         Real maxobj = lp.maxObj(i);

         // max -3 x
         // s.t. 5 x <= 8
         if (locnt == 0 && LT(maxobj, 0.0, epsZero()))
         {
            if (lp.lower(i) <= -infinity)
               return UNBOUNDED;           // LP is unbounded

            // fix column to lower bound
            lp.changeUpper(i, lp.lower(i));
         }
         // max  3 x
         // s.t. 5 x >= 8
         else if (upcnt == 0 && GT(maxobj, 0.0, epsZero()))
         {
            if (lp.upper(i) >= infinity)
               return UNBOUNDED;           // LP is unbounded

            // fix column to upper bound
            lp.changeLower(i, lp.upper(i));
         }
         else if (isZero(maxobj, epsZero()))
         {
            //TODO free variable with zero objective can be removed, or?
            //see simpleCols

            upcnt += (lp.upper(i) <  infinity) ? 1 : 0;
            locnt += (lp.lower(i) > -infinity) ? 1 : 0;

            if (locnt == 0)
            {
               // make variable free
               lp.changeUpper(i, infinity);
               chgBnds++;
#ifndef NDEBUG
               for(j = 0; j < col.size(); ++j )
               {
                  if (col.value(j) < 0.0)
                     assert(lp.rhs(col.index(j)) >= infinity); // lp.changeRhs(col.index(j), infinity);
                  else
                     assert(lp.lhs(col.index(j)) <= -infinity); // lp.changeLhs(col.index(j), -infinity);
               }
#endif // NDEBUG
            }
            if (upcnt == 0)
            {
               // make variable free
               lp.changeLower(i, -infinity);
               chgBnds++;
#ifndef NDEBUG
               for(j = 0; j < col.size(); ++j )
               {
                  if (col.value(j) > 0.0)
                     assert(lp.rhs(col.index(j)) >= infinity); // lp.changeRhs(col.index(j), infinity);
                  else
                     assert(lp.lhs(col.index(j)) <= -infinity); // lp.changeLhs(col.index(j), -infinity);
               }
#endif // NDEBUG
            }
         }
      }
      // remove fixed variables
      if (EQrel(lp.lower(i), lp.upper(i), deltaBnd()))
      {
         fixColumn(lp, i);

         rem[i] = -1;
         remCols++;
         remNzos += col.size();
      }
   }
   if (remCols > 0)
      removeCols(lp, rem, remCols);

   if (remCols + remNzos + chgBnds > 0)
   {
      again      = true;
      m_remCols += remCols;
      m_remNzos += remNzos;
      m_chgBnds += chgBnds;

      MSG_INFO2( spxout << "IREDSM16 redundant col simplifier removed "
                           << remCols << " cols, "
                           << remNzos << " nzos, changed "
                           << chgBnds << " col bounds" 
                           << std::endl; )
   }
   return OKAY;
}

/* Handle rows -------------------------------------------------------
 */
SPxSimplifier::Result SPxRedundantSM::simpleRows(SPxLP& lp, bool& again)
{
   DataArray<int> rem(lp.nRows());
   DataArray<int> tmp(lp.nRows());
   int            remRows = 0;
   int            remNzos = 0;

   for(int i = 0; i < lp.nRows(); ++i)
   {
      const SVector& row = lp.rowVector(i);
      rem[i]             = 0;

      // infeasible range row
      if (LTrel(lp.rhs(i), lp.lhs(i), deltaBnd()))
      {
         MSG_INFO3( spxout << "IREDSM17 infeasible row " << i 
                              <<"  lhs= " << lp.lhs(i) 
                              << " rhs= " << lp.rhs(i) 
                              << std::endl; )
         return INFEASIBLE;
      }
      // empty row ?
      if (row.size() == 0)
      {
         MSG_INFO3( spxout << "IREDSM18 empty row " << i; )

         if (LT(lp.rhs(i), 0.0, deltaBnd()) || GT(lp.lhs(i), 0.0, deltaBnd()))
         {
            MSG_INFO3( spxout << " infeasible lhs= " << lp.lhs(i) 
                                 << " rhs= " << lp.rhs(i) << std::endl; )
            return INFEASIBLE;
         }         
         MSG_INFO3( spxout << " removed" << std::endl; )

         rem[i] = -1;
         remRows++;
         continue;
      }
      // unconstraint constraint ?
      if (lp.rhs(i) >= infinity && lp.lhs(i) <= -infinity)
      {
         MSG_INFO3( spxout << "IREDSM19 unconstraint row " << i 
                              << " removed" << std::endl; )

         rem[i] = -1;
         remRows++;
         remNzos += row.size();
         continue;
      }
      // row singleton
      if (row.size() == 1)
      {
         Real x = row.value(0);
         int  j = row.index(0);
         Real up;
         Real lo;

         MSG_INFO3( spxout << "IREDSM20 row singleton " << i 
                              << " x= " << x 
                              << " lhs= " << lp.lhs(i) 
                              << " rhs= " << lp.rhs(i); )

         if (GT(x, 0.0, epsZero()))           // x > 0
         {
            up = (lp.rhs(i) >=  infinity) ?  infinity : (lp.rhs(i) / x);
            lo = (lp.lhs(i) <= -infinity) ? -infinity : (lp.lhs(i) / x);
         }
         else if (LT(x, 0.0, epsZero()))      // x < 0
         {
            lo = (lp.rhs(i) >=  infinity) ? -infinity : (lp.rhs(i) / x);
            up = (lp.lhs(i) <= -infinity) ?  infinity : (lp.lhs(i) / x);
         }
         else if (LT(lp.rhs(i), 0.0, deltaBnd()) || GT(lp.lhs(i), 0.0, deltaBnd()))  
         {
            // x == 0 rhs/lhs != 0
            MSG_INFO3( spxout << " infeasible" << std::endl; )

            return INFEASIBLE;
         }
         else                     // x == 0
         {
            lo = lp.lower(j);
            up = lp.upper(j);
         }
         rem[i] = -1;
         remRows++;
         remNzos++;
         
         if (isZero(lo, epsZero()))
            lo = 0.0;
         if (isZero(up, epsZero()))
            up = 0.0;
         
         ASSERT_WARN( "WREDSM36", LErel(lp.lower(j), lp.upper(j)) );

         MSG_INFO3( spxout << " removed lo= " << lo
                              << " up= " << up
                              << " lower= " << lp.lower(j)
                              << " upper= " << lp.upper(j)
                              << std::endl; )

         if (LT(up, lp.upper(j), epsZero()))
            lp.changeUpper(j, up);
         if (GT(lo, lp.lower(j), epsZero()))
            lp.changeLower(j, lo);

         assert(LErel(lp.lower(j), lp.upper(j), epsZero()));
      }
   }
   assert(remRows != 0 || remNzos == 0);

   if (remRows > 0)
   {
      removeRows(lp, rem, remRows);

      again      = true;
      m_remRows += remRows;
      m_remNzos += remNzos;

      MSG_INFO2( spxout << "IREDSM21 simple row simplifier removed "
                           << remRows << " rows, "
                           << remNzos << " nzos"
                           << std::endl; )
   }
   return OKAY;
}

/* Handle columns -----------------------------------------------------
 */
SPxSimplifier::Result SPxRedundantSM::simpleCols(SPxLP& lp, bool& again)
{
   DataArray<int> rem(lp.nCols());
   DataArray<int> tmp(lp.nCols());
   int            remCols = 0;
   int            remNzos = 0;

   for(int i = 0; i < lp.nCols(); ++i)
   {
      const SVector& col = lp.colVector(i);
      rem[i]             = 0;

      // Empty column ? 
      if (col.size() == 0)
      {
         MSG_INFO3( spxout << "IREDSM22 empty column " << i 
                              << " maxObj= " << lp.maxObj(i)
                              << " lower= " << lp.lower(i)
                              << " upper= " << lp.upper(i); )

         if (GT(lp.maxObj(i), 0.0, epsZero()))
         {
            if (lp.upper(i) >= infinity)
            {
               MSG_INFO3( spxout << " unbounded" << std::endl; )

               return UNBOUNDED;
            }
            m_pval.add(m_cperm[i], lp.upper(i));
         }
         else if (LT(lp.maxObj(i), 0.0, epsZero()))
         {
            if (lp.lower(i) <= -infinity)
            {
               MSG_INFO3( spxout << " unbounded" << std::endl; )

               return UNBOUNDED;
            }
            m_pval.add(m_cperm[i], lp.lower(i));
         }
         else 
         {
            assert(isZero(lp.maxObj(i), epsZero()));
            // any value within the bounds is ok
            if (lp.lower(i) > -infinity)
               m_pval.add(m_cperm[i], lp.lower(i));
            else if (lp.upper(i) < infinity)
               m_pval.add(m_cperm[i], lp.upper(i));
            else
               m_pval.add(m_cperm[i], 0.0);
         }
         MSG_INFO3( spxout << " removed" << std::endl; )

         rem[i] = -1;
         remCols++;
         continue;
      }

      // infeasible bounds ?
      if (GTrel(lp.lower(i), lp.upper(i), deltaBnd()))
      {
         MSG_INFO3( spxout << "IREDSM23 infeasible bounds column " << i 
                              << " lower= " << lp.lower(i)
                              << " upper= " << lp.upper(i)
                              << std::endl; )
         return INFEASIBLE;
      }
      // Fixed column ?
      if (EQrel(lp.lower(i), lp.upper(i), deltaBnd()))
      {
         fixColumn(lp, i);

         rem[i] = -1;
         remCols++;
         remNzos += col.size();
      }
      // a lot of this is also done in redundantRows.
      else if (col.size() == 1)
      {
         Real x = col.value(0);
         int  j = col.index(0);
         
         ASSERT_WARN( "WREDSM37", isNotZero(x) );

         if (x > 0.0)
         {
            // max -3 x
            // s.t. 5 x <= 8
            // l <= x
            if (lp.lhs(j) <= -infinity && lp.rhs(j) < infinity && LE(lp.maxObj(i), 0.0, epsZero()) 
               && lp.lower(i) > -infinity)
            {
               lp.changeUpper(i, lp.lower(i));
               fixColumn(lp, i);
               rem[i] = -1;
               remCols++;
               remNzos++;
               continue;
            }
            // max  3 x
            // s.t. 5 x >= 8
            // x <= u
            if (lp.lhs(j) > -infinity && lp.rhs(j) >= infinity && GE(lp.maxObj(i), 0.0, epsZero()) 
               && lp.upper(i) < infinity)
            {
               lp.changeLower(i, lp.upper(i));
               fixColumn(lp, i);
               rem[i] = -1;
               remCols++;
               remNzos++;
               continue;
            }
         }
      }
#if 0 // a lot of this is allready done in redundantRows.

      // Column singleton without objective
      else if (col.size() == 1 && isZero(lp.maxObj(i), epsilonSimplifier()))
      {
         Real x = col.value(0);
         int  j = col.index(0);
         
         Real           up;
         Real           lo;


         if (GT(x, 0.0, epsilonSimplifier()))
         {
            if (lp.lower(i) > -infinity)
               // TODO here should be checked if rhs < infinity
               up = lp.rhs(j) - lp.lower(i) * x;
            else
               up = infinity;
            
            if (lp.upper(i) < infinity)
               lo = lp.lhs(j) - lp.upper(i) * x;
            else
               lo = -infinity;
         }
         else if (LT(x, 0.0, epsilonSimplifier()))
         {
            if (lp.lower(i) > -infinity)
               lo = lp.lhs(j) - lp.lower(i) * x;
            else
               lo = -infinity;
            
            if (lp.upper(i) < infinity)
               up = lp.rhs(j) - lp.upper(i) * x;
            else
               up = infinity;
         }
         else
         {
            up = lp.rhs(j);
            lo = lp.lhs(j);
         }
         
         if (isZero(lo, epsilonSimplifier()))
            lo = 0.0;
         if (isZero(up, epsilonSimplifier()))
            up = 0.0;
         
         lp.changeRange(j, lo, up);
         rem[i] = -1;
         num++;
      }
#endif // 0
   }
   assert(remCols != 0 || remNzos == 0);

   if (remCols > 0)
   {
      removeCols(lp, rem, remCols);

      again      = true;
      m_remCols += remCols;
      m_remNzos += remNzos;

      MSG_INFO2( spxout << "IREDSM24 simple col simplifier removed "
                           << remCols << " cols, "
                           << remNzos << " nzos"
                           << std::endl; )
   }
   return OKAY;
}

void SPxRedundantSM::removeRows(SPxLP& lp, DataArray<int>& rem, int num) 
{
   int i;

   assert(rem.size() == lp.nRows());

   lp.removeRows(rem.get_ptr());

   DataArray<int> tmp(m_rperm);
      
   for(i = 0; i < rem.size(); ++i)
      if (rem[i] >= 0)
         m_rperm[rem[i]] = tmp[i];
   
   for(i = lp.nRows(); i < lp.nRows() + num; ++i)
      m_rperm[i] = -1;
   
   assert(lp.isConsistent());
}

void SPxRedundantSM::removeCols(SPxLP& lp, DataArray<int>& rem, int num) 
{
   int i;

   assert(rem.size() == lp.nCols());

   lp.removeCols(rem.get_ptr());
      
   DataArray<int> tmp(m_cperm);
      
   for(i = 0; i < rem.size(); ++i)
      if (rem[i] >= 0)
         m_cperm[rem[i]] = tmp[i];
   
   for(i = lp.nCols(); i < lp.nCols() + num; ++i)
      m_cperm[i] = -1;
#if 0
     for(i = 0; i < cperm.size(); ++i)
      std::cout << i << " -> " << cperm[i] << " " << tmp[i] << " " << rem[i] << "\n";
#endif
   assert(lp.isConsistent());
}

SPxSimplifier::Result SPxRedundantSM::simplify(SPxLP& lp, Real eps, Real delta)
{
   METHOD("SPxRedundantSM::simplify()");

   Result ret = OKAY;
   bool   again;
   bool   rcagain = false;
   bool   rragain = false;
   int    i;
   
   m_timeUsed.reset();
   m_timeUsed.start();

   m_epsilon = eps;
   m_delta   = delta;

   m_prim.reDim(lp.nCols());
   m_dual.reDim(lp.nRows());

   m_cperm.reSize(lp.nCols());
   m_rperm.reSize(lp.nRows());

   for(i = 0; i < lp.nRows(); ++i)
      m_rperm[i] = i;

   for(i = 0; i < lp.nCols(); ++i)
      m_cperm[i] = i;

   //std::cout << lp << std::endl;

   do
   {
      again = false;
 
      if (ret == OKAY)
         ret = simpleRows(lp, again);

      if (ret == OKAY)
         ret = simpleCols(lp, again);

      assert(ret == OKAY || !again);

      //std::cout << lp << std::endl;
   }
   while(again);

   //std::cout << lp << std::endl;
#if 1
   if (ret == OKAY)
      ret = redundantCols(lp, rcagain);

   if (ret == OKAY)
      ret = redundantRows(lp, rragain);

   again = (ret == OKAY) && (rcagain || rragain);

   // This has to be a loop, otherwise we could end up with
   // empty rows.
   while(again)
   {
      again = false;

      if (ret == OKAY)
         ret = simpleRows(lp, again);

      if (ret == OKAY)
         ret = simpleCols(lp, again);

      assert(ret == OKAY || !again);
   }
#endif
   MSG_INFO1( spxout << "IREDSM25 redundant simplifier removed "
                        << m_remRows << " rows, "
                        << m_remNzos << " nzos, changed "
                        << m_chgBnds << " col bounds " 
                        << m_chgLRhs << " row bounds,"
                        << std::endl; )

   if (lp.nCols() == 0 && lp.nRows() == 0)
   {
      MSG_INFO1( spxout << "IREDSM26 simplifier removed all rows and columns" << std::endl; )
      ret = VANISHED;
   }
   m_timeUsed.stop();

   return ret;
}

const Vector& SPxRedundantSM::unsimplifiedPrimal(const Vector& x)
{
   assert(x.dim()                 <= m_cperm.size());
   assert(x.dim() + m_pval.size() == m_prim.dim());

   int i;

   for(i = 0; i < x.dim(); ++i)
      m_prim[m_cperm[i]] = x[i];

   assert(m_cperm.size() == x.dim() || m_cperm[i] < 0);

   for(i = 0; i < m_pval.size(); ++i)
      m_prim[m_pval.index(i)] = m_pval.value(i);

   return m_prim;
}

const Vector& SPxRedundantSM::unsimplifiedDual(const Vector& pi)
{
   MSG_ERROR( spxout << "EREDSM27 SPxRedundantSM::unsimplifiedDual() not implemented\n"; )

   assert(false);

   return pi;
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
