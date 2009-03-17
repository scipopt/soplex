/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the class library                   */
/*       SoPlex --- the Sequential object-oriented simPlex.                  */
/*                                                                           */
/*    Copyright (C) 1997-1999 Roland Wunderling                              */
/*                  1997-2009 Konrad-Zuse-Zentrum                            */
/*                            fuer Informationstechnik Berlin                */
/*                                                                           */
/*  SoPlex is distributed under the terms of the ZIB Academic Licence.       */
/*                                                                           */
/*  You should have received a copy of the ZIB Academic License              */
/*  along with SoPlex; see the file COPYING. If not email to soplex@zib.de.  */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
#pragma ident "@(#) $Id: spxmainsm.cpp,v 1.20 2009/03/17 23:25:26 bzfgleix Exp $"

//#define DEBUGGING 1

#include <iostream>

#include "spxmainsm.h"
#include "array.h"
#include "dataarray.h"
#include "sorter.h"
#include "spxout.h"
#include <sstream>
#include <iostream>
#include <fstream>

#define FREE_BOUNDS             1
#define FREE_LHS_RHS            1
#define FREE_CONSTRAINT         1
#define EMPTY_CONSTRAINT        1
#define ROW_SINGLETON           1
#define FORCE_CONSTRAINT        1
#define EMPTY_COLUMN            1
#define FREE_ZERO_OBJ_VARIABLE  1
#define ZERO_OBJ_COL_SINGLETON  1
#define FREE_COL_SINGLETON      1
#define DOUBLETON_EQUATION      1
#define DOMINATED_COLUMN        1
#define WEAKLY_DOMINATED_COLUMN 1
#define FIX_VARIABLE            1

#define EXTREMES                1
#define ROWS                    1 
#define COLS                    1
#define DUAL                    1
#define DUPLICATE_ROWS          1
#define DUPLICATE_COLS          1

namespace soplex
{
void SPxMainSM::FreeConstraintPS::execute(DVector& x, DVector& y, DVector& s, DVector&,
                                          DataArray<SPxSolver::VarStatus>&,
                                          DataArray<SPxSolver::VarStatus>& rStatus) const
{
   // primal:
   Real slack = 0.0;

   for (int k = 0; k < m_row.size(); ++k)
      slack += m_row.value(k) * x[m_row.index(k)];
   
   s[m_i] = slack;

   // dual:
   y[m_i] = 0.0;
   
   // basis:
   rStatus[m_i] = SPxSolver::BASIC;
}

void SPxMainSM::EmptyConstraintPS::execute(DVector&, DVector& y, DVector& s, DVector&,
                                           DataArray<SPxSolver::VarStatus>&,
                                           DataArray<SPxSolver::VarStatus>& rStatus) const
{
   // primal: 
   s[m_i] = 0.0;

   // dual:
   y[m_i] = 0.0;
   
   // basis:
   rStatus[m_i] = SPxSolver::BASIC;
}
   
void SPxMainSM::RowSingletonPS::execute(DVector& x, DVector& y, DVector& s, DVector& r,
                                        DataArray<SPxSolver::VarStatus>& cStatus,
                                        DataArray<SPxSolver::VarStatus>& rStatus) const
{
   Real aij = m_col[m_i];
   
   // primal:
   s[m_i] = aij * x[m_j];
   
   // dual & basis:
   Real val = m_obj;
   
   for(int k = 0; k < m_col.size(); ++k)
      if (m_col.index(k) != m_i) 
         val -= m_col.value(k) * y[m_col.index(k)];
  
   if (cStatus[m_j] == SPxSolver::BASIC)
   {
      assert(isZero(r[m_j], eps()));
      
      y[m_i] = val / aij;
      r[m_j] = 0.0;
      
      rStatus[m_i] = SPxSolver::BASIC;
   }
   else if (cStatus[m_j] == SPxSolver::ON_LOWER && m_strictLo)
   {
      y[m_i] = val / aij;
      r[m_j] = 0.0;

      cStatus[m_j] = SPxSolver::BASIC;
      rStatus[m_i] = (aij > 0) ? SPxSolver::ON_LOWER : SPxSolver::ON_UPPER;
   }
   else if (cStatus[m_j] == SPxSolver::ON_UPPER && m_strictUp)
   {
      y[m_i] = val / aij;
      r[m_j] = 0.0;

      cStatus[m_j] = SPxSolver::BASIC;
      rStatus[m_i] = (aij > 0) ? SPxSolver::ON_UPPER : SPxSolver::ON_LOWER;
   }
   else if (cStatus[m_j] == SPxSolver::FIXED &&
            (( m_maxSense && ((r[m_j] > 0 && m_strictUp) || (r[m_j] < 0 && m_strictLo))) ||
             (!m_maxSense && ((r[m_j] > 0 && m_strictLo) || (r[m_j] < 0 && m_strictUp)))))
   {
      y[m_i] = val / aij;
      r[m_j] = 0.0;
      
      cStatus[m_j] = SPxSolver::BASIC;
      if (m_strictLo && m_strictUp)
         rStatus[m_i] = SPxSolver::FIXED;
      else if (m_strictLo)
         rStatus[m_i] = (aij > 0) ? SPxSolver::ON_LOWER : SPxSolver::ON_UPPER;
      else if (m_strictUp)
         rStatus[m_i] = (aij > 0) ? SPxSolver::ON_UPPER : SPxSolver::ON_LOWER;
   }
   else
   {
      y[m_i] = 0.0;
      r[m_j] = val;
       
      rStatus[m_i] = SPxSolver::BASIC;

      if (cStatus[m_j] == SPxSolver::FIXED && (!m_strictLo || !m_strictUp))
      {
         if (m_strictLo)
            rStatus[m_i] = (aij > 0) ? SPxSolver::ON_LOWER : SPxSolver::ON_UPPER;
         else 
            rStatus[m_i] = (aij > 0) ? SPxSolver::ON_UPPER : SPxSolver::ON_LOWER;
      }
   }
}

void SPxMainSM::ForceConstraintPS::execute(DVector&, DVector& y, DVector& s, DVector& r,
                                           DataArray<SPxSolver::VarStatus>& cStatus,
                                           DataArray<SPxSolver::VarStatus>& rStatus) const
{
   // primal:
   s[m_i] = m_lRhs;

   // dual:
   Real lo = -infinity;
   Real up =  infinity;

   for(int k = 0; k < m_row.size(); ++k)
   {
      const SVector& col = m_cols[k];
      Real aij           = m_row.value(k);
      Real bound         = m_objs[m_row.index(k)];
      
      ASSERT_WARN( "WMAISM72", isNotZero(aij) );
             
      for(int l = 0; l < col.size(); ++l)
         if (col.index(l) != m_i)
            bound -= col.value(l) * y[col.index(l)]; 
      
      bound /= aij;

      if (m_maxSense)
      {
         if (m_lhsFixed && bound < up)
            up = bound; 
         else if (!m_lhsFixed && bound > lo)
            lo = bound;             
      }
      else
      {
         if (m_lhsFixed && bound > lo)
            lo = bound;
         else if (!m_lhsFixed && bound < up)
            up = bound;
      }
   }
   
   assert(LE(lo, up));
  
   if (lo > -infinity)
      y[m_i] = lo;
   else if (up < infinity)
      y[m_i] = up;
   else
      y[m_i] = 0.0;

   if (isZero(y[m_i], eps()))
      y[m_i] = 0.0;
   
   for(int k = 0; k < m_row.size(); ++k)
      r[m_row.index(k)] -= m_row.value(k) * y[m_i];

   // basis:
   rStatus[m_i] = SPxSolver::BASIC;
   for(int k = 0; k < m_row.size(); ++k)
   {
      if (m_fixed[k])
         cStatus[m_row.index(k)] = SPxSolver::FIXED;
      else if (m_lhsFixed)
         cStatus[m_row.index(k)] = (m_row.value(k) > 0) ? SPxSolver::ON_LOWER : SPxSolver::ON_UPPER;
      else
         cStatus[m_row.index(k)] = (m_row.value(k) > 0) ? SPxSolver::ON_UPPER : SPxSolver::ON_LOWER;
   }
}

void SPxMainSM::FixVariablePS::execute(DVector& x, DVector& y, DVector& s, DVector& r,
                                       DataArray<SPxSolver::VarStatus>& cStatus,
                                       DataArray<SPxSolver::VarStatus>&) const
{
   // primal:
   x[m_j] = m_val;
   
   for(int k = 0; k < m_col.size(); ++k)
      s[m_col.index(k)] += m_col.value(k) * x[m_j];
   
   // dual:
   Real val = m_obj;

   for(int k = 0; k < m_col.size(); ++k)
      val -= m_col.value(k) * y[m_col.index(k)];

   r[m_j] = val;

   // basis:
   cStatus[m_j] = SPxSolver::FIXED;
}

void SPxMainSM::FixBoundsPS::execute(DVector&, DVector&, DVector&, DVector&,
                                     DataArray<SPxSolver::VarStatus>& cStatus,
                                     DataArray<SPxSolver::VarStatus>&) const
{
   // basis:
   cStatus[m_j] = m_status;
}

void SPxMainSM::FreeZeroObjVariablePS::execute(DVector& x, DVector& y, DVector& s, DVector& r,
                                               DataArray<SPxSolver::VarStatus>& cStatus,
                                               DataArray<SPxSolver::VarStatus>& rStatus) const
{
   // primal:
   int      domIdx = -1;
   DSVector slack(m_col.size());
   
   if (m_loFree)
   {
      Real minRowUp = infinity;
     
      for(int k = 0; k < m_rows.size(); ++k)
      {
	 Real           val = 0.0;
	 const SVector& row = m_rows[k];       
         
         for(int l = 0; l < row.size(); ++l)
            if (row.index(l) != m_j)
	       val += row.value(l) * x[row.index(l)];
        
         Real scale = maxAbs(m_lRhs[k], val);  
         
         if (scale < 1.0)
            scale = 1.0;
         
         Real z = (m_lRhs[k] / scale) - (val / scale);

         if (isZero(z))
            z = 0.0;
         
         Real up = z * scale / row[m_j];
	 slack.add(k, val);
         
         if (up < minRowUp)
         {
            minRowUp = up; 
            domIdx   = k;
         }
      }
      
      if (m_bnd < minRowUp)
      {
         x[m_j] = m_bnd;
         domIdx = -1;
      }
      else
         x[m_j] = minRowUp;
   }
   else
   {
      Real maxRowLo = -infinity;
     
      for(int k = 0; k < m_rows.size(); ++k)
      {
	 Real           val = 0.0;
         const SVector& row = m_rows[k]; 
         
         for(int l = 0; l < row.size(); ++l)
            if (row.index(l) != m_j)
               val += row.value(l) * x[row.index(l)];
         
         Real scale = maxAbs(m_lRhs[k], val);  
         
         if (scale < 1.0)
            scale = 1.0;
         
         Real z = (m_lRhs[k] / scale) - (val / scale);

         if (isZero(z))
            z = 0.0;
         
         Real lo = z * scale / row[m_j];
	 slack.add(k, val);
	 
         if (lo > maxRowLo)
         {
            maxRowLo = lo; 
            domIdx   = k;
         }
      }
      
      if (m_bnd > maxRowLo)
      {
         x[m_j] = m_bnd;
         domIdx = -1;
      }
      else
         x[m_j] = maxRowLo;
   }

   for(int k = 0; k < m_col.size(); ++k)
      s[m_col.index(k)] = slack[k] + m_col.value(k) * x[m_j];
   
   // dual:
   r[m_j] = 0.0;

   for(int k = 0; k < m_col.size(); ++k)
      y[m_col.index(k)] = 0.0;
   
   // basis:
   for(int k = 0; k < m_col.size(); ++k)
   {
      if (k != domIdx)         
         rStatus[m_col.index(k)] = SPxSolver::BASIC;
      
      else
      {
	 cStatus[m_j] = SPxSolver::BASIC;
         if (m_loFree)
            rStatus[m_col.index(k)] = (m_col.value(k) > 0) ? SPxSolver::ON_UPPER : SPxSolver::ON_LOWER;
         else 
            rStatus[m_col.index(k)] = (m_col.value(k) > 0) ? SPxSolver::ON_LOWER : SPxSolver::ON_UPPER;
      }
   }  
   if (domIdx == -1)
   {
      if (m_loFree)
         cStatus[m_j] = SPxSolver::ON_UPPER;
      else
         cStatus[m_j] = SPxSolver::ON_LOWER;
   }
}

void SPxMainSM::ZeroObjColSingletonPS::execute(DVector& x, DVector& y, DVector& s, DVector& r,
                                               DataArray<SPxSolver::VarStatus>& cStatus,
                                               DataArray<SPxSolver::VarStatus>& rStatus) const
{
   // primal & basis:
   Real aij = m_row[m_j];

   if (isZero(s[m_i], 1e-6))
      s[m_i] = 0.0;
   
   Real scale1 = maxAbs(m_lhs, s[m_i]);
   Real scale2 = maxAbs(m_rhs, s[m_i]);
   
   if (scale1 < 1.0)
      scale1 = 1.0;
   if (scale2 < 1.0)
      scale2 = 1.0;
            
   Real z1 = (m_lhs / scale1) - (s[m_i] / scale1);
   Real z2 = (m_rhs / scale2) - (s[m_i] / scale2);

   if (isZero(z1))
      z1 = 0.0;
   if (isZero(z2))
      z2 = 0.0;
         
   Real lo = (aij > 0) ? z1 * scale1 / aij : z2 * scale2 / aij;
   Real up = (aij > 0) ? z2 * scale2 / aij : z1 * scale1 / aij;
   
   if (isZero(lo, eps()))
      lo = 0.0;
   if (isZero(up, eps()))
      up = 0.0;   
  
   assert(LErel(lo, up));
   assert((LErel(m_lower, lo, 1e-10) && LErel(lo, m_upper, 1e-10)) ||
          (LErel(m_lower, up, 1e-10) && LErel(up, m_upper, 1e-10)) ||
          (LErel(lo, m_lower, 1e-10) && LErel(m_upper, up, 1e-10)));
   ASSERT_WARN( "WMAISM01", isNotZero(aij) );
   
   if (rStatus[m_i] == SPxSolver::ON_LOWER)
   {
      if (EQrel(m_lower, m_upper))
      {
         x[m_j]       = m_lower;
         cStatus[m_j] = SPxSolver::FIXED;
      }
      else if (aij > 0)
      {
         x[m_j]       = m_upper;
         cStatus[m_j] = SPxSolver::ON_UPPER;
      }
      else if (aij < 0)
      {
         x[m_j]       = m_lower;
         cStatus[m_j] = SPxSolver::ON_LOWER;
      }
      else
         throw SPxInternalCodeException("XMAISM01 This should never happen.");
   }
   else if (rStatus[m_i] == SPxSolver::ON_UPPER)
   {
      if (EQrel(m_lower, m_upper))
      {
         x[m_j]       = m_lower;
         cStatus[m_j] = SPxSolver::FIXED;
      }
      else if (aij > 0)
      {
         x[m_j]       = m_lower;
         cStatus[m_j] = SPxSolver::ON_LOWER;
      }
      else if (aij < 0)
      {
         x[m_j]       = m_upper;
         cStatus[m_j] = SPxSolver::ON_UPPER;
      }
      else
         throw SPxInternalCodeException("XMAISM02 This should never happen.");
   }
   else if (rStatus[m_i] == SPxSolver::FIXED)
   {
      assert(EQrel(m_lower, m_upper));
 
      x[m_j]        = m_lower;
      cStatus[m_j]  = SPxSolver::FIXED;
   }
   else if (rStatus[m_i] == SPxSolver::BASIC)
   {
      if (GErel(m_lower, lo, eps()) && LErel(m_lower, up, eps()))
      {
         x[m_j]       = m_lower;
         cStatus[m_j] = EQrel(m_lower, m_upper) ? SPxSolver::FIXED : SPxSolver::ON_LOWER;
      }
      else if (GErel(m_upper, lo, eps()) && LErel(m_upper, up, eps()))
      {
         x[m_j]       = m_upper;
         cStatus[m_j] = EQrel(m_lower, m_upper) ? SPxSolver::FIXED : SPxSolver::ON_UPPER;
      }
      else
      {
         // make m_i non-basic and m_j basic
         if (EQrel(lo, up))
         {
            x[m_j]       = lo;
            cStatus[m_j] = SPxSolver::BASIC;
            rStatus[m_i] = SPxSolver::FIXED;
         }
         else if (GE(lo, m_lower, eps()) && LE(lo, m_upper, eps()))
         {
            x[m_j]       = lo;
            cStatus[m_j] = SPxSolver::BASIC;
            rStatus[m_i] = (aij > 0) ? SPxSolver::ON_LOWER : SPxSolver::ON_UPPER;
         }
         else if (GE(up, m_lower, eps()) && LE(up, m_upper, eps()))
         {
            x[m_j]       = up;
            cStatus[m_j] = SPxSolver::BASIC;
            rStatus[m_i] = (aij > 0) ? SPxSolver::ON_UPPER : SPxSolver::ON_LOWER;
         }
         else
            throw SPxInternalCodeException("XMAISM03 This should never happen.");
      }
   }
   else
      throw SPxInternalCodeException("XMAISM04 This should never happen.");

   s[m_i] += aij * x[m_j];

   // dual:
   r[m_j] = -1.0 * aij * y[m_i];

   assert(cStatus[m_j] != SPxSolver::BASIC || isZero(r[m_j], eps()));
}

void SPxMainSM::FreeColSingletonPS::execute(DVector& x, DVector& y, DVector& s, DVector& r,
                                            DataArray<SPxSolver::VarStatus>& cStatus,
                                            DataArray<SPxSolver::VarStatus>& rStatus) const
{
   // primal:
   Real val = 0.0;
   Real aij = m_row[m_j];

   for(int k = 0; k < m_row.size(); ++k)
      if (m_row.index(k) != m_j)
 	 val += m_row.value(k) * x[m_row.index(k)];
  
   Real scale = maxAbs(m_lRhs, val);  
   
   if (scale < 1.0)
      scale = 1.0;
   
   Real z = (m_lRhs / scale) - (val / scale);

   if (isZero(z))
      z = 0.0;
 
   x[m_j] = z * scale / aij;
   s[m_i] = m_lRhs;
   
   // dual:
   y[m_i] = m_obj / aij;
   r[m_j] = 0.0;
   
   // basis:
   cStatus[m_j] = SPxSolver::BASIC;
      
   if (m_eqCons)
      rStatus[m_i] = SPxSolver::FIXED;
   else if (m_onLhs)
      rStatus[m_i] = SPxSolver::ON_LOWER;
   else
      rStatus[m_i] = SPxSolver::ON_UPPER;
}

void SPxMainSM::DoubletonEquationPS::execute(DVector&, DVector& y, DVector&, DVector& r,
                                             DataArray<SPxSolver::VarStatus>& cStatus,
                                             DataArray<SPxSolver::VarStatus>&) const
{
   // dual: 
   if ((cStatus[m_k]  != SPxSolver::BASIC) &&
       ((cStatus[m_k] == SPxSolver::ON_LOWER && m_strictLo) ||
        (cStatus[m_k] == SPxSolver::ON_UPPER && m_strictUp) ||
        (cStatus[m_k] == SPxSolver::FIXED    && 
         (( m_maxSense && ((r[m_j] > 0 && m_strictUp) || (r[m_j] < 0 && m_strictLo))) ||
          (!m_maxSense && ((r[m_j] > 0 && m_strictLo) || (r[m_j] < 0 && m_strictUp)))))))
   {
      Real val  = m_kObj;
      Real aik  = m_col[m_i];
      
      for(int k = 0; k < m_col.size(); ++k)
         if (m_col.index(k) != m_i)
            val -= m_col.value(k) * y[m_col.index(k)];
      
      y[m_i] = val / aik;
      r[m_k] = 0.0;

      r[m_j] = m_jObj - val * m_aij / aik;

      ASSERT_WARN( "WMAISM73", isNotZero(m_aij * aik) );
     
      // basis:
      if (cStatus[m_k] == SPxSolver::ON_LOWER)
      {
         if (m_jFixed)
            cStatus[m_j] = SPxSolver::FIXED;
         else 
            cStatus[m_j] = (m_aij * aik > 0) ? SPxSolver::ON_UPPER : SPxSolver::ON_LOWER;
      }
      else if (cStatus[m_k] == SPxSolver::ON_UPPER)
      {
         if (m_jFixed)
            cStatus[m_j] = SPxSolver::FIXED;
         else 
            cStatus[m_j] = (m_aij * aik > 0) ? SPxSolver::ON_LOWER : SPxSolver::ON_UPPER;
      }
      else
      {
         if (m_jFixed)
            cStatus[m_j] = SPxSolver::FIXED;
         else if (m_strictLo)
            cStatus[m_j] = (m_aij * aik > 0) ? SPxSolver::ON_UPPER : SPxSolver::ON_LOWER;
         else
            cStatus[m_j] = (m_aij * aik > 0) ? SPxSolver::ON_LOWER : SPxSolver::ON_UPPER;
      }
      
      cStatus[m_k] = SPxSolver::BASIC;
   }
}

void SPxMainSM::DuplicateRowsPS::execute(DVector&, DVector& y, DVector& s, DVector&,
                                         DataArray<SPxSolver::VarStatus>&,
                                         DataArray<SPxSolver::VarStatus>& rStatus) const
{
   // primal:
   for(int k = 0; k < m_scale.size(); ++k)
      if (m_scale.index(k) != m_i)
         s[m_scale.index(k)] = s[m_i] / m_scale.value(k); 

   // dual & basis:
   for(int k = 0; k < m_scale.size(); ++k)
   {
      int i = m_scale.index(k);

      ASSERT_WARN( "WMAISM02", isNotZero(m_scale.value(k)) );
      
      if (rStatus[m_i] == SPxSolver::FIXED && i == m_maxLhsIdx && i == m_minRhsIdx)
      {
         y[i]   = y[m_i] * m_scale.value(k);
         y[m_i] = 0.0;
         
         rStatus[i]   = SPxSolver::FIXED;
         rStatus[m_i] = SPxSolver::BASIC;
      }
      else if (i == m_maxLhsIdx &&
               ((rStatus[m_i] == SPxSolver::ON_LOWER) ||
                (rStatus[m_i] == SPxSolver::FIXED    &&
                 ((m_maxSense && y[m_i] < 0) || (!m_maxSense && y[m_i] > 0))))) 
      {
         y[i]   = y[m_i] * m_scale.value(k);
         y[m_i] = 0.0;

         rStatus[i]   = (m_scale.value(k) > 0) ? SPxSolver::ON_LOWER : SPxSolver::ON_UPPER;
         rStatus[m_i] = SPxSolver::BASIC;
      }
      else if (i == m_minRhsIdx &&
               ((rStatus[m_i] == SPxSolver::ON_UPPER) ||
                (rStatus[m_i] == SPxSolver::FIXED    &&
                 ((m_maxSense && y[m_i] > 0) || (!m_maxSense && y[m_i] < 0)))))
      {   
         y[i]   = y[m_i] * m_scale.value(k);
         y[m_i] = 0.0;

         rStatus[i]   = (m_scale.value(k) > 0) ? SPxSolver::ON_UPPER : SPxSolver::ON_LOWER;
         rStatus[m_i] = SPxSolver::BASIC;
      }
      else if (i != m_i)
      {
         y[i]       = 0.0;
         rStatus[i] = SPxSolver::BASIC;
      }
   }   
      
   // nothing to do for the reduced cost values
}

void SPxMainSM::DuplicateColsPS::execute(DVector& x, 
                                         DVector&, 
                                         DVector&, 
                                         DVector& r,
                                         DataArray<SPxSolver::VarStatus>& cStatus,
                                         DataArray<SPxSolver::VarStatus>&) const
{
   // primal & basis:
   ASSERT_WARN( "WMAISM03", isNotZero(m_scale) );

   if (cStatus[m_k] == SPxSolver::ON_LOWER)
   {
      x[m_k] = m_loK;
         
      if (m_scale > 0)
      {
         x[m_j]       = m_loJ;
         cStatus[m_j] = EQrel(m_loJ, m_upJ) ? SPxSolver::FIXED : SPxSolver::ON_LOWER;
      }
      else
      {
         x[m_j]       = m_upJ;
         cStatus[m_j] = EQrel(m_loJ, m_upJ) ? SPxSolver::FIXED : SPxSolver::ON_UPPER;
      }
   }
   else if (cStatus[m_k] == SPxSolver::ON_UPPER)
   {
      x[m_k] = m_upK;
         
      if (m_scale > 0)
      {
         x[m_j]       = m_upJ;
         cStatus[m_j] = EQrel(m_loJ, m_upJ) ? SPxSolver::FIXED : SPxSolver::ON_UPPER;
      }
      else
      {
         x[m_j]       = m_loJ;
         cStatus[m_j] = EQrel(m_loJ, m_upJ) ? SPxSolver::FIXED : SPxSolver::ON_LOWER;
      }
   }
   else if (cStatus[m_k] == SPxSolver::FIXED)
   {
      // => x[m_k] and x[m_j] are also fixed before the corresponding preprocessing step
      x[m_j]       = m_loJ;
      cStatus[m_j] = SPxSolver::FIXED;
   }
   else if (cStatus[m_k] == SPxSolver::ZERO)
   {
      if (m_loK <= -infinity && m_upK >= infinity)
      {
         x[m_j] = 0.0;
         
         if (isZero(m_loJ) && isZero(m_upJ))
            cStatus[m_j] = SPxSolver::FIXED;
         else if (isZero(m_loJ))
            cStatus[m_j] = SPxSolver::ON_LOWER;
         else if (isZero(m_upJ))
            cStatus[m_j] = SPxSolver::ON_UPPER;
         else if (m_loJ <= -infinity && m_upJ >= infinity)
            cStatus[m_j] = SPxSolver::ZERO;
         else
            throw SPxInternalCodeException("XMAISM05 This should never happen.");
      }
      else if (m_loJ <= -infinity && m_upJ >= infinity)
      {
         x[m_j]       = 0.0;
         cStatus[m_j] = SPxSolver::ZERO;
         
         if (isZero(m_loK) && isZero(m_upK))
            cStatus[m_k] = SPxSolver::FIXED;
         else if (isZero(m_loK))
            cStatus[m_k] = SPxSolver::ON_LOWER;
         else if (isZero(m_upK))
            cStatus[m_k] = SPxSolver::ON_UPPER;
         else if (m_loK <= -infinity && m_upK >= infinity)
            cStatus[m_k] = SPxSolver::ZERO;
         else
            throw SPxInternalCodeException("XMAISM06 This should never happen.");
      }
      else if (m_loK > -infinity)
      {
         x[m_j]       = x[m_k] - m_scale * m_loK;
         x[m_k]       = m_loK;
         cStatus[m_k] = SPxSolver::ON_LOWER;
         
         if (EQrel(x[m_j], m_loJ))
            cStatus[m_j] = EQrel(m_loJ, m_upJ) ? SPxSolver::FIXED : SPxSolver::ON_LOWER;
         else if (EQrel(x[m_j], m_upJ))
            cStatus[m_j] = EQrel(m_loJ, m_upJ) ? SPxSolver::FIXED : SPxSolver::ON_UPPER;
         else
            throw SPxInternalCodeException("XMAISM07 This should never happen.");
      }
      else if (m_upK < infinity)
      {
         x[m_j]       = x[m_k] - m_scale * m_upK;
         x[m_k]       = m_upK;
         cStatus[m_k] = SPxSolver::ON_UPPER;
         
         if (EQrel(x[m_j], m_loJ))
            cStatus[m_j] = EQrel(m_loJ, m_upJ) ? SPxSolver::FIXED : SPxSolver::ON_LOWER;
         else if (EQrel(x[m_j], m_upJ))
            cStatus[m_j] = EQrel(m_loJ, m_upJ) ? SPxSolver::FIXED : SPxSolver::ON_UPPER;
         else
            throw SPxInternalCodeException("XMAISM08 This should never happen.");
      }
      else
         throw SPxInternalCodeException("XMAISM09 This should never happen.");
   }
   else if (cStatus[m_k] == SPxSolver::BASIC)
   {
      Real scale1 = maxAbs(x[m_k], m_loK);  
      Real scale2 = maxAbs(x[m_k], m_upK);  
      
      if (scale1 < 1.0)
         scale1 = 1.0;
      if (scale2 < 1.0)
         scale2 = 1.0;
         
      Real z1 = (x[m_k] / scale1) - (m_loK / scale1);
      Real z2 = (x[m_k] / scale2) - (m_upK / scale2);

      if (isZero(z1))
         z1 = 0.0;
      if (isZero(z2))
         z2 = 0.0;
 
      if (m_loJ <= -infinity && m_upJ >= infinity && m_loK <= -infinity && m_upK >= infinity)
      {
         x[m_j]       = 0.0;
         cStatus[m_j] = SPxSolver::ZERO;
      }
      else if (m_loJ > -infinity                             &&
               GErel(x[m_k] - m_scale * m_loJ, m_loK, eps()) &&
               LErel(x[m_k] - m_scale * m_loJ, m_upK, eps()))
      {
         x[m_j]       = m_loJ;
         x[m_k]      -= m_scale * x[m_j];         
         cStatus[m_j] = EQrel(m_loJ, m_upJ) ? SPxSolver::FIXED : SPxSolver::ON_LOWER;
      }
      else if (m_upJ < infinity                              &&
               GErel(x[m_k] - m_scale * m_upJ, m_loK, eps()) &&
               LErel(x[m_k] - m_scale * m_upJ, m_upK, eps()))
      {
         x[m_j]       = m_upJ;
         x[m_k]      -= m_scale * x[m_j];   
         cStatus[m_j] = EQrel(m_loJ, m_upJ) ? SPxSolver::FIXED : SPxSolver::ON_UPPER;
      }
      else if (m_loK > -infinity                          &&
               GErel(z1 * scale1 / m_scale, m_loJ, eps()) && 
               LErel(z1 * scale1 / m_scale, m_upJ, eps()))
      {
         x[m_j]       = z1 * scale1 / m_scale;         
         x[m_k]       = m_loK;
         cStatus[m_j] = SPxSolver::BASIC;
         cStatus[m_k] = EQrel(m_loK, m_upK) ? SPxSolver::FIXED : SPxSolver::ON_LOWER;
      }
      else if (m_upK < infinity                           &&
               GErel(z2 * scale2 / m_scale, m_loJ, eps()) &&  
               LErel(z2 * scale2 / m_scale, m_upJ, eps()))
      {
         x[m_j]       = z2 * scale2 / m_scale;                
         x[m_k]       = m_upK;
         cStatus[m_j] = SPxSolver::BASIC;
         cStatus[m_k] = EQrel(m_loK, m_upK) ? SPxSolver::FIXED : SPxSolver::ON_UPPER;
      }
      else
         throw SPxInternalCodeException("XMAISM10 This should never happen.");
   }

   // dual:
   r[m_j] = m_scale * r[m_k];
}

void SPxMainSM::handleExtremes(SPxLP& lp)
{
   METHOD( "SPxMainSM::handleExtremes" );

   // This method handles extreme value of the given LP by
   //
   // 1. setting numbers of very small absolute values to zero and
   // 2. setting numbers of very large absolute values to -infinity or +infinity, respectively.

   Real maxVal  = infinity / 5.0;
   int  remRows = 0;
   int  remNzos = 0;
   int  chgBnds = 0;
   int  chgLRhs = 0;
   int  objCnt  = 0;
   
   for(int i = lp.nRows()-1; i >= 0; --i)
   {
      // lhs
      Real lhs = lp.lhs(i);
 
      if (lhs != 0.0 && isZero(lhs, epsZero()))
      {
         lp.changeLhs(i, 0.0);
         ++chgLRhs;
      }
      else if (lhs > -infinity && lhs < -maxVal)
      {
	 lp.changeLhs(i, -infinity);
         ++chgLRhs;
      }
      else if (lhs <  infinity && lhs >  maxVal)
      {
         lp.changeLhs(i,  infinity);
         ++chgLRhs;
      }

      // rhs
      Real rhs = lp.rhs(i);

      if (rhs != 0.0 && isZero(rhs, epsZero()))
      {
         lp.changeRhs(i, 0.0);
         ++chgLRhs;
      }
      else if (rhs > -infinity && rhs < -maxVal)
      {
         lp.changeRhs(i, -infinity);
         ++chgLRhs;
      }
      else if (rhs <  infinity && rhs >  maxVal)
      {
         lp.changeRhs(i,  infinity);
         ++chgLRhs;
      }

      if (lp.lhs(i) <= -infinity && lp.rhs(i) >= infinity)
      {
         m_hist.append(new FreeConstraintPS(lp, *this, i));

         removeRow(lp, i);
         ++remRows;
         
         ++m_stat[FREE_ROW];
      }
   }
   
   for(int j = 0; j < lp.nCols(); ++j)
   {
      // lower bound
      Real lo = lp.lower(j);

      if (lo != 0.0 && isZero(lo, epsZero()))
      {
         lp.changeLower(j, 0.0);
         ++chgBnds;
      }
      else if (lo > -infinity && lo < -maxVal)
      {
         lp.changeLower(j, -infinity);
         ++chgBnds;
      }
      else if (lo <  infinity && lo >  maxVal)
      {
         lp.changeLower(j,  infinity);
         ++chgBnds;
      }

      // upper bound
      Real up = lp.upper(j);

      if (up != 0.0 && isZero(up, epsZero()))
      {
         lp.changeUpper(j, 0.0);
         ++chgBnds;
      }
      else if (up > -infinity && up < -maxVal)
      {
         lp.changeUpper(j, -infinity);
         ++chgBnds;
      }
      else if (up <  infinity && up >  maxVal)
      {
         lp.changeUpper(j,  infinity);
         ++chgBnds;
      }

      // fixed columns will be eliminated later
      if (NE(lo, up))
      {
         lo = fabs(lo);
         up = fabs(up);

         Real absBnd = (lo > up) ? lo : up;

         if (absBnd < 1.0)
            absBnd = 1.0;
         
         // non-zeros
         SVector& col = lp.colVector_w(j);
         int        i = 0;
         
         while(i < col.size())
         {
            Real aij = fabs(col.value(i));
            
            if (isZero(aij, epsZero()) || isZero(aij * absBnd, deltaBnd()))
            {
               SVector& row = lp.rowVector_w(col.index(i));

               // this changes col.size()
               row.remove(row.number(j));
               col.remove(i);           

               MSG_INFO3( spxout << "IMAISM04 aij=" << aij 
                                 << " removed, absBnd=" << absBnd 
                                 << std::endl; )
               ++remNzos;
            }         
            else
            {
               if (aij > maxVal)
                  MSG_WARNING( spxout << "WMAISM05 Warning! Big value " << aij << std::endl; )

               ++i;
            }
         }
      }

      // objective
      Real obj = lp.obj(j);

      if (obj != 0.0 && isZero(obj, epsZero()))
      {
         lp.changeObj(j, 0.0);
         ++objCnt;
      }
      else if (obj > -infinity && obj < -maxVal)
      {
         lp.changeObj(j, -infinity);
         ++objCnt;
      }
      else if (obj <  infinity && obj >  maxVal)
      {
         lp.changeObj(j,  infinity);
         ++objCnt;
      }
   }

   if (remRows + remNzos + chgLRhs + chgBnds + objCnt > 0)
   {
      m_remRows += remRows;
      m_remNzos += remNzos;
      m_chgLRhs += chgLRhs;
      m_chgBnds += chgBnds;

      MSG_INFO2( spxout << "IMAISM06 Main simplifier (extremes) removed "
                        << remRows << " rows, "
                        << remNzos << " non-zeros, "
                        << chgBnds << " col bounds, "
                        << chgLRhs << " row bounds, "
                        << objCnt  << " objective coefficients" << std::endl; )
   }   
   assert(lp.isConsistent());
}

SPxSimplifier::Result SPxMainSM::removeEmpty(SPxLP& lp)
{
   METHOD( "SPxMainSM::removeEmpty" );
   
   // This method removes empty rows and columns from the LP.
  
   int remRows = 0;
   int remCols = 0;
   
   for(int i = lp.nRows()-1; i >= 0; --i)
   {
      const SVector& row = lp.rowVector(i);       
      
      if (row.size() == 0)
      {
         MSG_INFO3( spxout << "IMAISM07 row " << i
                           << ": empty ->"; )

         if (LT(lp.rhs(i), 0.0, deltaBnd()) || GT(lp.lhs(i), 0.0, deltaBnd()))
         {
            MSG_INFO3( spxout << " infeasible lhs=" << lp.lhs(i)
                              << " rhs=" << lp.rhs(i) << std::endl; )
            return INFEASIBLE;
         }         
         MSG_INFO3( spxout << " removed" << std::endl; )

         m_hist.append(new EmptyConstraintPS(rIdx(i)));
         
         ++remRows;
         removeRow(lp, i);

         ++m_stat[EMPTY_ROW];
      }
   }
   
   for(int j = lp.nCols()-1; j >= 0; --j)
   {
      const SVector& col = lp.colVector(j);
            
      if (col.size() == 0)
      {
	 MSG_INFO3( spxout << "IMAISM08 col " << j 
                           << ": empty -> maxObj=" << lp.maxObj(j)
                           << " lower=" << lp.lower(j)
                           << " upper=" << lp.upper(j); )

         Real val;

         if (GT(lp.maxObj(j), 0.0, epsZero()))
         {
            if (lp.upper(j) >= infinity)
            {
               MSG_INFO3( spxout << " unbounded" << std::endl; )
               return UNBOUNDED;
            }
            val = lp.upper(j);
         }
         else if (LT(lp.maxObj(j), 0.0, epsZero()))
         {
            if (lp.lower(j) <= -infinity)
            {
               MSG_INFO3( spxout << " unbounded" << std::endl; )
               return UNBOUNDED;
            }
            val = lp.lower(j);
         }
         else 
         {
            ASSERT_WARN( "WMAISM09", isZero(lp.maxObj(j), epsZero()) );
            // any value within the bounds is ok
            if (lp.lower(j) > -infinity)
               val = lp.lower(j);
            else if (lp.upper(j) < infinity)
               val = lp.upper(j);
            else
               val = 0.0;
         }
         MSG_INFO3( spxout << " removed" << std::endl; )

         m_hist.append(new FixBoundsPS(lp, *this, j, val));  
         m_hist.append(new FixVariablePS(lp, *this, j, val));
            
         ++remCols;
         removeCol(lp, j);

         ++m_stat[EMPTY_COL];
      }
   }
   
   if (remRows + remCols > 0)
   {
      m_remRows += remRows;
      m_remCols += remCols;
   
      MSG_INFO2( spxout << "IMAISM10 Main simplifier (empty rows/colums) removed "
                        << remRows << " rows, "
                        << remCols << " cols"
                        << std::endl; )

   }
   return OKAY;
}
 
SPxSimplifier::Result SPxMainSM::simplifyRows(SPxLP& lp, bool& again)
{
   METHOD( "SPxMainSM::simplifyRows" );
  
   // This method simplifies the rows of the LP.
   //
   // The following operations are done:
   // 1. detect implied free variables
   // 2. detect implied free constraints
   // 3. detect infeasible constraints
   // 4. remove unconstrained constraints
   // 5. remove empty constraints
   // 6. remove row singletons and tighten the corresponding variable bounds if necessary  
   // 7. detect forcing rows and fix the corresponding variables
   
   int remRows = 0;
   int remNzos = 0;
   int chgLRhs = 0;
   int chgBnds = 0;
   
   for(int i = lp.nRows()-1; i >= 0; --i)
   {
      const SVector& row = lp.rowVector(i);   
      
      // compute bounds on constraint value
      int  lhsCnt = 0;
      int  rhsCnt = 0;
      Real lhsBnd = 0.0;
      Real rhsBnd = 0.0;
      
      for(int k = 0; k < row.size(); ++k)
      {
         Real aij = row.value(k);
         int  j   = row.index(k);
         
         ASSERT_WARN( "WMAISM11", isNotZero(aij) );
         
         if (aij > 0.0)
         {
            if (lp.lower(j) <= -infinity)
               ++lhsCnt;
            else
               lhsBnd += aij * lp.lower(j);
            
            if (lp.upper(j) >= infinity)
               ++rhsCnt;
            else
               rhsBnd += aij * lp.upper(j);
         }
         else if (aij < 0.0)
         {
            if (lp.lower(j) <= -infinity)
               ++rhsCnt;
            else
               rhsBnd += aij * lp.lower(j);

            if (lp.upper(j) >= infinity)
               ++lhsCnt;
            else
               lhsBnd += aij * lp.upper(j);
         }
      }
    
#if FREE_BOUNDS
      // 1. detect implied free variables
      if (rhsCnt <= 1 || lhsCnt <= 1)
      {
         for(int k = 0; k < row.size(); ++k)
         {
            Real aij = row.value(k);
            int  j   = row.index(k);
            
            ASSERT_WARN( "WMAISM12", isNotZero(aij) );
                     
            if (aij > 0.0)
            {
               if (lp.lhs(i) > -infinity && lp.lower(j) > -infinity && rhsCnt <= 1 && NErel(lp.lhs(i), rhsBnd, deltaBnd()))
               {
                  Real lo    = -infinity;
                  Real scale = maxAbs(lp.lhs(i), rhsBnd);
                  
                  if (scale < 1.0)
                     scale = 1.0;
                  
                  Real z = (lp.lhs(i) / scale) - (rhsBnd / scale);

                  if (isZero(z, epsZero()))
                     z = 0.0;
                  
                  assert(rhsCnt > 0 || lp.upper(j) < infinity);
                  
                  if (rhsCnt == 0)
                     lo = lp.upper(j) + z * scale / aij;
                  else if (lp.upper(j) >= infinity)
                     lo = z * scale / aij;

                  if (isZero(lo, epsZero()))
                     lo = 0.0;

                  if (GErel(lo, lp.lower(j), deltaBnd()))
                  {
                     MSG_INFO3( spxout << "IMAISM13 row " << i
                                       << ": redundant lower bound on x" << j
                                       << " -> lower=" << lo
                                       << " (" << lp.lower(j)
                                       << ")" << std::endl; )
                 
                     ++lhsCnt;
                     lhsBnd -= aij * lp.lower(j);

                     lp.changeLower(j, -infinity);
                     ++chgBnds;    
                  }
               }
               if (lp.rhs(i) < infinity && lp.upper(j) < infinity && lhsCnt <= 1 && NErel(lp.rhs(i), lhsBnd, deltaBnd()))
               {
                  Real up    = infinity;
                  Real scale = maxAbs(lp.rhs(i), lhsBnd);
                  
                  if (scale < 1.0)
                     scale = 1.0;
                  
                  Real z = (lp.rhs(i) / scale) - (lhsBnd / scale);
                  
                  if (isZero(z, epsZero()))
                     z = 0.0;
 
                  assert(lhsCnt > 0 || lp.lower(j) > -infinity);
                 
                  if (lhsCnt == 0)
                     up = lp.lower(j) + z * scale / aij;
                  else if (lp.lower(j) <= -infinity)
                     up = z * scale / aij;
                 
                  if (isZero(up, epsZero()))
                     up = 0.0;

                  if (LErel(up, lp.upper(j), deltaBnd()))
                  {
                     MSG_INFO3( spxout << "IMAISM14 row " << i
                                       << ": redundant upper bound on x" << j
                                       << " -> upper=" << up
                                       << " (" << lp.upper(j)
                                       << ")" << std::endl; )
                  
                     ++rhsCnt;
                     rhsBnd -= aij * lp.upper(j);

                     lp.changeUpper(j, infinity);
                     ++chgBnds;                  
                  }
               }
            }
            else if (aij < 0.0)
            {
               if (lp.lhs(i) > -infinity && lp.upper(j) < infinity && rhsCnt <= 1 && NErel(lp.lhs(i), rhsBnd, deltaBnd()))
               {
                  Real up    = infinity;
                  Real scale = maxAbs(lp.lhs(i), rhsBnd);
                  
                  if (scale < 1.0)
                     scale = 1.0;
                  
                  Real z = (lp.lhs(i) / scale) - (rhsBnd / scale);
                  
                  if (isZero(z, epsZero()))
                     z = 0.0;
 
                  assert(rhsCnt > 0 || lp.lower(j) > -infinity);
                 
                  if (rhsCnt == 0)
                     up = lp.lower(j) + z * scale / aij;
                  else if (lp.lower(j) <= -infinity)
                     up = z * scale / aij;
                  
                  if (isZero(up, epsZero()))
                     up = 0.0;

                  if (LErel(up, lp.upper(j), deltaBnd()))
                  {
                     MSG_INFO3( spxout << "IMAISM15 row " << i
                                       << ": redundant upper bound on x" << j
                                       << " -> upper=" << up
                                       << " (" << lp.upper(j)
                                       << ")" << std::endl; )
                        
                     ++lhsCnt;
                     lhsBnd -= aij * lp.upper(j);
                     
                     lp.changeUpper(j, infinity);
                     ++chgBnds;                                          
                  }
               }
               if (lp.rhs(i) < infinity && lp.lower(j) > -infinity && lhsCnt <= 1 && NErel(lp.rhs(i), lhsBnd, deltaBnd()))
               {
                  Real lo    = -infinity;
                  Real scale = maxAbs(lp.rhs(i), lhsBnd);
         
                  if (scale < 1.0)
                     scale = 1.0;
                  
                  Real z = (lp.rhs(i) / scale) - (lhsBnd / scale);

                  if (isZero(z, epsZero()))
                     z = 0.0;
 
                  assert(lhsCnt > 0 || lp.upper(j) < infinity);

                  if (lhsCnt == 0)
                     lo = lp.upper(j) + z * scale / aij;
                  else if (lp.upper(j) >= infinity)
                     lo = z * scale / aij;
                 
                  if (isZero(lo, epsZero()))
                     lo = 0.0;

                  if (GErel(lo, lp.lower(j)))
                  {
                     MSG_INFO3( spxout << "IMAISM16 row " << i
                                       << ": redundant lower bound on x" << j
                                       << " -> lower=" << lo
                                       << " (" << lp.lower(j)
                                       << ")" << std::endl; )

                     ++rhsCnt;
                     rhsBnd -= aij * lp.lower(j);

                     lp.changeLower(j, -infinity);
                     ++chgBnds;                     
                  }
               }
            }
         }
      }
#endif

#if FREE_LHS_RHS
      // 2. detect implied free constraints
      if (lp.lhs(i) > -infinity && lhsCnt == 0 && GErel(lhsBnd, lp.lhs(i), deltaBnd()))
      {
         MSG_INFO3( spxout << "IMAISM17 row " << i 
                           << ": redundant lhs -> lhsBnd=" << lhsBnd
                           << " lhs=" << lp.lhs(i)
                           << std::endl; )
            
         lp.changeLhs(i, -infinity);
         ++chgLRhs;
      }
      if (lp.rhs(i) <  infinity && rhsCnt == 0 && LErel(rhsBnd, lp.rhs(i), deltaBnd()))
      {
         MSG_INFO3( spxout << "IMAISM18 row " << i 
                           << ": redundant rhs -> rhsBnd=" << rhsBnd
                           << " rhs=" << lp.rhs(i)
                           << std::endl; )
            
         lp.changeRhs(i, infinity);
         ++chgLRhs;
      }
#endif
      
      // 3. infeasible constraint
      if (LTrel(lp.rhs(i), lp.lhs(i), deltaBnd())                 ||
          (LTrel(rhsBnd,   lp.lhs(i), deltaBnd()) && rhsCnt == 0) ||
          (GTrel(lhsBnd,   lp.rhs(i), deltaBnd()) && lhsCnt == 0))
      {
         MSG_INFO3( spxout << "IMAISM19 row " << std::setprecision(20) << i 
                           << ": infeasible -> lhs=" << lp.lhs(i)
                           << " rhs=" << lp.rhs(i)
                           << " lhsBnd=" << lhsBnd
                           << " rhsBnd=" << rhsBnd
                           << std::endl; )
         return INFEASIBLE;
      }

#if FREE_CONSTRAINT      
      // 4. unconstrained constraint
      if (lp.lhs(i) <= -infinity && lp.rhs(i) >= infinity)
      {
         MSG_INFO3( spxout << "IMAISM20 row " << i 
                           << ": unconstrained -> removed" << std::endl; )

         m_hist.append(new FreeConstraintPS(lp, *this, i));
         
 	 ++remRows;
         remNzos += row.size();
         removeRow(lp, i);
         
         ++m_stat[FREE_ROW];
         
         continue;
      }
#endif

#if EMPTY_CONSTRAINT
      // 5. empty constraint
      if (row.size() == 0)
      {
         MSG_INFO3( spxout << "IMAISM21 row " << i
                           << ": empty ->"; )
            
         if (LT(lp.rhs(i), 0.0, deltaBnd()) || GT(lp.lhs(i), 0.0, deltaBnd()))
         {
            MSG_INFO3( spxout << " infeasible lhs=" << lp.lhs(i)
                              << " rhs=" << lp.rhs(i) << std::endl; )
            return INFEASIBLE;
         }         
         MSG_INFO3( spxout << " removed" << std::endl; )
         
         m_hist.append(new EmptyConstraintPS(rIdx(i)));   
        
         ++remRows;
         removeRow(lp, i);
         
         ++m_stat[EMPTY_ROW];

         continue;
      }
#endif

#if ROW_SINGLETON
      // 6. row singleton
      if (row.size() == 1)
      {
         Real aij = row.value(0);
         int  j   = row.index(0);
         Real lo  = -infinity;
         Real up  =  infinity;
         
         MSG_INFO3( spxout << "IMAISM22 row " << i 
                           << ": singleton -> val=" << aij 
                           << " lhs=" << lp.lhs(i)
                           << " rhs=" << lp.rhs(i); )
    
         if (GT(aij, 0.0, epsZero()))           // aij > 0
         {
            lo = (lp.lhs(i) <= -infinity) ? -infinity : lp.lhs(i) / aij;
            up = (lp.rhs(i) >=  infinity) ?  infinity : lp.rhs(i) / aij;
         }
         else if (LT(aij, 0.0, epsZero()))      // aij < 0
         {
            lo = (lp.rhs(i) >=  infinity) ? -infinity : lp.rhs(i) / aij;
            up = (lp.lhs(i) <= -infinity) ?  infinity : lp.lhs(i) / aij;
         }
         else if (LT(lp.rhs(i), 0.0, deltaBnd()) || GT(lp.lhs(i), 0.0, deltaBnd()))  
         {
            // aij == 0, rhs < 0 or lhs > 0
            MSG_INFO3( spxout << " infeasible" << std::endl; )

            return INFEASIBLE;
         }

         if (isZero(lo, epsZero()))
            lo = 0.0;

         if (isZero(up, epsZero()))
            up = 0.0;

         MSG_INFO3( spxout << " removed, lower=" << lo
                           << " (" << lp.lower(j)
                           << ") upper=" << up
                           << " (" << lp.upper(j)
                           << ")" << std::endl; )
            
         bool stricterUp = false;
         bool stricterLo = false;

         Real oldLo = lp.lower(j);
         Real oldUp = lp.upper(j);
         
         if (LTrel(up, lp.upper(j), deltaBnd()))
         {
            lp.changeUpper(j, up);
            stricterUp = true;
         }
         if (GTrel(lo, lp.lower(j), deltaBnd()))
         {
            lp.changeLower(j, lo);
            stricterLo = true;
         }

         m_hist.append(new RowSingletonPS(lp, *this, i, j, stricterLo, stricterUp, lp.lower(j), lp.upper(j), oldLo, oldUp));
         
         ++remRows;
         ++remNzos;
         removeRow(lp, i);
         
         ++m_stat[SINGLETON_ROW];
         
         continue;
      }
#endif

#if FORCE_CONSTRAINT
      // 7. forcing constraint (postsolving)
      // fix variables to obtain the upper bound on constraint value
      if (rhsCnt == 0 && EQrel(rhsBnd, lp.lhs(i), deltaBnd()))
      {
         MSG_INFO3( spxout << "IMAISM24 row " << i
                           << ": forcing constraint fix on lhs ->"
                           << " lhs=" << lp.lhs(i)
                           << " rhsBnd=" << rhsBnd
                           << std::endl; )

         for(int k = 0; k < row.size(); ++k)
         {
            Real aij = row.value(k);
            int  j   = row.index(k);
            
            ASSERT_WARN( "WMAISM25", isNotZero(aij, epsZero()) );
            
            if (aij > 0.0)
               lp.changeLower(j, lp.upper(j));
            else
               lp.changeUpper(j, lp.lower(j));
         }
         
         m_hist.append(new ForceConstraintPS(lp, *this, i, true));

         ++remRows;
         remNzos += row.size();
         removeRow(lp, i);
       
         ++m_stat[FORCE_ROW];

         continue;
      }
      // fix variables to obtain the lower bound on constraint value
      if (lhsCnt == 0 && EQrel(lhsBnd, lp.rhs(i), deltaBnd()))
      {
         MSG_INFO3( spxout << "IMAISM26 row " << i
                           << ": forcing constraint fix on rhs ->"
                           << " rhs=" << lp.rhs(i)
                           << " lhsBnd=" << lhsBnd
                           << std::endl; )

         for(int k = 0; k < row.size(); ++k)
         {
            Real aij   = row.value(k);
            int  j     = row.index(k);
            
            ASSERT_WARN( "WMAISM27", isNotZero(aij, epsZero()) );
            
            if (aij > 0.0)
               lp.changeUpper(j, lp.lower(j));
            else
               lp.changeLower(j, lp.upper(j));
         }
         m_hist.append(new ForceConstraintPS(lp, *this, i, false));
            
         ++remRows;
         remNzos += row.size();
         removeRow(lp, i);

         ++m_stat[FORCE_ROW];

         continue;
      }
#endif
   }

   assert(remRows > 0 || remNzos == 0);

   if (remRows + chgLRhs + chgBnds > 0)
   {
      again      = true;
      m_remRows += remRows;
      m_remNzos += remNzos;
      m_chgLRhs += chgLRhs;
      m_chgBnds += chgBnds;
   
      MSG_INFO2( spxout << "IMAISM28 Main simplifier (rows) removed "
                        << remRows << " rows, "
                        << remNzos << " non-zeros, "
                        << chgBnds << " col bounds, " 
                        << chgLRhs << " row bounds"
                        << std::endl; )

   }
   return OKAY;
}

SPxSimplifier::Result SPxMainSM::simplifyCols(SPxLP& lp, bool& again)
{  
   METHOD( "SPxMainSM::simplifyCols" );
  
   // This method simplifies the columns of the LP.
   //
   // The following operations are done:
   // 1. detect empty columns and fix corresponding variables
   // 2. detect variables that are unconstrained from below or above
   //    and fix corresponding variables or remove involved constraints
   // 3. fix variables
   // 4. use column singleton variables with zero objective to adjust constraint bounds 
   // 5. free column singleton combined with doubleton equation are
   //    used to make the column singleton variable free
   // 6. substitute (implied) free column singletons

   int remRows = 0;
   int remCols = 0;
   int remNzos = 0;
   int chgBnds = 0;
   
   for(int j = lp.nCols()-1; j >= 0; --j)
   {
       const SVector& col = lp.colVector(j);

      // infeasible bounds
      if (GTrel(lp.lower(j), lp.upper(j), deltaBnd()))
      {
         MSG_INFO3( spxout << "IMAISM29 col " << j
                           << ": infeasible bounds on x" << j
                           << " -> lower=" << lp.lower(j)
                           << " upper=" << lp.upper(j)
                           << std::endl; )
         return INFEASIBLE;
      }

      // 1. empty column
      if (col.size() == 0)
      {
#if EMPTY_COLUMN
	 MSG_INFO3( spxout << "IMAISM30 col " << j 
                           << ": empty -> maxObj=" << lp.maxObj(j)
                           << " lower=" << lp.lower(j)
                           << " upper=" << lp.upper(j); )

         Real val;

         if (GT(lp.maxObj(j), 0.0, epsZero()))
         {
            if (lp.upper(j) >= infinity)
            {
               MSG_INFO3( spxout << " unbounded" << std::endl; )
               return UNBOUNDED;
            }
            val = lp.upper(j);
         }
         else if (LT(lp.maxObj(j), 0.0, epsZero()))
         {
            if (lp.lower(j) <= -infinity)
            {
               MSG_INFO3( spxout << " unbounded" << std::endl; )
               return UNBOUNDED;
            }
            val = lp.lower(j);
         }
         else 
         {
            assert(isZero(lp.maxObj(j), epsZero()));
            // any value within the bounds is ok
            if (lp.lower(j) > -infinity)
               val = lp.lower(j);
            else if (lp.upper(j) < infinity)
               val = lp.upper(j);
            else
               val = 0.0;
         }
         MSG_INFO3( spxout << " removed" << std::endl; )

         m_hist.append(new FixBoundsPS(lp, *this, j, val));  
         m_hist.append(new FixVariablePS(lp, *this, j, val));
           
         ++remCols;
         removeCol(lp, j);

         ++m_stat[EMPTY_COL];

         continue;
#endif
      }
      
      if (NErel(lp.lower(j), lp.upper(j), deltaBnd()))
      {
         // will be set to false if any constraint implies a bound on the variable
         bool loFree = true;
         bool upFree = true;
      
         // 1. fix and remove variables
         for(int k = 0; k < col.size(); ++k)
         {            
            if (!loFree && !upFree)
               break;
            
            int i = col.index(k);
            
            // warn since this unhandled case may slip through unnoticed otherwise
            ASSERT_WARN( "WMAISM31", isNotZero(col.value(k), epsZero()) );

            if (col.value(k) > 0.0)
            {
               if (lp.rhs(i) <  infinity)
                  upFree = false;
               
               if (lp.lhs(i) > -infinity)
                  loFree = false;
            }
            else if (col.value(k) < 0.0)
            {
               if (lp.rhs(i) <  infinity)
                  loFree = false;
               
               if (lp.lhs(i) > -infinity)
                  upFree = false;
            }
         }
        
         // 2. detect variables that are unconstrained from below or above
         // max  3 x
         // s.t. 5 x >= 8
         if (GT(lp.maxObj(j), 0.0, epsZero()) && upFree)
         {
#if FIX_VARIABLE
            MSG_INFO3( spxout << "IMAISM32 col " << j
                              << ": x" << j
                              << " unconstrained above ->"; )
               
            if (lp.upper(j) >= infinity)
            {
               MSG_INFO3( spxout << " unbounded" << std::endl; )

               return UNBOUNDED;
            }         
            MSG_INFO3( spxout << " fixed at upper=" << lp.upper(j) << std::endl; )

            m_hist.append(new FixBoundsPS(lp, *this, j, lp.upper(j)));
            lp.changeLower(j, lp.upper(j));
         }
         // max -3 x
         // s.t. 5 x <= 8
         else if (LT(lp.maxObj(j), 0.0, epsZero()) && loFree)
         {
            MSG_INFO3( spxout << "IMAISM33 col " << j
                              << ": x" << j
                              << " unconstrained below ->"; )

            if (lp.lower(j) <= -infinity)
            {
               MSG_INFO3( spxout << " unbounded" << std::endl; )

               return UNBOUNDED;
            }
            MSG_INFO3( spxout << " fixed at lower=" << lp.lower(j) << std::endl; )
            
            m_hist.append(new FixBoundsPS(lp, *this, j, lp.lower(j)));
            lp.changeUpper(j, lp.lower(j));
#endif
         }
         else if (isZero(lp.maxObj(j), epsZero()))
         {
#if FREE_ZERO_OBJ_VARIABLE
            if (loFree && lp.lower(j) <= -infinity)
            {               
               MSG_INFO3( spxout << "IMAISM34 col " << j
                                 << ": x" << j
                                 << " unconstrained below with zero objective (" << lp.maxObj(j)
                                 << ")" << std::endl; )
         
               // variable j can be removed together with all constraints with variable j     
               for(int k = 0; k < col.size(); ++k)
               {
                  int l   =  col.index(k);
		  remNzos += lp.rowVector(l).size();
               }
               m_hist.append(new FreeZeroObjVariablePS(lp, *this, j, true));
               
               ++remCols;
               remRows += col.size();

	       for(int k = col.size()-1; k >= 0; --k)
                  removeRow(lp, col.index(k));
               
               removeCol(lp, j);
             
               ++m_stat[FREE_ZOBJ_COL];

               continue;
            }
            if (upFree && lp.upper(j) >= infinity)
            {
               MSG_INFO3( spxout << "IMAISM35 col " << j
                                 << ": x" << j
                                 << " unconstrained above with zero objective (" << lp.maxObj(j)
                                 << ")" << std::endl; )
         
               // variable j can be removed together with all constraints with variable j
               for(int k = 0; k < col.size(); ++k)
               {
                  int l   =  col.index(k);                 
                  remNzos += lp.rowVector(l).size();
               }
               m_hist.append(new FreeZeroObjVariablePS(lp, *this, j, false));
               
               ++remCols;
               remRows += col.size();

               for(int k = col.size()-1; k >= 0; --k)
                  removeRow(lp, col.index(k));
               
               removeCol(lp, j);
               
               ++m_stat[FREE_ZOBJ_COL];
               
               continue;
            }  
#endif
         }
      }
      
#if FIX_VARIABLE
      // 3. fix variable
      if (EQrel(lp.lower(j), lp.upper(j), deltaBnd()))
      {
         MSG_INFO3( spxout << "IMAISM36 col " << j
                           << ": x" << j
                           << " fixed -> lower=" << lp.lower(j)
                           << " upper=" << lp.upper(j) << std::endl; )         

         fixColumn(lp, j);
         
         ++remCols;
         remNzos += col.size();
         removeCol(lp, j);
        
         ++m_stat[FIX_COL];

         continue;
      }
#endif

      // handle column singletons
      if (col.size() == 1)
      {
         Real aij = col.value(0);
         int  i   = col.index(0);

         // 4. column singleton with zero objective
         if (isZero(lp.maxObj(j), epsZero()))
         {
#if ZERO_OBJ_COL_SINGLETON
            MSG_INFO3( spxout << "IMAISM37 col " << j
                              << ": singleton in row " << i
                              << " with zero objective"; )

            Real lhs = -infinity;
            Real rhs = +infinity;
            
            if (GT(aij, 0.0, epsZero()))
            {
               if (lp.lhs(i) > -infinity && lp.upper(j) <  infinity)
                  lhs = lp.lhs(i) - aij * lp.upper(j); 
               if (lp.rhs(i) <  infinity && lp.lower(j) > -infinity)
                  rhs = lp.rhs(i) - aij * lp.lower(j);
            }
            else if (LT(aij, 0.0, epsZero()))
            {
               if (lp.lhs(i) > -infinity && lp.lower(j) > -infinity)
                  lhs = lp.lhs(i) - aij * lp.lower(j);            
               if (lp.rhs(i) <  infinity && lp.upper(j) <  infinity)
                  rhs = lp.rhs(i) - aij * lp.upper(j);
            }
            else
            {
               lhs = lp.lhs(i);
               rhs = lp.rhs(i);
            }
         
            if (isZero(lhs, epsZero()))
               lhs = 0.0;
            if (isZero(rhs, epsZero()))
               rhs = 0.0;
            
            MSG_INFO3( spxout << " removed -> lhs=" << lhs
                              << " (" << lp.lhs(i)
                              << ") rhs=" << rhs
                              << " (" << lp.rhs(i)
                              << ")" << std::endl; )

            m_hist.append(new ZeroObjColSingletonPS(lp, *this, j, i));
            
            lp.changeRange(i, lhs, rhs);
 
            ++remCols;
            ++remNzos;
            removeCol(lp, j);
            
            ++m_stat[ZOBJ_SINGLETON_COL];

            if (lp.lhs(i) <= -infinity && lp.rhs(i) >= infinity)
            {
               m_hist.append(new FreeConstraintPS(lp, *this, i));
               
               ++remRows;
               removeRow(lp, i);

               ++m_stat[FREE_ROW];
            }

            continue;
#endif
         }
         
         // 5. not free column singleton combined with doubleton equation
         else if (EQrel(lp.lhs(i), lp.rhs(i), deltaBnd())             &&
                  lp.rowVector(i).size() == 2                         &&
                  (lp.lower(j) > -infinity || lp.upper(j) < infinity))
         {
#if DOUBLETON_EQUATION
            MSG_INFO3( spxout << "IMAISM38 col " << j
                              << ": singleton in row " << i
                              << " with doubleton equation ->"; )

            Real lhs = lp.lhs(i);
           
            const SVector& row = lp.rowVector(i);
            
            Real aik;
            int  k;

            if (row.index(0) == j)
            {
               aik = row.value(1);
               k   = row.index(1);
            }
            else if (row.index(1) == j)
            {
               aik = row.value(0);
               k   = row.index(0);
            }
            else
               throw SPxInternalCodeException("XMAISM11 This should never happen.");

            ASSERT_WARN( "WMAISM39", isNotZero(aik, epsZero()) );
           
            Real lo, up;
            Real oldLower = lp.lower(k);
            Real oldUpper = lp.upper(k);
            
            Real scale1 = maxAbs(lhs, aij * lp.upper(j));
            Real scale2 = maxAbs(lhs, aij * lp.lower(j));

            if (scale1 < 1.0)
               scale1 = 1.0;
            if (scale2 < 1.0)
               scale2 = 1.0;
            
            Real z1 = (lhs / scale1) - (aij * lp.upper(j) / scale1);
            Real z2 = (lhs / scale2) - (aij * lp.lower(j) / scale2);
            
            if (isZero(z1, epsZero()))
               z1 = 0.0;
            if (isZero(z2, epsZero()))
               z2 = 0.0;
         
            if (GT(aij * aik, 0.0, epsZero()))
            {
               lo = (lp.upper(j) >=  infinity) ? -infinity : z1 * scale1 / aik;
               up = (lp.lower(j) <= -infinity) ?  infinity : z2 * scale2 / aik;
            }            
            else if (LT(aij * aik, 0.0, epsZero()))
            {
               lo = (lp.lower(j) <= -infinity) ? -infinity : z2 * scale2 / aik;
               up = (lp.upper(j) >=  infinity) ?  infinity : z1 * scale1 / aik;
            }
            else
               throw SPxInternalCodeException("XMAISM12 This should never happen.");
            
            if (GTrel(lo, lp.lower(k), epsZero()))
               lp.changeLower(k, lo);
            
            if (LTrel(up, lp.upper(k), epsZero()))
               lp.changeUpper(k, up);
            
            MSG_INFO3( spxout << " made free, bounds on x" << k
                              << ": lower=" << lp.lower(k)
                              << " (" << oldLower
                              << ") upper=" << lp.upper(k)
                              << " (" << oldUpper
                              << ")" << std::endl; )
           
            // infeasible bounds
            if (GTrel(lp.lower(k), lp.upper(k), deltaBnd()))
            {
               MSG_INFO3( spxout << "new bounds are infeasible"
                                 << std::endl; )
               return INFEASIBLE;
            }
            
            m_hist.append(new DoubletonEquationPS(lp, *this, j, k, i, oldLower, oldUpper));

            if (lp.lower(j) > -infinity && lp.upper(j) < infinity)
               chgBnds += 2;
            else
               ++chgBnds;

            lp.changeBounds(j, -infinity, infinity);

            ++m_stat[DOUBLETON_ROW];
#endif
         }
         
         // 6. (implied) free column singleton
         if (lp.lower(j) <= -infinity && lp.upper(j) >= infinity)
         {
#if FREE_COL_SINGLETON
            Real slackVal = lp.lhs(i);
           
            // constraint i is an inequality constraint -> transform into equation type
            if (NErel(lp.lhs(i), lp.rhs(i), deltaBnd()))
            {
               MSG_INFO3( spxout << "IMAISM40 col " << j
                                 << ": free singleton in inequality constraint" << std::endl; )
               
               // do nothing if constraint i is unconstrained
               if (lp.lhs(i) <= -infinity && lp.rhs(i) >= infinity)
                  continue;
               
               // introduce slack variable to obtain equality constraint
               Real sMaxObj = lp.maxObj(j) / aij; // after substituting variable j in objective
               Real sLo     = lp.lhs(i);
               Real sUp     = lp.rhs(i);
               
               if (GT(sMaxObj, 0.0, epsZero()))
               {
                  if (sUp >= infinity)
                  {
                     MSG_INFO3( spxout << " -> problem unbounded" << std::endl; )
                     return UNBOUNDED;
                  }
                  slackVal = sUp;
               }
               else if (LT(sMaxObj, 0.0, epsZero()))
               {
                  if (sLo <= -infinity)
                  {
                     MSG_INFO3( spxout << " -> problem unbounded" << std::endl; )
                     return UNBOUNDED;
                  }
                  slackVal = sLo;
               }
               else 
               {
                  assert(isZero(sMaxObj, epsZero()));
                  // any value within the bounds is ok
                  if (sLo > -infinity)
                     slackVal = sLo;
                  else if (sUp < infinity)
                     slackVal = sUp;
                  else
                     throw SPxInternalCodeException("XMAISM13 This should never happen.");
               }
            }
           
            m_hist.append(new FreeColSingletonPS(lp, *this, j, i, slackVal));
            
            MSG_INFO3( spxout << "IMAISM41 col " << j
                              << ": free singleton removed" << std::endl; )
            
            const SVector& row = lp.rowVector(i);
            
            for (int h = 0; h < row.size(); ++h) 
            {
               int k = row.index(h);
               
               if (k != j)
               {
                  Real objAdd = -1.0 * lp.obj(j) * row.value(h) / aij;
                  lp.changeObj(k, lp.obj(k) + objAdd);
               }
            }
            
            ++remRows;
            ++remCols;
            remNzos += row.size();
            removeRow(lp, i);
            removeCol(lp, j);

            ++m_stat[FREE_SINGLETON_COL];
#endif
         }
      }
   }
   
   if (remCols + remRows > 0)
   {
      again      = true;
      m_remRows += remRows;
      m_remCols += remCols;
      m_remNzos += remNzos;
      m_chgBnds += chgBnds;

      MSG_INFO2( spxout << "IMAISM42 Main simplifier (columns) removed "
                        << remRows << " rows, "
                        << remCols << " cols, "
                        << remNzos << " non-zeros, "
                        << chgBnds << " col bounds" 
                        << std::endl; )
   }
   return OKAY;
}

SPxSimplifier::Result SPxMainSM::simplifyDual(SPxLP& lp, bool& again)
{
   METHOD( "SPxMainSM::simplifyDual" );
   
   // This method simplifies LP using the following dual structures:
   //
   // 1. dominated columns
   // 2. weakly dominated columns
   // 
   // For constructing the dual variables, it is assumed that the objective sense is max

   int remRows = 0;
   int remCols = 0;
   int remNzos = 0;

   DataArray<bool> colSingleton(lp.nCols());
   DVector         dualVarLo(lp.nRows());
   DVector         dualVarUp(lp.nRows());
   DVector         dualConsLo(lp.nCols());
   DVector         dualConsUp(lp.nCols());

   // init
   for(int i = lp.nRows()-1; i >= 0; --i)
   {
      // check for unconstrained constraints
      if (lp.lhs(i) <= -infinity && lp.rhs(i) >= infinity) 
      {
         MSG_INFO3( spxout << "IMAISM43 row " << i 
                           << ": unconstrained" << std::endl; )

	 m_hist.append(new FreeConstraintPS(lp, *this, i));

         ++remRows;
         remNzos += lp.rowVector(i).size();
         removeRow(lp, i);
         
         ++m_stat[FREE_ROW];
         
         continue;
      }
      
      // corresponds to maximization sense
      dualVarLo[i] = (lp.lhs(i) <= -infinity) ? 0.0 : -infinity;
      dualVarUp[i] = (lp.rhs(i) >=  infinity) ? 0.0 :  infinity;  
   }
   
   // compute bounds on the dual variables using column singletons
   for(int j = 0; j < lp.nCols(); ++j) 
   {
      if (lp.colVector(j).size() == 1)
      {
         int  i   = lp.colVector(j).index(0);
         Real aij = lp.colVector(j).value(0);
         
         ASSERT_WARN( "WMAISM44", isNotZero(aij, epsZero()) );
         
         Real bound = lp.maxObj(j) / aij;
         
         if (aij > 0)
         {
            if (lp.lower(j) <= -infinity && bound < dualVarUp[i])
               dualVarUp[i] = bound;
            if (lp.upper(j) >=  infinity && bound > dualVarLo[i])
               dualVarLo[i] = bound;
         }
         else if (aij < 0)
         {
            if (lp.lower(j) <= -infinity && bound > dualVarLo[i])
               dualVarLo[i] = bound;
            if (lp.upper(j) >=  infinity && bound < dualVarUp[i])
               dualVarUp[i] = bound;
         }             
      }
   }
   
   // compute bounds on the dual constraints
   for(int j = 0; j < lp.nCols(); ++j)
   {
      dualConsLo[j] = dualConsUp[j] = 0.0;
      
      const SVector& col = lp.colVector(j);
      
      for(int k = 0; k < col.size(); ++k)
      {
         if (dualConsLo[j] <= -infinity && dualConsUp[j] >= infinity)
            break;

         Real aij = col.value(k);
         int  i   = col.index(k);
         
         ASSERT_WARN( "WMAISM45", isNotZero(aij, epsZero()) );
         
         if (aij > 0)
         {
            if (dualVarLo[i] <= -infinity)
               dualConsLo[j] = -infinity;
            else 
               dualConsLo[j] += aij * dualVarLo[i];
            
            if (dualVarUp[i] >= infinity)
               dualConsUp[j] = infinity;
            else
               dualConsUp[j] += aij * dualVarUp[i];
         }
         else if (aij < 0)
         {
            if (dualVarLo[i] <= -infinity)
               dualConsUp[j] = infinity;
            else 
               dualConsUp[j] += aij * dualVarLo[i];
            
            if (dualVarUp[i] >= infinity)
               dualConsLo[j] = -infinity;
            else
               dualConsLo[j] += aij * dualVarUp[i];
         }
      }
   }
   
   for(int j = lp.nCols()-1; j >= 0; --j)
   {
      if (lp.colVector(j).size() <= 1)
         continue;
      
      // dual infeasibility checks
      if (LTrel(dualConsUp[j], dualConsLo[j], deltaBnd()))
      {
         MSG_INFO3( spxout << "IMAISM46 col " << j
                           << ": dual infeasible -> dual lhs bound=" << dualConsLo[j] 
                           << " dual rhs bound=" << dualConsUp[j] << std::endl; )
         return DUAL_INFEASIBLE;
      }      

      Real obj = lp.maxObj(j); 

      // 1. dominated column
      // Is the problem really unbounded in the cases below ??? Or is only dual infeasiblity be shown
      if (GTrel(obj, dualConsUp[j], deltaBnd()))
      {
#if DOMINATED_COLUMN
         MSG_INFO3( spxout << "IMAISM47 col " << j
                           << ": dominated -> maxObj=" << obj
                           << " dual rhs bound=" << dualConsUp[j] << std::endl; )
           
         if (lp.upper(j) >= infinity)
         {
            MSG_INFO2( spxout << " unbounded" << std::endl; )
            return UNBOUNDED;
         }
        
         MSG_INFO3( spxout << " fixed at upper=" << lp.upper(j) << std::endl; )
         
         m_hist.append(new FixBoundsPS(lp, *this, j, lp.upper(j)));
         lp.changeLower(j, lp.upper(j));

         ++m_stat[DOMINATED_COL];
#endif
      }
      else if (LTrel(obj, dualConsLo[j], deltaBnd()))
      {
#if DOMINATED_COLUMN
         MSG_INFO3( spxout << "IMAISM48 col " << j
                           << ": dominated -> maxObj=" << obj
                           << " dual lhs bound=" << dualConsLo[j] << std::endl; )
           
         if (lp.lower(j) <= -infinity)
         {
            MSG_INFO2( spxout << " unbounded" << std::endl; )
            return UNBOUNDED;
         }
         
         MSG_INFO3( spxout << " fixed at lower=" << lp.lower(j) << std::endl; )
         
         m_hist.append(new FixBoundsPS(lp, *this, j, lp.lower(j)));
         lp.changeUpper(j, lp.lower(j));

         ++m_stat[DOMINATED_COL];
#endif
      }

      // 2. weakly dominated column (no postsolving)
      else if (lp.upper(j) < infinity && EQrel(obj, dualConsUp[j], deltaBnd()))
      {
#if WEAKLY_DOMINATED_COLUMN
         MSG_INFO3( spxout << "IMAISM49 col " << j
                           << ": weakly dominated -> maxObj=" << obj
                           << " dual rhs bound=" << dualConsUp[j] << std::endl; )

         m_hist.append(new FixBoundsPS(lp, *this, j, lp.upper(j)));
         lp.changeLower(j, lp.upper(j));

         ++m_stat[WEAKLY_DOMINATED_COL];
#endif
      }        
      else if (lp.lower(j) > -infinity && EQrel(obj, dualConsLo[j], deltaBnd()))
      {
#if WEAKLY_DOMINATED_COLUMN
         MSG_INFO3( spxout << "IMAISM50 col " << j
                           << ": weakly dominated -> maxObj=" << obj
                           << " dual lhs bound=" << dualConsLo[j] << std::endl; )
         
         m_hist.append(new FixBoundsPS(lp, *this, j, lp.lower(j)));
         lp.changeUpper(j, lp.lower(j));

         ++m_stat[WEAKLY_DOMINATED_COL];
#endif
      }

      // fix column
      if (EQrel(lp.lower(j), lp.upper(j), deltaBnd()))
      {
#if FIX_VARIABLE
         fixColumn(lp, j);
         
         ++remCols;
         remNzos += lp.colVector(j).size();
         removeCol(lp, j);

         ++m_stat[FIX_COL];
#endif
      }
   }
   
   assert(remRows > 0 || remCols > 0 || remNzos == 0);

   if (remCols + remRows > 0)
   {
      again      = true;
      m_remRows += remRows;
      m_remCols += remCols;
      m_remNzos += remNzos;

      MSG_INFO2( spxout << "IMAISM51 Main simplifier (dual) removed "
                        << remRows << " rows, "
                        << remCols << " cols, "
                        << remNzos << " non-zeros"
                        << std::endl; )
   }
   return OKAY;
}

SPxSimplifier::Result SPxMainSM::duplicateRows(SPxLP& lp, bool& again)
{  
   METHOD( "SPxMainSM::duplicateRows" );

   // This method simplifies the LP by removing duplicate rows
   // Duplicates are detected using the algorithm of Bixby and Wagner [1987]
   
   // Possible extension: use generalized definition of duplicate rows according to Andersen and Andersen
   // However: the resulting sparsification is often very small since the involved rows are usually very sparse
   
   int remRows = 0;
   int remNzos = 0;

   // remove empty rows and columns
   SPxSimplifier::Result ret = removeEmpty(lp);
   if (ret != OKAY)
      return ret;
 
   if (lp.nRows() < 2)
      return OKAY;
    
   DataArray<int>    pClass(lp.nRows());           // class of parallel rows 
   DataArray<int>    classSize(lp.nRows());        // size of each class
   DataArray<double> scale(lp.nRows());            // scaling factor for each row
   int*              idxMem = 0;
   try
   {   
      spx_alloc(idxMem, lp.nRows());
      IdxSet idxSet(lp.nRows(), idxMem);           // set of feasible indices for new pClass
 
      // init
      pClass[0]    = 0;
      scale[0]     = 0.0;
      classSize[0] = lp.nRows();
      
      for(int i = 1; i < lp.nRows(); ++i)
      {
         pClass[i] = 0;
         scale[i]  = 0.0;
         classSize[i] = 0;
         idxSet.addIdx(i);
      }
      
      // stores parallel classes with non-zero colum entry     
      Array<DSVector> classSet(lp.nRows());
      Real oldVal;
   
      // main loop
      for(int j = 0; j < lp.nCols(); ++j)
      {
         const SVector& col = lp.colVector(j);
      
         for(int k = 0; k < col.size(); ++k)
         {
            Real aij = col.value(k);
            int  i   = col.index(k);
	
            ASSERT_WARN( "WMAISM52", isNotZero(aij, epsZero()) );
         
            if (scale[i] == 0.0)
               scale[i] = aij;

            classSet[pClass[i]].add(i, aij / scale[i]);
            if (--classSize[pClass[i]] == 0)
               idxSet.addIdx(pClass[i]);
         }

         // update each parallel class with non-zero column entry
         for(int m = 0; m < col.size(); ++m)
         {
            int k = pClass[col.index(m)];

            if (classSet[k].size() > 0)
            {            
               // sort classSet[k] w.r.t. scaled column values
               ElementCompare compare;

               if (classSet[k].size() > 1)
                  sorter_qsort(classSet[k].mem()+1, classSet[k].size(), compare);
           
               // use new index first
               int classIdx = idxSet.index(0);
               idxSet.remove(0);

               for(int l = 0; l < classSet[k].size(); ++l)
               {
                  if (l != 0 && NErel(classSet[k].value(l), oldVal, epsZero()))
                  {
                     classIdx = idxSet.index(0);
                     idxSet.remove(0);
                  }

                  pClass[classSet[k].index(l)] = classIdx; 
                  ++classSize[classIdx];
                  
                  oldVal = classSet[k].value(l);
               }
            
               classSet[k].clear();
            }
         }
      }  
   }
   catch(std::bad_alloc& x)
   {
      spx_free(idxMem);
      throw x;
   }
   catch(SPxMemoryException& x)
   {
      spx_free(idxMem);
      throw x;
   }
   spx_free(idxMem);
   
   // arrange duplicate rows using bucket sort w.r.t. their pClass values
   Array<DSVector> dupRows(lp.nRows());
   DataArray<bool> remRow(lp.nRows());
   
   for(int k = 0; k < dupRows.size(); ++k)
   {
      remRow[k] = false;
      dupRows[pClass[k]].add(k, 0.0);
   }
   
   // store rhs and lhs changes for combined update
   bool doChangeRanges = false;
   DVector newLhsVec(lp.lhs());
   DVector newRhsVec(lp.rhs());
   
   for(int k = 0; k < dupRows.size(); ++k)
   {
      if (dupRows[k].size() > 1)
      {
         MSG_INFO3( spxout << "IMAISM53 " << dupRows[k].size()
                           << " duplicate rows found" << std::endl; )
 
         m_stat[DUPLICATE_ROW] += dupRows[k].size()-1;

         // index of one non-column singleton row in dupRows[k]
         int  rowIdx    = -1;
         int  maxLhsIdx = -1;
         int  minRhsIdx = -1;
         Real maxLhs    = -infinity;
         Real minRhs    = +infinity;
     
         // determine strictest bounds on constraint
         for(int l = 0; l < dupRows[k].size(); ++l)
         {
            int i = dupRows[k].index(l);
               
            ASSERT_WARN( "WMAISM54", isNotZero(scale[i], epsZero()) );

            if (rowIdx == -1)
            {
               rowIdx = i;
               maxLhs = lp.lhs(rowIdx);
               minRhs = lp.rhs(rowIdx);
            }               
            else 
            {
               Real scaledLhs, scaledRhs;
               Real factor = scale[rowIdx] / scale[i];
                  
               if (factor > 0)
               {
                  scaledLhs = (lp.lhs(i) <= -infinity) ? -infinity : lp.lhs(i) * factor;
                  scaledRhs = (lp.rhs(i) >=  infinity) ?  infinity : lp.rhs(i) * factor;
               }
               if (factor < 0)
               {
                  scaledLhs = (lp.rhs(i) >=  infinity) ? -infinity : lp.rhs(i) * factor;
                  scaledRhs = (lp.lhs(i) <= -infinity) ?  infinity : lp.lhs(i) * factor;
               }
               if (scaledLhs > maxLhs)
               {
                  maxLhs    = scaledLhs;
                  maxLhsIdx = i;
               }
               if (scaledRhs < minRhs)
               {
                  minRhs    = scaledRhs;
                  minRhsIdx = i;
               }

               remRow[i] = true;
            }
         }
                  
         if (rowIdx != -1)
         {
            m_hist.append(new DuplicateRowsPS(lp, *this, rowIdx, maxLhsIdx, minRhsIdx, dupRows[k], scale));
            
            if (maxLhs > lp.lhs(rowIdx) || minRhs < lp.rhs(rowIdx))
            {
               // modify lhs and rhs of constraint rowIdx
               Real newLhs = (maxLhs > lp.lhs(rowIdx)) ? maxLhs : lp.lhs(rowIdx);  
               Real newRhs = (minRhs < lp.rhs(rowIdx)) ? minRhs : lp.rhs(rowIdx);
                     
               doChangeRanges = true;
               
               newLhsVec[rowIdx] = newLhs;
               newRhsVec[rowIdx] = newRhs;
                              
               if (LTrel(newRhs, newLhs, deltaBnd()))
               {
                  MSG_INFO3( spxout << "IMAISM55 duplicate rows yield infeasible bounds:"
                                    << " lhs=" << newLhs
                                    << " rhs=" << newRhs << std::endl; )
                  return INFEASIBLE;
               }
            }
         }
      }
   }

   // change ranges for all modified constraints by one single call (more efficient)
   if (doChangeRanges)
   {
      lp.changeRange(newLhsVec, newRhsVec);
   }
   
   // remove all rows by one single method call (more efficient)
   const int nRowsOld = lp.nRows();
   int* perm = new int[nRowsOld];
   
   for(int i = 0; i < nRowsOld; ++i)
   {
      if (remRow[i])
      {
         perm[i] = -1;
         ++remRows;
         remNzos += lp.rowVector(i).size();
      }
      else
         perm[i] = 0;
   }
   lp.removeRows(perm);
   
   for(int i = 0; i < nRowsOld; ++i)
   {
      if (perm[i] >= 0)
         m_rIdx[perm[i]] = m_rIdx[i];
   }   
   delete[] perm;   
   
   if (remRows + remNzos > 0)
   {
      again      = true;
      m_remRows += remRows;
      m_remNzos += remNzos;

      MSG_INFO2( spxout << "IMAISM56 Main simplifier (duplicate rows) removed "
                        << remRows << " rows, "
                        << remNzos << " non-zeros"
                        << std::endl; )
   }
   return OKAY;
}

SPxSimplifier::Result SPxMainSM::duplicateCols(SPxLP& lp, bool& again)
{  
   METHOD( "SPxMainSM::duplicateCols" );
   
   // This method simplifies the LP by removing duplicate columns
   // Duplicates are detected using the algorithm of Bixby and Wagner [1987]

   int remCols = 0;
   int remNzos = 0;
   
   // remove empty rows and columns
   SPxSimplifier::Result ret = removeEmpty(lp);
   if (ret != OKAY)
      return ret; 
   
   if (lp.nCols() < 2)
      return OKAY;
    
   DataArray<int>    pClass(lp.nCols());          // class of parallel columns
   DataArray<int>    classSize(lp.nCols());       // size of each class
   DataArray<double> scale(lp.nCols());           // scaling factor for each column
   int*              idxMem = 0;
   try
   {   
      spx_alloc(idxMem, lp.nCols());
      IdxSet idxSet(lp.nCols(), idxMem);  // set of feasible indices for new pClass
   
      // init
      pClass[0]    = 0;
      scale[0]     = 0.0;
      classSize[0] = lp.nCols();
      
      for(int j = 1; j < lp.nCols(); ++j)
      {
         pClass[j] = 0;
         scale[j]  = 0.0;
         classSize[j] = 0;
         idxSet.addIdx(j);
      }
      
      // stores parallel classes with non-zero row entry
      Array<DSVector> classSet(lp.nCols()); 

      Real oldVal;
   
      // main loop
      for(int i = 0; i < lp.nRows(); ++i)
      {
         const SVector& row = lp.rowVector(i);
      
         for(int k = 0; k < row.size(); ++k)
         {
            Real aij = row.value(k);
            int  j   = row.index(k);
            
            ASSERT_WARN( "WMAISM57", isNotZero(aij, epsZero()) );
            
            if (scale[j] == 0.0)
               scale[j] = aij;
            
            classSet[pClass[j]].add(j, aij / scale[j]);
            if (--classSize[pClass[j]] == 0)
               idxSet.addIdx(pClass[j]);
         }
            
         // update each parallel class with non-zero row entry
         for(int m = 0; m < row.size(); ++m)
         {
            int k = pClass[row.index(m)];
         
            if (classSet[k].size() > 0)
            {
               // sort classSet[k] w.r.t. scaled row values
               ElementCompare compare;
            
               if (classSet[k].size() > 1)
                  sorter_qsort(classSet[k].mem()+1, classSet[k].size(), compare); 
 
               // use new index first
               int classIdx = idxSet.index(0);
               idxSet.remove(0);
               
               for(int l = 0; l < classSet[k].size(); ++l)
               {
                  if (l != 0 && NErel(classSet[k].value(l), oldVal, epsZero()))
                  {
                     // start new parallel class
                     classIdx = idxSet.index(0);
                     idxSet.remove(0);
                  }
                  
                  pClass[classSet[k].index(l)] = classIdx; 
                  ++classSize[classIdx];

                  oldVal = classSet[k].value(l);
               }
            
               classSet[k].clear();
            }
         }
      } 
   }
   catch(std::bad_alloc& x)
   {
      spx_free(idxMem);
      throw x;
   }
   catch(SPxMemoryException& x)
   {
      spx_free(idxMem);
      throw x;
   }
    
   spx_free(idxMem);
  
   // arrange duplicate columns w.r.t. their pClass values
   Array<DSVector> dupCols(lp.nCols());
   
   for(int k = 0; k < dupCols.size(); ++k)
      dupCols[pClass[k]].add(k, 0.0);
   
   DataArray<bool> remCol(lp.nCols());
   
   for(int j = 0; j < lp.nCols(); ++j)
      remCol[j] = false;
   
   for(int k = 0; k < dupCols.size(); ++k)
   {
      if (dupCols[k].size() > 1)
      {
         MSG_INFO3( spxout << "IMAISM58 " << dupCols[k].size()
                           << " duplicate columns found" << std::endl; )
   
         for(int l = 0; l < dupCols[k].size(); ++l)
         {
            for(int m = 0; m < dupCols[k].size(); ++m)
            {               
               int j1  = dupCols[k].index(l);
               int j2  = dupCols[k].index(m);
                     
               if (l != m && !remCol[j1] && !remCol[j2])
               {
                  Real cj1 = lp.maxObj(j1);
                  Real cj2 = lp.maxObj(j2);
                  
                  // A.j1 = factor * A.j2
                  Real factor = scale[j1] / scale[j2];
                  Real objDif = cj1 - cj2 * scale[j1] / scale[j2];
                                    
                  ASSERT_WARN( "WMAISM59", isNotZero(factor, epsZero()) );
		  
                  if (isZero(objDif, epsZero()))
                  {
                     // case 1: objectives also duplicate
                     // variable substitution xj2' := xj2 + factor * xj1 <=> xj2 = -factor * xj1 + xj2'
                     m_hist.append(new DuplicateColsPS(lp, *this, j1, j2, factor));
                     
                     // update bounds of remaining column j2 (new column j2')
                     if (factor > 0)
                     {
                        if (lp.lower(j2) <= -infinity || lp.lower(j1) <= -infinity)
                           lp.changeLower(j2, -infinity);
                        else 
                           lp.changeLower(j2, lp.lower(j2) + factor * lp.lower(j1));
                        
                        if (lp.upper(j2) >= infinity || lp.upper(j1) >= infinity)
                           lp.changeUpper(j2, infinity);
                        else
                           lp.changeUpper(j2, lp.upper(j2) + factor * lp.upper(j1));
                     }
                     else if (factor < 0)
                     {
                        if (lp.lower(j2) <= -infinity || lp.upper(j1) >= infinity)
                           lp.changeLower(j2, -infinity);
                        else 
                           lp.changeLower(j2, lp.lower(j2) + factor * lp.upper(j1));
                        
                        if (lp.upper(j2) >= infinity || lp.lower(j1) <= -infinity)
                           lp.changeUpper(j2, infinity);
                        else
                           lp.changeUpper(j2, lp.upper(j2) + factor * lp.lower(j1));
                     }
                     
                     MSG_INFO3( spxout << "IMAISM60 two duplicate columns " << j1
                                       << ", " << j2
                                       << " replaced by one" << std::endl; )

                     remCol[j1] = true;

                     ++m_stat[SUB_DUPLICATE_COL];
                  }
                  else
                  {                     
                     // case 2: objectives not duplicate
                     // considered for maximization sense
                     if (lp.lower(j2) <= -infinity)
                     {
                        if (factor > 0 && objDif > 0)
                        {
                           if (lp.upper(j1) >= infinity)
                           {
                              MSG_INFO3( spxout << "IMAISM75 LP unbounded" << std::endl; )
                              return UNBOUNDED;
                           }
                           
                           // fix j1 at upper bound
                           MSG_INFO3( spxout << "IMAISM61 two duplicate columns " << j1
                                             << ", " << j2
                                             << " first one fixed at upper bound=" << lp.upper(j1) << std::endl; )
                           
                           m_hist.append(new FixBoundsPS(lp, *this, j1, lp.upper(j1)));
                           lp.changeLower(j1, lp.upper(j1));
                        }
                        else if (factor < 0 && objDif < 0)
                        {
                           if (lp.lower(j1) <= -infinity)
                           {
                              MSG_INFO3( spxout << "IMAISM76 LP unbounded" << std::endl; )
                              return UNBOUNDED;
                           }
                        
                           // fix j1 at lower bound
                           MSG_INFO3( spxout << "IMAISM62 two duplicate columns " << j1
                                             << ", " << j2
                                             << " first one fixed at lower bound=" << lp.lower(j1) << std::endl; )
                           
                           m_hist.append(new FixBoundsPS(lp, *this, j1, lp.lower(j1)));
                           lp.changeUpper(j1, lp.lower(j1));
                        }
                     }
                     else if (lp.upper(j2) >= infinity)
                     {
                        // fix j1 at upper bound
                        if (factor < 0 && objDif > 0)
                        {
                           if (lp.upper(j1) >= infinity)
                           {
                              MSG_INFO3( spxout << "IMAISM77 LP unbounded" << std::endl; )
                              return UNBOUNDED;
                           }
                        
                           // fix j1 at upper bound
                           MSG_INFO3( spxout << "IMAISM63 two duplicate columns " << j1
                                             << ", " << j2
                                             << " first one fixed at upper bound=" << lp.upper(j1) << std::endl; )
                           
                           m_hist.append(new FixBoundsPS(lp, *this, j1, lp.upper(j1)));
                           lp.changeLower(j1, lp.upper(j1));
                        }

                        // fix j1 at lower bound
                        else if (factor > 0 && objDif < 0)
                        {
                           if (lp.lower(j1) <= -infinity)
                           {
                              MSG_INFO3( spxout << "IMAISM78 LP unbounded" << std::endl; )
                              return UNBOUNDED;
                           }
                        
                           // fix j1 at lower bound
                           MSG_INFO3( spxout << "IMAISM64 two duplicate columns " << j1
                                             << ", " << j2
                                             << " first one fixed at lower bound=" << lp.lower(j1) << std::endl; )
                           
                           m_hist.append(new FixBoundsPS(lp, *this, j1, lp.lower(j1)));
                           lp.changeUpper(j1, lp.lower(j1));
                        } 
                     }
                     if (EQrel(lp.lower(j1), lp.upper(j1), deltaBnd()))
                     {
                        fixColumn(lp, j1);
                        remCol[j1] = true;

                        ++m_stat[FIX_DUPLICATE_COL];
                     }
                  }
               }
	    }
         }
      }     
   }
   
   // remove all columns by one single method call (more efficient)
   const int nColsOld = lp.nCols();
   int* perm = new int[nColsOld];
   
   for(int j = 0; j < nColsOld; ++j)
   {
      if (remCol[j])
      {
         perm[j] = -1;
         ++remCols;
         remNzos += lp.colVector(j).size();
      }
      else
         perm[j] = 0;
   }
   lp.removeCols(perm);
   
   for(int j = 0; j < nColsOld; ++j)
   {
      if (perm[j] >= 0)
         m_cIdx[perm[j]] = m_cIdx[j];
   }   
   delete[] perm;

   assert(remCols > 0 || remNzos == 0);

   if (remCols > 0)
   {
      again      = true;
      m_remCols += remCols;
      m_remNzos += remNzos;

      MSG_INFO2( spxout << "IMAISM65 Main simplifier (duplicate columns) removed "
                        << remCols << " cols, "
                        << remNzos << " non-zeros"
                        << std::endl; )   
   }
   return OKAY;
}

void SPxMainSM::fixColumn(SPxLP& lp, int j)
{
   METHOD( "SPxMainSM::fixColumn" );

   assert(EQrel(lp.lower(j), lp.upper(j), deltaBnd()));

   Real lo            = lp.lower(j);
   const SVector& col = lp.colVector(j);

   assert(NE(lo, infinity) && NE(lo, -infinity));

   MSG_INFO3( spxout << "IMAISM66 fix variable x" << j
                     << ": lower=" << lp.lower(j)
                     << " upper=" << lp.upper(j)
                     << std::endl; )
   
   if (isNotZero(lo, epsZero()))
   {
      for(int k = 0; k < col.size(); ++k)
      {
         int i = col.index(k);
         
         if (lp.rhs(i) < infinity)
         {
            Real y     = lo * col.value(k);
            Real scale = maxAbs(lp.rhs(i), y);
            
            if (scale < 1.0)
               scale = 1.0;
         
            Real rhs = (lp.rhs(i) / scale) - (y / scale);
                  
            if (isZero(rhs, epsZero()))
               rhs = 0.0;
            else
               rhs *= scale;

            MSG_INFO3( spxout << "IMAISM67 row " << i 
                              << ": rhs=" << rhs 
                              << " (" << lp.rhs(i) 
                              << ") aij=" << col.value(k) 
                              << std::endl; )

            lp.changeRhs(i, rhs);
         }
         if (lp.lhs(i) > -infinity)
         {
            Real y     = lo * col.value(k);
            Real scale = maxAbs(lp.lhs(i), y);
         
            if (scale < 1.0)
               scale = 1.0;
           
            Real lhs = (lp.lhs(i) / scale) - (y / scale);

            if (isZero(lhs, epsZero()))
               lhs = 0.0;
            else
               lhs *= scale;
                  
            MSG_INFO3( spxout << "IMAISM68 row " << i 
                              << ": lhs=" << lhs 
                              << " (" << lp.lhs(i) 
                              << ") aij=" << col.value(k) 
                              << std::endl; )

            lp.changeLhs(i, lhs);
         }
      }
   }
   
   m_hist.append(new FixVariablePS(lp, *this, j, lp.lower(j)));
}

SPxSimplifier::Result SPxMainSM::simplify(SPxLP& lp, Real eps, Real delta)
{
   METHOD( "SPxMainSM::simplify()" );
      
   m_timeUsed.reset();
   m_timeUsed.start();

   m_remRows = 0;
   m_remCols = 0;
   m_remNzos = 0;
   m_chgBnds = 0;
   m_chgLRhs = 0;
   
   Result ret   = OKAY;
   bool   again = true;

   m_prim.reDim(lp.nCols());
   m_slack.reDim(lp.nRows());
   m_dual.reDim(lp.nRows());
   m_redCost.reDim(lp.nCols());

   m_cBasisStat.reSize(lp.nCols());
   m_rBasisStat.reSize(lp.nRows());

   m_cIdx.reSize(lp.nCols());
   m_rIdx.reSize(lp.nRows());
   
   m_hist.reSize(0);
   m_postsolved = false;
   
   m_epsilon = eps;
   m_delta   = delta;
   
   for(int i = 0; i < lp.nRows(); ++i)
      m_rIdx[i] = i;
   
   for(int j = 0; j < lp.nCols(); ++j)
      m_cIdx[j] = j;

   m_stat.reSize(15);
   
   for(int k = 0; k < m_stat.size(); ++k)
      m_stat[k] = 0;
     
   // round extreme values (set all values smaller than eps to zero and all values bigger than infinity/5 to infinity)
#if EXTREMES
   handleExtremes(lp);
#endif

   // main presolving loop
   while(again && ret == OKAY) 
   {
      again = false;

#if ROWS      
      if (ret == OKAY)
         ret = simplifyRows(lp, again);
#endif

#if COLS
      if (ret == OKAY)
         ret = simplifyCols(lp, again);
#endif

#if DUAL
      if (ret == OKAY)
         ret = simplifyDual(lp, again);
#endif

#if DUPLICATE_ROWS
      if (ret == OKAY)
         ret = duplicateRows(lp, again);
#endif
      
#if DUPLICATE_COLS
      if (ret == OKAY)
         ret = duplicateCols(lp, again);
#endif
   }
   
   // preprocessing detected infeasibility or unboundness
   if (ret != OKAY)
      return ret;
   
   MSG_INFO1( spxout << "IMAISM69 Main simplifier removed "
                     << m_remRows << " rows, "
                     << m_remCols << " columns, "
                     << m_remNzos << " nonzeros, "
                     << m_chgBnds << " col bounds, " 
                     << m_chgLRhs << " row bounds"
                     << std::endl; )
      
   MSG_INFO1( spxout << "IMAISM74 Reduced LP has "
                     << lp.nRows() << " rows "
                     << lp.nCols() << " columns "
                     << lp.nNzos() << " nonzeros"
                     << std::endl; )   

   if (lp.nCols() == 0 && lp.nRows() == 0)
   {
      MSG_INFO1( spxout << "IMAISM70 Main simplifier removed all rows and columns" << std::endl; )
      ret = VANISHED;
   }
   
   MSG_INFO2( spxout << "\nIMAISM71 Main simplifier performed:\n"
                     << m_stat[EMPTY_ROW]            << " empty rows\n"
                     << m_stat[FREE_ROW]             << " free rows\n"
                     << m_stat[SINGLETON_ROW]        << " singleton rows\n"
                     << m_stat[FORCE_ROW]            << " forcing rows\n"
                     << m_stat[EMPTY_COL]            << " empty columns\n"
                     << m_stat[FIX_COL]              << " fixed columns\n"
                     << m_stat[FREE_ZOBJ_COL]        << " free columns with zero objective\n"
                     << m_stat[ZOBJ_SINGLETON_COL]   << " singleton columns with zero objective\n"
                     << m_stat[DOUBLETON_ROW]        << " singleton columns combined with a doubleton equation\n"
                     << m_stat[FREE_SINGLETON_COL]   << " free singleton columns\n"
                     << m_stat[DOMINATED_COL]        << " dominated columns\n"
                     << m_stat[WEAKLY_DOMINATED_COL] << " weakly dominated columns\n"
                     << m_stat[DUPLICATE_ROW]        << " duplicate rows\n"
                     << m_stat[FIX_DUPLICATE_COL]    << " duplicate columns (fixed)\n"
                     << m_stat[SUB_DUPLICATE_COL]    << " duplicate columns (substituted)\n"
                     << std::endl; );

   m_timeUsed.stop();
      
   return ret;
}

void SPxMainSM::unsimplify(const Vector& x, const Vector& y, const Vector& s, const Vector& r,
                           const SPxSolver::VarStatus rows[], const SPxSolver::VarStatus cols[])
{
   assert(x.dim() <= m_prim.dim());
   assert(y.dim() <= m_dual.dim());
   assert(x.dim() == r.dim());
   assert(y.dim() == s.dim());
      
   METHOD( "SPxMainSM::unsimplify()" );
   
   // assign values of variables in reduced LP
   for(int j = 0; j < x.dim(); ++j)
   {
      m_prim[cIdx(j)]       = isZero(x[j], epsZero()) ? 0.0 : x[j];
      m_redCost[cIdx(j)]    = isZero(r[j], epsZero()) ? 0.0 : r[j];
      m_cBasisStat[cIdx(j)] = cols[j];
   }
   for(int i = 0; i < y.dim(); ++i)
   {
      m_dual[rIdx(i)]       = isZero(y[i], epsZero()) ? 0.0 : y[i];
      m_slack[rIdx(i)]      = isZero(s[i], epsZero()) ? 0.0 : s[i];
      m_rBasisStat[rIdx(i)] = rows[i];
   }
   
   // undo preprocessing
   for(int k = m_hist.size()-1; k >= 0; --k)
   {
      const PostStep* psPtr = m_hist[k];
      
      psPtr->execute(m_prim, m_dual, m_slack, m_redCost, m_cBasisStat, m_rBasisStat);
      
      delete psPtr;
      psPtr = 0;
   }

   m_hist.clear();
   m_postsolved = true;
}

} //namespace soplex

//-----------------------------------------------------------------------------
//Emacs Local Variables:
//Emacs mode:c++
//Emacs c-basic-offset:3
//Emacs tab-width:8
//Emacs indent-tabs-mode:nil
//Emacs End:
//-----------------------------------------------------------------------------

