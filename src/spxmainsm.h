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
#pragma ident "@(#) $Id: spxmainsm.h,v 1.16 2009/08/11 12:48:40 bzfgleix Exp $"

/**@file  spxmainsm.h
 * @brief General methods in LP preprocessing.
 */
#ifndef _SPXMAINSM_H_
#define _SPXMAINSM_H_

#include <assert.h>

#include "spxdefines.h"
#include "spxsimplifier.h"
#include "array.h"
#include "exceptions.h"

namespace soplex
{
//---------------------------------------------------------------------
//  class SPxMainSM
//---------------------------------------------------------------------

/**@brief   LP simplifier for removing uneccessary row/columns.
   @ingroup Algo

   This #SPxSimplifier is mainly based on the paper "Presolving in
   linear programming" by E. Andersen and K. Andersen (Mathematical
   Programming, 1995).  It implements all proposed methods and some
   other preprocessing techniques for removing redundant rows and
   columns and bounds.  Also infeasibility and unboundness may be
   detected.

   Removed are:
   - empty rows / columns
   - unconstraint rows
   - row singletons
   - forcing rows
   - zero objective column singletons
   - (implied) free column singletons
   - doubleton equations combined with a column singleton
   - (implicitly) fixed columns
   - redundant lhs / rhs
   - redundant variable bounds
   - variables that are free in one direction
   - (weakly) dominated columns
   - duplicate rows / columns
*/
class SPxMainSM : public SPxSimplifier
{
private:
   //---------------------------------------------------------------------
   //  class PostsolveStep
   //---------------------------------------------------------------------
   
   /**@brief   Base class for postsolving operations.
      @ingroup Algo
      
      Class #PostStep is an abstract base class providing the
      interface for operations in the postsolving process.
   */
   class PostStep
   {
   private:
      const Real m_eps;
      
   public:
      /// constructor.
      PostStep()
         : m_eps(1e-6)
      {}
      /// copy constructor.
      PostStep(const PostStep& old)
         : m_eps(old.m_eps)
      {}
      /// assignment operator
      PostStep& operator=(const PostStep& /*rhs*/)
      {
         return *this;
      }
      /// destructor.
      virtual ~PostStep()
      {}
      /// clone function for polymorphism
      virtual PostStep* clone() const = 0;
      /// executes the postsolving.
      virtual void execute(DVector& x, DVector& y, DVector& s, DVector& r,
                           DataArray<SPxSolver::VarStatus>& cBasis, DataArray<SPxSolver::VarStatus>& rBasis) const = 0;

   protected:
      ///
      Real eps() const
      {
         return m_eps;
      }
   };

   /**@brief   Postsolves unconstraint constraints.
      @ingroup Algo
   */
   class FreeConstraintPS : public PostStep
   {
   private:
      const int m_i;
      DSVector  m_row;
   
   public:
      ///
      FreeConstraintPS(const SPxLP& lp, const SPxMainSM& simplifier, int i)
         : m_i(simplifier.rIdx(i))
         , m_row(lp.rowVector(i).size())
      {
         const SVector& row = lp.rowVector(i);
         
         for(int k = 0; k < row.size(); ++k)
            m_row.add(simplifier.cIdx(row.index(k)), row.value(k));
      }
      /// copy constructor
      FreeConstraintPS(const FreeConstraintPS& old)
         : PostStep(old)
         , m_i(old.m_i)
         , m_row(old.m_row)
      {}
      /// assignment operator
      FreeConstraintPS& operator=( const FreeConstraintPS& rhs)
      {
         if(this != &rhs)
         {
            PostStep::operator=(rhs);
            m_row = rhs.m_row;
         }

         return *this;
      }
      ///
      virtual void execute(DVector& x, DVector& y, DVector& s, DVector& r,
                           DataArray<SPxSolver::VarStatus>& cBasis, DataArray<SPxSolver::VarStatus>& rBasis) const;
      /// clone function for polymorphism
      inline virtual PostStep* clone() const
      {
         return new FreeConstraintPS(*this);
      }
   };

   /**@brief   Postsolves empty constraints.
      @ingroup Algo
   */
   class EmptyConstraintPS : public PostStep
   {
   private:
      const int m_i;
      
   public:
      ///
      EmptyConstraintPS(int i)
         : m_i(i)
      {}
      /// copy constructor
      EmptyConstraintPS(const EmptyConstraintPS& old)
         : PostStep(old)
         , m_i(old.m_i)
      {}
      /// assignment operator
      EmptyConstraintPS& operator=( const EmptyConstraintPS& rhs)
      {
         if(this != &rhs)
         {
            PostStep::operator=(rhs);
         }

         return *this;
      }
      ///
      virtual void execute(DVector& x, DVector& y, DVector& s, DVector& r,
                           DataArray<SPxSolver::VarStatus>& cBasis, DataArray<SPxSolver::VarStatus>& rBasis) const;
      /// clone function for polymorphism
      inline virtual PostStep* clone() const
      {
         return new EmptyConstraintPS(*this);
      }
   };
   
   /**@brief   Postsolves row singletons.
      @ingroup Algo
   */
   class RowSingletonPS : public PostStep
   {
   private:
      const int  m_i;
      const int  m_j;
      const bool m_strictLo;
      const bool m_strictUp;
      const bool m_maxSense;
      const Real m_obj;
      DSVector   m_col;
      const Real m_newLo;
      const Real m_newUp;
      const Real m_oldLo;
      const Real m_oldUp;
      
   public:
      ///
      RowSingletonPS(const SPxLP& lp, const SPxMainSM& simplifier, int i, int j, bool strictLo, bool strictUp,
                     Real newLo, Real newUp, Real oldLo, Real oldUp)
         : m_i(simplifier.rIdx(i))
         , m_j(simplifier.cIdx(j))
         , m_strictLo(strictLo)
         , m_strictUp(strictUp)
         , m_maxSense(lp.spxSense() == SPxLP::MAXIMIZE)
         , m_obj(lp.obj(j))
         , m_col(lp.colVector(j).size())
         , m_newLo(newLo)
         , m_newUp(newUp)
         , m_oldLo(oldLo)
         , m_oldUp(oldUp)
      {
         const SVector& col = lp.colVector(j);
         
         for(int k = 0; k < col.size(); ++k)
            m_col.add(simplifier.rIdx(col.index(k)), col.value(k));
      }
      /// copy constructor
      RowSingletonPS(const RowSingletonPS& old)
         : PostStep(old)
         , m_i(old.m_i)
         , m_j(old.m_j)
         , m_strictLo(old.m_strictLo)
         , m_strictUp(old.m_strictUp)
         , m_maxSense(old.m_maxSense)
         , m_obj(old.m_obj)
         , m_col(old.m_col)
         , m_newLo(old.m_newLo)
         , m_newUp(old.m_newUp)
         , m_oldLo(old.m_oldLo)
         , m_oldUp(old.m_oldUp)
      {}
      /// assignment operator
      RowSingletonPS& operator=( const RowSingletonPS& rhs)
      {
         if(this != &rhs)
         {
            PostStep::operator=(rhs);
            m_col = rhs.m_col;
         }

         return *this;
      }
      /// clone function for polymorphism
      inline virtual PostStep* clone() const
      {
         return new RowSingletonPS(*this);
      }
      ///
      virtual void execute(DVector& x, DVector& y, DVector& s, DVector& r,
                           DataArray<SPxSolver::VarStatus>& cBasis, DataArray<SPxSolver::VarStatus>& rBasis) const;
   };
 
   /**@brief   Postsolves forcing constraints.
      @ingroup Algo
   */
   class ForceConstraintPS : public PostStep
   {
   private:
      const int       m_i;
      const Real      m_lRhs;
      DSVector        m_row;
      DSVector        m_objs;
      DataArray<bool> m_fixed;
      Array<DSVector> m_cols;
      const bool      m_lhsFixed;
      const bool      m_maxSense;
      
   public:
      ///
      ForceConstraintPS(const SPxLP& lp, const SPxMainSM& simplifier, int i, bool lhsFixed)
         : m_i(simplifier.rIdx(i))
         , m_lRhs(lhsFixed ? lp.lhs(i) : lp.rhs(i))
         , m_row(lp.rowVector(i).size())
         , m_objs(lp.rowVector(i).size())
         , m_fixed(lp.rowVector(i).size())
         , m_cols(lp.rowVector(i).size())
         , m_lhsFixed(lhsFixed)
         , m_maxSense(lp.spxSense() == SPxLP::MAXIMIZE)
      {
         const SVector& row = lp.rowVector(i);
         
         for(int k = 0; k < row.size(); ++k)
         {
            int j = simplifier.cIdx(row.index(k));
            m_row.add(j, row.value(k));
            m_objs.add(j, lp.obj(row.index(k)));
            m_fixed[k] = EQrel(lp.lower(row.index(k)), lp.upper(row.index(k)));
            
            const SVector& col = lp.colVector(row.index(k));
            m_cols[k].setMax(col.size());
            
            for(int l = 0; l < col.size(); ++l)
               m_cols[k].add(simplifier.rIdx(col.index(l)), col.value(l));
         }
      }
      /// copy constructor
      ForceConstraintPS(const ForceConstraintPS& old)
         : PostStep(old)
         , m_i(old.m_i)
         , m_lRhs(old.m_lRhs)
         , m_row(old.m_row)
         , m_objs(old.m_objs)
         , m_fixed(old.m_fixed)
         , m_cols(old.m_cols)
         , m_lhsFixed(old.m_lhsFixed)
         , m_maxSense(old.m_maxSense)
      {}
      /// assignment operator
      ForceConstraintPS& operator=( const ForceConstraintPS& rhs)
      {
         if(this != &rhs)
         {
            PostStep::operator=(rhs);
            m_row = rhs.m_row;
            m_objs = rhs.m_objs;
            m_fixed = rhs.m_fixed;
            m_cols = rhs.m_cols;
         }

         return *this;
      }
      /// clone function for polymorphism
      inline virtual PostStep* clone() const
      {
         return new ForceConstraintPS(*this);
      }
      ///
      virtual void execute(DVector& x, DVector& y, DVector& s, DVector& r,
                           DataArray<SPxSolver::VarStatus>& cBasis, DataArray<SPxSolver::VarStatus>& rBasis) const;
   };
   
   /**@brief   Postsolves variable fixing.
      @ingroup Algo
   */
   class FixVariablePS : public PostStep
   {
   private:
      const int  m_j;
      const Real m_val;
      const Real m_obj;
      DSVector   m_col;
      
   public:
      ///
      FixVariablePS(const SPxLP& lp, const SPxMainSM& simplifier, int j, const Real val)
         : m_j(simplifier.cIdx(j))
         , m_val(val)
         , m_obj(lp.obj(j))
         , m_col(lp.colVector(j).size())
      {
         const SVector& col = lp.colVector(j);
         
         for(int k = 0; k < col.size(); ++k)
            m_col.add(simplifier.rIdx(col.index(k)), col.value(k));
      }
      /// copy constructor
      FixVariablePS(const FixVariablePS& old)
         : PostStep(old)
         , m_j(old.m_j)
         , m_val(old.m_val)
         , m_obj(old.m_obj)
         , m_col(old.m_col)
      {}
      /// assignment operator
      FixVariablePS& operator=( const FixVariablePS& rhs)
      {
         if(this != &rhs)
         {
            PostStep::operator=(rhs);
            m_col = rhs.m_col;
         }

         return *this;
      }
      /// clone function for polymorphism
      inline virtual PostStep* clone() const
      {
         return new FixVariablePS(*this);
      }
      ///
      virtual void execute(DVector& x, DVector& y, DVector& s, DVector& r,
                           DataArray<SPxSolver::VarStatus>& cBasis, DataArray<SPxSolver::VarStatus>& rBasis) const;
   };

   /**@brief   Postsolves variable bound fixing.
      @ingroup Algo
   */
   class FixBoundsPS : public PostStep
   {
   private:
      const int            m_j;
      SPxSolver::VarStatus m_status;
      
   public:
      ///
      FixBoundsPS(const SPxLP& lp, const SPxMainSM& simplifier, int j, Real val)
         : m_j(simplifier.cIdx(j))
      { 
         if (EQrel(lp.lower(j), lp.upper(j), eps()))
            m_status = SPxSolver::FIXED;
         else if (EQrel(val, lp.lower(j), eps()))
            m_status = SPxSolver::ON_LOWER;
         else if (EQrel(val, lp.upper(j), eps()))
            m_status = SPxSolver::ON_UPPER;
         else if (lp.lower(j) <= -infinity && lp.upper(j) >= infinity)
            m_status = SPxSolver::ZERO;
         else
         {
            throw SPxInternalCodeException("XMAISM14 This should never happen.");
         }
      }
      /// copy constructor
      FixBoundsPS(const FixBoundsPS& old)
         : PostStep(old)
         , m_j(old.m_j)
         , m_status(old.m_status)
      {}
      /// assignment operator
      FixBoundsPS& operator=( const FixBoundsPS& rhs)
      {
         if(this != &rhs)
         {
            PostStep::operator=(rhs);
            m_status = rhs.m_status;
         }

         return *this;
      }
      /// clone function for polymorphism
      inline virtual PostStep* clone() const
      {
         return new FixBoundsPS(*this);
      }
      ///
      virtual void execute(DVector& x, DVector& y, DVector& s, DVector& r,
                           DataArray<SPxSolver::VarStatus>& cBasis, DataArray<SPxSolver::VarStatus>& rBasis) const;
   };

   /**@brief Postsolves the case when constraints are removed due to a
             variable with zero objective that is free in one direction.
      @ingroup Algo
   */
   class FreeZeroObjVariablePS : public PostStep
   {
   private:
      const int       m_j;
      const Real      m_bnd;
      DSVector        m_col;
      DSVector        m_lRhs;
      Array<DSVector> m_rows;
      const bool      m_loFree;
      
   public:
      ///
      FreeZeroObjVariablePS(const SPxLP& lp, const SPxMainSM& simplifier, int j, bool loFree)
         : m_j(simplifier.cIdx(j))
         , m_bnd(loFree ? lp.upper(j) : lp.lower(j))
         , m_col(lp.colVector(j).size())
         , m_lRhs(lp.colVector(j).size())
         , m_rows(lp.colVector(j).size())
         , m_loFree(loFree)
      {
         const SVector& col = lp.colVector(j);
         
         for(int k = 0; k < col.size(); ++k)
         {
            int i = simplifier.rIdx(col.index(k));

            m_col.add(i, col.value(k));

	    assert(isNotZero(col.value(k)));
            
            if ((m_loFree  && col.value(k) > 0) || 
                (!m_loFree && col.value(k) < 0))
               m_lRhs.add(i, lp.rhs(col.index(k)));
            else
               m_lRhs.add(i, lp.lhs(col.index(k)));
            
            const SVector& row = lp.rowVector(col.index(k));
            m_rows[k].setMax(row.size());
            
            for(int l = 0; l < row.size(); ++l)
               m_rows[k].add(simplifier.cIdx(row.index(l)), row.value(l));
         }
      }
      /// copy constructor
      FreeZeroObjVariablePS(const FreeZeroObjVariablePS& old)
         : PostStep(old)
         , m_j(old.m_j)
         , m_bnd(old.m_bnd)
         , m_col(old.m_col)
         , m_lRhs(old.m_lRhs)
         , m_rows(old.m_rows)
         , m_loFree(old.m_loFree)
      {}
      /// assignment operator
      FreeZeroObjVariablePS& operator=( const FreeZeroObjVariablePS& rhs)
      {
         if(this != &rhs)
         {
            PostStep::operator=(rhs);
            m_col = rhs.m_col;
            m_lRhs = rhs.m_lRhs;
            m_rows = rhs.m_rows;
         }

         return *this;
      }
      /// clone function for polymorphism
      inline virtual PostStep* clone() const
      {
         return new FreeZeroObjVariablePS(*this);
      }
      ///
      virtual void execute(DVector& x, DVector& y, DVector& s, DVector& r,
                           DataArray<SPxSolver::VarStatus>& cBasis, DataArray<SPxSolver::VarStatus>& rBasis) const;
   };

   /**@brief   Postsolves column singletons with zero objective.
      @ingroup Algo
   */
   class ZeroObjColSingletonPS : public PostStep
   {
   private:
      const int  m_j;
      const int  m_i;
      const Real m_lhs;
      const Real m_rhs;
      const Real m_lower;
      const Real m_upper;
      DSVector   m_row;
      
    public:
      ///
      ZeroObjColSingletonPS(const SPxLP& lp, const SPxMainSM& simplifier, int j, int i)
         : m_j(simplifier.cIdx(j))
         , m_i(simplifier.rIdx(i))
         , m_lhs(lp.lhs(i))
         , m_rhs(lp.rhs(i))
         , m_lower(lp.lower(j))
         , m_upper(lp.upper(j))
         , m_row(lp.rowVector(i).size())
      {
         const SVector& row = lp.rowVector(i);
      
         for(int k = 0; k < row.size(); ++k)
            m_row.add(simplifier.cIdx(row.index(k)), row.value(k));
      }
      /// copy constructor
      ZeroObjColSingletonPS(const ZeroObjColSingletonPS& old)
         : PostStep(old)
         , m_j(old.m_j)
         , m_i(old.m_i)
         , m_lhs(old.m_lhs)
         , m_rhs(old.m_rhs)
         , m_lower(old.m_lower)
         , m_upper(old.m_upper)
         , m_row(old.m_row)
      {}
      /// assignment operator
      ZeroObjColSingletonPS& operator=( const ZeroObjColSingletonPS& rhs)
      {
         if(this != &rhs)
         {
            PostStep::operator=(rhs);
            m_row = rhs.m_row;
         }

         return *this;
      }
      /// clone function for polymorphism
      inline virtual PostStep* clone() const
      {
         return new ZeroObjColSingletonPS(*this);
      }
      ///
      virtual void execute(DVector& x, DVector& y, DVector& s, DVector& r,
                           DataArray<SPxSolver::VarStatus>& cBasis, DataArray<SPxSolver::VarStatus>& rBasis) const;
   };
   
   /**@brief   Postsolves free column singletons.
      @ingroup Algo
   */
   class FreeColSingletonPS : public PostStep
   {
   private:
      const int  m_j;
      const int  m_i;
      const Real m_obj;
      const Real m_lRhs;
      const bool m_onLhs;
      const bool m_eqCons;
      DSVector   m_row;
      
   public:
      ///
      FreeColSingletonPS(const SPxLP& lp, const SPxMainSM& simplifier, int j, int i, Real slackVal)
         : m_j(simplifier.cIdx(j))
         , m_i(simplifier.rIdx(i))
         , m_obj(lp.obj(j))
         , m_lRhs(slackVal)
         , m_onLhs(slackVal == lp.lhs(i))
         , m_eqCons(EQrel(lp.lhs(i), lp.rhs(i)))
         , m_row(lp.rowVector(i).size())
      {
         const SVector& row = lp.rowVector(i);
      
         for(int k = 0; k < row.size(); ++k)
            m_row.add(simplifier.cIdx(row.index(k)), row.value(k));
      }
      /// copy constructor
      FreeColSingletonPS(const FreeColSingletonPS& old)
         : PostStep(old)
         , m_j(old.m_j)
         , m_i(old.m_i)
         , m_obj(old.m_obj)
         , m_lRhs(old.m_lRhs)
         , m_onLhs(old.m_onLhs)
         , m_eqCons(old.m_eqCons)
         , m_row(old.m_row)
      {}
      /// assignment operator
      FreeColSingletonPS& operator=( const FreeColSingletonPS& rhs)
      {
         if(this != &rhs)
         {
            PostStep::operator=(rhs);
            m_row = rhs.m_row;
         }

         return *this;
      }
      /// clone function for polymorphism
      inline virtual PostStep* clone() const
      {
         return new FreeColSingletonPS(*this);
      }
      ///
      virtual void execute(DVector& x, DVector& y, DVector& s, DVector& r,
                           DataArray<SPxSolver::VarStatus>& cBasis, DataArray<SPxSolver::VarStatus>& rBasis) const;
   };

   /**@brief   Postsolves doubleton equations combined with a column singleton.
      @ingroup Algo
   */
   class DoubletonEquationPS : public PostStep
   {
   private:
      const int  m_j;
      const int  m_k;
      const int  m_i;
      const bool m_maxSense;
      const bool m_jFixed;
      const Real m_jObj;
      const Real m_kObj;
      const Real m_aij;
      const bool m_strictLo;
      const bool m_strictUp;
      const Real m_newLo;
      const Real m_newUp;
      const Real m_oldLo;
      const Real m_oldUp;
      DSVector   m_col;
      
   public:
      ///
      DoubletonEquationPS(const SPxLP& lp, const SPxMainSM& simplifier, int j, int k, int i, Real oldLo, Real oldUp)
         : m_j(simplifier.cIdx(j))
         , m_k(simplifier.cIdx(k))
         , m_i(simplifier.rIdx(i))
         , m_maxSense(lp.spxSense() == SPxLP::MAXIMIZE)
         , m_jFixed(EQrel(lp.lower(j), lp.upper(j)))
         , m_jObj(lp.obj(j))
         , m_kObj(lp.obj(k))
         , m_aij(lp.colVector(j).value(0))
         , m_strictLo(lp.lower(k) > oldLo)
         , m_strictUp(lp.upper(k) < oldUp)
         , m_newLo(lp.lower(k))
         , m_newUp(lp.upper(k))
         , m_oldLo(oldLo)
         , m_oldUp(oldUp)
         , m_col(lp.colVector(k).size())
      {
         const SVector& col = lp.colVector(k);
         
         for(int l = 0; l < col.size(); ++l)
            m_col.add(simplifier.rIdx(col.index(l)), col.value(l));
      }
      /// copy constructor
      DoubletonEquationPS(const DoubletonEquationPS& old)
         : PostStep(old)
         , m_j(old.m_j)
         , m_k(old.m_k)
         , m_i(old.m_i)
         , m_maxSense(old.m_maxSense)
         , m_jFixed(old.m_jFixed)
         , m_jObj(old.m_jObj)
         , m_kObj(old.m_kObj)
         , m_aij(old.m_aij)
         , m_strictLo(old.m_strictLo)
         , m_strictUp(old.m_strictUp)
         , m_newLo(old.m_newLo)
         , m_newUp(old.m_newUp)
         , m_oldLo(old.m_oldLo)
         , m_oldUp(old.m_oldUp)
         , m_col(old.m_col)
      {}
      /// assignment operator
      DoubletonEquationPS& operator=( const DoubletonEquationPS& rhs)
      {
         if(this != &rhs)
         {
            PostStep::operator=(rhs);
            m_col = rhs.m_col;
         }

         return *this;
      }
      /// clone function for polymorphism
      inline virtual PostStep* clone() const
      {
         return new DoubletonEquationPS(*this);
      }
      ///
      virtual void execute(DVector& x, DVector& y, DVector& s, DVector& r,
                           DataArray<SPxSolver::VarStatus>& cBasis, DataArray<SPxSolver::VarStatus>& rBasis) const;
   };

   /**@brief   Postsolves duplicate rows.
      @ingroup Algo
   */
   class DuplicateRowsPS : public PostStep
   {
   private:
      const int       m_i;
      const int       m_maxLhsIdx;
      const int       m_minRhsIdx;
      const bool      m_maxSense;
      DSVector        m_scale;
      
   public:
      DuplicateRowsPS(const SPxLP& lp, const SPxMainSM& simplifier, int i,
                      int maxLhsIdx, int minRhsIdx, const DSVector& dupRows, const DataArray<double> scale)
         : m_i(simplifier.rIdx(i))
         , m_maxLhsIdx((maxLhsIdx == -1) ? -1 : simplifier.rIdx(maxLhsIdx))
         , m_minRhsIdx((minRhsIdx == -1) ? -1 : simplifier.rIdx(minRhsIdx))
         , m_maxSense(lp.spxSense() == SPxLP::MAXIMIZE)
         , m_scale(dupRows.size())
      {
         Real rowScale = scale[i];
         
         for(int k = 0; k < dupRows.size(); ++k)
            m_scale.add(simplifier.rIdx(dupRows.index(k)), rowScale / scale[dupRows.index(k)]);    
      }
      /// copy constructor
      DuplicateRowsPS(const DuplicateRowsPS& old)
         : PostStep(old)
         , m_i(old.m_i)
         , m_maxLhsIdx(old.m_maxLhsIdx)
         , m_minRhsIdx(old.m_minRhsIdx)
         , m_maxSense(old.m_maxSense)  
         , m_scale(old.m_scale)
      {}
      /// assignment operator
      DuplicateRowsPS& operator=( const DuplicateRowsPS& rhs)
      {
         if(this != &rhs)
         {
            PostStep::operator=(rhs);
            m_scale = rhs.m_scale;
         }

         return *this;
      }
      /// clone function for polymorphism
      inline virtual PostStep* clone() const
      {
         return new DuplicateRowsPS(*this);
      }
      virtual void execute(DVector& x, DVector& y, DVector& s, DVector& r,
                           DataArray<SPxSolver::VarStatus>& cBasis, DataArray<SPxSolver::VarStatus>& rBasis) const;
   };

   /**@brief   Postsolves duplicate columns.
      @ingroup Algo
   */
   class DuplicateColsPS : public PostStep
   {
   private:
      const int            m_j;
      const int            m_k;
      const Real           m_loJ;
      const Real           m_upJ;
      const Real           m_loK;
      const Real           m_upK;
      const Real           m_scale;
      
   public:
      DuplicateColsPS(const SPxLP& lp, const SPxMainSM& simplifier, int j, int k, Real scale)
         : m_j(simplifier.cIdx(j))
         , m_k(simplifier.cIdx(k))
         , m_loJ(lp.lower(j))
         , m_upJ(lp.upper(j))
         , m_loK(lp.lower(k))
         , m_upK(lp.upper(k))
         , m_scale(scale)
      {}
      /// copy constructor
      DuplicateColsPS(const DuplicateColsPS& old)
         : PostStep(old)
         , m_j(old.m_j)
         , m_k(old.m_k)
         , m_loJ(old.m_loJ)
         , m_upJ(old.m_upJ)
         , m_loK(old.m_loK)
         , m_upK (old.m_upK)
         , m_scale (old.m_scale)
      {}
      /// assignment operator
      DuplicateColsPS& operator=( const DuplicateColsPS& rhs)
      {
         if(this != &rhs)
         {
            PostStep::operator=(rhs);
         }

         return *this;
      }
      /// clone function for polymorphism
      inline virtual PostStep* clone() const
      {
         return new DuplicateColsPS(*this);
      }
      virtual void execute(DVector& x, DVector& y, DVector& s, DVector& r,
                           DataArray<SPxSolver::VarStatus>& cBasis, DataArray<SPxSolver::VarStatus>& rBasis) const;
   };

   // friends 
   friend class FreeConstraintPS;
   friend class EmptyConstraintPS;
   friend class RowSingletonPS;
   friend class ForceConstraintPS;
   friend class FixVariablePS;
   friend class FixBoundsPS;
   friend class FreeZeroObjVariablePS;
   friend class ZeroObjColSingletonPS;
   friend class FreeColSingletonPS;
   friend class DoubletonEquationPS;
   friend class DuplicateRowsPS;
   friend class DuplicateColsPS;   
  
private:
   //------------------------------------
   //**@name Types */
   //@{
   /// Different simplification steps.
   enum SimpleStep
   {
      EMPTY_ROW            =  0,
      FREE_ROW             =  1,
      SINGLETON_ROW        =  2,
      FORCE_ROW            =  3,
      EMPTY_COL            =  4,
      FIX_COL              =  5,
      FREE_ZOBJ_COL        =  6,
      ZOBJ_SINGLETON_COL   =  7,
      DOUBLETON_ROW        =  8,
      FREE_SINGLETON_COL   =  9,
      DOMINATED_COL        = 10,
      WEAKLY_DOMINATED_COL = 11,
      DUPLICATE_ROW        = 12,
      FIX_DUPLICATE_COL    = 13,
      SUB_DUPLICATE_COL    = 14
   };
   //@}

   //------------------------------------
   //**@name Data */
   //@{
   ///
   DVector                         m_prim;       ///< unsimplified primal solution vector.
   DVector                         m_slack;      ///< unsimplified slack vector.
   DVector                         m_dual;       ///< unsimplified dual solution vector.
   DVector                         m_redCost;    ///< unsimplified reduced cost vector.   
   DataArray<SPxSolver::VarStatus> m_cBasisStat; ///< basis status of columns.
   DataArray<SPxSolver::VarStatus> m_rBasisStat; ///< basis status of rows.
   DataArray<int>                  m_cIdx;       ///< column index vector in original LP.
   DataArray<int>                  m_rIdx;       ///< row index vector in original LP.
   DataArray<PostStep*>            m_hist;       ///< vector of presolve history.
   bool                            m_postsolved; ///< status of postsolving.
   Real                            m_epsilon;    ///< epsilon zero.
   Real                            m_delta;      ///< maximum bound violation.
   DataArray<int>                  m_stat;       ///< preprocessing history.
   //@}
     
private:
   //------------------------------------
   //**@name Private helpers */
   //@{
   /// handles extreme values by setting them to zero or infinity.
   void handleExtremes(SPxLP& lp);
   
   /// removed empty rows and empty columns.
   Result removeEmpty(SPxLP& lp);
   
   /// performs simplification steps on the rows of the LP.
   Result simplifyRows(SPxLP& lp, bool& again);
   
   /// performs simplification steps on the columns of the LP.
   Result simplifyCols(SPxLP& lp, bool& again);
   
   /// performs simplification steps on the LP based on dual concepts.
   Result simplifyDual(SPxLP& lp, bool& again);
   
   /// removes duplicate rows.
   Result duplicateRows(SPxLP& lp, bool& again);
   
   /// removes duplicate columns
   Result duplicateCols(SPxLP& lp, bool& again);
   
   /// handles the fixing of a variable.
   void fixColumn(SPxLP& lp, int i);
   
   /// removes a row in the LP.
   void removeRow(SPxLP& lp, int i)
   {
      m_rIdx[i] = m_rIdx[lp.nRows()-1];
      lp.removeRow(i);
   }
   /// removes a column in the LP.
   void removeCol(SPxLP& lp, int j)
   {
      m_cIdx[j] = m_cIdx[lp.nCols()-1];
      lp.removeCol(j);
   }
   /// returns for a given row index of the (reduced) LP the corresponding row index in the unsimplified LP.
   int rIdx(int i) const
   {
      return m_rIdx[i];
   }
   /// returns for a given column index of the (reduced) LP the corresponding column index in the unsimplified LP.
   int cIdx(int j) const
   {
      return m_cIdx[j];
   }
   ///
   Real epsZero() const
   {
      return m_epsilon;
   }
   ///
   Real deltaBnd() const
   {
      return m_delta;
   }
   //@}

public:
   //------------------------------------
   //**@name Constructors / destructors */
   //@{
   /// default constructor.
   SPxMainSM() 
      : SPxSimplifier("MainSM")
      , m_stat(15)
   {}   
   /// copy constructor.
   SPxMainSM(const SPxMainSM& old) 
      : SPxSimplifier(old)
      , m_prim(old.m_prim)
      , m_slack(old.m_slack)
      , m_dual(old.m_dual)
      , m_redCost(old.m_redCost)
      , m_cBasisStat(old.m_cBasisStat)
      , m_rBasisStat(old.m_rBasisStat)
      , m_cIdx(old.m_cIdx)
      , m_rIdx(old.m_rIdx)
      , m_postsolved(old.m_postsolved)
      , m_epsilon(old.m_epsilon)
      , m_delta(old.m_delta)
      , m_stat(old.m_stat)
   {
      // copy pointers in m_hist
      m_hist.reSize(0);
      for(int k = 0; k < old.m_hist.size(); ++k)
      {
         if(old.m_hist[k] != 0)
            m_hist.append(old.m_hist[k]->clone());
         else
            m_hist.append(0);
      }
   }
   /// assignment operator
   SPxMainSM& operator=( const SPxMainSM& rhs)
   {
      if(this != &rhs)
      {
         SPxSimplifier::operator=(rhs);
         m_prim = rhs.m_prim;
         m_slack = rhs.m_slack;
         m_dual = rhs.m_dual;
         m_redCost = rhs.m_redCost;
         m_cBasisStat = rhs.m_cBasisStat;
         m_rBasisStat = rhs.m_rBasisStat;
         m_cIdx = rhs.m_cIdx;
         m_rIdx = rhs.m_rIdx;
         m_postsolved = rhs.m_postsolved;
         m_epsilon = rhs.m_epsilon;
         m_delta = rhs.m_delta;
         m_stat = rhs.m_stat;

         // delete pointers in m_hist
         for(int k = 0; k < m_hist.size(); ++k)
         {
            delete m_hist[k];
            m_hist[k] = 0;
         }

         m_hist.clear();

         // copy pointers in m_hist
         for(int k = 0; k < rhs.m_hist.size(); ++k)
         {
            if(rhs.m_hist[k] != 0)
               m_hist.append(rhs.m_hist[k]->clone());
            else
               m_hist.append(0);
         }
      }
      
      return *this;
   }   
   /// destructor.
   virtual ~SPxMainSM()
   {
      // delete pointers in m_hist
      for(int k = 0; k < m_hist.size(); ++k)
      {
         delete m_hist[k];
         m_hist[k] = 0;
      }
   }  
   /// clone function for polymorphism
   inline virtual SPxSimplifier* clone() const
   {
      return new SPxMainSM(*this);
   }
   //@}

   //------------------------------------
   //**@name LP simplification */
   //@{
   /// simplifies LP. 
   virtual Result simplify(SPxLP& lp, Real eps, Real delta);

   /// reconstructs an optimal solution for the unsimplified LP.
   virtual void unsimplify(const Vector& x, const Vector& y, const Vector& s, const Vector& r,
                           const SPxSolver::VarStatus rows[], const SPxSolver::VarStatus cols[]);

   /// specifies whether an optimal solution has already been unsimplified.
   virtual bool isUnsimplified() const
   {
      return m_postsolved;
   }
   /// returns a reference to the unsimplified primal solution.
   virtual const Vector& unsimplifiedPrimal()
   {
      assert(m_postsolved);
      return m_prim;
   }
   /// returns a reference to the unsimplified dual solution.
   virtual const Vector& unsimplifiedDual()
   {
      assert(m_postsolved);
      return m_dual;
   }
   /// returns a reference to the unsimplified slack values.
   virtual const Vector& unsimplifiedSlacks()
   {
      assert(m_postsolved);
      return m_slack;
   }
   /// returns a reference to the unsimplified reduced costs.
   virtual const Vector& unsimplifiedRedCost()
   {
      assert(m_postsolved);
      return m_redCost;
   }
   /// gets basis status for a single row.
   virtual SPxSolver::VarStatus getBasisRowStatus(int i) const
   {
      assert(m_postsolved);
      return m_rBasisStat[i];
   }
   /// gets basis status for a single column.
   virtual SPxSolver::VarStatus getBasisColStatus(int j) const
   {
      assert(m_postsolved);
      return m_cBasisStat[j];
   }
   /// get optimal basis.
   virtual void getBasis(SPxSolver::VarStatus rows[], SPxSolver::VarStatus cols[]) const
   {
      assert(m_postsolved);
      
      for(int i = 0; i < m_rBasisStat.size(); ++i)
         rows[i] = m_rBasisStat[i];

      for(int j = 0; j < m_rBasisStat.size(); ++j)
         cols[j] = m_cBasisStat[j];
   }   
   //@}

private:
   //------------------------------------
   //**@name Types */
   //@{
   /// 
   struct ElementCompare
   {
   public: 
      ElementCompare() {}
      
      int operator()(const SVector::Element& e1, const SVector::Element& e2) const
      {
	 if (EQ(e1.val, e2.val))
            return 0;
         if (e1.val < e2.val)
	    return -1;
	 else // (e1.val > e2.val)
	    return 1;
      }
   };
   //@}
};

} // namespace soplex
#endif // _SPXMAINSM_H_

//-----------------------------------------------------------------------------
//Emacs Local Variables:
//Emacs mode:c++
//Emacs c-basic-offset:3
//Emacs tab-width:8
//Emacs indent-tabs-mode:nil
//Emacs End:
//-----------------------------------------------------------------------------
