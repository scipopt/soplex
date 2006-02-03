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
#pragma ident "@(#) $Id: spxmainsm.h,v 1.2 2006/02/03 12:21:12 bzftuchs Exp $"

/**@file  spxmainsm.h
 * @brief General methods in LP preprocessing.
 */
#ifndef _SPXMAINSM_H_
#define _SPXMAINSM_H_

#include <assert.h>

#include "spxdefines.h"
#include "spxsimplifier.h"
#include "array.h"

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
   columns and bounds.  As Andersen and Andersen solely consider LPs
   with equality constraints, their methods have been slightly adapted
   for arbitrary constraint ranges.  Also infeasibility and
   unboundness may be detected.

   Removed are:
   - empty rows / columns
   - unconstraint rows
   - row singletons
   - forcing equality rows
   - zero objective column singletons
   - (implied) free column singletons
   - doubleton equations combined with a column singleton
   - (implicit) fixed columns
   - rows with all fixed variables due to implied bounds
   - redundant rhs / lhs
   - redundant column bounds
   - (weakly) dominated columns
   - columns with redundant bounds
   - dublicate rows / columns
*/
class SPxMainSM : public SPxSimplifier
{
private:
   //------------------------------------
   //**@name Types */
   //@{
   /// 
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

   class PostStep;
   
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
   DataArray<int>                  m_cIdx;       ///< removed column index vector in original LP.
   DataArray<int>                  m_rIdx;       ///< removed row index vector in original LP.
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
   
   /// removed empty rows and empty constraints.
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
      m_rIdx[i] = lp.rId(lp.nRows()-1).getIdx();
      lp.removeRow(i);
   }
   /// removes a column in the LP.
   void removeCol(SPxLP& lp, int j)
   {
      m_cIdx[j] = lp.cId(lp.nCols()-1).getIdx();
      lp.removeCol(j);
   }
   /// returns for a given row index of the given (reduced) LP the corresponding row index in the unsimplified LP.
   int rIdx(const SPxLP& lp, int i) const
   {
      return (m_rIdx[i] != -1) ? m_rIdx[i] : lp.rId(i).getIdx();
   }
   /// returns for a given column index of the given (reduced) LP the corresponding column index in the unsimplified LP.
   int cIdx(const SPxLP& lp, int j) const
   {
      return (m_cIdx[j] != -1) ? m_cIdx[j] : lp.cId(j).getIdx();
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
   {}   
   /// destructor.
   virtual ~SPxMainSM()
   {}  
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
         : m_eps(1e-7)
      {}
      /// destructor.
      virtual ~PostStep()
      {}
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

   class FreeConstraintPS : public PostStep
   {
   private:
      const int m_i;
      DSVector  m_row;
   
   public:
      ///
      FreeConstraintPS(const SPxLP& lp, const SPxMainSM& simplifier, int i)
         : m_i(simplifier.rIdx(lp, i))
         , m_row(lp.rowVector(i).size())
      {
         const SVector& row = lp.rowVector(i);
         
         for(int k = 0; k < row.size(); ++k)
            m_row.add(simplifier.cIdx(lp, row.index(k)), row.value(k));
      }
      ///
      virtual void execute(DVector& x, DVector& y, DVector& s, DVector& r,
                           DataArray<SPxSolver::VarStatus>& cBasis, DataArray<SPxSolver::VarStatus>& rBasis) const;
   };

   class EmptyConstraintPS : public PostStep
   {
   private:
      const int m_i;
      
   public:
      ///
      EmptyConstraintPS(int i)
         : m_i(i)
      {}
      ///
      virtual void execute(DVector& x, DVector& y, DVector& s, DVector& r,
                           DataArray<SPxSolver::VarStatus>& cBasis, DataArray<SPxSolver::VarStatus>& rBasis) const;
   };
   
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
      
      Real m_newLo;
      Real m_newUp;
      Real m_oldLo;
      Real m_oldUp;
      
   public:
      ///
      RowSingletonPS(const SPxLP& lp, const SPxMainSM& simplifier, int i, int j, bool strictLo, bool strictUp,
                     Real newLo, Real newUp, Real oldLo, Real oldUp)
         : m_i(simplifier.rIdx(lp, i))
         , m_j(simplifier.cIdx(lp, j))
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
            m_col.add(simplifier.rIdx(lp, col.index(k)), col.value(k));
      }
      ///
      virtual void execute(DVector& x, DVector& y, DVector& s, DVector& r,
                           DataArray<SPxSolver::VarStatus>& cBasis, DataArray<SPxSolver::VarStatus>& rBasis) const;
   };
 
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
         : m_i(simplifier.rIdx(lp, i))
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
            int j = simplifier.cIdx(lp, row.index(k));
            m_row.add(j, row.value(k));
            m_objs.add(j, lp.obj(row.index(k)));
            m_fixed[k] = EQrel(lp.lower(row.index(k)), lp.upper(row.index(k)));
            
            const SVector& col = lp.colVector(row.index(k));
            m_cols[k].setMax(col.size());
            
            for(int l = 0; l < col.size(); ++l)
               m_cols[k].add(simplifier.rIdx(lp, col.index(l)), col.value(l));
         }
      }
      ///
      virtual void execute(DVector& x, DVector& y, DVector& s, DVector& r,
                           DataArray<SPxSolver::VarStatus>& cBasis, DataArray<SPxSolver::VarStatus>& rBasis) const;
   };
   
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
         : m_j(simplifier.cIdx(lp, j))
         , m_val(val)
         , m_obj(lp.obj(j))
         , m_col(lp.colVector(j).size())
      {
         const SVector& col = lp.colVector(j);
         
         for(int k = 0; k < col.size(); ++k)
            m_col.add(simplifier.rIdx(lp, col.index(k)), col.value(k));
      }
      ///
      virtual void execute(DVector& x, DVector& y, DVector& s, DVector& r,
                           DataArray<SPxSolver::VarStatus>& cBasis, DataArray<SPxSolver::VarStatus>& rBasis) const;
   };

   class FixBoundsPS : public PostStep
   {
   private:
      const int            m_j;
      SPxSolver::VarStatus m_status;
      
   public:
      ///
      FixBoundsPS(const SPxLP& lp, const SPxMainSM& simplifier, int j, Real val)
         : m_j(simplifier.cIdx(lp, j))
      { 
         if (EQrel(lp.lower(j), lp.upper(j)))
            m_status = SPxSolver::FIXED;
         else if (EQrel(val, lp.lower(j)))
            m_status = SPxSolver::ON_LOWER;
         else if (EQrel(val, lp.upper(j)))
            m_status = SPxSolver::ON_UPPER;
         else if (lp.lower(j) <= -infinity && lp.upper(j) >= infinity)
            m_status = SPxSolver::ZERO;
         else
            abort();
      }
      ///
      virtual void execute(DVector& x, DVector& y, DVector& s, DVector& r,
                           DataArray<SPxSolver::VarStatus>& cBasis, DataArray<SPxSolver::VarStatus>& rBasis) const;
   };

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
         : m_j(simplifier.cIdx(lp, j))
         , m_bnd(loFree ? lp.upper(j) : lp.lower(j))
         , m_col(lp.colVector(j).size())
         , m_lRhs(lp.colVector(j).size())
         , m_rows(lp.colVector(j).size())
         , m_loFree(loFree)
      {
         const SVector& col = lp.colVector(j);
         
         for(int k = 0; k < col.size(); ++k)
         {
            int i = simplifier.rIdx(lp, col.index(k));

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
               m_rows[k].add(simplifier.cIdx(lp, row.index(l)), row.value(l));
         }
      }
      ///
      virtual void execute(DVector& x, DVector& y, DVector& s, DVector& r,
                           DataArray<SPxSolver::VarStatus>& cBasis, DataArray<SPxSolver::VarStatus>& rBasis) const;
   };

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
         : m_j(simplifier.cIdx(lp, j))
         , m_i(simplifier.rIdx(lp, i))
         , m_lhs(lp.lhs(i))
         , m_rhs(lp.rhs(i))
         , m_lower(lp.lower(j))
         , m_upper(lp.upper(j))
         , m_row(lp.rowVector(i).size())
      {
         const SVector& row = lp.rowVector(i);
      
         for(int k = 0; k < row.size(); ++k)
            m_row.add(simplifier.cIdx(lp, row.index(k)), row.value(k));
      }
      ///
      virtual void execute(DVector& x, DVector& y, DVector& s, DVector& r,
                           DataArray<SPxSolver::VarStatus>& cBasis, DataArray<SPxSolver::VarStatus>& rBasis) const;
   };
   
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
         : m_j(simplifier.cIdx(lp, j))
         , m_i(simplifier.rIdx(lp, i))
         , m_obj(lp.obj(j))
         , m_lRhs(slackVal)
         , m_onLhs(slackVal == lp.lhs(i))
         , m_eqCons(EQrel(lp.lhs(i), lp.rhs(i)))
         , m_row(lp.rowVector(i).size())
      {
         const SVector& row = lp.rowVector(i);
      
         for(int k = 0; k < row.size(); ++k)
            m_row.add(simplifier.cIdx(lp, row.index(k)), row.value(k));
      }
      ///
      virtual void execute(DVector& x, DVector& y, DVector& s, DVector& r,
                           DataArray<SPxSolver::VarStatus>& cBasis, DataArray<SPxSolver::VarStatus>& rBasis) const;
   };

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
         : m_j(simplifier.cIdx(lp, j))
         , m_k(simplifier.cIdx(lp, k))
         , m_i(simplifier.rIdx(lp, i))
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
            m_col.add(simplifier.rIdx(lp, col.index(l)), col.value(l));
      }
      ///
      virtual void execute(DVector& x, DVector& y, DVector& s, DVector& r,
                           DataArray<SPxSolver::VarStatus>& cBasis, DataArray<SPxSolver::VarStatus>& rBasis) const;
   };

   class ForceColumnPS : public PostStep
   {
   private:
      const int       m_j;
      const Real      m_lo;
      const Real      m_up;
      const Real      m_maxObj;
      DSVector        m_col;
      DSVector        m_lhs;
      DSVector        m_rhs;
      DSVector        m_y;
      Array<DSVector> m_rows;
      
   public:
      ///
      ForceColumnPS(const SPxLP& lp, const SPxMainSM& simplifier, int j, const DSVector& y)
         : m_j(simplifier.cIdx(lp, j))
         , m_lo(lp.lower(j))
         , m_up(lp.upper(j))
         , m_maxObj(lp.maxObj(j))
         , m_col(lp.colVector(j).size())
         , m_lhs(lp.colVector(j).size())
         , m_rhs(lp.colVector(j).size())
         , m_y(y)
         , m_rows(lp.colVector(j).size())
      {
         const SVector& col = lp.colVector(j);
         
         for(int k = 0; k < col.size(); ++k)
         {
            int i = simplifier.rIdx(lp, col.index(k));
            
            m_col.add(i, col.value(k));

	    assert(isNotZero(col.value(k)));
            
            m_lhs.add(i, lp.lhs(col.index(k)));
            m_rhs.add(i, lp.rhs(col.index(k)));
            
            const SVector& row = lp.rowVector(col.index(k));
            m_rows[k].setMax(row.size());
            
            for(int l = 0; l < row.size(); ++l)
               m_rows[k].add(simplifier.cIdx(lp, row.index(l)), row.value(l));
         }
      }
      ///
      virtual void execute(DVector& x, DVector& y, DVector& s, DVector& r,
                           DataArray<SPxSolver::VarStatus>& cBasis, DataArray<SPxSolver::VarStatus>& rBasis) const;
   };

   class DuplicateRowsPS : public PostStep
   {
   private:
      const int       m_i;
      const int       m_maxLhsIdx;
      const int       m_minRhsIdx;
      const bool      m_maxSense;
      DSVector        m_scale;
      DSVector        m_obj;
      Array<DSVector> m_cols;
      
   public:
      DuplicateRowsPS(const SPxLP& lp, const SPxMainSM& simplifier, int i,
                      int maxLhsIdx, int minRhsIdx, const DSVector& dupRows, const DataArray<double> scale)
         : m_i(simplifier.rIdx(lp, i))
         , m_maxLhsIdx((maxLhsIdx == -1) ? -1 : simplifier.rIdx(lp, maxLhsIdx))
         , m_minRhsIdx((minRhsIdx == -1) ? -1 : simplifier.rIdx(lp, minRhsIdx))
         , m_maxSense(lp.spxSense() == SPxLP::MAXIMIZE)
         , m_scale(dupRows.size())
         , m_obj(lp.rowVector(i).size())
         , m_cols(lp.rowVector(i).size())
      {
         Real rowScale = scale[i];
         
         for(int k = 0; k < dupRows.size(); ++k)
            m_scale.add(simplifier.rIdx(lp, dupRows.index(k)), rowScale / scale[dupRows.index(k)]);
         
         const SVector& row = lp.rowVector(i);
         
         for(int k = 0; k < row.size(); ++k)
         {
            m_obj.add(simplifier.cIdx(lp, row.index(k)), lp.obj(row.index(k)));
            
            const SVector& col = lp.colVector(row.index(k));
            m_cols[k].setMax(col.size());
            
            for(int l = 0; l < col.size(); ++l)
               m_cols[k].add(simplifier.rIdx(lp, col.index(l)), col.value(l));
         }         
      }
      virtual void execute(DVector& x, DVector& y, DVector& s, DVector& r,
                           DataArray<SPxSolver::VarStatus>& cBasis, DataArray<SPxSolver::VarStatus>& rBasis) const;
   };

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
         : m_j(simplifier.cIdx(lp, j))
         , m_k(simplifier.cIdx(lp, k))
         , m_loJ(lp.lower(j))
         , m_upJ(lp.upper(j))
         , m_loK(lp.lower(k))
         , m_upK(lp.upper(k))
         , m_scale(scale)
      {}
      virtual void execute(DVector& x, DVector& y, DVector& s, DVector& r,
                           DataArray<SPxSolver::VarStatus>& cBasis, DataArray<SPxSolver::VarStatus>& rBasis) const;
   };

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
