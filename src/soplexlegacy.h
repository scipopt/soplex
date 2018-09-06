/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the class library                   */
/*       SoPlexBase --- the Sequential object-oriented simPlex.                  */
/*                                                                           */
/*    Copyright (C) 1996-2017 Konrad-Zuse-Zentrum                            */
/*                            fuer Informationstechnik Berlin                */
/*                                                                           */
/*  SoPlexBase is distributed under the terms of the ZIB Academic Licence.       */
/*                                                                           */
/*  You should have received a copy of the ZIB Academic License              */
/*  along with SoPlexBase; see the file COPYING. If not email to soplex@zib.de.  */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/**@file  soplexlegacy.h
 * @brief preconfigured \ref soplex::SoPlexBaseLegacy "SoPlexBaseLegacy" LP-solver.
 */
#define _SOPLEXLEGACY_H_ // flipped the place to remove the legacy parts from the compiliation
#ifndef _SOPLEXLEGACY_H_


#include <assert.h>
#include <string.h>

#include "spxsolver.h"
#include "spxscaler.h"
#include "spxsimplifier.h"
#include "spxsteeppr.h"
#include "spxboundflippingrt.h"
#include "spxfileio.h"
#include "spxweightst.h"
#include "slufactor.h"

namespace soplex
{

/**@class SoPlexBaseLegacy
   @brief   Preconfigured SoPlexBaseLegacy LP-solver.
   @ingroup Algo
*/
class SoPlexBaseLegacy : public SPxLP
{
protected:

   //-------------------------
   //**@name Protected data */
   //@{
   //   SPxWeightST st;  ///< weight starter
   SLUFactor       m_slu;        ///< LU Factorisation
   SPxSolverBase       m_solver;     ///< solver
   SPxScaler*      m_postScaler; ///< post-scaler
   SPxSimplifier*  m_simplifier; ///< simplifier
   bool            m_vanished;   ///< did the presolver solve the problem ?
   bool            m_freePostScaler;  ///< true iff m_postScaler should be freed inside of this object
   bool            m_freeSimplifier;  ///< true iff m_simplifier should be freed inside of this object
   DataArray<SPxSolverBase::VarStatus> m_colsbasisstatus;
   DataArray<SPxSolverBase::VarStatus> m_rowsbasisstatus;
   //@}

public:

   //---------------------------------------
   //**@name Construction / destruction */
   //@{
   /// default construtor.
   explicit SoPlexBaseLegacy(
      SPxOut&                   outstream,
      SPxSolverBase::Type           type = SPxSolverBase::LEAVE,
      SPxSolverBase::Representation rep  = SPxSolverBase::COLUMN );
   virtual ~SoPlexBaseLegacy();
   /// assignment operator.
   SoPlexBaseLegacy& operator=(const SoPlexBaseLegacy& rhs);
   /// copy constructor.
   SoPlexBaseLegacy(const SoPlexBaseLegacy&);
   //@}

   /// message handler
//   mutable SPxOut spxout;

   //---------------------------------------
   //**@name Access / modification */
   //@{

   /// set verbosity
//   virtual void setVerbosity(const SPxOut::Verbosity v)
//   {
//      spxout.setVerbosity(v);
//      setOutstream(spxout);
//   }

   /// set update type for factorization.
   virtual void setUtype(SLUFactor::UpdateType tp)
   {
      m_slu.setUtype(tp);
   }
   /// return current Pricing.
   inline SPxSolverBase::Pricing pricing() const
   {
      return m_solver.pricing();
   }
   /// set FULL or PARTIAL pricing.
   virtual void setPricing(SPxSolverBase::Pricing pr)
   {
      m_solver.setPricing(pr);
   }
   /// return current Type.
   inline SPxSolverBase::Type type() const
   {
      return m_solver.type();
   }
   /// return current basis representation.
   inline SPxSolverBase::Representation rep() const
   {
      return m_solver.rep();
   }
   /// set LEAVE or ENTER algorithm.
   virtual void setType(SPxSolverBase::Type tp)
   {
      m_solver.setType(tp);
   }
   /// set ROW or COLUMN representation.
   virtual void setRep (SPxSolverBase::Representation p_rep)
   {
      m_solver.setRep(p_rep);
   }
   /// setup postscaler to use. If \p destroy is true, \p scaler will be freed in destructor.
   virtual void setPostScaler(SPxScaler* scaler, const bool destroy = false);
   /// setup simplifier to use. If \p destroy is true, \p simpli will be freed in destructor.
   virtual void setSimplifier(SPxSimplifier* simpli, const bool destroy = false);
   /// has a simplifier been set?
   inline bool has_simplifier() const
   {
      return m_simplifier != 0;
   }
   /// has a postscaler been set?
   inline bool has_postscaler() const
   {
      return m_postScaler != 0;
   }
   /// setup pricer to use.
   virtual void setPricer(SPxPricer* pricer, const bool destroy = false)
   {
      m_solver.setPricer(pricer, destroy);
   }
   /// setup ratio-tester to use.
   virtual void setTester(SPxRatioTester* tester, const bool destroy = false)
   {
      m_solver.setTester(tester, destroy);
   }
   /// setup starting basis generator to use.
   virtual void setStarter(SPxStarter* starter, const bool destroy = false)
   {
      m_solver.setStarter(starter, destroy);
   }
   /// @throw SPxStatusException if simplifier loaded, this is not yet implemented
   /// set starting basis
   virtual void setBasis(SPxSolverBase::VarStatus rows[], SPxSolverBase::VarStatus cols[])
   {
      if (has_simplifier())
      {
         MSG_ERROR( std::cerr << "ESOLVR04 setting starting basis with presolving not yet implemented" << std::endl; )
            throw SPxStatusException("XSOLVR04 setting starting basis with presolving not yet implemented");
      }

      m_colsbasisstatus.reSize(nCols());
      for (int i = 0; i < nCols(); i++)
         m_colsbasisstatus[i] = cols[i];

      m_rowsbasisstatus.reSize(nRows());
      for (int i = 0; i < nRows(); i++)
         m_rowsbasisstatus[i] = rows[i];
   }
   /// clear starting basis
   virtual void clearBasis()
   {
      m_colsbasisstatus.clear();
      m_rowsbasisstatus.clear();
      m_solver.reLoad();
   }
   /// set time limit.
   virtual void setTerminationTime(Real time = infinity)
   {
      m_solver.setTerminationTime(time);
   }
   /// return time limit.
   inline Real terminationTime() const
   {
      return m_solver.terminationTime();
   }
   /// set iteration limit.
   virtual void setTerminationIter(int iter = -1)
   {
      m_solver.setTerminationIter(iter);
   }
   /// return iteration limit.
   inline int terminationIter() const
   {
      return m_solver.terminationIter();
   }
   /// set objective limit.
   virtual void setTerminationValue(Real val = infinity)
   {
      m_solver.setTerminationValue(val);
   }
   /// return objective limit.
   inline Real terminationValue() const
   {
      return m_solver.terminationValue();
   }
   /// allowed primal feasibility tolerance.
   virtual Real feastol() const
   {
      return m_solver.feastol();
   }
   /// allowed optimality, i.e., dual feasibility tolerance.
   virtual Real opttol() const
   {
      return m_solver.opttol();
   }
   /// guaranteed primal and dual bound violation for optimal solution, returning the maximum of feastol() and opttol(), i.e., the less tight tolerance.
   virtual Real delta() const
   {
      return m_solver.delta();
   }
   /// set parameter \p feastol.
   virtual void setFeastol(Real d)
   {
      m_solver.setFeastol(d);
   }
   /// set parameter \p opttol.
   virtual void setOpttol(Real d)
   {
      m_solver.setOpttol(d);
   }
   /// set parameter \p delta, i.e., set \p feastol and \p opttol to same value.
   virtual void setDelta(Real d)
   {
      m_solver.setDelta(d);
   }
   //@}

   //---------------------------------------
   //**@name Solving and solution query */
   //@{
   /// @throw SPxStatusException if no problem loaded
   virtual SPxSolverBase::Status solve();
   ///
   virtual Real objValue() const;
   ///
   virtual SPxSolverBase::Status getPrimal(Vector& vector) const;
   ///
   virtual SPxSolverBase::Status getSlacks(Vector& vector) const;
   ///
   virtual SPxSolverBase::Status getDual(Vector& vector) const;
   ///
   virtual SPxSolverBase::Status getRedCost(Vector& vector) const;

   /// gets basis status for a single row.
   SPxSolverBase::VarStatus getBasisRowStatus(int row) const;

   /// gets basis status for a single column.
   SPxSolverBase::VarStatus getBasisColStatus(int col) const;

   /// get current basis, and return solver status.
   SPxSolverBase::Status getBasis(SPxSolverBase::VarStatus rows[], SPxSolverBase::VarStatus cols[]) const;

   const char* getColName(
      int            idx,
      const NameSet* cnames,
      char*          buf)
   {
      assert(buf != 0);
      assert(idx >= 0);
      assert(idx < nCols());

      if (cnames != 0)
      {
         DataKey key = cId(idx);

         if (cnames->has(key))
            return (*cnames)[key];
      }
      spxSnprintf(buf, SPX_MAXSTRLEN, "x%d", idx);

      return buf;
   }

   const char* getRowName(
      int            idx,
      const NameSet* rnames,
      char*          buf)
   {
      assert(buf != 0);
      assert(idx >= 0);
      assert(idx < nRows());

      if (rnames != 0)
      {
         DataKey key = rId(idx);

         if (rnames->has(key))
            return (*rnames)[key];
      }
      spxSnprintf(buf, SPX_MAXSTRLEN, "C%d", idx);

      return buf;
   }

   /// @throw SPxStatusException if simplifier loaded, this is not yet
   /// implemented
   virtual SPxSolverBase::Status getPrimalray(Vector& vector) const;

   /// @throw SPxStatusException if simplifier loaded, this is not yet
   /// implemented
   virtual SPxSolverBase::Status getDualfarkas(Vector& vector) const;

   /// get violation of constraints.
   virtual void qualConstraintViolation(Real& maxviol, Real& sumviol) const;
   /// get violations of bounds.
   virtual void qualBoundViolation(Real& maxviol, Real& sumviol) const;
#if 0
   /// get the residuum |Ax-b|.
   virtual void qualSlackViolation(Real& maxviol, Real& sumviol) const;
   /// get violation of optimality criterion.
   virtual void qualRedCostViolation(Real& maxviol, Real& sumviol) const;
#endif
   /// time spent in factorizations
   virtual Real getFactorTime() const
   {
      return m_vanished ? REAL(0.0) : m_slu.getFactorTime();
   }
   /// number of factorizations performed
   virtual int getFactorCount() const
   {
      return m_vanished ? 0 : m_slu.getFactorCount();
   }
   /// time spent in solves
   virtual Real getSolveTime() const
   {
      return m_vanished ? REAL(0.0) : m_slu.getSolveTime();
   }
   /// number of solves performed
   virtual int getSolveCount() const
   {
      return m_vanished ? 0 : m_slu.getSolveCount();
   }
   ///
   virtual int iteration() const
   {
      return m_vanished ? 0 : m_solver.basis().iteration();
   }
   ///
   virtual bool terminate()
   {
      return m_solver.terminate();
   }
   /// returns the current status
   virtual SPxSolverBase::Status status() const
   {
      if (m_vanished)
         return SPxSolverBase::OPTIMAL;

      return m_solver.status();
   }
   //@}

   //---------------------------------------
   //**@name I/O */
   //@{

   /** Load basis from \p filename in MPS format. If \p rowNames and \p
    *  colNames are \c NULL, default names are used for the constraints and
    *  variables.
    */
   virtual bool readBasisFile(const char* filename,
      const NameSet* rowNames, const NameSet* colNames);

   /** Write basis to \p filename in MPS format. If \p rowNames and \p
    *  colNames are \c NULL, default names are used for the constraints and
    *  variables.
    */
   virtual bool writeBasisFile(const char* filename,
      const NameSet* rowNames, const NameSet* colNames);

   /** Write LP, basis and parameter settings of the current SPxSolverBase object
    *  (i.e. after simplifying and scaling).
    *  LP is written in MPS format to "\p filename".mps, basis is written in
    *  "\p filename".bas, and parameters are written to "\p filename".set.
    *  If \p rowNames and \p colNames are \c NULL, default names are used for
    *  the constraints and variables.
    */
   virtual bool writeState(const char* filename,
      const NameSet* rowNames = NULL, const NameSet* colNames = NULL) const;

   /// returns statistical information in form of a string.
   std::string statistics() const
   {
      return m_solver.statistics();
   }
   //@}

private:

   //------------------------------------
   //**@name Private helpers */
   //@{
   /// undoes preprocessing such that the unsimplified solution values and basis is available
   void unsimplify() const;
   //@}
};
} // namespace soplex
#endif // _SOPLEX_H_
