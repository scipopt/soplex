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

#include <assert.h>
#include <iostream>
#include <sstream>

#include "soplex/spxdefines.h"
#include "soplex.h"
#include "soplex/spxpricer.h"
#include "soplex/spxratiotester.h"
#include "soplex/spxstarter.h"
#include "soplex/spxout.h"
#include "soplex/timerfactory.h"

namespace soplex
{
/// definition of signature to prevent instantiation before def. error
template <>
void SPxSolverBase<Real>::resetClockStats();

template <>
bool SPxSolverBase<Real>::isConsistent() const;

template <>
void SPxSolverBase<Real>::setPricing(Pricing pr);

template <>
void SPxSolverBase<Real>::reDim();

template <>
void SPxSolverBase<Real>::clear();

template <>
void SPxSolverBase<Real>::clearUpdateVecs(void);

template <>
void SPxSolverBase<Real>::setFeastol(Real d);

template <>
void SPxSolverBase<Real>::setOpttol(Real d);

template <>
bool SPxSolverBase<Real>::read(std::istream& in, NameSet* rowNames,
                               NameSet* colNames, DIdxSet* intVars)
{
   if(initialized)
   {
      clear();
      unInit();

      if(thepricer)
         thepricer->clear();

      if(theratiotester)
         theratiotester->clear();
   }

   this->unLoad();

   if(!SPxLP::read(in, rowNames, colNames, intVars))
      return false;

   this->theLP = this;

   return true;
}

template <>
void SPxSolverBase<Real>::reLoad()
{
   forceRecompNonbasicValue();
   unInit();
   this->unLoad();
   this->theLP = this;
   m_status = SPxSolverBase<Real>::UNKNOWN;

   if(thepricer)
      thepricer->clear();

   if(theratiotester)
      theratiotester->clear();
}

template <>
void SPxSolverBase<Real>::loadLP(const SPxLP& lp, bool initSlackBasis)
{
   clear();
   unInit();
   this->unLoad();
   resetClockStats();

   if(thepricer)
      thepricer->clear();

   if(theratiotester)
      theratiotester->clear();

   SPxLP::operator=(lp);
   reDim();
   SPxBasisBase<Real>::load(this, initSlackBasis);
}

template <>
void SPxSolverBase<Real>::setBasisSolver(SLinSolver* slu, const bool destroy)
{
   // we need to set the outstream before we load the solver to ensure that the basis
   // can be initialized with this pointer in loadSolver()
   assert(spxout != 0);
   slu->spxout = spxout;
   SPxBasisBase<Real>::loadBasisSolver(slu, destroy);
}

template <>
void SPxSolverBase<Real>::loadBasis(const typename SPxBasisBase<Real>::Desc& p_desc)
{
   unInit();

   if(SPxBasisBase<Real>::status() == SPxBasisBase<Real>::NO_PROBLEM)
   {
      SPxBasisBase<Real>::load(this, false);
   }

   setBasisStatus(SPxBasisBase<Real>::REGULAR);
   SPxBasisBase<Real>::loadDesc(p_desc);
}

template <>
void SPxSolverBase<Real>::setPricer(SPxPricer<Real>* x, const bool destroy)
{

   assert(!freePricer || thepricer != 0);

   if(freePricer)
   {
      delete thepricer;
      thepricer = 0;
   }

   if(x != 0 && x != thepricer)
   {
      setPricing(FULL);

      if(isInitialized())
         x->load(this);
      else
         x->clear();
   }

   if(thepricer && thepricer != x)
      thepricer->clear();

   thepricer = x;

   freePricer = destroy;
}

template <>
void SPxSolverBase<Real>::setTester(SPxRatioTester<Real>* x, const bool destroy)
{
   assert(!freeRatioTester || theratiotester != 0);

   if(freeRatioTester)
   {
      delete theratiotester;
      theratiotester = 0;
   }

   theratiotester = x;

   // set the solver pointer inside the ratiotester
   if(theratiotester != 0)
   {
      if(isInitialized())
         theratiotester->load(this);
      else
         theratiotester->clear();
   }

   freeRatioTester = destroy;
}

template <>
void SPxSolverBase<Real>::setStarter(SPxStarter<Real>* x, const bool destroy)
{

   assert(!freeStarter || thestarter != 0);

   if(freeStarter)
   {
      delete thestarter;
      thestarter = 0;
   }

   thestarter = x;

   freeStarter = destroy;
}

template <>
void SPxSolverBase<Real>::setType(Type tp)
{

   if(theType != tp)
   {
      theType = tp;

      forceRecompNonbasicValue();

      unInit();
#if 0
      else
      {
         if(!matrixIsSetup)
         {
            SPxBasisBase<Real>::load(this);
            // SPxBasisBase<Real>::load(desc());
            // not needed, because load(this) allready loads descriptor
         }

         factorized = false;
         m_numCycle = 0;
#endif
         MSG_INFO3((*spxout), (*spxout) << "Switching to "
                   << static_cast<const char*>((tp == LEAVE)
                                               ? "leaving" : "entering")
                   << " algorithm" << std::endl;)
      }
   }

   template <>
   void SPxSolverBase<Real>::initRep(Representation p_rep)
   {

      Real tmpfeastol = feastol();
      Real tmpopttol = opttol();

      theRep = p_rep;

      if(theRep == COLUMN)
      {
         thevectors   = this->colSet();
         thecovectors = this->rowSet();
         theFrhs      = &primRhs;
         theFvec      = &primVec;
         theCoPrhs    = &dualRhs;
         theCoPvec    = &dualVec;
         thePvec      = &addVec;
         theRPvec     = theCoPvec;
         theCPvec     = thePvec;
         theUbound    = &theUCbound;
         theLbound    = &theLCbound;
         theCoUbound  = &theURbound;
         theCoLbound  = &theLRbound;
      }
      else
      {
         assert(theRep == ROW);

         thevectors   = this->rowSet();
         thecovectors = this->colSet();
         theFrhs      = &dualRhs;
         theFvec      = &dualVec;
         theCoPrhs    = &primRhs;
         theCoPvec    = &primVec;
         thePvec      = &addVec;
         theRPvec     = thePvec;
         theCPvec     = theCoPvec;
         theUbound    = &theURbound;
         theLbound    = &theLRbound;
         theCoUbound  = &theUCbound;
         theCoLbound  = &theLCbound;
      }

      unInit();
      reDim();

      forceRecompNonbasicValue();

      setFeastol(tmpfeastol);
      setOpttol(tmpopttol);

      SPxBasisBase<Real>::setRep();

      if(SPxBasisBase<Real>::status() > SPxBasisBase<Real>::NO_PROBLEM)
         SPxBasisBase<Real>::loadDesc(this->desc());

      if(thepricer && thepricer->solver() == this)
         thepricer->setRep(p_rep);
   }

   template <>
   void SPxSolverBase<Real>::setRep(Representation p_rep)
   {

      if(p_rep != theRep)
         initRep(p_rep);
   }

   // needed for strongbranching. use carefully
   template <>
   void SPxSolverBase<Real>::reinitializeVecs()
   {

      initialized = true;

      if(type() == ENTER)
      {
         if(rep() == COLUMN)
            setPrimalBounds();
         else
            setDualRowBounds();

         setEnterBounds();
         computeEnterCoPrhs();
      }
      else
      {
         if(rep() == ROW)
            setPrimalBounds();
         else
            setDualColBounds();

         setLeaveBounds();
         computeLeaveCoPrhs();
      }

      SPxBasisBase<Real>::coSolve(*theCoPvec, *theCoPrhs);
      computePvec();
      computeFrhs();
      SPxBasisBase<Real>::solve(*theFvec, *theFrhs);

      theShift  = 0.0;
      lastShift = 0.0;

      if(type() == ENTER)
      {
         computeCoTest();
         computeTest();
      }
      else
      {
         computeFtest();
      }

      assert((testBounds(), 1));
   }

   template <>
   void SPxSolverBase<Real>::resetClockStats()
   {
      nClckSkipsLeft = 0;
      nCallsToTimelim = 0;
      theCumulativeTime = 0.0;
   }

   template <>
   void SPxSolverBase<Real>::init()
   {

      assert(thepricer      != 0);
      assert(theratiotester != 0);

      if(!initialized)
      {
         initialized = true;
         clearUpdateVecs();
         reDim();

         if(SPxBasisBase<Real>::status() <= SPxBasisBase<Real>::NO_PROBLEM || this->solver() != this)
            SPxBasisBase<Real>::load(this);

         initialized = false;
      }

      if(!this->matrixIsSetup)
         SPxBasisBase<Real>::loadDesc(this->desc());

      // Inna/Tobi: don't "upgrade" a singular basis to a regular one
      if(SPxBasisBase<Real>::status() == SPxBasisBase<Real>::SINGULAR)
         return;

      // catch pathological case for LPs with zero constraints
      if(dim() == 0)
      {
         this->factorized = true;
      }

      // we better factorize explicitly before solving
      if(!this->factorized)
      {
         try
         {
            SPxBasisBase<Real>::factorize();
         }
         catch(const SPxException&)
         {
            // reload inital slack basis in case the factorization failed
            assert(SPxBasisBase<Real>::status() <= SPxBasisBase<Real>::SINGULAR);
            SPxBasisBase<Real>::restoreInitialBasis();
            SPxBasisBase<Real>::factorize();
            assert(this->factorized);
         }
      }

      m_numCycle = 0;

      if(type() == ENTER)
      {
         if(rep() == COLUMN)
         {
            setPrimalBounds();
            setBasisStatus(SPxBasisBase<Real>::PRIMAL);
         }
         else
         {
            setDualRowBounds();
            setBasisStatus(SPxBasisBase<Real>::DUAL);
         }

         setEnterBounds();
         computeEnterCoPrhs();
         // prepare support vectors for sparse pricing
         infeasibilities.setMax(dim());
         infeasibilitiesCo.setMax(coDim());
         isInfeasible.reSize(dim());
         isInfeasibleCo.reSize(coDim());
         theratiotester->setDelta(entertol());
      }
      else
      {
         if(rep() == ROW)
         {
            setPrimalBounds();
            setBasisStatus(SPxBasisBase<Real>::PRIMAL);
         }
         else
         {
            setDualColBounds();
            setBasisStatus(SPxBasisBase<Real>::DUAL);
         }

         setLeaveBounds();
         computeLeaveCoPrhs();
         // prepare support vectors for sparse pricing
         infeasibilities.setMax(dim());
         isInfeasible.reSize(dim());
         theratiotester->setDelta(leavetol());
      }

      SPxBasisBase<Real>::coSolve(*theCoPvec, *theCoPrhs);
      computePvec();
      computeFrhs();
      SPxBasisBase<Real>::solve(*theFvec, *theFrhs);

      theShift = 0.0;

      if(type() == ENTER)
      {
         shiftFvec();
         lastShift = theShift + entertol();

         computeCoTest();
         computeTest();
      }
      else
      {
         shiftPvec();
         lastShift = theShift + leavetol();

         computeFtest();
      }

      if(!initialized)
      {
         // if(thepricer->solver() != this)
         thepricer->load(this);
         // if(theratiotester->solver() != this)
         theratiotester->load(this);
         initialized = true;
      }
   }

   template <>
   void SPxSolverBase<Real>::setPricing(Pricing pr)
   {
      thePricing = pr;

      if(initialized && type() == ENTER)
      {
         computePvec();
         computeCoTest();
         computeTest();
      }
   }

   template <>
   void SPxSolverBase<Real>::setDecompStatus(DecompStatus decomp_stat)
   {
      if(decomp_stat == FINDSTARTBASIS)
         getStartingDecompBasis = true;
      else
         getStartingDecompBasis = false;
   }

   /*
     The following method resizes all vectors and arrays of |SoPlex|
     (excluding inherited vectors).
   */
   template <>
   void SPxSolverBase<Real>::reDim()
   {

      int newsize = SPxLP::nCols() > SPxLP::nRows() ? SPxLP::nCols() : SPxLP::nRows();

      if(newsize > unitVecs.size())
      {
         unitVecs.reSize(newsize);

         while(newsize-- > 0)
            unitVecs[newsize] = UnitVector(newsize);
      }

      if(isInitialized())
      {
         theFrhs->reDim(dim());
         theFvec->reDim(dim());
         thePvec->reDim(coDim());

         theCoPrhs->reDim(dim());
         theCoPvec->reDim(dim());

         theTest.reDim(coDim());
         theCoTest.reDim(dim());

         theURbound.reDim(SPxLP::nRows());
         theLRbound.reDim(SPxLP::nRows());
         theUCbound.reDim(SPxLP::nCols());
         theLCbound.reDim(SPxLP::nCols());
         theUBbound.reDim(dim());
         theLBbound.reDim(dim());
      }
   }

   template <>
   void SPxSolverBase<Real>::clear()
   {
      unitVecs.reSize(0);

      dualRhs.clear();
      dualVec.clear();
      primRhs.clear();
      primVec.clear();
      addVec.clear();
      theURbound.clear();
      theLRbound.clear();
      theUCbound.clear();
      theLCbound.clear();
      theTest.clear();
      theCoTest.clear();

      forceRecompNonbasicValue();
      unInit();
      SPxLP::clear();
      setBasisStatus(SPxBasisBase<Real>::NO_PROBLEM);

      // clear the basis only when theLP is present, because LP data (nrows, ncols) is used in reDim()
      if(this->theLP != 0)
         SPxBasisBase<Real>::reDim();

      infeasibilities.clear();
      infeasibilitiesCo.clear();
      isInfeasible.clear();
      isInfeasibleCo.clear();
   }

   template <>
   void SPxSolverBase<Real>::unscaleLPandReloadBasis()
   {
      SPxLPBase<Real>::unscaleLP();
      SPxBasisBase<Real>::invalidate();
      unInit();
      init();
   }

   template <>
   void SPxSolverBase<Real>::clearUpdateVecs(void)
   {
      theFvec->clearUpdate();
      thePvec->clearUpdate();
      theCoPvec->clearUpdate();
      solveVector2 = 0;
      solveVector3 = 0;
      coSolveVector2 = 0;
      coSolveVector3 = 0;
   }

   /*
     When the basis matrix factorization is recomputed from scratch,
     we also recompute the vectors.
   */
   template <>
   void SPxSolverBase<Real>::factorize()
   {

      MSG_INFO3((*spxout), (*spxout) << " --- refactorizing basis matrix" << std::endl;)

      try
      {
         SPxBasisBase<Real>::factorize();
      }
      catch(const SPxStatusException&)
      {
         assert(SPxBasisBase<Real>::status() == SPxBasisBase<Real>::SINGULAR);
         m_status = SINGULAR;
         std::stringstream s;
         s << "Basis is singular (numerical troubles, feastol = " << feastol() << ", opttol = " << opttol()
           << ")";
         throw SPxStatusException(s.str());
      }

      if(!initialized)
      {
         init();  // not sure if init() is neccessary here
         // we must not go on here because not all vectors (e.g. fVec) may be set up correctly
         return;
      }

      if(SPxBasisBase<Real>::status() >= SPxBasisBase<Real>::REGULAR)
      {
#ifndef NDEBUG
         DVector ftmp(fVec());
         DVector ptmp(pVec());
         DVector ctmp(coPvec());
#endif  // NDEBUG

         if(type() == LEAVE)
         {
            /* we have to recompute theFrhs, because roundoff errors can occur during updating, especially when
             * columns/rows with large bounds are present
             */
            computeFrhs();
            SPxBasisBase<Real>::solve(*theFvec, *theFrhs);
            SPxBasisBase<Real>::coSolve(*theCoPvec, *theCoPrhs);

#ifndef NDEBUG
            ftmp -= fVec();
            ptmp -= pVec();
            ctmp -= coPvec();

            if(ftmp.length() > DEFAULT_BND_VIOL)
            {
               MSG_DEBUG(std::cout << "DSOLVE21 fVec:   " << ftmp.length() << std::endl;)
               ftmp = fVec();
               multBaseWith(ftmp);
               ftmp -= fRhs();

               if(ftmp.length() > DEFAULT_BND_VIOL)
                  MSG_INFO1((*spxout), (*spxout) << "ESOLVE29 " << iteration() << ": fVec error = "
                            << ftmp.length() << " exceeding DEFAULT_BND_VIOL = " << DEFAULT_BND_VIOL << std::endl;)
               }

            if(ctmp.length() > DEFAULT_BND_VIOL)
            {
               MSG_DEBUG(std::cout << "DSOLVE23 coPvec: " << ctmp.length() << std::endl;)
               ctmp = coPvec();
               multWithBase(ctmp);
               ctmp -= coPrhs();

               if(ctmp.length() > DEFAULT_BND_VIOL)
                  MSG_INFO1((*spxout), (*spxout) << "ESOLVE30 " << iteration() << ": coPvec error = "
                            << ctmp.length() << " exceeding DEFAULT_BND_VIOL = " << DEFAULT_BND_VIOL << std::endl;)
               }

            if(ptmp.length() > DEFAULT_BND_VIOL)
            {
               MSG_DEBUG(std::cout << "DSOLVE24 pVec:   " << ptmp.length() << std::endl;)
            }

#endif  // NDEBUG

            computeFtest();
         }
         else
         {
            assert(type() == ENTER);

            SPxBasisBase<Real>::coSolve(*theCoPvec, *theCoPrhs);
            computeCoTest();

            if(pricing() == FULL)
            {
               /* to save time only recompute the row activities (in row rep) when we are already nearly optimal to
                * avoid missing any violations from previous updates */
               if(rep() == ROW && m_pricingViolCo < entertol() && m_pricingViol < entertol())
                  computePvec();

               /* was deactivated, but this leads to warnings in testVecs() */
               computeTest();
            }
         }
      }

#ifdef ENABLE_ADDITIONAL_CHECKS

      /* moved this test after the computation of fTest and coTest below, since these vectors might not be set up at top, e.g. for an initial basis */
      if(SPxBasisBase<Real>::status() > SPxBasisBase<Real>::SINGULAR)
         testVecs();

#endif
   }

   /* We compute how much the current solution violates (primal or dual) feasibility. In the
      row/enter or column/leave algorithm the maximum violation of dual feasibility is
      computed. In the row/leave or column/enter algorithm the primal feasibility is checked.
      Additionally, the violation from pricing is taken into account. */
   template <>
   Real SPxSolverBase<Real>::maxInfeas() const
   {
      Real inf = 0.0;

      if(type() == ENTER)
      {
         if(m_pricingViolUpToDate && m_pricingViolCoUpToDate)
            inf = m_pricingViol + m_pricingViolCo;

         for(int i = 0; i < dim(); i++)
         {
            if((*theFvec)[i] > theUBbound[i])
               inf = MAXIMUM(inf, (*theFvec)[i] - theUBbound[i]);
            else if((*theFvec)[i] < theLBbound[i])
               inf = MAXIMUM(inf, theLBbound[i] - (*theFvec)[i]);
         }
      }
      else
      {
         assert(type() == LEAVE);

         if(m_pricingViolUpToDate)
            inf = m_pricingViol;

         for(int i = 0; i < dim(); i++)
         {
            if((*theCoPvec)[i] > (*theCoUbound)[i])
               inf = MAXIMUM(inf, (*theCoPvec)[i] - (*theCoUbound)[i]);
            else if((*theCoPvec)[i] < (*theCoLbound)[i])
               inf = MAXIMUM(inf, (*theCoLbound)[i] - (*theCoPvec)[i]);
         }

         for(int i = 0; i < coDim(); i++)
         {
            if((*thePvec)[i] > (*theUbound)[i])
               inf = MAXIMUM(inf, (*thePvec)[i] - (*theUbound)[i]);
            else if((*thePvec)[i] < (*theLbound)[i])
               inf = MAXIMUM(inf, (*theLbound)[i] - (*thePvec)[i]);
         }
      }

      return inf;
   }

   /* check for (dual) violations above tol and immediately return false w/o checking the remaining values
      This method is useful for verifying whether an objective limit can be used as termination criterion */
   template <>
   bool SPxSolverBase<Real>::noViols(Real tol) const
   {
      assert(tol >= 0.0);

      if(type() == ENTER)
      {
         for(int i = 0; i < dim(); i++)
         {
            if((*theFvec)[i] - theUBbound[i] > tol)
               return false;

            if(theLBbound[i] - (*theFvec)[i] > tol)
               return false;
         }
      }
      else
      {
         assert(type() == LEAVE);

         for(int i = 0; i < dim(); i++)
         {
            if((*theCoPvec)[i] - (*theCoUbound)[i] > tol)
               return false;

            if((*theCoLbound)[i] - (*theCoPvec)[i] > tol)
               return false;
         }

         for(int i = 0; i < coDim(); i++)
         {
            if((*thePvec)[i] - (*theUbound)[i] > tol)
               return false;

            if((*theLbound)[i] - (*thePvec)[i] > tol)
               return false;
         }
      }

      return true;
   }

   template <>
   Real SPxSolverBase<Real>::nonbasicValue()
   {
      int i;
      StableSum<Real> val;
      const typename SPxBasisBase<Real>::Desc& ds = this->desc();

#ifndef ENABLE_ADDITIONAL_CHECKS

      // if the value is available we don't need to recompute it
      if(m_nonbasicValueUpToDate)
         return m_nonbasicValue;

#endif

      if(rep() == COLUMN)
      {
         if(type() == LEAVE)
         {
            for(i = this->nCols() - 1; i >= 0; --i)
            {
               switch(ds.colStatus(i))
               {
               case SPxBasisBase<Real>::Desc::P_ON_UPPER :
                  val += theUCbound[i] * SPxLP::upper(i);
                  //@ val += maxObj(i) * SPxLP::upper(i);
                  break;

               case SPxBasisBase<Real>::Desc::P_ON_LOWER :
                  val += theLCbound[i] * SPxLP::lower(i);
                  //@ val += maxObj(i) * SPxLP::lower(i);
                  break;

               case SPxBasisBase<Real>::Desc::P_FIXED :
                  val += this->maxObj(i) * SPxLP::lower(i);
                  break;

               default:
                  break;
               }
            }

            for(i = this->nRows() - 1; i >= 0; --i)
            {
               switch(ds.rowStatus(i))
               {
               case SPxBasisBase<Real>::Desc::P_ON_UPPER :
                  val += theLRbound[i] * SPxLP::rhs(i);
                  break;

               case SPxBasisBase<Real>::Desc::P_ON_LOWER :
                  val += theURbound[i] * SPxLP::lhs(i);
                  break;

               case SPxBasisBase<Real>::Desc::P_FIXED :
                  val += this->maxRowObj(i) * SPxLP::lhs(i);
                  break;

               default:
                  break;
               }
            }
         }
         else
         {
            assert(type() == ENTER);

            for(i = this->nCols() - 1; i >= 0; --i)
            {
               switch(ds.colStatus(i))
               {
               case SPxBasisBase<Real>::Desc::P_ON_UPPER :
                  val += this->maxObj(i) * theUCbound[i];
                  break;

               case SPxBasisBase<Real>::Desc::P_ON_LOWER :
                  val += this->maxObj(i) * theLCbound[i];
                  break;

               case SPxBasisBase<Real>::Desc::P_FIXED :
                  assert(theLCbound[i] == theUCbound[i]);
                  val += this->maxObj(i) * theLCbound[i];
                  break;

               default:
                  break;
               }
            }

            for(i = this->nRows() - 1; i >= 0; --i)
            {
               switch(ds.rowStatus(i))
               {
               case SPxBasisBase<Real>::Desc::P_ON_UPPER :
                  val += this->maxRowObj(i) * theLRbound[i];
                  break;

               case SPxBasisBase<Real>::Desc::P_ON_LOWER :
                  val += this->maxRowObj(i) * theURbound[i];
                  break;

               case SPxBasisBase<Real>::Desc::P_FIXED :
                  val += this->maxRowObj(i) * theURbound[i];
                  break;

               default:
                  break;
               }
            }
         }
      }
      else
      {
         assert(rep() == ROW);
         assert(type() == ENTER);

         for(i = this->nCols() - 1; i >= 0; --i)
         {
            switch(ds.colStatus(i))
            {
            case SPxBasisBase<Real>::Desc::D_ON_UPPER :
               val += theUCbound[i] * this->lower(i);
               break;

            case SPxBasisBase<Real>::Desc::D_ON_LOWER :
               val += theLCbound[i] * this->upper(i);
               break;

            case SPxBasisBase<Real>::Desc::D_ON_BOTH :
               val += theLCbound[i] * this->upper(i);
               val += theUCbound[i] * this->lower(i);
               break;

            default:
               break;
            }
         }

         for(i = this->nRows() - 1; i >= 0; --i)
         {
            switch(ds.rowStatus(i))
            {
            case SPxBasisBase<Real>::Desc::D_ON_UPPER :
               val += theURbound[i] * this->lhs(i);
               break;

            case SPxBasisBase<Real>::Desc::D_ON_LOWER :
               val += theLRbound[i] * this->rhs(i);
               break;

            case SPxBasisBase<Real>::Desc::D_ON_BOTH :
               val += theLRbound[i] * this->rhs(i);
               val += theURbound[i] * this->lhs(i);
               break;

            default:
               break;
            }
         }
      }

#ifdef ENABLE_ADDITIONAL_CHECKS

      if(m_nonbasicValueUpToDate && NE(m_nonbasicValue, val))
      {
         MSG_ERROR(std::cerr << "stored nonbasic value: " << m_nonbasicValue
                   << ", correct nonbasic value: " << val
                   << ", violation: " << val - m_nonbasicValue << std::endl;)
         assert(EQrel(m_nonbasicValue, val, 1e-12));
      }

#endif

      if(!m_nonbasicValueUpToDate)
      {
         m_nonbasicValue = val;
         m_nonbasicValueUpToDate = true;
      }

      return val;
   }

   template <>
   Real SPxSolverBase<Real>::value()
   {
      assert(isInitialized());

      Real x;

      // calling value() without having a suitable status is an error.
      if(!isInitialized())
         return infinity;

      if(rep() == ROW)
      {
         if(type() == LEAVE)
            x = SPxLP::spxSense() * (coPvec() * fRhs()); // the contribution of maxRowObj() is missing
         else
            x = SPxLP::spxSense() * (nonbasicValue() + (coPvec() * fRhs()));
      }
      else
         x = SPxLP::spxSense() * (nonbasicValue() + fVec() * coPrhs());

      return x + this->objOffset();
   }

   template <>
   bool SPxSolverBase<Real>::updateNonbasicValue(Real objChange)
   {
      if(m_nonbasicValueUpToDate)
         m_nonbasicValue += objChange;

      MSG_DEBUG(std::cout
                << "Iteration: " << iteration()
                << ": updated objValue: " << objChange
                << ", new value: " << m_nonbasicValue
                << ", correct value: " << nonbasicValue()
                << std::endl;
               )

      return m_nonbasicValueUpToDate;
   }



   template <>
   void SPxSolverBase<Real>::setFeastol(Real d)
   {

      if(d <= 0.0)
         throw SPxInterfaceException("XSOLVE30 Cannot set feastol less than or equal to zero.");

      if(theRep == COLUMN)
         m_entertol = d;
      else
         m_leavetol = d;
   }

   template <>
   void SPxSolverBase<Real>::setOpttol(Real d)
   {

      if(d <= 0.0)
         throw SPxInterfaceException("XSOLVE31 Cannot set opttol less than or equal to zero.");

      if(theRep == COLUMN)
         m_leavetol = d;
      else
         m_entertol = d;
   }

   template <>
   void SPxSolverBase<Real>::setDelta(Real d)
   {

      if(d <= 0.0)
         throw SPxInterfaceException("XSOLVE32 Cannot set delta less than or equal to zero.");

      m_entertol = d;
      m_leavetol = d;
   }

   template <>
   void SPxSolverBase<Real>::hyperPricing(bool h)
   {
      hyperPricingEnter = h;
      hyperPricingLeave = h;

      if(h)
      {
         updateViols.setMax(dim());
         updateViolsCo.setMax(coDim());
      }
   }

   template <>
   SPxSolverBase<Real>::SPxSolverBase(
      Type            p_type,
      Representation  p_rep,
      Timer::TYPE     ttype)
      : theType(p_type)
      , thePricing(FULL)
      , theRep(p_rep)
      , polishObj(POLISH_OFF)
      , theTime(0)
      , timerType(ttype)
      , theCumulativeTime(0.0)
      , maxIters(-1)
      , maxTime(infinity)
      , nClckSkipsLeft(0)
      , nCallsToTimelim(0)
      , objLimit(infinity)
      , m_status(UNKNOWN)
      , m_nonbasicValue(0.0)
      , m_nonbasicValueUpToDate(false)
      , m_pricingViol(0.0)
      , m_pricingViolUpToDate(false)
      , m_pricingViolCo(0.0)
      , m_pricingViolCoUpToDate(false)
      , m_numViol(0)
      , theShift(0)
      , m_maxCycle(100)
      , m_numCycle(0)
      , initialized(false)
      , solveVector2(0)
      , solveVector3(0)
      , coSolveVector2(0)
      , coSolveVector3(0)
      , freePricer(false)
      , freeRatioTester(false)
      , freeStarter(false)
      , displayLine(0)
      , displayFreq(200)
      , sparsePricingFactor(SPARSITYFACTOR)
      , getStartingDecompBasis(false)
      , computeDegeneracy(false)
      , degenCompIterOffset(0)
      , fullPerturbation(false)
      , printBasisMetric(0)
      , unitVecs(0)
      , primVec(0, Param::epsilon())
      , dualVec(0, Param::epsilon())
      , addVec(0, Param::epsilon())
      , thepricer(0)
      , theratiotester(0)
      , thestarter(0)
      , boundrange(0.0)
      , siderange(0.0)
      , objrange(0.0)
      , infeasibilities(0)
      , infeasibilitiesCo(0)
      , isInfeasible(0)
      , isInfeasibleCo(0)
      , sparsePricingLeave(false)
      , sparsePricingEnter(false)
      , sparsePricingEnterCo(false)
      , hyperPricingLeave(true)
      , hyperPricingEnter(true)
      , remainingRoundsLeave(0)
      , remainingRoundsEnter(0)
      , remainingRoundsEnterCo(0)
      , weights(0)
      , coWeights(0)
      , weightsAreSetup(false)
      , integerVariables(0)
   {
      theTime = TimerFactory::createTimer(timerType);

      setDelta(DEFAULT_BND_VIOL);

      this->theLP = this;
      initRep(p_rep);

      // info: SPxBasisBase is not consistent in this moment.
      //assert(SPxSolverBase<Real>::isConsistent());
   }

   template <>
   SPxSolverBase<Real>::~SPxSolverBase()
   {
      assert(!freePricer || thepricer != 0);
      assert(!freeRatioTester || theratiotester != 0);
      assert(!freeStarter || thestarter != 0);

      if(freePricer)
      {
         delete thepricer;
         thepricer = 0;
      }

      if(freeRatioTester)
      {
         delete theratiotester;
         theratiotester = 0;
      }

      if(freeStarter)
      {
         delete thestarter;
         thestarter = 0;
      }

      // free timer
      assert(theTime);
      theTime->~Timer();
      spx_free(theTime);
   }


   template <>
   SPxSolverBase<Real>& SPxSolverBase<Real>::operator=(const SPxSolverBase<Real>& base)
   {
      if(this != &base)
      {
         SPxLP::operator=(base);
         SPxBasisBase<Real>::operator=(base);
         theType = base.theType;
         thePricing = base.thePricing;
         theRep = base.theRep;
         polishObj = base.polishObj;
         timerType = base.timerType;
         maxIters = base.maxIters;
         maxTime = base.maxTime;
         objLimit = base.objLimit;
         m_status = base.m_status;
         m_nonbasicValue = base.m_nonbasicValue;
         m_nonbasicValueUpToDate = base.m_nonbasicValueUpToDate;
         m_pricingViol = base.m_pricingViol;
         m_pricingViolUpToDate = base.m_pricingViolUpToDate;
         m_pricingViolCo = base.m_pricingViolCo;
         m_pricingViolCoUpToDate = base.m_pricingViolCoUpToDate;
         m_numViol = base.m_numViol;
         m_entertol = base.m_entertol;
         m_leavetol = base.m_leavetol;
         theShift = base.theShift;
         lastShift = base.lastShift;
         m_maxCycle = base.m_maxCycle;
         m_numCycle = base.m_numCycle;
         initialized = base.initialized;
         instableLeaveNum = base.instableLeaveNum;
         instableLeave = base.instableLeave;
         instableLeaveVal = base.instableLeaveVal;
         instableEnterId = base.instableEnterId;
         instableEnter = base.instableEnter;
         instableEnterVal = base.instableEnterVal;
         displayLine = base.displayLine;
         displayFreq = base.displayFreq;
         sparsePricingFactor = base.sparsePricingFactor;
         getStartingDecompBasis = base.getStartingDecompBasis;
         computeDegeneracy = base.computeDegeneracy;
         degenCompIterOffset = base.degenCompIterOffset;
         decompIterationLimit = base.decompIterationLimit;
         fullPerturbation = base.fullPerturbation;
         printBasisMetric = base.printBasisMetric;
         unitVecs = base.unitVecs;
         primRhs = base.primRhs;
         primVec = base.primVec;
         dualRhs = base.dualRhs;
         dualVec = base.dualVec;
         addVec = base.addVec;
         theURbound = base.theURbound;
         theLRbound = base.theLRbound;
         theUCbound = base.theUCbound;
         theLCbound = base.theLCbound;
         theUBbound = base.theUBbound;
         theLBbound = base.theLBbound;
         theCoTest = base.theCoTest;
         theTest = base.theTest;
         primalRay = base.primalRay;
         dualFarkas = base.dualFarkas;
         leaveCount = base.leaveCount;
         enterCount = base.enterCount;
         theCumulativeTime = base.theCumulativeTime;
         primalCount = base.primalCount;
         polishCount = base.polishCount;
         boundflips = base.boundflips;
         totalboundflips = base.totalboundflips;
         enterCycles = base.enterCycles;
         leaveCycles = base.leaveCycles;
         enterDegenCand = base.enterDegenCand;
         leaveDegenCand = base.leaveDegenCand;
         primalDegenSum = base.primalDegenSum;
         boundrange = base.boundrange;
         siderange = base.siderange;
         objrange = base.objrange;
         infeasibilities = base.infeasibilities;
         infeasibilitiesCo = base.infeasibilitiesCo;
         isInfeasible = base.isInfeasible;
         isInfeasibleCo = base.isInfeasibleCo;
         sparsePricingLeave = base.sparsePricingLeave;
         sparsePricingEnter = base.sparsePricingEnter;
         sparsePricingEnterCo = base.sparsePricingEnterCo;
         sparsePricingFactor = base.sparsePricingFactor;
         hyperPricingLeave = base.hyperPricingLeave;
         hyperPricingEnter = base.hyperPricingEnter;
         remainingRoundsLeave = base.remainingRoundsLeave;
         remainingRoundsEnter = base.remainingRoundsEnter;
         remainingRoundsEnterCo = base.remainingRoundsEnterCo;
         weights = base.weights;
         coWeights = base.coWeights;
         weightsAreSetup = base.weightsAreSetup;
         spxout = base.spxout;
         integerVariables = base.integerVariables;

         if(base.theRep == COLUMN)
         {
            thevectors   = this->colSet();
            thecovectors = this->rowSet();
            theFrhs      = &primRhs;
            theFvec      = &primVec;
            theCoPrhs    = &dualRhs;
            theCoPvec    = &dualVec;
            thePvec      = &addVec;
            theRPvec     = theCoPvec;
            theCPvec     = thePvec;
            theUbound    = &theUCbound;
            theLbound    = &theLCbound;
            theCoUbound  = &theURbound;
            theCoLbound  = &theLRbound;
         }
         else
         {
            assert(base.theRep == ROW);

            thevectors   = this->rowSet();
            thecovectors = this->colSet();
            theFrhs      = &dualRhs;
            theFvec      = &dualVec;
            theCoPrhs    = &primRhs;
            theCoPvec    = &primVec;
            thePvec      = &addVec;
            theRPvec     = thePvec;
            theCPvec     = theCoPvec;
            theUbound    = &theURbound;
            theLbound    = &theLRbound;
            theCoUbound  = &theUCbound;
            theCoLbound  = &theLCbound;
         }

         SPxBasisBase<Real>::theLP = this;

         assert(!freePricer || thepricer != 0);
         assert(!freeRatioTester || theratiotester != 0);
         assert(!freeStarter || thestarter != 0);

         // thepricer
         if(freePricer)
         {
            delete thepricer;
            thepricer = 0;
         }

         if(base.thepricer == 0)
         {
            thepricer = 0;
            freePricer = false;
         }
         else
         {
            thepricer = base.thepricer->clone();
            freePricer = true;
            thepricer->load(this);
         }

         // theratiotester
         if(freeRatioTester)
         {
            delete theratiotester;
            theratiotester = 0;
         }

         if(base.theratiotester == 0)
         {
            theratiotester = 0;
            freeRatioTester = false;
         }
         else
         {
            theratiotester = base.theratiotester->clone();
            freeRatioTester = true;
            theratiotester->load(this);
         }

         // thestarter
         if(freeStarter)
         {
            delete thestarter;
            thestarter = 0;
         }

         if(base.thestarter == 0)
         {
            thestarter = 0;
            freeStarter = false;
         }
         else
         {
            thestarter = base.thestarter->clone();
            freeStarter = true;
         }

         assert(SPxSolverBase<Real>::isConsistent());
      }

      return *this;
   }


   template <class R>
   SPxSolverBase<R>::SPxSolverBase(const SPxSolverBase<R>& base)
      : SPxLP(base)
      , SPxBasisBase<R>(this->basSe)
      , theType(base.theType)
      , thePricing(base.thePricing)
      , theRep(base.theRep)
      , polishObj(base.polishObj)
      , timerType(base.timerType)
      , theCumulativeTime(base.theCumulativeTime)
      , maxIters(base.maxIters)
      , maxTime(base.maxTime)
      , nClckSkipsLeft(base.nClckSkipsLeft)
      , nCallsToTimelim(base.nCallsToTimelim)
      , objLimit(base.objLimit)
      , m_status(base.m_status)
      , m_nonbasicValue(base.m_nonbasicValue)
      , m_nonbasicValueUpToDate(base.m_nonbasicValueUpToDate)
      , m_pricingViol(base.m_pricingViol)
      , m_pricingViolUpToDate(base.m_pricingViolUpToDate)
      , m_pricingViolCo(base.m_pricingViolCo)
      , m_pricingViolCoUpToDate(base.m_pricingViolCoUpToDate)
      , m_numViol(base.m_numViol)
      , m_entertol(base.m_entertol)
      , m_leavetol(base.m_leavetol)
      , theShift(base.theShift)
      , lastShift(base.lastShift)
      , m_maxCycle(base.m_maxCycle)
      , m_numCycle(base.m_numCycle)
      , initialized(base.initialized)
      , solveVector2(0)
      , solveVector2rhs(base.solveVector2rhs)
      , solveVector3(0)
      , solveVector3rhs(base.solveVector3rhs)
      , coSolveVector2(0)
      , coSolveVector2rhs(base.coSolveVector2rhs)
      , coSolveVector3(0)
      , coSolveVector3rhs(base.coSolveVector3rhs)
      , instableLeaveNum(base.instableLeaveNum)
      , instableLeave(base.instableLeave)
      , instableLeaveVal(base.instableLeaveVal)
      , instableEnterId(base.instableEnterId)
      , instableEnter(base.instableEnter)
      , instableEnterVal(base.instableEnterVal)
      , displayLine(base.displayLine)
      , displayFreq(base.displayFreq)
      , sparsePricingFactor(base.sparsePricingFactor)
      , getStartingDecompBasis(base.getStartingDecompBasis)
      , computeDegeneracy(base.computeDegeneracy)
      , degenCompIterOffset(base.degenCompIterOffset)
      , decompIterationLimit(base.decompIterationLimit)
      , fullPerturbation(base.fullPerturbation)
      , printBasisMetric(base.printBasisMetric)
      , unitVecs(base.unitVecs)
      , primRhs(base.primRhs)
      , primVec(base.primVec)
      , dualRhs(base.dualRhs)
      , dualVec(base.dualVec)
      , addVec(base.addVec)
      , theURbound(base.theURbound)
      , theLRbound(base.theLRbound)
      , theUCbound(base.theUCbound)
      , theLCbound(base.theLCbound)
      , theUBbound(base.theUBbound)
      , theLBbound(base.theLBbound)
      , theCoTest(base.theCoTest)
      , theTest(base.theTest)
      , primalRay(base.primalRay)
      , dualFarkas(base.dualFarkas)
      , leaveCount(base.leaveCount)
      , enterCount(base.enterCount)
      , primalCount(base.primalCount)
      , polishCount(base.polishCount)
      , boundflips(base.boundflips)
      , totalboundflips(base.totalboundflips)
      , enterCycles(base.enterCycles)
      , leaveCycles(base.leaveCycles)
      , enterDegenCand(base.enterDegenCand)
      , leaveDegenCand(base.leaveDegenCand)
      , primalDegenSum(base.primalDegenSum)
      , dualDegenSum(base.dualDegenSum)
      , boundrange(base.boundrange)
      , siderange(base.siderange)
      , objrange(base.objrange)
      , infeasibilities(base.infeasibilities)
      , infeasibilitiesCo(base.infeasibilitiesCo)
      , isInfeasible(base.isInfeasible)
      , isInfeasibleCo(base.isInfeasibleCo)
      , sparsePricingLeave(base.sparsePricingLeave)
      , sparsePricingEnter(base.sparsePricingEnter)
      , sparsePricingEnterCo(base.sparsePricingEnterCo)
      , hyperPricingLeave(base.hyperPricingLeave)
      , hyperPricingEnter(base.hyperPricingEnter)
      , remainingRoundsLeave(base.remainingRoundsLeave)
      , remainingRoundsEnter(base.remainingRoundsEnter)
      , remainingRoundsEnterCo(base.remainingRoundsEnterCo)
      , weights(base.weights)
      , coWeights(base.coWeights)
      , weightsAreSetup(base.weightsAreSetup)
      , spxout(base.spxout)
      , integerVariables(base.integerVariables)
   {
      theTime = TimerFactory::createTimer(timerType);

      if(base.theRep == COLUMN)
      {
         thevectors   = this->colSet();
         thecovectors = this->rowSet();
         theFrhs      = &primRhs;
         theFvec      = &primVec;
         theCoPrhs    = &dualRhs;
         theCoPvec    = &dualVec;
         thePvec      = &addVec;
         theRPvec     = theCoPvec;
         theCPvec     = thePvec;
         theUbound    = &theUCbound;
         theLbound    = &theLCbound;
         theCoUbound  = &theURbound;
         theCoLbound  = &theLRbound;
      }
      else
      {
         assert(base.theRep == ROW);

         thevectors   = this->rowSet();
         thecovectors = this->colSet();
         theFrhs      = &dualRhs;
         theFvec      = &dualVec;
         theCoPrhs    = &primRhs;
         theCoPvec    = &primVec;
         thePvec      = &addVec;
         theRPvec     = thePvec;
         theCPvec     = theCoPvec;
         theUbound    = &theURbound;
         theLbound    = &theLRbound;
         theCoUbound  = &theUCbound;
         theCoLbound  = &theLCbound;
      }

      SPxBasisBase<R>::theLP = this;

      if(base.thepricer == 0)
      {
         thepricer = 0;
         freePricer = false;
      }
      else
      {
         thepricer = base.thepricer->clone();
         freePricer = true;
         thepricer->clear();
         thepricer->load(this);
      }

      if(base.theratiotester == 0)
      {
         theratiotester = 0;
         freeRatioTester = false;
      }
      else
      {
         theratiotester = base.theratiotester->clone();
         freeRatioTester = true;
         theratiotester->clear();
         theratiotester->load(this);
      }

      if(base.thestarter == 0)
      {
         thestarter = 0;
         freeStarter = false;
      }
      else
      {
         thestarter = base.thestarter->clone();
         freeStarter = true;
      }

      assert(SPxSolverBase<R>::isConsistent());
   }

   template <>
   bool SPxSolverBase<Real>::isConsistent() const
   {
#ifdef ENABLE_CONSISTENCY_CHECKS

      if(epsilon() < 0)
         return MSGinconsistent("SPxSolverBase");

      if(primVec.delta().getEpsilon() != dualVec.delta().getEpsilon())
         return MSGinconsistent("SPxSolverBase");

      if(dualVec.delta().getEpsilon() != addVec.delta().getEpsilon())
         return MSGinconsistent("SPxSolverBase");

      if(unitVecs.size() < SPxLP::this->nCols() || unitVecs.size() < SPxLP::this->nRows())
         return MSGinconsistent("SPxSolverBase");

      if(initialized)
      {
         if(theFrhs->dim() != dim())
            return MSGinconsistent("SPxSolverBase");

         if(theFvec->dim() != dim())
            return MSGinconsistent("SPxSolverBase");

         if(theCoPrhs->dim() != dim())
            return MSGinconsistent("SPxSolverBase");

         if(thePvec->dim() != coDim())
            return MSGinconsistent("SPxSolverBase");

         if(theCoPvec->dim() != dim())
            return MSGinconsistent("SPxSolverBase");

         if(theTest.dim() != coDim())
            return MSGinconsistent("SPxSolverBase");

         if(theCoTest.dim() != dim())
            return MSGinconsistent("SPxSolverBase");

         if(theURbound.dim() != SPxLP::this->nRows())
            return MSGinconsistent("SPxSolverBase");

         if(theLRbound.dim() != SPxLP::this->nRows())
            return MSGinconsistent("SPxSolverBase");

         if(theUCbound.dim() != SPxLP::this->nCols())
            return MSGinconsistent("SPxSolverBase");

         if(theLCbound.dim() != SPxLP::this->nCols())
            return MSGinconsistent("SPxSolverBase");

         if(theUBbound.dim() != dim())
            return MSGinconsistent("SPxSolverBase");

         if(theLBbound.dim() != dim())
            return MSGinconsistent("SPxSolverBase");
      }

      if(rep() == COLUMN)
      {
         if(thecovectors !=
               reinterpret_cast<const SVSet*>(static_cast<const LPRowSet*>(this))
               || thevectors !=
               reinterpret_cast<const SVSet*>(static_cast<const LPColSet*>(this))
               || theFrhs != &primRhs ||
               theFvec != &primVec ||
               theCoPrhs != &dualRhs ||
               theCoPvec != &dualVec ||
               thePvec != &addVec ||
               theRPvec != theCoPvec ||
               theCPvec != thePvec ||
               theUbound != &theUCbound ||
               theLbound != &theLCbound ||
               theCoUbound != &theURbound ||
               theCoLbound != &theLRbound)
            return MSGinconsistent("SPxSolverBase");
      }
      else
      {
         if(thecovectors
               != reinterpret_cast<const SVSet*>(static_cast<const LPColSet*>(this))
               || thevectors
               != reinterpret_cast<const SVSet*>(static_cast<const LPRowSet*>(this))
               || theFrhs != &dualRhs ||
               theFvec != &dualVec ||
               theCoPrhs != &primRhs ||
               theCoPvec != &primVec ||
               thePvec != &addVec ||
               theRPvec != thePvec ||
               theCPvec != theCoPvec ||
               theUbound != &theURbound ||
               theLbound != &theLRbound ||
               theCoUbound != &theUCbound ||
               theCoLbound != &theLCbound)
            return MSGinconsistent("SPxSolverBase");
      }

      return SPxLP::isConsistent()
             && primRhs.isConsistent()
             && primVec.isConsistent()
             && dualRhs.isConsistent()
             && dualVec.isConsistent()
             && addVec.isConsistent()
             && theTest.isConsistent()
             && theCoTest.isConsistent()
             && theURbound.isConsistent()
             && theLRbound.isConsistent()
             && theUCbound.isConsistent()
             && theLCbound.isConsistent()
             && SPxBasisBase<Real>::isConsistent()
             ;
#else
      return true;
#endif
   }


   template <>
   void SPxSolverBase<Real>::setTerminationTime(Real p_time)
   {
      if(p_time < 0.0)
         p_time = 0.0;

      maxTime = p_time;
   }

   template <>
   Real SPxSolverBase<Real>::terminationTime() const
   {
      return maxTime;
   }

   template <>
   void SPxSolverBase<Real>::setTerminationIter(int p_iteration)
   {
      if(p_iteration < 0)
         p_iteration = -1;

      maxIters = p_iteration;
   }

   template <>
   int SPxSolverBase<Real>::terminationIter() const
   {
      return maxIters;
   }

   // returns whether current time limit is reached; call to time() may be skipped unless \p forceCheck is true
   template <>
   bool SPxSolverBase<Real>::isTimeLimitReached(const bool forceCheck)
   {
      // always update the number of calls, since the user might set a time limit later in the solving process
      ++nCallsToTimelim;

      // check if a time limit is actually set
      if(maxTime >= infinity)
         return false;

      // check if the expensive system call to update the time should be skipped again
      if(forceCheck || nCallsToTimelim < NINITCALLS ||  nClckSkipsLeft <= 0)
      {
         Real currtime = time();

         if(currtime >= maxTime)
            return true;

         // determine the number of times the clock can be skipped again.
         int nClckSkips = MAXNCLCKSKIPS;
         Real avgtimeinterval = (currtime + cumulativeTime()) / (Real)(nCallsToTimelim);

         // it would not be safe to skip the clock so many times since we are approaching the time limit
         if(SAFETYFACTOR * (maxTime - currtime) / (avgtimeinterval + 1e-6) < nClckSkips)
            nClckSkips = 0;

         nClckSkipsLeft = nClckSkips;
      }
      else
         --nClckSkipsLeft;

      return false;
   }


   /**@todo A first version for the termination value is
    *       implemented. Currently we check if no bound violations (shifting)
    *       is present. It might be even possible to use this termination
    *       value in case of bound violations (shifting) but in this case it
    *       is quite difficult to determine if we already reached the limit.
    */
   template <>
   void SPxSolverBase<Real>::setTerminationValue(Real p_value)
   {
      objLimit = p_value;
   }

   template <>
   Real SPxSolverBase<Real>::terminationValue() const
   {
      return objLimit;
   }

   template <>
   typename SPxSolverBase<Real>::VarStatus
   SPxSolverBase<Real>::basisStatusToVarStatus(typename SPxBasisBase<Real>::Desc::Status stat) const
   {
      VarStatus vstat;

      switch(stat)
      {
      case SPxBasisBase<Real>::Desc::P_ON_LOWER:
         vstat = ON_LOWER;
         break;

      case SPxBasisBase<Real>::Desc::P_ON_UPPER:
         vstat = ON_UPPER;
         break;

      case SPxBasisBase<Real>::Desc::P_FIXED:
         vstat = FIXED;
         break;

      case SPxBasisBase<Real>::Desc::P_FREE:
         vstat = ZERO;
         break;

      case SPxBasisBase<Real>::Desc::D_ON_UPPER:
      case SPxBasisBase<Real>::Desc::D_ON_LOWER:
      case SPxBasisBase<Real>::Desc::D_ON_BOTH:
      case SPxBasisBase<Real>::Desc::D_UNDEFINED:
      case SPxBasisBase<Real>::Desc::D_FREE:
         vstat = BASIC;
         break;

      default:
         MSG_ERROR(std::cerr << "ESOLVE26 ERROR: unknown basis status (" << static_cast<int>(stat) << ")"
                   << std::endl;)
         throw SPxInternalCodeException("XSOLVE22 This should never happen.");
      }

      return vstat;
   }

   template <>
   typename SPxBasisBase<Real>::Desc::Status
   SPxSolverBase<Real>::varStatusToBasisStatusRow(int row, SPxSolverBase<Real>::VarStatus stat) const
   {
      typename SPxBasisBase<Real>::Desc::Status rstat;

      switch(stat)
      {
      case FIXED :
         assert(EQ(this->rhs(row), this->lhs(row), feastol()));
         rstat = SPxBasisBase<Real>::Desc::P_FIXED;
         break;

      case ON_UPPER :
         assert(this->rhs(row) < infinity);
         rstat = this->lhs(row) < this->rhs(row)
                 ? SPxBasisBase<Real>::Desc::P_ON_UPPER
                 : SPxBasisBase<Real>::Desc::P_FIXED;
         break;

      case ON_LOWER :
         assert(this->lhs(row) > -infinity);
         rstat = this->lhs(row) < this->rhs(row)
                 ? SPxBasisBase<Real>::Desc::P_ON_LOWER
                 : SPxBasisBase<Real>::Desc::P_FIXED;
         break;

      case ZERO :

         /* A 'free' row (i.e., infinite lower & upper bounds) does not really make sense. The user
          * might (think to) know better, e.g., when temporarily turning off a row. We therefore apply
          * the same adjustment as in the column case in varStatusToBasisStatusCol(). */
         if(this->lhs(row) <= -infinity && this->rhs(row) >= infinity)
            rstat = SPxBasisBase<Real>::Desc::P_FREE;
         else
         {
            if(this->lhs(row) == this->rhs(row))
            {
               assert(this->rhs(row) < infinity);
               rstat = SPxBasisBase<Real>::Desc::P_FIXED;
            }
            else
            {
               if(this->lhs(row) > -infinity)
                  rstat = SPxBasisBase<Real>::Desc::P_ON_LOWER;
               else
               {
                  assert(this->rhs(row) < infinity);
                  rstat = SPxBasisBase<Real>::Desc::P_ON_UPPER;
               }
            }
         }

         break;

      case BASIC :
         rstat = this->dualRowStatus(row);
         break;

      default:
         MSG_ERROR(std::cerr << "ESOLVE27 ERROR: unknown VarStatus (" << int(stat) << ")"
                   << std::endl;)
         throw SPxInternalCodeException("XSOLVE23 This should never happen.");
      }

      return rstat;
   }

   template <>
   typename SPxBasisBase<Real>::Desc::Status
   SPxSolverBase<Real>::varStatusToBasisStatusCol(int col, SPxSolverBase<Real>::VarStatus stat) const
   {
      typename SPxBasisBase<Real>::Desc::Status cstat;

      switch(stat)
      {
      case FIXED :
         if(this->upper(col) == this->lower(col))
            cstat = SPxBasisBase<Real>::Desc::P_FIXED;
         else if(this->maxObj(col) > 0.0)
            cstat = SPxBasisBase<Real>::Desc::P_ON_UPPER;
         else
            cstat = SPxBasisBase<Real>::Desc::P_ON_LOWER;

         break;

      case ON_UPPER :
         assert(this->upper(col) < infinity);
         cstat = this->lower(col) < this->upper(col)
                 ? SPxBasisBase<Real>::Desc::P_ON_UPPER
                 : SPxBasisBase<Real>::Desc::P_FIXED;
         break;

      case ON_LOWER :
         assert(this->lower(col) > -infinity);
         cstat = this->lower(col) < this->upper(col)
                 ? SPxBasisBase<Real>::Desc::P_ON_LOWER
                 : SPxBasisBase<Real>::Desc::P_FIXED;
         break;

      case ZERO :

         /* In this case the upper and lower bounds on the variable should be infinite. The bounds
          * might, however, have changed and we try to recover from this by changing the status to
          * 'resonable' settings. A solve has to find the correct values afterwards. Note that the
          * approach below is consistent with changesoplex.cpp (e.g., changeUpperStatus() and
          * changeLowerStatus() ). */
         if(this->lower(col) <= -infinity && this->upper(col) >= infinity)
            cstat = SPxBasisBase<Real>::Desc::P_FREE;
         else
         {
            if(this->lower(col) == this->upper(col))
            {
               assert(this->upper(col) < infinity);
               cstat = SPxBasisBase<Real>::Desc::P_FIXED;
            }
            else
            {
               if(this->lower(col) > -infinity)
                  cstat = SPxBasisBase<Real>::Desc::P_ON_LOWER;
               else
               {
                  assert(this->upper(col) < infinity);
                  cstat = SPxBasisBase<Real>::Desc::P_ON_UPPER;
               }
            }
         }

         break;

      case BASIC :
         cstat = this->dualColStatus(col);
         break;

      default:
         MSG_ERROR(std::cerr << "ESOLVE28 ERROR: unknown VarStatus (" << int(stat) << ")"
                   << std::endl;)
         throw SPxInternalCodeException("XSOLVE24 This should never happen.");
      }

      return cstat;
   }

   template <>
   typename SPxSolverBase<Real>::VarStatus SPxSolverBase<Real>::getBasisRowStatus(int row) const
   {
      assert(0 <= row && row < this->nRows());
      return basisStatusToVarStatus(this->desc().rowStatus(row));
   }

   template <>
   typename SPxSolverBase<Real>::VarStatus SPxSolverBase<Real>::getBasisColStatus(int col) const
   {
      assert(0 <= col && col < this->nCols());
      return basisStatusToVarStatus(this->desc().colStatus(col));
   }

   template <>
   typename SPxSolverBase<Real>::Status SPxSolverBase<Real>::getBasis(VarStatus row[], VarStatus col[],
         const int rowsSize, const int colsSize) const
   {
      const typename SPxBasisBase<Real>::Desc& d = this->desc();
      int i;

      assert(rowsSize < 0 || rowsSize >= this->nRows());
      assert(colsSize < 0 || colsSize >= this->nCols());

      if(col)
         for(i = this->nCols() - 1; i >= 0; --i)
            col[i] = basisStatusToVarStatus(d.colStatus(i));

      if(row)
         for(i = this->nRows() - 1; i >= 0; --i)
            row[i] = basisStatusToVarStatus(d.rowStatus(i));

      return status();
   }

   template <>
   bool SPxSolverBase<Real>::isBasisValid(DataArray<VarStatus> p_rows, DataArray<VarStatus> p_cols)
   {

      int basisdim;

      if(p_rows.size() != this->nRows() || p_cols.size() != this->nCols())
         return false;

      basisdim = 0;

      for(int row = this->nRows() - 1; row >= 0; --row)
      {
         if(p_rows[row] == UNDEFINED)
            return false;
         // row is basic
         else if(p_rows[row] == BASIC)
         {
            basisdim++;
         }
         // row is nonbasic
         else
         {
            if((p_rows[row] == FIXED && this->lhs(row) != this->rhs(row))
                  || (p_rows[row] == ON_UPPER && this->rhs(row) >= infinity)
                  || (p_rows[row] == ON_LOWER && this->lhs(row) <= -infinity))
               return false;
         }
      }

      for(int col = this->nCols() - 1; col >= 0; --col)
      {
         if(p_cols[col] == UNDEFINED)
            return false;
         // col is basic
         else if(p_cols[col] == BASIC)
         {
            basisdim++;
         }
         // col is nonbasic
         else
         {
            if((p_cols[col] == FIXED && this->lower(col) != this->upper(col))
                  || (p_cols[col] == ON_UPPER && this->upper(col) >= infinity)
                  || (p_cols[col] == ON_LOWER && this->lower(col) <= -infinity))
               return false;
         }
      }

      if(basisdim != dim())
         return false;

      // basis valid
      return true;
   }

   template <>
   void SPxSolverBase<Real>::setBasis(const VarStatus p_rows[], const VarStatus p_cols[])
   {
      if(SPxBasisBase<Real>::status() == SPxBasisBase<Real>::NO_PROBLEM)
         SPxBasisBase<Real>::load(this, false);

      typename SPxBasisBase<Real>::Desc ds = this->desc();
      int i;

      for(i = 0; i < this->nRows(); i++)
         ds.rowStatus(i) = varStatusToBasisStatusRow(i, p_rows[i]);

      for(i = 0; i < this->nCols(); i++)
         ds.colStatus(i) = varStatusToBasisStatusCol(i, p_cols[i]);

      loadBasis(ds);
      forceRecompNonbasicValue();
   }

   // NOTE: This only works for the row representation. Need to update to account for column representation.
   // The degenvec differs relative to the algorithm being used.
   // For the primal simplex, degenvec is the primal solution values.
   // For the dual simplex, the degenvec is the feasvec (ROW) and pVec (COLUMN).
   template <>
   Real SPxSolverBase<Real>::getDegeneracyLevel(Vector degenvec)
   {
      int numDegenerate = 0;
      Real degeneracyLevel = 0;

      // iterating over all columns in the basis matrix
      // this identifies the basis indices and those that have a zero dual multiplier (rows) or zero reduced cost (cols).
      if(rep() == ROW)
      {
         for(int i = 0; i < this->nCols();
               ++i)   // @todo Check the use of numColsReal for the reduced problem.
         {
            // degeneracy in the dual simplex exists if there are rows with a zero dual multiplier or columns with a zero
            // reduced costs. This requirement is regardless of the objective sense.
            if(isZero(degenvec[i], feastol()))
               numDegenerate++;
         }

         if(type() == ENTER)     // dual simplex
            degeneracyLevel = Real(numDegenerate) / this->nCols();
         else                    // primal simplex
         {
            assert(type() == LEAVE);
            Real degenVars = (numDegenerate > (this->nCols() - this->nRows())) ? Real(numDegenerate -
                             (this->nCols() - this->nRows())) : 0.0;
            degeneracyLevel = degenVars / this->nRows();
         }
      }
      else
      {
         assert(rep() == COLUMN);

         for(int i = 0; i < this->nCols(); i++)
         {
            if(type() == LEAVE)     // dual simplex
            {
               if(isZero(this->maxObj()[i] - degenvec[i], feastol()))
                  numDegenerate++;
            }
            else                    // primal simplex
            {
               assert(type() == ENTER);

               if(isZero(degenvec[i], feastol()))
                  numDegenerate++;
            }
         }


         if(type() == LEAVE)     // dual simplex
         {
            Real degenVars = this->nRows() > numDegenerate ? Real(this->nRows() - numDegenerate) : 0.0;
            degeneracyLevel = degenVars / this->nCols();
         }
         else                    // primal simplex
         {
            assert(type() == ENTER);
            Real degenVars = (numDegenerate > (this->nCols() - this->nRows())) ? Real(numDegenerate -
                             (this->nCols() - this->nRows())) : 0.0;
            degeneracyLevel = degenVars / this->nRows();
         }
      }

      return degeneracyLevel;
   }

   template <>
   void SPxSolverBase<Real>::getNdualNorms(int& nnormsRow, int& nnormsCol) const
   {
      nnormsRow = 0;
      nnormsCol = 0;

      if(weightsAreSetup)
      {
         if(type() == SPxSolverBase<Real>::LEAVE && rep() == SPxSolverBase<Real>::COLUMN)
         {
            nnormsRow = coWeights.dim();
            nnormsCol = 0;

            assert(nnormsRow == dim());
         }
         else if(type() == SPxSolverBase<Real>::ENTER && rep() == SPxSolverBase<Real>::ROW)
         {
            nnormsRow = weights.dim();
            nnormsCol = coWeights.dim();

            assert(nnormsRow == coDim());
            assert(nnormsCol == dim());
         }
      }
   }

   template <>
   bool SPxSolverBase<Real>::getDualNorms(int& nnormsRow, int& nnormsCol, Real * norms) const
   {
      nnormsRow = 0;
      nnormsCol = 0;

      if(!weightsAreSetup)
         return false;

      if(type() == SPxSolverBase<Real>::LEAVE && rep() == SPxSolverBase<Real>::COLUMN)
      {
         nnormsCol = 0;
         nnormsRow = coWeights.dim();

         assert(nnormsRow == dim());

         for(int i = 0; i < nnormsRow; ++i)
            norms[i] = coWeights[i];
      }
      else if(type() == SPxSolverBase<Real>::ENTER && rep() == SPxSolverBase<Real>::ROW)
      {
         nnormsRow = weights.dim();
         nnormsCol = coWeights.dim();

         assert(nnormsCol == dim());
         assert(nnormsRow == coDim());

         for(int i = 0; i < nnormsRow; ++i)
            norms[i] = weights[i];

         for(int i = 0; i < nnormsCol; ++i)
            norms[nnormsRow + i] = coWeights[i];
      }
      else
         return false;

      return true;
   }

   template <>
   bool SPxSolverBase<Real>::setDualNorms(int nnormsRow, int nnormsCol, Real * norms)
   {
      weightsAreSetup = false;

      if(type() == SPxSolverBase<Real>::LEAVE && rep() == SPxSolverBase<Real>::COLUMN)
      {
         coWeights.reDim(dim(), false);
         assert(coWeights.dim() >= nnormsRow);

         for(int i = 0; i < nnormsRow; ++i)
            coWeights[i] = norms[i];

         weightsAreSetup = true;
      }
      else if(type() == SPxSolverBase<Real>::ENTER && rep() == SPxSolverBase<Real>::ROW)
      {
         weights.reDim(coDim(), false);
         coWeights.reDim(dim(), false);
         assert(weights.dim() >= nnormsRow);
         assert(coWeights.dim() >= nnormsCol);

         for(int i = 0; i < nnormsRow; ++i)
            weights[i] = norms[i];

         for(int i = 0; i < nnormsCol; ++i)
            coWeights[i] = norms[nnormsRow + i];

         weightsAreSetup = true;
      }
      else
         return false;

      return true;
   }

   template <>
   void SPxSolverBase<Real>::setIntegralityInformation(int ncols, int* intInfo)
   {
      assert(ncols == this->nCols() || (ncols == 0 && intInfo == NULL));

      integerVariables.reSize(ncols);

      for(int i = 0; i < ncols; ++i)
      {
         integerVariables[i] = intInfo[i];
      }
   }



   //
   // Auxiliary functions.
   //

   // Pretty-printing of variable status.
   template <class R>
   std::ostream& operator<<(std::ostream & os,
                            const typename SPxSolverBase<Real>::VarStatus & status)
   {
      switch(status)
      {
      case SPxSolverBase<Real>::BASIC:
         os << "BASIC";
         break;

      case SPxSolverBase<Real>::FIXED:
         os << "FIXED";
         break;

      case SPxSolverBase<Real>::ON_LOWER:
         os << "ON_LOWER";
         break;

      case SPxSolverBase<Real>::ON_UPPER:
         os << "ON_UPPER";
         break;

      case SPxSolverBase<Real>::ZERO:
         os << "ZERO";
         break;

      case SPxSolverBase<Real>::UNDEFINED:
         os << "UNDEFINED";
         break;

      default:
         os << "?invalid?";
         break;
      }

      return os;
   }

   // Pretty-printing of solver status.
   template <class R>
   std::ostream& operator<<(std::ostream & os,
                            const typename SPxSolverBase<R>::Status & status)
   {
      switch(status)
      {
      case SPxSolverBase<R>::ERROR:
         os << "ERROR";
         break;

      case SPxSolverBase<R>::NO_RATIOTESTER:
         os << "NO_RATIOTESTER";
         break;

      case SPxSolverBase<R>::NO_PRICER:
         os << "NO_PRICER";
         break;

      case SPxSolverBase<R>::NO_SOLVER:
         os << "NO_SOLVER";
         break;

      case SPxSolverBase<R>::NOT_INIT:
         os << "NOT_INIT";
         break;

      case SPxSolverBase<R>::ABORT_CYCLING:
         os << "ABORT_CYCLING";
         break;

      case SPxSolverBase<R>::ABORT_TIME:
         os << "ABORT_TIME";
         break;

      case SPxSolverBase<R>::ABORT_ITER:
         os << "ABORT_ITER";
         break;

      case SPxSolverBase<R>::ABORT_VALUE:
         os << "ABORT_VALUE";
         break;

      case SPxSolverBase<R>::SINGULAR:
         os << "SINGULAR";
         break;

      case SPxSolverBase<R>::NO_PROBLEM:
         os << "NO_PROBLEM";
         break;

      case SPxSolverBase<R>::REGULAR:
         os << "REGULAR";
         break;

      case SPxSolverBase<R>::RUNNING:
         os << "RUNNING";
         break;

      case SPxSolverBase<R>::UNKNOWN:
         os << "UNKNOWN";
         break;

      case SPxSolverBase<R>::OPTIMAL:
         os << "OPTIMAL";
         break;

      case SPxSolverBase<R>::UNBOUNDED:
         os << "UNBOUNDED";
         break;

      case SPxSolverBase<R>::INFEASIBLE:
         os << "INFEASIBLE";
         break;

      default:
         os << "?other?";
         break;
      }

      return os;
   }

   // Pretty-printing of algorithm.
   template <class R>
   std::ostream& operator<<(std::ostream & os,
                            const typename SPxSolverBase<R>::Type & status)
   {
      switch(status)
      {
      case SPxSolverBase<R>::ENTER:
         os << "ENTER";
         break;

      case SPxSolverBase<R>::LEAVE:
         os << "LEAVE";
         break;

      default:
         os << "?other?";
         break;
      }

      return os;
   }

   // Pretty-printing of representation.
   template <class R>
   std::ostream& operator<<(std::ostream & os,
                            const typename SPxSolverBase<R>::Representation & status)
   {
      switch(status)
      {
      case SPxSolverBase<R>::ROW:
         os << "ROW";
         break;

      case SPxSolverBase<R>::COLUMN:
         os << "COLUMN";
         break;

      default:
         os << "?other?";
         break;
      }

      return os;
   }


} // namespace soplex
