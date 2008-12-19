/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the class library                   */
/*       SoPlex --- the Sequential object-oriented simPlex.                  */
/*                                                                           */
/*    Copyright (C) 1997-1999 Roland Wunderling                              */
/*                  1997-2008 Konrad-Zuse-Zentrum                            */
/*                            fuer Informationstechnik Berlin                */
/*                                                                           */
/*  SoPlex is distributed under the terms of the ZIB Academic Licence.       */
/*                                                                           */
/*  You should have received a copy of the ZIB Academic License              */
/*  along with SoPlex; see the file COPYING. If not email to soplex@zib.de.  */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
#pragma ident "@(#) $Id: spxsolver.cpp,v 1.37 2008/12/19 10:01:51 bzfgleix Exp $"

//#define DEBUGGING 1

#include <assert.h>
#include <iostream>
//#include <fstream>

#include "spxdefines.h"
#include "soplex.h"
#include "spxpricer.h"
#include "spxratiotester.h"
#include "spxstarter.h"
#include "spxout.h"
#include "exceptions.h"


namespace soplex
{
#define MAXIMUM(x,y)        ((x)>(y) ? (x) : (y))

bool SPxSolver::read(std::istream& in, NameSet* rowNames, 
                  NameSet* colNames, DIdxSet* intVars)
{
   METHOD( "SPxSolver::read()" );
   clear();
   unInit();
   unLoad();

   if (thepricer)
      thepricer->clear();

   if (theratiotester)
      theratiotester->clear();

   if (!SPxLP::read(in, rowNames, colNames, intVars))
      return false;

   SPxBasis::load(this);

   return true;
}

void SPxSolver::reLoad()
{
   METHOD( "SPxSolver::reLoad()" );
   unInit();
   unLoad();
   theLP = this;
   if (thepricer)
      thepricer->clear();
   if (theratiotester)
      theratiotester->clear();
}

void SPxSolver::loadLP(const SPxLP& lp)
{
   METHOD( "SPxSolver::loadLP()" );
   clear();
   unInit();
   unLoad();
   if (thepricer)
      thepricer->clear();
   if (theratiotester)
      theratiotester->clear();
   SPxLP::operator=(lp);
   reDim();
   SPxBasis::load(this);
}

void SPxSolver::setSolver(SLinSolver* slu)
{
   METHOD( "SPxSolver::setSolver()" );
   SPxBasis::loadSolver(slu);
}

void SPxSolver::loadBasis(const SPxBasis::Desc& p_desc)
{
   METHOD( "SPxSolver::loadBasis()" );
   unInit();
   if (SPxBasis::status() == SPxBasis::NO_PROBLEM)
      SPxBasis::load(this);
   SPxBasis::loadDesc(p_desc);
}

void SPxSolver::setPricer(SPxPricer* x)
{
   METHOD( "SPxSolver::setPricer()" );
   if (x != 0)
   {
      setPricing(FULL);
      if (isInitialized() && x != thepricer)
         x->load(this);
      else
         x->clear();
   }
   if (thepricer && thepricer != x)
      thepricer->clear();
   thepricer = x;
}

void SPxSolver::setTester(SPxRatioTester* x)
{
   METHOD( "SPxSolver::setTester()" );
   if (x)
   {
      if (isInitialized() && x != theratiotester)
         x->load(this);
      else
         x->clear();
   }
   if (theratiotester && theratiotester != x)
      theratiotester->clear();
   theratiotester = x;
}

void SPxSolver::setStarter(SPxStarter* x)
{
   METHOD( "SPxSolver::setStarter()" );
   thestarter = x;
}

void SPxSolver::setType(Type tp)
{
   METHOD( "SPxSolver::setType()" );

   if (theType != tp)
   {
      theType = tp;

      unInit();
#if 0
      else
      {
      if (!matrixIsSetup)
      {
         SPxBasis::load(this);
         // SPxBasis::load(desc());
         // not needed, because load(this) allready loads descriptor
      }
      factorized = false;
      m_numCycle = 0;
#endif
      MSG_DEBUG( spxout << "DSOLVE20 switching to " 
                        << static_cast<const char*>((tp == LEAVE)
                           ? "leaving" : "entering")
                        << " algorithm" << std::endl; )
   }
}

void SPxSolver::setRep(Representation p_rep)
{
   METHOD( "SPxSolver::setRep()" );
   if (p_rep == COLUMN)
   {
      thevectors   = colSet();
      thecovectors = rowSet(); 
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
      assert(p_rep == ROW);

      thevectors   = rowSet(); 
      thecovectors = colSet();
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
   theRep = p_rep;
   unInit();
   reDim();

   if (SPxBasis::status() > SPxBasis::NO_PROBLEM)
   {
      SPxBasis::setRep();
      SPxBasis::loadDesc(desc());
   }
   if (thepricer && thepricer->solver() == this)
      thepricer->setRep(p_rep);
}

// needed for strongbranching. use carefully
void SPxSolver::reinitializeVecs()
{
   METHOD( "SPxSolver::reinitializeVecs" );
   
   initialized = true;

   if (type() == ENTER)
   {
      if (rep() == COLUMN)
         setPrimalBounds();
      else
         setDualRowBounds();

      setEnterBounds();
      computeEnterCoPrhs();
   }
   else
   {
      if (rep() == ROW)
         setPrimalBounds();
      else
         setDualColBounds();

      setLeaveBounds();
      computeLeaveCoPrhs();
   }

   SPxBasis::coSolve(*theCoPvec, *theCoPrhs);
   computePvec();
   computeFrhs();
   SPxBasis::solve(*theFvec, *theFrhs);

   theShift  = 0.0;
   lastShift = 0.0;

   if (type() == ENTER)
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

void SPxSolver::init()
{
   METHOD( "SPxSolver::init()" );

   assert(thepricer      != 0);
   assert(theratiotester != 0);

   if (!initialized)
   {
      initialized = true;
      reDim();
      if (SPxBasis::status() <= SPxBasis::NO_PROBLEM || solver() != this)
         SPxBasis::load(this);
      initialized = false;
   }
   if (!matrixIsSetup)
      SPxBasis::loadDesc(desc());

   // Inna/Tobi: don't "upgrade" a singular basis to a regular one
   if( SPxBasis::status() == SPxBasis::SINGULAR )
      return;

   //factorized = false;
   m_numCycle = 0;

   if (type() == ENTER)
   {
      if (rep() == COLUMN)
      {
         setPrimalBounds();
         setBasisStatus(SPxBasis::PRIMAL);
      }
      else
      {
         setDualRowBounds();
         setBasisStatus(SPxBasis::DUAL);
      }
      setEnterBounds();
      computeEnterCoPrhs();
   }
   else
   {
      if (rep() == ROW)
      {
         setPrimalBounds();
         setBasisStatus(SPxBasis::PRIMAL);
      }
      else
      {
         setDualColBounds();
         setBasisStatus(SPxBasis::DUAL);
      }
      setLeaveBounds();
      computeLeaveCoPrhs();
   }

   SPxBasis::coSolve(*theCoPvec, *theCoPrhs);
   computePvec();
   computeFrhs();
   SPxBasis::solve(*theFvec, *theFrhs);

   theShift = 0.0;

   if (type() == ENTER)
   {
      shiftFvec();
      computeCoTest();
      computeTest();
   }
   else
   {
      shiftPvec();
      computeFtest();
   }
   lastShift = theShift + delta();

   if (!initialized)
   {
      // if(thepricer->solver() != this)
      thepricer->load(this);
      // if(theratiotester->solver() != this)
      theratiotester->load(this);
      initialized = true;
   }
}

void SPxSolver::setPricing(Pricing pr)
{
   METHOD( "SPxSolver::setPricing()" );
   thePricing = pr;
   if (initialized && type() == ENTER)
   {
      computePvec();
      computeCoTest();
      computeTest();
   }
}

/*
    The following method resizes all vectors and arrays of |SoPlex|
    (excluding inherited vectors).
 */
void SPxSolver::reDim()
{
   METHOD( "SPxSolver::reDim()" );

   int newdim = (rep() == ROW) ? SPxLP::nCols() : SPxLP::nRows();

   if (newdim > unitVecs.size())
   {
      unitVecs.reSize(newdim);

      while(newdim-- > 0)
         unitVecs[newdim] = UnitVector( newdim );
   }

   if (isInitialized())
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

void SPxSolver::clear()
{
   METHOD( "SPxSolver::clear()" );
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

   unInit();
   SPxLP::clear();
   setBasisStatus(SPxBasis::NO_PROBLEM);
   SPxBasis::reDim();
}

void SPxSolver::clearUpdateVecs(void)
{
   METHOD( "SPxSolver::clearUpdateVecs()" );
   theFvec->clearUpdate();
   thePvec->clearUpdate();
   theCoPvec->clearUpdate();
   solveVector2 = 0;
   coSolveVector2 = 0;
}

/*
    When the basis matrix factorization is recomputed from scratch, 
    we also recompute the vectors.
 */
void SPxSolver::factorize()
{
   METHOD( "SPxSolver::factorize()" );

   MSG_INFO1( spxout << "ISOLVE01 " 
                     << "iteration = "    << std::setw(8) << basis().iteration() 
                     << "\tlastUpdate = " << std::setw(4) << basis().lastUpdate()
                     << "\tvalue = "      << (isInitialized() ? value() : 0.0) 
                     << std::endl; )


   SPxBasis::factorize();

   if (SPxBasis::status() >= SPxBasis::REGULAR)
   {
      // #undef       NDEBUG
#ifndef NDEBUG
      DVector ftmp(fVec());
      DVector ptmp(pVec());
      DVector ctmp(coPvec());
      /* move this test after the computation of fTest and coTest below, since these vectors might not be set up, e.g., for an initial basis. */
      // testVecs();
#endif  // NDEBUG

      if (type() == LEAVE)
      {
         SPxBasis::solve (*theFvec, *theFrhs);
         SPxBasis::coSolve(*theCoPvec, *theCoPrhs);
      }

#ifndef NDEBUG
      ftmp -= fVec();
      ptmp -= pVec();
      ctmp -= coPvec();
      if (ftmp.length() > delta())
      {
         MSG_DEBUG( spxout << "DSOLVE21 fVec:   " << ftmp.length() << std::endl; )
         ftmp = fVec();
         multBaseWith(ftmp);
         ftmp -= fRhs();
         if (ftmp.length() > delta())
            MSG_ERROR( spxout << "ESOLVE29 " << iteration() << ": fVec error = " 
                            << ftmp.length() << std::endl; )
      }
      if (ctmp.length() > delta())
      {
         MSG_DEBUG( spxout << "DSOLVE23 coPvec: " << ctmp.length() << std::endl; )
         ctmp = coPvec();
         multWithBase(ctmp);
         ctmp -= coPrhs();
         if (ctmp.length() > delta())
            MSG_ERROR( spxout << "ESOLVE30 " << iteration() << ": coPvec error = " 
                            << ctmp.length() << std::endl; )
      }
      if (ptmp.length() > delta())
      {
         MSG_DEBUG( spxout << "DSOLVE24 pVec:   " << ptmp.length() << std::endl; )
      }
#endif  // NDEBUG

      if (type() == ENTER)
      {
         computeCoTest();
//         if(pricing() == FULL)
//         {
//             computePvec();
//             computeTest();
//         }
      }
      else
      {
         computeFtest();
//          computePvec();
      }
   }

#ifndef NDEBUG
   /* see note near top of function */
   testVecs();
#endif  // NDEBUG

   if (SPxBasis::status() == SPxBasis::SINGULAR)
   {
      m_status = SINGULAR;
      throw SPxStatusException("XSOLVE21 Basis is singular (numerical troubles)");
   }
}

/* We compute how much the current solution violates (primal or dual) feasibility. In the
   row/enter or column/leave algorithm the maximum violation of dual feasibility is
   computed. In the row/leave or column/enter algorithm the primal feasibility is checked. */
Real SPxSolver::maxInfeas() const
{
   METHOD( "SPxSolver::maxInfeas()" );
   Real inf = 0.0;

   if (type() == ENTER)
   {
      for (int i = 0; i < dim(); i++)
      {
         if ((*theFvec)[i] > theUBbound[i])
            inf = MAXIMUM(inf, (*theFvec)[i] - theUBbound[i]);
         if (theLBbound[i] > (*theFvec)[i])
            inf = MAXIMUM(inf, theLBbound[i] - (*theFvec)[i]);
      }
   }
   else
   {
      assert(type() == LEAVE);

      for (int i = 0; i < dim(); i++)
      {
         if ((*theCoPvec)[i] > (*theCoUbound)[i])
            inf = MAXIMUM(inf, (*theCoPvec)[i] - (*theCoUbound)[i]);
         if ((*theCoLbound)[i] > (*theCoPvec)[i])
            inf = MAXIMUM(inf, (*theCoLbound)[i] - (*theCoPvec)[i]);
      }
      for (int i = 0; i < coDim(); i++)
      {
         if ((*thePvec)[i] > (*theUbound)[i])
            inf = MAXIMUM(inf, (*thePvec)[i] - (*theUbound)[i]);
         else if ((*thePvec)[i] < (*theLbound)[i])
            inf = MAXIMUM(inf, (*theLbound)[i] - (*thePvec)[i]);
      }
   }

   return inf;
}

Real SPxSolver::nonbasicValue() const
{
   METHOD( "SPxSolver::nonbasicValue()" );

   int i;
   Real val = 0;
   const SPxBasis::Desc& ds = desc();

   if (rep() == COLUMN)
   {
      if (type() == LEAVE)
      {
         for (i = nCols() - 1; i >= 0; --i)
         {
            switch (ds.colStatus(i))
            {
            case SPxBasis::Desc::P_ON_UPPER :
               val += theUCbound[i] * SPxLP::upper(i);
               //@ val += maxObj(i) * SPxLP::upper(i);
               break;
            case SPxBasis::Desc::P_ON_LOWER :
               val += theLCbound[i] * SPxLP::lower(i);
               //@ val += maxObj(i) * SPxLP::lower(i);
               break;
            case SPxBasis::Desc::P_ON_UPPER + SPxBasis::Desc::P_ON_LOWER :
               val += maxObj(i) * SPxLP::lower(i);
               break;
            default:
               break;
            }
         }
         for (i = nRows() - 1; i >= 0; --i)
         {
            switch (ds.rowStatus(i))
            {
            case SPxBasis::Desc::P_ON_UPPER :
               val += theLRbound[i] * SPxLP::rhs(i);
               break;
            case SPxBasis::Desc::P_ON_LOWER :
               val += theURbound[i] * SPxLP::lhs(i);
               break;
            default:
               break;
            }
         }
      }
      else
      {
         assert(type() == ENTER);
         for (i = nCols() - 1; i >= 0; --i)
         {
            switch (ds.colStatus(i))
            {
            case SPxBasis::Desc::P_ON_UPPER :
               val += maxObj(i) * theUCbound[i];
               break;
            case SPxBasis::Desc::P_ON_LOWER :
               val += maxObj(i) * theLCbound[i];
               break;
            case SPxBasis::Desc::P_ON_UPPER + SPxBasis::Desc::P_ON_LOWER :
               assert(theLCbound[i] == theUCbound[i]);
               val += maxObj(i) * theLCbound[i];
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
      for (i = nCols() - 1; i >= 0; --i)
      {
         switch (ds.colStatus(i))
         {
         case SPxBasis::Desc::D_ON_UPPER :
            val += theUCbound[i] * lower(i);
            break;
         case SPxBasis::Desc::D_ON_LOWER :
            val += theLCbound[i] * upper(i);
            break;
         case SPxBasis::Desc::D_ON_BOTH :
            val += theLCbound[i] * upper(i);
            val += theUCbound[i] * lower(i);
            break;
         default:
            break;
         }
      }
      for (i = nRows() - 1; i >= 0; --i)
      {
         switch (ds.rowStatus(i))
         {
         case SPxBasis::Desc::D_ON_UPPER :
            val += theURbound[i] * lhs(i);
            break;
         case SPxBasis::Desc::D_ON_LOWER :
            val += theLRbound[i] * rhs(i);
            break;
         case SPxBasis::Desc::D_ON_BOTH :
            val += theLRbound[i] * rhs(i);
            val += theURbound[i] * lhs(i);
            break;
         default:
            break;
         }
      }
   }

   return val;
}

Real SPxSolver::value() const
{
   METHOD( "SPxSolver::value()" );
   Real x;

   assert(isInitialized());

   // calling value() without having a suitable status is an error.
   if (!isInitialized())
      return infinity;

   if (rep() == ROW)
   {
      if (type() == LEAVE)
         x = SPxLP::spxSense() * (coPvec() * fRhs());
      else
         x = SPxLP::spxSense() * (nonbasicValue() + (coPvec() * fRhs()));
   }
   else
      x = SPxLP::spxSense() * (nonbasicValue() + fVec() * coPrhs());

   return x;
}

void SPxSolver::setDelta(Real d)
{
   METHOD( "SPxSolver::setDelta()" );
   theDelta = d;
}

SPxSolver::SPxSolver(
   Type            p_type, 
   Representation  p_rep,
   SPxPricer*      pric, 
   SPxRatioTester* rt,
   SPxStarter*     start)
   : theType (p_type)
   , thePricing(FULL)
   , maxIters (-1)
   , maxTime (infinity)
   , objLimit(infinity)
   , m_status(UNKNOWN)
   , theShift (0)
   , m_maxCycle(100)
   , m_numCycle(0)
   , initialized (false)
   , solveVector2 (0)
   , coSolveVector2(0)
   , unitVecs (0)
   , primVec (0, Param::epsilon())
   , dualVec (0, Param::epsilon())
   , addVec (0, Param::epsilon())
   , thepricer (pric)
   , theratiotester(rt)
   , thestarter (start)
{
   METHOD( "SPxSolver::SPxSolver()" );
   setRep (p_rep);
   setDelta (DEFAULT_BND_VIOL);
   theLP = this;
}

#ifndef NO_CONSISTENCY_CHECKS
bool SPxSolver::isConsistent() const
{
   METHOD( "SPxSolver::isConsistent()" );
   if (epsilon() < 0)
      return MSGinconsistent("SPxSolver");

   if (primVec.delta().getEpsilon() != dualVec.delta().getEpsilon())
      return MSGinconsistent("SPxSolver");

   if (dualVec.delta().getEpsilon() != addVec.delta().getEpsilon())
      return MSGinconsistent("SPxSolver");

   if (unitVecs.size() < ((rep() == ROW) ? SPxLP::nCols() : SPxLP::nRows()))
      return MSGinconsistent("SPxSolver");

   if (initialized)
   {
      if (theFrhs->dim() != dim())
         return MSGinconsistent("SPxSolver");
      if (theFvec->dim() != dim())
         return MSGinconsistent("SPxSolver");

      if (theCoPrhs->dim() != dim())
         return MSGinconsistent("SPxSolver");
      if (thePvec->dim() != coDim())
         return MSGinconsistent("SPxSolver");
      if (theCoPvec->dim() != dim())
         return MSGinconsistent("SPxSolver");

      if (theTest.dim() != coDim())
         return MSGinconsistent("SPxSolver");
      if (theCoTest.dim() != dim())
         return MSGinconsistent("SPxSolver");

      if (theURbound.dim() != SPxLP::nRows())
         return MSGinconsistent("SPxSolver");
      if (theLRbound.dim() != SPxLP::nRows())
         return MSGinconsistent("SPxSolver");
      if (theUCbound.dim() != SPxLP::nCols())
         return MSGinconsistent("SPxSolver");
      if (theLCbound.dim() != SPxLP::nCols())
         return MSGinconsistent("SPxSolver");
      if (theUBbound.dim() != dim())
         return MSGinconsistent("SPxSolver");
      if (theLBbound.dim() != dim())
         return MSGinconsistent("SPxSolver");
   }

   if (rep() == COLUMN)
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
         return MSGinconsistent("SPxSolver");
   }
   else
   {
      if (thecovectors 
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
         return MSGinconsistent("SPxSolver");
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
          && SPxBasis::isConsistent()
         ;
}
#endif


void SPxSolver::setTerminationTime(Real p_time)
{
   METHOD( "SPxSolver::setTerminationTime()" );
   if( p_time < 0.0 )
      p_time = infinity;
   maxTime = p_time;
}

Real SPxSolver::terminationTime() const
{
   METHOD( "SPxSolver::terminationTime()" );
   return maxTime;
}

void SPxSolver::setTerminationIter(int p_iteration)
{
   METHOD( "SPxSolver::setTerminationIter()" );
   if( p_iteration < 0 )
      p_iteration = -1;
   maxIters = p_iteration;
}

int SPxSolver::terminationIter() const
{
   METHOD( "SPxSolver::terminationIter()" );
   return maxIters;
}

/**@todo A first version for the termination value is
 *       implemented. Currently we check if no bound violations (shifting)
 *       is present. It might be even possible to use this termination
 *       value in case of bound violations (shifting) but in this case it
 *       is quite difficult to determine if we already reached the limit.
 */
void SPxSolver::setTerminationValue(Real p_value)
{
   METHOD( "SPxSolver::setTerminationValue()" );
   objLimit = p_value;
}

Real SPxSolver::terminationValue() const
{
   METHOD( "SPxSolver::terminationValue()" );
   return objLimit;
}
   
SPxSolver::VarStatus
SPxSolver::basisStatusToVarStatus( SPxBasis::Desc::Status stat ) const
{
   METHOD( "SPxSolver::VarStatus()" );
   VarStatus vstat;

   switch( stat )
   {
   case SPxBasis::Desc::P_ON_LOWER:
      vstat = ON_LOWER;
      break;
   case SPxBasis::Desc::P_ON_UPPER:
      vstat = ON_UPPER;
      break;
   case SPxBasis::Desc::P_FIXED:
      vstat = FIXED;
      break;
   case SPxBasis::Desc::P_FREE:
      vstat = ZERO;
      break;
   case SPxBasis::Desc::D_ON_UPPER:
   case SPxBasis::Desc::D_ON_LOWER:
   case SPxBasis::Desc::D_ON_BOTH:
   case SPxBasis::Desc::D_UNDEFINED:
   case SPxBasis::Desc::D_FREE:
      vstat = BASIC;
      break;
   default:
      MSG_ERROR( spxout << "ESOLVE26 ERROR: unknown basis status (" << stat << ")" 
                        << std::endl; )
      throw SPxInternalCodeException("XSOLVE22 This should never happen.");
   }
   return vstat;
}

SPxBasis::Desc::Status
SPxSolver::varStatusToBasisStatusRow( int row, SPxSolver::VarStatus stat ) const
{
   METHOD( "SPxSolver::varStatusToBasisStatusRow()" );
   SPxBasis::Desc::Status rstat;

   switch( stat )
   {
   case FIXED :
      assert(rhs(row) == lhs(row));
      rstat = SPxBasis::Desc::P_FIXED;
      break;
   case ON_UPPER :
      assert(rhs(row) < infinity);
      rstat = lhs(row) < rhs(row)
         ? SPxBasis::Desc::P_ON_UPPER
         : SPxBasis::Desc::P_FIXED;
      break;
   case ON_LOWER :
      assert(lhs(row) > -infinity);
      rstat = lhs(row) < rhs(row)
         ? SPxBasis::Desc::P_ON_LOWER
         : SPxBasis::Desc::P_FIXED;
      break;
   case ZERO :
      /* A 'free' row (i.e., infinite lower & upper bounds) does not really make sense. The user
       * might (think to) know better, e.g., when temporarily turning off a row. We therefore apply
       * the same adjustment as in the column case in varStatusToBasisStatusCol(). */
      if (lhs(row) <= -infinity && rhs(row) >= infinity)
         rstat = SPxBasis::Desc::P_FREE;
      else
      {
         if ( lhs(row) == rhs(row) )
         {
            assert( rhs(row) < infinity );
            rstat = SPxBasis::Desc::P_FIXED;
         }
         else
         {
            if ( lhs(row) > -infinity )
               rstat = SPxBasis::Desc::P_ON_LOWER;
            else
            {
               assert( rhs(row) < infinity );
               rstat = SPxBasis::Desc::P_ON_UPPER;
            }
         }
      }
      break;
   case BASIC :
      rstat = dualRowStatus(row);
      break;
   default:
      MSG_ERROR( spxout << "ESOLVE27 ERROR: unknown VarStatus (" << int(stat) << ")" 
                        << std::endl; )
      throw SPxInternalCodeException("XSOLVE23 This should never happen.");
   }
   return rstat;
}

SPxBasis::Desc::Status 
SPxSolver::varStatusToBasisStatusCol( int col, SPxSolver::VarStatus stat ) const
{
   METHOD( "SPxSolver::varStatusToBasisStatusCol()" );
   SPxBasis::Desc::Status cstat;

   switch( stat )
   {
   case FIXED :
      assert(upper(col) == lower(col));
      cstat = SPxBasis::Desc::P_FIXED;
      break;
   case ON_UPPER :
      assert(upper(col) < infinity);
      cstat = lower(col) < upper(col)
         ? SPxBasis::Desc::P_ON_UPPER
         : SPxBasis::Desc::P_FIXED;
      break;
   case ON_LOWER :
      assert(lower(col) > -infinity);
      cstat = lower(col) < upper(col)
         ? SPxBasis::Desc::P_ON_LOWER
         : SPxBasis::Desc::P_FIXED;
      break;
   case ZERO :
      /* In this case the upper and lower bounds on the variable should be infinite. The bounds
       * might, however, have changed and we try to recover from this by changing the status to
       * 'resonable' settings. A solve has to find the correct values afterwards. Note that the
       * approach below is consistent with changesoplex.cpp (e.g., changeUpperStatus() and
       * changeLowerStatus() ). */
      if (lower(col) <= -infinity && upper(col) >= infinity)
         cstat = SPxBasis::Desc::P_FREE;
      else
      {
         if ( lower(col) == upper(col) )
         {
            assert( upper(col) < infinity );
            cstat = SPxBasis::Desc::P_FIXED;
         }
         else
         {
            if ( lower(col) > -infinity )
               cstat = SPxBasis::Desc::P_ON_LOWER;
            else
            {
               assert( upper(col) < infinity );
               cstat = SPxBasis::Desc::P_ON_UPPER;
            }
         }
      }
      break;
   case BASIC :
      cstat = dualColStatus(col);
      break;
   default:
      MSG_ERROR( spxout << "ESOLVE28 ERROR: unknown VarStatus (" << int(stat) << ")" 
                        << std::endl; )
      throw SPxInternalCodeException("XSOLVE24 This should never happen.");
   }
   return cstat;
}

SPxSolver::VarStatus SPxSolver::getBasisRowStatus( int row ) const
{
   METHOD( "SPxSolver::VarStatus()" );
   assert( 0 <= row && row < nRows() );
   return basisStatusToVarStatus( desc().rowStatus( row ) );
}

SPxSolver::VarStatus SPxSolver::getBasisColStatus( int col ) const
{
   METHOD( "SPxSolver::VarStatus()" );
   assert( 0 <= col && col < nCols() );
   return basisStatusToVarStatus( desc().colStatus( col ) );
}

SPxSolver::Status SPxSolver::getBasis(VarStatus row[], VarStatus col[]) const
{
   METHOD( "SPxSolver::Status()" );
   const SPxBasis::Desc& d = desc();
   int i;

   if (col)
      for (i = nCols() - 1; i >= 0; --i)
         col[i] = basisStatusToVarStatus( d.colStatus(i) );

   if (row)
      for (i = nRows() - 1; i >= 0; --i)
         row[i] = basisStatusToVarStatus( d.rowStatus(i) );

   return status();
}

void SPxSolver::setBasis(const VarStatus p_rows[], const VarStatus p_cols[])
{
   METHOD( "SPxSolver::setBasis()" );
   if (SPxBasis::status() == SPxBasis::NO_PROBLEM)
      SPxBasis::load(this);

   SPxBasis::Desc ds = desc();
   int i;

   for(i = 0; i < nRows(); i++)
      ds.rowStatus(i) = varStatusToBasisStatusRow( i, p_rows[i] );

   for(i = 0; i < nCols(); i++)
      ds.colStatus(i) = varStatusToBasisStatusCol( i, p_cols[i] );

   loadBasis(ds);
}

//
// Auxiliary functions.
//

// Pretty-printing of variable status.
std::ostream& operator<<( std::ostream& os,
                          const SPxSolver::VarStatus& status )
{
   switch( status )
      {
      case SPxSolver::BASIC:
         os << "BASIC";
         break;
      case SPxSolver::FIXED:
         os << "FIXED";
         break;
      case SPxSolver::ON_LOWER:
         os << "ON_LOWER";
         break;
      case SPxSolver::ON_UPPER:
         os << "ON_UPPER";
         break;
      case SPxSolver::ZERO:
         os << "ZERO";
         break;
      default:
         os << "?invalid?";
         break;
      }
   return os;
}

// Pretty-printing of solver status.
std::ostream& operator<<( std::ostream& os,
                          const SPxSolver::Status& status )
{
   switch ( status )
      {
      case SPxSolver::ERROR:
         os << "ERROR";
         break;
      case SPxSolver::NO_RATIOTESTER:
         os << "NO_RATIOTESTER";
         break;
      case SPxSolver::NO_PRICER:
         os << "NO_PRICER";
         break;
      case SPxSolver::NO_SOLVER:
         os << "NO_SOLVER";
         break;
      case SPxSolver::NOT_INIT:
         os << "NOT_INIT";
         break;
      case SPxSolver::ABORT_CYCLING:
         os << "ABORT_CYCLING";
         break;
      case SPxSolver::ABORT_TIME:
         os << "ABORT_TIME";
         break;
      case SPxSolver::ABORT_ITER:
         os << "ABORT_ITER";
         break;
      case SPxSolver::ABORT_VALUE:
         os << "ABORT_VALUE";
         break;
      case SPxSolver::SINGULAR:
         os << "SINGULAR";
         break;
      case SPxSolver::NO_PROBLEM:
         os << "NO_PROBLEM";
         break;
      case SPxSolver::REGULAR:
         os << "REGULAR";
         break;
      case SPxSolver::RUNNING:
         os << "RUNNING";
         break;
      case SPxSolver::UNKNOWN:
         os << "UNKNOWN";
         break;
      case SPxSolver::OPTIMAL:
         os << "OPTIMAL";
         break;
      case SPxSolver::UNBOUNDED:
         os << "UNBOUNDED";
         break;
      case SPxSolver::INFEASIBLE:
         os << "INFEASIBLE";
         break;
      default:
         os << "?other?";
         break;
      }
   return os;
}

// Pretty-printing of algorithm.
std::ostream& operator<<( std::ostream& os,
                          const SPxSolver::Type& status )
{
   switch ( status )
      {
      case SPxSolver::ENTER:
         os << "ENTER";
         break;
      case SPxSolver::LEAVE:
         os << "LEAVE";
         break;
      default:
         os << "?other?";
         break;
      }
   return os;
}

// Pretty-printing of representation.
std::ostream& operator<<( std::ostream& os,
                          const SPxSolver::Representation& status )
{
   switch ( status )
      {
      case SPxSolver::ROW:
         os << "ROW";
         break;
      case SPxSolver::COLUMN:
         os << "COLUMN";
         break;
      default:
         os << "?other?";
         break;
      }
   return os;
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
