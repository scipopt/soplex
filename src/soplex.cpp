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
#pragma ident "@(#) $Id: soplex.cpp,v 1.27 2002/01/12 11:41:25 bzfkocht Exp $"

#include <assert.h>
#include <iostream>
#include <fstream>

#include "soplex.h"
#include "spxpricer.h"
#include "spxratiotester.h"
#include "spxstarter.h"
#include "spxsimplifier.h"

namespace soplex
{
#define MAXIMUM(x,y)        ((x)>(y) ? (x) : (y))

bool SoPlex::read(std::istream& in, NameSet* rowNames, NameSet* colNames)
{
   clear();
   unInit();
   unLoad();

   if (thepricer)
      thepricer->clear();

   if (theratiotester)
      theratiotester->clear();

   if (!SPxLP::read(in, rowNames, colNames))
      return false;

   SPxBasis::load(this);

   int tmp = coDim() / (20 * dim()) + 1;
   coVecDim = coDim() / tmp + 1;

   return true;
}

void SoPlex::reLoad()
{
   unInit();
   unLoad();
   theLP = this;
   if (thepricer)
      thepricer->clear();
   if (theratiotester)
      theratiotester->clear();
}

void SoPlex::loadLP(const SPxLP& lp)
{
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

void SoPlex::setSolver(SLinSolver* slu)
{
   SPxBasis::load(slu);
}

void SoPlex::loadBasis(const SPxBasis::Desc& p_desc)
{
   unInit();
   if (SPxBasis::status() == SPxBasis::NO_PROBLEM)
      SPxBasis::load(this);
   SPxBasis::load(p_desc);
}

void SoPlex::setPricer(SPxPricer* x)
{
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

void SoPlex::setTester(SPxRatioTester* x)
{
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

void SoPlex::setStarter(SPxStarter* x)
{
   thestarter = x;
}

void SoPlex::setSimplifier(SPxSimplifier* x)
{
   thesimplifier = x;
}


void SoPlex::setType(Type tp)
{
   if (isInitialized() && theType != tp)
   {
      theType = tp;
      init();
   }
   else
   {
      theType = tp;

      if (!matrixIsSetup)
      {
         SPxBasis::load(this);
         SPxBasis::load(desc());
      }
      factorized = false;
      m_numCycle = 0;
   }
   if ((thepricer != 0) && (thepricer->solver() == this))
      thepricer->setType(tp);
   if ((theratiotester != 0) && (theratiotester->solver() == this))
      theratiotester->setType(tp);

   std::cout << "switching to " 
             << static_cast<const char*>
                ((tp == LEAVE) ? "leaving" : "entering")
             << " algorithm" 
             << std::endl;
}

void SoPlex::setRep(int rp)
{
   Representation l_rep = Representation(rp);
   if (l_rep == COLUMN)
   {
      thevectors = colset();
      thecovectors = rowset();
      theFrhs = &primRhs;
      theFvec = &primVec;
      theCoPrhs = &dualRhs;
      theCoPvec = &dualVec;
      thePvec = &addVec;
      theRPvec = theCoPvec;
      theCPvec = thePvec;
      theUbound = &theUCbound;
      theLbound = &theLCbound;
      theCoUbound = &theURbound;
      theCoLbound = &theLRbound;
   }
   else
   {
      assert(l_rep == ROW);

      thevectors = rowset();
      thecovectors = colset();
      theFrhs = &dualRhs;
      theFvec = &dualVec;
      theCoPrhs = &primRhs;
      theCoPvec = &primVec;
      thePvec = &addVec;
      theRPvec = thePvec;
      theCPvec = theCoPvec;
      theUbound = &theURbound;
      theLbound = &theLRbound;
      theCoUbound = &theUCbound;
      theCoLbound = &theLCbound;
   }
   therep = l_rep;
   unInit();
   reDim();

   if (SPxBasis::status() > SPxBasis::NO_PROBLEM)
   {
      SPxBasis::setRep();
      SPxBasis::load(desc());
   }

   if (thepricer && thepricer->solver() == this)
      thepricer->setRep(l_rep);
}

void SoPlex::init()
{
   assert(thepricer != 0);
   assert(theratiotester != 0);

   if (!initialized)
   {
      initialized = true;
      reDim();
      if (SPxBasis::status() <= NO_PROBLEM || solver() != this)
         SPxBasis::load(this);
      initialized = false;
   }
   if (!matrixIsSetup)
      SPxBasis::load(desc());
   factorized = false;
   m_numCycle = 0;

   if (type() == ENTER)
   {
      if (rep() == COLUMN)
      {
         setPrimalBounds();
         setStatus(SPxBasis::PRIMAL);
      }
      else
      {
         setDualRowBounds();
         setStatus(SPxBasis::DUAL);
      }
      setEnterBounds();
      computeEnterCoPrhs();
   }
   else
   {
      if (rep() == ROW)
      {
         setPrimalBounds();
         setStatus(SPxBasis::PRIMAL);
      }
      else
      {
         setDualColBounds();
         setStatus(SPxBasis::DUAL);
      }
      setLeaveBounds();
      computeLeaveCoPrhs();
   }

   SPxBasis::coSolve(*theCoPvec, *theCoPrhs);
   computePvec();

   computeFrhs();
   SPxBasis::solve(*theFvec, *theFrhs);

   theShift = 0;
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

int SoPlex::sortLP(int pe, int nPes)
{
   int n = 0;
   int i;
   for (i = pe; i < thecovectors->num(); i += nPes)
   {
      n += (*thecovectors)[i].size();
      (const_cast<SVector*>(&((*thecovectors)[i])))->sort();
   }
   return n;
}

void SoPlex::setPricing(Pricing pr)
{
   thePricing = pr;
   if (initialized && type() == ENTER)
   {
      computePvec();
      computeCoTest();
      computeTest();
   }
}

void SoPlex::splitLP(int pe, int nPes)
{
   assert(pe   >= 0);
   assert(nPes > 0);

   int i, j, n;
   int start;
   int nnes = 0;
   int end = 0;
   int nVecs = subcovectors.size();

   for (n = pe; n < nVecs; n += nPes)
   {
      subcovectors[n].reSize(dim());
      start = nnes = 0;
      assert(nVecs > 0);
      for (j = (n * nNZEs) / nVecs; start < coDim(); ++start)
      {
         if (nnes >= j)
            break;
         nnes += vector(start).size();
      }
      end = start;
      for (j = ((n + 1) * nNZEs) / nVecs; end < coDim(); ++end)
      {
         if (nnes >= j)
            break;
         nnes += vector(end).size();
      }
      for (i = 0; i < thecovectors->num(); i++)
      {
         const SVector& vec = (*thecovectors)[i];
         int first = -1;
         for (j = 0; j < vec.size(); ++j)
         {
            if (vec.index(j) >= start && first < 0)
               first = j;
            if (vec.index(j) >= end)
               break;
         }
         if (first < 0)
            first = j;
         //  subcovectors[n][i] = SubSVector(&vec, first, j - first);
         subcovectors[n][i].assign(&vec, first, j - first);
      }
   }

#ifdef DEBUG
   if (pe == 0)
   {
      for (i = 0; i < dim(); i++)
      {
         int sum = 0;
         for (j = 0; j < nVecs; ++j)
            sum += subcovectors[j][i].size();
         if (sum != (*thecovectors)[i].size())
            std::cerr << pe << ": " << sum << "<->"
            << (*thecovectors)[i].size() << std::endl;
      }
   }
#endif // DEBUG
}

void SoPlex::splitLP()
{
   subcovectors.reSize(coDim() / coVecDim + 1);
   if (subcovectors.size() > 1)
   {
      nNZEs = sortLP (0, 1);
      splitLP(0, 1);
   }
}


/*
    The following method resizes all vectors and arrays of |SoPlex|
    (excluding inherited vectors).
 */
void SoPlex::reDim()
{
   int newdim = (rep() == ROW) ? SPxLP::nCols() : SPxLP::nRows();

   if (dim() > 0 && coDim() > 0)
   {
      int tmp = coDim() / (20 * dim()) + 1;
      coVecDim = coDim() / tmp + 1;
   }

   if (newdim > unitVecs.size())
   {
      unitVecs.reSize (newdim);
      while (newdim-- > 0)
         unitVecs[newdim] = newdim;
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

bool SoPlex::readFile(const char* filename)
{
   std::ifstream file(filename);

   if (!file)
      return false;
   else
      return  read(file);
}

void SoPlex::dumpFile(const char* filename) const
{
   std::ofstream file(filename);

   if (file.good())
      file << *this;
}

void SoPlex::clear()
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
   if (thesimplifier)
      thesimplifier->unload();

   unInit();
   setStatus(NO_PROBLEM);
   SPxLP::clear();
}

void SoPlex::clearUpdateVecs(void)
{
   theFvec->clearUpdate();
   thePvec->clearUpdate();
   theCoPvec->clearUpdate();
   solveVector2 = 0;
   coSolveVector2 = 0;
}

/*
    When the basis matrix factorization is recomputed from scratch, we also
    recompute the vectors.
 */
void SoPlex::factorize()
{
   SPxBasis::factorize();

   if (SPxBasis::status() >= REGULAR)
   {
      // #undef       NDEBUG
#ifndef NDEBUG
      DVector ftmp(fVec());
      DVector ptmp(pVec());
      DVector ctmp(coPvec());
      testVecs();
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
         std::cerr << std::endl;
         std::cerr << "fVec:      ";
         std::cerr << ftmp.length() << std::endl;
         std::cerr << std::endl;
         ftmp = fVec();
         multBaseWith(ftmp);
         ftmp -= fRhs();
         if (ftmp.length() > delta())
            std::cerr << iteration() << ": fVec error = " << ftmp.length() << std::endl;
      }
      if (ctmp.length() > delta())
      {
         std::cerr << std::endl;
         std::cerr << "coPvec:    ";
         std::cerr << ctmp.length() << std::endl;
         std::cerr << std::endl;
         ctmp = coPvec();
         multWithBase(ctmp);
         ctmp -= coPrhs();
         if (ctmp.length() > delta())
            std::cerr << iteration() << ": coPvec error = " << ctmp.length() << std::endl;
      }
      if (ptmp.length() > delta())
      {
         std::cerr << std::endl;
         std::cerr << "pVec:      ";
         std::cerr << ptmp.length() << std::endl;
         std::cerr << std::endl;
      }
#endif  // NDEBUG

      if (type() == ENTER)
      {
         computeCoTest();
         /*
                     if(pricing() == FULL)
                     {
                         computePvec();
                         computeTest();
                     }
         */
      }
      else
      {
         computeFtest();
         //          computePvec();
      }

   }
}

double SoPlex::maxInfeas() const
{
   int i;
   double inf = 0;

   if (type() == ENTER)
   {
      for (i = dim() - 1; i >= 0; --i)
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
      for (i = dim() - 1; i >= 0; --i)
      {
         if ((*theCoPvec)[i] > (*theCoUbound)[i])
            inf = MAXIMUM(inf, (*theCoPvec)[i] - (*theCoUbound)[i]);
         if ((*theCoLbound)[i] > (*theCoPvec)[i])
            inf = MAXIMUM(inf, (*theCoLbound)[i] - (*theCoPvec)[i]);
      }
      for (i = coDim() - 1; i >= 0; --i)
      {
         if ((*thePvec)[i] > (*theUbound)[i])
            inf = MAXIMUM(inf, (*thePvec)[i] - (*theUbound)[i]);
         else if ((*thePvec)[i] < (*theLbound)[i])
            inf = MAXIMUM(inf, (*theLbound)[i] - (*thePvec)[i]);
      }
   }

   return inf;
}

double SoPlex::nonbasicValue() const
{

   int i;
   double val = 0;
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

double SoPlex::value() const
{
   double x;

   if (!isInitialized())      
      /**@todo patch suggests returning SPxLP::infinity instead of initializing */
      (const_cast<SoPlex*>(this))->init();

   if (rep() == ROW)
   {
      if (type() == LEAVE)
         x = SPxLP::spxSense() * (coPvec() * fRhs());
      else
         x = SPxLP::spxSense() * (nonbasicValue() + (coPvec() * fRhs()));
   }
   else
      x = SPxLP::spxSense() * (nonbasicValue() + fVec() * coPrhs());

   if (thesimplifier)
      return thesimplifier->value(x);
   return x;
}

void SoPlex::setDelta(double d)
{
   thedelta = d;
}

void SoPlex::setEpsilon(double eps)
{
   primVec.delta().epsilon = eps;
   dualVec.delta().epsilon = eps;
   addVec.delta().epsilon = eps;
}

SoPlex::SoPlex(Type p_type, Representation p_rep,
                SPxPricer* pric, SPxRatioTester* rt,
                SPxStarter* start, SPxSimplifier* simple)
   : theType (p_type)
   , thePricing(FULL)
   , maxIters (-1)
   , maxTime (-1)
   , maxValue(infinity)
   , theShift (0)
   , m_maxCycle(100)
   , m_numCycle(0)
   , initialized (false)
   , solveVector2 (0)
   , coSolveVector2(0)
   , cacheProductFactor(4.0)
   , unitVecs (0)
   , primVec (0, 1e-16)
   , dualVec (0, 1e-16)
   , addVec (0, 1e-16)
   , thepricer (pric)
   , theratiotester(rt)
   , thestarter (start)
   , thesimplifier (simple)
{
   setRep (p_rep);
   setDelta (1e-6);
   setEpsilon(1e-16);
   theLP = this;
   coVecDim = 400;
}

/* We forbid the copyconstructor and the assignment operator
 * untill we are sure the work and we have an idea what exactly
 * they should be used to. Why would somebody copy an SoPlex Object?
 */
#if 0 
SoPlex::SoPlex(const SoPlex& old)
   : SPxLP (old)
   , SPxBasis (old)
   , theType (old.theType)
   , thePricing (old.thePricing)
   , therep (old.therep)  /// ??? siehe unten
   , maxIters (old.maxIters)
   , maxTime (old.maxTime)
   , maxValue(old.maxValue)
   , theShift (old.theShift)
   , m_maxCycle (old.m_maxCycle)
   , m_numCycle (old.m_numCycle)
   , initialized (old.initialized)
   , solveVector2 (0)
   , coSolveVector2 (0)
   , cacheProductFactor(old.cacheProductFactor)
   , unitVecs (old.unitVecs)
   , primRhs (old.primRhs)
   , primVec (old.primVec)
   , dualRhs (old.dualRhs)
   , dualVec (old.dualVec)
   , addVec (old.addVec)
   , theURbound (old.theURbound)
   , theLRbound (old.theLRbound)
   , theUCbound (old.theUCbound)
   , theLCbound (old.theLCbound)
   , theUBbound (old.theUBbound)
   , theLBbound (old.theLBbound)
   , theCoTest (old.theCoTest)
   , theTest (old.theTest)
   , thepricer (old.thepricer)
   , theratiotester (old.theratiotester)
   , thestarter (old.thestarter)
   , thesimplifier (old.thesimplifier)
{
   setRep (old.rep());
   setDelta(old.thedelta);
   coVecDim = 400;
   theLP = this;
}

SoPlex& SoPlex::operator=(const SoPlex& old)
{
   *(static_cast<SPxLP*>(this)) = old;
   *(static_cast<SPxBasis*>(this)) = old;

   therep = old.therep;
   unitVecs = old.unitVecs;
   theCoTest = old.theCoTest;
   theTest = old.theTest;
   theType = old.theType;
   thePricing = old.thePricing;
   primRhs = old.primRhs;
   primVec = old.primVec;
   dualRhs = old.dualRhs;
   dualVec = old.dualVec;
   addVec = old.addVec;
   theURbound = old.theURbound;
   theLRbound = old.theLRbound;
   theUCbound = old.theUCbound;
   theLCbound = old.theLCbound;
   theUBbound = old.theUBbound;
   theLBbound = old.theLBbound;
   m_maxCycle = old.m_maxCycle;
   m_numCycle = old.m_numCycle;
   theShift = old.theShift;
   initialized = old.initialized;
   cacheProductFactor = old.cacheProductFactor;
   thepricer = old.thepricer;
   theratiotester = old.theratiotester;
   thestarter = old.thestarter;
   thesimplifier = old.thesimplifier;
   solveVector2 = 0;
   coSolveVector2 = 0;

   setRep (old.rep());
   setDelta(old.thedelta);

   theLP = this;
   return *this;
}
#endif // no copy constructor and assignment operator

#define Inconsistent(file, line) \
{                                                                             \
   std::cout << file << "(" << line << ") ERROR: Inconsistency SoPlex\n";     \
   return 0;                                                                  \
}

int SoPlex::isConsistent() const
{
   if (epsilon() < 0)
      Inconsistent(__FILE__,__LINE__)

   if (primVec.delta().epsilon != dualVec.delta().epsilon)
      Inconsistent(__FILE__,__LINE__)
   if (dualVec.delta().epsilon != addVec.delta().epsilon)
      Inconsistent(__FILE__,__LINE__)

   if (unitVecs.size() < ((rep() == ROW) ? SPxLP::nCols() : SPxLP::nRows()))
      Inconsistent(__FILE__,__LINE__)

   if (initialized)
   {
      if (theFrhs->dim() != dim())
         Inconsistent(__FILE__,__LINE__)
      if (theFvec->dim() != dim())
         Inconsistent(__FILE__,__LINE__)

      if (theCoPrhs->dim() != dim())
         Inconsistent(__FILE__,__LINE__)
      if (thePvec->dim() != coDim())
         Inconsistent(__FILE__,__LINE__)
      if (theCoPvec->dim() != dim())
         Inconsistent(__FILE__,__LINE__)

      if (theTest.dim() != coDim())
         Inconsistent(__FILE__,__LINE__)
      if (theCoTest.dim() != dim())
         Inconsistent(__FILE__,__LINE__)

      if (theURbound.dim() != SPxLP::nRows())
         Inconsistent(__FILE__,__LINE__)
      if (theLRbound.dim() != SPxLP::nRows())
         Inconsistent(__FILE__,__LINE__)
      if (theUCbound.dim() != SPxLP::nCols())
         Inconsistent(__FILE__,__LINE__)
      if (theLCbound.dim() != SPxLP::nCols())
         Inconsistent(__FILE__,__LINE__)
      if (theUBbound.dim() != dim())
         Inconsistent(__FILE__,__LINE__)
      if (theLBbound.dim() != dim())
         Inconsistent(__FILE__,__LINE__)
   }

   if (rep() == COLUMN)
   {
      if
      (
         thecovectors != reinterpret_cast<const SVSet*>(static_cast<const LPRowSet*>(this)) ||
         thevectors != reinterpret_cast<const SVSet*>(static_cast<const LPColSet*>(this)) ||
         theFrhs != &primRhs ||
         theFvec != &primVec ||
         theCoPrhs != &dualRhs ||
         theCoPvec != &dualVec ||
         thePvec != &addVec ||
         theRPvec != theCoPvec ||
         theCPvec != thePvec ||
         theUbound != &theUCbound ||
         theLbound != &theLCbound ||
         theCoUbound != &theURbound ||
         theCoLbound != &theLRbound
     )
         Inconsistent(__FILE__,__LINE__)
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
         theCoLbound != &theLCbound
     )
         Inconsistent(__FILE__,__LINE__)
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

int SoPlex::nofNZEs() const
{
   int n = 0;
   for (int i = nCols(); --i >= 0;)
      n += colVector(i).size();
   return n;
}

void SoPlex::setTermination(double p_time, int p_iteration, double p_value)
{
   maxTime  = p_time;
   maxIters = p_iteration;
   maxValue = p_value;
}

void SoPlex::getTermination(double* p_time, int* p_iteration, double* p_value)
const
{
   if (p_time != 0)
      *p_time = maxTime;
   if (p_iteration != 0)
      *p_iteration = maxIters;
   if (p_value != 0)
      *p_value = maxValue;
}

SoPlex::Status SoPlex::getBasis(signed char row[], signed char col[]) const
{
   const SPxBasis::Desc& d = desc();
   int i;
   if (col)
   {
      for (i = nCols() - 1; i >= 0; --i)
         switch (d.colStatus(i))
         {
         case SPxBasis::Desc::P_ON_LOWER:
            col[i] = ON_LOWER;
            break;
         case SPxBasis::Desc::P_ON_UPPER:
            col[i] = ON_UPPER;
            break;
         case SPxBasis::Desc::P_FIXED:
            col[i] = FIXED;
            break;
         case SPxBasis::Desc::P_FREE:
            col[i] = ZERO;
            break;
         case SPxBasis::Desc::D_ON_UPPER:
         case SPxBasis::Desc::D_ON_LOWER:
         case SPxBasis::Desc::D_ON_BOTH:
         case SPxBasis::Desc::D_UNDEFINED:
         case SPxBasis::Desc::D_FREE:
            col[i] = BASIC;
            break;
         default:
            std::cout << std::endl << "ERROR: "
            << __FILE__ << ":" << __LINE__ << std::endl << std::endl;
            break;
         }
   }

   if (row)
   {
      for (i = nRows() - 1; i >= 0; --i)
         switch (d.rowStatus(i))
         {
         case SPxBasis::Desc::P_ON_LOWER:
            row[i] = ON_LOWER;
            break;
         case SPxBasis::Desc::P_ON_UPPER:
            row[i] = ON_UPPER;
            break;
         case SPxBasis::Desc::P_FIXED:
            row[i] = FIXED;
            break;
         case SPxBasis::Desc::P_FREE:
            row[i] = ZERO;
            break;
         case SPxBasis::Desc::D_ON_UPPER:
         case SPxBasis::Desc::D_ON_LOWER:
         case SPxBasis::Desc::D_ON_BOTH:
         case SPxBasis::Desc::D_UNDEFINED:
         case SPxBasis::Desc::D_FREE:
            row[i] = BASIC;
            break;
         default:
            std::cout << std::endl << "ERROR: "
            << __FILE__ << ":" << __LINE__ << std::endl << std::endl;
            break;
         }
   }
   return status();
}

void SoPlex::setBasis(const signed char p_rows[], const signed char p_cols[])
{
   if (SPxBasis::status() == SPxBasis::NO_PROBLEM)
      SPxBasis::load(this);

   SPxBasis::Desc ds = desc();
   int i;

   for (i = nRows() - 1; i >= 0; i--)
   {
      switch (p_rows[i])
      {
      case FIXED :
         assert(rhs(i) == lhs(i));
         ds.rowStatus(i) = SPxBasis::Desc::P_FIXED;
         break;
      case ON_UPPER :
         assert(rhs(i) < SPxLP::infinity);
         ds.rowStatus(i) = SPxBasis::Desc::P_ON_UPPER;
         break;
      case ON_LOWER :
         assert(lhs(i) > -SPxLP::infinity);
         ds.rowStatus(i) = SPxBasis::Desc::P_ON_LOWER;
         break;
      case ZERO :
         assert(lhs(i) <= -SPxLP::infinity && rhs(i) >= SPxLP::infinity);
         ds.rowStatus(i) = SPxBasis::Desc::P_FREE;
         break;
      case BASIC :
         ds.rowStatus(i) = dualRowStatus(i);
         break;
      default:
         std::cout << std::endl << "ERROR: " << __FILE__ << ":" << __LINE__
         << std::endl << std::endl;
      }
   }

   for (i = nCols() - 1; i >= 0; i--)
   {
      switch (p_cols[i])
      {
      case FIXED :
         assert(upper(i) == lower(i));
         ds.colStatus(i) = SPxBasis::Desc::P_FIXED;
         break;
      case ON_UPPER :
         assert(upper(i) < SPxLP::infinity);
         ds.colStatus(i) = SPxBasis::Desc::P_ON_UPPER;
         break;
      case ON_LOWER :
         assert(lower(i) > -SPxLP::infinity);
         ds.colStatus(i) = SPxBasis::Desc::P_ON_LOWER;
         break;
      case ZERO :
         assert(lower(i) <= -SPxLP::infinity && upper(i) >= SPxLP::infinity);
         ds.colStatus(i) = SPxBasis::Desc::P_FREE;
         break;
      case BASIC :
         ds.colStatus(i) = dualColStatus(i);
         break;
      default:
         std::cout << std::endl << "ERROR: " << __FILE__ << ":" << __LINE__
         << std::endl << std::endl;
      }
   }
   loadBasis(ds);
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
