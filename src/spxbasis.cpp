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
#pragma ident "@(#) $Id: spxbasis.cpp,v 1.17 2002/01/19 18:59:16 bzfkocht Exp $"

#include <assert.h>
#include <iostream>
#include <math.h>

#include "real.h"
#include "spxbasis.h"
#include "didxset.h"
#include "dvector.h"
#include "soplex.h"
#include "mpsinput.h"
#include "message.h"

namespace soplex
{
static Real minStab;
#define EPS     minStab
//#define       EPS     1e-6

SPxBasis::Desc::Status
SPxBasis::dualStatus(const SPxLP::SPxColId& id) const
{
   return dualColStatus(static_cast<SPxLP*>(theLP)->number(id));
}

SPxBasis::Desc::Status
SPxBasis::dualStatus(const SPxLP::SPxRowId& id) const
{
   return dualRowStatus((static_cast<SPxLP*>(theLP))->number(id));
}

SPxBasis::Desc::Status
SPxBasis::dualRowStatus(int i) const
{
   assert(theLP != 0);

   if (theLP->rhs(i) < SPxLP::infinity)
   {
      if (theLP->lhs(i) > -SPxLP::infinity)
      {
         if (theLP->lhs(i) == theLP->rhs(i))
            return Desc::D_FREE;
         else
            return Desc::D_ON_BOTH;
      }
      else
         return Desc::D_ON_LOWER;
   }
   else if (theLP->lhs(i) > -SPxLP::infinity)
      return Desc::D_ON_UPPER;
   else
      return Desc::D_UNDEFINED;
}

SPxBasis::Desc::Status
SPxBasis::dualColStatus(int i) const
{
   assert(theLP != 0);

   if (theLP->SPxLP::upper(i) < SPxLP::infinity)
   {
      if (theLP->SPxLP::lower(i) > -SPxLP::infinity)
      {
         if (theLP->SPxLP::lower(i) == theLP->SPxLP::upper(i))
            return Desc::D_FREE;
         else
            return Desc::D_ON_BOTH;
      }
      else
         return Desc::D_ON_LOWER;
   }
   else if (theLP->SPxLP::lower(i) > -SPxLP::infinity)
      return Desc::D_ON_UPPER;
   else
      return Desc::D_UNDEFINED;
}

void SPxBasis::loadMatrixVecs()
{
   assert(theLP != 0);

   int i;
   nzCount = 0;
   for (i = theLP->dim() - 1; i >= 0; --i)
   {
      matrix[i] = &theLP->vector(baseId(i));
      nzCount += matrix[i]->size();
   }
   matrixIsSetup = true;
   factorized = false;
}

/*
    Loading a #Desc# into the basis can be done more efficiently, by
    explicitely programming both cases, for the rowwise and for the columnwise
    representation. This implementation hides this distingtion in the use of
    methods #isBasic()# and #vector()#.
 */
void SPxBasis::load(const Desc& ds)
{
   assert(status() > NO_PROBLEM);
   assert(theLP != 0);
   assert(ds.nRows() == theLP->nRows());
   assert(ds.nCols() == theLP->nCols());

   SPxLP::Id none;
   int i, j;

   lastin = none;
   lastout = none;
   lastidx = -1;
   iterCount = 0;
   updateCount = 0;

   if (&ds != &thedesc)
   {
      thedesc = ds;
      setRep();
   }

   nzCount = 0;
   for (j = i = 0; i < theLP->nRows(); ++i)
   {
      if (theLP->isBasic(thedesc.rowStatus(i)))
      {
         SPxLP::SPxRowId id = theLP->SPxLP::rId(i);
         theBaseId[j] = id;
         matrix[j] = &theLP->vector(id);
         nzCount += matrix[j++]->size();
         /*
         std::cerr << "\tR" << theLP->number(id);
         if(j % 8 == 0)
             std::cerr << std::endl;
          */
      }
   }

   for (i = 0; i < theLP->nCols(); ++i)
   {
      if (theLP->isBasic(thedesc.colStatus(i)))
      {
         SPxLP::SPxColId id = theLP->SPxLP::cId(i);
         theBaseId[j] = id;
         matrix[j] = &theLP->vector(id);
         nzCount += matrix[j++]->size();
         /*
         std::cerr << "\tC" << theLP->number(id);
         if(j % 8 == 0)
             std::cerr << std::endl;
          */
      }
   }

   assert(j == matrix.size());

   matrixIsSetup = true;
   factorized = false;
   if (factor != 0)
      factor->clear();
}

void SPxBasis::setRep()
{
   assert(theLP != 0);

   reDim();
   minStab = 1e-4;

   if (theLP->rep() == SoPlex::ROW)
   {
      thedesc.stat = & thedesc.rowstat;
      thedesc.costat = & thedesc.colstat;
   }
   else
   {
      thedesc.stat = & thedesc.colstat;
      thedesc.costat = & thedesc.rowstat;
   }
}

void SPxBasis::load(SoPlex* lp)
{
   assert(lp != 0);
   theLP = lp;

   setRep();

   addedRows(lp->nRows());
   addedCols(lp->nCols());

   setStatus(REGULAR);

   load(thedesc);
}

void SPxBasis::load(SLinSolver* p_solver)
{
   factor = p_solver;
   factorized = false;
   factor->clear();
}

/**@todo This routine is untested.
 */
void SPxBasis::readBasis(std::istream& is, NameSet& rn, NameSet& cn)
{
   assert(theLP != 0);

   int  i;
   Desc l_desc(thedesc);

   for(i = 0; i < theLP->nRows(); i++)
      l_desc.rowstat[i] = dualRowStatus(i);

   for(i = 0; i < theLP->nCols(); i++)
      l_desc.colstat[i] = Desc::P_ON_LOWER;

   MPSInput mps(is);

   if (mps.readLine() && (mps.field0() != 0) && !strcmp(mps.field0(), "NAME"))
   {
      while(mps.readLine())
      {
         int c = -1;
         int r = -1;

         if (!strcmp(mps.field0(), "ENDATA"))
         {
            mps.setSection(MPSInput::ENDATA);
            break;
         }
         if ((mps.field1() == 0) || (mps.field2() == 0))
            break;

         if ((c = cn.number(mps.field2())) < 0)
            break;

         if (*mps.field1() == 'X')
            if ((mps.field3() == 0) || ((r = rn.number(mps.field3())) < 0))
               break;

         if (!strcmp(mps.field1(), "XU"))
         {
            l_desc.colstat[c] = dualColStatus(c);
            l_desc.rowstat[r] = Desc::P_ON_UPPER;
         }
         else if (!strcmp(mps.field1(), "XL"))
         {
            l_desc.colstat[c] = dualColStatus(c);
            l_desc.rowstat[r] = Desc::P_ON_LOWER;
         }
         else if (!strcmp(mps.field1(), "UL"))
         {
            l_desc.colstat[c] = Desc::P_ON_UPPER;
         }
         else if (!strcmp(mps.field1(), "LL"))
         {
            l_desc.colstat[c] = Desc::P_ON_LOWER;
         }
         else
         {
            break;
         }
      }
   }
   if (!mps.hasError())
   {
      if (mps.section() == MPSInput::ENDATA)
         load(l_desc);
      else
         mps.syntaxError();
   }
}

/*      \SubSection{Pivoting Methods}
 */
int SPxBasis::doFactorize()
{
   if (nonzeroFactor < 0)
      return (updateCount >= -nonzeroFactor);

   Real newFac = nzFac + factor->memory();
   Real neu = (newFac + lastFill * nzCount) / (updateCount + 1);
   Real alt = (nzFac + lastFill * nzCount) / updateCount;

   return (updateCount >= maxUpdates || neu > alt);
}

void SPxBasis::change
(
   int i,
   SPxLP::Id& id,
   const SVector* enterVec,
   const SSVector* eta
)
{
   assert(!id.isValid() || (enterVec != 0));

   assert(factor != 0);
   lastidx = i;
   lastin = id;


   if (id.isValid() && i >= 0)
   {
      /*
      if(iterCount < 80)
   {
          std::cerr << i << ":\t";
          if(theBaseId[i].isSPxRowId())
              std::cerr << 'R';
          else
              std::cerr << 'C';
          std::cerr << theLP->number(theBaseId[i]) << "\t->  ";
          if(id.isSPxRowId())
              std::cerr << 'R';
          else
              std::cerr << 'C';
          std::cerr << theLP->number(id) << "\n";
   }
       */
      assert(enterVec != 0);

      nzCount = nzCount - matrix[i]->size() + enterVec->size();
      matrix[i] = enterVec;
      lastout = theBaseId[i];
      theBaseId[i] = id;

      ++iterCount;
      ++updateCount;
      Real newFac = nzFac + factor->memory();
      if (doFactorize())
         factorize();
      else
      {
         // Real    s = factor->stability();
         factor->change(i, *enterVec, eta);
         if (factor->status() != SLinSolver::OK
              || factor->stability() < EPS)
         {
            // std::cerr << s << " -> " << factor->stability() << '\t';
            factorize();
         }
         else
            nzFac = newFac;
      }
   }
   else
      lastout = id;
}

void SPxBasis::factorize()
{
   assert(factor != 0);

   if (!matrixIsSetup)
      load(thedesc);

   updateCount = 0;
   factorized = true;
   switch (factor->load(matrix.get_ptr(), matrix.size()))
   {
   case SLinSolver::OK :
      if (status() == SINGULAR)
         setStatus(REGULAR);
      minStab = factor->stability();
      if (minStab > 1e-4)
         minStab *= 0.001;
      if (minStab > 1e-5)
         minStab *= 0.01;
      if (minStab > 1e-6)
         minStab *= 0.1;
      break;
   case SLinSolver::SINGULAR :
      setStatus(SINGULAR);
      break;
   default :
      std::cerr << "ERROR: unknown status of factorization.\n";
      abort();
      // factorized = false;
   }
   lastFill = Real(factor->memory()) * nonzeroFactor / Real(nzCount);
   nzFac = 0;
}

Vector& SPxBasis::multWithBase(Vector& x) const
{
   assert(status() > SINGULAR);
   assert(theLP->dim() == x.dim());

   int i;
   DVector tmp(x);

   if (!matrixIsSetup)
      (const_cast<SPxBasis*>(this))->load(thedesc);

   for (i = x.dim() - 1; i >= 0; --i)
      x[i] = *(matrix[i]) * tmp;

   return x;
}

Vector& SPxBasis::multBaseWith(Vector& x) const
{
   assert(status() > SINGULAR);
   assert(theLP->dim() == x.dim());

   int i;
   DVector tmp(x);

   if (!matrixIsSetup)
      (const_cast<SPxBasis*>(this))->load(thedesc);

   x.clear();
   for (i = x.dim() - 1; i >= 0; --i)
   {
      if (tmp[i])
         x.multAdd(tmp[i], *(matrix[i]));
   }

   return x;
}

bool SPxBasis::isConsistent() const
{
   int primals = 0;
   int i;

   if (status() > NO_PROBLEM)
   {
      if (theLP == 0)
         return MSGinconsistent("SPxBasis");

      if (theBaseId.size() != theLP->dim() || matrix.size() != theLP->dim())
         return MSGinconsistent("SPxBasis");

      if (thedesc.nCols() != theLP->nCols() 
         || thedesc.nRows() != theLP->nRows())
         return MSGinconsistent("SPxBasis");

      for (i = thedesc.nRows() - 1; i >= 0; --i)
      {
         if (thedesc.rowStatus(i) >= 0)
         {
            if (thedesc.rowStatus(i) != dualRowStatus(i))
               return MSGinconsistent("SPxBasis");
         }
         else
            ++primals;
      }
      
      for (i = thedesc.nCols() - 1; i >= 0; --i)
      {
         if (thedesc.colStatus(i) >= 0)
         {
            if (thedesc.colStatus(i) != dualColStatus(i))
               return MSGinconsistent("SPxBasis");
         }
         else
            ++primals;
      }
      if (primals != thedesc.nCols())
         return MSGinconsistent("SPxBasis");
   }
   return thedesc.isConsistent()
          && theBaseId.isConsistent()
          && matrix.isConsistent()
          && factor->isConsistent();
}

SPxBasis::SPxBasis()
   : theLP (0)
   , matrixIsSetup (false)
   , factor (0)
   , maxUpdates (1000)
   , nonzeroFactor (10)
   , nzCount (1)
   , thestatus (NO_PROBLEM)
{}
} // namespace soplex

//-----------------------------------------------------------------------------
//Emacs Local Variables:
//Emacs mode:c++
//Emacs c-basic-offset:3
//Emacs tab-width:8
//Emacs indent-tabs-mode:nil
//Emacs End:
//-----------------------------------------------------------------------------
