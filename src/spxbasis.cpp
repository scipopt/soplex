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
#pragma ident "@(#) $Id: spxbasis.cpp,v 1.9 2001/12/10 15:46:49 bzfkocht Exp $"



//@ -----------------------------------------------------------------------------
/*      \Section{Complex Method Implementations}
 */
#include <assert.h>
#include <iostream>
#include <math.h>

#include "spxbasis.h"



#include "didxset.h"
#include "dvector.h"
#include "soplex.h"

namespace soplex
{




static double minStab;
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
   if (factor)
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

/**@todo Reimplement with new MPS Reader.
 */
void SPxBasis::readBasis(std::istream& is, NameSet& rn, NameSet& cn)
{
   std::cout << "Not implemented yet!" << std::endl;
#if 0
   int i, r, c;
   char* f1;
   char* f2;
   char* f3;
   char* f4;
   char* f5;
   char* f6;
   Desc l_desc(thedesc);

   assert(theLP != 0);

   for (i = theLP->nRows() - 1; i >= 0; --i)
      l_desc.rowstat[i] = dualRowStatus(i);

   for (i = theLP->nCols() - 1; i >= 0; --i)
      l_desc.colstat[i] = Desc::P_ON_LOWER;

   if (!SPxLP::readLine(is, f1, f2, f3, f4, f5, f6))
   {
      std::cerr << "ERROR: incorrect basis file format\n";
      return;
   }
   if (strncmp(f1, "NAME", 4) != 0)
   {
      std::cerr << "ERROR: incorrect basis file format\n";
      return;
   }

   do
   {
      if (SPxLP::readLine(is, f1, f2, f3, f4, f5, f6))
      {
         if (strncmp(f1, "ENDATA", 6) == 0)
            break;
         std::cerr << "ERROR: incorrect basis file format\n";
         return;
      }

      if (f2 == 0)
      {
         std::cerr << "ERROR: incorrect basis file format\n";
         return;
      }

      c = cn.number(f2);
      if (c < 0)
      {
         std::cerr << "ERROR: incorrect basis file format\n";
         return;
      }

      if (strncmp(f1, "XU", 2))
      {
         r = rn.number(f3);
         if (r < 0)
         {
            std::cerr << "ERROR: incorrect basis file format\n";
            return;
         }
         l_desc.colstat[c] = dualColStatus(c);
         l_desc.rowstat[r] = Desc::P_ON_UPPER;
      }
      else if (strncmp(f1, "XL", 2))
      {
         r = rn.number(f3);
         if (r < 0)
         {
            std::cerr << "ERROR: incorrect basis file format\n";
            return;
         }
         l_desc.colstat[c] = dualColStatus(c);
         l_desc.rowstat[r] = Desc::P_ON_LOWER;
      }
      else if (strncmp(f1, "UL", 2))
      {
         l_desc.colstat[c] = Desc::P_ON_UPPER;
      }
      else if (strncmp(f1, "LL", 2))
      {
         l_desc.colstat[c] = Desc::P_ON_LOWER;
      }
      else
      {
         std::cerr << "ERROR: incorrect basis file format\n";
         return;
      }
      //      std::cerr << f1 << '\t' << f2 << '\t' << f3 << std::endl;
   }

   while (is.good());

   load(l_desc);
   return;
#endif
}

//@ -----------------------------------------------------------------------------
/*      \SubSection{Pivoting Methods}
 */
int SPxBasis::doFactorize()
{
   if (nonzeroFactor < 0)
      return (updateCount >= -nonzeroFactor);

   double newFac = nzFac + factor->memory();
   double neu = (newFac + lastFill * nzCount) / (updateCount + 1);
   double alt = (nzFac + lastFill * nzCount) / updateCount;

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
      nzCount = nzCount - matrix[i]->size() + enterVec->size();
      matrix[i] = enterVec;
      lastout = theBaseId[i];
      theBaseId[i] = id;

      ++iterCount;
      ++updateCount;
      double newFac = nzFac + factor->memory();
      if (doFactorize())
         factorize();
      else
      {
         // double    s = factor->stability();
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
   lastFill = double(factor->memory()) * nonzeroFactor / double(nzCount);
   nzFac = 0;
}


//@ -----------------------------------------------------------------------------
/*      \SubSection{Linear Algebra}
 */
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

//@ -----------------------------------------------------------------------------
/*      \SubSubSection{Consistency Check}
 */
#define inconsistent                                            \
{                                                               \
std::cout << "Inconsistency detected in class SPxBasis\n";      \
return 0;                                                  \
}

int SPxBasis::isConsistent() const
{
   int primals = 0;
   int i;

   if (status() > NO_PROBLEM)
   {
      if (theLP == 0)
         inconsistent

         if
         (
            theBaseId.size() != theLP->dim() ||
            matrix.size() != theLP->dim()
        )
            inconsistent

            if
            (
               thedesc.nCols() != theLP->nCols() ||
               thedesc.nRows() != theLP->nRows()
           )
               inconsistent

               for (i = thedesc.nRows() - 1; i >= 0; --i)
               {
                  if (thedesc.rowStatus(i) >= 0)
                  {
                     if (thedesc.rowStatus(i) != dualRowStatus(i))
                        inconsistent;
                  }
                  else
                     ++primals;
               }

      for (i = thedesc.nCols() - 1; i >= 0; --i)
      {
         if (thedesc.colStatus(i) >= 0)
         {
            if (thedesc.colStatus(i) != dualColStatus(i))
               inconsistent;
         }
         else
            ++primals;
      }

      if (primals != thedesc.nCols())
         inconsistent
      }

   return thedesc.isConsistent()
          && theBaseId.isConsistent()
          && matrix.isConsistent()
          && factor->isConsistent();
}

SPxBasis::SPxBasis()
   : theLP (0)
   , thestatus (NO_PROBLEM)
   , matrixIsSetup (false)
   , factor (0)
   , maxUpdates (1000)
   , nonzeroFactor (10)
   , nzCount (1)
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
