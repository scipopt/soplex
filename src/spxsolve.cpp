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
#pragma ident "@(#) $Id: spxsolve.cpp,v 1.7 2001/11/25 14:58:29 bzfkocht Exp $"

#include <assert.h>
#include <iostream>

#include "soplex.h"
#include "spxpricer.h"
#include "spxratiotester.h"
#include "spxstarter.h"
#include "spxsimplifier.h"

namespace soplex
{
LPSolver::Status SoPlex::solve()
{
   Id enterId;
   int leaveNum;

   if (dim() <= 0 && coDim() <= 0)          // no problem loaded
      return LPSolver::ERROR;

   if (slinSolver() == 0)             // linear system solver is required.
      return LPSolver::ERROR;

   if (thesimplifier)
   {
      // if (thesimplifier->loadedLP() != this)
      thesimplifier->load(this);

      switch (thesimplifier->simplify())
      {
      case 1:
         setStatus(SPxBasis::UNBOUNDED);
         return LPSolver::UNBOUNDED;
      case - 1:
         setStatus(SPxBasis::INFEASIBLE);
         return LPSolver::INFEASIBLE;
      default:
         break;
      }
   }
   if (thepricer == 0)                                // pricer is required.
      return LPSolver::ERROR;

   theTime.reset();
   theTime.start();

   m_numCycle = 0;

   splitLP();

   if (!isInitialized())
   {
      /*
      if(SPxBasis::status() <= NO_PROBLEM)
          SPxBasis::load(this);
       */
      if (thestarter)                            // no basis and no starter.
         thestarter->generate(*this);              // generate start basis.
      init();
   }
   thepricer->setEpsilon(delta());
   setType(type());

#ifdef  DEBUG
   std::cerr << "starting value = " << value() << '\n';
   std::cerr << "starting shift = " << shift() << '\n';
   {
      int i;

      std::cout << "column status':\t";
      for (i = 0; i < desc().nCols(); ++i)
      {
         switch (desc().colStatus(i))
         {
         case SPxBasis::Desc::P_ON_LOWER + SPxBasis::Desc::P_ON_UPPER:
            std::cout << "X ";
            break;
         case SPxBasis::Desc::P_ON_LOWER:
            std::cout << "L ";
            break;
         case SPxBasis::Desc::P_ON_UPPER:
            std::cout << "U ";
            break;
         case SPxBasis::Desc::P_FREE:
            std::cout << "F ";
            break;
         case SPxBasis::Desc::D_ON_LOWER + SPxBasis::Desc::D_ON_UPPER:
            std::cout << "x ";
            break;
         case SPxBasis::Desc::D_ON_LOWER:
            std::cout << "l ";
            break;
         case SPxBasis::Desc::D_ON_UPPER:
            std::cout << "u ";
            break;
         case SPxBasis::Desc::D_FREE:
            std::cout << "f ";
            break;
         default:
            std::cout << ". ";
            break;
         }
      }
      std::cout << "\n";

      std::cout << "row status':\t";
      for (i = 0; i < desc().nRows(); ++i)
      {
         switch (desc().rowStatus(i))
         {
         case SPxBasis::Desc::P_ON_LOWER + SPxBasis::Desc::P_ON_UPPER:
            std::cout << "X ";
            break;
         case SPxBasis::Desc::P_ON_LOWER:
            std::cout << "L ";
            break;
         case SPxBasis::Desc::P_ON_UPPER:
            std::cout << "U ";
            break;
         case SPxBasis::Desc::P_FREE:
            std::cout << "F ";
            break;
         case SPxBasis::Desc::D_ON_LOWER + SPxBasis::Desc::D_ON_UPPER:
            std::cout << "x ";
            break;
         case SPxBasis::Desc::D_ON_LOWER:
            std::cout << "l ";
            break;
         case SPxBasis::Desc::D_ON_UPPER:
            std::cout << "u ";
            break;
         case SPxBasis::Desc::D_FREE:
            std::cout << "f ";
            break;
         default:
            std::cout << ". ";
            break;
         }
      }
      std::cout << "\n" << flush;
   }
#endif  // DEBUG

   if (SPxBasis::status() == SPxBasis::OPTIMAL)
      setStatus(SPxBasis::REGULAR);
   int stop = terminate();
   leaveCount = enterCount = 0;
   while (!stop)
   {
      if (type() == ENTER)
      {
#ifndef NDEBUG
         std::cout << "Enter iteration: " << iteration()
         << "\tValue = " << value()
         << "\tShift = " << shift() << std::endl;
#endif  // NDEBUG
         do
         {
            enterId = thepricer->selectEnter();
            if (!enterId.isValid())
               break;
            enter(enterId);
            assert((testBounds(), 1));
            thepricer->entered4(lastEntered(), lastIndex());
            stop = terminate();
            clearUpdateVecs();
            enterCount += (lastIndex() >= 0);
            //@ assert(isConsistent());
         }

         while (!stop);
         if (!stop)
         {
            if (shift() <= epsilon())
            {
               factorize();
               unShift();
               if (maxInfeas() + shift() <= delta())
               {
                  setStatus(SPxBasis::OPTIMAL);
                  break;
               }
            }
            setType(LEAVE);
         }
      }

      else
      {
         assert(type() == LEAVE);
#ifndef NDEBUG
         std::cout << "Leave Iteration: " << iteration()
         << "\tValue = " << value()
         << "\tShift = " << shift() << std::endl;
#endif  // NDEBUG
         do
         {
            leaveNum = thepricer->selectLeave();
            if (leaveNum < 0)
               break;
            leave(leaveNum);
            assert((testBounds(), 1));
            thepricer->left4(lastIndex(), lastLeft());
            stop = terminate();
            clearUpdateVecs();
            leaveCount += lastEntered().isValid();
            //@ assert(isConsistent());
         }

         while (!stop);
         if (!stop)
         {
            if (shift() <= epsilon())
            {
               factorize();
               unShift();
               if (maxInfeas() + shift() <= delta())
               {
                  setStatus(SPxBasis::OPTIMAL);
                  break;
               }
            }
            setType(ENTER);
         }
      }
   }

   theTime.stop();
   return status();
}

void SoPlex::testVecs()
{
   int i;
   DVector tmp(dim());

   tmp = *theCoPvec;
   multWithBase(tmp);
   tmp -= *theCoPrhs;
   if (tmp.length() > delta())
   {
      std::cerr << iteration() << ":\tcoP error = " << tmp.length();
      tmp.clear();
      SPxBasis::coSolve(tmp, *theCoPrhs);
      multWithBase(tmp);
      tmp -= *theCoPrhs;
      std::cerr << "\t[" << tmp.length() << "]\t(";
      tmp.clear();
      SPxBasis::coSolve(tmp, *theCoPrhs);
      tmp -= *theCoPvec;
      std::cerr << tmp.length() << ")\n";
   }

   tmp = *theFvec;
   multBaseWith(tmp);
   tmp -= *theFrhs;
   if (tmp.length() > delta())
   {
      std::cerr << iteration() << ":\t  F error = " << tmp.length() << "\t(";
      tmp.clear();
      SPxBasis::solve(tmp, *theFrhs);
      tmp -= *theFvec;
      std::cerr << tmp.length() << ")\n";
   }

   if (type() == ENTER)
   {
      for (i = 0; i < dim(); ++i)
      {
         if (theCoTest[i] < 0 && isCoBasic(i))
            std::cerr << "this shalt not be!\n";
      }
      for (i = 0; i < coDim(); ++i)
      {
         if (theTest[i] < 0 && isBasic(i))
            std::cerr << "this shalt not be!\n";
      }
   }
}

int SoPlex::terminate()
{
#ifndef NDEBUG
   testVecs();
#endif  // NDEBUG

   int redo = dim();
   if (redo < 1000)
      redo = 1000;

   if (iteration() > 10 && iteration() % redo == 0)
   {
#ifndef NDEBUG
      DVector cr(*theCoPrhs);
      DVector fr(*theFrhs);
#endif  // !NDEBUG

      if (type() == ENTER)
         computeEnterCoPrhs();
      else
         computeLeaveCoPrhs();
      computeFrhs();

#ifndef NDEBUG
      cr -= *theCoPrhs;
      fr -= *theFrhs;
      if (cr.length() > delta())
         std::cerr << "unexpected change of coPrhs " << cr.length() << std::endl;
      if (fr.length() > delta())
         std::cerr << "unexpected change of   Frhs " << fr.length() << std::endl;
#endif  // !NDEBUG

      if (updateCount > 1)
         factorize();
      SPxBasis::coSolve(*theCoPvec, *theCoPrhs);
      SPxBasis::solve (*theFvec, *theFrhs);

      if (pricing() == FULL)
      {
         computePvec();
         if (type() == ENTER)
            computeTest();
      }
      if (shift() > 0)
         unShift();
   }

   if (maxIters >= 0 && iterations() >= maxIters)
      return 1;
   if (maxTime >= 0 && time() >= maxTime)
      return 1;

   return (SPxBasis::status() >= SPxBasis::OPTIMAL
            || SPxBasis::status() <= SPxBasis::SINGULAR);
}


//@ ----------------------------------------------------------------------------
/*      \SubSection{Accessing Results}
 */
LPSolver::Status SoPlex::getPrimal (Vector& p_vector) const
{
   if (!isInitialized())
      const_cast<SoPlex*>(this)->init();

   if (rep() == ROW)
      p_vector = coPvec();

   else
   {
      int i;
      const SPxBasis::Desc& ds = desc();
      for (i = nCols() - 1; i >= 0; --i)
      {
         switch (ds.colStatus(i))
         {
         case SPxBasis::Desc::P_ON_LOWER :
            p_vector[i] = SPxLP::lower(i);
            break;
         case SPxBasis::Desc::P_ON_UPPER :
         case SPxBasis::Desc::P_FIXED :
            p_vector[i] = SPxLP::upper(i);
            break;
         case SPxBasis::Desc::P_FREE :
            p_vector[i] = 0;
            break;
         case SPxBasis::Desc::D_FREE :
         case SPxBasis::Desc::D_ON_UPPER :
         case SPxBasis::Desc::D_ON_LOWER :
         case SPxBasis::Desc::D_ON_BOTH :
         case SPxBasis::Desc::D_UNDEFINED :
            break;
         default:
            abort();
         }
      }
      for (i = dim() - 1; i >= 0; --i)
      {
         if (baseId(i).isSPxColId())
            p_vector[ number(SPxColId(baseId(i))) ] = fVec()[i];
      }
   }

   return status();
}

LPSolver::Status SoPlex::getDual (Vector& p_vector) const
{
   if (!isInitialized())
      const_cast<SoPlex*>(this)->init();

   if (rep() == ROW)
   {
      int i;
      p_vector.clear ();
      for (i = nCols() - 1; i >= 0; --i)
      {
         if (baseId(i).isSPxRowId())
            p_vector[ number(SPxRowId(baseId(i))) ] = fVec()[i];
      }
   }

   else
      p_vector = coPvec();

   p_vector *= spxSense();

   return status();
}

LPSolver::Status SoPlex::getRdCost (Vector& p_vector) const
{
   if (!isInitialized())
      const_cast<SoPlex*>(this)->init();

   if (rep() == ROW)
   {
      int i;
      p_vector.clear();
      if (spxSense() == SPxLP::MINIMIZE)
      {
         for (i = dim() - 1; i >= 0; --i)
         {
            if (baseId(i).isSPxColId())
               p_vector[ number(SPxColId(baseId(i))) ] = -fVec()[i];
         }
      }
      else
      {
         for (i = dim() - 1; i >= 0; --i)
         {
            if (baseId(i).isSPxColId())
               p_vector[ number(SPxColId(baseId(i))) ] = fVec()[i];
         }
      }
   }

   else
   {
      p_vector = maxObj();
      p_vector -= pVec();
      if (spxSense() == SPxLP::MINIMIZE)
         p_vector *= -1;
   }

   return status();
}

LPSolver::Status SoPlex::getSlacks (Vector& p_vector) const
{
   if (!isInitialized())
      const_cast<SoPlex*>(this)->init();

   if (rep() == COLUMN)
   {
      int i;
      const SPxBasis::Desc& ds = desc();
      for (i = nRows() - 1; i >= 0; --i)
      {
         switch (ds.rowStatus(i))
         {
         case SPxBasis::Desc::P_ON_LOWER :
            p_vector[i] = lhs(i);
            break;
         case SPxBasis::Desc::P_ON_UPPER :
         case SPxBasis::Desc::P_FIXED :
            p_vector[i] = rhs(i);
            break;
         case SPxBasis::Desc::P_FREE :
            p_vector[i] = 0;
            break;
         case SPxBasis::Desc::D_FREE :
         case SPxBasis::Desc::D_ON_UPPER :
         case SPxBasis::Desc::D_ON_LOWER :
         case SPxBasis::Desc::D_ON_BOTH :
         case SPxBasis::Desc::D_UNDEFINED :
            break;
         default:
            abort();
         }
      }
      for (i = dim() - 1; i >= 0; --i)
      {
         if (baseId(i).isSPxRowId())
            p_vector[ number(SPxRowId(baseId(i))) ] = -(*theFvec)[i];
      }
   }

   else
      p_vector = pVec();

   return status();
}

LPSolver::Status SoPlex::status() const
{
   switch (SPxBasis::status())
   {
   case SPxBasis::NO_PROBLEM :
      return LPSolver::UNKNOWN;
   case SPxBasis::SINGULAR :
      return LPSolver::ERROR;
   case SPxBasis::REGULAR :
      return LPSolver::UNKNOWN;
   case SPxBasis::DUAL :
      return LPSolver::DUAL;
   case SPxBasis::PRIMAL :
      return LPSolver::PRIMAL;
   case SPxBasis::OPTIMAL :
      return LPSolver::SOLVED;
   case SPxBasis::UNBOUNDED :
      return LPSolver::UNBOUNDED;
   case SPxBasis::INFEASIBLE :
      return LPSolver::INFEASIBLE;
   default:
      return LPSolver::ERROR;
   }
}

LPSolver::Status SoPlex::getResult(
   double* p_value,
   Vector* p_primal,
   Vector* p_slacks,
   Vector* p_dual,
   Vector* reduCosts) const
{
   if (p_value)
      *p_value = this->value();
   if (p_primal)
      getPrimal(*p_primal);
   if (p_slacks)
      getSlacks(*p_slacks);
   if (p_dual)
      getDual(*p_dual);
   if (reduCosts)
      getRdCost(*reduCosts);
   return status();
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
