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
#pragma ident "@(#) $Id: spxsolve.cpp,v 1.16 2002/01/10 13:34:49 bzfpfend Exp $"

#include <assert.h>
#include <iostream>

#include "soplex.h"
#include "spxpricer.h"
#include "spxratiotester.h"
#include "spxstarter.h"
#include "spxsimplifier.h"

//#define DEBUG 1

#ifdef DEBUG
#undef NDEBUG
#endif // DEBUG

namespace soplex
{
SoPlex::Status SoPlex::solve()
{
   Id enterId;
   int leaveNum;

   if (dim() <= 0 && coDim() <= 0)          // no problem loaded
      return ERROR;

   if (slinSolver() == 0)             // linear system solver is required.
      return ERROR;

   if (thesimplifier)
   {
      // if (thesimplifier->loadedLP() != this)
      thesimplifier->load(this);

      switch (thesimplifier->simplify())
      {
      case 1:
         setStatus(SPxBasis::UNBOUNDED);
         return UNBOUNDED;
      case - 1:
         setStatus(SPxBasis::INFEASIBLE);
         return INFEASIBLE;
      default:
         break;
      }
   }
   if (thepricer == 0)                                // pricer is required.
      return ERROR;

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
      if (thestarter != 0)                         // no basis and no starter.
         thestarter->generate(*this);              // generate start basis.
      init();
   }
   thepricer->setEpsilon(delta());
   setType(type());

#ifdef DEBUG
   std::cout << "starting value = " << value() << '\n';
   std::cout << "starting shift = " << shift() << '\n';
   {
      int i;

      std::cout << "column status:\t";

      for (i = 0; i < desc().nCols(); ++i)
         std::cout << desc().colStatus(i);

      std::cout << "\nrow status:\t";

      for (i = 0; i < desc().nRows(); ++i)
         std::cout << desc().rowStatus(i);

      std::cout << std::endl;
   }
#endif  // DEBUG

   if (SPxBasis::status() == SPxBasis::OPTIMAL)
      setStatus(SPxBasis::REGULAR);
   bool stop = terminate();
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
#ifdef DEBUG
            std::cout << "iteration: " << iteration() 
                      << " value: " << value()
                      << " shift: " << shift()
                      << " epsilon: " << epsilon() << std::endl;
#endif // DEBUG

            if (shift() <= epsilon())
            {
               factorize();
               unShift();

#ifdef DEBUG
               std::cout << "maxInfeas: " << maxInfeas()
                         << " shift: " << shift()
                         << " delta: " << delta() << std::endl;
#endif // DEBUG

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
#ifdef DEBUG
            std::cout << "iteration: " << iteration() 
                      << " value: " << value()
                      << " shift: " << shift()
                      << " epsilon: " << epsilon() << std::endl;
#endif // DEBUG

            if (shift() <= epsilon())
            {
               factorize();
               unShift();

#ifdef DEBUG
               std::cout << "maxInfeas: " << maxInfeas()
                         << " shift: " << shift()
                         << " delta: " << delta() << std::endl;
#endif // DEBUG
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
      std::cout << iteration() << ":\tcoP error = " << tmp.length();
      tmp.clear();
      SPxBasis::coSolve(tmp, *theCoPrhs);
      multWithBase(tmp);
      tmp -= *theCoPrhs;
      std::cout << "\t[" << tmp.length() << "]\t(";
      tmp.clear();
      SPxBasis::coSolve(tmp, *theCoPrhs);
      tmp -= *theCoPvec;
      std::cout << tmp.length() << ")\n";
   }

   tmp = *theFvec;
   multBaseWith(tmp);
   tmp -= *theFrhs;
   if (tmp.length() > delta())
   {
      std::cout << iteration() << ":\t  F error = " << tmp.length() << "\t(";
      tmp.clear();
      SPxBasis::solve(tmp, *theFrhs);
      tmp -= *theFvec;
      std::cout << tmp.length() << ")\n";
   }

   if (type() == ENTER)
   {
      for (i = 0; i < dim(); ++i)
      {
         if (theCoTest[i] < 0 && isCoBasic(i))
            std::cout << "testVecs: theCoTest: this shalt not be!\n";
      }
      for (i = 0; i < coDim(); ++i)
      {
         if (theTest[i] < 0 && isBasic(i))
            std::cout << "testVecs: theTest: this shalt not be!\n";
      }
   }
}

bool SoPlex::terminate()
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
         std::cout << "unexpected change of coPrhs " 
                   << cr.length() << std::endl;
      if (fr.length() > delta())
         std::cout << "unexpected change of   Frhs " 
                   << fr.length() << std::endl;
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
   {
#ifndef NDEBUG
      std::cout << "Maximum number of iterations reached" << std::endl;
#endif // !NDEBUG

      return true;
   }
   if (maxTime >= 0 && time() >= maxTime)
   {
#ifndef NDEBUG
      std::cout << "Timelimit reached" << std::endl;
#endif // !NDEBUG

      return true;   
   }
   if (maxValue < SPxLP::infinity)
   {
      // This code is *NOT* tested.

      double objective = value();

      if (  ((spxSense() == SPxLP::MAXIMIZE) && (objective > maxValue))
         || ((spxSense() == SPxLP::MINIMIZE) && (objective < maxValue)))
      {
#ifndef NDEBUG
         std::cout << "Objective value limit reached" << std::endl;
#endif // !NDEBUG

         return true;
      }
   }
   return SPxBasis::status() >= SPxBasis::OPTIMAL
      || SPxBasis::status() <= SPxBasis::SINGULAR;
}

SoPlex::Status SoPlex::getPrimal (Vector& p_vector) const
{
   if (!isInitialized())
      /**@todo patch suggests returning ERROR instead of initializing */
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

SoPlex::Status SoPlex::getDual (Vector& p_vector) const
{
   if (!isInitialized())
      /**@todo patch suggests returning ERROR instead of initializing */
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

SoPlex::Status SoPlex::getRdCost (Vector& p_vector) const
{
   if (!isInitialized())
      /**@todo patch suggests returning ERROR instead of initializing */
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

SoPlex::Status SoPlex::getSlacks (Vector& p_vector) const
{
   if (!isInitialized())
      /**@todo patch suggests returning ERROR instead of initializing */
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

SoPlex::Status SoPlex::status() const
{
   switch (SPxBasis::status())
   {
   case SPxBasis::NO_PROBLEM :
      return UNKNOWN;
   case SPxBasis::SINGULAR :
      return ERROR;
   case SPxBasis::REGULAR :
      return UNKNOWN;
   case SPxBasis::DUAL :
      return DUAL;
   case SPxBasis::PRIMAL :
      return PRIMAL;
   case SPxBasis::OPTIMAL :
      return SOLVED;
   case SPxBasis::UNBOUNDED :
      return UNBOUNDED;
   case SPxBasis::INFEASIBLE :
      return INFEASIBLE;
   default:
      return ERROR;
   }
}

SoPlex::Status SoPlex::getResult(
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
