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
#pragma ident "@(#) $Id: leave.cpp,v 1.10 2002/01/19 18:59:15 bzfkocht Exp $"

/* Updating the Basis for Leaving Variables
 */
#include        <assert.h>
#include        <stdio.h>
#include "real.h"
#include        "soplex.h"


#include "spxratiotester.h"

namespace soplex
{
static const Real reject_leave_tol = 1e-8;

/*
    Vector |fTest| gives the feasibility test of all basic variables. For its
    compution |fVec|, |theUBbound| and |theLBbound| must be setup correctly.
    Values of |fTest| $<0$ represent infeasible variables, which are eligable
    for leaving the basis in the simplex loop.
 */
void SoPlex::computeFtest()
{
   assert(type() == LEAVE);
   Vector& ftest = theCoTest;                  // |== fTest()|
   assert(&ftest == &fTest());

   for (int i = dim() - 1; i >= 0; --i)
   {
      ftest[i] = ((*theFvec)[i] > theUBbound[i])
         ? theUBbound[i] - (*theFvec)[i]
         : (*theFvec)[i] - theLBbound[i];
   }
}

void SoPlex::updateFtest()
{
   const IdxSet& idx = theFvec->idx();
   Vector& ftest = theCoTest;      // |== fTest()|
   assert(&ftest == &fTest());

   assert(type() == LEAVE);
   for (int j = idx.size() - 1; j >= 0; --j)
   {
      int i = idx.index(j);

      ftest[i] = ((*theFvec)[i] > theUBbound[i])
         ? theUBbound[i] - (*theFvec)[i]
         : (*theFvec)[i] - theLBbound[i];
   }
}


/* compute statistics on leaveing variable 
   Compute a set of statistical values on the variable selected for leaving the
   basis.
 */
void SoPlex::getLeaveVals
(
   int leaveIdx,
   SPxBasis::Desc::Status& leaveStat,
   Id& leaveId,
   Real& leaveMax,
   Real& leavebound,
   int& leaveNum
)
{
   SPxBasis::Desc& ds = desc();
   leaveId = baseId(leaveIdx);

   if (leaveId.isSPxRowId())
   {
      leaveNum = number(SPxRowId(leaveId));
      leaveStat = ds.rowStatus(leaveNum);
      //@ std::cerr << "R" << leaveNum << ":" << int(leaveStat);

      assert(isBasic(leaveStat));
      switch (leaveStat)
      {
      case SPxBasis::Desc::P_ON_UPPER :
         if (SPxLP::lhs(leaveNum) > -SPxLP::infinity)
            ds.rowStatus(leaveNum) = SPxBasis::Desc::D_ON_BOTH;
         else
            ds.rowStatus(leaveNum) = SPxBasis::Desc::D_ON_LOWER;
         leavebound = 0;
         leaveMax = -SPxLP::infinity;
         break;
      case SPxBasis::Desc::P_ON_LOWER :
         if (SPxLP::rhs(leaveNum) < SPxLP::infinity)
            ds.rowStatus(leaveNum) = SPxBasis::Desc::D_ON_BOTH;
         else
            ds.rowStatus(leaveNum) = SPxBasis::Desc::D_ON_UPPER;
         leavebound = 0;
         leaveMax = SPxLP::infinity;
         break;
      case SPxBasis::Desc::P_FREE :
         abort();

      case SPxBasis::Desc::D_FREE :
         ds.rowStatus(leaveNum) = SPxBasis::Desc::P_FIXED;
         assert(lhs(leaveNum) == rhs(leaveNum));
         leavebound = -rhs(leaveNum);
         if ((*theFvec)[leaveIdx] < theLBbound[leaveIdx])
            leaveMax = SPxLP::infinity;
         else
            leaveMax = -SPxLP::infinity;
         break;
      case SPxBasis::Desc::D_ON_LOWER :
         ds.rowStatus(leaveNum) = SPxBasis::Desc::P_ON_UPPER;
         leavebound = -rhs(leaveNum);                // slack !!
         leaveMax = SPxLP::infinity;
         break;
      case SPxBasis::Desc::D_ON_UPPER :
         ds.rowStatus(leaveNum) = SPxBasis::Desc::P_ON_LOWER;
         leavebound = -lhs(leaveNum);                // slack !!
         leaveMax = -SPxLP::infinity;
         break;
      case SPxBasis::Desc::D_ON_BOTH :
         if ((*theFvec)[leaveIdx] > theLBbound[leaveIdx])
         {
            ds.rowStatus(leaveNum) = SPxBasis::Desc::P_ON_LOWER;
            theLRbound[leaveNum] = -SPxLP::infinity;
            leavebound = -lhs(leaveNum);            // slack !!
            leaveMax = -SPxLP::infinity;
         }
         else
         {
            ds.rowStatus(leaveNum) = SPxBasis::Desc::P_ON_UPPER;
            theURbound[leaveNum] = SPxLP::infinity;
            leavebound = -rhs(leaveNum);            // slack !!
            leaveMax = SPxLP::infinity;
         }
         break;

      default:
         abort();
      }
   }

   else
   {
      assert(leaveId.isSPxColId());
      leaveNum = number(SPxColId(leaveId));
      leaveStat = ds.colStatus(leaveNum);
      //@ std::cerr << "C" << leaveNum << ":" << int(leaveStat);

      assert(isBasic(leaveStat));
      switch (leaveStat)
      {
      case SPxBasis::Desc::P_ON_UPPER :
         if (SPxLP::lower(leaveNum) > -SPxLP::infinity)
            ds.colStatus(leaveNum) = SPxBasis::Desc::D_ON_BOTH;
         else
            ds.colStatus(leaveNum) = SPxBasis::Desc::D_ON_LOWER;
         leavebound = 0;
         leaveMax = -SPxLP::infinity;
         break;
      case SPxBasis::Desc::P_ON_LOWER :
         if (SPxLP::upper(leaveNum) < SPxLP::infinity)
            ds.colStatus(leaveNum) = SPxBasis::Desc::D_ON_BOTH;
         else
            ds.colStatus(leaveNum) = SPxBasis::Desc::D_ON_UPPER;
         leavebound = 0;
         leaveMax = SPxLP::infinity;
         break;
      case SPxBasis::Desc::P_FREE :
         if ((*theFvec)[leaveIdx] < theLBbound[leaveIdx])
         {
            leavebound = theLBbound[leaveIdx];
            leaveMax = -SPxLP::infinity;
         }
         else
         {
            leavebound = theUBbound[leaveIdx];
            leaveMax = SPxLP::infinity;
         }
         ds.colStatus(leaveNum) = SPxBasis::Desc::D_UNDEFINED;
         break;

      case SPxBasis::Desc::D_FREE :
         assert(SPxLP::upper(leaveNum) == SPxLP::lower(leaveNum));
         ds.colStatus(leaveNum) = SPxBasis::Desc::P_FIXED;
         leavebound = SPxLP::upper(leaveNum);
         if ((*theFvec)[leaveIdx] < theLBbound[leaveIdx])
            leaveMax = SPxLP::infinity;
         else
            leaveMax = -SPxLP::infinity;
         break;
      case SPxBasis::Desc::D_ON_LOWER :
         ds.colStatus(leaveNum) = SPxBasis::Desc::P_ON_UPPER;
         leavebound = SPxLP::upper(leaveNum);
         leaveMax = -SPxLP::infinity;
         break;
      case SPxBasis::Desc::D_ON_UPPER :
         ds.colStatus(leaveNum) = SPxBasis::Desc::P_ON_LOWER;
         leavebound = SPxLP::lower(leaveNum);
         leaveMax = SPxLP::infinity;
         break;
      case SPxBasis::Desc::D_ON_BOTH :
         if ((*theFvec)[leaveIdx] > theUBbound[leaveIdx])
         {
            leaveMax = -SPxLP::infinity;
            leavebound = SPxLP::upper(leaveNum);
            theLCbound[leaveNum] = -SPxLP::infinity;
            ds.colStatus(leaveNum) = SPxBasis::Desc::P_ON_UPPER;
         }
         else
         {
            leaveMax = SPxLP::infinity;
            leavebound = SPxLP::lower(leaveNum);
            theUCbound[leaveNum] = SPxLP::infinity;
            ds.colStatus(leaveNum) = SPxBasis::Desc::P_ON_LOWER;
         }
         break;

      default:
         abort();
      }
   }
}

void SoPlex::getLeaveVals2(
   Real leaveMax,
   Id enterId,
   Real& enterBound,
   Real& newUBbound,
   Real& newLBbound,
   Real& newCoPrhs
)
{
   SPxBasis::Desc& ds = desc();

   enterBound = 0;
   if (enterId.isSPxRowId())
   {
      int idx = number(SPxRowId(enterId));
      //@ std::cerr << "\t->\tC" << idx << ": " << int(ds.rowStatus(idx)) << '\n';
      switch (ds.rowStatus(idx))
      {
      case SPxBasis::Desc::D_FREE :
         assert(rep() == ROW);
         if (thePvec->delta()[idx] * leaveMax < 0)
            newCoPrhs = theLRbound[idx];
         else
            newCoPrhs = theURbound[idx];
         newUBbound = SPxLP::infinity;
         newLBbound = -SPxLP::infinity;
         ds.rowStatus(idx) = SPxBasis::Desc::P_FIXED;
         break;
      case SPxBasis::Desc::D_ON_UPPER :
         assert(rep() == ROW);
         newUBbound = 0;
         newLBbound = -SPxLP::infinity;
         ds.rowStatus(idx) = SPxBasis::Desc::P_ON_LOWER;
         newCoPrhs = theLRbound[idx];
         break;
      case SPxBasis::Desc::D_ON_LOWER :
         assert(rep() == ROW);
         newUBbound = SPxLP::infinity;
         newLBbound = 0;
         ds.rowStatus(idx) = SPxBasis::Desc::P_ON_UPPER;
         newCoPrhs = theURbound[idx];
         break;
      case SPxBasis::Desc::D_ON_BOTH :
         assert(rep() == ROW);
         if (leaveMax * thePvec->delta()[idx] < 0)
         {
            newUBbound = 0;
            newLBbound = -SPxLP::infinity;
            ds.rowStatus(idx) = SPxBasis::Desc::P_ON_LOWER;
            newCoPrhs = theLRbound[idx];
         }
         else
         {
            newUBbound = SPxLP::infinity;
            newLBbound = 0;
            ds.rowStatus(idx) = SPxBasis::Desc::P_ON_UPPER;
            newCoPrhs = theURbound[idx];
         }
         break;

      case SPxBasis::Desc::P_ON_UPPER :
         if (lhs(idx) > -SPxLP::infinity)
         {
            ds.rowStatus(idx) = SPxBasis::Desc::D_ON_BOTH;
            theURbound[idx] = theLRbound[idx];
            //@             theURbound[idx]   = 0;
         }

         else
            ds.rowStatus(idx) = SPxBasis::Desc::D_ON_LOWER;
         newCoPrhs = theLRbound[idx];        // slack !!
         newUBbound = -lhs(idx);
         newLBbound = -rhs(idx);
         enterBound = -rhs(idx);
         break;
      case SPxBasis::Desc::P_ON_LOWER :
         if (rhs(idx) < SPxLP::infinity)
         {
            ds.rowStatus(idx) = SPxBasis::Desc::D_ON_BOTH;
            theLRbound[idx] = theURbound[idx];
            //@             theLRbound[idx]   = 0;
         }

         else
            ds.rowStatus(idx) = SPxBasis::Desc::D_ON_UPPER;
         newCoPrhs = theURbound[idx];        // slack !!
         newLBbound = -rhs(idx);
         newUBbound = -lhs(idx);
         enterBound = -lhs(idx);
         break;
      case SPxBasis::Desc::P_FREE :
         abort();
#if 0
         ds.rowStatus(idx) = SPxBasis::Desc::D_UNDEFINED;
         std::cerr << __FILE__ << __LINE__ << "ERROR: not yet debugged!\n";

         newCoPrhs = theURbound[idx];        // slack !!
         newUBbound = SPxLP::infinity;
         newLBbound = -SPxLP::infinity;
         enterBound = 0;
#endif
         break;
      case SPxBasis::Desc::P_FIXED :
         abort();

      default:
         abort();
         break;
      }
   }

   else
   {
      assert(enterId.isSPxColId());
      int idx = number(SPxColId(enterId));
      //@ std::cerr << "\t->\tC" << idx << ": " << int(ds.colStatus(idx)) << '\n';
      switch (ds.colStatus(idx))
      {
      case SPxBasis::Desc::D_ON_UPPER :
         assert(rep() == ROW);
         newUBbound = 0;
         newLBbound = -SPxLP::infinity;
         ds.colStatus(idx) = SPxBasis::Desc::P_ON_LOWER;
         newCoPrhs = theLCbound[idx];
         break;
      case SPxBasis::Desc::D_ON_LOWER :
         assert(rep() == ROW);
         newUBbound = SPxLP::infinity;
         newLBbound = 0;
         ds.colStatus(idx) = SPxBasis::Desc::P_ON_UPPER;
         newCoPrhs = theUCbound[idx];
         break;
      case SPxBasis::Desc::D_FREE :
         assert(rep() == ROW);
         newUBbound = SPxLP::infinity;
         newLBbound = -SPxLP::infinity;
         newCoPrhs = theLCbound[idx];
         ds.colStatus(idx) = SPxBasis::Desc::P_FIXED;
         break;
      case SPxBasis::Desc::D_ON_BOTH :
         assert(rep() == ROW);
         if (leaveMax * theCoPvec->delta()[idx] < 0)
         {
            newUBbound = 0;
            newLBbound = -SPxLP::infinity;
            ds.colStatus(idx) = SPxBasis::Desc::P_ON_LOWER;
            newCoPrhs = theLCbound[idx];
         }
         else
         {
            newUBbound = SPxLP::infinity;
            newLBbound = 0;
            ds.colStatus(idx) = SPxBasis::Desc::P_ON_UPPER;
            newCoPrhs = theUCbound[idx];
         }
         break;

      case SPxBasis::Desc::P_ON_UPPER :
         if (SPxLP::lower(idx) > -SPxLP::infinity)
         {
            ds.colStatus(idx) = SPxBasis::Desc::D_ON_BOTH;
            theLCbound[idx] = theUCbound[idx];
            //@             theLCbound[idx]   = object[i];
         }

         else
            ds.colStatus(idx) = SPxBasis::Desc::D_ON_LOWER;
         newCoPrhs = theUCbound[idx];
         newUBbound = SPxLP::upper(idx);
         newLBbound = SPxLP::lower(idx);
         enterBound = SPxLP::upper(idx);
         break;
      case SPxBasis::Desc::P_ON_LOWER :
         if (SPxLP::upper(idx) < SPxLP::infinity)
         {
            ds.colStatus(idx) = SPxBasis::Desc::D_ON_BOTH;
            theUCbound[idx] = theLCbound[idx];
            //@             theUCbound[idx]   = object[i];
         }

         else
            ds.colStatus(idx) = SPxBasis::Desc::D_ON_UPPER;
         newCoPrhs = theLCbound[idx];
         newUBbound = SPxLP::upper(idx);
         newLBbound = SPxLP::lower(idx);
         enterBound = SPxLP::lower(idx);
         break;
      case SPxBasis::Desc::P_FREE :
         if (thePvec->delta()[idx] * leaveMax > 0)
            newCoPrhs = theUCbound[idx];
         else
            newCoPrhs = theLCbound[idx];
         ds.colStatus(idx) = SPxBasis::Desc::D_UNDEFINED;
         newUBbound = SPxLP::upper(idx);
         newLBbound = SPxLP::lower(idx);
         enterBound = 0;
         break;
      case SPxBasis::Desc::P_FIXED :
         abort();

      default:
         abort();
      }
   }

}

void SoPlex::rejectLeave(
   int leaveNum,
   Id leaveId,
   SPxBasis::Desc::Status leaveStat,
   const SVector* //newVec
)
{
   SPxBasis::Desc& ds = desc();
   if (leaveId.isSPxRowId())
   {
      if (leaveStat == SPxBasis::Desc::D_ON_BOTH)
      {
         if (ds.rowStatus(leaveNum) == SPxBasis::Desc::P_ON_LOWER)
            theLRbound[leaveNum] = theURbound[leaveNum];
         else
            theURbound[leaveNum] = theLRbound[leaveNum];
      }
      ds.rowStatus(leaveNum) = leaveStat;
   }
   else
   {
      if (leaveStat == SPxBasis::Desc::D_ON_BOTH)
      {
         if (ds.colStatus(leaveNum) == SPxBasis::Desc::P_ON_UPPER)
            theLCbound[leaveNum] = theUCbound[leaveNum];
         else
            theUCbound[leaveNum] = theLCbound[leaveNum];
      }
      ds.colStatus(leaveNum) = leaveStat;
   }
}


int SoPlex::leave(int leaveIdx)
{
   assert(leaveIdx < dim() && leaveIdx >= 0);
   assert(type() == LEAVE);
   assert(initialized);

   /*
       Before performing the actual basis update, we must determine, how this
       is to be accomplished.
    */
   if (theCoPvec->delta().isSetup() && theCoPvec->delta().size() == 0)
   {
      coSolve(theCoPvec->delta(), unitVecs[leaveIdx]);
   }
#ifndef NDEBUG
   else
   {
      SSVector tmp(dim());
      tmp.clear();
      coSolve(tmp, unitVecs[leaveIdx]);
      tmp -= theCoPvec->delta();
      if (tmp.length() > delta())
         std::cerr << basis().iteration() << ": coPvec.delta error = "
         << tmp.length() << std::endl;
   }
   if (!theCoPvec->isConsistent())
      std::cerr << "fuck\n";
#endif  // NDEBUG
   setupPupdate();

   assert(thePvec->isConsistent());
   assert(theCoPvec->isConsistent());

   SPxBasis::Desc::Status leaveStat;      // status of leaving var
   Id leaveId;        // id of leaving var
   Real leaveMax;       // maximium lambda of leaving var
   Real leavebound;     // current fVec value of leaving var
   int leaveNum;       // number of leaveId in bounds
   getLeaveVals(leaveIdx, leaveStat, leaveId, leaveMax, leavebound, leaveNum);

   if (m_numCycle > m_maxCycle)
   {
      if (leaveMax > 0)
         perturbMaxLeave();
      else
         perturbMinLeave();
      //@ m_numCycle /= 2;
   }
   //@ testBounds();


   for(;;)
   {
      Real enterVal = leaveMax;
      Id enterId = theratiotester->selectEnter(enterVal);

      /*
          No variable could be selected to enter the basis and even the leaving
          variable is unbounded --- this is a failure.
       */
      if (!enterId.isValid())
      {
         Id none;
         change(leaveIdx, none, 0);
         if (enterVal != leaveMax)
         {
            // std::cerr << "rejecting leave\n";
            rejectLeave(leaveNum, leaveId, leaveStat);
            theCoTest[leaveIdx] *= 0.01;            // #== fTest()#
            theCoTest[leaveIdx] -= 2 * delta();       // #== fTest()#
            return 1;
         }
         if (rep() != COLUMN)
            setStatus(SPxBasis::UNBOUNDED);
         else
            setStatus(SPxBasis::INFEASIBLE);
         return 0;
      }


      /*
          If an entering variable has been found, a regular basis update is to be
          performed.
       */
      else if (enterId != baseId(leaveIdx))
      {
         const SVector& newVector = *enterVector(enterId);


         // update feasibility vectors
         if (solveVector2)
            SPxBasis::solve4update (theFvec->delta(), *solveVector2,
                                     newVector, *solveVector2rhs);
         else
            SPxBasis::solve4update (theFvec->delta(), newVector);

#ifndef NDEBUG
         {
            SSVector tmp(dim());
            SPxBasis::solve(tmp, newVector);
            tmp -= fVec().delta();
            if (tmp.length() > delta())
               std::cerr << "\t(" << tmp.length() << ")\n";
         }
#endif  // NDEBUG


         if (fabs(theFvec->delta()[leaveIdx]) < reject_leave_tol)
         {
            Id none;
            change(leaveIdx, none, 0);
            theFvec->delta().clear();
            rejectLeave(leaveNum, leaveId, leaveStat, &newVector);
            std::cerr << "rejecting leave\n";
            // factorize();
            theCoTest[leaveIdx] *= 0.01;            // #== fTest()#
            return 1;
         }

         //      process leaving variable
         if (leavebound > epsilon() || leavebound < -epsilon())
            theFrhs->multAdd(-leavebound, baseVec(leaveIdx));



         //      process entering variable
         Real enterBound;
         Real newUBbound;
         Real newLBbound;
         Real newCoPrhs;

         getLeaveVals2(leaveMax, enterId,
                        enterBound, newUBbound, newLBbound, newCoPrhs);

         theUBbound[leaveIdx] = newUBbound;
         theLBbound[leaveIdx] = newLBbound;
         (*theCoPrhs)[leaveIdx] = newCoPrhs;

         if (enterBound > epsilon() || enterBound < -epsilon())
            theFrhs->multAdd(enterBound, newVector);

         // update pricing vectors
         theCoPvec->value() = enterVal;
         thePvec->value() = enterVal;
         if (enterVal > epsilon() || enterVal < -epsilon())
            doPupdate();


         // update feasibility vector
         theFvec->value() = -((*theFvec)[leaveIdx] - leavebound)
                            / theFvec->delta()[leaveIdx];
         theFvec->update();
         (*theFvec)[leaveIdx] = enterBound - theFvec->value();
         updateFtest();


         //  change basis matrix
         change(leaveIdx, enterId, &newVector, &(theFvec->delta()));
      }


      /*
          No entering vector has been selected from the basis. However, if the
          shift amount for |coPvec| is bounded, we are in the case, that the
          entering variable is moved from one bound to its other, before any of
          the basis feasibility variables reaches their bound. This may only
          happen in primal/columnwise case with upper and lower bounds on
          variables.
       */
      else
      {
         assert(rep() == ROW);
         SPxBasis::Desc& ds = desc();


         Id none;
         change(leaveIdx, none, 0);

         if (leaveStat == SPxBasis::Desc::P_ON_UPPER)
         {
            if (leaveId.isSPxRowId())
            {
               ds.rowStatus(leaveNum) = SPxBasis::Desc::P_ON_LOWER;
               (*theCoPrhs)[leaveIdx] = theLRbound[leaveNum];
            }
            else
            {
               ds.colStatus(leaveNum) = SPxBasis::Desc::P_ON_LOWER;
               (*theCoPrhs)[leaveIdx] = theLCbound[leaveNum];
            }
            theUBbound[leaveIdx] = 0;
            theLBbound[leaveIdx] = -SPxLP::infinity;
         }
         else
         {
            if (leaveId.isSPxRowId())
            {
               ds.rowStatus(leaveNum) = SPxBasis::Desc::P_ON_UPPER;
               (*theCoPrhs)[leaveIdx] = theURbound[leaveNum];
            }
            else
            {
               ds.colStatus(leaveNum) = SPxBasis::Desc::P_ON_UPPER;
               (*theCoPrhs)[leaveIdx] = theUCbound[leaveNum];
            }
            theUBbound[leaveIdx] = SPxLP::infinity;
            theLBbound[leaveIdx] = 0;
         }


         // update copricing vector
         theCoPvec->value() = enterVal;
         thePvec->value() = enterVal;
         if (enterVal > epsilon() || enterVal < -epsilon())
            doPupdate();


         // update feasibility vectors
         theFvec->value() = 0;
         theCoTest[leaveIdx] *= -1;
      }

      if ((leaveMax > delta() && enterVal <= delta())
           || (leaveMax < -delta() && enterVal >= -delta()))
      {
         m_numCycle += ((theUBbound[leaveIdx] < SPxLP::infinity ||
                        theLBbound[leaveIdx] > -SPxLP::infinity)
                        && leaveStat != SPxBasis::Desc::P_FREE
                        && leaveStat != SPxBasis::Desc::D_FREE);
      }
      else
         m_numCycle /= 2;

#ifndef NDEBUG
      {
         DVector tmp = fVec();
         multBaseWith(tmp);
         tmp -= fRhs();
         if (tmp.length() > delta())
         {
            std::cerr << '\t' << basis().iteration()
            << ": fVec error = " << tmp.length();
            SPxBasis::solve(tmp, fRhs());
            tmp -= fVec();
            std::cerr << "\t(" << tmp.length() << ")\n";
         }
      }
#endif  // NDEBUG

      return 1;
   }
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
