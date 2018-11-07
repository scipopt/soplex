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

/* Updating the Basis for Leaving Variables
 */
#include <assert.h>
#include <stdio.h>

#include "soplex/spxdefines.h"
#include "soplex/spxpricer.h"
#include "soplex/spxsolver.h"
#include "soplex/spxratiotester.h"
#include "soplex/spxout.h"
#include "soplex/exceptions.h"

namespace soplex
{
static const Real reject_leave_tol = 1e-10; // = LOWSTAB as defined in spxfastrt.cpp

/*
  Vector |fTest| gives the feasibility test of all basic variables. For its
  computation |fVec|, |theUBbound| and |theLBbound| must be setup correctly.
  Values of |fTest| $<0$ represent infeasible variables, which are eligible
  for leaving the basis in the simplex loop.
*/
template <>
void SPxSolverBase<Real>::computeFtest()
{

   assert(type() == LEAVE);

   Real theeps = entertol();
   m_pricingViolUpToDate = true;
   m_pricingViolCoUpToDate = true;
   m_pricingViol = 0;
   m_pricingViolCo = 0;
   infeasibilities.clear();
   int ninfeasibilities = 0;
   int sparsitythreshold = (int)(sparsePricingFactor * dim());

   for(int i = 0; i < dim(); ++i)
   {
      theCoTest[i] = ((*theFvec)[i] > theUBbound[i])
                     ? theUBbound[i] - (*theFvec)[i]
                     : (*theFvec)[i] - theLBbound[i];

      if(remainingRoundsLeave == 0)
      {
         if(theCoTest[i] < -theeps)
         {
            m_pricingViol -= theCoTest[i];
            infeasibilities.addIdx(i);
            isInfeasible[i] = SPxPricer<Real>::VIOLATED;
            ++ninfeasibilities;
         }
         else
            isInfeasible[i] = SPxPricer<Real>::NOT_VIOLATED;

         if(ninfeasibilities > sparsitythreshold)
         {
            MSG_INFO2((*spxout), (*spxout) << " --- using dense pricing"
                      << std::endl;)
            remainingRoundsLeave = DENSEROUNDS;
            sparsePricingLeave = false;
            ninfeasibilities = 0;
         }
      }
      else if(theCoTest[i] < -theeps)
         m_pricingViol -= theCoTest[i];
   }

   if(ninfeasibilities == 0 && !sparsePricingLeave)
   {
      --remainingRoundsLeave;
   }
   else if(ninfeasibilities <= sparsitythreshold && !sparsePricingLeave)
   {
      MSG_INFO2((*spxout),
                std::streamsize prec = spxout->precision();

                if(hyperPricingLeave)
                (*spxout) << " --- using hypersparse pricing, ";
                else
                   (*spxout) << " --- using sparse pricing, ";
                   (*spxout) << "sparsity: "
                   << std::setw(6) << std::fixed << std::setprecision(4)
                   << (Real) ninfeasibilities / dim()
                   << std::scientific << std::setprecision(int(prec))
                   << std::endl;
                  )
            sparsePricingLeave = true;
   }
}

template <>
void SPxSolverBase<Real>::updateFtest()
{
   const IdxSet& idx = theFvec->idx();
   Vector& ftest = theCoTest;      // |== fTest()|
   assert(&ftest == &fTest());

   assert(type() == LEAVE);

   updateViols.clear();
   Real theeps = entertol();

   for(int j = idx.size() - 1; j >= 0; --j)
   {
      int i = idx.index(j);

      if(m_pricingViolUpToDate && ftest[i] < -theeps)
         m_pricingViol += ftest[i];

      ftest[i] = ((*theFvec)[i] > theUBbound[i])
                 ? theUBbound[i] - (*theFvec)[i]
                 : (*theFvec)[i] - theLBbound[i];


      if(sparsePricingLeave && ftest[i] < -theeps)
      {
         assert(remainingRoundsLeave == 0);

         if(m_pricingViolUpToDate)
            m_pricingViol -= ftest[i];

         if(isInfeasible[i] == SPxPricer<Real>::NOT_VIOLATED)
         {
            // this can cause problems - we cannot keep on adding indeces to infeasibilities,
            // because they are not deleted in hyper mode...
            //             if( !hyperPricingLeave )
            infeasibilities.addIdx(i);
            isInfeasible[i] = SPxPricer<Real>::VIOLATED;
         }

         if(hyperPricingLeave)
            updateViols.addIdx(i);
      }
      else if(m_pricingViolUpToDate && ftest[i] < -theeps)
         m_pricingViol -= ftest[i];

   }

   // if boundflips were performed, we need to update these indices as well
   if(boundflips > 0)
   {
      Real eps = epsilon();

      for(int j = 0; j < solveVector3->size(); ++j)
      {
         if(spxAbs(solveVector3->value(j)) > eps)
         {
            int i = solveVector3->index(j);

            if(m_pricingViolUpToDate && ftest[i] < -theeps)
               m_pricingViol += ftest[i];

            ftest[i] = ((*theFvec)[i] > theUBbound[i]) ? theUBbound[i] - (*theFvec)[i] :
                       (*theFvec)[i] - theLBbound[i];

            if(sparsePricingLeave && ftest[i] < -theeps)
            {
               assert(remainingRoundsLeave == 0);

               if(m_pricingViolUpToDate)
                  m_pricingViol -= ftest[i];

               if(!isInfeasible[i])
               {
                  infeasibilities.addIdx(i);
                  isInfeasible[i] = true;
               }
            }
            else if(m_pricingViolUpToDate && ftest[i] < -theeps)
               m_pricingViol -= ftest[i];
         }
      }
   }
}


/* compute statistics on leaving variable
   Compute a set of statistical values on the variable selected for leaving the
   basis.
*/
template <>
void SPxSolverBase<Real>::getLeaveVals(
   int leaveIdx,
   typename SPxBasisBase<Real>::Desc::Status& leaveStat,
   SPxId& leaveId,
   Real& leaveMax,
   Real& leavebound,
   int& leaveNum,
   Real& objChange)
{
   typename SPxBasisBase<Real>::Desc& ds = this->desc();
   leaveId = this->baseId(leaveIdx);

   if(leaveId.isSPxRowId())
   {
      leaveNum = this->number(SPxRowId(leaveId));
      leaveStat = ds.rowStatus(leaveNum);

      assert(isBasic(leaveStat));

      switch(leaveStat)
      {
      case SPxBasisBase<Real>::Desc::P_ON_UPPER :
         assert(rep() == ROW);
         ds.rowStatus(leaveNum) = this->dualRowStatus(leaveNum);
         leavebound = 0;
         leaveMax = -infinity;
         break;

      case SPxBasisBase<Real>::Desc::P_ON_LOWER :
         assert(rep() == ROW);
         ds.rowStatus(leaveNum) = this->dualRowStatus(leaveNum);
         leavebound = 0;
         leaveMax = infinity;
         break;

      case SPxBasisBase<Real>::Desc::P_FREE :
         assert(rep() == ROW);
         throw SPxInternalCodeException("XLEAVE01 This should never happen.");

      case SPxBasisBase<Real>::Desc::D_FREE :
         assert(rep() == COLUMN);
         ds.rowStatus(leaveNum) = SPxBasisBase<Real>::Desc::P_FIXED;
         assert(lhs(leaveNum) == rhs(leaveNum));
         leavebound = -this->rhs(leaveNum);

         if((*theFvec)[leaveIdx] < theLBbound[leaveIdx])
            leaveMax = infinity;
         else
            leaveMax = -infinity;

         break;

      case SPxBasisBase<Real>::Desc::D_ON_LOWER :
         assert(rep() == COLUMN);
         ds.rowStatus(leaveNum) = SPxBasisBase<Real>::Desc::P_ON_UPPER;
         leavebound = -this->rhs(leaveNum);                // slack !!
         leaveMax = infinity;
         objChange += theLRbound[leaveNum] * this->rhs(leaveNum);
         break;

      case SPxBasisBase<Real>::Desc::D_ON_UPPER :
         assert(rep() == COLUMN);
         ds.rowStatus(leaveNum) = SPxBasisBase<Real>::Desc::P_ON_LOWER;
         leavebound = -this->lhs(leaveNum);                // slack !!
         leaveMax = -infinity;
         objChange += theURbound[leaveNum] * this->lhs(leaveNum);
         break;

      case SPxBasisBase<Real>::Desc::D_ON_BOTH :
         assert(rep() == COLUMN);

         if((*theFvec)[leaveIdx] > theLBbound[leaveIdx])
         {
            ds.rowStatus(leaveNum) = SPxBasisBase<Real>::Desc::P_ON_LOWER;
            theLRbound[leaveNum] = -infinity;
            leavebound = -this->lhs(leaveNum);            // slack !!
            leaveMax = -infinity;
            objChange += theURbound[leaveNum] * this->lhs(leaveNum);
         }
         else
         {
            ds.rowStatus(leaveNum) = SPxBasisBase<Real>::Desc::P_ON_UPPER;
            theURbound[leaveNum] = infinity;
            leavebound = -this->rhs(leaveNum);            // slack !!
            leaveMax = infinity;
            objChange += theLRbound[leaveNum] * this->rhs(leaveNum);
         }

         break;

      default:
         throw SPxInternalCodeException("XLEAVE02 This should never happen.");
      }

      MSG_DEBUG(std::cout << "DLEAVE51 SPxSolverBase<Real>::getLeaveVals() : row " << leaveNum
                << ": " << leaveStat
                << " -> " << ds.rowStatus(leaveNum)
                << " objChange: " << objChange
                << std::endl;)
   }

   else
   {
      assert(leaveId.isSPxColId());
      leaveNum = this->number(SPxColId(leaveId));
      leaveStat = ds.colStatus(leaveNum);

      assert(isBasic(leaveStat));

      switch(leaveStat)
      {
      case SPxBasisBase<Real>::Desc::P_ON_UPPER :
         assert(rep() == ROW);
         ds.colStatus(leaveNum) = this->dualColStatus(leaveNum);
         leavebound = 0;
         leaveMax = -infinity;
         break;

      case SPxBasisBase<Real>::Desc::P_ON_LOWER :
         assert(rep() == ROW);
         ds.colStatus(leaveNum) = this->dualColStatus(leaveNum);
         leavebound = 0;
         leaveMax = infinity;
         break;

      case SPxBasisBase<Real>::Desc::P_FREE :
         assert(rep() == ROW);
         ds.colStatus(leaveNum) = this->dualColStatus(leaveNum);

         if((*theFvec)[leaveIdx] < theLBbound[leaveIdx])
         {
            leavebound = theLBbound[leaveIdx];
            leaveMax = -infinity;
         }
         else
         {
            leavebound = theUBbound[leaveIdx];
            leaveMax = infinity;
         }

         break;

      case SPxBasisBase<Real>::Desc::D_FREE :
         assert(rep() == COLUMN);
         assert(SPxLP::upper(leaveNum) == SPxLP::lower(leaveNum));
         ds.colStatus(leaveNum) = SPxBasisBase<Real>::Desc::P_FIXED;
         leavebound = SPxLP::upper(leaveNum);
         objChange += this->maxObj(leaveNum) * leavebound;

         if((*theFvec)[leaveIdx] < theLBbound[leaveIdx])
            leaveMax = infinity;
         else
            leaveMax = -infinity;

         break;

      case SPxBasisBase<Real>::Desc::D_ON_LOWER :
         assert(rep() == COLUMN);
         ds.colStatus(leaveNum) = SPxBasisBase<Real>::Desc::P_ON_UPPER;
         leavebound = SPxLP::upper(leaveNum);
         objChange += theUCbound[leaveNum] * leavebound;
         leaveMax = -infinity;
         break;

      case SPxBasisBase<Real>::Desc::D_ON_UPPER :
         assert(rep() == COLUMN);
         ds.colStatus(leaveNum) = SPxBasisBase<Real>::Desc::P_ON_LOWER;
         leavebound = SPxLP::lower(leaveNum);
         objChange += theLCbound[leaveNum] * leavebound;
         leaveMax = infinity;
         break;

      case SPxBasisBase<Real>::Desc::D_ON_BOTH :
         assert(rep() == COLUMN);

         if((*theFvec)[leaveIdx] > theUBbound[leaveIdx])
         {
            leaveMax = -infinity;
            leavebound = SPxLP::upper(leaveNum);
            objChange += theUCbound[leaveNum] * leavebound;
            theLCbound[leaveNum] = -infinity;
            ds.colStatus(leaveNum) = SPxBasisBase<Real>::Desc::P_ON_UPPER;
         }
         else
         {
            leaveMax = infinity;
            leavebound = SPxLP::lower(leaveNum);
            objChange += theLCbound[leaveNum] * leavebound;
            theUCbound[leaveNum] = infinity;
            ds.colStatus(leaveNum) = SPxBasisBase<Real>::Desc::P_ON_LOWER;
         }

         break;

      default:
         throw SPxInternalCodeException("XLEAVE03 This should never happen.");
      }

      MSG_DEBUG(std::cout << "DLEAVE52 SPxSolverBase<Real>::getLeaveVals() : col " << leaveNum
                << ": " << leaveStat
                << " -> " << ds.colStatus(leaveNum)
                << " objChange: " << objChange
                << std::endl;)
   }
}

template <>
void SPxSolverBase<Real>::getLeaveVals2(
   Real leaveMax,
   SPxId enterId,
   Real& enterBound,
   Real& newUBbound,
   Real& newLBbound,
   Real& newCoPrhs,
   Real& objChange
)
{
   typename SPxBasisBase<Real>::Desc& ds = this->desc();

   enterBound = 0;

   if(enterId.isSPxRowId())
   {
      int idx = this->number(SPxRowId(enterId));
      typename SPxBasisBase<Real>::Desc::Status enterStat = ds.rowStatus(idx);

      switch(enterStat)
      {
      case SPxBasisBase<Real>::Desc::D_FREE :
         assert(rep() == ROW);

         if(thePvec->delta()[idx] * leaveMax < 0)
            newCoPrhs = theLRbound[idx];
         else
            newCoPrhs = theURbound[idx];

         newUBbound = infinity;
         newLBbound = -infinity;
         ds.rowStatus(idx) = SPxBasisBase<Real>::Desc::P_FIXED;
         break;

      case SPxBasisBase<Real>::Desc::D_ON_UPPER :
         assert(rep() == ROW);
         newUBbound = 0;
         newLBbound = -infinity;
         ds.rowStatus(idx) = SPxBasisBase<Real>::Desc::P_ON_LOWER;
         newCoPrhs = theLRbound[idx];
         break;

      case SPxBasisBase<Real>::Desc::D_ON_LOWER :
         assert(rep() == ROW);
         newUBbound = infinity;
         newLBbound = 0;
         ds.rowStatus(idx) = SPxBasisBase<Real>::Desc::P_ON_UPPER;
         newCoPrhs = theURbound[idx];
         break;

      case SPxBasisBase<Real>::Desc::D_ON_BOTH :
         assert(rep() == ROW);

         if(leaveMax * thePvec->delta()[idx] < 0)
         {
            newUBbound = 0;
            newLBbound = -infinity;
            ds.rowStatus(idx) = SPxBasisBase<Real>::Desc::P_ON_LOWER;
            newCoPrhs = theLRbound[idx];
         }
         else
         {
            newUBbound = infinity;
            newLBbound = 0;
            ds.rowStatus(idx) = SPxBasisBase<Real>::Desc::P_ON_UPPER;
            newCoPrhs = theURbound[idx];
         }

         break;

      case SPxBasisBase<Real>::Desc::P_ON_UPPER :
         assert(rep() == COLUMN);
         ds.rowStatus(idx) = this->dualRowStatus(idx);

         if(this->lhs(idx) > -infinity)
            theURbound[idx] = theLRbound[idx];

         newCoPrhs = theLRbound[idx];        // slack !!
         newUBbound = -this->lhs(idx);
         newLBbound = -this->rhs(idx);
         enterBound = -this->rhs(idx);
         objChange -= newCoPrhs * this->rhs(idx);
         break;

      case SPxBasisBase<Real>::Desc::P_ON_LOWER :
         assert(rep() == COLUMN);
         ds.rowStatus(idx) = this->dualRowStatus(idx);

         if(this->rhs(idx) < infinity)
            theLRbound[idx] = theURbound[idx];

         newCoPrhs = theURbound[idx];        // slack !!
         newLBbound = -this->rhs(idx);
         newUBbound = -this->lhs(idx);
         enterBound = -this->lhs(idx);
         objChange -= newCoPrhs * this->lhs(idx);
         break;

      case SPxBasisBase<Real>::Desc::P_FREE :
         assert(rep() == COLUMN);
#if 1
         throw SPxInternalCodeException("XLEAVE04 This should never happen.");
#else
         MSG_ERROR(std::cerr << "ELEAVE53 ERROR: not yet debugged!" << std::endl;)
         ds.rowStatus(idx) = this->dualRowStatus(idx);
         newCoPrhs = theURbound[idx];        // slack !!
         newUBbound = infinity;
         newLBbound = -infinity;
         enterBound = 0;
#endif
         break;

      case SPxBasisBase<Real>::Desc::P_FIXED :
         assert(rep() == COLUMN);
         MSG_ERROR(std::cerr << "ELEAVE54 "
                   << "ERROR! Tried to put a fixed row variable into the basis: "
                   << "idx="   << idx
                   << ", lhs=" << this->lhs(idx)
                   << ", rhs=" << this->rhs(idx) << std::endl;)
         throw SPxInternalCodeException("XLEAVE05 This should never happen.");

      default:
         throw SPxInternalCodeException("XLEAVE06 This should never happen.");
      }

      MSG_DEBUG(std::cout << "DLEAVE55 SPxSolverBase<Real>::getLeaveVals2(): row " << idx
                << ": " << enterStat
                << " -> " << ds.rowStatus(idx)
                << " objChange: " << objChange
                << std::endl;)
   }

   else
   {
      assert(enterId.isSPxColId());
      int idx = this->number(SPxColId(enterId));
      typename SPxBasisBase<Real>::Desc::Status enterStat = ds.colStatus(idx);

      switch(enterStat)
      {
      case SPxBasisBase<Real>::Desc::D_ON_UPPER :
         assert(rep() == ROW);
         newUBbound = 0;
         newLBbound = -infinity;
         ds.colStatus(idx) = SPxBasisBase<Real>::Desc::P_ON_LOWER;
         newCoPrhs = theLCbound[idx];
         break;

      case SPxBasisBase<Real>::Desc::D_ON_LOWER :
         assert(rep() == ROW);
         newUBbound = infinity;
         newLBbound = 0;
         ds.colStatus(idx) = SPxBasisBase<Real>::Desc::P_ON_UPPER;
         newCoPrhs = theUCbound[idx];
         break;

      case SPxBasisBase<Real>::Desc::D_FREE :
         assert(rep() == ROW);
         newUBbound = infinity;
         newLBbound = -infinity;
         newCoPrhs = theLCbound[idx];
         ds.colStatus(idx) = SPxBasisBase<Real>::Desc::P_FIXED;
         break;

      case SPxBasisBase<Real>::Desc::D_ON_BOTH :
         assert(rep() == ROW);

         if(leaveMax * theCoPvec->delta()[idx] < 0)
         {
            newUBbound = 0;
            newLBbound = -infinity;
            ds.colStatus(idx) = SPxBasisBase<Real>::Desc::P_ON_LOWER;
            newCoPrhs = theLCbound[idx];
         }
         else
         {
            newUBbound = infinity;
            newLBbound = 0;
            ds.colStatus(idx) = SPxBasisBase<Real>::Desc::P_ON_UPPER;
            newCoPrhs = theUCbound[idx];
         }

         break;

      case SPxBasisBase<Real>::Desc::P_ON_UPPER :
         assert(rep() == COLUMN);
         ds.colStatus(idx) = this->dualColStatus(idx);

         if(SPxLP::lower(idx) > -infinity)
            theLCbound[idx] = theUCbound[idx];

         newCoPrhs = theUCbound[idx];
         newUBbound = SPxLP::upper(idx);
         newLBbound = SPxLP::lower(idx);
         enterBound = SPxLP::upper(idx);
         objChange -= newCoPrhs * enterBound;
         break;

      case SPxBasisBase<Real>::Desc::P_ON_LOWER :
         assert(rep() == COLUMN);
         ds.colStatus(idx) = this->dualColStatus(idx);

         if(SPxLP::upper(idx) < infinity)
            theUCbound[idx] = theLCbound[idx];

         newCoPrhs = theLCbound[idx];
         newUBbound = SPxLP::upper(idx);
         newLBbound = SPxLP::lower(idx);
         enterBound = SPxLP::lower(idx);
         objChange -= newCoPrhs * enterBound;
         break;

      case SPxBasisBase<Real>::Desc::P_FREE :
         assert(rep() == COLUMN);
         ds.colStatus(idx) = this->dualColStatus(idx);

         if(thePvec->delta()[idx] * leaveMax > 0)
            newCoPrhs = theUCbound[idx];
         else
            newCoPrhs = theLCbound[idx];

         newUBbound = SPxLP::upper(idx);
         newLBbound = SPxLP::lower(idx);
         enterBound = 0;
         break;

      case SPxBasisBase<Real>::Desc::P_FIXED :
         assert(rep() == COLUMN);
         MSG_ERROR(std::cerr << "ELEAVE56 "
                   << "ERROR! Tried to put a fixed column variable into the basis. "
                   << "idx="     << idx
                   << ", lower=" << this->lower(idx)
                   << ", upper=" << this->upper(idx) << std::endl;)
         throw SPxInternalCodeException("XLEAVE07 This should never happen.");

      default:
         throw SPxInternalCodeException("XLEAVE08 This should never happen.");
      }

      MSG_DEBUG(std::cout << "DLEAVE57 SPxSolverBase<Real>::getLeaveVals2(): col " << idx
                << ": " << enterStat
                << " -> " << ds.colStatus(idx)
                << " objChange: " << objChange
                << std::endl;)
   }

}

template <>
void SPxSolverBase<Real>::rejectLeave(
   int leaveNum,
   SPxId leaveId,
   typename SPxBasisBase<Real>::Desc::Status leaveStat,
   const SVector* //newVec
)
{
   typename SPxBasisBase<Real>::Desc& ds = this->desc();

   if(leaveId.isSPxRowId())
   {
      MSG_DEBUG(std::cout << "DLEAVE58 rejectLeave()  : row " << leaveNum
                << ": " << ds.rowStatus(leaveNum)
                << " -> " << leaveStat << std::endl;)

      if(leaveStat == SPxBasisBase<Real>::Desc::D_ON_BOTH)
      {
         if(ds.rowStatus(leaveNum) == SPxBasisBase<Real>::Desc::P_ON_LOWER)
            theLRbound[leaveNum] = theURbound[leaveNum];
         else
            theURbound[leaveNum] = theLRbound[leaveNum];
      }

      ds.rowStatus(leaveNum) = leaveStat;
   }
   else
   {
      MSG_DEBUG(std::cout << "DLEAVE59 rejectLeave()  : col " << leaveNum
                << ": " << ds.colStatus(leaveNum)
                << " -> " << leaveStat << std::endl;)

      if(leaveStat == SPxBasisBase<Real>::Desc::D_ON_BOTH)
      {
         if(ds.colStatus(leaveNum) == SPxBasisBase<Real>::Desc::P_ON_UPPER)
            theLCbound[leaveNum] = theUCbound[leaveNum];
         else
            theUCbound[leaveNum] = theLCbound[leaveNum];
      }

      ds.colStatus(leaveNum) = leaveStat;
   }
}


template <>
void SPxSolverBase<Real>::computePrimalray4Row(Real direction)
{
   Real sign = (direction > 0 ? 1.0 : -1.0);

   primalRay.clear();
   primalRay.setMax(coPvec().delta().size());

   for(int i = 0; i < coPvec().delta().size(); ++i)
      primalRay.add(coPvec().delta().index(i), sign * coPvec().delta().value(i));
}

template <>
void SPxSolverBase<Real>::computeDualfarkas4Col(Real direction)
{
   Real sign = (direction > 0 ? -1.0 : 1.0);

   dualFarkas.clear();
   dualFarkas.setMax(coPvec().delta().size());

   for(int i = 0; i < coPvec().delta().size(); ++i)
      dualFarkas.add(coPvec().delta().index(i), sign * coPvec().delta().value(i));
}

template <>
bool SPxSolverBase<Real>::leave(int leaveIdx, bool polish)
{
   assert(leaveIdx < dim() && leaveIdx >= 0);
   assert(type() == LEAVE);
   assert(initialized);

   bool instable = instableLeave;
   assert(!instable || instableLeaveNum >= 0);

   /*
     Before performing the actual basis update, we must determine, how this
     is to be accomplished.
     When using steepest edge pricing this solve is already performed by the pricer
   */
   if(theCoPvec->delta().isSetup() && theCoPvec->delta().size() == 0)
   {
      this->coSolve(theCoPvec->delta(), unitVecs[leaveIdx]);
   }

#ifdef ENABLE_ADDITIONAL_CHECKS
   else
   {
      SSVector tmp(dim(), epsilon());
      tmp.clear();
      coSolve(tmp, unitVecs[leaveIdx]);
      tmp -= theCoPvec->delta();

      if(tmp.length() > leavetol())
      {
         // This happens very frequently and does usually not hurt, so print
         // these warnings only with verbose level INFO2 and higher.
         MSG_INFO2((*spxout), (*spxout) << "WLEAVE60 iteration=" << basis().iteration()
                   << ": coPvec.delta error = " << tmp.length()
                   << std::endl;)
      }
   }

#endif  // ENABLE_ADDITIONAL_CHECKS

   setupPupdate();

   assert(thePvec->isConsistent());
   assert(theCoPvec->isConsistent());

   typename SPxBasisBase<Real>::Desc::Status leaveStat;      // status of leaving var
   SPxId leaveId;        // id of leaving var
   SPxId none;           // invalid id used if leave fails
   Real leaveMax;       // maximium lambda of leaving var
   Real leavebound;     // current fVec value of leaving var
   int  leaveNum;       // number of leaveId in bounds
   Real objChange = 0.0; // amount of change in the objective function

   getLeaveVals(leaveIdx, leaveStat, leaveId, leaveMax, leavebound, leaveNum, objChange);

   if(!polish && m_numCycle > m_maxCycle)
   {
      if(leaveMax > 0)
         perturbMaxLeave();
      else
         perturbMinLeave();

      //@ m_numCycle /= 2;
      // perturbation invalidates the currently stored nonbasic value
      forceRecompNonbasicValue();
   }

   //@ testBounds();

   Real enterVal = leaveMax;
   boundflips = 0;
   Real oldShift = theShift;
   SPxId enterId = theratiotester->selectEnter(enterVal, leaveIdx, polish);

   if(NE(theShift, oldShift))
   {
      MSG_DEBUG(std::cout << "DLEAVE71 trigger recomputation of nonbasic value due to shifts in ratiotest"
                << std::endl;)
      forceRecompNonbasicValue();
   }

   assert(!enterId.isValid() || !isBasic(enterId));

   instableLeaveNum = -1;
   instableLeave = false;

   /*
     No variable could be selected to enter the basis and even the leaving
     variable is unbounded.
   */
   if(!enterId.isValid())
   {
      /* the following line originally was below in "rejecting leave" case;
         we need it in the unbounded/infeasible case, too, to have the
         correct basis size */
      rejectLeave(leaveNum, leaveId, leaveStat);
      this->change(-1, none, 0);
      objChange = 0.0; // the nonbasicValue is not supposed to be updated in this case

      if(polish)
         return false;

      if(NE(enterVal, leaveMax))
      {
         MSG_DEBUG(std::cout << "DLEAVE61 rejecting leave A (leaveIdx=" << leaveIdx
                   << ", theCoTest=" << theCoTest[leaveIdx] << ")"
                   << std::endl;)

         /* In the LEAVE algorithm, when for a selected leaving variable we find only
            an instable entering variable, then the basis change is not conducted.
            Instead, we save the leaving variable's index in instableLeaveNum and scale
            theCoTest[leaveIdx] down by some factor, hoping to find a different leaving
            variable with a stable entering variable.
            If this fails, however, and no more leaving variable is found, we have to
            perform the instable basis change using instableLeaveNum. In this (and only
            in this) case, the flag instableLeave is set to true.

            enterVal != leaveMax is the case that selectEnter has found only an instable entering
            variable. We store this leaving variable for later -- if we are not already in the
            instable case: then we continue and conclude unboundedness/infeasibility */
         if(!instable)
         {
            instableLeaveNum = leaveIdx;

            // Note: These changes do not survive a refactorization
            instableLeaveVal = theCoTest[leaveIdx];
            theCoTest[leaveIdx] = instableLeaveVal / 10.0;

            return true;
         }
      }

      if(this->lastUpdate() > 1)
      {
         MSG_INFO3((*spxout), (*spxout) << "ILEAVE01 factorization triggered in "
                   << "leave() for feasibility test" << std::endl;)

         try
         {
            factorize();
         }
         catch(const SPxStatusException& E)
         {
            // don't exit immediately but handle the singularity correctly
            assert(SPxBasisBase<Real>::status() == SPxBasisBase<Real>::SINGULAR);
            MSG_INFO3((*spxout), (*spxout) << "Caught exception in factorization: " << E.what() << std::endl;)
         }

         /* after a factorization, the leaving column/row might not be infeasible or suboptimal anymore, hence we do
          * not try to call leave(leaveIdx), but rather return to the main solving loop and call the pricer again
          */
         return true;
      }

      /* do not exit with status infeasible or unbounded if there is only a very small violation */
      if(!recomputedVectors && spxAbs(enterVal) < leavetol())
      {
         MSG_INFO3((*spxout), (*spxout) << "ILEAVE11 clean up step to reduce numerical errors" << std::endl;)

         computeFrhs();
         SPxBasisBase<Real>::solve(*theFvec, *theFrhs);
         computeFtest();

         /* only do this once per solve */
         recomputedVectors = true;

         return true;
      }

      MSG_INFO3((*spxout), (*spxout) << "ILEAVE02 unboundedness/infeasibility found "
                << "in leave()" << std::endl;)

      if(rep() != COLUMN)
      {
         computePrimalray4Row(enterVal);
         setBasisStatus(SPxBasisBase<Real>::UNBOUNDED);
      }
      else
      {
         computeDualfarkas4Col(enterVal);
         setBasisStatus(SPxBasisBase<Real>::INFEASIBLE);
      }

      return false;
   }
   else
   {
      /*
        If an entering variable has been found, a regular basis update is to
        be performed.
      */
      if(enterId != this->baseId((leaveIdx)))
      {
         const SVector& newVector = *enterVector(enterId);

         // update feasibility vectors
         if(solveVector2 != NULL && solveVector3 != NULL)
         {
            assert(solveVector2->isConsistent());
            assert(solveVector2rhs->isSetup());
            assert(solveVector3->isConsistent());
            assert(solveVector3rhs->isSetup());
            assert(boundflips > 0);
            SPxBasisBase<Real>::solve4update(theFvec->delta(),
                                             *solveVector2,
                                             *solveVector3,
                                             newVector,
                                             *solveVector2rhs,
                                             *solveVector3rhs);

            // perform update of basic solution
            primVec -= (*solveVector3);
            MSG_DEBUG(std::cout << "ILBFRT02 breakpoints passed / bounds flipped = " << boundflips << std::endl;
                     )
            totalboundflips += boundflips;
         }
         else if(solveVector2 != NULL)
         {
            assert(solveVector2->isConsistent());
            assert(solveVector2rhs->isSetup());

            SPxBasisBase<Real>::solve4update(theFvec->delta(),
                                             *solveVector2,
                                             newVector,
                                             *solveVector2rhs);
         }
         else if(solveVector3 != NULL)
         {
            assert(solveVector3->isConsistent());
            assert(solveVector3rhs->isSetup());
            assert(boundflips > 0);
            SPxBasisBase<Real>::solve4update(theFvec->delta(),
                                             *solveVector3,
                                             newVector,
                                             *solveVector3rhs);

            // perform update of basic solution
            primVec -= (*solveVector3);
            MSG_DEBUG(std::cout << "ILBFRT02 breakpoints passed / bounds flipped = " << boundflips << std::endl;
                     )
            totalboundflips += boundflips;
         }
         else
            SPxBasisBase<Real>::solve4update(theFvec->delta(), newVector);

#ifdef ENABLE_ADDITIONAL_CHECKS
         {
            SSVector tmp(dim(), epsilon());
            SPxBasisBase<Real>::solve(tmp, newVector);
            tmp -= fVec().delta();

            if(tmp.length() > entertol())
            {
               // This happens very frequently and does usually not hurt, so print
               // these warnings only with verbose level INFO2 and higher.
               MSG_INFO2((*spxout), (*spxout) << "WLEAVE62\t(" << tmp.length() << ")\n";)
            }
         }
#endif  // ENABLE_ADDITIONAL_CHECKS


         if(spxAbs(theFvec->delta()[leaveIdx]) < reject_leave_tol)
         {
            if(instable)
            {
               /* We are in the case that for all leaving variables only instable entering
                  variables were found: Thus, above we already accepted such an instable
                  entering variable. Now even this seems to be impossible, thus we conclude
                  unboundedness/infeasibility. */
               MSG_INFO3((*spxout), (*spxout) << "ILEAVE03 unboundedness/infeasibility found "
                         << "in leave()" << std::endl;)

               rejectLeave(leaveNum, leaveId, leaveStat);
               this->change(-1, none, 0);
               objChange = 0.0; // the nonbasicValue is not supposed to be updated in this case

               /**@todo if shift() is not zero we must not conclude unboundedness */
               if(rep() == ROW)
               {
                  computePrimalray4Row(enterVal);
                  setBasisStatus(SPxBasisBase<Real>::UNBOUNDED);
               }
               else
               {
                  computeDualfarkas4Col(enterVal);
                  setBasisStatus(SPxBasisBase<Real>::INFEASIBLE);
               }

               return false;
            }
            else
            {
               theFvec->delta().clear();
               rejectLeave(leaveNum, leaveId, leaveStat, &newVector);
               this->change(-1, none, 0);
               objChange = 0.0; // the nonbasicValue is not supposed to be updated in this case

               MSG_DEBUG(std::cout << "DLEAVE63 rejecting leave B (leaveIdx=" << leaveIdx
                         << ", theCoTest=" << theCoTest[leaveIdx]
                         << ")" << std::endl;)

               // Note: These changes do not survive a refactorization
               theCoTest[leaveIdx] *= 0.01;

               return true;
            }
         }

         //      process leaving variable
         if(leavebound > epsilon() || leavebound < -epsilon())
            theFrhs->multAdd(-leavebound, this->baseVec(leaveIdx));

         //      process entering variable
         Real enterBound;
         Real newUBbound;
         Real newLBbound;
         Real newCoPrhs;

         try
         {
            getLeaveVals2(leaveMax, enterId, enterBound, newUBbound, newLBbound, newCoPrhs, objChange);
         }
         catch(const SPxException& F)
         {
            rejectLeave(leaveNum, leaveId, leaveStat);
            this->change(-1, none, 0);
            objChange = 0.0; // the nonbasicValue is not supposed to be updated in this case
            throw F;
         }

         theUBbound[leaveIdx] = newUBbound;
         theLBbound[leaveIdx] = newLBbound;
         (*theCoPrhs)[leaveIdx] = newCoPrhs;

         if(enterBound > epsilon() || enterBound < -epsilon())
            theFrhs->multAdd(enterBound, newVector);

         // update pricing vectors
         theCoPvec->value() = enterVal;
         thePvec->value() = enterVal;

         if(enterVal > epsilon() || enterVal < -epsilon())
            doPupdate();

         // update feasibility vector
         theFvec->value() = -((*theFvec)[leaveIdx] - leavebound)
                            / theFvec->delta()[leaveIdx];
         theFvec->update();
         (*theFvec)[leaveIdx] = enterBound - theFvec->value();
         updateFtest();

         // update objective funtion value
         updateNonbasicValue(objChange);

         //  change basis matrix
         this->change(leaveIdx, enterId, &newVector, &(theFvec->delta()));
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
         // @todo update obj function value here!!!
         assert(rep() == ROW);
         typename SPxBasisBase<Real>::Desc& ds = this->desc();

         this->change(leaveIdx, none, 0);

         if(leaveStat == SPxBasisBase<Real>::Desc::P_ON_UPPER)
         {
            if(leaveId.isSPxRowId())
            {
               ds.rowStatus(leaveNum) = SPxBasisBase<Real>::Desc::P_ON_LOWER;
               (*theCoPrhs)[leaveIdx] = theLRbound[leaveNum];
            }
            else
            {
               ds.colStatus(leaveNum) = SPxBasisBase<Real>::Desc::P_ON_LOWER;
               (*theCoPrhs)[leaveIdx] = theLCbound[leaveNum];
            }

            theUBbound[leaveIdx] = 0;
            theLBbound[leaveIdx] = -infinity;
         }
         else
         {
            assert(leaveStat == SPxBasisBase<Real>::Desc::P_ON_LOWER);

            if(leaveId.isSPxRowId())
            {
               ds.rowStatus(leaveNum) = SPxBasisBase<Real>::Desc::P_ON_UPPER;
               (*theCoPrhs)[leaveIdx] = theURbound[leaveNum];
            }
            else
            {
               ds.colStatus(leaveNum) = SPxBasisBase<Real>::Desc::P_ON_UPPER;
               (*theCoPrhs)[leaveIdx] = theUCbound[leaveNum];
            }

            theUBbound[leaveIdx] = infinity;
            theLBbound[leaveIdx] = 0;
         }

         // update copricing vector
         theCoPvec->value() = enterVal;
         thePvec->value() = enterVal;

         if(enterVal > epsilon() || enterVal < -epsilon())
            doPupdate();

         // update feasibility vectors
         theFvec->value() = 0;
         assert(theCoTest[leaveIdx] < 0.0);
         m_pricingViol += theCoTest[leaveIdx];
         theCoTest[leaveIdx] *= -1;
      }

      if((leaveMax > entertol() && enterVal <= entertol()) || (leaveMax < -entertol()
            && enterVal >= -entertol()))
      {
         if((theUBbound[leaveIdx] < infinity || theLBbound[leaveIdx] > -infinity)
               && leaveStat != SPxBasisBase<Real>::Desc::P_FREE
               && leaveStat != SPxBasisBase<Real>::Desc::D_FREE)
         {
            m_numCycle++;
            leaveCycles++;
         }
      }
      else
         m_numCycle /= 2;

#ifdef ENABLE_ADDITIONAL_CHECKS
      {
         DVector tmp = fVec();
         multBaseWith(tmp);
         tmp -= fRhs();

         if(tmp.length() > entertol())
         {
            // This happens very frequently and does usually not hurt, so print
            // these warnings only with verbose level INFO2 and higher.
            MSG_INFO2((*spxout), (*spxout) << "WLEAVE64\t" << basis().iteration()
                      << ": fVec error = " << tmp.length() << std::endl;)
            SPxBasisBase<Real>::solve(tmp, fRhs());
            tmp -= fVec();
            MSG_INFO2((*spxout), (*spxout) << "WLEAVE65\t(" << tmp.length() << ")\n";)
         }
      }
#endif  // ENABLE_ADDITIONAL_CHECKS

      return true;
   }
}
} // namespace soplex
