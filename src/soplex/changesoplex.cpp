/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the class library                   */
/*       SoPlex --- the Sequential object-oriented simPlex.                  */
/*                                                                           */
/*    Copyright (C) 1996-2019 Konrad-Zuse-Zentrum                            */
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

#include "soplex/spxdefines.h"
#include "soplex/spxsolver.h"
#include "soplex/spxpricer.h"
#include "soplex/spxratiotester.h"
#include "soplex/exceptions.h"

namespace soplex
{

template <>
void SPxSolverBase<Real>::addedRows(int n)
{

   if(n > 0)
   {
      SPxLP::addedRows(n);

      unInit();
      reDim();

      if(SPxBasisBase<Real>::status() > SPxBasisBase<Real>::NO_PROBLEM)
         SPxBasisBase<Real>::addedRows(n);
   }

   /* we must not assert consistency here, since addedCols() might be still necessary to obtain a consistent basis */
}

template <>
void SPxSolverBase<Real>::addedCols(int n)
{

   if(n > 0)
   {
      SPxLP::addedCols(n);

      unInit();
      reDim();

      if(SPxBasisBase<Real>::status() > SPxBasisBase<Real>::NO_PROBLEM)
         SPxBasisBase<Real>::addedCols(n);
   }

   /* we must not assert consistency here, since addedRows() might be still necessary to obtain a consistent basis */
}

template <>
void SPxSolverBase<Real>::doRemoveRow(int i)
{

   SPxLP::doRemoveRow(i);

   unInit();

   if(SPxBasisBase<Real>::status() > SPxBasisBase<Real>::NO_PROBLEM)
   {
      removedRow(i);

      switch(SPxBasisBase<Real>::status())
      {
      case SPxBasisBase<Real>::DUAL:
      case SPxBasisBase<Real>::INFEASIBLE:
         setBasisStatus(SPxBasisBase<Real>::REGULAR);
         break;

      case SPxBasisBase<Real>::OPTIMAL:
         setBasisStatus(SPxBasisBase<Real>::PRIMAL);
         break;

      default:
         break;
      }
   }
}

template <>
void SPxSolverBase<Real>::doRemoveRows(int perm[])
{

   SPxLP::doRemoveRows(perm);

   unInit();

   if(SPxBasisBase<Real>::status() > SPxBasisBase<Real>::NO_PROBLEM)
   {
      removedRows(perm);

      switch(SPxBasisBase<Real>::status())
      {
      case SPxBasisBase<Real>::DUAL:
      case SPxBasisBase<Real>::INFEASIBLE:
         setBasisStatus(SPxBasisBase<Real>::REGULAR);
         break;

      case SPxBasisBase<Real>::OPTIMAL:
         setBasisStatus(SPxBasisBase<Real>::PRIMAL);
         break;

      default:
         break;
      }
   }
}

template <>
void SPxSolverBase<Real>::doRemoveCol(int i)
{
   forceRecompNonbasicValue();

   SPxLP::doRemoveCol(i);

   unInit();

   if(SPxBasisBase<Real>::status() > SPxBasisBase<Real>::NO_PROBLEM)
   {
      removedCol(i);

      switch(SPxBasisBase<Real>::status())
      {
      case SPxBasisBase<Real>::PRIMAL:
      case SPxBasisBase<Real>::UNBOUNDED:
         setBasisStatus(SPxBasisBase<Real>::REGULAR);
         break;

      case SPxBasisBase<Real>::OPTIMAL:
         setBasisStatus(SPxBasisBase<Real>::DUAL);
         break;

      default:
         break;
      }
   }
}

template <>
void SPxSolverBase<Real>::doRemoveCols(int perm[])
{
   forceRecompNonbasicValue();

   SPxLP::doRemoveCols(perm);

   unInit();

   if(SPxBasisBase<Real>::status() > SPxBasisBase<Real>::NO_PROBLEM)
   {
      removedCols(perm);

      switch(SPxBasisBase<Real>::status())
      {
      case SPxBasisBase<Real>::PRIMAL:
      case SPxBasisBase<Real>::UNBOUNDED:
         setBasisStatus(SPxBasisBase<Real>::REGULAR);
         break;

      case SPxBasisBase<Real>::OPTIMAL:
         setBasisStatus(SPxBasisBase<Real>::DUAL);
         break;

      default:
         break;
      }
   }
}

template <>
void SPxSolverBase<Real>::changeObj(const Vector& newObj, bool scale)
{
   forceRecompNonbasicValue();

   SPxLP::changeObj(newObj, scale);

   /**@todo Factorization remains valid, we do not need a reDim()
    *       pricing vectors should be recomputed.
    */
   unInit();
}

template <>
void SPxSolverBase<Real>::changeObj(int i, const Real& newVal, bool scale)
{
   forceRecompNonbasicValue();

   SPxLP::changeObj(i, newVal, scale);


   /**@todo Factorization remains valid, we do not need a reDim()
    *       pricing vectors should be recomputed.
    */
   unInit();
}

template <>
void SPxSolverBase<Real>::changeMaxObj(const Vector& newObj, bool scale)
{
   forceRecompNonbasicValue();

   SPxLP::changeMaxObj(newObj, scale);

   /**@todo Factorization remains valid, we do not need a reDim()
    *       pricing vectors should be recomputed.
    */
   unInit();
}

template <>
void SPxSolverBase<Real>::changeMaxObj(int i, const Real& newVal, bool scale)
{
   forceRecompNonbasicValue();

   SPxLP::changeMaxObj(i, newVal, scale);

   /**@todo Factorization remains valid, we do not need a reDim()
    *       pricing vectors should be recomputed.
    */
   unInit();
}

template <>
void SPxSolverBase<Real>::changeRowObj(const Vector& newObj, bool scale)
{
   forceRecompNonbasicValue();

   SPxLP::changeRowObj(newObj, scale);

   /**@todo Factorization remains valid, we do not need a reDim()
    *       pricing vectors should be recomputed.
    */
   unInit();
}

template <>
void SPxSolverBase<Real>::changeRowObj(int i, const Real& newVal, bool scale)
{
   forceRecompNonbasicValue();

   SPxLP::changeRowObj(i, newVal, scale);

   /**@todo Factorization remains valid, we do not need a reDim()
    *       pricing vectors should be recomputed.
    */
   unInit();
}

template <>
void SPxSolverBase<Real>::changeLowerStatus(int i, Real newLower, Real oldLower)
{
   SPxBasisBase<Real>::Desc::Status& stat      = desc().colStatus(i);
   Real                    currUpper = upper(i);
   Real                    objChange = 0.0;

   MSG_DEBUG(std::cout << "DCHANG01 changeLowerStatus(): col " << i
             << "[" << newLower << ":" << currUpper << "] " << stat;)

   switch(stat)
   {
   case SPxBasisBase<Real>::Desc::P_ON_LOWER:
      if(newLower <= -infinity)
      {
         if(currUpper >= infinity)
         {
            stat = SPxBasisBase<Real>::Desc::P_FREE;

            if(m_nonbasicValueUpToDate && rep() == COLUMN)
               objChange = -theLCbound[i] * oldLower;
         }
         else
         {
            stat = SPxBasisBase<Real>::Desc::P_ON_UPPER;

            if(m_nonbasicValueUpToDate && rep() == COLUMN)
               objChange = (theUCbound[i] * currUpper) - (theLCbound[i] * oldLower);
         }
      }
      else if(EQ(newLower, currUpper))
      {
         stat = SPxBasisBase<Real>::Desc::P_FIXED;

         if(m_nonbasicValueUpToDate && rep() == COLUMN)
            objChange = maxObj(i) * (newLower - oldLower);
      }
      else if(m_nonbasicValueUpToDate && rep() == COLUMN)
         objChange = theLCbound[i] * (newLower - oldLower);

      break;

   case SPxBasisBase<Real>::Desc::P_ON_UPPER:
      if(EQ(newLower, currUpper))
         stat = SPxBasisBase<Real>::Desc::P_FIXED;

      break;

   case SPxBasisBase<Real>::Desc::P_FREE:
      if(newLower > -infinity)
      {
         stat = SPxBasisBase<Real>::Desc::P_ON_LOWER;

         if(m_nonbasicValueUpToDate && rep() == COLUMN)
            objChange = theLCbound[i] * newLower;
      }

      break;

   case SPxBasisBase<Real>::Desc::P_FIXED:
      if(NE(newLower, currUpper))
      {
         stat = SPxBasisBase<Real>::Desc::P_ON_UPPER;

         if(isInitialized())
            theUCbound[i] = maxObj(i);
      }

      break;

   case SPxBasisBase<Real>::Desc::D_FREE:
   case SPxBasisBase<Real>::Desc::D_ON_UPPER:
   case SPxBasisBase<Real>::Desc::D_ON_LOWER:
   case SPxBasisBase<Real>::Desc::D_ON_BOTH:
   case SPxBasisBase<Real>::Desc::D_UNDEFINED:
      if(rep() == ROW && theShift > 0.0)
         forceRecompNonbasicValue();

      stat = dualColStatus(i);
      break;

   default:
      throw SPxInternalCodeException("XCHANG01 This should never happen.");
   }

   MSG_DEBUG(std::cout << " -> " << stat << std::endl;)

   // we only need to update the nonbasic value in column representation (see nonbasicValue() for comparison/explanation)
   if(rep() == COLUMN)
      updateNonbasicValue(objChange);
}

template <>
void SPxSolverBase<Real>::changeLower(const Vector& newLower, bool scale)
{
   // we better recompute the nonbasic value when changing all lower bounds
   forceRecompNonbasicValue();

   SPxLP::changeLower(newLower, scale);

   if(SPxBasisBase<Real>::status() > SPxBasisBase<Real>::NO_PROBLEM)
   {
      for(int i = 0; i < newLower.dim(); ++i)
         changeLowerStatus(i, lower(i));

      unInit();
   }
}

template <>
void SPxSolverBase<Real>::changeLower(int i, const Real& newLower, bool scale)
{
   if(newLower != (scale ? lowerUnscaled(i) : lower(i)))
   {
      Real oldLower = lower(i);
      // This has to be done before calling changeLowerStatus() because that is calling
      // basis.dualColStatus() which calls lower() and needs the changed value.
      SPxLP::changeLower(i, newLower, scale);

      if(SPxBasisBase<Real>::status() > SPxBasisBase<Real>::NO_PROBLEM)
      {
         changeLowerStatus(i, lower(i), oldLower);
         unInit();
      }
   }
}

template <>
void SPxSolverBase<Real>::changeUpperStatus(int i, Real newUpper, Real oldUpper)
{
   SPxBasisBase<Real>::Desc::Status& stat      = desc().colStatus(i);
   Real                    currLower = lower(i);
   Real                    objChange = 0.0;

   MSG_DEBUG(std::cout << "DCHANG02 changeUpperStatus(): col " << i
             << "[" << currLower << ":" << newUpper << "] " << stat;)

   switch(stat)
   {
   case SPxBasisBase<Real>::Desc::P_ON_LOWER:
      if(newUpper == currLower)
         stat = SPxBasisBase<Real>::Desc::P_FIXED;

      break;

   case SPxBasisBase<Real>::Desc::P_ON_UPPER:
      if(newUpper >= infinity)
      {
         if(currLower <= -infinity)
         {
            stat = SPxBasisBase<Real>::Desc::P_FREE;

            if(m_nonbasicValueUpToDate && rep() == COLUMN)
               objChange = -theUCbound[i] * oldUpper;
         }
         else
         {
            stat = SPxBasisBase<Real>::Desc::P_ON_LOWER;

            if(m_nonbasicValueUpToDate && rep() == COLUMN)
               objChange = (theLCbound[i] * currLower) - (theUCbound[i] * oldUpper);
         }
      }
      else if(EQ(newUpper, currLower))
      {
         stat = SPxBasisBase<Real>::Desc::P_FIXED;

         if(m_nonbasicValueUpToDate && rep() == COLUMN)
            objChange = maxObj(i) * (newUpper - oldUpper);
      }
      else if(m_nonbasicValueUpToDate && rep() == COLUMN)
         objChange = theUCbound[i] * (newUpper - oldUpper);

      break;

   case SPxBasisBase<Real>::Desc::P_FREE:
      if(newUpper < infinity)
      {
         stat = SPxBasisBase<Real>::Desc::P_ON_UPPER;

         if(m_nonbasicValueUpToDate && rep() == COLUMN)
            objChange = theUCbound[i] * newUpper;
      }

      break;

   case SPxBasisBase<Real>::Desc::P_FIXED:
      if(NE(newUpper, currLower))
      {
         stat = SPxBasisBase<Real>::Desc::P_ON_LOWER;

         if(isInitialized())
            theLCbound[i] = maxObj(i);
      }

      break;

   case SPxBasisBase<Real>::Desc::D_FREE:
   case SPxBasisBase<Real>::Desc::D_ON_UPPER:
   case SPxBasisBase<Real>::Desc::D_ON_LOWER:
   case SPxBasisBase<Real>::Desc::D_ON_BOTH:
   case SPxBasisBase<Real>::Desc::D_UNDEFINED:
      if(rep() == ROW && theShift > 0.0)
         forceRecompNonbasicValue();

      stat = dualColStatus(i);
      break;

   default:
      throw SPxInternalCodeException("XCHANG02 This should never happen.");
   }

   MSG_DEBUG(std::cout << " -> " << stat << std::endl;);

   // we only need to update the nonbasic value in column representation (see nonbasicValue() for comparison/explanation)
   if(rep() == COLUMN)
      updateNonbasicValue(objChange);
}

template <>
void SPxSolverBase<Real>::changeUpper(const Vector& newUpper, bool scale)
{
   // we better recompute the nonbasic value when changing all upper bounds
   forceRecompNonbasicValue();

   SPxLP::changeUpper(newUpper, scale);

   if(SPxBasisBase<Real>::status() > SPxBasisBase<Real>::NO_PROBLEM)
   {
      for(int i = 0; i < newUpper.dim(); ++i)
         changeUpperStatus(i, upper(i));

      unInit();
   }
}

template <>
void SPxSolverBase<Real>::changeUpper(int i, const Real& newUpper, bool scale)
{
   if(newUpper != (scale ? upperUnscaled(i) : upper(i)))
   {
      Real oldUpper = upper(i);
      SPxLP::changeUpper(i, newUpper, scale);

      if(SPxBasisBase<Real>::status() > SPxBasisBase<Real>::NO_PROBLEM)
      {
         changeUpperStatus(i, upper(i), oldUpper);
         unInit();
      }
   }
}

template <>
void SPxSolverBase<Real>::changeBounds(const Vector& newLower, const Vector& newUpper, bool scale)
{
   changeLower(newLower, scale);
   changeUpper(newUpper, scale);
}

template <>
void SPxSolverBase<Real>::changeBounds(int i, const Real& newLower, const Real& newUpper,
                                       bool scale)
{
   changeLower(i, newLower, scale);
   changeUpper(i, newUpper, scale);
}

template <>
void SPxSolverBase<Real>::changeLhsStatus(int i, Real newLhs, Real oldLhs)
{
   SPxBasisBase<Real>::Desc::Status& stat      = desc().rowStatus(i);
   Real                    currRhs   = rhs(i);
   Real                    objChange = 0.0;

   MSG_DEBUG(std::cout << "DCHANG03 changeLhsStatus()  : row " << i
             << ": " << stat;)

   switch(stat)
   {
   case SPxBasisBase<Real>::Desc::P_ON_LOWER:
      if(newLhs <= -infinity)
      {
         if(currRhs >= infinity)
         {
            stat = SPxBasisBase<Real>::Desc::P_FREE;

            if(m_nonbasicValueUpToDate && rep() == COLUMN)
               objChange = -theURbound[i] * oldLhs;
         }
         else
         {
            stat = SPxBasisBase<Real>::Desc::P_ON_UPPER;

            if(m_nonbasicValueUpToDate && rep() == COLUMN)
               objChange = (theLRbound[i] * currRhs) - (theURbound[i] * oldLhs);
         }
      }
      else if(EQ(newLhs, currRhs))
      {
         stat = SPxBasisBase<Real>::Desc::P_FIXED;

         if(m_nonbasicValueUpToDate && rep() == COLUMN)
            objChange = maxRowObj(i) * (newLhs - oldLhs);
      }
      else if(m_nonbasicValueUpToDate && rep() == COLUMN)
         objChange = theURbound[i] * (newLhs - oldLhs);

      break;

   case SPxBasisBase<Real>::Desc::P_ON_UPPER:
      if(EQ(newLhs, currRhs))
         stat = SPxBasisBase<Real>::Desc::P_FIXED;

      break;

   case SPxBasisBase<Real>::Desc::P_FREE:
      if(newLhs > -infinity)
      {
         stat = SPxBasisBase<Real>::Desc::P_ON_LOWER;

         if(m_nonbasicValueUpToDate && rep() == COLUMN)
            objChange = theURbound[i] * newLhs;
      }

      break;

   case SPxBasisBase<Real>::Desc::P_FIXED:
      if(NE(newLhs, currRhs))
      {
         stat = SPxBasisBase<Real>::Desc::P_ON_UPPER;

         if(isInitialized())
            theLRbound[i] = maxRowObj(i);
      }

      break;

   case SPxBasisBase<Real>::Desc::D_FREE:
   case SPxBasisBase<Real>::Desc::D_ON_UPPER:
   case SPxBasisBase<Real>::Desc::D_ON_LOWER:
   case SPxBasisBase<Real>::Desc::D_ON_BOTH:
   case SPxBasisBase<Real>::Desc::D_UNDEFINED:
      if(rep() == ROW && theShift > 0.0)
         forceRecompNonbasicValue();

      stat = dualRowStatus(i);
      break;

   default:
      throw SPxInternalCodeException("XCHANG03 This should never happen.");
   }

   MSG_DEBUG(std::cout << " -> " << stat << std::endl;)

   // we only need to update the nonbasic value in column representation (see nonbasicValue() for comparison/explanation)
   if(rep() == COLUMN)
      updateNonbasicValue(objChange);
}

template <>
void SPxSolverBase<Real>::changeLhs(const Vector& newLhs, bool scale)
{
   // we better recompute the nonbasic value when changing all lhs
   forceRecompNonbasicValue();

   SPxLP::changeLhs(newLhs, scale);

   if(SPxBasisBase<Real>::status() > SPxBasisBase<Real>::NO_PROBLEM)
   {
      for(int i = 0; i < nRows(); ++i)
         changeLhsStatus(i, lhs(i));

      unInit();
   }
}

template <>
void SPxSolverBase<Real>::changeLhs(int i, const Real& newLhs, bool scale)
{
   if(newLhs != (scale ? lhsUnscaled(i) : lhs(i)))
   {
      Real oldLhs = lhs(i);
      SPxLP::changeLhs(i, newLhs, scale);

      if(SPxBasisBase<Real>::status() > SPxBasisBase<Real>::NO_PROBLEM)
      {
         changeLhsStatus(i, lhs(i), oldLhs);
         unInit();
      }
   }
}

template <>
void SPxSolverBase<Real>::changeRhsStatus(int i, Real newRhs, Real oldRhs)
{
   SPxBasisBase<Real>::Desc::Status& stat      = desc().rowStatus(i);
   Real                    currLhs   = lhs(i);
   Real                    objChange = 0.0;

   MSG_DEBUG(std::cout << "DCHANG04 changeRhsStatus()  : row " << i
             << ": " << stat;)

   switch(stat)
   {
   case SPxBasisBase<Real>::Desc::P_ON_UPPER:
      if(newRhs >= infinity)
      {
         if(currLhs <= -infinity)
         {
            stat = SPxBasisBase<Real>::Desc::P_FREE;

            if(m_nonbasicValueUpToDate && rep() == COLUMN)
               objChange = -theLRbound[i] * oldRhs;
         }
         else
         {
            stat = SPxBasisBase<Real>::Desc::P_ON_LOWER;

            if(m_nonbasicValueUpToDate && rep() == COLUMN)
               objChange = (theURbound[i] * currLhs) - (theLRbound[i] * oldRhs);
         }
      }
      else if(EQ(newRhs, currLhs))
      {
         stat = SPxBasisBase<Real>::Desc::P_FIXED;

         if(m_nonbasicValueUpToDate && rep() == COLUMN)
            objChange = maxRowObj(i) * (newRhs - oldRhs);
      }
      else if(m_nonbasicValueUpToDate && rep() == COLUMN)
         objChange = theLRbound[i] * (newRhs - oldRhs);

      break;

   case SPxBasisBase<Real>::Desc::P_ON_LOWER:
      if(EQ(newRhs, currLhs))
         stat = SPxBasisBase<Real>::Desc::P_FIXED;

      break;

   case SPxBasisBase<Real>::Desc::P_FREE:
      if(newRhs < infinity)
      {
         stat = SPxBasisBase<Real>::Desc::P_ON_UPPER;

         if(m_nonbasicValueUpToDate && rep() == COLUMN)
            objChange = theLRbound[i] * newRhs;
      }

      break;

   case SPxBasisBase<Real>::Desc::P_FIXED:
      if(NE(newRhs, currLhs))
      {
         stat = SPxBasisBase<Real>::Desc::P_ON_LOWER;

         if(isInitialized())
            theURbound[i] = maxRowObj(i);
      }

      break;

   case SPxBasisBase<Real>::Desc::D_FREE:
   case SPxBasisBase<Real>::Desc::D_ON_UPPER:
   case SPxBasisBase<Real>::Desc::D_ON_LOWER:
   case SPxBasisBase<Real>::Desc::D_ON_BOTH:
   case SPxBasisBase<Real>::Desc::D_UNDEFINED:
      if(rep() == ROW && theShift > 0.0)
         forceRecompNonbasicValue();

      stat = dualRowStatus(i);
      break;

   default:
      throw SPxInternalCodeException("XCHANG04 This should never happen.");
   }

   MSG_DEBUG(std::cout << " -> " << stat << std::endl;)

   // we only need to update the nonbasic value in column representation (see nonbasicValue() for comparison/explanation)
   if(rep() == COLUMN)
      updateNonbasicValue(objChange);
}


template <>
void SPxSolverBase<Real>::changeRhs(const Vector& newRhs, bool scale)
{
   // we better recompute the nonbasic value when changing all rhs
   forceRecompNonbasicValue();

   SPxLP::changeRhs(newRhs, scale);

   if(SPxBasisBase<Real>::status() > SPxBasisBase<Real>::NO_PROBLEM)
   {
      for(int i = 0; i < nRows(); ++i)
         changeRhsStatus(i, rhs(i));

      unInit();
   }
}

template <>
void SPxSolverBase<Real>::changeRhs(int i, const Real& newRhs, bool scale)
{
   if(newRhs != (scale ? rhsUnscaled(i) : rhs(i)))
   {
      Real oldRhs = rhs(i);
      SPxLP::changeRhs(i, newRhs, scale);

      if(SPxBasisBase<Real>::status() > SPxBasisBase<Real>::NO_PROBLEM)
      {
         changeRhsStatus(i, rhs(i), oldRhs);
         unInit();
      }
   }
}

template <>
void SPxSolverBase<Real>::changeRange(const Vector& newLhs, const Vector& newRhs, bool scale)
{
   // we better recompute the nonbasic value when changing all ranges
   forceRecompNonbasicValue();

   SPxLP::changeLhs(newLhs, scale);
   SPxLP::changeRhs(newRhs, scale);

   if(SPxBasisBase<Real>::status() > SPxBasisBase<Real>::NO_PROBLEM)
   {
      for(int i = nRows() - 1; i >= 0; --i)
      {
         changeLhsStatus(i, lhs(i));
         changeRhsStatus(i, rhs(i));
      }

      unInit();
   }
}

template <>
void SPxSolverBase<Real>::changeRange(int i, const Real& newLhs, const Real& newRhs, bool scale)
{
   Real oldLhs = lhs(i);
   Real oldRhs = rhs(i);

   SPxLP::changeLhs(i, newLhs, scale);
   SPxLP::changeRhs(i, newRhs, scale);

   if(SPxBasisBase<Real>::status() > SPxBasisBase<Real>::NO_PROBLEM)
   {
      changeLhsStatus(i, lhs(i), oldLhs);
      changeRhsStatus(i, rhs(i), oldRhs);
      unInit();
   }
}

template <>
void SPxSolverBase<Real>::changeRow(int i, const LPRow& newRow, bool scale)
{
   forceRecompNonbasicValue();

   SPxLP::changeRow(i, newRow, scale);

   if(SPxBasisBase<Real>::status() > SPxBasisBase<Real>::NO_PROBLEM)
      SPxBasisBase<Real>::changedRow(i);

   unInit();
}

template <>
void SPxSolverBase<Real>::changeCol(int i, const LPCol& newCol, bool scale)
{
   if(i < 0)
      return;

   forceRecompNonbasicValue();

   SPxLP::changeCol(i, newCol, scale);

   if(SPxBasisBase<Real>::status() > SPxBasisBase<Real>::NO_PROBLEM)
      SPxBasisBase<Real>::changedCol(i);

   unInit();
}

template <>
void SPxSolverBase<Real>::changeElement(int i, int j, const Real& val, bool scale)
{
   if(i < 0 || j < 0)
      return;

   forceRecompNonbasicValue();

   SPxLP::changeElement(i, j, val, scale);

   if(SPxBasisBase<Real>::status() > SPxBasisBase<Real>::NO_PROBLEM)
      SPxBasisBase<Real>::changedElement(i, j);

   unInit();
}

template <>
void SPxSolverBase<Real>::changeSense(SPxSense sns)
{

   SPxLP::changeSense(sns);
   unInit();
}
} // namespace soplex
