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
#pragma ident "@(#) $Id: spxchangebasis.cpp,v 1.5 2001/11/29 14:00:25 bzfkocht Exp $"


//@ -----------------------------------------------------------------------------
/*      \SubSection{Manipulation Methods}
 */
#include <iostream>
#include <math.h>

#include "spxbasis.h"



#include "soplex.h"

namespace soplex
{





//@ -----------------------------------------------------------------------------
void SPxBasis::reDim()
{
   assert(theLP != 0);

   thedesc.reSize (theLP->nRows(), theLP->nCols());

   if (theLP->dim() != matrix.size())
   {
      matrix.reSize (theLP->dim());
      theBaseId.reSize(theLP->dim());
      matrixIsSetup = false;
      factorized = false;
   }
}

void SPxBasis::addedRows(int n)
{
   assert(theLP != 0);

   reDim();

   if (theLP->rep() == SoPlex::COLUMN)
   {
      for (int i = theLP->nRows() - n; i < theLP->nRows(); ++i)
      {
         thedesc.rowStatus(i) = dualRowStatus(i);
         baseId(i) = theLP->SPxLP::rId(i);
      }
      matrixIsSetup = false;
   }
   else
   {
      assert(theLP->rep() == SoPlex::ROW);
      for (int i = theLP->nRows() - n; i < theLP->nRows(); ++i)
         thedesc.rowStatus(i) = dualRowStatus(i);
   }

   if (status() > NO_PROBLEM && matrixIsSetup)
      loadMatrixVecs();
}

void SPxBasis::removedRow(int i)
{
   assert(status() > NO_PROBLEM);
   assert(theLP != 0);

   if (theLP->rep() == SoPlex::ROW)
   {
      if (theLP->isBasic(thedesc.rowStatus(i)))
      {
         setStatus(NO_PROBLEM);
         factorized = false;
         // std::cout << "Are you sure, you wanna do that?\n";
      }

   }
   else
   {
      assert(theLP->rep() == SoPlex::COLUMN);
      factorized = false;
      if (!theLP->isBasic(thedesc.rowStatus(i)))
      {
         setStatus(NO_PROBLEM);
         // std::cout << "Are you sure, you wanna do that?\n";
      }

      else if (status() > NO_PROBLEM && matrixIsSetup)
      {
         for (int j = theLP->dim(); j >= 0; --j)
         {
            SoPlex::Id id = baseId(j);
            if (id.isSPxRowId()
                 && theLP->number(SPxLP::SPxRowId(id)) < 0)
            {
               baseId(j) = baseId(theLP->dim());
               if (j < theLP->dim())
                  matrix[j] = &theLP->vector(baseId(j));
               break;
            }
         }
      }
   }

   thedesc.rowStatus(i) = thedesc.rowStatus(theLP->nRows());
   reDim();
}

void SPxBasis::removedRows(int perm[])
{
   assert(status() > NO_PROBLEM);
   assert(theLP != 0);

   int i;
   int n = thedesc.nRows();

   if (theLP->rep() == SoPlex::ROW)
   {
      for (i = 0; i < n; ++i)
      {
         if (perm[i] != i)
         {
            if (perm[i] < 0)               // row got removed
            {
               if (theLP->isBasic(thedesc.rowStatus(i)))
               {
                  setStatus(NO_PROBLEM);
                  factorized = false;
                  // std::cout << "Are you sure, you wanna do that?\n";
               }

            }
            else                            // row was moved
               thedesc.rowStatus(perm[i]) = thedesc.rowStatus(i);
         }
      }
   }
   else
   {
      assert(theLP->rep() == SoPlex::COLUMN);
      factorized = matrixIsSetup = false;
      for (i = 0; i < n; ++i)
      {
         if (perm[i] != i)
         {
            if (perm[i] < 0)               // row got removed
            {
               if (!theLP->isBasic(thedesc.rowStatus(i)))
                  setStatus(NO_PROBLEM);
            }
            else                            // row was moved
               thedesc.rowStatus(perm[i]) = thedesc.rowStatus(i);
         }
      }
   }

   reDim();
}


static SPxBasis::Desc::Status
primalColStatus(int i, const SPxLP* theLP)
{
   assert(theLP != 0);

   if (theLP->upper(i) < SPxLP::infinity)
   {
      if (theLP->lower(i) > -SPxLP::infinity)
      {
         if (theLP->lower(i) == theLP->SPxLP::upper(i))
            return SPxBasis::Desc::P_FIXED;
         /*
             else
                 return (-theLP->lower(i) < theLP->upper(i))
                             ? SPxBasis::Desc::P_ON_LOWER
                          : SPxBasis::Desc::P_ON_UPPER;
         */
         else if (theLP->maxObj(i) == 0)
            return (-theLP->lower(i) < theLP->upper(i))
                   ? SPxBasis::Desc::P_ON_LOWER
                : SPxBasis::Desc::P_ON_UPPER;
         else
            return (theLP->maxObj(i) < 0)
                   ? SPxBasis::Desc::P_ON_LOWER
                : SPxBasis::Desc::P_ON_UPPER;
      }
      else
         return SPxBasis::Desc::P_ON_UPPER;
   }
   else if (theLP->lower(i) > -SPxLP::infinity)
      return SPxBasis::Desc::P_ON_LOWER;
   else
      return SPxBasis::Desc::P_FREE;
}


void SPxBasis::addedCols(int n)
{
   assert(theLP != 0);

   reDim();

   if (theLP->rep() == SoPlex::ROW)
   {
      for (int i = theLP->nCols() - n; i < theLP->nCols(); ++i)
      {
         thedesc.colStatus(i) = primalColStatus(i, theLP);
         baseId(i) = theLP->SPxLP::cId(i);
      }
   }
   else
   {
      assert(theLP->rep() == SoPlex::COLUMN);
      for (int i = theLP->nCols() - n; i < theLP->nCols(); ++i)
         thedesc.colStatus(i) = primalColStatus(i, theLP);
   }

   if (status() > NO_PROBLEM && matrixIsSetup)
      loadMatrixVecs();
}

void SPxBasis::removedCol(int i)
{
   assert(status() > NO_PROBLEM);
   assert(theLP != 0);

   if (theLP->rep() == SoPlex::COLUMN)
   {
      if (theLP->isBasic(thedesc.colStatus(i)))
         setStatus(NO_PROBLEM);
   }
   else
   {
      assert(theLP->rep() == SoPlex::ROW);
      factorized = false;
      if (!theLP->isBasic(thedesc.colStatus(i)))
         setStatus(NO_PROBLEM);
      else if (status() > NO_PROBLEM)
      {
         for (int j = theLP->dim(); j >= 0; --j)
         {
            SoPlex::Id id = baseId(j);
            if (id.isSPxColId()
                 && theLP->number(SPxLP::SPxColId(id)) < 0)
            {
               baseId(j) = baseId(theLP->dim());
               if (j < theLP->dim())
                  matrix[j] = &theLP->vector(baseId(j));
               break;
            }
         }
      }
   }

   thedesc.colStatus(i) = thedesc.colStatus(theLP->nCols());
   reDim();
}

void SPxBasis::removedCols(int perm[])
{
   assert(status() > NO_PROBLEM);
   assert(theLP != 0);

   int i;
   int n = thedesc.nCols();

   if (theLP->rep() == SoPlex::COLUMN)
   {
      for (i = 0; i < n; ++i)
      {
         if (perm[i] < 0)           // column got removed
         {
            if (theLP->isBasic(thedesc.colStatus(i)))
               setStatus(NO_PROBLEM);
         }
         else                        // column was potentially moved
            thedesc.colStatus(perm[i]) = thedesc.colStatus(i);
      }
   }
   else
   {
      assert(theLP->rep() == SoPlex::ROW);
      factorized = matrixIsSetup = false;
      for (i = 0; i < n; ++i)
      {
         if (perm[i] != i)
         {
            if (perm[i] < 0)               // column got removed
            {
               if (!theLP->isBasic(thedesc.colStatus(i)))
                  setStatus(NO_PROBLEM);
            }
            else                            // column was moved
               thedesc.colStatus(perm[i]) = thedesc.colStatus(i);
         }
      }
   }

   reDim();
}

/**@todo is this correctly implemented?
 */
void SPxBasis::changedRow(int /*row*/)
{
   factorized = matrixIsSetup = false;
}

/**@todo is this correctly implemented?
 */
void SPxBasis::changedCol(int /*col*/)
{
   factorized = matrixIsSetup = false;
}

/**@todo is this correctly implemented?
 */
void SPxBasis::changedElement(int /*row*/, int /*col*/)
{
   factorized = matrixIsSetup = false;
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
