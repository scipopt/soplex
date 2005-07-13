/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the class library                   */
/*       SoPlex --- the Sequential object-oriented simPlex.                  */
/*                                                                           */
/*    Copyright (C) 1997-1999 Roland Wunderling                              */
/*                  1997-2002 Konrad-Zuse-Zentrum                            */
/*                            fuer Informationstechnik Berlin                */
/*                                                                           */
/*  SoPlex is distributed under the terms of the ZIB Academic Licence.       */
/*                                                                           */
/*  You should have received a copy of the ZIB Academic License              */
/*  along with SoPlex; see the file COPYING. If not email to soplex@zib.de.  */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
#pragma ident "@(#) $Id: spxchangebasis.cpp,v 1.21 2005/07/13 19:05:32 bzforlow Exp $"

//#define DEBUGGING 1

#include <iostream>

#include "spxdefines.h"
#include "spxbasis.h"
#include "spxsolver.h"
#include "spxout.h"

namespace soplex
{

void SPxBasis::reDim()
{
   METHOD( "SPxBasis::reDim()" );

   assert(theLP != 0);

   DEBUG({ s_spxout << "SPxBasis::reDim():"
                    << " matrixIsSetup=" << matrixIsSetup
                    << " fatorized=" << factorized
                    << std::endl; });

   thedesc.reSize (theLP->nRows(), theLP->nCols());

   if (theLP->dim() != matrix.size())
   {
      VERBOSE3({ s_spxout << "basis redimensioning invalidates factorization" << std::endl; });

      matrix.reSize (theLP->dim());
      theBaseId.reSize(theLP->dim());
      matrixIsSetup = false;
      factorized    = false;
   }

   DEBUG({ s_spxout << "SPxBasis::reDim(): -->"
                    << " matrixIsSetup=" << matrixIsSetup
                    << " fatorized=" << factorized
                    << std::endl; });

   assert( matrix.size()    >= theLP->dim() );
   assert( theBaseId.size() >= theLP->dim() );
}

void SPxBasis::addedRows(int n)
{
   METHOD( "SPxBasis::addedRows()" );

   assert(theLP != 0);

   if( n > 0 )
   {
      reDim();
      
      if (theLP->rep() == SPxSolver::COLUMN)
      {
         /* I think, after adding rows in column representation,
            reDim should set these bools to false. */
         assert( !matrixIsSetup && !factorized );

         for (int i = theLP->nRows() - n; i < theLP->nRows(); ++i)
         {
            thedesc.rowStatus(i) = dualRowStatus(i);
            baseId(i) = theLP->SPxLP::rId(i);
         }

         /* ??? I think, this cannot happen. */
         /* if matrix was set up, load new basis vectors to the matrix */
         if (status() > NO_PROBLEM && matrixIsSetup)
            loadMatrixVecs();
      }
      else
      {
         assert(theLP->rep() == SPxSolver::ROW);

         for (int i = theLP->nRows() - n; i < theLP->nRows(); ++i)
            thedesc.rowStatus(i) = dualRowStatus(i);
      }

      /* update basis status */
      switch (status())
      {
      case PRIMAL:
      case UNBOUNDED:
         setStatus(REGULAR);
         break;
      case OPTIMAL:
      case INFEASIBLE:
         setStatus(DUAL);
         break;
      case NO_PROBLEM:
      case SINGULAR:
      case REGULAR:
      case DUAL:
         break;
      default:
         ERROR( s_spxout << "Unknown basis status!" << std::endl; )
         assert(false);
      }
   }
}

void SPxBasis::removedRow(int i)
{
   METHOD( "SPxBasis::removedRow()" );

   assert(status() >  NO_PROBLEM);
   assert(theLP    != 0);

   if (theLP->rep() == SPxSolver::ROW)
   {
      if (theLP->isBasic(thedesc.rowStatus(i)))
      {
         setStatus(NO_PROBLEM);
         factorized = false;

         DEBUG( s_spxout << "Are you sure, you wanna do that?\n"; );
      }
   }
   else
   {
      assert(theLP->rep() == SPxSolver::COLUMN);
      factorized = false;
      if (!theLP->isBasic(thedesc.rowStatus(i)))
      {
         setStatus(NO_PROBLEM);
         DEBUG( s_spxout << "Are you sure, you wanna do that?\n"; );
      }
      else if (status() > NO_PROBLEM && matrixIsSetup)
      {
         for (int j = theLP->dim(); j >= 0; --j)
         {
            SPxId id = baseId(j);

            if (id.isSPxRowId() && theLP->number(SPxRowId(id)) < 0)
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

void SPxBasis::removedRows(const int perm[])
{
   METHOD( "SPxBasis::removedRows()" );
   assert(status() > NO_PROBLEM);
   assert(theLP != 0);

   int i;
   int n = thedesc.nRows();

   if (theLP->rep() == SPxSolver::ROW)
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
                  factorized = matrixIsSetup = false;
                  DEBUG( s_spxout << "Are you sure, you wanna do that?\n"; );
               }
            }
            else                            // row was moved
               thedesc.rowStatus(perm[i]) = thedesc.rowStatus(i);
         }
      }
   }
   else
   {
      assert(theLP->rep() == SPxSolver::COLUMN);

      factorized    = false;
      matrixIsSetup = false;

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

   if (theLP->upper(i) < infinity)
   {
      if (theLP->lower(i) > -infinity)
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
   else if (theLP->lower(i) > -infinity)
      return SPxBasis::Desc::P_ON_LOWER;
   else
      return SPxBasis::Desc::P_FREE;
}


void SPxBasis::addedCols(int n)
{
   METHOD( "SPxBasis::addedCols()" );
   assert(theLP != 0);

   if( n > 0 )
   {
      reDim();
      
      if (theLP->rep() == SPxSolver::ROW)
      {
         /* I think, after adding columns in row representation,
            reDim should set these bools to false. */
         assert( !matrixIsSetup && !factorized );

         for (int i = theLP->nCols() - n; i < theLP->nCols(); ++i)
         {
            thedesc.colStatus(i) = primalColStatus(i, theLP);
            baseId(i) = theLP->SPxLP::cId(i);
         }

         /* ??? I think, this cannot happen. */
         /* if matrix was set up, load new basis vectors to the matrix */
         if (status() > NO_PROBLEM && matrixIsSetup)
            loadMatrixVecs();
      }
      else
      {
         assert(theLP->rep() == SPxSolver::COLUMN);
         for (int i = theLP->nCols() - n; i < theLP->nCols(); ++i)
            thedesc.colStatus(i) = primalColStatus(i, theLP);
      }
         
      switch (status())
      {
      case DUAL:
      case INFEASIBLE:
         setStatus(REGULAR);
         break;
      case OPTIMAL:
      case UNBOUNDED:
         setStatus(PRIMAL);
         break;
      case NO_PROBLEM:
      case SINGULAR:
      case REGULAR:
      case PRIMAL:
         break;
      default:
         ERROR( s_spxout << "Unknown basis status!" << std::endl; )
         assert(false);
      }
   }
}

void SPxBasis::removedCol(int i)
{
   METHOD( "SPxBasis::removedCol()" );
   assert(status() > NO_PROBLEM);
   assert(theLP != 0);

   if (theLP->rep() == SPxSolver::COLUMN)
   {
      if (theLP->isBasic(thedesc.colStatus(i)))
         setStatus(NO_PROBLEM);
   }
   else
   {
      assert(theLP->rep() == SPxSolver::ROW);
      factorized = false;
      if (!theLP->isBasic(thedesc.colStatus(i)))
         setStatus(NO_PROBLEM);
      else if (status() > NO_PROBLEM)
      {
         assert( matrixIsSetup );
         for (int j = theLP->dim(); j >= 0; --j)
         {
            SPxId id = baseId(j);
            if (id.isSPxColId()
                 && theLP->number(SPxColId(id)) < 0)
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

void SPxBasis::removedCols(const int perm[])
{
   METHOD( "SPxBasis::removedCols()" );
   assert(status() > NO_PROBLEM);
   assert(theLP != 0);

   int i;
   int n = thedesc.nCols();

   if (theLP->rep() == SPxSolver::COLUMN)
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
      assert(theLP->rep() == SPxSolver::ROW);
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
void SPxBasis::invalidate()
{
   METHOD( "SPxBasis::invalidate()" );

   VERBOSE3({ s_spxout << "explicit invalidation of factorization" << std::endl; });

   factorized    = false;
   matrixIsSetup = false;
}

/**@todo is this correctly implemented?
 */
void SPxBasis::changedRow(int /*row*/)
{
   METHOD( "SPxBasis::changedRow()" );
   invalidate();
}

/**@todo is this correctly implemented?
 */
void SPxBasis::changedCol(int /*col*/)
{
   METHOD( "SPxBasis::changedCol()" );
   invalidate();
}

/**@todo is this correctly implemented?
 */
void SPxBasis::changedElement(int /*row*/, int /*col*/)
{
   METHOD( "SPxBasis::changedElement()" );
   invalidate();
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

