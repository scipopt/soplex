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
#pragma ident "@(#) $Id: cachelpsolver.cpp,v 1.1 2001/11/06 16:18:31 bzfkocht Exp $"

/*      \Section{Complex Methods}
 */

/*  Import system include files
 */
#include <assert.h>
#include <iostream>


/*  and class header files
 */
#include "cachelpsolver.h"

namespace soplex
{


//@ ----------------------------------------------------------------------------
const Vector& CacheLPSolver::primal() const
{
   if (thePrimal.dim() != nofCols())
      thePrimal.reDim(nofCols());

   getPrimal(thePrimal);

   return thePrimal ;
}

const Vector& CacheLPSolver::slacks() const
{
   if (theSlack.dim() != nofRows())
      theSlack.reDim( nofRows());

   getSlacks(theSlack);

   return theSlack;
}

const Vector& CacheLPSolver::dual() const
{
   if (theDual.dim() != nofRows())
      theDual.reDim(nofRows());

   getDual(theDual);
   
   return theDual;
}

const Vector& CacheLPSolver::rdCost() const
{
   if (theRdCost.dim() != nofCols())
      theRdCost.reDim(nofCols());
  
   getRdCost(theRdCost);

   return theRdCost;
}

const Vector& CacheLPSolver::obj() const
{
   if (theObj.dim() != nofCols())
      theObj.reDim( nofCols());

   getObj(theObj);

   return theObj;
}

void CacheLPSolver::getRowIds(RowId ids[]) const
{
   for (int i = nofRows(); --i >= 0;)
      ids[i] = rowId(i);
}

void CacheLPSolver::getColIds(ColId ids[]) const
{
   for (int i = nofCols(); --i >= 0;)
      ids[i] = colId(i);
}

const signed char* CacheLPSolver::rowBasis() const
{
   if (rowBase.size() != nofRows())
      rowBase.reSize(nofRows());  

   getBasis(rowBase.get_ptr(), 0);

   return rowBase.get_const_ptr();
}

const signed char* CacheLPSolver::colBasis() const
{
   if (colBase.size() != nofCols())
      colBase.reSize(nofCols());

   getBasis(0, colBase.get_ptr());

   return colBase.get_const_ptr();
}
} // namespace soplex

//-----------------------------------------------------------------------------
//Emacs Local Variables:
//Emacs c-basic-offset:3
//Emacs tab-width:8
//Emacs indent-tabs-mode:nil
//Emacs End:
//-----------------------------------------------------------------------------


