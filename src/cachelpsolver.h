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
#pragma ident "@(#) $Id: cachelpsolver.h,v 1.4 2001/11/15 22:35:13 bzfkocht Exp $"

/**@file  cachelpsolver.h
 * @brief LP solver with result caching.
 */
#ifndef _CACHELPSOLVER_H_
#define _CACHELPSOLVER_H_

#include <assert.h>

#include "lpsolver.h"
#include "dvector.h"
#include "dataarray.h"

namespace soplex
{
/**@brief   LP solver with result caching.
   @ingroup Algorithmic 

   Class CacheLPSolver provides some caching for some of the data access
   methods of LPSolver. E.g. when multiple successive calls to method
   primal() occur, only the first call the primal solution vector is really
   constructed. The following calls will return the previously constructed
   vector.
 */
class CacheLPSolver : public LPSolver
{
private:
   mutable DVector theObj;
   mutable DVector thePrimal;
   mutable DVector theRdCost;
   mutable DataArray < signed char > colBase;

   mutable DVector theDual;
   mutable DVector theSlack;
   mutable DataArray < signed char > rowBase;

public:
   /// #i-th value of objective vector.
   virtual double obj(int i) const = 0;

   /// #id-th value of objective vector.
   virtual double obj(ColId id) const = 0;

   /// return const objective vector.
   virtual const Vector& obj() const;

   /// return const solution vector for primal variables if available.
   virtual const Vector& primal() const;

   /// return const vector of slack variables if available.
   virtual const Vector& slacks() const;

   /// return const solution vector for dual variables if available.
   virtual const Vector& dual() const;

   /// return const vector of reduced costs if available.
   virtual const Vector& rdCost() const;

   /// return row ids if available.
   virtual void getRowIds(RowId ids[]) const;

   /// return column ids if available.
   virtual void getColIds(ColId ids[]) const;

   /// return status for primal vars if available.
   virtual const signed char* rowBasis() const;

   /// return status for dual vars if available.
   virtual const signed char* colBasis() const;
};

} // namespace soplex
#endif // _CACHELPSOLVER_H_

//-----------------------------------------------------------------------------
//Emacs Local Variables:
//Emacs mode:c++
//Emacs c-basic-offset:3
//Emacs tab-width:8
//Emacs indent-tabs-mode:nil
//Emacs End:
//-----------------------------------------------------------------------------
