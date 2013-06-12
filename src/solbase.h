/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the class library                   */
/*       SoPlex --- the Sequential object-oriented simPlex.                  */
/*                                                                           */
/*    Copyright (C) 1996      Roland Wunderling                              */
/*                  1996-2013 Konrad-Zuse-Zentrum                            */
/*                            fuer Informationstechnik Berlin                */
/*                                                                           */
/*  SoPlex is distributed under the terms of the ZIB Academic Licence.       */
/*                                                                           */
/*  You should have received a copy of the ZIB Academic License              */
/*  along with SoPlex; see the file COPYING. If not email to soplex@zib.de.  */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/**@file  solbase.h
 * @brief Class for storing a primal-dual solution with basis information
 */
#ifndef _SOLBASE_H_
#define _SOLBASE_H_

#include <assert.h>
#include <string.h>
#include <math.h>
#include <iostream>

#include "basevectors.h"
#include "spxsolver.h" // needed for basis information

namespace soplex
{
/**@class   SolBase
 * @brief   Class for storing a primal-dual solution with basis information
 * @ingroup Algo
 */
template< class R >
class SolBase
{
   friend class SoPlex2;

public:
   /// is a primal feasible solution available?
   bool hasPrimal() const
   {
      return _hasPrimal;
   }

   /// gets the primal solution vector if available; returns true on success
   bool getPrimal(VectorBase<R>& vector) const
   {
      if( _hasPrimal )
         vector = _primal;

      return _hasPrimal;
   }

   /// gets the vector of slack values if available; returns true on success
   bool getSlacks(VectorBase<R>& vector) const
   {
      if( _hasPrimal )
         vector = _slacks;

      return _hasPrimal;
   }

   /// is a primal unbounded ray available?
   bool hasPrimalray() const
   {
      return _hasPrimalray;
   }

   /// gets the primal unbounded ray if available; returns true on success
   bool getPrimalray(VectorBase<R>& vector) const
   {
      if( _hasPrimalray )
         vector = _primalray;

      return _hasPrimalray;
   }

   /// is a dual solution available?
   bool hasDual() const
   {
      return _hasDual;
   }

   /// gets the dual solution vector if available; returns true on success
   bool getDual(VectorBase<R>& vector) const
   {
      if( _hasDual )
         vector = _dual;

      return _hasDual;
   }

   /// gets the vector of reduced cost values if available; returns true on success
   bool getRedcost(VectorBase<R>& vector) const
   {
      if( _hasDual )
         vector = _redcost;

      return _hasDual;
   }

   /// is a dual farkas ray available?
   bool hasDualfarkas() const
   {
      return _hasDualfarkas;
   }

   /// gets the Farkas proof if available; returns true on success
   bool getDualfarkas(VectorBase<R>& vector) const
   {
      if( _hasDualfarkas )
         vector = _dualfarkas;

      return _hasDualfarkas;
   }

   /// is basis information available?
   bool hasBasis() const
   {
      return _hasBasis;
   }

   /// gets basis information if available; returns true on success
   bool getBasis(DataArray<SPxSolver::VarStatus> rows, DataArray<SPxSolver::VarStatus> cols) const
   {
      if( _hasBasis )
      {
         rows = _basisStatusRows;
         cols = _basisStatusCols;
      }

      return _hasBasis;
   }

private:
   DVectorBase<R> _primal;
   DVectorBase<R> _slacks;
   DVectorBase<R> _primalray;
   DVectorBase<R> _dual;
   DVectorBase<R> _redcost;
   DVectorBase<R> _dualfarkas;

   DataArray< SPxSolver::VarStatus > _basisStatusRows;
   DataArray< SPxSolver::VarStatus > _basisStatusCols;

   unsigned int _hasPrimal:1;
   unsigned int _hasPrimalray:1;
   unsigned int _hasDual:1;
   unsigned int _hasDualfarkas:1;
   unsigned int _hasBasis:1;

   /// default constructor only for friends
   SolBase<R>()
   {
      _invalidate();
   }

   /// invalidate solution
   void _invalidate()
   {
      _hasPrimal = false;
      _hasPrimalray = false;
      _hasDual = false;
      _hasDualfarkas = false;
      _hasBasis = false;
   }
};
} // namespace soplex
#endif // _SOLBASE_H_

// ---------------------------------------------------------------------------------------------------------------------
// Emacs Local Variables:
// Emacs mode:c++
// Emacs c-basic-offset:3
// Emacs tab-width:8
// Emacs indent-tabs-mode:nil
// Emacs End:
// ---------------------------------------------------------------------------------------------------------------------
