/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the class library                   */
/*       SoPlex --- the Sequential object-oriented simPlex.                  */
/*                                                                           */
/*    Copyright (C) 1996      Roland Wunderling                              */
/*                  1996-2014 Konrad-Zuse-Zentrum                            */
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
   friend class SoPlex;
   template < class S > friend class SolBase;

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
   bool hasPrimalRay() const
   {
      return _hasPrimalRay;
   }

   /// gets the primal unbounded ray if available; returns true on success
   bool getPrimalRay(VectorBase<R>& vector) const
   {
      if( _hasPrimalRay )
         vector = _primalRay;

      return _hasPrimalRay;
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
   bool getRedCost(VectorBase<R>& vector) const
   {
      if( _hasDual )
         vector = _redCost;

      return _hasDual;
   }

   /// is a dual farkas ray available?
   bool hasDualFarkas() const
   {
      return _hasDualFarkas;
   }

   /// gets the Farkas proof if available; returns true on success
   bool getDualFarkas(VectorBase<R>& vector) const
   {
      if( _hasDualFarkas )
         vector = _dualFarkas;

      return _hasDualFarkas;
   }

   /// invalidate solution
   void invalidate()
   {
      _hasPrimal = false;
      _hasPrimalRay = false;
      _hasDual = false;
      _hasDualFarkas = false;
   }

private:
   DVectorBase<R> _primal;
   DVectorBase<R> _slacks;
   DVectorBase<R> _primalRay;
   DVectorBase<R> _dual;
   DVectorBase<R> _redCost;
   DVectorBase<R> _dualFarkas;

   R _primalObjVal;
   R _dualObjVal;

   unsigned int _hasPrimal:1;
   unsigned int _hasPrimalRay:1;
   unsigned int _hasDual:1;
   unsigned int _hasDualFarkas:1;

   /// default constructor only for friends
   SolBase<R>()
   {
      invalidate();
   }

   /// assignment operator only for friends
   SolBase<R>& operator=(const SolBase<R>& sol)
   {
      if( this != &sol )
      {

         _hasPrimal = sol._hasPrimal;
         if( _hasPrimal )
         {
            _primal = sol._primal;
            _slacks = sol._slacks;
            _primalObjVal = sol._primalObjVal;
         }

         _hasPrimalRay = sol._hasPrimalRay;
         if( _hasPrimalRay )
            _primalRay = sol._primalRay;

         _hasDual = sol._hasDual;
         if( _hasDual )
         {
            _dual = sol._dual;
            _redCost = sol._redCost;
            _dualObjVal = sol._dualObjVal;
         }

         _hasDualFarkas = sol._hasDualFarkas;
         if( _hasDualFarkas )
            _dualFarkas = sol._dualFarkas;
      }

      return *this;
   }

   /// assignment operator only for friends
   template < class S >
   SolBase<R>& operator=(const SolBase<S>& sol)
   {
      if( (SolBase<S>*)this != &sol )
      {

         _hasPrimal = sol._hasPrimal;
         if( _hasPrimal )
         {
            _primal = sol._primal;
            _slacks = sol._slacks;
            _primalObjVal = R(sol._primalObjVal);
         }

         _hasPrimalRay = sol._hasPrimalRay;
         if( _hasPrimalRay )
            _primalRay = sol._primalRay;

         _hasDual = sol._hasDual;
         if( _hasDual )
         {
            _dual = sol._dual;
            _redCost = sol._redCost;
            _dualObjVal = R(sol._dualObjVal);
         }

         _hasDualFarkas = sol._hasDualFarkas;
         if( _hasDualFarkas )
            _dualFarkas = sol._dualFarkas;
      }

      return *this;
   }
};
} // namespace soplex
#endif // _SOLBASE_H_
