/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the class library                   */
/*       SoPlex --- the Sequential object-oriented simPlex.                  */
/*                                                                           */
/*    Copyright (C) 1996-2020 Konrad-Zuse-Zentrum                            */
/*                            fuer Informationstechnik Berlin                */
/*                                                                           */
/*  SoPlex is distributed under the terms of the ZIB Academic Licence.       */
/*                                                                           */
/*  You should have received a copy of the ZIB Academic License              */
/*  along with SoPlex; see the file COPYING. If not email to soplex@zib.de.  */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/**@file  spxequilisc.h
 * @brief LP equilibrium scaling.
 */
#ifndef _SPXEQUILISC_H_
#define _SPXEQUILISC_H_

#include <assert.h>

#include "soplex/spxdefines.h"
#include "soplex/spxscaler.h"

namespace soplex
{
/**@brief Equilibrium row/column scaling.
   @ingroup Algo

   This SPxScaler implementation performs equilibrium scaling of the
   LPs rows and columns.
*/
template <class R>
class SPxEquiliSC : public SPxScaler<R>
{
public:
   /// compute equilibrium scaling vector rounded to power of two
   static void computeEquiExpVec(const SVSetBase<R>* vecset, const DataArray<int>& coScaleExp,
                                 DataArray<int>& scaleExp);

   /// compute equilibrium scaling vector rounded to power of two
   static void computeEquiExpVec(const SVSetBase<R>* vecset, const std::vector<R>& coScaleVal,
                                 DataArray<int>& scaleExp);

   /// compute equilibrium scaling rounded to power of 2 for existing R scaling factors (preRowscale, preColscale)
   static void computePostequiExpVecs(const SPxLPBase<R>& lp, const std::vector<R>& preRowscale,
                                      const std::vector<R>& preColscale,
                                      DataArray<int>& rowscaleExp, DataArray<int>& colscaleExp);
   //-------------------------------------
   /**@name Construction / destruction */
   ///@{
   /// default constructor (this scaler makes no use of inherited member m_colFirst)
   explicit SPxEquiliSC(bool doBoth = true);
   /// copy constructor
   SPxEquiliSC(const SPxEquiliSC& old);
   /// assignment operator
   SPxEquiliSC& operator=(const SPxEquiliSC&);
   /// destructor
   virtual ~SPxEquiliSC()
   {}
   /// clone function for polymorphism
   inline virtual SPxScaler<R>* clone() const override
   {
      return new SPxEquiliSC<R>(*this);
   }
   ///@}

   //-------------------------------------
   /**@name Scaling */
   ///@{
   /// Scale the loaded SPxLP.
   virtual void scale(SPxLPBase<R>& lp, bool persistent = false) override;
   ///@}
};
} // namespace soplex

#include "spxequilisc.hpp"

#endif // _SPXEQUILISC_H_
