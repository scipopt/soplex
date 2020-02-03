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

/**@file  spxleastsqsc.h
 * @brief LP least squares scaling.
 */
#ifndef _SPXLEASTSQSC_H_
#define _SPXLEASTSQSC_H_

#include <assert.h>
#include "soplex/spxsolver.h"
#include "soplex/spxscaler.h"
#include "soplex/spxlp.h"


namespace soplex
{

/**@brief Least squares scaling.
   @ingroup Algo

   This SPxScaler implementation performs least squares scaling as suggested by Curtis and Reid in:
   On the Automatic Scaling of Matrices for Gaussian Elimination (1972).
*/
template <class R>
class SPxLeastSqSC : public SPxScaler<R>
{
public:
   //-------------------------------------
   /**@name Construction / destruction */
   ///@{
   /// default constructor (this scaler makes no use of inherited member m_colFirst)
   explicit SPxLeastSqSC();
   /// copy constructor
   SPxLeastSqSC(const SPxLeastSqSC& old);
   /// assignment operator
   SPxLeastSqSC& operator=(const SPxLeastSqSC&);
   /// destructor
   virtual ~SPxLeastSqSC()
   {}
   /// clone function for polymorphism
   inline virtual SPxScaler<R>* clone() const override
   {
      return new SPxLeastSqSC(*this);
   }
   ///@}

   //-------------------------------------
   /**@name Access / modification */
   ///@{
   /// set Real param (conjugate gradient accuracy)
   virtual void setRealParam(R param, const char* name) override;
   /// set int param (maximal conjugate gradient rounds)
   virtual void setIntParam(int param, const char* name) override;
   ///@}

   //-------------------------------------
   /**@name Scaling */
   ///@{
   /// Scale the loaded SPxLP.
   virtual void scale(SPxLPBase<R>& lp, bool persistent = true) override;


private:
   R acrcydivisor = R(1000.0);
   int maxrounds = 20;

};
} // namespace soplex

// For the general templated files
#include "spxleastsqsc.hpp"

#endif // _SPXLEASTSQSC_H_

