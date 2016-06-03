/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the class library                   */
/*       SoPlex --- the Sequential object-oriented simPlex.                  */
/*                                                                           */
/*    Copyright (C) 1996-2016 Konrad-Zuse-Zentrum                            */
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
#include "spxsolver.h"
#include "spxscaler.h"
#include "spxlp.h"

#define MAX_ROUNDS 50
#define ACCURACY_DIVISOR 1000.0

namespace soplex
{
/**@brief Least squares scaling.
   @ingroup Algo

   This SPxScaler implementation performs least squares scaling as suggested by Curtis and Reid in:
   On the Automatic Scaling of Matrices for Gaussian Elimination (1972).
*/
class SPxLeastSqSC : public SPxScaler
{
public:

   //-------------------------------------
   /**@name Construction / destruction */
   //@{
   /// default constructor (this scaler makes no use of inherited member m_colFirst)
   explicit SPxLeastSqSC();
   /// copy constructor
   SPxLeastSqSC(const SPxLeastSqSC& old);
   /// assignment operator
   SPxLeastSqSC& operator=(const SPxLeastSqSC& );
   /// destructor
   virtual ~SPxLeastSqSC()
   {}
   /// clone function for polymorphism
   inline virtual SPxScaler* clone() const
   {
      return new SPxLeastSqSC(*this);
   }
   //@}

   //-------------------------------------
   /**@name Access / modification */
   //@{
   /// set real param (conjugate gradient accuracy)
   virtual void setRealParam(Real param, const char* name);
   /// set int param (maximal conjugate gradient rounds)
   virtual void setIntParam(int param, const char* name);
   //@}

   //-------------------------------------
   /**@name Scaling */
   //@{
   /// Scale the loaded SPxLP.
   virtual void scale(SPxLP& lp);


protected:

   //-------------------------------------
   /**@name Protected helpers */
   //@{
   /// Does nothing but returning \p maxi.
   virtual Real computeScale(Real /*mini*/, Real maxi) const;
   //@}

   Real acrcydivisor = ACCURACY_DIVISOR;
   int maxrounds = MAX_ROUNDS;

};
} // namespace soplex
#endif // _SPXLEASTSQSC_H_
