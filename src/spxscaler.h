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
#pragma ident "@(#) $Id: spxscaler.h,v 1.4 2003/01/10 12:46:14 bzfkocht Exp $"

/**@file  spxscaler.h
 * @brief LP scaling base class.
 */
#ifndef _SPXSCALER_H_
#define _SPXSCALER_H_

#include <assert.h>

#include "spxdefines.h"
#include "dataarray.h"
#include "spxlp.h"

namespace soplex
{
/**@brief   LP scaler abstract base class.
   @ingroup Algo

   Instances of classes derived from #SPxScaler may be loaded to #SoPlex in
   order to scale LPs before solving them. #SoPlex# will #load()# itself to
   the #SPxScaler and then call #scale(). Generally any #SPxLP can be
   loaded to a #SPxScaler for #scale()%ing it. The scaling can
   be undone by calling #unscale().
*/
class SPxScaler
{
protected:
   const char*        m_name;      ///< Name of the scaler
   DataArray < Real > m_colscale;  ///< column scaleing factors
   DataArray < Real > m_rowscale;  ///< row scaleing factors
   bool               m_colFirst;  ///< do column scaling first 
   bool               m_doBoth;    ///< do columns and rows

   /// setup scale array for the LP.
   virtual void setup(SPxLP& lp);
   ///
   virtual Real computeColscale(const SVector& col) const = 0;
   ///
   virtual Real computeRowscale(const SVector& row) const = 0;

public:
   friend std::ostream& operator<<(std::ostream& s, const SPxScaler& sc);

   /// constructor
   explicit SPxScaler(const char* name, bool colFirst = false, bool doBoth = true);
   /// copy constructor
   SPxScaler(const SPxScaler& old);
   /// destructor.
   virtual ~SPxScaler();
   /// assignment operator
   SPxScaler& operator=(const SPxScaler& rhs);

   /// get name of scaler.
   virtual const char* getName() const;
   /// set scaling order.
   virtual void setOrder(bool colFirst); 
   /// set wether column and row scaling should be performed.
   virtual void setBoth(bool both); 

   /// Scale #SPxLP. 
   virtual void scale(SPxLP& lp);
   /// Unscale dense primal solution vector given in \p x. 
   virtual void unscalePrimal(Vector& x) const;
   /// Unscale dense dual solution vector given in \p pi. 
   virtual void unscaleDual(Vector& pi) const;

   /// absolute smallest column scaling factor
   virtual Real minAbsColscale() const;
   /// absolute biggest column scaling factor
   virtual Real maxAbsColscale() const;
   /// absolute smallest row scaling factor
   virtual Real minAbsRowscale() const;
   /// absolute biggest row scaling factor
   virtual Real maxAbsRowscale() const;
   /// maximum ratio between absolute biggest and smallest element in any column.
   virtual Real maxColRatio(const SPxLP& lp) const;
   /// maximum ratio between absolute biggest and smallest element in any row.
   virtual Real maxRowRatio(const SPxLP& lp) const;

   /// consistency check
   virtual bool isConsistent() const;
};
} // namespace soplex
#endif // _SPXSCALER_H_

//-----------------------------------------------------------------------------
//Emacs Local Variables:
//Emacs mode:c++
//Emacs c-basic-offset:3
//Emacs tab-width:8
//Emacs indent-tabs-mode:nil
//Emacs End:
//-----------------------------------------------------------------------------

