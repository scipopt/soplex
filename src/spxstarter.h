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
#pragma ident "@(#) $Id: spxstarter.h,v 1.9 2005/09/16 12:42:36 bzfhille Exp $"


/**@file  spxstarter.h
 * @brief #SoPlex start basis generation base class.
 */
#ifndef _SPXDSTARTER_H_
#define _SPXDSTARTER_H_

#include <assert.h>

#include "spxdefines.h"
#include "spxsolver.h"

namespace soplex
{

/**@brief   #SoPlex start basis generation base class.
   @ingroup Algo
   
   #SPxStarter is the virtual base class for classes generating a starter basis
   for the Simplex solver #SoPlex. When a #SPxStarter object has been loaded
   to a #SoPlex solver, the latter will call method #generate() in order to
   have a start basis generated. Implementations of method #generate() must
   terminate by #SPxSolver::load()%ing the generated basis to #SoPlex. Loaded
   basises must be nonsingular.
*/
class SPxStarter
{
protected:
   const char* m_name;

public:
   /// constructor
   explicit SPxStarter(const char* name)
      : m_name(name)
   {}
   /// destructor.
   virtual ~SPxStarter()
   {
      m_name = 0;
   }
   /// get name of starter.
   virtual const char* getName() const
   {
      return m_name;
   }
   /// generates start basis for loaded basis.
   virtual void generate(SPxSolver& base) = 0;
#ifndef NO_CONSISTENCY_CHECKS
   /// checks consistency.
   virtual bool isConsistent() const;
#endif
};


} // namespace soplex
#endif // _SPXDSTARTER_H_

//-----------------------------------------------------------------------------
//Emacs Local Variables:
//Emacs mode:c++
//Emacs c-basic-offset:3
//Emacs tab-width:8
//Emacs indent-tabs-mode:nil
//Emacs End:
//-----------------------------------------------------------------------------
