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
#pragma ident "@(#) $Id: dsvector.cpp,v 1.6 2001/11/13 21:55:17 bzfkocht Exp $"

/*      \Section{Complex Methods}
 */

/*  Import system include files
 */
#include <assert.h>
#include <iostream>


/*  and class header files
 */
#include "dsvector.h"
#include "spxalloc.h"

namespace soplex
{



//@ ----------------------------------------------------------------------------
/*      \SubSection{Consistency}
 */
void DSVector::allocMem(int len)
{
   spx_alloc(theelem, len);
   setMem(len, theelem);
}

void DSVector::setMax(int newmax)
{
   int siz = size();
   int len = ((newmax < siz) ? siz : newmax) + 1;

   spx_realloc(theelem, len);
   setMem (len, theelem);
   set_size( siz );
}

DSVector& DSVector::operator=(const Vector& vec)
{
   clear();
   setMax(vec.dim());
   SVector::operator=(vec);
   return *this;
}

DSVector& DSVector::assign(const Vector& vec, double eps)
{
   clear();
   setMax(vec.dim());
   SVector::assign(vec, eps);
   return *this;
}

DSVector::DSVector(const SVector& old)
{
   allocMem(old.size() + 1);
   SVector::operator= ( old );
}

DSVector::DSVector(const DSVector& old):
   SVector()
{
   allocMem(old.size() + 1);
   SVector::operator= ( old );
}

DSVector::DSVector(int n)
{
   allocMem((n < 1) ? 2 : n + 1);
}

DSVector::DSVector(const Vector& vec)
{
   allocMem((vec.dim() < 1) ? 2 : vec.dim() + 1);
   *this = vec;
}

DSVector::~DSVector()
{
   spx_free(theelem);
}

int DSVector::isConsistent() const
{
   if (theelem && m_elem - 1 != theelem)
   {
      std::cout << "ERROR: Inconsistency detected in class DSVector\n";
      return 0;
   }
   return 1;
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


