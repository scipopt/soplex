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
#pragma ident "@(#) $Id: dsvector.cpp,v 1.10 2002/01/19 13:06:29 bzfkocht Exp $"

#include <assert.h>
#include <iostream>

#include "dsvector.h"
#include "spxalloc.h"
#include "message.h"

namespace soplex
{
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

DSVector::DSVector(const SSVector& old)
{
   allocMem(old.size() + 1);
   SVector::operator= ( old );
}

DSVector::DSVector(const DSVector& old)
   : SVector()
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

bool DSVector::isConsistent() const
{
   if ((theelem != 0) && (mem() != theelem))
      return MSGinconsistent("DSVector");

   return true;
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


