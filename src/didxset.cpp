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
#pragma ident "@(#) $Id: didxset.cpp,v 1.9 2002/01/31 08:19:26 bzfkocht Exp $"

#include "didxset.h"
#include "spxalloc.h"

namespace soplex
{

void DIdxSet::setMax(int newmax)
{ 
   assert(idx   != 0);
   assert(max() >  0);

   len = (newmax < size()) ? size() : newmax;
   len = (len < 1) ? 1 : len;

   assert(len > 0);

   spx_realloc(idx, len);
}

DIdxSet::DIdxSet(const IdxSet& old)
   : IdxSet(0, 0)
{
   len = old.size();
   len = (len < 1) ? 1 : len;
   spx_alloc(idx, len);

   IdxSet::operator= ( old );
}

DIdxSet::DIdxSet(const DIdxSet& old)
   : IdxSet(0, 0)
{
   len = old.size();
   len = (len < 1) ? 1 : len;
   spx_alloc(idx, len);

   IdxSet::operator= ( old );
}

DIdxSet::DIdxSet(int n)
   : IdxSet()
{
   len = (n < 1) ? 1 : n;
   spx_alloc(idx, len);
}

DIdxSet::~DIdxSet()
{
   spx_free(idx);
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
