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
#pragma ident "@(#) $Id: didxset.cpp,v 1.3 2001/11/07 17:31:15 bzfbleya Exp $"


/*      \Section{Complex Members}
 */
#include <iostream>
#include "didxset.h"

namespace soplex
{

void DIdxSet::setMax(int newmax)
{
   len = (newmax < size()) ? size() : newmax;
   len = (len < 1) ? 1 : len;
   idx = reinterpret_cast<int*>(realloc(idx, len * sizeof(int)));
   if (idx == 0)
   {
      std::cerr << "ERROR: DIdxSet could not reallocate memory\n";
      exit(-1);
   }
}

DIdxSet::DIdxSet(const IdxSet& old)
   : IdxSet(0, 0)
{
   len = old.size();
   idx = reinterpret_cast<int*>(malloc(len * sizeof(int)));
   if (idx == 0)
   {
      std::cerr << "ERROR: DIdxSet could not allocate memory\n";
      exit(-1);
   }
   IdxSet::operator= ( old );
}

DIdxSet::DIdxSet(const DIdxSet& old)
   : IdxSet(0, 0)
{
   len = old.size();
   idx = reinterpret_cast<int*>(malloc(len * sizeof(int)));
   if (idx == 0)
   {
      std::cerr << "ERROR: DIdxSet could not allocate memory\n";
      exit(-1);
   }
   IdxSet::operator= ( old );
}

DIdxSet::DIdxSet(int n)
   : IdxSet()
{
   len = (n < 1) ? 1 : n;
   idx = reinterpret_cast<int*>(malloc(len * sizeof(int)));
   if (idx == 0)
   {
      std::cerr << "ERROR: DIdxSet could not allocate memory\n";
      exit(-1);
   }
}

DIdxSet::~DIdxSet()
{
   free(idx);
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
