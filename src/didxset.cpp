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
#pragma ident "@(#) $Id: didxset.cpp,v 1.2 2001/11/06 23:31:01 bzfkocht Exp $"


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
   idx = (int*)realloc(idx, len * sizeof(int));
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
   idx = (int*)malloc(len * sizeof(int));
   if (idx == 0)
   {
      std::cerr << "ERROR: DIdxSet could not allocate memory\n";
      exit(-1);
   }
   *(IdxSet*)this = (IdxSet&)old;
}

DIdxSet::DIdxSet(const DIdxSet& old)
   : IdxSet(0, 0)
{
   len = old.size();
   idx = (int*)malloc(len * sizeof(int));
   if (idx == 0)
   {
      std::cerr << "ERROR: DIdxSet could not allocate memory\n";
      exit(-1);
   }
   *(IdxSet*)this = (IdxSet&)old;
}

DIdxSet::DIdxSet(int n)
   : IdxSet()
{
   len = (n < 1) ? 1 : n;
   idx = (int*)malloc(len * sizeof(int));
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
