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
#pragma ident "@(#) $Id: spxalloc.h,v 1.2 2001/11/13 21:55:19 bzfkocht Exp $"

#ifndef _SPXALLOC_H_
#define _SPXALLOC_H_

#include <iostream>
#include <stdlib.h>
#include <assert.h>

namespace soplex
{
///
template <class T>
inline void spx_alloc(T& p, int n)
{
   assert(n >= 0);

   if (n == 0)
      n = 1;
   
   p = reinterpret_cast<T>(malloc(sizeof(*p) * n));

   if (0 == p)
   {
      std::cerr << "malloc: Out of memory" << std::endl;
      abort();
   }
}

/// 
template <class T>
inline void spx_realloc(T& p, int n)
{
   assert(p != 0);
   assert(n >= 0);

   p = reinterpret_cast<T>(realloc(p, sizeof(*p) * n));

   if (0 == p)
   {
      std::cerr << "realloc: Out of memory" << std::endl;
      abort();
   }
}

///
template <class T>
inline void spx_free(T& p)
{
   assert(p != 0);

   free(p);
   
   p = 0;
}

} // namespace soplex
#endif // _SPXALLOC_H_

//-----------------------------------------------------------------------------
//Emacs Local Variables:
//Emacs mode:c++
//Emacs c-basic-offset:3
//Emacs tab-width:8
//Emacs indent-tabs-mode:nil
//Emacs End:
//-----------------------------------------------------------------------------



