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
#pragma ident "@(#) $Id: clumembers.h,v 1.3 2001/11/21 09:30:12 bzfkocht Exp $"

#ifndef _CLUMEMBERS_H_
#define _CLUMEMBERS_H_

#include "clutypes.h"

namespace soplex
{

typedef struct CLUFactor
{
   int thedim;                 /* dimension of factorized matrix   */
   int stat;                   /* Status indicator. */
   int nzCnt;                  /* number of nonzeros in U      */
   double initMaxabs;     /* maximum abs number in initail Matrix */
   double maxabs;         /* maximum abs number in L and U        */

   double rowMemMult;     /* factor of minimum Memory * #of nonzeros */
   double colMemMult;     /* factor of minimum Memory * #of nonzeros */
   double lMemMult;       /* factor of minimum Memory * #of nonzeros */

   Perm row;            /* row permutation matrices */
   Perm col;            /* column permutation matrices */

   L l;              /* L matrix */
   double *diag;          /* Array of pivot elements          */
   U u;              /* U matrix */

   double *work;          /* Working array: must always be left as 0! */
   double *work2;         /* Working array: must always be left as 0! */
} CLUFactor;

} // namespace soplex
#endif // _CLUMEMBERS_H_

//-----------------------------------------------------------------------------
//Emacs Local Variables:
//Emacs mode:c++
//Emacs c-basic-offset:3
//Emacs tab-width:8
//Emacs indent-tabs-mode:nil
//Emacs End:
//-----------------------------------------------------------------------------
