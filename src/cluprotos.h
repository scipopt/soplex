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
#pragma ident "@(#) $Id: cluprotos.h,v 1.9 2001/12/04 18:25:56 bzfkocht Exp $"


#ifndef _CLUPROTOS_H_
#define _CLUPROTOS_H_


#define WITH_L_ROWS             1

#if 0
#define CLU_OK                  0
#define CLU_INSTABLE            1
#define CLU_SINGULAR            2
#define CLU_OUT_OF_MEMORY       8

#include "svector.h"
#include "clutypes.h"
#endif
namespace soplex
{

#if 0

int updateCLUFactor (CLUFactor*, int, double*, const int*, int);
int updateCLUFactorNoClear(CLUFactor*, int, const double*, const int*, int);


int forestUpdateCLUFactor (CLUFactor* fac, int col, double* work, int n, int *nonz);

int factor(CLUFactor*,
            SVector**,       /* Array of column vector pointers   */
            double,          /* pivoting threshold                */
            double           /* epsilon for zero detection        */
         );
void solveLright(CLUFactor* fac, double* vec);
void remaxRow(CLUFactor*, int, int);
void remaxCol(CLUFactor*, int, int);
void packRows(CLUFactor* fac);
int CLUFactorIsConsistent (const CLUFactor*);
int makeLvec(CLUFactor*, int, int);
void dumpCLUFactor(const CLUFactor* fac);
#endif

} // namespace soplex
#endif // _CLUPROTOS_H_

//-----------------------------------------------------------------------------
//Emacs Local Variables:
//Emacs mode:c++
//Emacs c-basic-offset:3
//Emacs tab-width:8
//Emacs indent-tabs-mode:nil
//Emacs End:
//-----------------------------------------------------------------------------
