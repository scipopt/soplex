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
#pragma ident "@(#) $Id: vector_c.h,v 1.1 2001/11/06 16:18:33 bzfkocht Exp $"

#ifndef _VECTOR_C_H_
#define _VECTOR_C_H_

namespace soplex
{
void Vector_MultAddSSVector(double* vec, double x,
                             int n, const int *idx,
                             const double *val);
void Vector_MultAddSVector(double* vec, double x,
                            int n, void *elem);
void Vector_MultAddVector(double x, int n,
                           double* v, const double* w);
void Vector_Set0toSSVector(double*zero, int n,
                            const int* idx, const double* val);
double MultiplyVectorSSVector(const double* dense, int n,
                               const int* idx, const double* val);
double MultiplyVectorVector(const double* v1, int dim,
                             const double* v2);

} // namespace soplex
#endif /* _VECTOR_C_H_ */

//-----------------------------------------------------------------------------
//Emacs Local Variables:
//Emacs c-basic-offset:3
//Emacs tab-width:8
//Emacs indent-tabs-mode:nil
//Emacs End:
//-----------------------------------------------------------------------------
