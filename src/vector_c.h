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
#pragma ident "@(#) $Id: vector_c.h,v 1.4 2001/11/23 14:34:55 bzfbleya Exp $"

/**@file  vector_v.h
 * @brief Algebraic functions for vectors.
 */
#ifndef _VECTOR_C_H_
#define _VECTOR_C_H_

namespace soplex
{
/**@name C-implementation of vector algebra and utils 
   @ingroup Algebra 
   
   Several algebraic and utility methods of vector are implemented in
   C-code for efficiency.
*/
//@{

/**@brief   C-implementation of Vector::multAdd
   @ingroup Algebra 
*/
void Vector_MultAddSSVector(double* vec, double x,
                            int n, const int *idx,
                            const double *val);
/**@brief   C-implementation of Vector::multAdd
   @ingroup Algebra 
*/
void Vector_MultAddSVector(double* vec, double x,
                           int n, void *elem);
/**@brief   C-implementation of Vector::multAdd
   @ingroup Algebra 
*/
void Vector_MultAddVector(double x, int n,
                          double* v, const double* w);
/**@brief   C-implementation to assign 0-vector to SSVector
   @ingroup Algebra 
*/
void Vector_Set0toSSVector(double*zero, int n,
                           const int* idx, const double* val);
/**@brief   C-implementation of Vector::operator*
   @ingroup Algebra 
*/
double MultiplyVectorSSVector(const double* dense, int n,
                              const int* idx, const double* val);
/**@brief   C-implementation of Vector::operator*
   @ingroup Algebra 
*/
double MultiplyVectorVector(const double* v1, int dim,
                            const double* v2);
/**@brief   C-implementation of UpdateVector::update
   @ingroup Algebra 
*/
void UpdateUpdateVector(double*, double, int, const int*, const double*);

//@}

} // namespace soplex
#endif /* _VECTOR_C_H_ */

//-----------------------------------------------------------------------------
//Emacs Local Variables:
//Emacs mode:c++
//Emacs c-basic-offset:3
//Emacs tab-width:8
//Emacs indent-tabs-mode:nil
//Emacs End:
//-----------------------------------------------------------------------------
