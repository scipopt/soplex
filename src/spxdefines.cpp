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
#pragma ident "@(#) $Id: spxdefines.cpp,v 1.4 2002/08/27 07:20:37 bzfkocht Exp $"

/**@file  spxdefines.cpp
 * @brief Debugging, floating point type and parameter definitions.
 */
#include "spxdefines.h"

namespace soplex
{

#if defined(TRACE_METHOD)

int TraceMethod::s_indent = 0;

#endif //TRACE_METHOD

const Real infinity                 = DEFAULT_INFINITY;

Real Param::s_epsilon               = DEFAULT_EPS_ZERO;
Real Param::s_epsilon_factorization = DEFAULT_EPS_FACTOR;
Real Param::s_epsilon_update        = DEFAULT_EPS_UPDATE;
int  Param::s_verbose               = 1;
   
void Param::setEpsilon(Real eps)
{
   s_epsilon = eps;
}

void Param::setVerbose(int p_verbose)
{
   s_verbose = p_verbose;
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

