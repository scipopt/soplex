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
#pragma ident "@(#) $Id: spxgeometsc.cpp,v 1.1 2003/01/05 19:03:16 bzfkocht Exp $"

/**@file  spxgeometsc.cpp
 * @brief Geometric mean row/column scaling.
 */
#include <assert.h>

#include "spxgeometsc.h"

namespace soplex
{
static const char* makename(bool colFirst, bool doBoth)
{
   const char* name;

   if (doBoth)
      name = colFirst ? "CR-Geometric" : "RC-Geometric";
   else
      name = colFirst ? "C-Geometric" : "R-Geometric";

   return name;
}

SPxGeometSC::SPxGeometSC(bool colFirst, bool doBoth)
   : SPxScaler(makename(colFirst, doBoth), colFirst, doBoth)
{}


Real SPxGeometSC::computeColscale(const SVector& col) const
{
   Real mini = col.minAbs();
   Real maxi = col.maxAbs();

   return sqrt(mini * maxi);
}

Real SPxGeometSC::computeRowscale(const SVector& row) const
{
   Real mini = row.minAbs();
   Real maxi = row.maxAbs();

   return sqrt(mini * maxi);
}

void SPxGeometSC::scale(SPxLP& lp) 
{
   VERBOSE2({ std::cout << "Geometric scaling LP" << std::endl; });   

   VERBOSE3({ std::cout << "\tNon zero min= " << lp.minAbsNzo()
                        << " max= " << lp.maxAbsNzo()
                        << std::endl; });
   
   SPxScaler::scale(lp);

   VERBOSE3({ std::cout << "\tRow scaling min= " << minAbsRowscale()
                        << " max= " << maxAbsRowscale()
                        << std::endl
                        << "\tCol scaling min= " << minAbsColscale()
                        << " max= " << maxAbsColscale()
                        << std::endl
                        << "\tNon zero min= " << lp.minAbsNzo()
                        << " max= " << lp.maxAbsNzo()
                        << std::endl; });
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



