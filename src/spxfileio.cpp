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
#pragma ident "@(#) $Id: spxfileio.cpp,v 1.6 2006/08/25 18:17:18 bzforlow Exp $"

//#define DEBUGGING 1

#include <assert.h>

#include "spxdefines.h"
#include "spxsolver.h"
#include "spxfileio.h"

namespace soplex
{
bool SPxSolver::readBasisFile(
   const char*    filename, 
   const NameSet& rowNames,
   const NameSet& colNames)
{
   METHOD( "SPxSolver::readBasisFile()" );

   spxifstream file(filename);

   if (!file)
      return false;
 
   return readBasis(file, rowNames, colNames);
}

bool SPxSolver::writeBasisFile(
   const char*    filename, 
   const NameSet& rowNames,
   const NameSet& colNames)
   const
{
   METHOD( "SPxSolver::writeBasisFile()" );
   std::ofstream file(filename);

   if (!file)
      return false;
 
   writeBasis(file, rowNames, colNames);

   return true;
}

bool SPxLP::readFile( 
   const char* filename, 
   NameSet*    rowNames,
   NameSet*    colNames, 
   DIdxSet*    intVars)
{
   METHOD( "SPxLP::readFile()" );

   spxifstream file(filename);

   if (!file)
      return false;

   return read(file, rowNames, colNames, intVars);
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
