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
#pragma ident "@(#) $Id: spxfileio.cpp,v 1.2 2002/12/14 07:48:30 bzfkocht Exp $"

//#define DEBUGGING 1
#define USE_GZSTREAM  1

#include <assert.h>
#include <iostream>
#include <fstream>

#if USE_GZSTREAM
#include "gzstream.h"
#endif // USE_GZSTREAM

#include "spxdefines.h"
#include "soplex.h"

namespace soplex
{

bool SoPlex::readBasisFile(
   const char*    filename, 
   const NameSet& rowNames,
   const NameSet& colNames)
{
   METHOD( "SoPlex::readBasisFile()" );

#if USE_GZSTREAM
   gzstream::igzstream file(filename);
#else
   std::ifstream file(filename);
#endif // USE_GZSTREAM

   if (!file)
      return false;
 
   return readBasis(file, rowNames, colNames);
}

bool SoPlex::writeBasisFile(
   const char*    filename, 
   const NameSet& rowNames,
   const NameSet& colNames)
{
   METHOD( "SoPlex::writeBasisFile()" );
   std::ofstream file(filename);

   if (!file)
      return false;
 
   writeBasis(file, rowNames, colNames);

   return true;
}


bool SoPlex::readFile( 
   const char* filename, 
   NameSet*    rowNames,
   NameSet*    colNames, 
   DIdxSet*    intVars)
{
   METHOD( "SoPlex::readFile()" );

#if USE_GZSTREAM
   gzstream::igzstream file(filename);
#else
   std::ifstream file(filename);
#endif // USE_GZSTREAM

   if (!file)
      return false;

   return read(file, rowNames, colNames, intVars);
}

void SoPlex::dumpFile(const char* filename) const
{
   METHOD( "SoPlex::dumpFile()" );
   std::ofstream file(filename);

   if (file.good())
      file << *this;
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
