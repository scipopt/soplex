/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the class library                   */
/*       SoPlex --- the Sequential object-oriented simPlex.                  */
/*                                                                           */
/*    Copyright (C) 1997-1999 Roland Wunderling                              */
/*                  1997-2007 Konrad-Zuse-Zentrum                            */
/*                            fuer Informationstechnik Berlin                */
/*                                                                           */
/*  SoPlex is distributed under the terms of the ZIB Academic Licence.       */
/*                                                                           */
/*  You should have received a copy of the ZIB Academic License              */
/*  along with SoPlex; see the file COPYING. If not email to soplex@zib.de.  */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
#pragma ident "@(#) $Id: spxio.cpp,v 1.26 2007/08/27 15:35:11 bzfberth Exp $"


//#define DEBUGGING 1

#include <iostream>
#include <assert.h>

#include "spxdefines.h"
#include "spxlp.h"

#include "dvector.h"
#include "dataarray.h"
#include "lprow.h"
#include "lpcol.h"
#include "lprowset.h"
#include "lpcolset.h"
#include "nameset.h"
#include "spxout.h"

namespace soplex
{
/**@param is       input stream. 
 * @param rowNames contains after the call the names of the constraints 
 *                 (rows) in the same order as the rows in the LP.
 *                 Constraints without a name (only possible with LPF 
 *                 files) are automatically assigned a name.
 *                 May be NULL if the names are not needed.
 * @param colNames contains after the call the names of the variables
 *                 (columns) in the same order as the columns in the LP.
 *                 May be NULL if the names are not needed.
 * @param intVars  contains after the call the indices of those variables
 *                 that where marked as being integer in the file.
 *                 May be NULL if the information is not needed.
 * @todo  Make sure the Ids in the NameSet%s are the same as in the LP.
 */
bool SPxLP::read(
   std::istream& is, 
   NameSet* rowNames,
   NameSet* colNames,
   DIdxSet* intVars)
{
   bool ok;
   char c;

   is.get(c);
   is.putback(c);

   /* MPS starts either with a comment mark '*' or with the keyword
    * 'NAME' at the first column.
    * LPF starts either with blanks, a comment mark '\' or with
    * the keyword "MAX" or "MIN" in upper or lower case.
    * There is no possible valid LPF file starting with a '*' or 'N'.
    */
   ok = ((c == '*') || (c == 'N'))
      ? readMPS(is, rowNames, colNames, intVars)
      : readLPF(is, rowNames, colNames, intVars);

   MSG_DEBUG( spxout << "DSPXIO01\n" << *this; );

   return ok;
}

// LP output operator with default row and column names
std::ostream& operator<<(std::ostream& s, const SPxLP& lp)
{
   lp.writeLPF(s, NULL, NULL, NULL);
   return s;
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
