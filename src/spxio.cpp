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
#pragma ident "@(#) $Id: spxio.cpp,v 1.6 2002/01/11 21:05:30 bzfkocht Exp $"


#include <iostream>
#include <stdio.h>
#include <assert.h>

#include "spxlp.h"

#include "dvector.h"
#include "dataarray.h"
#include "lprow.h"
#include "lpcol.h"
#include "lprowset.h"
#include "lpcolset.h"
#include "nameset.h"

namespace soplex
{
bool SPxLP::read(
   std::istream& is, 
   NameSet* rowNames,
   NameSet* colNames)
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
      ? readMPS(is, rowNames, colNames)
      : readLPF(is, rowNames, colNames);

   //   std::cout << *this;

   return ok;
}

static void dumpRows(std::ostream& s, const SPxLP& lp)
{
   int i;

   s << "\nSubject To\n";
   for (i = 0; i < lp.nRows(); ++i)
   {
      s << "  C" << i << ": ";
      double low;
      double up;
      low = lp.lhs(i);
      up = lp.rhs(i);
      if (low > -SPxLP::infinity && up < SPxLP::infinity && low != up)
      {
         s << low << " <= " << lp.rowVector(i) << " <= " << up << '\n';
      }
      else if (low == up)
         s << lp.rowVector(i) << " = " << up << '\n';
      else if (low <= -SPxLP::infinity)
         s << lp.rowVector(i) << " <= " << up << '\n';
      else
         s << lp.rowVector(i) <<">= " << low << '\n';
   }
}

static void dumpBounds(std::ostream& s, const SPxLP& lp)
{
   int i;

   s << "Bounds\n";
   for (i = 0; i < lp.nCols(); ++i)
   {
      double up;
      double low;
      up = lp.upper(i);
      low = lp.lower(i);
      if (low == up)
         s << "  x" << i << " = " << up << '\n';
      else if (low > -SPxLP::infinity)
      {
         if (up < SPxLP::infinity)
            s << "  " << low << " <= x" << i
            << " <= " << up << '\n';
         else if (low != 0)
            s << "  x" << i <<">= " << low << '\n';
      }
      else if (up < SPxLP::infinity)
         s << "  x" << i << " <= " << up << '\n';
      else
         s << "  x" << i << " FREE\n";
   }
}

std::ostream& operator<<(std::ostream& s, const SPxLP& lp)
{
   int i, j;
   int sns = lp.spxSense();

   s << ((sns == SPxLP::MINIMIZE) ? "Minimize\n" : "Maximize\n");

   s << "  obj: ";
   for (i = j = 0; i < lp.nCols(); ++i)
   {
      double obj = lp.obj(i);
      if (obj != 0)
      {
         if (j)
         {
            if (obj < 0)
               s << " - " << -obj << " x" << i;
            else
               s << " + " << obj << " x" << i;
         }
         else
            s << obj << " x" << i;
         j++;
         if (j % 5 == 0)
            s << "\n\t";
      }
   }
   dumpRows(s, lp);
   dumpBounds(s, lp);

   s << "End\n";
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
