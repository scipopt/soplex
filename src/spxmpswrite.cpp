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
#pragma ident "@(#) $Id: spxmpswrite.cpp,v 1.1 2002/03/06 10:28:52 bzfkocht Exp $"

/**@file  spxmpswrite.cpp
 * @brief Write LP as MPS format file.
 */
//#define DEBUGGING 1

#include <assert.h>
#include <stdio.h>
#include <iostream>
#include <iomanip>

#include "spxdefines.h"
#include "spxlp.h"

namespace soplex
{
static void writeRecord(
   std::ostream&  os, 
   const char*    indicator,
   const char*    name,
   const char*    name1  = 0,
   const Real     value1 = 0.0,
   const char*    name2  = 0,
   const Real     value2 = 0.0) 
{
#if 0
   //scientific std::ios::showpoint | 
   os.setf(std::ios::scientific | std::ios::left);

   os << " " 
      << std::setw(2) <<  ((indicator == 0) ? "" : indicator)
      << std::setw(8) <<  ((name == 0)      ? "" : name);

   if (name1 != 0)
   {
      os << "  "
         << std::setw(8) << name1
         << "  "
         << std::setw(12) << std::setprecision(9) << value1;

      if (name2 != 0)
      {
      os << "   "
         << std::setw(8) << name2
         << "  "
         << std::setw(12) << std::setprecision(9) << value2;
      }
   }
   os << std::endl;
#endif   
   char buf[81];

   sprintf(buf, " %-2.2s %-8.8s",
      (indicator == 0) ? "" : indicator,
      (name == 0)      ? "" : name);

   os << buf;

   if (name1 != 0)
   {
      sprintf(buf, "  %-8.8s  %12.9g", name1, value1);
      os << buf;

      if (name2 != 0)
      {
         sprintf(buf, "   %-8.8s  %12.9g", name2, value2);
         os << buf;
      }
   }
   os << std::endl;
}

static Real getRHS(Real left, Real right)
{
   Real rhsval;

   ///@todo Ranges are not handle by this.
   if (left > -infinity)
      rhsval = left;
   else if (right <  infinity)
      rhsval = right;
   else
      ABORT();

   return rhsval;
}

/**@todo Ranges are not yet supported.
 */
void SPxLP::writeMPS(
   std::ostream&  p_output, 
   const NameSet* p_rnames,          ///< row names.
   const NameSet* p_cnames,          ///< column names.
   const DIdxSet* /*p_intvars*/)     ///< integer variables.
   const
{
   METHOD("writeMPS");

   assert(p_rnames != 0);
   assert(p_cnames != 0);

   const char* indicator;

   int i;
   int k;

   // --- NAME Section ---
   p_output << "NAME          MPSDATA" << std::endl;

   // --- ROWS Setction ---
   p_output << "ROWS" << std::endl;

   for(i = 0; i < nRows(); i++)
   {
      if (lhs(i) == rhs(i))
         indicator = "E";
      else if ((lhs(i) > -infinity) && (rhs(i) < infinity))
         indicator = "R";
      else if (lhs(i) > -infinity)
         indicator = "G";
      else if (rhs(i) <  infinity)
         indicator = "L";
      else
         ABORT();

      writeRecord(p_output, indicator, (*p_rnames)[rId(i)]); 
   }
   writeRecord(p_output, "N", "MINIMIZE"); 
   
   // --- COLUMNS Section ---
   p_output << "COLUMNS" << std::endl;

   for(i = 0; i < nCols(); i++)
   {
      const SVector& col = colVector(i);
      int colsize2       = (col.size() / 2) * 2;

      assert(colsize2 % 2 == 0);

      for(k = 0; k < colsize2; k += 2)
         writeRecord(p_output, 0, (*p_cnames)[cId(i)], 
            (*p_rnames)[rId(col.index(k))    ], col.value(k), 
            (*p_rnames)[rId(col.index(k + 1))], col.value(k + 1));

      if (colsize2 != col.size())
         writeRecord(p_output, 0, (*p_cnames)[cId(i)], 
            (*p_rnames)[rId(col.index(k))], col.value(k));

      if (isNotZero(maxObj(i)))
         writeRecord(p_output, 0, (*p_cnames)[cId(i)], "MINIMIZE", -maxObj(i));
   }
   // --- RHS Section ---
   p_output << "RHS" << std::endl;

   i = 0;
   while(i < nRows())
   {
      Real rhsval1;
      Real rhsval2;

      for(; i < nRows(); i++)
         if ((rhsval1 = getRHS(lhs(i), rhs(i))) != 0.0)
            break;

      if (i < nRows())
      {
         for(k = i + 1; k < nRows(); k++)
            if ((rhsval2 = getRHS(lhs(k), rhs(k))) != 0.0)
               break;

         if (k < nRows())
            writeRecord(p_output, 0, "RHS", 
               (*p_rnames)[rId(i)], rhsval1, 
               (*p_rnames)[rId(k)], rhsval2);
         else
            writeRecord(p_output, 0, "RHS", 
               (*p_rnames)[rId(i)], rhsval1); 

         i = k + 1;
      }
   }

   // --- RANGES Section ---
   // --- BOUNDS Section ---
   p_output << "BOUNDS" << std::endl;

   for(i = 0; i < nCols(); i++)
   {
      if (lower(i) == upper(i))
      {
         writeRecord(p_output, "FX", "BOUND", (*p_cnames)[cId(i)], lower(i));
         continue;
      }
      if ((lower(i) <= -infinity) && (upper(i) >= infinity))
      {
         writeRecord(p_output, "FR", "BOUND", (*p_cnames)[cId(i)]);
         continue;
      }
      if (lower(i) != 0.0)
      {
         if (lower(i) > -infinity)
            writeRecord(p_output, "LO", "BOUND", 
               (*p_cnames)[cId(i)], lower(i));
         else
            writeRecord(p_output, "MI", "BOUND", (*p_cnames)[cId(i)]);
      }
      if (upper(i) < infinity)
      {
         writeRecord(p_output, "UP", "BOUND", (*p_cnames)[cId(i)], upper(i));
         continue;
      }
   }   
   // --- ENDATA Section ---
   p_output << "ENDATA" << std::endl;   
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
