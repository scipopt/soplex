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
#pragma ident "@(#) $Id: spxmps.cpp,v 1.8 2001/11/21 09:30:15 bzfkocht Exp $"


#include <iostream>
#include <stdio.h>
#include <ctype.h>
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





//@ -----------------------------------------------------------------------------

#define LINE    81
static int line;
static char objName[LINE];

enum MPS_Section
{
   ROWS,
   COLUMNS,
   RHS,
   RANGES,
   BOUNDS,
   ENDATA
};

#define TRACE        0

#define DATA_FIELDS   6

#define POS_INDIC     1
#define POS_NAME      4
#define POS_NAME1    14
#define POS_VALUE1   24
#define POS_NAME2    39
#define POS_VALUE2   49
#define POS_END      61

#define SRT_INDIC     0
#define SRT_NAME      1
#define SRT_NAME1     2
#define SRT_VALUE1    3
#define SRT_NAME2     4
#define SRT_VALUE2    5

int SPxLP::readLine(
   std::istream& is,
   char*& f1,
   char*& f2,
   char*& f3,
   char*& f4,
   char*& f5,
   char*& f6)
{
   static const int beg [] =
      {
         1, 4, 14, 24, 39, 49, -1
      };
   static const int end [] =
      {
         3, 12, 22, 36, 47, 61, -1
      };

   static char buf [1024];

   char* srt[DATA_FIELDS];
   int length;
   int i;
   int k;

   /* CONSTCOND */
   assert((sizeof(beg) / sizeof(int)) == DATA_FIELDS + 1);
   /* CONSTCOND */
   assert((sizeof(end) / sizeof(int)) == DATA_FIELDS + 1);

   while (is.getline(buf, sizeof(buf)))
   {
      line++;

      /* Maximale inhaltsvolle Zeilenlaenge
       */
      buf[POS_END] = '\0';

      /* Kommentare, Leerzeilen
       */
      if ((buf[0] == '*') || (buf[0] == '\n'))
         continue;

#if TRACE
      std::cerr << '[' << buf << ']' << std::endl;
#endif
      /* Neue Sektion
       */
      if (buf[0] != ' ')
      {
         f1 = &buf[0];
         f2 = 0;
         f3 = 0;
         f4 = 0;
         f5 = 0;
         f6 = 0;

         return 1;
      }
      length = static_cast<int>(strlen(buf));   // -1

      /* Split Sections
       */
      for (i = 0; beg[i] > 0; i++)
      {
         if (beg[i] >= length)
         {
            buf[beg[i]] = '\0';
            srt[i] = &buf[beg[i]];
         }
         else
         {
            /* Skip leading blanks, tabs, etc.
             */
            for (k = beg[i]; (k < end[i]) && isspace(buf[k]); k++)
              ;
            /* Remember start of data in the field.
             */
            srt[i] = &buf[k];

            /* Skip over field-data.
             */
            for (; (k < end[i]) && isgraph(buf[k]); k++)
              ;
            /* Null out the rest. (And put an ending null at the end.)
            */
            for (; k <= end[i]; k++)
               buf[k] = '\0';
         }
      }
#if TRACE > 1
      for (i = beg[1]; i <= end[1]; i++)
         std::cout << i << ':' << buf[i] << std::endl;
#endif
      if (*srt[SRT_NAME1] == '$')
      {
         *srt[SRT_NAME1] = '\0';
         *srt[SRT_VALUE1] = '\0';
         *srt[SRT_NAME2] = '\0';
         *srt[SRT_VALUE2] = '\0';
      }

      if (*srt[SRT_NAME2] == '$')
      {
         *srt[SRT_NAME2] = '\0';
         *srt[SRT_VALUE2] = '\0';
      }

#if TRACE
      (void)fprintf(stderr, "[%s]\n", srt[SRT_NAME]);
      (void)fprintf(stderr, "[%s]\n", srt[SRT_NAME1]);
      (void)fprintf(stderr, "[%s]\n", srt[SRT_NAME2]);
      (void)fprintf(stderr, "[%s]\n", srt[SRT_VALUE1]);
      (void)fprintf(stderr, "[%s]\n", srt[SRT_VALUE2]);
#endif
      f1 = srt[0];
      f2 = srt[1];
      f3 = srt[2];
      f4 = srt[3];
      f5 = srt[4];
      f6 = srt[5];

      return 0;
   }
   return -1;
}

inline int SPxLP__readLine
(
   std::istream& is,
   char*& f1,
   char*& f2,
   char*& f3,
   char*& f4,
   char*& f5,
   char*& f6
)
{
   return SPxLP::readLine(is, f1, f2, f3, f4, f5, f6);
}

/**@todo suspicious: readRows only finds the keyword ROWS, 
         readObj reads all rows (which should be done by readRows?)
*/
static MPS_Section readRows(
   std::istream& is,
   LPRowSet& /*rowSet*/,
   NameSet& /*rowNames*/,
   LPColSet& /*colSet*/,
   NameSet& /*colNames*/
)
{
   char* f1;
   char* f2;
   char* f3;
   char* f4;
   char* f5;
   char* f6;
   MPS_Section next = ENDATA;

   do
   {
      if (SPxLP__readLine(is, f1, f2, f3, f4, f5, f6))
      {
         if (strncmp(f1, "ROWS", 4) == 0)
         {
            next = ROWS;
            break;
         }
         else if (strncmp(f1, "NAME", 4) == 0)
            continue;
         else
            std::cerr << "illegal section order in MPS file\n";
      }
   }
   while (is.good());
   return next;
}

static MPS_Section readObj(
   std::istream& is,
   LPRowSet& rowSet,
   NameSet& rowNames,
   LPColSet& /*colSet*/,
   NameSet& /*colNames*/
)
{
   char* f1;
   char* f2;
   char* f3;
   char* f4;
   char* f5;
   char* f6;
   MPS_Section next = ENDATA;

   objName[0] = 0;
   LPRow row;
   do
   {
      if (SPxLP__readLine(is, f1, f2, f3, f4, f5, f6))
      {
         if (strncmp(f1, "COLUMNS", 7) == 0)
         {
            next = COLUMNS;
            break;
         }
         else
            std::cerr << "illegal section order in MPS file\n";
      }
      else
      {
         switch (*f1)
         {
         case 'G':
            assert(!rowNames.has(f2));
            rowNames.add(f2);
            row.lhs() = 0;
            row.rhs() = SPxLP::infinity;
            rowSet.add(row);
            break;
         case 'E':
            assert(!rowNames.has(f2));
            rowNames.add(f2);
            row.lhs() = 0;
            row.rhs() = 0;
            rowSet.add(row);
            break;
         case 'L':
            assert(!rowNames.has(f2));
            rowNames.add(f2);
            row.lhs() = -SPxLP::infinity;
            row.rhs() = 0;
            rowSet.add(row);
            break;
         case 'N':
            if (!objName[0])
               strcpy(objName, f2);
            break;
         default:
            abort();
         }
#ifndef    NDEBUG
         if (*f1 != 'N')
         {
            int num = rowNames.number(f2);
            assert(num == rowSet.num() - 1);
         }
#endif  // NDEBUG
      }
   }
   while (f1);
   return next;
}


static MPS_Section readCols(
   std::istream& is,
   LPRowSet& rowSet,
   NameSet& rowNames,
   LPColSet& colSet,
   NameSet& colNames
)
{
   char* f1;
   char* f2;
   char* f3;
   char* f4;
   char* f5;
   char* f6;
   MPS_Section next = ENDATA;

   double val;
   int idx;
   int first = 1;
   int colNum;
   LPCol col(rowSet.num());
   char colStr[LINE];
   colStr[0] = 0;
   col.obj() = 0;
   col.lower() = 0;
   col.upper() = SPxLP::infinity;
   col.colVector().clear();
   do
   {
      if (SPxLP__readLine(is, f1, f2, f3, f4, f5, f6))
      {
         if (strncmp(f1, "RHS", 3) == 0)
         {
            next = RHS;
            break;
         }
         else
            std::cerr << "illegal section order in MPS file\n";
      }
      if (strcmp(colStr, f2) != 0)
      {
         strcpy(colStr, f2);
         colNames.add(f2);
         if (!first)
            colSet.add(col);
         col.colVector().clear();
         col.obj() = first = 0;
         colNum = colNames.number(f2);
      }
      if (f4)
      {
         if (strcmp(f3, objName) == 0)
            col.obj() = atof(f4);
         else if ((idx = rowNames.number(f3)) >= 0)
         {
            val = atof(f4);
            if (val != 0)
               col.colVector().add(idx, val);
         }
      }
      if (f6)
      {
         if (strcmp(f5, objName) == 0)
            col.obj() = atof(f6);
         else if ((idx = rowNames.number(f5)) >= 0)
         {
            val = atof(f6);
            if (val != 0)
               col.colVector().add(idx, val);
         }
      }
   }
   while (f2);
   if (!first)
      colSet.add(col);
   return next;
}

static MPS_Section readRhs(
   std::istream& is,
   LPRowSet& rowSet,
   NameSet& rowNames,
   LPColSet& /*colSet*/,
   NameSet& /*colNames*/
)
{
   char* f1;
   char* f2;
   char* f3;
   char* f4;
   char* f5;
   char* f6;
   MPS_Section next = ENDATA;
   char rhsStr[LINE];
   int rhsSet = 0;
   rhsStr[0] = 0;

   do
   {
      if (SPxLP__readLine(is, f1, f2, f3, f4, f5, f6))
      {
         if (strncmp(f1, "RANGES", 6) == 0)
         {
            next = RANGES;
            break;
         }
         if (strncmp(f1, "BOUNDS", 6) == 0)
         {
            next = BOUNDS;
            break;
         }
         if (strncmp(f1, "ENDATA", 7) == 0)
         {
            next = ENDATA;
            break;
         }
      }
      if (!rhsSet)
      {
         rhsSet = 1;
         if (f2)
            strcpy(rhsStr, f2);
      }
      if ((f2 && strcmp(f2, rhsStr) == 0)
           || (f2 == 0 && *rhsStr == 0))
      {
         if (f4)
         {
            int idx = rowNames.number(f3);
            if (idx >= 0)
            {
               if (rowSet.rhs(idx) < SPxLP::infinity)
                  rowSet.rhs(idx) = atof(f4);
               if (rowSet.lhs(idx) > -SPxLP::infinity)
                  rowSet.lhs(idx) = atof(f4);
            }
            else
               std::cerr << "ignoring RHS for unknown row " << f3 << std::endl;
         }
         if (f6)
         {
            int idx = rowNames.number(f5);
            if (idx >= 0)
            {
               if (rowSet.rhs(idx) < SPxLP::infinity)
                  rowSet.rhs(idx) = atof(f6);
               if (rowSet.lhs(idx) > -SPxLP::infinity)
                  rowSet.lhs(idx) = atof(f6);
            }
            else
               std::cerr << "ignoring RHS for unknown row " << f5 << std::endl;
         }
      }
   }
   while (f4);
   return next;
}

static MPS_Section readRanges(
   std::istream& is,
   LPRowSet& rowSet,
   NameSet& rowNames,
   LPColSet& /*colSet*/,
   NameSet& /*colNames*/
)
{
   char* f1;
   char* f2;
   char* f3;
   char* f4;
   char* f5;
   char* f6;
   MPS_Section next = ENDATA;

   char rngStr[LINE];
   int rngSet;
   rngSet = 0;
   rngStr[0] = 0;
   do
   {
      if (SPxLP__readLine(is, f1, f2, f3, f4, f5, f6))
      {
         if (strncmp(f1, "BOUNDS", 6) == 0)
         {
            next = BOUNDS;
            break;
         }
         if (strncmp(f1, "ENDATA", 7) == 0)
         {
            next = ENDATA;
            break;
         }
      }
      if (!rngSet)
      {
         rngSet = 1;
         if (f2)
            strcpy(rngStr, f2);
      }
      if ((f2 && strcmp(f2, rngStr) == 0)
           || (f2 == 0 && *rngStr == 0))
      {
         if (f4)
         {
            int idx = rowNames.number(f3);
            if (idx >= 0)
            {
               double val = atof(f4);
               if (val >= 0)
               {
                  if (rowSet.lhs(idx) > -SPxLP::infinity)
                     rowSet.rhs(idx) = rowSet.lhs(idx) + val;
                  else
                     rowSet.lhs(idx) = rowSet.rhs(idx) - val;
               }
               else
               {
                  assert(rowSet.rhs(idx) == rowSet.lhs(idx));
                  rowSet.lhs(idx) += val;
               }
            }
            else
               std::cerr << "ignoring RANGE for unknown row " << f3 << std::endl;
         }
         if (f6)
         {
            int idx = rowNames.number(f5);
            if (idx >= 0)
            {
               double val = atof(f6);
               if (val >= 0)
               {
                  if (rowSet.lhs(idx) > -SPxLP::infinity)
                     rowSet.rhs(idx) = rowSet.lhs(idx) + val;
                  else
                     rowSet.lhs(idx) = rowSet.rhs(idx) - val;
               }
               else
               {
                  assert(rowSet.rhs(idx) == rowSet.lhs(idx));
                  rowSet.lhs(idx) += val;
               }
            }
            else
               std::cerr << "ignoring RANGE for unknown row " << f5 << std::endl;
         }
      }
   }
   while (f4);
   return next;
}

static MPS_Section readBounds(
   std::istream& is,
   LPRowSet& /*rowSet*/,
   NameSet& /*rowNames*/,
   LPColSet& colSet,
   NameSet& colNames
)
{
   char* f1;
   char* f2;
   char* f3;
   char* f4;
   char* f5;
   char* f6;
   MPS_Section next = ENDATA;

   char bndStr[LINE];
   int bndSet;
   bndSet = 0;
   bndStr[0] = 0;

   do
   {
      if (SPxLP__readLine(is, f1, f2, f3, f4, f5, f6))
      {
         if (strncmp(f1, "ENDATA", 6) == 0)
         {
            next = ENDATA;
            break;
         }
      }
      if (!bndSet)
      {
         bndSet = 1;
         if (f2)
            strcpy(bndStr, f2);
      }
      if ((f2 && strcmp(f2, bndStr) == 0)
           || (f2 == 0 && *bndStr == 0))
      {
         int idx = colNames.number(f3);
         switch (*f1)
         {
         case 'L':
            colSet.lower(idx) = atof(f4);
            break;
         case 'U':
            colSet.upper(idx) = atof(f4);
            break;
         case 'F':
            if (f1[1] == 'X')
            {
               colSet.lower(idx) = atof(f4);
               colSet.upper(idx) = atof(f4);
            }
            else
            {
               colSet.lower(idx) = -SPxLP::infinity;
               colSet.upper(idx) = SPxLP::infinity;
            }
            break;
         case 'M':
            colSet.lower(idx) = -SPxLP::infinity;
            break;
         case 'P':
            colSet.upper(idx) = SPxLP::infinity;
            break;
         default:
            std::cerr << "ignoring unknown bound type for column " << f3 << std::endl;
            break;
         }
      }
   }
   while (f1);
   return next;
}

void SPxLP::readMPS(std::istream& is, NameSet* rn, NameSet* cn)
{
   LPRowSet& rowSet = *this; //(static_cast<LPRowSet*>(this));
   LPColSet& colSet = *this; //(static_cast<LPColSet*>(this));
   NameSet* rowNames = rn;
   NameSet* colNames = cn;
   NameSet _rowNames(1000, 10000);
   NameSet _colNames(1000, 10000);

   if (!rowNames)
      rowNames = &_rowNames;
   if (!colNames)
      colNames = &_colNames;
   rowNames->clear();
   colNames->clear();

   colSet.memRemax(10000);
   colSet.reMax(1000);

   clear();

   MPS_Section next;
   line = 0;
   next = readRows(is, rowSet, *rowNames, colSet, *colNames);

   if (next == ROWS)
      next = readObj(is, rowSet, *rowNames, colSet, *colNames);

   addedRows(rowSet.num());

   if (next == COLUMNS)
      next = readCols(is, rowSet, *rowNames, colSet, *colNames);

   if (next == RHS)
      next = readRhs(is, rowSet, *rowNames, colSet, *colNames);

   if (next == RANGES)
      next = readRanges(is, rowSet, *rowNames, colSet, *colNames);

   if (next == BOUNDS)
      next = readBounds(is, rowSet, *rowNames, colSet, *colNames);

   assert(next == ENDATA);
   changeSense(MINIMIZE);

   added2Set(
      *(reinterpret_cast<SVSet*>(static_cast<LPRowSet*>(this))), 
      *(reinterpret_cast<SVSet*>(static_cast<LPColSet*>(this))), 
      colSet.num());
   addedCols(colSet.num());
   assert(isConsistent());
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
