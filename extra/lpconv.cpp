/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the class library                   */
/*       SoPlex --- the Sequential object-oriented simPlex.                  */
/*                                                                           */
/*  Copyright (c) 1996-2025 Zuse Institute Berlin (ZIB)                      */
/*                                                                           */
/*  Licensed under the Apache License, Version 2.0 (the "License");          */
/*  you may not use this file except in compliance with the License.         */
/*  You may obtain a copy of the License at                                  */
/*                                                                           */
/*      http://www.apache.org/licenses/LICENSE-2.0                           */
/*                                                                           */
/*  Unless required by applicable law or agreed to in writing, software      */
/*  distributed under the License is distributed on an "AS IS" BASIS,        */
/*  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied. */
/*  See the License for the specific language governing permissions and      */
/*  limitations under the License.                                           */
/*                                                                           */
/*  You should have received a copy of the Apache-2.0 license                */
/*  along with SoPlex; see the file LICENSE. If not email to soplex@zib.de.  */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

#include <assert.h>
#include <iostream>
#include <fstream>

#include "spxdefines.h"
#include "spxlp.h"
#include "spxfileio.h"

using namespace soplex;

// write constraint line in latte format
//
static void write_row(
   const SVector& row,
   int            cols,
   double         rhs,
   double         dir,
   std::ofstream& ofile)
{
   int idx = 0;

   ofile << rhs * dir;

   //
   for(int k = 0; k < row.size(); ++k)
   {
      for(; idx < row.index(k); idx++)
         ofile << " 0";
      assert(idx == row.index(k));
      ofile << " " << -row.value(k) * dir;
      idx++;
   }
   while(idx < cols)
   {
      ofile << " 0";
      idx++;
   }
   ofile << std::endl;
}

// Write Latte format: http://www.math.ucdavis.edu/~latte/manual/node8.html
//
static void write_latte(const SPxLP& lp, std::ofstream& ofile)
{
   int lowers  = 0;
   int uppers  = 0;
   int nonnegs = 0;
   int equals  = 0;
   int ranges  = 0;

   for(int i = 0; i < lp.nCols(); ++i)
   {
      if (isZero(lp.lower(i)))
         ++nonnegs;
      else if (lp.lower(i) > -infinity)
         ++lowers;

      if (lp.upper(i) < infinity)
         ++uppers;
   }
   for(int i = 0; i < lp.nRows(); ++i)
   {
      switch(lp.rowType(i))
      {
      case LPRow::EQUAL :
         equals++;
         break;
      case LPRow::RANGE :
         ranges++;
         break;
      default :
         break;
      }
   }
   ofile << lp.nRows() + ranges + lowers + uppers << " " << lp.nCols() << std::endl;

   for(int i = 0; i < lp.nRows(); ++i)
   {
      switch(lp.rowType(i))
      {
      case LPRow::EQUAL :
         ++equals;
         write_row(lp.rowVector(i), lp.nCols(), lp.rhs(i), 1, ofile);
         break;
      case LPRow::LESS_EQUAL :
         write_row(lp.rowVector(i), lp.nCols(), lp.rhs(i), 1, ofile);
         break;
      case LPRow::GREATER_EQUAL :
         write_row(lp.rowVector(i), lp.nCols(), lp.lhs(i), -1, ofile);
         break;
      case LPRow::RANGE :
         write_row(lp.rowVector(i), lp.nCols(), lp.lhs(i), -1, ofile);
         write_row(lp.rowVector(i), lp.nCols(), lp.rhs(i),  1, ofile);
         break;
      default:
         abort();
      }
   }
   for(int i = 0; i < lp.nCols(); i++)
   {
      if (isNotZero(lp.lower(i)) && lp.lower(i) > -infinity)
      {
         ofile << -lp.lower(i);

         for(int k = 0; k < lp.nCols(); ++k)
            ofile << " " << (k == i) ? 1 : 0;

         ofile << std::endl;
      }

      if (lp.upper(i) < infinity)
      {
         ofile << lp.upper(i);

         for(int k = 0; k < lp.nCols(); ++k)
            ofile << " " << (k == i) ? -1 : 0;

         ofile << std::endl;
      }
   }

   // Write linearity line if neccessary
   if (equals > 0)
   {
      ofile << "linearity " << equals;

      for(int i = 0; i < lp.nRows(); ++i)
         if (lp.rowType(i) == LPRow::EQUAL)
            ofile << " " << i + 1;

      ofile << std::endl;
   }

   // Write nonnegative line if neccessary
   if (nonnegs > 0)
   {
      ofile << "nonnegative " << nonnegs;

      for(int i = 0; i < lp.nCols(); ++i)
         if (isZero(lp.lower(i)))
            ofile << " " << i + 1;

      ofile << std::endl;
   }
}

static void read_latte(
   SPxLP& lp,
   spxifstream& ifile,
   NameSet& rownames,
   NameSet& colnames,
   DIdxSet& intvars)
{
   int rows;
   int cols;

   ifile >> rows;
   ifile >> cols;

   std::cout << "Reading Latte file with " << rows << " rows and " << cols << " columns\n";

   for(int r = 0; r < rows; ++r)
   {
      DSVector vec(cols);
      double  rhs;

      ifile >> rhs;

      for(int c = 0; c < cols; ++c)
      {
         double val;

         ifile >> val;

         if (isNotZero(val))
            vec.add(c, -val);
      }
      lp.addRow(LPRow(vec, LPRow::LESS_EQUAL, rhs));
   }
   char keyword[256];

   ifile >> keyword;

   if (!strcmp(keyword, "linearity"))
   {
      std::cout << "Found linearity line\n";

      int count;

      ifile >> count;

      for(int n = 0; n < count; ++n)
      {
         int c;

         ifile >> c;

         lp.changeLhs(c - 1, lp.rhs(c - 1));
      }
   }
   for(int c = 0; c < cols; ++c)
      lp.changeLower(c, -infinity);

   if (!strcmp(keyword, "nonnegative"))
   {
      std::cout << "Found nonnegative line\n";

      int count;

      ifile >> count;

      for(int n = 0; n < count; ++n)
      {
         int c;

         ifile >> c;

         lp.changeLower(c - 1, 0.0);
      }
   }
   for(int c = 0; c < cols; ++c)
   {
      char tmp[32];
      sprintf(tmp, "x%d", c + 1);
      colnames.add(tmp);
      intvars.add(c);
   }
   for(int r = 0; r < rows; ++r)
   {
      char tmp[32];
      sprintf(tmp, "c%d", r + 1);
      rownames.add(tmp);
   }
}

int main(int argc, char **argv)
{
   const char* banner =
   "****************************************************************************\n"
   "*                                                                          *\n"
   "*       LPConv --- Convert LPF to MPS format.                              *\n"
   "*                  Release 1.0.2                                           *\n"
   "*    Copyright (c) 2007-2020 Zuse Institute Berlin (ZIB)                   *\n"
   "*                            fuer Informationstechnik Berlin               *\n"
   "*                                                                          *\n"
   "*  LPConv is distributed under the terms of the Apache 2.0 Licence.        *\n"
   "*  You should have received a copy of the Apache-2.0 license               *\n"
   "*  along with SoPlex; see the file LICENSE. If not email to soplex@zib.de. *\n"
   "*                                                                          *\n"
   "****************************************************************************\n"
   ;

   const char* usage =
   "[options] input-file output-file\n\n"
   "          input-file can be either in MPS or LPF format\n\n"
   "options:  (*) indicates default\n"
   " -i        Latte input format\n"
   " -l        Latte output format\n"
   " -m        MPS output format\n"
   " -vLevel   set verbosity Level [0-3], default 1\n"
   " -V        show program version\n"
   " -h        show this help\n"
   ;

   enum { LPF, MPS, LATTE } output = LPF;
   bool latte_input = false;
   int  verbose = 1;
   int  optidx;

   for(optidx = 1; optidx < argc; optidx++)
   {
      if (*argv[optidx] != '-')
         break;

      switch(argv[optidx][1])
      {
      case 'i' :
         latte_input = true;
         break;
      case 'l' :
         output = LATTE;
         break;
      case 'm' :
         output = MPS;
         break;
      case 'v' :
         verbose = atoi(&argv[optidx][2]);
         break;
      case 'V' :
         std::cout << banner << std::endl;
         exit(0);
      case 'h' :
      case '?' :
         std::cout << banner << std::endl;
         /*FALLTHROUGH*/
      default :
         std::cerr << "usage: " << argv[0] << " " << usage << std::endl;
         exit(0);
      }
   }
   if ((argc - optidx) < 2)
   {
      std::cerr << "usage: " << argv[0] << " " << usage << std::endl;
      exit(0);
   }
   const char* inpfile  = argv[optidx];
   const char* outfile  = argv[optidx + 1];

   SPxLP       lp;
   NameSet     rownames;
   NameSet     colnames;
   DIdxSet     intvars;

   spxifstream ifile(inpfile);

   if (!ifile)
   {
      std::cerr << "Can't open file: " << inpfile << std::endl;
      exit(1);
   }

   if (latte_input)
      read_latte(lp, ifile, rownames, colnames, intvars);
   else
   {
      if (!lp.read(ifile, &rownames, &colnames, &intvars))
      {
         std::cerr << "Error while reading file: " << inpfile << std::endl;
         exit(1);
      }
   }

   std::ofstream ofile(outfile);

   if (!ofile)
   {
      std::cerr << "Can't open file: " << outfile << std::endl;
      exit(1);
   }
   switch(output)
   {
   case MPS :
      lp.writeMPS(ofile, &rownames, &colnames, &intvars);
      break;
   case LPF :
      lp.writeLPF(ofile, &rownames, &colnames, &intvars);
      break;
   case LATTE :
      write_latte(lp, ofile);
      break;
   default :
      abort();
   }
   return 0;
}
