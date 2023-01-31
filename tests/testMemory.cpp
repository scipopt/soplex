/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the class library                   */
/*       SoPlex --- the Sequential object-oriented simPlex.                  */
/*                                                                           */
/*  Copyright (c) 1996-2023 Zuse Institute Berlin (ZIB)                      */
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

/* Test a for a bug in the memory allocation, which results in a dramatic increase of the used
 * memory - at least to date (May 2009).
 *
 * The reason is in svset.
 */

#include "spxsolver.h"
#include "slufactor.h"
#include "spxfastrt.h"
#include "spxsteeppr.h"

using namespace soplex;

#define NCOLS 100
#define MAXNROWS 1500

int main(int argc, const char* const argv[])
{
   SPxSolver work(SPxSolver::LEAVE, SPxSolver::COLUMN);

   /* no starter, no simplifier, no scaler */
   SLUFactor solver;
   SPxFastRT tester;
   SPxSteepPR pricer;
   work.setBasisSolver(&solver);
   work.setTester(&tester);
   work.setPricer(&pricer);
   work.changeSense(work.MINIMIZE);

   /* create columns */
   for (int j = 0; j < NCOLS; ++j)
      work.addCol( LPCol(1.0, DSVector(), 1.0, 0.0) );

   /* now add rows in rounds */
   for (int i = 0; i < MAXNROWS; ++i)
   {
      DSVector row;

      /* create row with random coefficients, but dense */
      for (int j = 0; j < NCOLS; ++j)
	 row.add(j, drand48());

      work.addRow( LPRow(row, LPRow::GREATER_EQUAL, 1.0) );
      assert( work.isConsistent() );

      /* solve the problem */
      // Param::setVerbose(5);

      /* uncomment to yield faster way to desaster: */
      work.solve();

      std::cout << "n = " << NCOLS << " m = " << i << "   value: " << std::setw(10) << work.value() << "   iters: " << work.iterations() << std::endl;
   }

   // work.writeLPF(std::cout, NULL, NULL);
   // std::cout << std::endl;

   return 0;
}
