/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the class library                   */
/*       SoPlex --- the Sequential object-oriented simPlex.                  */
/*                                                                           */
/*  Copyright (c) 1996-2026 Zuse Institute Berlin (ZIB)                      */
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

/* Test a for a bug in changesoplex.cpp (version 1.3.3), which produces an assert if columns or rows
 *  are added and then changeElement is called.
 */

#include "spxsolver.h"
#include "slufactor.h"
#include "spxfastrt.h"
#include "spxsteeppr.h"

using namespace soplex;

int main(int argc, const char* const argv[])
{
   SPxSolver work(SPxSolver::LEAVE, SPxSolver::COLUMN);

   /* no starter, no simplifier, no scaler */
   SLUFactor<Real> solver;
   SPxFastRT tester;
   SPxSteepPR pricer;
   work.setBasisSolver(&solver);
   work.setTester(&tester);
   work.setPricer(&pricer);
   work.changeSense(work.MAXIMIZE);

   /* create some columns */
   work.addCol( LPCol(1.0, DSVector(), 4.0, 0.0) );
   work.addCol( LPCol(1.0, DSVector(), 3.0, 0.0) );
   work.addCol( LPCol(1.0, DSVector(), 2.0, 0.0) );
   work.addCol( LPCol(1.0, DSVector(), 1.0, 0.0) );

   /* add rows */
   DSVector row;
   row.add(0, 1.0);
   row.add(1, 1.0);
   row.add(2, 1.0);
   row.add(3, 1.0);
   work.addRow( LPRow(row, LPRow::EQUAL, 1.0) );

   assert( work.isConsistent() );
   std::cout << "Initial problem:" << std::endl << std::endl;
   work.writeLPF(std::cout, NULL, NULL);
   std::cout << std::endl;

   work.changeElement(work.rowId(0), work.colId(0), 2.0);

   assert( work.isConsistent() );

   std::cout << "Problem after changing element (1,1):" << std::endl << std::endl;
   work.writeLPF(std::cout, NULL, NULL);
   std::cout << std::endl;

   return 0;
}
