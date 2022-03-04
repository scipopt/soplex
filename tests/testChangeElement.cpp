/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the class library                   */
/*       SoPlex --- the Sequential object-oriented simPlex.                  */
/*                                                                           */
/*    Copyright (C) 1996-2022 Konrad-Zuse-Zentrum                            */
/*                            fuer Informationstechnik Berlin                */
/*                                                                           */
/*  SoPlex is distributed under the terms of the ZIB Academic Licence.       */
/*                                                                           */
/*  You should have received a copy of the ZIB Academic License              */
/*  along with SoPlex; see the file COPYING. If not email to soplex@zib.de.  */
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
