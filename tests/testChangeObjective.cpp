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

/* Test a for a bug in soplex (version 1.4), which avoids the change of the objective. */

#include <sstream>
#include <soplex.h>

using namespace soplex;

const char lpstr[] =
  "Maximize               \n"
  " obj: x0 + x1          \n"
  "Subject To             \n"
  " c1: x0 + x1 <= 1      \n"
  " c2: x0 + x1 >= -1     \n"
  " c3: x0 - x1 <= 1      \n"
  " c4: x0 - x1 >= -1     \n"
  "Bounds                 \n"
  " 0 <= x0               \n"
  " 0 <= x1               \n"
  "End                    \n";

int main()
{
   SoPlex lp;
   std::istringstream is(lpstr);
   lp.readLPF(is);

   lp.optimize();
   std::cerr << lp.objValue() << std::endl;

   lp.changeSense(SoPlex::MINIMIZE);
   lp.optimize();
   std::cerr << lp.objValue() << std::endl;

   return 0;
}
