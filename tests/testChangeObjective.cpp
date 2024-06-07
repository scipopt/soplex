/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the class library                   */
/*       SoPlex --- the Sequential object-oriented simPlex.                  */
/*                                                                           */
/*  Copyright (c) 1996-2024 Zuse Institute Berlin (ZIB)                      */
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
