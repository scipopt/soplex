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

/**@file  spxdefines.cpp
 * @brief Debugging, floating point type and parameter definitions.
 */
#include "assert.h"
#include "soplex/spxdefines.h"
#include "soplex/spxout.h"
#include "soplex/rational.h"

namespace soplex
{
// Overloaded EQ function
bool EQ(int a, int b)
{
   return (a == b);
}

SOPLEX_THREADLOCAL const Real infinity                 = SOPLEX_DEFAULT_INFINITY;

bool msginconsistent(const char* name, const char* file, int line)
{
   assert(name != nullptr);
   assert(file != nullptr);
   assert(line >= 0);

   SPX_MSG_ERROR(std::cerr << file << "(" << line << ") "
                 << "Inconsistency detected in " << name << std::endl;)

   return 0;
}

// Note: Tolerances methods are now inline in the templated TolerancesBase<R> class
// defined in spxdefines.h. The typedef `Tolerances = TolerancesBase<Real>` provides
// backward compatibility.

} // namespace soplex
