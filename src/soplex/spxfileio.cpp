/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the class library                   */
/*       SoPlex --- the Sequential object-oriented simPlex.                  */
/*                                                                           */
/*    Copyright (C) 1996-2017 Konrad-Zuse-Zentrum                            */
/*                            fuer Informationstechnik Berlin                */
/*                                                                           */
/*  SoPlex is distributed under the terms of the ZIB Academic Licence.       */
/*                                                                           */
/*  You should have received a copy of the ZIB Academic License              */
/*  along with SoPlex; see the file COPYING. If not email to soplex@zib.de.  */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

#include <assert.h>

#include "soplex/spxdefines.h"
#include "soplex/spxsolver.h"
#include "soplex/spxfileio.h"

namespace soplex
{
  template <>
  bool SPxSolver<Real>::readBasisFile(
                                   const char*    filename, 
                                   const NameSet* rowNames,
                                   const NameSet* colNames)
  {

    spxifstream file(filename);

    if (!file)
      return false;
 
    return this->readBasis(file, rowNames, colNames);
  }

  template <>
  bool SPxSolver<Real>::writeBasisFile
  ( const char*    filename, 
    const NameSet* rowNames,
    const NameSet* colNames,
    const bool cpxFormat ) const
  {
    std::ofstream file(filename);

    if (!file)
      return false;
 
    this->writeBasis(file, rowNames, colNames);

    return true;
  }

} // namespace soplex
