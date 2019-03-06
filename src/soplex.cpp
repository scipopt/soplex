/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the class library                   */
/*       SoPlex --- the Sequential object-oriented simPlex.                  */
/*                                                                           */
/*    Copyright (C) 1996-2018 Konrad-Zuse-Zentrum                            */
/*                            fuer Informationstechnik Berlin                */
/*                                                                           */
/*  SoPlex is distributed under the terms of the ZIB Academic Licence.       */
/*                                                                           */
/*  You should have received a copy of the ZIB Academic License              */
/*  along with SoPlex; see the file COPYING. If not email to soplex@zib.de.  */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/**@file  soplex.cpp
 * @brief Preconfigured SoPlex LP solver
 */

#include <assert.h>
#include "limits.h"
#include <iostream>

#ifndef _MSC_VER
#include <strings.h>
#endif

#include "soplex.h"
#include "soplex/spxfileio.h"
#include "soplex/statistics.h"
#include "soplex/mpsinput.h"

#ifdef _MSC_VER
#define strncasecmp _strnicmp
#endif

namespace soplex
{

  // template <class R>
	// bool SoPlexBase<R>::saveSettingsFile(const char* filename, const bool onlyChanged) const;

  // template <class R>
  // void SoPlexBase<R>::printShortStatistics(std::ostream& os);

  // template <class R>
  // SoPlexBase<R>::Settings::BoolParam::BoolParam();


  // template <class R>
  // R SoPlexBase<R>::realParam(const RealParam param) const;

  // template <class R>
  // SoPlexBase<R>::Settings::Settings();






} // namespace soplex
