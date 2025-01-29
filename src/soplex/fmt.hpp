/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the class library                   */
/*       SoPlex --- the Sequential object-oriented simPlex.                  */
/*                                                                           */
/*  Copyright 1996-2025 Zuse Institute Berlin                                */
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

#ifndef _SOPLEX_FMT_HPP_
#define _SOPLEX_FMT_HPP_

#ifndef FMT_HEADER_ONLY
#define FMT_HEADER_ONLY
#endif

/* if those macros are not defined and fmt includes windows.h
 * then many macros are defined that can interfere with standard C++ code
 */
#ifndef NOMINMAX
#define NOMINMAX
#define SOPLEX_DEFINED_NOMINMAX
#endif

#ifndef WIN32_LEAN_AND_MEAN
#define WIN32_LEAN_AND_MEAN
#define SOPLEX_DEFINED_WIN32_LEAN_AND_MEAN
#endif

#ifndef NOGDI
#define NOGDI
#define SOPLEX_DEFINED_NOGDI
#endif

/* to provide backwards compatibility use fmt of PaPILO <= 2.3 due to breaking changes in fmt 7 */
#if !defined(SOPLEX_WITH_PAPILO) || PAPILO_VERSION_MAJOR > 2 || (PAPILO_VERSION_MAJOR == 2 && PAPILO_VERSION_MINOR > 3)
#include "soplex/external/fmt/format.h"
#include "soplex/external/fmt/ostream.h"
#else
#include "papilo/external/fmt/format.h"
#include "papilo/external/fmt/ostream.h"
#endif

#ifdef SOPLEX_DEFINED_NOGDI
#undef NOGDI
#undef SOPLEX_DEFINED_NOGDI
#endif

#ifdef SOPLEX_DEFINED_NOMINMAX
#undef NOMINMAX
#undef SOPLEX_DEFINED_NOMINMAX
#endif

#ifdef SOPLEX_DEFINED_WIN32_LEAN_AND_MEAN
#undef WIN32_LEAN_AND_MEAN
#undef SOPLEX_DEFINED_WIN32_LEAN_AND_MEAN
#endif

#endif
