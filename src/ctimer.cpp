/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the class library                   */
/*       SoPlex --- the Sequential object-oriented simPlex.                  */
/*                                                                           */
/*    Copyright (C) 1997-1999 Roland Wunderling                              */
/*                  1997-2001 Konrad-Zuse-Zentrum                            */
/*                            fuer Informationstechnik Berlin                */
/*                                                                           */
/*  SoPlex is distributed under the terms of the ZIB Academic Licence.       */
/*                                                                           */
/*  You should have received a copy of the ZIB Academic License              */
/*  along with SoPlex; see the file COPYING. If not email to soplex@zib.de.  */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
#pragma ident "@(#) $Id: ctimer.cpp,v 1.1 2001/11/26 15:38:29 bzfbleya Exp $"



/* import system includes */
#include <assert.h>

#ifdef _WIN32

#include <time.h>

#else   /* !_WIN32 */

#include <sys/times.h>
#include <limits.h>
#ifdef HAVE_UNISTD_H
#include <unistd.h>
#endif /* HAVE_UNISTD_H */
#include <sys/param.h>

#endif  /* !_WIN32 */

/* import interface of timer */
#include "timer.h"

namespace soplex
{


/* determine TIMES_TICKS_PER_SEC for clock ticks delivered by times().
 * (don't use CLOCKS_PER_SEC since this is related to clock() only).
 */
#if defined(CLK_TCK)
# define TIMES_TICKS_PER_SEC CLK_TCK
#elif defined(_SC_CLK_TCK)
# define TIMES_TICKS_PER_SEC sysconf(_SC_CLK_TCK)
#elif defined(HZ)
# define TIMES_TICKS_PER_SEC HZ
#else /* !CLK_TCK && !_SC_CLK_TCK && !HZ */
# define TIMES_TICKS_PER_SEC 60
#endif /* !CLK_TCK && !_SC_CLK_TCK && !HZ */


/******************************************************************************
 *
 *      Implementation of Module |Timer|
 *
 ******************************************************************************
 */


/*
 * Private Functions:
 */


/* convert ticks to seconds */
double Timer_ticks2sec(clock_t ticks)
{
   clock_t sec, msec;

   assert(TIMES_TICKS_PER_SEC > 0);
   sec = ticks / TIMES_TICKS_PER_SEC;
   msec = ((ticks % TIMES_TICKS_PER_SEC) * 1000) / TIMES_TICKS_PER_SEC;
   return (static_cast<double>(sec) + (static_cast<double>(msec) / 1000.0));

   /* another way to compute seconds ...
      return (((double) ((ticks * 1000) / TIMES_TICKS_PER_SEC)) / 1000.0); */
}


/* get actual user, system and real time from system */
void Timer_getTicks(clock_t* usrTicks,
                    clock_t* sysTicks,
                    clock_t* realTicks)
{
#ifdef  _WIN32

   *usrTicks = *sysTicks = *realTicks = clock();

#else   /* !_WIN32 */

   struct tms now;
   clock_t ret = times(&now);

   /* ignore errors returned by times() */
   if (ret == -1)
      now.tms_utime = now.tms_stime = ret = 0;

   if (usrTicks)
      *usrTicks = now.tms_utime;

   if (sysTicks)
      *sysTicks = now.tms_stime;

   if (realTicks)
      *realTicks = ret;

#endif  /* !_WIN32 */
}


/*
 * Public Functions:
 */


/* return resolution of timer as 1/seconds */
long Timer_resolution(void)
{
   assert(TIMES_TICKS_PER_SEC > 0);
   return static_cast<long>(TIMES_TICKS_PER_SEC);
}


/******************************************************************************
 *
 *      End of Implementation of Module |Timer|
 *
 ******************************************************************************
 */
} // namespace soplex

//-----------------------------------------------------------------------------
//Emacs Local Variables:
//Emacs mode:c++
//Emacs c-basic-offset:3
//Emacs tab-width:8
//Emacs indent-tabs-mode:nil
//Emacs End:
//-----------------------------------------------------------------------------
