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
#pragma ident "@(#) $Id: ctimer.h,v 1.1 2001/11/26 15:38:29 bzfbleya Exp $"

/**@file  ctimer.h
 * @brief C-style timing data types and functions.
 */

#ifndef _C_TIMER_H_
#define _C_TIMER_H_

/* import system include files */

/* import user include files */

/* import system include files */
#include <sys/types.h>
#include <time.h>
#include <assert.h>

namespace soplex
{

#ifdef  __cplusplus
extern "C"
{
#endif

   /* define `inline' as appropriate */
#if !defined (__cplusplus)
#if !defined (inline)
#if defined (__GNUC__)
#define inline  __inline__
#else
#define inline
#endif  /* !__GNUC__ */
#endif  /* !inline */
#endif  /* Not C++  */


   /**@name    Timing routines
    * @ingroup Elementary
    * Here we provide c-style data types and functions used for timing.
    */
   /*\@{*/
   /**@brief timer status */
   enum _Timer_Status {
      Timer_RESET,   ///< reset
      Timer_STOPPED, ///< stopped
      Timer_RUNNING  ///< running
   };

   /**@struct _Timer_Struct
      @brief timer data structure */
   struct _Timer_Struct
   {
      enum _Timer_Status status;  ///< timer status
      clock_t uAccount;           ///< user time
      clock_t sAccount;           ///< system time
      clock_t rAccount;           ///< real time
   };


   /// convert ticks to seconds
   extern double Timer_ticks2sec(clock_t ticks);


   /// get actual user, system and real time from system
   extern void Timer_getTicks(clock_t* usrTicks,
                                 clock_t* sysTicks,
                                 clock_t* realTicks);

   /// abstract type of Timer
   typedef struct _Timer_Struct Timer_t;


   /// initialize timer
   /** set timing accounts to zero */
   static inline void Timer_reset(Timer_t *timer)
   {
      assert(timer != 0);
      timer->status = Timer_RESET;
      timer->uAccount = timer->rAccount = timer->sAccount = 0;
   }


   /// start timer
   /** resume accounting user, system and real time */
   static inline void Timer_start(Timer_t *timer)
   {
      assert(timer != 0);
      assert(timer->status == Timer_RESET
             || timer->status == Timer_STOPPED
             || timer->status == Timer_RUNNING);

      /* ignore start request if timer is runnning */
      if (timer->status != Timer_RUNNING)
      {
         clock_t uTime, sTime, rTime;
         Timer_getTicks(&uTime, &sTime, &rTime);
         timer->uAccount -= uTime;
         timer->sAccount -= sTime;
         timer->rAccount -= rTime;
         timer->status = Timer_RUNNING;
      }
   }


   /// stop timer
   /** return accounted user time */
   static inline double Timer_stop(Timer_t *timer)
   {
      assert(timer != 0);
      assert(timer->status == Timer_RESET
             || timer->status == Timer_STOPPED
             || timer->status == Timer_RUNNING);

      /* status remains unchanged if timer is not running */
      if (timer->status == Timer_RUNNING)
      {
         clock_t uTime, sTime, rTime;
         Timer_getTicks(&uTime, &sTime, &rTime);
         timer->uAccount += uTime;
         timer->sAccount += sTime;
         timer->rAccount += rTime;
         timer->status = Timer_STOPPED;
      }

      return (Timer_ticks2sec(timer->uAccount));
   }


   /// get user, system or real time accounted by timer
   /** null pointers for times are allowed.
    */
   static inline void Timer_getTimes(const Timer_t *timer,
                                     double *userTime,
                                     double *systemTime,
                                     double *realTime)
   {
      assert(timer != 0);
      assert(timer->status == Timer_RESET
             || timer->status == Timer_STOPPED
             || timer->status == Timer_RUNNING);

      if (timer->status != Timer_RUNNING)
      {
         if (userTime)
            *userTime = Timer_ticks2sec(timer->uAccount);
         if (systemTime)
            *systemTime = Timer_ticks2sec(timer->sAccount);
         if (realTime)
            *realTime = Timer_ticks2sec(timer->rAccount);
      }
      else
      {
         clock_t uTime, sTime, rTime;
         Timer_getTicks(&uTime, &sTime, &rTime);
         if (userTime)
            *userTime = Timer_ticks2sec(uTime + timer->uAccount);
         if (systemTime)
            *systemTime = Timer_ticks2sec(sTime + timer->sAccount);
         if (realTime)
            *realTime = Timer_ticks2sec(rTime + timer->rAccount);
      }

      assert(userTime ? *userTime >= 0.0 : 1);
      assert(systemTime ? *systemTime >= 0.0 : 1);
      assert(realTime ? *realTime >= 0.0 : 1);
   }


   /// return user time accounted by timer
   static inline double Timer_userTime(const Timer_t *timer)
   {
      double uTime;
      assert(timer != 0);
      assert(timer->status == Timer_RESET
             || timer->status == Timer_STOPPED
             || timer->status == Timer_RUNNING);

      Timer_getTimes(timer, &uTime, 0, 0);
      return (uTime);
   }


   /// return system time accounted by timer
   static inline double Timer_systemTime(const Timer_t *timer)
   {
      double sTime;
      assert(timer != 0);
      assert(timer->status == Timer_RESET
             || timer->status == Timer_STOPPED
             || timer->status == Timer_RUNNING);

      Timer_getTimes(timer, 0, &sTime, 0);
      return (sTime);
   }


   /// return real time accounted by timer
   static inline double Timer_realTime(const Timer_t *timer)
   {
      double rTime;
      assert(timer != 0);
      assert(timer->status == Timer_RESET
             || timer->status == Timer_STOPPED
             || timer->status == Timer_RUNNING);

      Timer_getTimes(timer, 0, 0, &rTime);
      return (rTime);
   }


   /// return resolution of timer as 1/seconds
   extern long Timer_resolution(void);

   /*\@}*/

#ifdef  __cplusplus
}
#endif


} // namespace soplex
#endif // _TIMER_H_

//-----------------------------------------------------------------------------
//Emacs Local Variables:
//Emacs mode:c++
//Emacs c-basic-offset:3
//Emacs tab-width:8
//Emacs indent-tabs-mode:nil
//Emacs End:
//-----------------------------------------------------------------------------
