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
#pragma ident "@(#) $Id: timer.h,v 1.5 2001/11/26 13:52:23 bzfbleya Exp $"

/**@file  timer.h
 * @brief Timer classes and functions.
 */

#ifndef _TIMER_H_
#define _TIMER_H_

/* import system include files */

/* import user include files */

/* import system include files */
#include <sys/types.h>
#include <time.h>

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
   //@{
   /**@brief timer status */
   enum _Timer_Status {
      Timer_RESET,   ///< reset
      Timer_STOPPED, ///< stopped
      Timer_RUNNING  ///< running
   };

   /**@brief timer data structure */
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

   //@}

#ifdef  __cplusplus
}
#endif


   /**@name    Timer
      @ingroup Elementary

    In C or C++ programs, the usual way to measure time intervalls, e.g.
    running times of some complex computations, is to call one of the provided
    system functions like %clock(), %time(), %times(), %gettimeofday(), %getrusage() etc.  
    By these functions one can gather
    information about the process' user and system time and the system clock
    (real time).
 
    Unfortunately, these functions are rather clumsy.  The programmer determines
    computation times by querying a (virtual) clock value at the beginning and
    another one at the end of some computation and converting the difference of
    these values into seconds.  Some functions impose some restrictions, for
    instance, the values of the ANSI C function %clock() are of high
    resolution but will wrap around after about 36 minutes (cpu time).  Most
    timing functions take some data structure as argument that has to be
    allocated before the call and from which the user has to pick up the
    information of interest after the call.  Problems can arise, when porting
    programs to other operating systems that do not support standards like POSIX
    etc.
 
    In order to simplify measuring computation times and to hide the
    system-dependencies involved, a concept of \em timers accounting the
    process' user, system and real time is implemented.  C and C++ interfaces
    are provided as a set of functions operating on timers and a timer class
    respectively.
 
    The idea is to provide a type Timer for objects that act like a stopwatch.
    Operations on such an objects include: start accounting time, stop
    accounting, read the actual time account and reset the objects time account
    to zero.
 
    After initialization, accounting for user, system and real time can be
    started by calling a function start(). Accounting is suspended by calling
    a function stop() and can be resumed at any time by calling start()
    again.
 
    The user, system or real time actually accounted by a timer can be accessed
    at any time by the methods shown in this code section:
 
    \verbatim
       double utime, stime, rtime;
         
       utime = timer.userTime();
       stime = timer.systemTime();
       rtime = timer.realTime();
         
       timer.getTimes(utime, stime rtime);
    \endverbatim
 
    For convenience, the actually accounted user time is returned by stop()
    too.  Function reset() re-initializes a timer clearing all time
    accounts.
 
    Function resolution() returns the smallest (non-zero) time intervall which
    is resolved by the underlying system function: res = 1/Timer_resolution().


    The process' user and system times are accessed by calling function %times(),
    which is declared in \c <sys/times.h>.  If OS supports POSIX
    compatibility through providing \c <sys/unistd.h>, set
    \c -DHAVE_UNISTD_H when compiling \c timer.c.  Ignore compiler
    warnings about missing prototypes of functions.
 */
class Timer
{
public:

   /// initialize timer, set timing accounts to zero.
   void reset();

   /// start timer, resume accounting user, system and real time.
   void start();

   /// stop timer, return accounted user time.
   double stop();

   /// return accounted user time.
   double userTime() const;

   /// return accounted system time.
   double systemTime() const;

   /// return accounted real time.
   double realTime() const;

   /// get accounted user, system or real time.
   void getTimes(double *userTime = 0,
                 double *systemTime = 0,
                 double *realTime = 0) const;

   /// return resolution of timer as 1/seconds.
   static long resolution();

   // Constructors and Destructor:

   ///
   Timer();

protected:
   // relay on C implementation of timer
   Timer_t timer;
};


// initialize timer, set time accounts to zero
inline void Timer::reset()
{
   Timer_reset(&timer);
}


// start timer, resume accounting user, system and real time
inline void Timer::start()
{
   Timer_start(&timer);
}


// stop timer, return accounted user time
inline double Timer::stop()
{
   return Timer_stop(&timer);
}


// return accounted user time
inline double Timer::userTime() const
{
   return Timer_userTime(&timer);
}


// return accounted system time
inline double Timer::systemTime() const
{
   return Timer_systemTime(&timer);
}


// return accounted real time
inline double Timer::realTime() const
{
   return Timer_realTime(&timer);
}


// get accounted user, system or real time
inline void Timer::getTimes(double *puserTime,
                            double *psystemTime,
                            double *prealTime) const
{
   Timer_getTimes(&timer, puserTime, psystemTime, prealTime);
}


// return resolution of timer as 1/seconds
inline long Timer::resolution()
{
   return Timer_resolution();
}


// default constructor (initializing timer)
inline Timer::Timer()
{
   Timer_reset(&timer);
}


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
