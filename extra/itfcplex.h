/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the class library                   */
/*       SoPlex --- the Sequential object-oriented simPlex.                  */
/*                                                                           */
/*  Copyright 1996-2022 Zuse Institute Berlin                                */
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

/**@file  itfcplex.h
 * @brief Simple CPLEX Interface to SoPlex
 *
 *  Some of the CPLEX interface functions are provided by this module. This
 *  makes it possible to use SoPlex as a replacement for CPLEX, when only a
 *  subset of its API is used.
 */
#ifndef _ITFCPLEX_H_
#define _ITFCPLEX_H_

    /** what to use as $\infty$ */
#define INFBOUND 1e+100 /* SoPlex::infinity ; */

    /** optimality status */
#define  CPX_OPTIMAL 1
    /** infeasibility status */
#define  CPX_INFEASIBLE 2
    /** unboundedness status */
#define  CPX_UNBOUNDED 3
    /** objectiv limit status */
#define  CPX_OBJ_LIM    4

#define  CPXERR_PRESLV_INForUNBD 1101
#define  CPX_DPRIIND_FULLSTEEP      3
#define  CPXALG_PRIMAL              1
#define  CPXALG_DUAL                2


#ifdef FOR_DOCXX
    /** Faked CPLEX pointer.
     */
typedef struct SoPlex* CPXLPptr ;
#endif

#ifndef __cplusplus
struct SoPlex ;
typedef struct SoPlex* CPXLPptr ;
#else
class SoPlex ;
typedef SoPlex* CPXLPptr ;
#endif

/* Function Declarations
 */
/** */
int lpread  (char *, int *, int *, int *, double **, double **,
             char **, int **, int **, int **, double **, double **,
             double **, char **, char **, char ***, char **,
             char ***, char **, int *, int *, int *, unsigned *,
             unsigned *);

/** */
int lpmread (char *, int *, int *, int *, double **, double **,
             char **, int **, int **, int **, double **, double **,
             double **, char **, char **, char ***, char **,
             char ***, char **, int *, int *, int *, unsigned *,
             unsigned *, char** );


/** */
CPXLPptr loadprob (char *, int, int, int, int, double *, double *,
                   char *, int *, int *, int *, double *, double *,
                   double *, double *, int *, int *, int *, int *,
                   int *, double *, char *, char *, char *, char *,
                   char *, char **, char *, char **, char *, char **,
                   char *, int, int, int, int, int, unsigned, unsigned,
                   unsigned);

/** */
CPXLPptr loadlp(char   *probname,
        int    mac,
        int    mar,
        int    objsen,
        double *objx,
        double *rhsx,
        char   *senx,
        int    *matbeg,
        int    *matcnt,
        int    *matind,
        double *matval,
        double *bdl,
        double *bdu,
        double *rngval,
        int    macsz,
        int    marsz,
        int    matsz) ;

/* @ManMemo:
void unloadprob (CPXLPptr*);
*/

/** */
void freeprob   (CPXLPptr*);

/** */
int setmlim  (int, int *, int *);
/** */
int setnlim  (int, int *, int *);
/** */
int setnzlim (int, int *, int *);
/** */
void getmlim  (int *i) ;
/** */
void getnlim  (int *i) ;
/** */
void getnzlim (int *i) ;
/** */
int openCPLEX ( void )  ;
/** */
int closeCPLEX( void ) ;
/** */
int setlogfile(FILE *fp) ;

/** */
int setdpriind  (int, int *, int *);

/** */
int dualopt (CPXLPptr);

/** */
int optimize (CPXLPptr);

/** */
int solution  (CPXLPptr, int *, double *, double *, double *,
               double *, double *);

/** */
int getrows  (CPXLPptr, int *, int *, int *, double *, int,
              int *, int, int);
/** */
int getmac   (CPXLPptr);
/** */
int getmar   (CPXLPptr);
/** */
int getmat   (CPXLPptr);
/** */
int getbase  (CPXLPptr, int *, int *);
/** */
int getsense (CPXLPptr, char *, int, int);
/** */
int getobj   (CPXLPptr, double *, int, int);
/** */
int getrhs   (CPXLPptr, double *, int, int);
/** */
int getbdl   (CPXLPptr, double *, int, int);
/** */
int getbdu   (CPXLPptr, double *, int, int);

/** */
int chgbds (CPXLPptr, int, int *, char *, double *);

/** */
int chgobj (CPXLPptr, int, int *, double *);

/** */
int addrows (CPXLPptr, int, int, int, double *, char *,
             int *, int *, double *, char **, char **);
/** */
int delrows (CPXLPptr cplex, int begin, int end ) ;

/** */
int addcols (CPXLPptr, int, int, double*,
      int*, int*, double*, double*, double*, char**) ;
/** */
int delcols (CPXLPptr cplex, int begin, int end ) ;

/** */
int delsetrows (CPXLPptr, int *);
/** */
int delsetcols (CPXLPptr, int *);

/** */
int lpwrite   (CPXLPptr, char *) ;

int savread (char *, int *, int *, int *, int *, double **,
             double **, char **, int **, int **, int **,
             double **, double **, double **, double **, int **,
             int **, int **, int **, int **, double **,
             char **, char **, char **, char **, char **,
             char ***, char **, char ***, char **, char ***,
             char **, int *, int *, int *, int *, int *,
             unsigned *, unsigned *, unsigned *, int **, int **) ;

/** */
int mpsread (char *, int *, int *, int *, int *, double **,
             double **, char **, int **, int **, int **,
             double **, double **, double **, double **, int **,
             int **, int **, int **, int **, double **,
             char **, char **, char **, char **, char **,
             char ***, char **, char ***, char **, char ***,
             char **, int *, int *, int *, int *, int *,
             unsigned *, unsigned *, unsigned *) ;

/** */
int loadbase(CPXLPptr, int *, int *) ;

/** */
int setobjulim  (double, double*, double*) ;

/** */
int setscr_ind  (int) ;

/** */
int settilim    (double, double*, double*) ;
/** */
void   gettilim    (double *) ;

/** */
int setcraind  ( int, int*, int* ) ;
/** */
int setadvind  ( int ) ;

/** */
int getobjsen   (CPXLPptr) ;

/** */
int setitlim    (int, int *, int *) ;

/** */
void getitlim    (int *) ;

/** */
int setobjllim  (double, double *, double *) ;

/** */
void getobjulim  (double *) ;

/** */
void getobjllim  (double *) ;

/** */
int getobjval   (CPXLPptr lp, double *z) ;

/** */
int getx        (CPXLPptr lp, double *x, int from, int upto) ;

/** */
int getpi       (CPXLPptr lp, double *pi, int from, int upto) ;

/** */
int getslack    (CPXLPptr lp, double *slack, int from, int upto) ;

/** */
int getdj       (CPXLPptr lp, double *dj, int from, int upto) ;

/** */
int getmethod   (CPXLPptr lp) ;

/** */
void setitfoind(int, int*, int*) ;

/** */
void getitfoind(int*) ;

/** */
void getdpriind(int*) ;

/** */
void getadvind(int*) ;

/** */
int getitc (CPXLPptr lp);

/** */
int getcols (CPXLPptr lp, int *nzcnt, int *cmatbeg, int *cmatind,
      double *cmatval, int cmatsz, int *surplus, int begin,
      int end);

#endif // _ITFCPLEX_H_
