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
#pragma ident "@(#) $Id: clutypes.h,v 1.2 2001/11/06 23:31:00 bzfkocht Exp $"

#ifndef _CLUTYPES_H_
#define _CLUTYPES_H_

namespace soplex
{

/*
 *      This file defines all C data structures involved in the high-speed
 *      C sparse LU factorization routine.
 */

/*      Data structures for saving the row and column permutations.
 */
typedef struct
{
   int *orig;          /* orig[p] original index from p */
   int *perm;          /* perm[i] permuted index from i */
} Perm;


/*      Double linked ring structure for garbage collection of column or
 *      row file in working matrix.
 */
typedef struct _dr_
{
   struct _dr_ *next;
   struct _dr_ *prev;
   int idx;
} Dring;

/*      Data structures for saving the working matrix and U factor.
 */
typedef struct
{
   struct Row
   {
      Dring list;           /* Double linked ringlist of vector
                                         indices in the order they appear
                                         in the row file */
      Dring *elem;          /* Array of ring elements.            */
      int size;           /* size of arrays val and idx         */
      int used;           /* used entries of arrays idx and val */
      double *val;           /* hold nonzero values                */
      int *idx;           /* hold nonzero indices               */
      int *start;         /* starting positions in val and idx  */
      int *len;           /* used nonzeros per row vectors      */
      int *max;           /* maximum available nonzeros per row:
                                         start[i] + max[i] == start[elem[i].next->idx] 
                                         len[i] <= max[i].
                                       */
   } row;
   struct Col
   {
      Dring list;           /* Double linked ringlist of vector
                                         indices in the order they appear
                                         in the column file */
      Dring *elem;          /* Array of ring elements.            */
      int size;           /* size of array idx                  */
      int used;           /* used entries of array idx          */
      int *idx;           /* hold nonzero indices               */
      double *val;           /* hold nonzero values: this is only initialized
                                         in the end of the factorization with DEFAULT
                                         updates.
                                       */
      int *start;         /* starting positions in val and idx  */
      int *len;           /* used nonzeros per column vector    */
      int *max;           /* maximum available nonzeros per colunn:
                                         start[i] + max[i] == start[elem[i].next->idx] 
                                         len[i] <= max[i].
                                       */
   }
   col;
   int lastColSing;            /* stage of last eliminated column singleton */
   int lastRowSing;            /* stage of last eliminated row singleton */
} U;


/*      Data structures for saving the working matrix and U factor.
 */
typedef struct
{
   int size;           /* size of arrays val and idx        */
   double *val;           /* values of L vectors               */
   int *idx;           /* indices of L vectors              */
   int startSize;      /* size of array start               */
   int firstUpdate;    /* number of first update L vector   */
   int firstUnused;    /* number of first unused L vector   */
   int *start;         /* starting positions in val and idx */
   int *row;           /* column indices of L vectors       */
   int updateType;     /* type of updates to be used.       */

   /*  The following arrays have length |firstUpdate|, since they keep
       rows of the L-vectors occuring during the factorization (without
       updates), only:
    */
   double *rval;          /* values of rows of L               */
   int *ridx;          /* indices of rows of L              */
   int *rbeg;          /* start of rows in rval and ridx    */
   int *rorig;         /* original row permutation          */
   int *rperm;         /* original row permutation          */
} L;


} // namespace soplex
#endif // _CLUTYPES_H_

//-----------------------------------------------------------------------------
//Emacs Local Variables:
//Emacs mode:c++
//Emacs c-basic-offset:3
//Emacs tab-width:8
//Emacs indent-tabs-mode:nil
//Emacs End:
//-----------------------------------------------------------------------------
