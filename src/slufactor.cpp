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
#pragma ident "@(#) $Id: slufactor.cpp,v 1.1 2001/11/06 16:18:31 bzfkocht Exp $"



#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stddef.h>
#include <assert.h>

#include "slufactor.h"


#include "cluprotos.h"
#include "cring.h"

namespace soplex
{
extern double verySparseFactor4right;
extern double verySparseFactor4left;


#define MINSTABILITY    1e-2

//@ ----------------------------------------------------------------------------

void SLUFactor::solve2right(Vector& x, Vector& b) const
{
soplex::solveRight((CLUFactor*)this, x.get_ptr(), b.get_ptr());
}

void SLUFactor::solve2right(Vector& x, SSVector& b) const
{
   vSolveRightNoNZ((CLUFactor*)this, x.get_ptr(), b.epsilon,
                    b.altValues(), b.altIndexMem(), b.size());
}

void SLUFactor::solve2right(SSVector& x, Vector& b) const
{
   x.clear();
soplex::solveRight((CLUFactor*)this, x.altValues(), b.get_ptr());
}

void SLUFactor::solve2right(SSVector& x, SSVector& b) const
{
   int n;
   int bs = b.size();
   x.clear();

   n = vSolveRight4update((CLUFactor*)this, x.epsilon, x.altValues(), x.altIndexMem(),
                           b.altValues(), b.altIndexMem(), bs, 0, 0, 0);
   if (n > 0)
   {
      x.setSize(n);
      x.forceSetup();
   }
   else
      x.unSetup();

   b.setSize(0);
   b.forceSetup();
}

void SLUFactor::solveRight (Vector& x, const Vector& b) const
{
   ((SLUFactor*)this)->vec = b;
   solve2right(x, ((SLUFactor*)this)->vec);
}

void SLUFactor::solveRight (Vector& x,
                             const SVector& b) const
{
   ((SLUFactor*)this)->vec.assign(b);
   solve2right(x, ((SLUFactor*)this)->vec);
}

void SLUFactor::solveRight (SSVector& x,
                             const Vector& b) const
{
   ((SLUFactor*)this)->vec = b;
   solve2right(x, ((SLUFactor*)this)->vec);
}

void SLUFactor::solveRight (SSVector& x,
                             const SVector& b) const
{
   ((SLUFactor*)this)->vec.assign(b);
   solve2right(x, ((SLUFactor*)this)->vec);
}

void SLUFactor::solveRight4update(SSVector& x,
                                   const SVector& b)
{
   int m, n, f;

   x.clear();
   ssvec = b;
   n = b.size();
   if (l.updateType == ETA)
   {
      m = vSolveRight4update(this, x.epsilon,
                              x.altValues(), x.altIndexMem(),
                              ssvec.altValues(), ssvec.altIndexMem(), n,
                              0, 0, 0);
      x.setSize(m);
      x.forceSetup();
      eta.setup_and_assign(x);
   }
   else
   {
      forest.clear();
      m = vSolveRight4update(this, x.epsilon,
                              x.altValues(), x.altIndexMem(),
                              ssvec.altValues(), ssvec.altIndexMem(), n,
                              forest.altValues(), &f, forest.altIndexMem());
      forest.setSize(f);
      forest.forceSetup();
      x.setSize(m);
      x.forceSetup();
   }
   ((SLUFactor*)this)->usetup = 1;
}

void SLUFactor::solve2right4update(SSVector& x,
                                    Vector& y,
                                    const SVector& b,
                                    SSVector& rhs)
{
   int m, n, f;

   Vector& svec = *(Vector*) & ((SLUFactor*)this)->ssvec;
   int* sidx = ssvec.altIndexMem();
   int rsize = rhs.size();
   int* ridx = rhs.altIndexMem();

   x.clear();
   y.clear();
   ((SLUFactor*)this)->usetup = 1;
   ((SLUFactor*)this)->ssvec = b;
   if (l.updateType == ETA)
   {
      n = b.size();
      m = vSolveRight4update2(this, x.epsilon,
                               x.altValues(), x.altIndexMem(), svec.get_ptr(),
                               sidx, n, y.get_ptr(),
                               rhs.epsilon, rhs.altValues(), ridx, rsize,
                               0, 0, 0);
      x.setSize(m);
      x.forceSetup();
      eta.setup_and_assign(x);
      //      ((SLUFactor*)this)->eta.setup_and_assign(x);
   }
   else
   {
      ((SLUFactor*)this)->forest.clear();
      n = ssvec.size();
      m = vSolveRight4update2(this, x.epsilon,
                               x.altValues(), x.altIndexMem(), svec.get_ptr(),
                               sidx, n, y.get_ptr(),
                               rhs.epsilon, rhs.altValues(), ridx, rsize,
                               forest.altValues(), &f, forest.altIndexMem());
      x.setSize(m);
      x.forceSetup();
      forest.setSize(f);
      ((SLUFactor*)this)->forest.forceSetup();
   }
}

//@ ----------------------------------------------------------------------------

void SLUFactor::solve2left (Vector& x, Vector& b) const
{
   x.clear();
soplex::solveLeft(x.get_ptr(), (CLUFactor*)this, b.get_ptr());
}

void SLUFactor::solve2left(Vector& x, SSVector& b) const
{
   x.clear();
   vSolveLeftNoNZ((CLUFactor*)this, b.epsilon, x.get_ptr(),
                   b.altValues(), b.altIndexMem(), b.size());
}

void SLUFactor::solve2left(SSVector& x, Vector& b) const
{
   int n;
   x.clear();
   n = soplex::solveLeftEps (x.altValues(), (CLUFactor*)this,
                        b.get_ptr(), x.altIndexMem(), x.epsilon);
   if (n)
   {
      x.setSize(n);
      x.forceSetup();
   }
}

void SLUFactor::solve2left(SSVector& x, SSVector& b) const
{
   int n;
   int bs = b.size();
   x.clear();

   n = vSolveLeft((CLUFactor*)this, x.epsilon, x.altValues(), x.altIndexMem(),
                   b.altValues(), b.altIndexMem(), bs);

   if (n > 0)
   {
      x.setSize(n);
      x.forceSetup();
   }
   else
      x.unSetup();
   b.setSize(0);
   b.forceSetup();
}

void SLUFactor::solveLeft(Vector& x,
                           const SVector& b) const
{
   ((SLUFactor*)this)->ssvec = b;
   solve2left(x, ((SLUFactor*)this)->ssvec);
}


void SLUFactor::solveLeft (Vector& x,
                            const Vector& b) const
{
   ((SLUFactor*)this)->vec = b;
   solve2left(x, ((SLUFactor*)this)->vec);
}

void SLUFactor::solveLeft (SSVector& x,
                            const Vector& b) const
{
   ((SLUFactor*)this)->vec = b;
   solve2left(x, ((SLUFactor*)this)->vec);
}

void SLUFactor::solveLeft (SSVector& x,
                            const SVector& b) const
{
   ((SLUFactor*)this)->ssvec.assign(b);
   SLUFactor::solve2left(x, ((SLUFactor*)this)->ssvec);
}

void SLUFactor::solveLeft (SSVector& x,
                            Vector& y,
                            const SVector& rhs1,
                            SSVector& rhs2) const
{
   int n;
   double* svec = ((SLUFactor*)this)->ssvec.altValues();
   int* sidx = ((SLUFactor*)this)->ssvec.altIndexMem();
   int rn = rhs2.size();
   int* ridx = rhs2.altIndexMem();

   x.clear();
   y.clear();
   ((SLUFactor*)this)->ssvec.assign(rhs1);
   n = ssvec.size();
   n = vSolveLeft2((CLUFactor*)this, x.epsilon,
                    x.altValues(), x.altIndexMem(), svec, sidx, n,
                    y.get_ptr(), rhs2.altValues(), ridx, rn);
   x.setSize(n);
   if (n > 0)
      x.forceSetup();
   else
      x.unSetup();
   rhs2.setSize(0);
   rhs2.forceSetup();
   ((SLUFactor*)this)->ssvec.setSize(0);
   ((SLUFactor*)this)->ssvec.forceSetup();
}


//@ ----------------------------------------------------------------------------
double SLUFactor::stability() const
{
   if (status() != OK)
      return 0;
   if (maxabs < initMaxabs)
      return 1;
   return initMaxabs / maxabs;
}

void SLUFactor::changeEta(int idx, SSVector& et)
{
   int es = et.size();
   stat = updateCLUFactor(this, idx, et.altValues(), et.altIndexMem(), es);
   et.setSize(0);
   et.forceSetup();
}

SLUFactor::Status SLUFactor::change(
   int idx,
   const SVector& subst,
   const SSVector* e
)
{
   if (usetup)
   {
      if (l.updateType)                      /// Forest-Tomlin updates
      {
         int fsize = forest.size();
         stat = forestUpdateCLUFactor(this, idx,
                                       forest.altValues(), fsize, forest.altIndexMem());
         forest.setSize(0);
         forest.forceSetup();
      }
      else                                    /// ETA updates
         changeEta(idx, eta);
   }
   else if (e)                                /// ETA updates
   {
      l.updateType = ETA;
      stat = updateCLUFactorNoClear(this, idx, e->values(), e->indexMem(), e->size());
   }
   else if (l.updateType)                     /// Forest-Tomlin updates
   {
      forest = subst;
      solveLright(this, forest.altValues());
      stat = forestUpdateCLUFactor(this, idx, forest.altValues(), 0, 0);
      forest.setSize(0);
      forest.forceSetup();
   }
   else                                        /// ETA updates
   {
      vec = subst;
      solve2right(eta, vec);
      changeEta(idx, eta);
   }
   usetup = 0;
#ifdef  DEBUG
   std::cerr << "\tupdated\t\tstability = " << stability() << std::endl;
#endif 
   return status();
}

//@ ----------------------------------------------------------------------------

void SLUFactor::clear()
{
   rowMemMult = 5;             /* factor of minimum Memory * #of nonzeros */
   colMemMult = 5;             /* factor of minimum Memory * #of nonzeros */
   lMemMult = 1;             /* factor of minimum Memory * #of nonzeros */

   l.firstUpdate = 0;
   l.firstUnused = 0;
   thedim = 1;

   epsilon = 1e-24;
   usetup = 0;
   maxabs = 1;
   initMaxabs = 1;
   minThreshold = 0.01;
   lastThreshold = minThreshold;
   minStability = MINSTABILITY;
   stat = UNLOADED;

   vec.clear();
   eta.clear();
   ssvec.clear();
   forest.clear();

   u.row.size = 100;
   u.col.size = 100;
   l.size = 100;
   l.startSize = 100;

   if (l.val)
   {
      Free(u.row.val);
      Free(u.row.idx);
      Free(u.col.idx);
      Free(l.val);
      Free(l.idx);
      Free(l.start);
      Free(l.row);
   }

   u.row.val = (double*)Malloc(u.row.size * sizeof(double));
   u.row.idx = (int *)Malloc(u.row.size * sizeof(int));
   u.col.idx = (int *)Malloc(u.col.size * sizeof(int));
   l.val = (double*)Malloc(l.size * sizeof(double));
   l.idx = (int *)Malloc(l.size * sizeof(int));
   l.start = (int *)Malloc(l.startSize * sizeof(int));
   l.row = (int *)Malloc(l.startSize * sizeof(int));
}

void SLUFactor::assign(const SLUFactor& old)
{
   thedim = old.thedim;
   rowMemMult = old.rowMemMult;
   colMemMult = old.colMemMult;
   lMemMult = old.lMemMult;
   stat = old.stat;
   epsilon = old.epsilon;
   minThreshold = old.minThreshold;
   minStability = old.minStability;
   maxabs = old.maxabs;
   initMaxabs = old.initMaxabs;

   row.perm = (int*)Malloc(thedim * sizeof(int));
   row.orig = (int*)Malloc(thedim * sizeof(int));
   col.perm = (int*)Malloc(thedim * sizeof(int));
   col.orig = (int*)Malloc(thedim * sizeof(int));
   diag = (double*)Malloc(thedim * sizeof(double));

   memcpy(row.perm, old.row.perm, thedim * sizeof(int));
   memcpy(row.orig, old.row.orig, thedim * sizeof(int));
   memcpy(col.perm, old.col.perm, thedim * sizeof(int));
   memcpy(col.orig, old.col.orig, thedim * sizeof(int));
   memcpy(diag, old.diag, thedim * sizeof(double));

   work = vec.get_ptr();

   u.row.size = old.u.row.size;
   u.row.used = old.u.row.used;
   u.row.elem = (Dring *)Malloc(thedim * sizeof(Dring));
   u.row.val = (double*)Malloc(u.row.size * sizeof(double));
   u.row.idx = (int *)Malloc(u.row.size * sizeof(int));
   u.row.start = (int *)Malloc((thedim + 1) * sizeof(int));
   u.row.len = (int *)Malloc((thedim + 1) * sizeof(int));
   u.row.max = (int *)Malloc((thedim + 1) * sizeof(int));

   memcpy(u.row.elem , old.u.row.elem, thedim * sizeof(Dring));
   memcpy(u.row.val , old.u.row.val, u.row.size * sizeof(double));
   memcpy(u.row.idx , old.u.row.idx, u.row.size * sizeof(int));
   memcpy(u.row.start , old.u.row.start, (thedim + 1) * sizeof(int));
   memcpy(u.row.len , old.u.row.len, (thedim + 1) * sizeof(int));
   memcpy(u.row.max , old.u.row.max, (thedim + 1) * sizeof(int));

   if (thedim && stat == OK)
   {           // make row list ok
      u.row.list = old.u.row.list;
      const Dring* oring = &old.u.row.list;
      Dring* ring = &u.row.list;
      for (; oring->next != &old.u.row.list; oring = oring->next, ring = ring->next)
      {
         ring->next = &(u.row.elem[oring->next->idx]);
         ring->next->prev = ring;
      }
      ring->next = &u.row.list;
      ring->next->prev = ring;
   }

   u.col.size = old.u.col.size;
   u.col.used = old.u.col.used;
   u.col.elem = (Dring *)Malloc(thedim * sizeof(Dring));
   u.col.idx = (int *)Malloc(u.col.size * sizeof(int));
   u.col.start = (int *)Malloc((thedim + 1) * sizeof(int));
   u.col.len = (int *)Malloc((thedim + 1) * sizeof(int));
   u.col.max = (int *)Malloc((thedim + 1) * sizeof(int));
   if (old.u.col.val)
   {
      u.col.val = (double*)Malloc(u.col.size * sizeof(double));
      memcpy(u.col.val, old.u.col.val, u.col.size * sizeof(double));
   }
   else
      u.col.val = 0;

   memcpy(u.col.elem , old.u.col.elem , thedim * sizeof(Dring));
   memcpy(u.col.idx , old.u.col.idx , u.col.size * sizeof(int));
   memcpy(u.col.start , old.u.col.start , (thedim + 1) * sizeof(int));
   memcpy(u.col.len , old.u.col.len , (thedim + 1) * sizeof(int));
   memcpy(u.col.max , old.u.col.max , (thedim + 1) * sizeof(int));

   if (thedim && stat == OK)
   {           // make col list ok
      u.col.list = old.u.col.list;
      const Dring* oring = &old.u.col.list;
      Dring* ring = &u.col.list;
      for (; oring->next != &old.u.col.list; oring = oring->next, ring = ring->next)
      {
         ring->next = &(u.col.elem[oring->next->idx]);
         ring->next->prev = ring;
      }
      ring->next = &u.col.list;
      ring->next->prev = ring;
   }

   nzCnt = old.nzCnt;
   uptype = old.uptype;
   l.size = old.l.size;
   l.startSize = old.l.startSize;
   l.firstUpdate = old.l.firstUpdate;
   l.firstUnused = old.l.firstUnused;
   l.updateType = old.l.updateType;
   l.val = (double*)Malloc(l.size * sizeof(double));
   l.idx = (int *)Malloc(l.size * sizeof(int));
   l.start = (int *)Malloc(l.startSize * sizeof(int));
   l.row = (int *)Malloc(l.startSize * sizeof(int));
   if (old.l.rbeg)
   {
      l.rval = (double*)Malloc(l.firstUpdate * sizeof(double));
      l.ridx = (int *)Malloc(l.firstUpdate * sizeof(int));
      l.rbeg = (int *)Malloc((thedim + 1) * sizeof(int));
      memcpy(l.rval, old.l.rval, l.firstUpdate * sizeof(double));
      memcpy(l.ridx, old.l.ridx, l.firstUpdate * sizeof(int));
      memcpy(l.rbeg, old.l.rbeg, (thedim + 1) * sizeof(int));
   }
   else
   {
      l.rval = 0;
      l.ridx = 0;
      l.rbeg = 0;
   }

   memcpy(l.val , old.l.val , l.size * sizeof(double));
   memcpy(l.idx , old.l.idx , l.size * sizeof(int));
   memcpy(l.start , old.l.start , l.startSize * sizeof(int));
   memcpy(l.row , old.l.row , l.startSize * sizeof(int));

   assert
   (
      row.perm != 0
      || row.orig != 0
      || col.perm != 0
      || col.orig != 0
      || diag != 0
      || u.row.elem != 0
      || u.row.val != 0
      || u.row.idx != 0
      || u.row.start != 0
      || u.row.len != 0
      || u.row.max != 0
      || u.col.elem != 0
      || u.col.idx != 0
      || u.col.start != 0
      || u.col.len != 0
      || u.col.max != 0
      || l.val != 0
      || l.idx != 0
      || l.start != 0
      || l.row != 0
  );
}

SLUFactor& SLUFactor::operator=(const SLUFactor& old)
{
   freeAll();
   vec = old.vec;
   ssvec = old.ssvec;
   assign(old);
   return *this;
}

SLUFactor::SLUFactor()
   : vec (1)
      , ssvec (1, 1e-16)
      , usetup (0)
      , uptype (FOREST_TOMLIN)
      , eta (1, 1e-16)
      , forest (1, 1e-16)
{
   nzCnt = 0;
   thedim = 1;

   row.perm = (int*)Malloc(thedim * sizeof(int));
   row.orig = (int*)Malloc(thedim * sizeof(int));
   col.perm = (int*)Malloc(thedim * sizeof(int));
   col.orig = (int*)Malloc(thedim * sizeof(int));

   diag = (double*)Malloc(thedim * sizeof(double));
   assert(diag && "ERROR: out of memory");
   work = vec.get_ptr();

   u.row.size = 1;
   u.row.used = 0;
   u.row.elem = (Dring *)Malloc(thedim * sizeof(Dring));
   u.row.val = (double*)Malloc(u.row.size * sizeof(double));
   u.row.idx = (int *)Malloc(u.row.size * sizeof(int));
   u.row.start = (int *)Malloc((thedim + 1) * sizeof(int));
   u.row.len = (int *)Malloc((thedim + 1) * sizeof(int));
   u.row.max = (int *)Malloc((thedim + 1) * sizeof(int));

   u.row.list.idx = thedim;
   u.row.start[thedim] = 0;
   u.row.max[thedim] = 0;
   u.row.len[thedim] = 0;

   u.col.size = 1;
   u.col.used = 0;
   u.col.elem = (Dring *)Malloc(thedim * sizeof(Dring));
   u.col.idx = (int *)Malloc(u.col.size * sizeof(int));
   u.col.start = (int *)Malloc((thedim + 1) * sizeof(int));
   u.col.len = (int *)Malloc((thedim + 1) * sizeof(int));
   u.col.max = (int *)Malloc((thedim + 1) * sizeof(int));
   u.col.val = 0;

   u.col.list.idx = thedim;
   u.col.start[thedim] = 0;
   u.col.max[thedim] = 0;
   u.col.len[thedim] = 0;

   l.size = 1;
   l.val = (double*)Malloc(l.size * sizeof(double));
   l.idx = (int *)Malloc(l.size * sizeof(int));
   l.startSize = 1;
   l.firstUpdate = 0;
   l.firstUnused = 0;
   l.start = (int *)Malloc(l.startSize * sizeof(int));
   l.row = (int *)Malloc(l.startSize * sizeof(int));
   l.rval = 0;
   l.ridx = 0;
   l.rbeg = 0;
   clear();

   assert
   (
      row.perm != 0
      || row.orig != 0
      || col.perm != 0
      || col.orig != 0
      || diag != 0
      || u.row.elem != 0
      || u.row.val != 0
      || u.row.idx != 0
      || u.row.start != 0
      || u.row.len != 0
      || u.row.max != 0
      || u.col.elem != 0
      || u.col.idx != 0
      || u.col.start != 0
      || u.col.len != 0
      || u.col.max != 0
      || l.val != 0
      || l.idx != 0
      || l.start != 0
      || l.row != 0
  );
}

void SLUFactor::freeAll()
{
   Free (row.perm);
   Free (row.orig);
   Free (col.perm);
   Free (col.orig);
   Free (u.row.elem);
   Free (u.row.val);
   Free (u.row.idx);
   Free (u.row.start);
   Free (u.row.len);
   Free (u.row.max);
   Free (u.col.elem);
   Free (u.col.idx);
   Free (u.col.start);
   Free (u.col.len);
   Free (u.col.max);
   Free (l.val);
   Free (l.idx);
   Free (l.start);
   Free (l.row);
   Free (diag);
   row.perm = 0;
   row.orig = 0;
   col.perm = 0;
   col.orig = 0;
   u.row.elem = 0;
   u.row.val = 0;
   u.row.idx = 0;
   u.row.start = 0;
   u.row.len = 0;
   u.row.max = 0;
   u.col.elem = 0;
   u.col.idx = 0;
   u.col.start = 0;
   u.col.len = 0;
   u.col.max = 0;
   l.val = 0;
   l.idx = 0;
   l.row = 0;
   l.start = 0;
   diag = 0;
   if (u.col.val)
   {
      free (u.col.val);
      u.col.val = 0;
   }
   if (l.rbeg)
   {
      Free(l.rval);
      Free(l.ridx);
      Free(l.rbeg);
      Free(l.rorig);
      Free(l.rperm);
      l.rval = 0;
      l.ridx = 0;
      l.rbeg = 0;
      l.rorig = 0;
      l.rperm = 0;
   }
}

SLUFactor::~SLUFactor()
{
   freeAll();
}

static double betterThreshold(double th)
{
   if (th < 0.1)
      th *= 10;
   else if (th < 0.9)
      th = (th + 1) / 2;
   else if (th < 0.999)
      th = 0.99999;
   return th;
}

SLUFactor::Status SLUFactor::load(const SVector* matrix[], int dm)
{
   assert(dm > 0);
   assert(matrix);
   double lastStability = stability();

   initDR(u.row.list);
   initDR(u.col.list);

   usetup = 0;
   l.updateType = uptype;
   l.firstUpdate = 0;
   l.firstUnused = 0;

   if (dm != thedim)
   {
      clear();

      thedim = dm;
      vec.reDim(thedim);
      ssvec.reDim(thedim);
      eta.reDim(thedim);
      forest.reDim(thedim);
      work = vec.get_ptr();

      row.perm = (int*)Realloc(row.perm, thedim * sizeof(int));
      row.orig = (int*)Realloc(row.orig, thedim * sizeof(int));
      col.perm = (int*)Realloc(col.perm, thedim * sizeof(int));
      col.orig = (int*)Realloc(col.orig, thedim * sizeof(int));
      diag = (double*)Realloc(diag, thedim * sizeof(double));

      u.row.elem = (Dring*)Realloc(u.row.elem, thedim * sizeof(Dring));
      u.row.len = (int *)Realloc(u.row.len, (thedim + 1) * sizeof(int));
      u.row.max = (int *)Realloc(u.row.max, (thedim + 1) * sizeof(int));
      u.row.start = (int *)Realloc(u.row.start, (thedim + 1) * sizeof(int));

      u.col.elem = (Dring*)Realloc(u.col.elem, thedim * sizeof(Dring));
      u.col.len = (int *)Realloc(u.col.len, (thedim + 1) * sizeof(int));
      u.col.max = (int *)Realloc(u.col.max, (thedim + 1) * sizeof(int));
      u.col.start = (int *)Realloc(u.col.start, (thedim + 1) * sizeof(int));

      l.startSize = thedim + MAXUPDATES;
      l.row = (int *)Realloc(l.row, (l.startSize) * sizeof(int));
      l.start = (int *)Realloc(l.start, (l.startSize) * sizeof(int));

      assert
      (
         row.perm != 0
         || row.orig != 0
         || col.perm != 0
         || col.orig != 0
         || diag != 0
         || u.row.elem != 0
         || u.row.start != 0
         || u.row.len != 0
         || u.row.max != 0
         || u.col.elem != 0
         || u.col.start != 0
         || u.col.len != 0
         || u.col.max != 0
         || l.row != 0
         || l.start != 0
     );
   }
   else if (lastStability > 2*minStability)
   {
      double last = minThreshold;
      double better = betterThreshold(last);
      while (better < lastThreshold)
      {
         last = better;
         better = betterThreshold(last);
      }
      minStability = 2 * MINSTABILITY;
      lastThreshold = last;
   }

   u.row.list.idx = thedim;
   u.row.start[thedim] = 0;
   u.row.max[thedim] = 0;
   u.row.len[thedim] = 0;

   u.col.list.idx = thedim;
   u.col.start[thedim] = 0;
   u.col.max[thedim] = 0;
   u.col.len[thedim] = 0;

   for (;;)
   {
      stat = OK;
      factor(this, (SVector**)matrix, lastThreshold, epsilon);
      if (stability() >= minStability)
         break;
      double x = lastThreshold;
      lastThreshold = betterThreshold(lastThreshold);
      if (x == lastThreshold)
         break;
      minStability /= 2;
   }
   /*
   std::cerr << "threshold = " << lastThreshold
        << "\tstability = " << stability()
        << "\tminStability = " << minStability << std::endl;
    */

#ifndef NDEBUG
   int i = 0;
   if (i)
   {
      FILE* fl = fopen("dump.lp", "w");
      std::cout << "Basis:\n";
      int j = 0;
      for (i = 0; i < dim(); ++i)
         j += matrix[i]->size();
      for (i = 0; i < dim(); ++i)
      {
         for (j = 0; j < matrix[i]->size(); ++j)
            fprintf(fl, "%8d  %8d  %16g\n",
                     i + 1, matrix[i]->index(j) + 1, matrix[i]->value(j));
      }
      fclose(fl);
      std::cout << "LU-Factors:\n";
      dump();
   }
   std::cerr << "threshold = " << lastThreshold << "\tstability = " << stability() << std::endl;
#endif

   assert(isConsistent());
   return Status(stat);
}


int SLUFactor::isConsistent() const
{
   return CLUFactorIsConsistent((CLUFactor*)this);
}

void SLUFactor::dump() const
{
soplex::dumpCLUFactor((CLUFactor*)this);
}
} // namespace soplex

//-----------------------------------------------------------------------------
//Emacs Local Variables:
//Emacs c-basic-offset:3
//Emacs tab-width:8
//Emacs indent-tabs-mode:nil
//Emacs End:
//-----------------------------------------------------------------------------
