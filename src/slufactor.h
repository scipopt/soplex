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
#pragma ident "@(#) $Id: slufactor.h,v 1.3 2001/11/07 17:31:20 bzfbleya Exp $"

#ifndef _SLUFACTOR_H_
#define _SLUFACTOR_H_


//@ ----------------------------------------------------------------------------
/*      \Section{Imports}
    Import required system include files
 */
#include <assert.h>


/*  and class header files
 */

#include "dvector.h"
#include "slinsolver.h"

//extern "C"
//{
#include "clutypes.h"
#include "clumembers.h"
//}


namespace soplex
{

#define MAXUPDATES      1000            // maximum nr. of factorization updates
// allowed before refactorization .


//@ ----------------------------------------------------------------------------
/*      \Section{Class Declaration}
 */
/** Sparse LU factorization.
 *  This is an implementation class for \Ref{SLinSolver} using sparse LU
 *  factorization.
 */
class SLUFactor : public SLinSolver, private CLUFactor
{
public:
   /// how to perform #change# method.
   enum UpdateType
   {       ///
      ETA = 0,
      ///
      FOREST_TOMLIN
   };

protected:
   void assign(const SLUFactor& old);
   void freeAll();
   void changeEta(int idx, SSVector& eta);

   DVector vec;
   SSVector ssvec;

   int usetup;         // 1 if update vector has been setup
   UpdateType uptype;
   SSVector eta;
   SSVector forest;
   double lastThreshold;

public:
   ///
   typedef SLinSolver::Status Status;

   /**@name Control Parameters */
   //@{
   /// minimum threshold to use.
   double minThreshold;

   /// #|x| < epsililon# is considered to be 0.
   double epsilon;

   /// minimum stability to acchieve by setting threshold.
   double minStability;

   ///
   UpdateType utype()
   {
      return uptype;
   }

   /** Set UpdateType.
    *  The new #UpdateType# becomes valid only after the next call to
    *  method #load()#.
    */
   void setUtype(UpdateType tp)
   {
      uptype = tp;
   }
   //@}

   /**@name derived from \Ref{SLinSolver}
    *  See documentation of \Ref{SLinSolver} for a documentation of these
    *  methods.
    */
   //@{
   ///
   void clear();
   ///
   int dim() const
   {
      return thedim;
   }
   ///
   int memory() const
   {
      return nzCnt + l.start[l.firstUnused];
   }
   ///
   Status status() const
   {
      return Status(stat);
   }
   ///
   double stability() const;

   ///
   Status load(const SVector* vec[], int dim);

   ///
   void solve2right(Vector& x, Vector& b);
   ///
   void solve2right(Vector& x, SSVector& b);
   ///
   void solve2right(SSVector& x, Vector& b);
   ///
   void solve2right(SSVector& x, SSVector& b);

   ///
   void solveRight (Vector& x,
                     const Vector& b);
   ///
   void solveRight (Vector& x,
                     const SVector& b);
   ///
   void solveRight (SSVector& x,
                     const Vector& b);
   ///
   void solveRight (SSVector& x,
                     const SVector& b);

   ///
   void solveRight4update(SSVector& x,
                           const SVector& b);
   ///
   void solve2right4update(SSVector& x,
                            Vector& two,
                            const SVector& b,
                            SSVector& rhs);

   ///
   void solve2left(Vector& x, Vector& b);
   ///
   void solve2left(Vector& x, SSVector& b);
   ///
   void solve2left(SSVector& x, Vector& b);
   ///
   void solve2left(SSVector& x, SSVector& b);

   ///
   void solveLeft (Vector& x,
                    const Vector& b);
   ///
   void solveLeft (Vector& x,
                    const SVector& b);
   ///
   void solveLeft (SSVector& x,
                    const Vector& b);
   ///
   void solveLeft (SSVector& x,
                    const SVector& b);

   ///
   void solveLeft (SSVector& x,
                    Vector& two,
                    const SVector& b,
                    SSVector& rhs2);

   ///
   Status change(int idx, const SVector& subst, const SSVector* eta = 0);
   //@}

   /** A zero vector.
    *  Return a zero #Vector# of the factorizations dimension. This may
    *  {\em temporarily} be used by other the caller in order to save
    *  memory (management overhead), but {\em must be reset to 0} when a
    *  method of #SLUFactor# is called.
    */
   Vector& zeroVec() //const
   {
      return vec; //((SLUFactor*)this)->vec;
   }

   ///
   void dump() const;
   ///
   int isConsistent() const;
   ///
   SLUFactor& operator=(const SLUFactor& old);
   ///
   SLUFactor(const SLUFactor& old)
      : SLinSolver( old )
         , vec (old.vec)
         , ssvec (old.ssvec)
         , eta (old.eta)
         , forest(old.forest)
   {
      assign(old);
   }
   ///
   SLUFactor();
   ///
   virtual ~SLUFactor();
};

} // namespace soplex
#endif

//-----------------------------------------------------------------------------
//Emacs Local Variables:
//Emacs mode:c++
//Emacs c-basic-offset:3
//Emacs tab-width:8
//Emacs indent-tabs-mode:nil
//Emacs End:
//-----------------------------------------------------------------------------
