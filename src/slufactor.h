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
#pragma ident "@(#) $Id: slufactor.h,v 1.6 2001/12/01 18:21:16 bzfbleya Exp $"

/**@file  slufactor.h
 * @brief Sparse LU factorization.
 */
#ifndef _SLUFACTOR_H_
#define _SLUFACTOR_H_

#include <assert.h>

#include "dvector.h"
#include "slinsolver.h"

#include "clutypes.h"
#include "clumembers.h"

namespace soplex
{

/// maximum nr. of factorization updates allowed before refactorization.
#define MAXUPDATES      1000     


/**@brief   Sparse LU factorization.
   @ingroup Algo
   
   This is an implementation class for #SLinSolver using sparse LU
   factorization.
*/
class SLUFactor : public SLinSolver, private CLUFactor
{
public:
   /**@todo document the two change methods ETA and FOREST_TOMLIN */
   /// how to perform #change method.
   enum UpdateType
   {
      ETA = 0,       ///< 
      FOREST_TOMLIN  ///<
   };

protected:
   /**@todo document these protected methods and attributes */
   void assign(const SLUFactor& old);
   void freeAll();
   void changeEta(int idx, SSVector& eta);

   DVector    vec;           ///<
   SSVector   ssvec;         ///<

   int        usetup;        ///< TRUE iff update vector has been setup
   UpdateType uptype;        ///< the current #UpdateType.
   SSVector   eta;           ///< 
   SSVector   forest;        ///<
   double     lastThreshold; ///<

public:
   typedef SLinSolver::Status Status;

   /**@name Control Parameters */
   //@{
   /// minimum threshold to use.
   double minThreshold;

   /// |x| < epsililon is considered to be 0.
   double epsilon;

   /// minimum stability to achieve by setting threshold.
   double minStability;

   /// returns the current update type #uptype.
   UpdateType utype()
   {
      return uptype;
   }

   /// sets update type.
   /** The new #UpdateType becomes valid only after the next call to
       method #load().
   */
   void setUtype(UpdateType tp)
   {
      uptype = tp;
   }
   //@}

   /**@todo should we document reimplemented derived methods again? */
   /**@name derived from SLinSolver
      See documentation of #SLinSolver for a documentation of these
      methods.
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

   /**@name Miscellaneous */
   //@{
   /// returns a zero vector.
   /** Returns a zero #Vector of the factorizations dimension. This may
       \em temporarily be used by other the caller in order to save
       memory (management overhead), but \em must \em be \em reset \em to \em 0
       when a method of #SLUFactor is called.
   */
   Vector& zeroVec() //const
   {
      return vec; //((SLUFactor*)this)->vec;
   }

   /// prints the LU factorization to stdout.
   void dump() const;

   /// consistency check.
   bool isConsistent() const;
   //@}

   /**@name Constructors / Destructors */
   //@{
   /// default constructor.
   SLUFactor();

   /// copy constructor.
   SLUFactor(const SLUFactor& old)
      : SLinSolver( old )
      , vec (old.vec)
      , ssvec (old.ssvec)
      , eta (old.eta)
      , forest(old.forest)
   {
      assign(old);
   }

   /// assignment operator.
   SLUFactor& operator=(const SLUFactor& old);

   /// destructor.
   virtual ~SLUFactor();
   //@}
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
