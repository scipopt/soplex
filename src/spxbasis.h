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
#pragma ident "@(#) $Id: spxbasis.h,v 1.3 2001/11/07 17:31:21 bzfbleya Exp $"



#ifndef _SPXBASIS_H_
#define _SPXBASIS_H_


/*      \Section{Imports}
 */
#include <assert.h>


#include "spxlp.h"
#include "svector.h"
#include "ssvector.h"
#include "dataarray.h"
#include "slinsolver.h"
#include "nameset.h"

namespace soplex
{





class SoPlex;

//@ -----------------------------------------------------------------------------
/*      \Section{Class SPxBasis}
 */

/** simplex basis.
    Consider the linear program as provided from class #SPxLP#
    \[
        \begin{array}{rl}
            \hbox{max}  & c^T x         \\
            \hbox{s.t.} & l_r \le Ax \le u_r    \\
                            & l_c \le x \le u_c
        \end{array}
    \]
    where $c, l_c, u_c, x \in {\bf R}^n$, $l_r, u_r \in {\bf R}^m$ and $A \in
    {\bf R}^{m \times n}$.  Solving this LP with the simplex algorithm requires
    the definition of a {\em basis}. Such can be defined as a set of column vectors
    or a set of row vectors building a nonsingular matrix. We will refer to the
    first case as the {\em columnwise representation} and the latter case will
    be called the {\em rowwise representation}. In both cases, a {\em basis} is
    a set of vectors forming a nonsigular matrix.  The dimension of the vectors
    is refferred to as the {\em basis' dimension}, whereas the number of vectors
    belonging to the LP is called the basis' {\em codimension}.
 
    Class #SPxBasis# is designed to represent a generic simplex basis, suitable
    for both representations. At any time the representation can be changed by
    calling method #setRep()#.
 
    Class #SPxBasis# provides methods for solving linear systems with the basis
    matrix. However, #SPxBasis# does not provide a linear solver by its own.
    Instead, a #SLinSolver# object must be #load#ed to a #SPxBasis# which will
    be called for solving linear systems.
*/
/*      \Section{Class Declaration}
 */
class SPxBasis
{
protected:
   SoPlex* theLP;
   
public:

   /**@name Basis status */
   //@{
   /** Basis status.
       Each #SPxBasis# is assigned a status flag, which can take on of the
       above values.
    */
   enum SPxStatus
   {
      /// No Problem has been loaded to the basis
      NO_PROBLEM = -2 ,
      /// Basis is singular
      SINGULAR = -1 ,
      /// Basis is not know to be dual nor primal feasible
      REGULAR = 0 ,
      /// Basis is dual feasible
      DUAL = 1 ,
      /// Basis is primal feasible
      PRIMAL = 2 ,
      /// Basis is optimal, i.e. dual and primal feasible
      OPTIMAL = 3 ,
      /// LP has been proven to be unbounded
      UNBOUNDED = 4 ,
      /// LP has been proven to be infeasible
      INFEASIBLE = 5
   };

private:
   SPxStatus thestatus;              // current status of the basis

public:
   /// current #SPxStatus#.
   SPxStatus status() const
   {
      return thestatus;
   }

   /// set basis #SPxStatus# to #stat#.
   void setStatus(SPxStatus stat)
   {
      thestatus = stat;
   }
   //@}



   /**@name Basis descriptor */
   //@{
   ///
   class Desc
   {
public:
      /** Status of a variable.
          A basis is described by assigning a #Status# to all of the LPs
          variables and covariables. This assignment is maintained by the
          basis #Desc#riptor.

          Variables and covariables may have a primal or dual #Status#. The
          first type specifies that a variable is set on a primal bound, while
          the later type indicates a dual variable to be set on a bound.
          If a row variable has a primal status, #P_ON_UPPER# say, this means
          that the upper bound of the inequality is set to be tight. Hence,
          in this case the upper bound must not be infinity.

          Equivalently, if the status of a variable is dual, #D_ON_UPPER# say,
          it means that the dual variable corresponding to the upper bound
          inequality of this variable is set to 0.

          For a column basis, primal #Status#'s correspond to nonbasic
          variables, while dual ones are basic. This is reversed for a row
          basis. We will now reveil in more detail the significance of
          variable #Status#'s.

          {\bf Primal Variables}\\
          Let consider a range inequality $l_r \le a^T x \le u_r$ or bounds on
          a variable $l_c \le x_c \le u_c$ The following table reveils what is
          implied if the corresponding variable or covariable is assigned to a
          primal #Status#:

      \begin{center}
      \begin{tabular}{lcl}
      $l_c \le x_c \le u_c$   & {\bf #Status#($x_i$)} & $l_r \le a^T x \le u_r$ \\
      \hline
      $x_c = u_c < \infty$    & #P_ON_UPPER#  & $a^T x = u_r < \infty$        \\
      $x_c = l_c > -\infty$   & #P_ON_LOWER#  & $a^T x = l_r > -\infty$\\
      $-\infty < l_c = x_c = u_c < \infty$
                              & #P_FIXED#     &
                                      $-\infty < l_r = a^T x = u_r < \infty$  \\
      $-\infty = l_i < x_i=0 < u_i = \infty$
                              & #P_FREE#      &
                                  $-\infty = l_r < a^T x = 0 < u_r = \infty$  \\
      \end{tabular}
      \end{center}

          Note that to determine if a variable with #Status stat# is set to
          its upper bound, on can compute the following test:
          #(-stat | -P_ON_UPPER)#. This will yield true even if the variable
          is fixed, i.e. sitting on both bounds at the same time.

          {\bf Dual Variables}\\
          In principle for implementing the Simplex algorithm it would suffice
          to use only one dual #Status#. However, for performance reasons it
          is advisable to introduce various dual #Status# types, reflecting
          the structure of the bounds. Given an upper bound $u$ and a lower
          bound $l$ of a constraint or variable, the following table indicates
          the setting of the dual #Status# of this variable.

          \begin{center}
          \begin{tabular}{cl}
              $l \le ... \le u$               & #Status#              \\
          \hline
              $-\infty < l \ne u < \infty$    & #D_ON_BOTH#   \\
              $-\infty < l \ne u = \infty$    & #D_ON_UPPER#  \\
              $-\infty = l \ne u < \infty$    & #D_ON_LOWER#  \\
              $-\infty < l  =  u < \infty$    & #D_FREE#              \\
              $-\infty = l \ne u = \infty$    & #D_UNDEFINED# \\
          \end{tabular}
          \end{center}

          Note that unbounded primal variables are reflected by an #D_UNDEFINED#
          dual variable, since no DUAL variables exist to them. To facilate
          the assignment of dual #Status#'s, class #SPxBasis# provides methods
          \Ref{dualStatus}, \Ref{dualColStatus} and \Ref{dualRowStatus}.
      */
      enum Status
      {
         /// primal variable is set to its lower bound
         P_ON_LOWER = -4,
         /// primal variable is set to its upper bound
         P_ON_UPPER = -2,
         /// primal variable is left free, but not unset
         P_FREE = -1,
         /// primal variable is fixed to both bounds
         P_FIXED = P_ON_UPPER + P_ON_LOWER,
         /// dual variable is left free, but not unset
         D_FREE = 1,
         /// dual variable is set to its upper bound
         D_ON_UPPER = 2,
         /// dual variable is set to its lower bound
         D_ON_LOWER = 4,
         /// dual variable has two bounds
         D_ON_BOTH = D_ON_LOWER + D_ON_UPPER ,
         /// primal or dual variable has no status
         D_UNDEFINED = 8
      };

private:
      friend class SPxBasis;
      DataArray < Status > rowstat;
      DataArray < Status > colstat;
      DataArray < Status > * stat;
      DataArray < Status > * costat;

public:
      /// number of columns.
      int nCols() const
      {
         return colstat.size();
      }
      /// number of rows.
      int nRows() const
      {
         return rowstat.size();
      }

      /// dimension.
      int dim() const
      {
         return stat->size();
      }
      /// codimension.
      int coDim() const
      {
         return costat->size();
      }


      ///
      Status& rowStatus(int i)
      {
         return rowstat[i];
      }
      /// status of row #i#.
      Status rowStatus(int i) const
      {
         return rowstat[i];
      }
      /// array of row #Status#'s.
      const Status* rowStatus(void) const
      {
         return rowstat.get_const_ptr();
      }

      ///
      Status& colStatus(int i)
      {
         return colstat[i];
      }
      /// status of column #i#.
      Status colStatus(int i) const
      {
         return colstat[i];
      }
      /// array of column #Status#'s.
      const Status* colStatus(void) const
      {
         return colstat.get_const_ptr();
      }

      ///
      Status& status(int i)
      {
         return (*stat)[i];
      }
      /// status of variable #i#.
      Status status(int i) const
      {
         return (*stat)[i];
      }
      /// array of variable #Status#'s.
      const Status* status(void) const
      {
         return stat->get_const_ptr();
      }
      ///
      Status& coStatus(int i)
      {
         return (*costat)[i];
      }
      /// status of covariable #i#.
      Status coStatus(int i) const
      {
         return (*costat)[i];
      }
      /// array of covariable #Status#'s.
      const Status* coStatus(void) const
      {
         return costat->get_const_ptr();
      }

      /// reset dimensions.
      void reSize(int rowDim, int colDim);

      ///
      int isConsistent() const;
   };


private:
   Desc thedesc;        // the basis' #Desc#riptor

public:
   ///
   const Desc& desc() const
   {
      return thedesc;
   }
   /// current #Desc# of basis.
   Desc& desc()
   {
      return thedesc;
   }

   /// dual #Status# for the #i#-th column variable of the loaded LP.
   Desc::Status dualColStatus(int i) const;
   /// dual #Status# for the #id#-th column variable of the loaded LP.
   Desc::Status dualStatus(const SPxLP::SPxColId& id) const;

   /// dual #Status# for the #i#-th row variable of the loaded LP.
   Desc::Status dualRowStatus(int i) const;
   /// dual #Status# for the #id#-th row variable of the loaded LP.
   Desc::Status dualStatus(const SPxLP::SPxRowId& id) const;

   /// dual #Status# for the #id#-th variable of the loaded LP.
   Desc::Status dualStatus(const SPxLP::Id& id) const
   {
      return id.isSPxRowId()
             ? dualStatus(SPxLP::SPxRowId(id))
          : dualStatus(SPxLP::SPxRowId(id));
   }

   /** Setup basis.
       Loads a #Desc# to the basis and sets up the basis matrix and
       all vectors accordingly.  The #Desc# must the same number of rows
       and columns as the currently loaded LP.
    */
   virtual void load(const Desc&);            // load basis

   /// read a file in MPS basis format from #in#.
   virtual void readBasis(std::istream& in, NameSet& rowNames, NameSet& colNames);
   //@}


protected:
   /*
       For storing the basis matrix we keep two arrays: Array #theBaseId#
       contains the #Id#s of the basis vectors, and #matrix# the pointers to
       the vectors themselfes. Method #loadMatrixVecs()# serves for loading
       #matrix# according to the #Id#s stored in #theBaseId#. This method must
       be called whenever there is a chance, that the vector pointers may have
       changed due to manipulations of the LP.
    */
   DataArray < SPxLP::Id > theBaseId;      // #Id#s of basic vectors
   DataArray < const SVector* > matrix;         // vectors of the basis matrix
   void loadMatrixVecs();
   int matrixIsSetup;

public:
   /**@name Control Parameters */
   //@{
   /** number of updates before refactorization.
       When a vector of the basis matrix is exchanged by a call to method
       #change()#, the LU factorization of the matrix is updated
       accordingly. However, after atmost #maxUpdates# updates of the
       factorization, it is recomputed in order to regain numerical
       stability and reduce fill in.
    */
   int maxUpdates;

   /** increase of nonzeros before refactorization.
       When the number of nonzeros in LU factorization exceeds
       #nonzeroFactor# times the number of nonzeros of a fresh on, the
       basis matrix is refactorized.
    */
   double nonzeroFactor;
   //@}


   /**@name Inquiry Methods */
   //@{
   ///
   SPxLP::Id& baseId(int i)
   {
      return theBaseId[i];
   }
   /// id of #i#-th basic vector.
   SPxLP::Id baseId(int i) const
   {
      return theBaseId[i];
   }
   /// #i#-th basic vector.
   const SVector& baseVec(int i) const
   {
      return *matrix[i];
   }

protected:
   /*
       The factorization of the matrix is stored in #factor# if #factorized != 0#.
       Otherwise #factor# is undefined.
    */
   SLinSolver* factor;                 // LU factorization of basis matrix
   int factorized;             // 1 if #factor = matrix#$^{-1}$

   /*
       Rank-1-updates to the basis may be performed via method #change#. In
       this case, the factorization is updated, and the following members are
       reset.
    */
   int iterCount;              // number of calls to #change()#
   // since last manipulation
   int updateCount;            // number of calls to #change()#
   // since last #factorize()#
   int nzCount;                // \# of nonzeros in basis matrix
   double nzFac;                  // current nzFactor
   double lastFill;               // fill occured during last factorization

   SPxLP::Id lastin;                 // #lastEntered()#
   SPxLP::Id lastout;                // #lastLeft()#
   int lastidx;                // #lastIndex()#

public:
   /// #Id# of last vector included to the basis.
   SPxLP::Id lastEntered() const
   {
      return lastin;
   }

   /// #Id# of last vector that left the basis.
   SPxLP::Id lastLeft() const
   {
      return lastout;
   }

   /// index in basis where last update was done.
   int lastIndex() const
   {
      return lastidx;
   }

   /// number of basis changes since last refactorization.
   int lastUpdate() const
   {
      return updateCount;
   }

   /// number of basis changes since last #load#.
   int iteration() const
   {
      return iterCount;
   }

   /**@name return loaded solver */
   SoPlex* solver() const
   {
      return theLP;
   }

   //@}

   /**@name Linear Algebra */
   //@{
   /** Basis-vector product.
       Depending on the representation, for a #SPxBasis B#,
       #B.multBaseWith(x)# computes
       \begin{description}
       \item{$x \leftarrow Bx$}    in the columnwise case and
       \item{$x \leftarrow x^TB$}  in the rowwise case.
       \end{description}
       Both can be seen uniformly as multiplying the basis matrix #B# with
       a vector #x# alligned the same way as the {\em vectors} of #B#.
    */
   Vector& multBaseWith(Vector& x) const;

   /** Vector-basis product.
       Depending on the representation, for a #SPxBasis B#,
       #B.multWithBase(x)# computes
       \begin{description}
       \item{$x \leftarrow x^TB$}  in the columnwise case and
       \item{$x \leftarrow Bx$}    in the rowwise case.
       \end{description}
       Both can be seen uniformly as multiplying the basis matrix #B# with
       a vector #x# alligned the same way as the {\em covectors} of #B#.
    */
   Vector& multWithBase(Vector& x) const;

   /// stability of basis matrix.
   double stability() const
   {
      return factor->stability();
   }

   /** Solve linear system with basis matrix.
       Depending on the representation, for a #SPxBasis B#,
       #B.solve(x)# computes
       \begin{description}
       \item{$x \leftarrow B^{-1}x$}       in the columnwise case and
       \item{$x \leftarrow x^TB^{-1}$}     in the rowwise case.
       \end{description}
       Both can be seen uniformly as solving a linear system with the basis
       matrix #B# and a right handside vector #x# aligned the same way as
       the {\em vectors} of #B#.
    */
    void solve2 (Vector& x, Vector& rhs)
    {
       if (!factorized) factorize();
       factor->solve2right(x, rhs);
    }
    ///
    void solve2 (Vector& x, SSVector& rhs)
    {
       if (!factorized) factorize();
       factor->solve2right(x, rhs);
    }
   ///
   void solve2 (SSVector& x, Vector& rhs)
   {
      if (!factorized) factorize();
      factor->solve2right(x, rhs);
   }
   ///
   void solve2 (SSVector& x, SSVector& rhs)
   {
      if (!factorized) factorize();
      factor->solve2right(x, rhs);
   }

    ///
    void solve (Vector& x, const Vector& rhs)
    {
       if (!factorized) factorize();
       factor->solveRight(x, rhs);
    }
    ///
    void solve (Vector& x, const SVector& rhs)
    {
       if (!factorized) factorize();
       factor->solveRight(x, rhs);
    }
    ///
    void solve (SSVector& x, const SVector& rhs)
    {
       if (!factorized) factorize();
       factor->solveRight(x, rhs);
    }
    ///
    void solve (SSVector& x, const Vector& rhs)
    {
       if (!factorized) factorize();
       factor->solveRight(x, rhs);
    }

   ///
   void solve4update(SSVector& x, const SVector& rhs)
   {
      if (!factorized) factorize();
      factor->solveRight4update(x, rhs);
   }
   ///
   void solve4update(SSVector& x, Vector& y,
                     const SVector& rhsx, SSVector& rhsy)
   {
      if (!factorized) factorize();
      factor->solve2right4update(x, y, rhsx, rhsy);
   }

   /** Cosolve linear system with basis matrix.
       Depending on the representation, for a #SPxBasis B#,
       #B.solve(x)# computes
       \begin{description}
       \item{$x \leftarrow x^TB^{-1}$}     in the columnwise case and
       \item{$x \leftarrow B^{-1}x$}       in the rowwise case.
       \end{description}
       Both can be seen uniformly as solving a linear system with the basis
       matrix #B# and a right handside vector #x# alligned the same way as
       the {\em covectors} of #B#.

       If #idx != 0# is given, upon return #idx# contains the indeces of
       the nonzeros of the result vector. #idx# must be allocated to fit
       enough indeces.
    */
   void coSolve2(Vector& x, Vector& rhs)
   {
      if (!factorized) factorize();
      factor->solve2left(x, rhs);
   }
   ///
   void coSolve2(Vector& x, SSVector& rhs)
   {
      if (!factorized) factorize();
      factor->solve2left(x, rhs);
   }
   ///
   void coSolve2(SSVector& x, Vector& rhs)
   {
      if (!factorized) factorize();
      factor->solve2left(x, rhs);
   }
   ///
   void coSolve2(SSVector& x, SSVector& rhs)
   {
      if (!factorized) factorize();
      factor->solve2left(x, rhs);
   }

   ///
   void coSolve(Vector& x, const Vector& rhs)
   {
      if (!factorized) factorize();
      factor->solveLeft(x, rhs);
   }
   ///
   void coSolve(Vector& x, const SVector& rhs)
   {
      if (!factorized) factorize();
      factor->solveLeft(x, rhs);
   }
   ///
   void coSolve(SSVector& x, const SVector& rhs)
   {
      if (!factorized) factorize();
      factor->solveLeft(x, rhs);
   }
   ///
   void coSolve(SSVector& x, const Vector& rhs)
   {
      if (!factorized) factorize();
      factor->solveLeft(x, rhs);
   }

   /// solve 2 systems in 1 call.
   void coSolve(SSVector& x, Vector& y,
                 const SVector& rhsx, SSVector& rhsy)
   {
      if (!factorized) factorize();
      factor->solveLeft(x, y, rhsx, rhsy);
   }
   //@}

   /**@name Modification.
       These methods must be called, after the loaded LP has been modified.
    */
   //@{
protected:
   /** Resize internal arrays.
       When a new LP is loaded, the basis matrix and vectors become invalid
       and possibly also of the wrong dimension. Hence, after loading an
       LP, #reDim()# is called to reset all arrays etc.~accoriding to the
       dimensions of the loaded LP.
    */
   void reDim();

public:
   ///
   void addedRows(int n);
   ///
   void removedRow(int i);
   ///
   void removedRows(int perm[]);
   ///
   void addedCols(int n);
   ///
   void removedCol(int i);
   ///
   void removedCols(int perm[]);
   ///
   void changedRow(int);
   ///
   void changedCol(int);
   ///
   void changedElement(int, int);
   //@}

   /**@name Miscellaneous */
   //@{
   /** perform basis update.
       Changes the $i$-th vector of the basis with the vector associated to
       #id#. This includes:
       \begin{itemize}
       \item       updating the factorization, or recomputing it from
                   scratch by calling \hbox{#factorize()#}
       \item       resetting #lastEnered()#
       \item       resetting #lastIndex()#
       \item       resetting #lastLeft()#
       \item       resetting #lastUpdate()#
       \item       resetting #iterations()#
       \end{itemize}
       The basis #Desc# is {\em not modified}, since #factor()#
       cannot know about how to set up the status of the involved variables
       correctly.\\
       A vector #vec# may be passed for a fast #ETA# update of the LU
       factorization associated to the basis. It must be initialized with
       the solution vector $x$ of the right linear system $Bx = b$ with the
       entering vector as right hand side vetor $b$, where $B$ denotes the
       basis matrix. This can be computed using method #solve(b)#.
       When using #FAST# updates, a vector #upd# may be passed for
       improved performance. It must be initialized by a call to
       #factor->solveRightUpdate# as described in \Ref{SLinSolver}. What
       implementation is hidden behind #FAST# updates, depends on the
       #SLinSolver# implementation class.
    */
   virtual void change(int i, SPxLP::Id& id,
                        const SVector* enterVec, const SSVector* eta = 0);

protected:
   /// refactorize basis instead of updating the factorization?.
   virtual int doFactorize();
   /// factorize the basis matrix.
   virtual void factorize();
   /// Set descriptor representation according to loaded LP
   void setRep();


public:
   /// setup linear solver to use.
   void load(SLinSolver* solver);
   /**
       Loads the #lp# to the basis. This involves resetting all counters to
       0 and setting up a regular default basis consisting of slacks,
       artificial variables or bounds.
    */
   void load(SoPlex* lp);
   ///
   void unLoad()
   {
      theLP = 0;
      setStatus(NO_PROBLEM);
   }


   ///
   int isConsistent() const;
   ///
   SPxBasis();
   virtual ~SPxBasis()
   {}
   //@}
}
;

} // namespace soplex
#endif // _SPXBASIS_H_

//-----------------------------------------------------------------------------
//Emacs Local Variables:
//Emacs mode:c++
//Emacs c-basic-offset:3
//Emacs tab-width:8
//Emacs indent-tabs-mode:nil
//Emacs End:
//-----------------------------------------------------------------------------
