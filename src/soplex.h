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
#pragma ident "@(#) $Id: soplex.h,v 1.9 2001/12/12 10:26:06 bzfkocht Exp $"

/**@file  soplex.h
 * @brief Sequential Objectoriented simPlex
 */
#ifndef _SOPLEX_H_
#define _SOPLEX_H_

#include <assert.h>

#include "cachelpsolver.h"
#include "timer.h"
#include "spxlp.h"
#include "spxbasis.h"
#include "array.h"
#include "random.h"
#include "unitvector.h"
#include "updatevector.h"
#include "subsvector.h"

namespace soplex
{
class SPxPricer;
class SPxRatioTester;
class SPxStarter;
class SPxSimplifier;

/**@brief   Sequential objectoriented simPlex.
   @ingroup Algo

   #SoPlex is an LP solver class using the revised Simplex algorithm. It
   provids two basis representations, namely a column basis and a row basis
   (see #Representation). For both representations, a primal and
   dual algorithm is available (see \Ref{Type}).
 
   In addition, #SoPlex can be custumized with various respects:
   - pricing algorithms using #SPxPricer
   - ratio test using class #SPxRatioTester
   - computation of a start basis using class #SPxStarter
   - preprocessing of the LP using class #SPxSimplifier
   - termination criteria by overriding 
 
   #SoPlex is derived from #SPxLP that is used to store the LP to be solved.
   Hence, the LPs solved with #SoPlex have the general format
 
   \f[
   \begin{array}{rl}
       \hbox{max}   & \mbox{maxObj}^T x                 \\
       \hbox{s.t.}  & \mbox{lhs} \le Ax \le \mbox{rhs}  \\
                    & \mbox{low} \le x  \le \mbox{up}
   \end{array}
   \f]
 
   Also, #SPxLP provide all manipulation methods for the LP. They allow
   #SoPlex to be used within cutting plane algorithms. (@see Programming)
*/
class SoPlex : public CacheLPSolver, public SPxLP, protected SPxBasis
{
   friend class SPxFastRT;

public:
   /**@name Data Types */
   //@{
   /// LP basis representation.
   /** Solving LPs with the Simplex algorithm requires the definition of a
    *  \em basis. A basis can be defined as a set of column vectors or a
    *  set of row vectors building a nonsingular matrix. We will refer to
    *  the first case as the \em columnwise representation and the latter
    *  case will be called the \em rowwise representation.
    *
    *  Type #Representation determines the representation of #SoPlex, i.e.
    *  a columnwise (#COLUMN == 1) or rowwise (#ROW == -1) one.
    */
   enum Representation  
   {
      ROW    = -1,  ///< rowwise representation.
      COLUMN =  1   ///< columnwise representation.
   };

   /// Algorithmic type.
   /** #SoPlex uses the reviesed Simplex algorithm to solve LPs.
    *  Mathematically, one distinguishes the \em primal from the
    *  \em dual algorihm. Algorithmically, these relate to the two
    *  types #ENTER or #LEAVE. How they relate, depends on the chosen
    *  basis representation. This is desribed by the following table:
    *
    *  <TABLE>
    *  <TR><TD>&nbsp;</TD><TD>ENTER </TD><TD>LEAVE </TD></TR>
    *  <TR><TD>ROW   </TD><TD>DUAL  </TD><TD>PRIMAL</TD></TR>
    *  <TR><TD>COLUMN</TD><TD>PRIMAL</TD><TD>DUAL  </TD></TR>
    *  </TABLE>
    */
   enum Type
   {
      /// Entering Simplex.
      /** The Simplex loop for the entering Simplex can be sketched
       *  as follows:
       *  - \em Pricing : Select a variable to #ENTER the basis.
       *  - \em Ratio-Test : Select variable to #LEAVE the
       *    basis such that the basis remains feasible.
       *  - Perform the basis update.
       */
      ENTER = -1,
      /// Leaving Simplex.
      /** The Simplex loop for the leaving Simplex can be sketched
       *  as follows:
       *  - \em Pricing: Select a variable to #LEAVE the basis.
       *  - \em Ratio-Test: Select variable to #ENTER the
       *    basis such that the basis remains priced.
       *  - Perform the basis update.
       */
      LEAVE = 1
   };

   /// Pricing type.
   /** In case of the #ENTER%ing Simplex algorithm, for performance
    *  reasons it may be advisable not to compute and maintain up to
    *  date vectors #pVec() and #test() and instead compute only some
    *  of its elements explicitely. This is constroled by the #Pricing type.
    */
   enum Pricing
   {
      /// Full pricing.
      /** If #FULL pricing in selected for the #ENTER%ing Simplex,
       *  vectors #pVec() and #test() are kept up to date by
       *  #SoPlex. An #SPxPricer only needs to select an #Id such
       *  that the #test() or #coTest() value is < 0.
       */
      FULL,
      /// Partial pricing.
      /** When #PARTIAL pricing in selected for the #ENTER%ing
       *  Simplex, vectors #pVec() and #test() are not set up and
       *  updated by #SoPlex. However, vectors #coPvec() and
       *  #coTest() are still kept up to date by #SoPlex.
       *  An #SPxPricer object needs to compute the values for
       *  #pVec() and #test() itself in order to select an
       *  appropriate pivot with #test() < 0. Methods #computePvec(i)
       *  and #computeTest(i) will assist the used to do so. Note,
       *  that it may be feasable for a pricer to return an #Id with
       *  #test() > 0; such will be rejected by #SoPlex.
       */
      PARTIAL  
   };
   //@}

private:
   Type           theType;     ///< entering or leaving algortihm.
   Pricing        thePricing;  ///< full or partial pricing.
   Representation therep;      ///< row or column representation.
   Timer          theTime;
   int            maxIters;    ///< maximum allowed iterations.
   double         maxTime;     ///< maximum allowed time.
   double         thedelta;
   double         theShift;    ///< shift of r/lhs or objective.
   double         lastShift;   ///< for forcing feasibility.
   int            m_maxCycle;  ///< maximum steps before cycling is detected.
   int            m_numCycle;  ///< actual number of degenerate steps so far.
   bool           initialized; ///< true, if all vectors are setup.

   Vector*        solveVector2;      ///< when 2 systems are to solve at a time
   SSVector*      solveVector2rhs;   ///< when 2 systems are to solve at a time
   Vector*        coSolveVector2;    ///< when 2 systems are to solve at a time
   SSVector*      coSolveVector2rhs; ///< when 2 systems are to solve at a time

protected:
   Array < UnitVector > unitVecs; ///< array of unit vectors
   const SVSet*   thevectors;
   const SVSet*   thecovectors;

   int            nNZEs;          ///< number of nonzero elements
   int            coVecDim;
   Array < Array < SubSVector > > subcovectors;

   DVector        primRhs;     ///< rhs vector for computing the primal vector
   UpdateVector   primVec;     ///< primal vector
   DVector        dualRhs;     ///< rhs vector for computing the dual vector
   UpdateVector   dualVec;     ///< dual vector
   UpdateVector   addVec;      ///< additional vector

   DVector        theURbound;  ///< Upper Row    Feasibility bound
   DVector        theLRbound;  ///< Lower Row    Feasibility bound
   DVector        theUCbound;  ///< Upper Column Feasibility bound
   DVector        theLCbound;  ///< Lower Column Feasibility bound

   /** In entering Simplex algorithm, the ratio test must ensure that all 
    *  \em basic variables remain within their feasibility bounds. To give fast
    *  acces to them, the bounds of basic variables are copied into the
    *  following two vectors.
    */
   DVector        theUBbound;  ///< Upper Basic Feasibility bound
   DVector        theLBbound;  ///< Lower Basic Feasibility bound

   DVector*       theFrhs;
   UpdateVector*  theFvec;

   DVector*       theCoPrhs;
   UpdateVector*  theCoPvec;
   UpdateVector*  thePvec;

   /** set on #theCoPvec or #thePvec */
   UpdateVector*  theRPvec;    ///< row pricing vector
   /** set on #thePvec or #theCoPvec */
   UpdateVector*  theCPvec;    /// column pricing vector

   // The following vectors serve for the virtualization of shift bounds
   DVector*       theUbound;      ///< Upper bound for vars
   DVector*       theLbound;      ///< Lower bound for vars
   DVector*       theCoUbound;    ///< Upper bound for covars
   DVector*       theCoLbound;    ///< Lower bound for covars

   // The following vectors serve for the virtualization of testing vectors
   DVector        theCoTest;
   DVector        theTest;

   int            leaveCount;     ///< number of LEAVE iterations
   int            enterCount;     ///< number of ENTER iterations

   SPxPricer*      thepricer;
   SPxRatioTester* theratiotester;
   SPxStarter*     thestarter;
   SPxSimplifier*  thesimplifier;

public:
   /// return the current basis representation.
   Representation rep() const
   {
      return therep;
   }

   /// return current #Type.
   Type type() const
   {
      return theType;
   }

   /// return current #Pricing.
   Pricing pricing() const
   {
      return thePricing;
   }

   /**@name Setup
    *  Before solving an LP with an instance of #SoPlex, 
    *  the following steps must be performed:
    *
    *  -# Load the LP by copying an external LP or reading it from an
    *     input stream.
    *  -# Setup the pricer to use by loading an #SPxPricer object
    *     (only neccessary, if not done in a previous call).
    *  -# Setup the ratio test method to use by loading an #SPxRatioTester 
    *     object (only neccessary, if not done in a previous call).
    *  -# Setup the linear system solver to use by loading an
    *     #SLinSolver object (only neccessary, if not done in a previous call).
    *  -# Optionally setup an LP simplifier by loading an SPxSimplifier object.
    *  -# Optionally setup an start basis generation method by loading an
    *     SPxStarter object.
    *  -# Optionally setup a start basis by loading a SPxBasis::Desc object.
    *  -# Optionally switch to another basis Representation by calling
    *     method #setRep().
    *  -# Optionally switch to another algorithm #Type by calling
    *     method #setType().
    *
    *  Now the solver is ready for execution. If the loaded LP is to be solved
    *  again from scratch, this can be done with method #reLoad(). Finally,
    *  #clear() removes the LP from the solver.
    */
   //@{
   /// read LP from input stream using operator syntax.
   friend std::istream& operator>>(std::istream& is, SoPlex& lp)
   {
      lp.read(is);
      return is;
   }
   /// read LP from input stream.
   virtual void read(std::istream& in, 
      NameSet* rowNames = 0, NameSet* colNames = 0);

   /// copy LP.
   void load(const SPxLP& LP);
   /// setup linear solver to use.
   void load(SLinSolver* slu);
   /// setup pricer to use.
   void load(SPxPricer*);
   /// setup ratio-tester to use.
   void load(SPxRatioTester*);
   /// setup starting basis generator to use.
   void load(SPxStarter*);
   /// setup simplifier to use.
   void load(SPxSimplifier*);
   /// set a start basis.
   void load(const SPxBasis::Desc&);

   /// set #ROW or #COLUMN representation.
   void setRep (int rep);
   /// set #LEAVE or #ENTER algorithm.
   void setType(Type tp);
   /// set #FULL or #PARTIAL pricing.
   void setPricing(Pricing pr);

   /// reload LP.
   virtual void reLoad();

   /// load LP from \p filename in MPS or LPF format.
   void readFile(char* filename);

   /// dump loaded LP to \p filename in LPF format.
   void dumpFile(char* filename) const;

   /// clear all data in solver.
   void clear();
   //@}

   /**@name Solving LPs */
   //@{
   /// solve loaded LP.
   /** Solves the loaded LP by processing the Simplex iteration until
    *  the termination criteria is fullfilled (see #terminate()). 
    *  The #SPxStatus of the solver will indicate the reason for termination.
    */
   LPSolver::Status solve();

   /// #Status of basis.
   LPSolver::Status status() const;

   /// current objective value.
   /**@return Objective value of the current solution vector 
    *         (see #getPrimal()).
    */
   virtual double value() const;

   /// get solution vector for primal variables.
   /** This method returns the #Status of the basis.
    *  If it is #REGULAR or better,
    *  the primal solution vector of the current basis will be copied
    *  to the argument \p vector. Hence, \p vector must be of dimension
    *  #nCols().
    */
   LPSolver::Status getPrimal (Vector& vector) const;

   /// get vector of slack variables.
   /** This method returns the #Status of the basis.
    *  If it is #REGULAR or better,
    *  the slack variables of the current basis will be copied
    *  to the argument \p vector. Hence, \p vector must be of dimension
    *  #nRows().
    *
    *  @warning Because #SoPlex supports range constraints as its
    *     default, slack variables are defined in a nonstandard way:
    *     Let \i x be the current solution vector and \i A the constraint
    *     matrix. Then the vector of slack variables is defined as
    *     \f$s = Ax\f$.
    */
   LPSolver::Status getSlacks (Vector& vector) const;

   /// get current solution vector for dual variables.
   /** This method returns the #Status of the basis.
    *  If it is #REGULAR or better,
    *  the vector of dual variables of the current basis will be copied
    *  to the argument \p vector. Hence, \p vector must be of dimension
    *  #nRows().
    *
    *  @warning Even though mathematically, each range constraint would
    *     account for two dual variables (one for each inequaility), only
    *     #nRows() dual variables are setup via the following
    *     construction: Given a range constraint, there are three possible
    *     situations:
    *     - None of its inequalities is tight: The dual variables
    *       for both are 0. However, when shifting (see below)
    *       occurs, it may be set to a value other than 0, which
    *       models a perturbed objective vector.
    *     - Both of its inequalities are tight: In this case the
    *       range constraint models an equality and we adopt the
    *       standard definition.
    *     - One of its inequalities is tight while the other is not:
    *       In this case only the dual variable for the tight
    *       constraint is given with the standard definition, while
    *       the other constraint is implicitely set to 0.
    */
   LPSolver::Status getDual (Vector& vector) const;

   /// get current solution vector for dual variables.
   /** This method returns the \Ref{Status} of the basis.
    *  If it is #REGULAR or better,
    *  the vector of reduced costs of the current basis will be copied
    *  to the argument \p vector. Hence, \p vector must be of dimension
    *  #nCols().
    *
    *  Let \i d denote the vector of dual variables, as defined above,
    *  and \i A the LPs constraint matrix. Then the reduced cost vector
    *  \i r is defined as \f$r^T = c^T - d^TA\f$.
    */
   LPSolver::Status getRdCost (Vector& vector) const;

   /// Termination criterion.
   /** This method is called in each Simplex iteration to determine, if
    *  the algorithm is to terminate. In this case a nonzero value is
    *  returned.
    *
    *  This method is declared virtual to allow for implementation of
    *  other stopping criteria or using it as callback method within the
    *  Simplex loop, by overriding the method in a derived class.
    *  However, all implementations must terminate with the
    *  statement \c return #SoPlex::terminate(), if no own termination
    *  criteria is encountered.
    *
    *  Note, that the Simplex loop stopped even when #terminate()
    *  returns 0, if the LP has been solved to optimality (i.e. no
    *  further pricing succeeds and no shift is present).
    */
   virtual int terminate ();
   //@}

   /**@name Control Parameters */
   //@{

   /// values \f$|x| < \epsilon\f$ are considered to be 0.
   double epsilon() const
   {
      return primVec.delta().epsilon;
   }
   /// set parameter \p epsilon.
   void setEpsilon(double eps);

   /// allowed bound violation for optimal Solution.
   /** When all vectors do not violate their bounds by more than \f$\delta\f$,
    *  the basis is considered optimal.
    */
   double delta() const
   {
      return thedelta;
   }
   /// set parameter \p delta.
   void setDelta(double d);

   /** #SoPlex consideres a Simplex step as degenerate, if the
    *  steplength does not exceed #epsilon. Cycling occurs, if only
    *  degenerate steps are taken. To prevent this situation, #SoPlex
    *  perturbs the problem such that nondegenerate steps are ensured.
    *
    *  maxCycle() controls, how agressive such perturbation is
    *  performed, since no more than #maxCycle() degenerate steps are
    *  accepted before perturbing the LP. The current number of consequtive
    *  degenerate steps is counted in variable numCycle().
    */
   /// maximum number of degenerate simplex steps before we detect cycling.
   int maxCycle() const 
   {
      return m_maxCycle;
   }
   /// actual number of degenerate simplex steps encountered so far.
   int numCycle() const 
   {
      return m_numCycle;
   }
   //@}

   /**@name LP Modification
    *  These methods are the overridden counterparts of the base class
    *  #SPxLP. See there for a more detailed documentation of these
    *  methods.
    */
   //@{
private:
   ///
   void localAddRows(int start);
   ///
   void localAddCols(int start);

protected:
   ///
   virtual void addedRows(int n);
   ///
   virtual void addedCols(int n);

   ///
   virtual void doRemoveRow(int i);
   ///
   virtual void doRemoveRows(int perm[]);
   ///
   virtual void doRemoveCol(int i);
   ///
   virtual void doRemoveCols(int perm[]);

public:
   ///
   virtual void changeObj(const Vector& newObj);
   ///
   virtual void changeObj(int i, double newVal);
   ///
   virtual void changeObj(SPxLP::SPxColId p_id, double p_newVal)
   {
      changeObj(number(p_id), p_newVal);
   }
   ///
   virtual void changeLower(const Vector& newLower);
   ///
   virtual void changeLower(int i, double newLower);
   ///
   virtual void changeLower(SPxLP::SPxColId p_id, double p_newLower)
   {
      changeLower(number(p_id), p_newLower);
   }
   ///
   virtual void changeUpper(const Vector& newUpper);
   ///
   virtual void changeUpper(int i, double newUpper);
   ///
   virtual void changeUpper(SPxLP::SPxColId p_id, double p_newUpper)
   {
      changeUpper(number(p_id), p_newUpper);
   }
   ///
   virtual void changeBounds(const Vector& newLower, const Vector& newUpper);
   ///
   virtual void changeBounds(int i, double newLower, double newUpper);
   ///
   virtual void changeBounds(
      SPxLP::SPxColId p_id, double p_newLower, double p_newUpper)
   {
      changeBounds(number(p_id), p_newLower, p_newUpper);
   }
   ///
   virtual void changeLhs(const Vector& newLhs);
   ///
   virtual void changeLhs(int i, double newLhs);
   ///
   virtual void changeLhs(SPxLP::SPxRowId p_id, double p_newLhs)
   {
      changeLhs(number(p_id), p_newLhs);
   }
   ///
   virtual void changeRhs(const Vector& newRhs);
   ///
   virtual void changeRhs(int i, double newRhs);
   ///
   virtual void changeRhs(SPxLP::SPxRowId p_id, double p_newRhs)
   {
      changeRhs(number(p_id), p_newRhs);
   }
   ///
   virtual void changeRange(const Vector& newLhs, const Vector& newRhs);
   ///
   virtual void changeRange(int i, double newLhs, double newRhs);
   ///
   virtual void changeRange(
      SPxLP::SPxRowId p_id, double p_newLhs, double p_newRhs)
   {
      changeRange(number(p_id), p_newLhs, p_newRhs);
   }
   ///
   virtual void changeRow(int i, const LPRow& newRow);
   ///
   virtual void changeRow(SPxLP::SPxRowId p_id, const LPRow& p_newRow)
   {
      changeRow(number(p_id), p_newRow);
   }
   ///
   virtual void changeCol(int i, const LPCol& newCol);
   ///
   virtual void changeCol(SPxLP::SPxColId p_id, const LPCol& p_newCol)
   {
      changeCol(number(p_id), p_newCol);
   }
   ///
   virtual void changeElement(int i, int j, double val);
   ///
   virtual void changeElement(
      SPxLP::SPxRowId rid, SPxLP::SPxColId cid, double val)
   {
      changeElement(number(rid), number(cid), val);
   }
   ///
   virtual void changeSense(SPxSense sns);
   //@}

   /// dimension of basis matrix.
   int dim() const
   {
      return thecovectors->num();
   }
   /// codimension.
   int coDim() const
   {
      return thevectors->num();
   }
   /// number of row \p p_id.
   int number(SPxRowId p_id) const
   {
      return SPxLP::number(p_id);
   }
   /// number of column \p p_id.
   int number(SPxColId p_id) const
   {
      return SPxLP::number(p_id);
   }
   /// number of column \p p_id.
   int number(Id p_id) const
   {
      return SPxLP::number(p_id);
   }

   /**@name Variables and Covariables
    *  Class #SPxLP introduces #Id%s to identify row or column data of
    *  an LP. #SoPlex uses this concept to access data with respect to the
    *  chosen representation.
    */
   //@{
   /// id of \p i 'th vector.
   /** The \p i 'th #id is the \p i 'th #SPxRowId for a rowwise and the
    *  \p i 'th #SPxColId for a columnwise basis represenation. Hence,
    *  0 <= i < #coDim().
    */
   Id id(int i) const
   {
      if (rep() == ROW)
      {
         SPxLP::SPxRowId rid = SPxLP::rId(i);
         return Id(rid);
      }
      else
      {
         SPxLP::SPxColId cid = SPxLP::cId(i);
         return Id(cid);
      }
   }

   /// id of \p i 'th covector.
   /** The \p i 'th #coId() is the \p i 'th #SPxColId for a rowwise and the
    *  \p i 'th #SPxRowId for a columnwise basis represenation. Hence,
    *  0 <= i < #dim().
    */
   Id coId(int i) const
   {
      if (rep() == ROW)
      {
         SPxLP::SPxColId cid = SPxLP::cId(i);
         return Id(cid);
      }
      else
      {
         SPxLP::SPxRowId rid = SPxLP::rId(i);
         return Id(rid);
      }
   }

   /// Is \p p_id an Id ?
   /** This method returns wheather or not \p p_id identifies a vector
    *  with respect to the chosen representation.
    */
   int isId(SPxLP::Id p_id) const
   {
      return p_id.info * therep > 0;
   }

   /// Is \p p_id a CoId.
   /** This method returns wheather or not \p p_id identifies a coVector
    *  with respect to the chosen representation.
    */
   int isCoId(SPxLP::Id p_id) const
   {
      return p_id.info * therep < 0;
   }
   //@}

   /**@name Vectors and Covectors */
   //@{

protected:
   /// ???
   int sortLP (int pe, int nPes);
   /// ???
   void splitLP(int pe, int nPes);
   /// ???
   virtual void splitLP();

public:
   /// \p i 'th vector.
   /**@return a reference to the \p i 'th, 0 <= i < #coDim(), vector of
    *         the loaded LP (with respect to the chosen representation).
    */
   const SVector& vector(int i) const
   {
      return (*thevectors)[i];
   }

   ///
   const SVector& vector(const SPxLP::SPxRowId& rid) const
   {
      assert(rid.isValid());
      return (rep() == ROW)
             ? (*thevectors)[number(rid)]
          : static_cast<const SVector&>(unitVecs[number(rid)]);
   }
   ///
   const SVector& vector(const SPxLP::SPxColId& cid) const
   {
      assert(cid.isValid());
      return (rep() == COLUMN)
             ? (*thevectors)[number(cid)]
          : static_cast<const SVector&>(unitVecs[number(cid)]);
   }

   /// vector associated to \p p_id.
   /**@return Returns a reference to the vector of the loaded LP corresponding
    *  to \p id (with respect to the chosen representation). If \p p_id is
    *  an id, a vector of the constraint matrix is returned, otherwise
    *  the corresponding unit vector (of the slack variable or bound
    *  inequality) is returned.
    *  @todo The implementation does not exactly look exactly it will do
    *        what is promised in the describtion.
    */
   const SVector& vector(const Id& p_id) const
   {
      assert(p_id.isValid());

      return p_id.isSPxRowId()
         ? vector(SPxLP::SPxRowId(p_id))
         : vector(SPxLP::SPxColId(p_id));
   }

   /// \p i 'th covector of LP.
   /**@return a reference to the \p i 'th, 0 <= i < #dim(), covector of
    *  the loaded LP (with respect to the chosen representation).
    */
   const SVector& coVector(int i) const
   {
      return (*thecovectors)[i];
   }
   ///
   const SVector& coVector(const SPxLP::SPxRowId& rid) const
   {
      assert(rid.isValid());
      return (rep() == COLUMN)
             ? (*thecovectors)[number(rid)]
          : static_cast<const SVector&>(unitVecs[number(rid)]);
   }
   ///
   const SVector& coVector(const SPxLP::SPxColId& cid) const
   {
      assert(cid.isValid());
      return (rep() == ROW)
             ? (*thecovectors)[number(cid)]
          : static_cast<const SVector&>(unitVecs[number(cid)]);
   }
   /// coVector associated to \p p_id.
   /**@return a reference to the covector of the loaded LP
    *  corresponding to \p p_id (with respect to the chosen
    *  representation). If \p p_id is a coid, a covector of the constraint
    *  matrix is returned, otherwise the corresponding unit vector is
    *  returned.
    */
   const SVector& coVector(const Id& p_id) const
   {
      assert(p_id.isValid());
      return p_id.isSPxRowId()
             ? coVector(SPxLP::SPxRowId(p_id))
          : coVector(SPxLP::SPxColId(p_id));
   }
   /// return \p i 'th unit vector.
   const SVector& unitVector(int i) const
   {
      return unitVecs[i];
   }
   //@}

   /**@name Variable status
    *  The Simplex basis assigns a #SPxBasis::Desc::Status to each
    *  variable and covariable. Depending on the representation, the status
    *  indicates that the corresponding vector is in the basis matrix or not.
    */
   //@{
   /// #Status of \p i 'th variable.
   SPxBasis::Desc::Status status(int i);

   /// #Status of \p i 'th covariable.
   SPxBasis::Desc::Status coStatus(int i);

   /// does \p stat describe a basic index ?
   int isBasic(SPxBasis::Desc::Status stat) const
   {
      return (stat * rep() > 0);
   }

   /// is the \p p_id 'th vector basic ?
   int isBasic(SPxLP::Id p_id) const
   {
      assert(p_id.isValid());
      return p_id.isSPxRowId()
             ? isBasic(SPxLP::SPxRowId(p_id))
          : isBasic(SPxLP::SPxColId(p_id));
   }

   /// is the \p rid 'th vector basic ?
   int isBasic(SPxLP::SPxRowId rid) const
   {
      return isBasic(desc().rowStatus(number(rid)));
   }

   /// is the \p cid 'th vector basic ?
   int isBasic(SPxLP::SPxColId cid) const
   {
      return isBasic(desc().colStatus(number(cid)));
   }

   /// is the \p i 'th row vector basic ?
   int isRowBasic(int i) const
   {
      return isBasic(desc().rowStatus(i));
   }

   /// is the \p i 'th column vector basic ?
   int isColBasic(int i) const
   {
      return isBasic(desc().colStatus(i));
   }

   /// is the \p i 'th vector basic ?
   int isBasic(int i) const
   {
      return isBasic(desc().status(i));
   }

   /// is the \p i 'th covector basic ?
   int isCoBasic(int i) const
   {
      return isBasic(desc().coStatus(i));
   }
   //@}

   /// feasibility vector.
   /** This method return the \em feasibility vector. If it satisfies its
    *  bound, the basis is called feasible (independently of the chosen
    *  representation). The feasibility vector has dimension #dim().
    *
    *  For the entering Simplex, #fVec is kept within its bounds. In
    *  contrast to this, the pricing of the leaving Simplex selects an
    *  element of #fVec, that violates its bounds.
    */
   UpdateVector& fVec() const
   {
      return *theFvec;
   }
   /// right-hand side vector for #fVec.
   /** The feasibility vector is computed by solving a linear system with the
    *  basis matrix. The right-hand side vector of this system is refferd to as
    *  \em feasibility \em, right-hand \em side \em vector #fRhs().
    *
    *  For a row basis, #fRhs() is the objective vector (ignoring shifts).
    *  For a column basis, it is the sum of all nonbasic vectors scaled by
    *  the factor of their bound.
    */
   const Vector& fRhs() const
   {
      return *theFrhs;
   }
   ///
   const Vector& ubBound() const
   {
      return theUBbound;
   }
   /// upper bound for #fVec.
   /** This method returns the upper bound for the feasibility vector.
    *  It may only be called for the #ENTER%ing Simplex.
    *  
    *  For the #ENTER%ing Simplex algorithms, the feasibility vector is
    *  maintained to fullfill its bounds. As #fVec itself, also its
    *  bounds depend on the chosen representation. Furhter, they may
    *  need to be shifted (see below).
    */
   Vector& ubBound()
   {
      return theUBbound;
   }
   ///
   const Vector& lbBound() const
   {
      return theLBbound;
   }
   /// lower bound for #fVec.
   /** This method returns the lower bound for the feasibility vector.
    *  It may only be called for the #ENTER%ing Simplex.
    *
    *  For the #ENTER%ing Simplex algorithms, the feasibility vector is
    *  maintained to fullfill its bounds. As #fVec itself, also its
    *  bound depend on the chosen representation. Further, they may
    *  need to be shifted (see below).
    */
   Vector& lbBound()
   {
      return theLBbound;
   }

   /// Violations of #fVec.
   /** For the leaving Simplex algorithm, pricing involves selecting a
    *  variable from #fVec that violates its bounds that is to leave
    *  the basis. When a #SPxPricer is called to select such a
    *  leaving variable, #fTest() contains the vector of violations:
    *  For #fTest()[i] < 0, the \c i 'th basic variable violates one of
    *  its bounds by the given value. Otherwise no bound is violated.
    */
   const Vector& fTest() const
   {
      assert(type() == LEAVE);
      return theCoTest;
   }

   /// copricing vector.
   /** The copricing vector #coPvec along with the pricing vector
    *  #pVec are used for pricing in the #ENTER%ing Simplex algorithm,
    *  i.e. one variable is selected, that violates its bounds. In
    *  contrast to this, the #LEAVE%ing Simplex algorithm keeps both
    *  vectors within their bounds.
    */
   UpdateVector& coPvec() const
   {
      return *theCoPvec;
   }

   /// Right-hand side vector for #coPvec.
   /** Vector #coPvec is computed by solving a linear system with the
    *  basis matrix and #coPrhs as the right-hand side vector. For
    *  column basis representation, #coPrhs is build up of the
    *  objective vector elements of all basic variables. For a row
    *  basis, it consists of the thight bounds of all basic
    *  constraints.
    */
   const Vector& coPrhs() const
   {
      return *theCoPrhs;
   }

   ///
   const Vector& ucBound() const
   {
      assert(theType == LEAVE);
      return *theCoUbound;
   }
   /// upper bound for #coPvec.
   /** This method returns the upper bound for #coPvec. It may only be
    *  called for the leaving Simplex algorithm.
    *
    *  For the leaving Simplex algorithms #coPvec is maintained to
    *  fullfill its bounds. As #coPvec itself, also its bound depend
    *  on the chosen representation. Further, they may need to be
    *  shifted (see below).
    */
   Vector& ucBound()
   {
      assert(theType == LEAVE);
      return *theCoUbound;
   }

   ///
   const Vector& lcBound() const
   {
      assert(theType == LEAVE);
      return *theCoLbound;
   }
   /// lower bound for #coPvec.
   /** This method returns the lower bound for #coPvec. It may only be
    *  called for the leaving Simplex algorithm.
    *
    *  For the leaving Simplex algorithms #coPvec is maintained to
    *  fullfill its bounds. As #coPvec itself, also its bound depend
    *  on the chosen representation. Further, they may need to be
    *  shifted (see below).
    */
   Vector& lcBound()
   {
      assert(theType == LEAVE);
      return *theCoLbound;
   }

   /// violations of #coPvec.
   /** In entering Simplex pricing selects checks vectors #coPvec()
    *  and #pVec()# for violation of its bounds. #coTest() contains
    *  the violations for #coPvec() which are indicated by a negative
    *  value. I.e. if #coTest(i) < 0, the \p i 'th element of #coPvec()
    *  is violated by #-coTest(i).
    */
   const Vector& coTest() const
   {
      assert(type() == ENTER);
      return theCoTest;
   }
   /// pricing vector.
   /** The pricing vector #pVec is the product of #coPvec with the
    *  constraint matrix. As #coPvec, also #pVec is maintained within
    *  its bound for the leaving Simplex algorithm, while the bounds
    *  are tested for the entering Simplex. #pVec is of dimension
    *  #dim(). Vector #pVec() is only up to date for #LEAVE%ing
    *  Simplex or #FULL pricing in #ENTER%ing Simplex.
    */
   UpdateVector& pVec() const
   {
      return *thePvec;
   }
   ///
   const Vector& upBound() const
   {
      assert(theType == LEAVE);
      return *theUbound;
   }
   /// upper bound for #pVec.
   /** This method returns the upper bound for #pVec. It may only be
    *  called for the leaving Simplex algorithm.
    *
    *  For the leaving Simplex algorithms #pVec is maintained to
    *  fullfill its bounds. As #pVec itself, also its bound depend
    *  on the chosen representation. Further, they may need to be
    *  shifted (see below).
    */
   Vector& upBound()
   {
      assert(theType == LEAVE);
      return *theUbound;
   }

   ///
   const Vector& lpBound() const
   {
      assert(theType == LEAVE);
      return *theLbound;
   }
   /// lower bound for #pVec.
   /** This method returns the lower bound for #pVec. It may only be
    *  called for the leaving Simplex algorithm.
    *
    *  For the leaving Simplex algorithms #pVec is maintained to
    *  fullfill its bounds. As #pVec itself, also its bound depend
    *  on the chosen representation. Further, they may need to be
    *  shifted (see below).
    */
   Vector& lpBound()
   {
      assert(theType == LEAVE);
      return *theLbound;
   }

   /// Violations of #pVec.
   /** In entering Simplex pricing selects checks vectors #coPvec()
    *  and #pVec() for violation of its bounds. Vector #test()
    *  contains the violations for #pVec(), i.e.~if #test(i) < 0,
    *  the i'th element of #pVec() is violated by #test(i).
    *  Vector #test() is only up to date for #FULL pricing.
    */
   const Vector& test() const
   {
      assert(type() == ENTER);
      return theTest;
   }

   /// compute and return #pVec(i).
   double computePvec(int i);
   /// compute entire #pVec().
   void computePvec();
   /// compute and return #test(i) in #ENTER%ing Simplex.
   double computeTest(int i);
   /// compute test vector in #ENTER%ing Simplex.
   void computeTest();

   /**@name Shifting
    *  The task of the ratio test (implemented in #SPxRatioTester classes)
    *  is to select a variable for the basis update, such that the basis
    *  remains priced (i.e. both, the pricing and copricing vectors satisfy
    *  their bounds) or feasible (i.e. the feasibility vector satisfies its
    *  bounds). However, this can lead to numerically instable basis matrices
    *  or -- after accumulation of various errors -- even to a singular basis
    *  matrix.
    *
    *  The key to overcome this problem is to allow the basis to become "a
    *  bit" infeasible or unpriced, in order provide a better choice for the
    *  ratio test to select a stable variable. This is equivalent to enlarging
    *  the bounds by a small amount. This is referred to as \em shifting.
    *
    *  These methods serve for shifting feasibility bounds, either in order
    *  to maintain numerical stability or initially for computation of
    *  phase 1. The sum of all shifts applied to any bound is stored in
    *  #theShift.
    *
    *  The following methods are used to shift individual bounds. They are
    *  mainly intended for stable implenentations of #SPxRatioTester.
    */
   //@{
   /// Perform initial shifting to optain an feasible or pricable basis.
   void shiftFvec();
   /// Perform initial shifting to optain an feasible or pricable basis.
   void shiftPvec();

   /// shift \p i 'th #ubBound to \p to.
   void shiftUBbound(int i, double to)
   {
      assert(theType == ENTER);
      theShift += to - theUBbound[i];
      theUBbound[i] = to;
   }
   /// shift \p i 'th #lbBound to \p to.
   void shiftLBbound(int i, double to)
   {
      assert(theType == ENTER);
      theShift += theLBbound[i] - to;
      theLBbound[i] = to;
   }
   /// shift \p i 'th #upBound to \p to.
   void shiftUPbound(int i, double to)
   {
      assert(theType == LEAVE);
      theShift += to - (*theUbound)[i];
      (*theUbound)[i] = to;
   }
   /// shift \p i 'th #lpBound# to \p to.
   void shiftLPbound(int i, double to)
   {
      assert(theType == LEAVE);
      theShift += (*theLbound)[i] - to;
      (*theLbound)[i] = to;
   }
   /// shift \p i 'th #ucBound# to \p to.
   void shiftUCbound(int i, double to)
   {
      assert(theType == LEAVE);
      theShift += to - (*theCoUbound)[i];
      (*theCoUbound)[i] = to;
   }
   /// shift \p i 'th #lcBound# to \p to.
   void shiftLCbound(int i, double to)
   {
      assert(theType == LEAVE);
      theShift += (*theCoLbound)[i] - to;
      (*theCoLbound)[i] = to;
   }
   ///
   void testBounds() const;

   /// total current shift amount.
   virtual double shift() const
   {
      return theShift;
   }
   /// remove shift as much as possible.
   virtual void unShift(void);

private:
   ///
   void perturbMin(
      const UpdateVector& vec, Vector& low, Vector& up, double eps,
      int start = 0, int incr = 1);
   ///
   void perturbMax(
      const UpdateVector& vec, Vector& low, Vector& up, double eps,
      int start = 0, int incr = 1);
   ///
   double perturbMin(const UpdateVector& uvec,
      Vector& low, Vector& up, double eps, double delta,
      const SPxBasis::Desc::Status* stat, int start, int incr);
   ///
   double perturbMax(const UpdateVector& uvec,
      Vector& low, Vector& up, double eps, double delta,
      const SPxBasis::Desc::Status* stat, int start, int incr);
   //@}

   /**@name The Simplex Loop
    *  We now present a set of methods that may be usefull when implementing
    *  own #SPxPricer or #SPxRatioTester classes. Here is, how
    *  #SoPlex will call methods from its loaded #SPxPricer and
    *  #SPxRatioTester.
    *  
    *  For the entering Simplex:
    *    -# #SPxPricer::selectEnter()
    *    -# #SPxRatioTester::selectLeave()
    *    -# #SPxPricer::entered4()
    *  
    *  For the leaving Simplex:
    *    -# #SPxPricer::selectLeave()
    *    -# #SPxRatioTester::selectEnter()
    *    -# #SPxPricer::left4()
    */
   //@{
public:
   /// Setup vectors to be solved within Simplex loop.
   /** Load vector \p y to be #solve%d with the basis matrix during the
    *  #LEAVE Simplex. The system will be solved after #SoPlex%'s call
    *  to #SPxRatioTester.  The system will be solved along with
    *  another system. Solving two linear system at a time has
    *  performance advantages over solving the two linear systems
    *  seperately.
    */
   void setup4solve(Vector* p_y, SSVector* p_rhs)
   {
      assert(type() == LEAVE);
      solveVector2    = p_y;
      solveVector2rhs = p_rhs;
   }
   /// Setup vectors to be cosolved within Simplex loop.
   /** Load vector \p y to be #coSolve%#d with the basis matrix during
    *  the #ENTER Simplex. The system will be solved after #SoPlex%'s
    *  call to #SPxRatioTester.  The system will be solved along
    *  with another system. Solving two linear system at a time has
    *  performance advantages over solving the two linear systems
    *  seperately.
    */
   void setup4coSolve(Vector* p_y, SSVector* p_rhs)
   {
      assert(type() == ENTER);
      coSolveVector2    = p_y;
      coSolveVector2rhs = p_rhs;
   }
   /// maximal infeasibility of basis
   /** This method is called for prooving optimality. Since it is
    *  possible, that some stable implementation of class
    *  #SPxRatioTester yielded a slightly infeasible (or unpriced)
    *  basis, this must be checked before terminating with an optimal
    *  solution.
    */
   virtual double maxInfeas() const;

   /// Return current basis.
   /**@note The basis can be used to solve linear systems or use
    *  any other of its (const) methods.  It is, however, encuraged
    *  to use methods #setup4solve() and #setup4coSolve() for solving
    *  systems, since this is likely to have perfomance advantages.
    */
   const SPxBasis& basis() const
   {
      return *this;
   }
   ///
   SPxBasis& basis()
   {
      return *this;
   }
   /// return loaded #SPxPricer.
   const SPxPricer* pricer() const
   {
      return thepricer;
   }
   /// return loaded #SLinSolver.
   const SLinSolver* slinSolver() const
   {
      return SPxBasis::factor;
   }
   /// return loaded #SPxRatioTester.
   const SPxRatioTester* ratiotester() const
   {
      return theratiotester;
   }
   /// return loaded #SPxStarter.
   const SPxStarter* starter() const
   {
      return thestarter;
   }
   /// return loaded #SPxSimplifier.
   const SPxSimplifier* simplifier() const
   {
      return thesimplifier;
   }

protected:
   /// Factorize basis matrix.
   virtual void factorize();

private:
   int leave(int i);
   int enter(Id& id);

   /// test coVector #i# with status #stat#.
   double coTest(int, SPxBasis::Desc::Status) const;
   /// compute coTest vector.
   void computeCoTest();
   /// recompute coTest vector.
   void updateCoTest();

   /// test vector #i# with status #stat#.
   double test(int i, SPxBasis::Desc::Status stat) const;
   /// recompute test vector.
   void updateTest();

   /// compute basis feasibility test vector.
   void computeFtest();
   /// update basis feasibility test vector.
   void updateFtest();
   //@}

   /**@name Parallelization
    *  In this section we present the methods, that are provided in order to
    *  allow a parallel version to be implemented as a derived class, thereby
    *  inheriting most of the code of #SoPlex#.
    *
    *  @par Initialization
    *  These methods are used to setup all the vectors used in the Simplex
    *  loop, that where described in the previous sectios.
    */
   //@{
protected:
   /// has the internal data been initialized?
   /** As long as an instance of #SoPlex is not #initialized, no member
    *  contains setup data. Initialization is performed via method
    *  #init().  Afterwards all data structures are kept up to date (even
    *  for all manipulation methods), until #unInit() is called. However,
    *  some manipulation methods call #unInit()# themselfs.
    */
   bool isInitialized() const
   {
      return initialized;
   }
public:
   /// intialize data structures.
   /** If #SoPlex is not #isInitialized(), method solve calls
    *  #init() to setup all vectors and internal data structures.
    *  Most of the other methods within this section are called by
    *  #init().
    *
    *  Derived classes should add the initialization of additional
    *  data structures by overriding this method. Don't forget,
    *  however to call #SoPlex::init().
    */
   virtual void init();

protected:
   /// unintialize data structures.
   virtual void unInit()
   {
      initialized = false;
   }

   /// reset dimensions of vectors according to loaded LP.
   virtual void reDim();
   /// compute feasibility vector from scratch.
   void computeFrhs();
   ///
   virtual void computeFrhsXtra();
   ///
   virtual void computeFrhs1(const Vector&, const Vector&);
   ///
   void computeFrhs2(const Vector&, const Vector&);
   /// compute #theCoPrhs for entering Simplex.
   virtual void computeEnterCoPrhs();
   ///
   void computeEnterCoPrhs4Row(int i, int n);
   ///
   void computeEnterCoPrhs4Col(int i, int n);
   /// compute #theCoPrhs for leaving Simplex.
   virtual void computeLeaveCoPrhs();
   ///
   void computeLeaveCoPrhs4Row(int i, int n);
   ///
   void computeLeaveCoPrhs4Col(int i, int n);

   /// Compute part of objective value.
   /** This method is called from #value() in order to compute the part of
    *  the objective value resulting form nonbasic variables for #COLUMN
    *  #Representation.
    */
   double nonbasicValue() const;

   /// Get pointer to the \p id 'th vector
   virtual const SVector* enterVector(const Id& p_id)
   {
      assert(p_id.isValid());
      return p_id.isSPxRowId() 
         ? &vector(SPxRowId(p_id)) : &vector(SPxColId(p_id));
   }
   ///
   virtual void getLeaveVals(int i,
      SPxBasis::Desc::Status& leaveStat, Id& leaveId,
      double& leaveMax, double& leavebound, int& leaveNum);
   ///
   virtual void getLeaveVals2(double leaveMax, Id enterId,
      double& enterBound, double& newUBbound,
      double& newLBbound, double& newCoPrhs);
   ///
   virtual void getEnterVals(Id id, double& enterTest,
      double& enterUB, double& enterLB, double& enterVal, double& enterMax,
      double& enterPric, SPxBasis::Desc::Status& enterStat, double& enterRO);
   ///
   virtual void getEnterVals2(int leaveIdx, 
      double enterMax, double& leaveBound);
   ///
   virtual void ungetEnterVal(Id enterId, SPxBasis::Desc::Status enterStat,
      double leaveVal, const SVector& vec);
   ///
   virtual void rejectEnter(Id enterId,
      double enterTest, SPxBasis::Desc::Status enterStat);
   ///
   virtual void rejectLeave(int leaveNum, Id leaveId,
      SPxBasis::Desc::Status leaveStat, const SVector* newVec = 0);
   ///
   virtual void setupPupdate(void);
   ///
   virtual void doPupdate(void);
   ///
   virtual void clearUpdateVecs(void);
   ///
   virtual void perturbMinEnter(void);
   /// perturb basis bounds.
   virtual void perturbMaxEnter(void);
   ///
   virtual void perturbMinLeave(void);
   /// perturb nonbasic bounds.
   virtual void perturbMaxLeave(void);

   //@}
   /*  The following methods serve for initializing the bounds for dual or
    *  primal Simplex algorithm of entering or leaving type.
    */
   //@{
   void clearDualBounds(SPxBasis::Desc::Status, double&, double&);

   void setDualColBounds();
   void setDualRowBounds();
   void setPrimalBounds();

   void setEnterBound4Col(int, int);
   void setEnterBound4Row(int, int);
   virtual void setEnterBounds();

   void setLeaveBound4Row(int i, int n);
   void setLeaveBound4Col(int i, int n);
   virtual void setLeaveBounds();
   //@}

public:
   /**@name Derived from #LPSolver 
    */
   //@{

   /// adjust conditions for termination.
   void setTermination(double value = LPSolver::infinity,
      double time = -1, int iteration = -1);

   /// get adjusted conditions for termination.
   virtual void getTermination(double* value = 0 ,
      double* time = 0, int* iteration = 0) const;

   /// get objective value of current solution.
   double objValue() const
   {
      return value();
   }

   /// get all results of last solve.
   LPSolver::Status getResult(double* value = 0, Vector* primal = 0,
      Vector* slacks = 0, Vector* dual = 0, Vector* reduCost = 0) const;

   /// set #LPSolver#s basis.
   void setBasis(const signed char rows[], const signed char cols[]);

   /// get current basis.
   LPSolver::Status getBasis(signed char rows[], signed char cols[]) const;

   /// get number of iterations of current solution.
   int iterations() const
   {
      return basis().iteration();
   }

   /// time spent in last call to method #solve().
   double time() const
   {
      return theTime.userTime();
   }
#if 0
   /// add \p row to #LPSolver%s LP.
   void addRow(const LPRow& row)
   {
      SPxLP::addRow(row);
   }
   /// add \p row to #LPSolver%s LP.
   void addRow(LPSolver::RowId& p_id, const LPRow& p_row)
   {
      SPxLP::addRow( *reinterpret_cast<SPxRowId*>(&p_id), p_row);
   }

   /// add all #LPRow%s of \p p_set to #LPSolver%s LP.
   void addRows(const LPRowSet& p_set)
   {
      SPxLP::addRows(p_set);
   }
   /// add all #LPRow%s of \p p_set to #LPSolver%s LP.
   void addRows(LPSolver::RowId p_id[], const LPRowSet& p_set)
   {
      SPxLP::addRows( reinterpret_cast<SPxRowId*>(p_id), p_set);
   }

   /// add \p p_col to #LPSolver%s LP.
   void addCol(const LPCol& col)
   {
      SPxLP::addCol(col);
   }
   /// add \p p_col to #LPSolver%s LP.
   void addCol(LPSolver::ColId& p_id, const LPCol& p_col)
   {
      SPxLP::addCol(*reinterpret_cast<SPxColId*>(&p_id), p_col);
   }

   /// add all #LPCol%s of \p set to #LPSolver%s LP.
   void addCols(const LPColSet& p_set)
   {
      SPxLP::addCols(p_set);
   }
   /// add all #LPCol%s of \p set to #LPSolver%s LP.
   void addCols(LPSolver::ColId p_id[], const LPColSet& p_set)
   {
      SPxLP::addCols(reinterpret_cast<SPxColId*>(p_id), p_set);
   }


   /// remove \p i 'th row.
   void removeRow(int i)
   {
      SPxLP::removeRow(i);
   }
   /// remove row with #RowId \p id.
   void removeRow(LPSolver::RowId p_id)
   {
      SPxLP::removeRow(*reinterpret_cast<SPxRowId*>(&p_id));
   }

   /// remove \p i 'th column.
   void removeCol(int i)
   {
      SPxLP::removeCol(i);
   }
   /// remove column with #LPSolver::ColId \p id.
   void removeCol(LPSolver::ColId p_id)
   {
      SPxLP::removeCol(*reinterpret_cast<SPxColId*>(&p_id));
   }

   /// remove \p n rows.
   void removeRows(LPSolver::RowId p_id[], int p_n, int p_perm[] = 0)
   {
      SPxLP::removeRows(reinterpret_cast<SPxRowId*>(p_id), p_n, p_perm);
   }
   /// remove \p n rows.
   void removeRows(int nums[], int n, int perm[] = 0)
   {
      SPxLP::removeRows(nums, n, perm);
   }
   /// remove multiple rows.
   void removeRows(int perm[])
   {
      SPxLP::removeRows(perm);
   }
   /// remove rows from \p start to \p end (including both).
   void removeRowRange(int start, int end, int perm[] = 0)
   {
      SPxLP::removeRowRange(start, end, perm);
   }

   /// remove \p n columns.
   void removeCols(LPSolver::ColId p_id[], int p_n, int p_perm[] = 0)
   {
      SPxLP::removeCols(reinterpret_cast<SPxColId*>(p_id), p_n, p_perm);
   }
   /// remove \p n columns.
   void removeCols(int nums[], int n, int perm[] = 0)
   {
      SPxLP::removeCols(nums, n, perm);
   }
   /// remove multiple columns.
   void removeCols(int perm[])
   {
      SPxLP::removeCols(perm);
   }
   /// remove columns from \p start to \p end (including both).
   void removeColRange(int start, int end, int perm[] = 0)
   {
      SPxLP::removeColRange(start, end, perm);
   }
#endif

   /// change objective value to variable with #ColId \p id.
   void changeObj(LPSolver::ColId p_id, double p_newVal)
   {
      changeObj(*reinterpret_cast<SPxColId*>(&p_id), p_newVal);
   }

   /// change lower bound of variable with #ColId \p id.
   void changeLower(LPSolver::ColId p_id, double p_newLower)
   {
      changeLower(*reinterpret_cast<SPxColId*>(&p_id), p_newLower);
   }

   /// change \p id 'th upper bound.
   void changeUpper(LPSolver::ColId p_id, double p_newUpper)
   {
      changeUpper(*reinterpret_cast<SPxColId*>(&p_id), p_newUpper);
   }

   /// change \p id 'th lower and upper bound.
   void changeBounds(LPSolver::ColId p_id, 
      double p_newLower, double p_newUpper)
   {
      changeBounds(*reinterpret_cast<SPxColId*>(&p_id), p_newLower, p_newUpper);
   }

   /// change #id#'th lhs value.
   void changeLhs(LPSolver::RowId p_id, double p_newLhs)
   {
      changeLhs(*reinterpret_cast<SPxRowId*>(&p_id), p_newLhs);
   }

   /// change #id#'th rhs value.
   void changeRhs(LPSolver::RowId p_id, double p_newRhs)
   {
      changeRhs(*reinterpret_cast<SPxRowId*>(&p_id), p_newRhs);
   }

   /// change #id#'th lhs and rhs value.
   void changeRange(LPSolver::RowId p_id, double p_newLhs, double p_newRhs)
   {
      changeRange(*reinterpret_cast<SPxRowId*>(&p_id), p_newLhs, p_newRhs);
   }

   /// change #id#'th row of LP.
   void changeRow(LPSolver::RowId p_id, const LPRow& p_newRow)
   {
      changeRow(*reinterpret_cast<SPxRowId*>(&p_id), p_newRow);
   }

   /// change #id#'th column of LP.
   void changeCol(LPSolver::ColId p_id, const LPCol& p_newCol)
   {
      changeCol(*reinterpret_cast<SPxColId*>(&p_id), p_newCol);
   }

   /// change LP element (#rid#, #cid#).
   void changeElement(LPSolver::RowId p_rid, LPSolver::ColId p_cid, double p_val)
   {
      changeElement(*reinterpret_cast<SPxRowId*>(&p_rid),
                    *reinterpret_cast<SPxColId*>(&p_cid),
                    p_val);
   }

   /// change optimization sense to #sns#.
   void changeSense(LPSolver::Sense p_sns)
   {
      changeSense(SPxSense(int(p_sns)));
   }


   /// get \p i 'th row.
   void getRow(int p_i, LPRow& p_row) const
   {
      SPxLP::getRow(p_i, p_row);
   }

   /// get #id#'th row.
   void getRow(LPSolver::RowId p_id, LPRow& p_row) const
   {
      SPxLP::getRow(*reinterpret_cast<SPxRowId*>(&p_id), p_row);
   }

   /// get rows #start# .. #end#.
   void getRows(int p_start, int p_end, LPRowSet& p_set) const
   {
      SPxLP::getRows(p_start, p_end, p_set);
   }

   /// return const \p i 'th row if available.
   const SVector& rowVector(int i) const
   {
      return SPxLP::rowVector(i);
   }

   /// return const #id#'th row if available.
   const SVector& rowVector(LPSolver::RowId p_id) const
   {
      return SPxLP::rowVector(*reinterpret_cast<SPxRowId*>(&p_id));
   }

   /// return const lp's rows if available.
   const LPRowSet& rows() const
   {
      return *lprowset();
   }

   /// get \p i 'th column.
   void getCol(int p_i, LPCol& p_column) const
   {
      SPxLP::getCol(p_i, p_column);
   }

   /// get #id#'th column.
   void getCol(LPSolver::ColId p_id, LPCol& p_column) const
   {
      SPxLP::getCol(*reinterpret_cast<SPxColId*>(&p_id), p_column);
   }

   /// get columns #start# .. #end#.
   void getCols(int p_start, int p_end, LPColSet& p_set) const
   {
      SPxLP::getCols(p_start, p_end, p_set);
   }

   /// return const \p i 'th col if available.
   const SVector& colVector(int i) const
   {
      return SPxLP::colVector(i);
   }

   /// return const #id#'th col if available.
   const SVector& colVector(LPSolver::ColId p_id) const
   {
      return SPxLP::colVector(*reinterpret_cast<SPxColId*>(&p_id));
   }

   /// return const lp's cols if available.
   const LPColSet& cols() const
   {
      return *lpcolset();
   }

   /// \p i 'th value of objective vector.
   double obj(int i) const
   {
      return SPxLP::obj(i);
   }

   /// #id#'th value of objective vector.
   double obj(LPSolver::ColId p_id) const
   {
      return SPxLP::obj(*reinterpret_cast<SPxColId*>(&p_id));
   }
   ///
   virtual const Vector& obj() const
   {
      return CacheLPSolver::obj();
   }

   /// copy objective vector to #obj#.
   void getObj(Vector& p_obj) const
   {
      SPxLP::getObj(p_obj);
   }

   /// \p i 'th lower bound.
   double lower(int i) const
   {
      return SPxLP::lower(i);
   }

   /// #id#'th lower bound.
   double lower(LPSolver::ColId p_id) const
   {
      return SPxLP::lower(*reinterpret_cast<SPxColId*>(&p_id));
   }

   /// copy lower bound vector to #low#.
   void getLower(Vector& lw) const
   {
      lw = SPxLP::lower();
   }

   /// return const lower bound vector.
   const Vector& lower() const
   {
      return SPxLP::lower();
   }


   /// \p i 'th upper bound.
   double upper(int i) const
   {
      return SPxLP::upper(i);
   }

   /// #id#'th upper bound.
   double upper(LPSolver::ColId p_id) const
   {
      return SPxLP::upper(*reinterpret_cast<SPxColId*>(&p_id));
   }

   /// copy upper bound vector to #up#.
   void getUpper(Vector& upp) const
   {
      upp = SPxLP::upper();
   }

   /// return const upper bound vector.
   const Vector& upper() const
   {
      return SPxLP::upper();
   }


   /// \p i 'th lhs value.
   double lhs(int i) const
   {
      return SPxLP::lhs(i);
   }

   /// #id#'th lhs value.
   double lhs(LPSolver::RowId p_id) const
   {
      return SPxLP::lhs(*reinterpret_cast<SPxRowId*>(&p_id));
   }

   /// copy lhs value vector to #lhs#.
   void getLhs(Vector& p_lhs) const
   {
      p_lhs = SPxLP::lhs();
   }

   /// return const lhs vector.
   const Vector& lhs() const
   {
      return SPxLP::lhs();
   }


   /// \p i 'th rhs value.
   double rhs(int i) const
   {
      return SPxLP::rhs(i);
   }

   /// #id#'th rhs value.
   double rhs(LPSolver::RowId p_id) const
   {
      return SPxLP::rhs(*reinterpret_cast<SPxRowId*>(&p_id));
   }

   /// copy rhs value vector to #rhs#.
   void getRhs(Vector& p_rhs) const
   {
      p_rhs = SPxLP::rhs();
   }

   /// return const rhs vector.
   const Vector& rhs() const
   {
      return SPxLP::rhs();
   }


   /// optimization sense.
   LPSolver::Sense sense() const
   {
      return LPSolver::Sense(spxSense());
   }


   /// number of rows of loaded LP.
   int nofRows() const
   {
      return nRows();
   }
   /// number of columns of loaded LP.
   int nofCols() const
   {
      return nCols();
   }
   /// number of nonzeros of loaded LP.
   int nofNZEs() const;

   /// #RowId# of \p i 'th inequality.
   LPSolver::RowId rowId(int i) const
   {
      SPxRowId p_id = rId(i);
      return *reinterpret_cast<LPSolver::RowId*>(&p_id);
   }
   /// #ColId# of \p i 'th column.
   LPSolver::ColId colId(int i) const
   {
      SPxColId p_id = cId(i);
      return *reinterpret_cast<LPSolver::ColId*>(&p_id);
   }

   /// number of row #id#.
   int number(LPSolver::RowId p_id) const
   {
      return number(*reinterpret_cast<SPxRowId*>(&p_id));
   }
   /// number of column #id#.
   int number(LPSolver::ColId p_id) const
   {
      return number(*reinterpret_cast<SPxColId*>(&p_id));
   }

   /// does #LPSolver# have row with #id#.
   int has(LPSolver::RowId p_id) const
   {
      return number(p_id) >= 0;
   }
   /// does #LPSolver# have column with #id#.
   int has(LPSolver::ColId p_id) const
   {
      return number(p_id) >= 0;
   }
   //@}

   /**@name Miscellaneous */
   //@{
   /// assignment operator.
   SoPlex& operator=(const SoPlex& base);

   /// copy constructor.
   SoPlex(const SoPlex& base);

   /// default constructor.
   SoPlex(Type type = LEAVE, Representation rep = ROW,
           SPxPricer* pric = 0, SPxRatioTester* rt = 0,
           SPxStarter* start = 0, SPxSimplifier* simple = 0
        );

   /// check consistency.
   int isConsistent() const;
   //@}

private:
   void testVecs();

   double cacheProductFactor;
};

} // namespace soplex
#endif // _SOPLEX_H_

//-----------------------------------------------------------------------------
//Emacs Local Variables:
//Emacs mode:c++
//Emacs c-basic-offset:3
//Emacs tab-width:8
//Emacs indent-tabs-mode:nil
//Emacs End:
//-----------------------------------------------------------------------------
