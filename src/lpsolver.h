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
#pragma ident "@(#) $Id: lpsolver.h,v 1.6 2001/12/12 10:26:06 bzfkocht Exp $"


/**@file  lpsolver.h
 * @brief Interface class for LP solver.
 */
#ifndef _LPSOLVER_H_
#define _LPSOLVER_H_

#include "svector.h"
#include "vector.h"
#include "lpcolset.h"
#include "lprowset.h"

namespace soplex
{
/**@brief   Interface class for LP solver.
   @ingroup Algebra

   Class #LPSolver provides a generic interface to simplex-type linear program
   solvers. After construction, an LP solver is available with an empty problem
   loaded to it. The solver can be loaded with a nontrival problem by means of
   methods #load() or #read(). The latter method reads an LP from a file,
   while the former passes an LP as argument.
   
   The problem loaded to an #LPSolver is considered to be of the form:
   \f[
     \begin{array}{rl}
     \hbox{min/max}  & c^T x                 \\
     \hbox{s.t.}     & l \le Ax \le r        \\
                     & w \le  x \le u
     \end{array}
   \f]
   Vector \f$c\f$ will be referred to as the objective Vector, \f$l\f$ and
   \f$r\f$ as the left and right hand side vectors, respectively, and \f$w\f$
   and \f$u\f$ as lower and upper bounds, respectively.
   
   The problem loaded to an #LPSolver can be solved with either the primal or
   dual simplex algorithm, by calling method #solve(). A termination criterion
   for the solution method can be set with method #setTermination().
   
   An #LPSolver has a #Status associated to it, indicating what the
   solver knows about the problem loaded to it. Most non-#const
   methods change the status of the solver.
   
   There are various methods available for changing the problem loaded to an
   #LPSolver. They allow to add, change or remove rows or columns of the
   constraint matrix or to modify the variable's upper and lower bound, the
   objective vector or the left- and right hand sides and senses of the
   constraint inequalities.
   
   The rows and columns of the LP loaded to an #LPSolver are numbered
   implicitly from 0 to #nofRows()-1 and 0 to #nofCols()-1,
   respectively. When rows or columns are added to the LP, they
   receive indices #nofRows(), #nofRows()+1, ... or #nofCols(),
   #nofCols()+1, ..., respectively. The indices of the first rows and
   columns in the LP remain unchanged. When adding rows or columns, no
   precautions with respect to memory management are required,
   i.e. all memory is reset to the required size automatically.
   
   When removing rows or columns from the LP loaded to an #LPSolver the
   remaining rows and columns are renumbered. However, all rows or columns with
   number lower than the first one removed keep their number unchanged.
   
   Class #LPSolver provides two means of dealing with the renumbering of the
   remaining rows or columns.  First, each removal method provides an optional
   parameter, where an array of indices (#int* #idx) can be passed. Upon return
   of the remove method, this number array contains the permutations due to the
   removal, i.e. #idx[i] is the new number of the row or column, that had
   number \p i \em before the last remove method was called. Hence, the number
   arrays must be at least of size #nofRows() or #nofCols(), respectively.
   
   The second concept for dealing with the renumbering of row or column indices
   in removal methods is the provision of #RowId%s and #ColId%s. An #LPSolver
   associates a unique #RowId to each row and a #ColId to each column of the
   loaded LP.  The Id of a row or column remains fixed throughout the time, the
   row or column belongs to the LP. In most methods, #RowId%s and #ColId%s may
   be used instead of indices.  */
class LPSolver
{
public:
   /**@name Datatypes */
   //@{
   /// solver status.
   /** This enumeration type describes the amount of information, the solver
       knows about its LP.
   */
   enum Status
   {
      UNKNOWN = 0,   ///< nothing known on loaded problem.
      UNBOUNDED,     ///< loaded problem is unbounded.
      INFEASIBLE,    ///< loaded problem is infeasible.
      PRIMAL,        ///< primal (not yet optimal) solution available.
      DUAL,          ///< dual (not yet optimal) solution available.
      SOLVED,        ///< loaded problem has been solved.
      ERROR          ///< an error occurred.
   };

   /// optimization sense.
   enum Sense
   {
      MAXIMIZE = 1,  ///< maximization problem.
      MINIMIZE = -1  ///< minimization problem.
   };

   /// unique id to access columns in an #LPSolver.
   class ColId : public DataKey
   {
   };
   //   int id;        ///< column id number.

   /// unique id to access rows in an #LPSolver.
   class RowId : public DataKey
   {
   };
   //   int id;        ///< row id number

   /// status of variables.
   /** A basis assigns a #VarStatus to each variable of the loaded
       LP. The names for the variable status come from a standard column
       basis representation. These are:
       - #ON_UPPER   if the variable is nonbasic and on its upper bound
       - #ON_LOWER   if the variable is nonbasic and on its upper bound
       - #FIXED      if the variable is nonbasic and its upper bound is equal
                     to its lower bound
       - #FIXED      if the variable is nonbasic and unbounded
       - #BASIC      if the variable is basic

       For slack variables the interpretation of #ON_UPPER means, that the
       upper bound of the inequality is tight.  Similarily, if the status
       of a slack variables is #ON_LOWER the lower bound of the inequality
       is tight.
    */
   enum VarStatus
   {
      ON_UPPER,   ///< variable set to its upper bound.
      ON_LOWER,   ///< variable set to its lower bound.
      FIXED,      ///< variable fixed to identical bounds.
      ZERO,       ///< free variable fixed to zero.
      BASIC       ///< variable is basic.
   };

   static const double infinity;   ///< value used as \f$\infty\f$.
   //@}


   /**@name Solving LPs */
   //@{
   /// solves current LP with the simplex method.
   virtual Status solve() = 0;

   /// sets current basis.
   /** Each variable is set to the status specified in the arrays \p rows
       and \p cols, which must be of adequate size.
       @see        #VarStatus
   */
   virtual void setBasis(const signed char rows[],
                         const signed char cols[]) = 0;

   /// adjusts conditions for termination.
   /**@todo the termination criterion "value" is not yet implemented! */
   virtual void setTermination(double value = infinity,
                               double time = -1,
                               int iteration = -1) = 0;

   /// gets adjusted conditions for termination.
   virtual void getTermination(double* value = 0,
                               double* time = 0,
                               int* iteration = 0) const = 0;
   //@}


   /**@name Accessing Computational Results */
   //@{
   /// returns objective value of current solution.
   virtual double objValue() const = 0;

   /// gets current solution vector for primal variables.
   virtual Status getPrimal(Vector& vector) const = 0;

   /// returns const solution vector for primal variables.
   virtual const Vector& primal() const = 0;

   /// gets current solution vector for dual variables.
   virtual Status getDual(Vector& vector) const = 0;

   /// returns const solution vector for dual variables.
   virtual const Vector& dual() const = 0;

   /// gets current vector of slack variables.
   virtual Status getSlacks(Vector& vector) const = 0;

   /// returns const vector of slack variables.
   virtual const Vector& slacks() const = 0;

   /// gets vector of reduced costs.
   virtual Status getRdCost(Vector& vector) const = 0;

   /// returns const vector of reduced costs.
   virtual const Vector& rdCost() const = 0;

   /// gets all results of last solve.
   virtual Status getResult(double* value = 0,
                            Vector* primal = 0,
                            Vector* slacks = 0,
                            Vector* dual = 0,
                            Vector* reduCost = 0) const = 0;

   /// Get current basis.
   /** The status information for each variable and row is copied into the
       arrays \p rows and \p cols, which must be of adequate size.
       @see    #VarStatus
   */
   virtual Status getBasis(signed char rows[],
                           signed char cols[]) const = 0;

   /// returns status for primal vars if available.
   virtual const signed char* rowBasis() const = 0;

   /// returns status for dual vars if available.
   virtual const signed char* colBasis() const = 0;

   /// gets number of iterations of current solution.
   virtual int iterations() const = 0;

   /// gets time for computing current solution.
   virtual double time() const = 0;
   //@}


   /**@name Loading LPs */
   //@{
   /// unloads current problem and re-initializes lpsolver.
   virtual void clear() = 0;

   /// loads LP from \p filename in MPS or LP format.
   virtual void readFile(char* filename) = 0;

   /// dumps loaded LP to \p filename in LP format.
   virtual void dumpFile(char* filename) const = 0;
   //@}
#if 0
   /**@name Adding Rows and Columns
      Each row and column in an #LPSolver is associated to #RowId and
      #ColId, respectively, which remains unchanged, as long as the row or
      column remains in the loaded LP. Theses #RowId%s and #ColId%s are
      assigned by the #LPSolver when rows and columns are added. Hence, all
      methods for adding rows and columns come with two signatures. One of
      them, provides a first parameter for returning the assigned #RowId(s)
      or #ColId(s), respectively.
   */
   //@{
   ///
   virtual void addRow(const LPRow& row) = 0;
   /// adds \p row to #LPSolver's LP.
   virtual void addRow(RowId& id, const LPRow& row) = 0;

   ///
   virtual void addRows(const LPRowSet& set) = 0;
   /// adds all #LPRow%s of \p set to #LPSolver's LP.
   virtual void addRows(RowId id[], const LPRowSet& set) = 0;

   ///
   virtual void addCol(const LPCol& col) = 0;
   /// adds \p col to #LPSolver's LP.
   virtual void addCol(ColId& id, const LPCol& col) = 0;

   ///
   virtual void addCols(const LPColSet& set) = 0;
   /// adds all #LPCol%s of \p set to #LPSolver's LP.
   virtual void addCols(ColId id[], const LPColSet& set) = 0;
   //@}


   /**@name Removing Rows and Columns
      Either single or multiple rows or columns may be removed in one method
      call. In general, both lead to renumbering of the remaining rows or
      columns.
      
      If only one row or column is removed at a time, the situation is simple.
      Either it is the last row or column, already, or the last one is moved
      to the  position (number) of the removed row or column.
      
      When multiple rows or columns are removed with one method
      invocation, the renumbering scheme is not specified. Instead,
      all such removal methods provide an additional parameter \p
      perm. If nonzero, \p perm must point to an array of #int%s of
      (at least) the size of the number of rows or columns. After
      termination, #perm[i] is the permuted number of the former \p i
      'th row or column, or < 0, if the \p i 'th row or column has
      been removed.  
    */
   //@{
   /// removes \p i 'th row.
   virtual void removeRow(int i) = 0;
   /// removes row with #RowId \p id.
   virtual void removeRow(RowId id) = 0;

   /// removes \p i 'th column.
   virtual void removeCol(int i) = 0;
   /// removes column with #ColId \p id.
   virtual void removeCol(ColId id) = 0;

   /// removes \p n rows given by #RowId%s \p id.
   virtual void removeRows(RowId id[], int n, int perm[] = 0) = 0;
   /// removes \p n rows given by their row numbers \p nums.
   virtual void removeRows(int nums[], int n, int perm[] = 0) = 0;

   /// removes multiple rows.
   /** Removes all #LPRow%s with a number \p i such that \p perm[i] < 0. Upon
       completion, \p perm[i] >= 0 indicates the new number where the \p i 'th
       #LPRow has been moved to due to this removal. Note, that \p perm
       must point to an array of at least #rowNumber() #int%s.
   */
   virtual void removeRows(int perm[]) = 0;
   /// removes rows from \p start to \p end (including both).
   virtual void removeRowRange(int start, int end, int perm[] = 0) = 0;

   /// removes \p n columns given by #ColId%s \p id.
   virtual void removeCols(ColId id[], int n, int perm[] = 0) = 0;
   /// removes \p n columns given by their column numbers \p nums.
   virtual void removeCols(int nums[], int n, int perm[] = 0) = 0;

   /// removes multiple columns.
   /** Removes all #LPCol%s with a number \p i such that \p perm[i] < 0. Upon
       completion, \p perm[i] >= 0 indicates the new number where the \p i 'th
       #LPCol has been moved to due to this removal. Note, that \p perm
       must point to an array of at least #colNumber() #int%s.
   */
   virtual void removeCols(int perm[]) = 0;
   /// removes columns from \p start to \p end (including both).
   virtual void removeColRange(int start, int end, int perm[] = 0) = 0;
   //@}
#endif

   /**@name Manipulating the LP */
   //@{
   /// makes \p newObj the new objective vector.
   virtual void changeObj(const Vector& newObj) = 0;

   /// changes \p i 'th objective value to \p newVal.
   virtual void changeObj(int i, double newVal) = 0;

   /// changes the objective value of column with #ColId \p id to \p newVal.
   virtual void changeObj(ColId id, double newVal) = 0;

   /// makes \p newLower the new lower bound vector.
   virtual void changeLower(const Vector& newLower) = 0;

   /// changes \p i 'th lower bound value to \p newLower.
   virtual void changeLower(int i, double newLower) = 0;

   /// changes the lower bound of column with #ColId \p id to \p newLower.
   virtual void changeLower(ColId id, double newLower) = 0;

   /// makes \p newUpper the new upper bound vector.
   virtual void changeUpper(const Vector& newUpper) = 0;

   /// changes \p i 'th upper bound value to \p newUpper.
   virtual void changeUpper(int i, double newUpper) = 0;

   /// changes the upper bound of column with #ColId \p id to \p newUpper.
   virtual void changeUpper(ColId id, double newUpper) = 0;

   /// makes \p newLower the new lower bound and \p newUpper 
   /// the new upper bound vector.
   virtual void changeBounds(
      const Vector& newLower, const Vector& newUpper) = 0;

   /// changes the bounds of column \p i to \p newLower and \p newUpper.
   virtual void changeBounds(int i, double newLower, double newUpper) = 0;

   /// changes the bounds of column with #ColId \p id to \p newLower 
   /// and \p newUpper.
   virtual void changeBounds(ColId id, double newLower, double newUpper) = 0;

   /// makes \p newLhs the new left hand side vector for constraints.
   virtual void changeLhs(const Vector& newLhs) = 0;

   /// changes the left hand side value of row \p i to \p newLhs.
   virtual void changeLhs(int i, double newLhs) = 0;

   /// changes the left hand side value of row with #RowId \p id to \p newLhs.
   virtual void changeLhs(RowId id, double newLhs) = 0;

   /// makes \p newRhs the new right hand side vector for constraints.
   virtual void changeRhs(const Vector& newRhs) = 0;

   /// changes the right hand side value of row \p i to \p newRhs.
   virtual void changeRhs(int i, double newRhs) = 0;

   /// changes the right hand side value of row with #RowId \p id to \p newRhs.
   virtual void changeRhs(RowId id, double newRhs) = 0;

   /// makes \p newLhs the new left and \p newRhs 
   /// the new right hand side vector for constraints.
   virtual void changeRange(const Vector& newLhs, const Vector& newRhs) = 0;

   /// changes the range values of row \p i to \p newLhs and \p newRhs.
   virtual void changeRange(int i, double newLhs, double newRhs) = 0;

   /// changes the range values of row with #Rowid \p id 
   /// to \p newLhs and \p newRhs.
   virtual void changeRange(RowId id, double newLhs, double newRhs) = 0;

   /// replaces row \p i with \p newRow.
   virtual void changeRow(int i, const LPRow& newRow) = 0;

   /// replaces row with #RowId \p id with \p newRow.
   virtual void changeRow(RowId id, const LPRow& newRow) = 0;

   /// replaces column \p i with \p newCol.
   virtual void changeCol(int i, const LPCol& newCol) = 0;

   /// replaces column with #ColId \p id with \p newCol.
   virtual void changeCol(ColId id, const LPCol& newCol) = 0;

   /// changes element (\p i, \p j) of coefficient matrix A 
   /// to new value \p val.
   virtual void changeElement(int i, int j, double val) = 0;

   /// changes element of coefficient matrix A corresponding 
   /// to #RowId \p rid and #ColId \p cid to new value \p val.
   virtual void changeElement(RowId rid, ColId cid, double val) = 0;

   /// change optimization sense to \p sns.
   virtual void changeSense(Sense sns) = 0;
   //@}

   /**@name Accessing Loaded LP */
   //@{
   /// gets \p i 'th row.
   virtual void getRow(int i, LPRow& row) const = 0;

   /// gets row with #RowId \p id.
   virtual void getRow(RowId id, LPRow& row) const = 0;

   /// gets rows \p start through \p end.
   virtual void getRows(int start, int end, LPRowSet& set) const = 0;

   /// returns the coefficient vector of \p i 'th row.
   virtual const SVector& rowVector(int i) const = 0;

   /// returns the coefficient vector of row with #RowId \p id.
   virtual const SVector& rowVector(RowId id) const = 0;

   /// returns the set of all rows of the LP.
   virtual const LPRowSet& rows() const = 0;

   /// gets \p i 'th column.
   virtual void getCol(int i, LPCol& column) const = 0;

   /// gets column with #ColId \p id.
   virtual void getCol(ColId id, LPCol& column) const = 0;

   /// gets columns \p start through \p end.
   virtual void getCols(int start, int end, LPColSet& set) const = 0;

   /// returns the coefficient vector of \p i 'th column.
   virtual const SVector& colVector(int i) const = 0;

   /// returns the coefficient vector of column with #ColId \p id.
   virtual const SVector& colVector(ColId id) const = 0;

   /// returns the set of all columns of the LP.
   virtual const LPColSet& cols() const = 0;

   /// returns the left hand side of row \p i.
   virtual double lhs(int i) const = 0;

   /// returns the left hand side of row with #RowId \p id.
   virtual double lhs(RowId id) const = 0;

   /// copies left hand side vector to \p lhs.
   virtual void getLhs(Vector& lhs) const = 0;

   /// returns left hand side vector.
   virtual const Vector& lhs() const = 0;

   /// returns the right hand side of row \p i.
   virtual double rhs(int i) const = 0;

   /// returns the right hand side of row with #RowId \p id.
   virtual double rhs(RowId id) const = 0;

   /// copies right hand side vector to \p rhs.
   virtual void getRhs(Vector& rhs) const = 0;

   /// returns right hand side vector.
   virtual const Vector& rhs() const = 0;

   /// returns the objective value of column \p i.
   virtual double obj(int i) const = 0;

   /// returns the objective value of column with #ColId \p id.
   virtual double obj(ColId id) const = 0;

   /// copies objective vector to \p obj.
   virtual void getObj(Vector& obj) const = 0;

   /// returns objective vector.
   virtual const Vector& obj() const = 0;

   /// returns the lower bound of column \p i.
   virtual double lower(int i) const = 0;

   /// returns the lower bound of column with #ColId \p id.
   virtual double lower(ColId id) const = 0;

   /// copies lower bound vector to \p low.
   virtual void getLower(Vector& low) const = 0;

   /// returns lower bound vector.
   virtual const Vector& lower() const = 0;

   /// returns the upper bound of column \p i.
   virtual double upper(int i) const = 0;

   /// returns the upper bound of column with #ColId \p id.
   virtual double upper(ColId id) const = 0;

   /// copies upper bound vector to \p up.
   virtual void getUpper(Vector& up) const = 0;

   /// returns upper bound vector.
   virtual const Vector& upper() const = 0;

   /// returns optimization sense.
   virtual Sense sense() const = 0;
   //@}

   /**@name Inquiry */
   //@{
   /// returns the solvers #Status.
   virtual Status status() const = 0;

   /// returns the number of columns of loaded LP.
   virtual int nofCols() const = 0;
   /// returns the number of rows of loaded LP.
   virtual int nofRows() const = 0;
   /// returns the number of nonzero coefficients of loaded LP.
   virtual int nofNZEs() const = 0;

   /// returns the number of row with #RowId \p id.
   virtual int number(RowId id) const = 0;
   /// returns the number of column with #ColId \p id.
   virtual int number(ColId id) const = 0;

   /// returns the #RowId of the \p i 'th inequality.
   virtual RowId rowId(int i) const = 0;
   /// returns the #ColId of the \p i 'th column.
   virtual ColId colId(int i) const = 0;

   /// test whether #LPSolver contains a row with #RowId \p id.
   virtual int has(RowId id) const
   {
      return number(id) >= 0;
   }
   /// test whether #LPSolver contains a column with #ColId \p id.
   virtual int has(ColId id) const
   {
      return number(id) >= 0;
   }

   /// gets the row ids.
   virtual void getRowIds(RowId ids[]) const = 0;

   /// gets the column ids.
   virtual void getColIds(ColId ids[]) const = 0;
   //@}

   /**@name Constructors / Destructors */
   //@{
   /// destructor.
   virtual ~LPSolver()
   {}
   //@}
};
} // namespace soplex
#endif // _LPSOLVER_H_

//-----------------------------------------------------------------------------
//Emacs Local Variables:
//Emacs mode:c++
//Emacs c-basic-offset:3
//Emacs tab-width:8
//Emacs indent-tabs-mode:nil
//Emacs End:
//-----------------------------------------------------------------------------
