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
#pragma ident "@(#) $Id: spxlp.h,v 1.2 2001/11/06 23:31:04 bzfkocht Exp $"

#ifndef _SPXLP_H_
#define _SPXLP_H_

/*
    \Section{Imports}
 */
#include <assert.h>
#include <iostream>

#include "dvector.h"
#include "svset.h"
#include "dataarray.h"
#include "lprow.h"
#include "lpcol.h"
#include "lprowset.h"
#include "lpcolset.h"
#include "nameset.h"

namespace soplex
{

/*
    \Section{Class Declaration}
 */
class SPxLP_Id;

/** Ids for LP columns.
    This class should be called #SPxLP::SPxColId#. However, since AT\&T cfront
    does not support subclasses as template arguments, we had to perform an own
    name mangling.
 
    Class #SPxLP_SPxRowId# provides #DataSet::Key#s for the column indices of an
    #SPxLP#.
 */
class SPxLP_SPxColId : public SVSet_Key
{
   friend class SPxLP;
   SPxLP_SPxColId(const SVSet::Key& key)
      : SVSet_Key(key)
   { }

public:
   ///
   SPxLP_SPxColId()
   { }

   ///

   SPxLP_SPxColId(const SPxLP_Id& key);
};

/** Ids for LP rows.
    This class should be called #SPxLP::SPxRowId#. However, since AT\&T cfront
    does not support subclasses as template arguments, we had to perform an own
    name mangling.
 
    Class #SPxLP_SPxColId# provides #DataSet::Key#s for the row indices of an
    #SPxLP#.
 */
class SPxLP_SPxRowId : public SVSet_Key
{
   friend class SPxLP;
   SPxLP_SPxRowId(const SVSet::Key& key)
      : SVSet_Key(key)
   { }

public:
   ///
   SPxLP_SPxRowId()
   { }

   ///

   SPxLP_SPxRowId(const SPxLP_Id& key);
};


/*
    \SubSubSection{SPxLP\_Ids}
    Ids use the 2 msb of an #int# to store the id's #Type#,
    #SPX_MAX_PE_LOG# for storing the PE number
    and the remaining for the index.
*/
/** Generic Ids for LP rows or columns.
    This class should be called #SPxLP::Id#. However, since AT\&T cfront does
    not support subclasses as template arguments, we had to perform an own name
    mangling.
 
    Both, #SPxLP_SPxColId#s and #SPxLP_SPxRowId#s may be treated uniformly as
    #SPxLP_Id#s:
 
    Rows and columns are numbered from 0 to #num()# and 0 to #dim()#
    respectively.  These numbers may be used to select individual rows or
    columns. However, these numbers may change if other rows or columns are
    added or removed.
 
    Further, each row or column of the problem matrix is assigned a
    #SPxLP_SPxRowId# or #SPxLP_SPxColId#, respectively. They are be used to select
    individual rows or columns just like numbers. In contrast to row and column
    numbers, ids remain unchanges for the time a row or column belongs to a
    #SPxLP#, no matter what other rows or columns are added to it or removed
    from it.
 */
class SPxLP_Id : public SVSet_Key
{
public:
   ///
   enum Type
   {
      /// row identifier.
      ROWID = -1,
      /// none.
      NONE = 0,
      /// column identifier.
      COLID = 1
   };

   ///
   Type type() const
   {
      return info ? (info < 0 ? ROWID : COLID) : NONE;
   }

   ///
   int isValid() const
   {
      return info != 0;
   }

   ///
   void inValidate()
   {
      info = 0;
   }

   ///
   int isSPxRowId() const
   {
      return info < 0;
   }

   ///
   int isSPxColId() const
   {
      return info > 0;
   }


   ///
   int operator==(const SPxLP_Id& id)
   {
      return (*(int*)this == *(int*)&id);
   }

   ///
   int operator!=(const SPxLP_Id& id)
   {
      return (*(int*)this != *(int*)&id);
   }

   ///
   SPxLP_Id& operator=(const SPxLP_Id& id)
   {
      *((int*)this) = *((int*)(&id));
      return *this;
   }

   ///
   SPxLP_Id& operator=(const SPxLP_SPxColId cid)
   {
      *(int*)this = *(int*) & cid;
      info = COLID * (cid.info + 1);
      return *this;
   }

   ///
   SPxLP_Id& operator=(const SPxLP_SPxRowId rid)
   {
      *(int*)this = *(int*) & rid;
      info = ROWID * (rid.info + 1);
      return *this;
   }



   ///
   SPxLP_Id()
   {
      info = NONE;
      idx = -1;
   }

   ///
   SPxLP_Id(const SPxLP_SPxColId& cid)
   {
      info = COLID * (cid.info + 1);
      idx = cid.idx;
   }

   ///
   SPxLP_Id(const SPxLP_SPxRowId& rid)
   {
      info = ROWID * (rid.info + 1);
      idx = rid.idx;
   }
};

inline SPxLP_SPxRowId::SPxLP_SPxRowId(const SPxLP_Id& key)
   : SVSet_Key(key)
{
   assert(!key.isSPxColId());
   info = info * SPxLP_Id::ROWID - 1;
}

inline SPxLP_SPxColId::SPxLP_SPxColId(const SPxLP_Id& key)
   : SVSet_Key(key)
{
   assert(!key.isSPxRowId());
   info = info * SPxLP_Id::COLID - 1;
}

class SoPlex;

//@ -----------------------------------------------------------------------------
/** saving LPs in a form suitable for SoPlex.
    Class #SPxLP# provides the data structures required for saving a linear program
    in the form
    \[
    \begin{array}{rl}
        \hbox{max}      & c^T x         \\
        \hbox{s.t.}     & l_r \le Ax \le u_r    \\
                        & l_c \le x \le u_c
    \end{array}
    \]
    suitable for solving with #SoPlex#. This includes:
    \begin{itemize}
    \item       SVSets for both, columns and rows
    \item       objective Vector
    \item       upper and lower bound Vectors for variables ($l_c$ and $u_c$)
    \item       upper and lower bound Vectors for inequalities ($l_r$ and $u_r$)
    \end{itemize}
 
    Note, that the optimization sense is not saved directly. Instead, the
    objective function are multiplied by -1 to transform the LP to our standard
    form maximizing the objective function. However, the sense of the loaded LP
    can be retreived with method #spxSense()#.
 
    Further, equality constraints are modelled by $l_r = u_r$. Analogously, fixed
    variables have $l_c = u_c$.
*/
class SPxLP : protected LPRowSet, protected LPColSet
{
   friend class SPxBasis;
   friend class SPxScale;
   friend int getmarsz (SoPlex*);
   friend int getmartz (SoPlex*);
   friend int SPxLP__readLine(
      std::istream& is,
      char*& f1,
      char*& f2,
      char*& f3,
      char*& f4,
      char*& f5,
      char*& f6
  );

   static int readLine(
      std::istream& is,
      char*& f1,
      char*& f2,
      char*& f3,
      char*& f4,
      char*& f5,
      char*& f6
  );

   SVector& colVector(int i)
   {
      return LPColSet::colVector(i);
   }
   SVector& rowVector(int i)
   {
      return LPRowSet::rowVector(i);
   }

protected:
   const LPRowSet* lprowset() const
   {
      LPRowSet* svs = (LPRowSet*)this;
      return svs;
   }
   const LPColSet* lpcolset() const
   {
      LPColSet* svs = (LPColSet*)this;
      return svs;
   }

   SVSet* rowset()
   {
      SVSet* svs = (SVSet*)(LPRowSet*)this;
      return svs;
   }
   SVSet* colset()
   {
      SVSet* svs = (SVSet*)(LPColSet*)this;
      return svs;
   }

   /*
       \SubSection{Data structures and layout}
       #SPxLP#s are saved as an #SVSet# for both, the columns and rows. Note that this
       is redundant but eases the access.
   */
public:
   /**@name Datatypes */
   //@{
   ///
   typedef SPxLP_SPxRowId SPxRowId;
   ///
   typedef SPxLP_SPxColId SPxColId;
   /** Generic Identifier for LP rows and columns.
       All #Id# classes should better be local classes of #SPxLP#. However,
       AT\&T's cfront compiler won't allow to do do so. :-(
    */
   typedef SPxLP_Id Id;

   /// Optimization SPxSense.
   enum SPxSense
   {
      MAXIMIZE = 1,
      MINIMIZE = -1
   };
   //@}

private:
   SPxSense thesense;

public:
   ///
   static const double infinity;               //@Memo: value used as $\infty$.

   ///
   int nRows() const
   {
      return LPRowSet::num();
   }

   ///
   int nCols() const
   {
      return LPColSet::num();
   }

   /// get #i#-th row.
   void getRow(int i, LPRow& row) const;

   /// get #id#-th row.
   void getRow(SPxRowId id, LPRow& row) const
   {
      getRow(number(id), row);
   }

   /// get rows #start# .. #end#.
   void getRows(int start, int end, LPRowSet& set) const;

   ///
   const SVector& rowVector(int i) const
   {
      return LPRowSet::rowVector(i);
   }
   ///
   const SVector& rowVector(SPxRowId& id) const
   {
      return LPRowSet::rowVector(id);
   }

   ///
   const Vector& rhs() const
   {
      return LPRowSet::rhs();
   }
   ///
   double rhs(int i) const
   {
      return LPRowSet::rhs(i);
   }
   ///
   double& rhs(int i)
   {
      return LPRowSet::rhs(i);
   }
   ///
   double rhs(SPxRowId& id) const
   {
      return LPRowSet::rhs(id);
   }

   ///
   const Vector& lhs() const
   {
      return LPRowSet::lhs();
   }
   ///
   double lhs(int i) const
   {
      return LPRowSet::lhs(i);
   }
   ///
   double& lhs(int i)
   {
      return LPRowSet::lhs(i);
   }
   ///
   double lhs(SPxRowId& id) const
   {
      return LPRowSet::lhs(id);
   }

   /// get #i#-th column.
   void getCol(int i, LPCol& column) const;

   /// get #id#-th column.
   void getCol(SPxColId id, LPCol& col) const
   {
      getCol(number(id), col);
   }

   /// get columns #start# .. #end#.
   void getCols(int start, int end, LPColSet& set) const;

   ///
   const SVector& colVector(int i) const
   {
      return LPColSet::colVector(i);
   }
   ///
   const SVector& colVector(SPxColId& id) const
   {
      return LPColSet::colVector(id);
   }

   ///
   void getObj(Vector& obj) const;

   ///
   double obj(int i) const
   {
      return spxSense() * maxObj(i);
   }
   ///
   double obj(SPxColId& id) const
   {
      return spxSense() * maxObj(id);
   }

   /** Objective vector for maximization problem.
       Methods #maxObj()# return the objective vector or its elements, after
       transformation to a maximization problem. Since this is how #SPxLP#
       internally stores any LP these methods are generally faster. The
       following condition holds: #obj() = spxSense() * maxObj()#.
    */
   const Vector& maxObj() const
   {
      return LPColSet::obj();
   }
   ///
   double maxObj(int i) const
   {
      return LPColSet::obj(i);
   }
   ///
   double& maxObj(int i)
   {
      return LPColSet::obj(i);
   }
   ///
   double maxObj(SPxColId& id) const
   {
      return LPColSet::obj(id);
   }

   ///
   const Vector& upper() const
   {
      return LPColSet::upper();
   }
   ///
   double upper(int i) const
   {
      return LPColSet::upper(i);
   }
   ///
   double& upper(int i)
   {
      return LPColSet::upper(i);
   }
   ///
   double upper(SPxColId& id) const
   {
      return LPColSet::upper(id);
   }

   ///
   const Vector& lower() const
   {
      return LPColSet::lower();
   }
   ///
   double lower(int i) const
   {
      return LPColSet::lower(i);
   }
   ///
   double& lower(int i)
   {
      return LPColSet::lower(i);
   }
   ///
   double lower(SPxColId& id) const
   {
      return LPColSet::lower(id);
   }

   ///
   SPxSense spxSense() const
   {
      return thesense;
   }

   ///
   int number(const SPxRowId& id) const
   {
      return LPRowSet::number(id);
   }

   ///
   int number(const SPxColId& id) const
   {
      return LPColSet::number(id);
   }

   ///
   int number(const Id& id) const
   {
      return (id.type() == Id::COLID)
             ? LPColSet::number(id)
          : LPRowSet::number(id);
   }

   ///
   SPxRowId rId(int n) const
   {
      return LPRowSet::key(n);
   }

   ///
   SPxColId cId(int n) const
   {
      return LPColSet::key(n);
   }

private:
   void doAddRow (const LPRow& row);
   void doAddRows(const LPRowSet& set);
   void doAddCol (const LPCol& col);
   void doAddCols(const LPColSet& set);

public:
   /**@name Extension */
   //@{
   ///
   void addRow(const LPRow& row)
   {
      doAddRow(row);
   }
   /// add #row# to #LPRowSet#.
   void addRow(SPxRowId& id, const LPRow& row)
   {
      addRow(row);
      id = rId(nRows() - 1);
   }

   ///
   void addRows(const LPRowSet& set)
   {
      doAddRows(set);
   }
   /// add all #LPRow#s of #set# to #LPRowSet#.
   void addRows(SPxRowId id[], const LPRowSet& set);

   ///
   void addCol(const LPCol& col)
   {
      doAddCol(col);
   }
   /// add #col# to #LPColSet#.
   void addCol(SPxColId& id, const LPCol& col)
   {
      addCol(col);
      id = cId(nCols() - 1);
   }

   ///
   void addCols(const LPColSet& set)
   {
      doAddCols(set);
   }
   /// add all #LPCol#s of #set# to #LPColSet#.
   void addCols(SPxColId id[], const LPColSet& set);

protected:
   /// called after the last #n# rows have just been added.
   virtual void addedRows(int n)
   {
      (void)n;
   }
   /// called after the last #n# columns have just been added.
   virtual void addedCols(int n)
   {
      (void)n;
   }

   void added2Set(SVSet& set, const SVSet& add, int n);
   //@}


public:
   /**@name Shrinking */
   //@{
   /// remove #i#-th row.
   void removeRow(int i)
   {
      doRemoveRow(i);
   }
   /// remove #id#-th row.
   void removeRow(SPxRowId id)
   {
      removeRow(number(id));
   }

   /// remove multiple columns.
   void removeCols(int perm[])
   {
      doRemoveCols(perm);
   }
   /** Remove multiple rows.
       This method removes all #LPRow#s from the #SPxLP# with an
       index #i# such that #perm[i] < 0#. Upon completion, #perm[i] >= 0#
       indicates the new index where the #i#-th #LPRow# has been moved to
       due to this removal. Note, that #perm# must point to an array of at
       least #nRows()# #int#s.
    */
   void removeRows(int perm[])
   {
      doRemoveRows(perm);
   }

   ///
   void removeRows(SPxRowId id[], int n, int perm[] = 0);
   /** Remove #n# #LPRow#s.
       Removing multiple rows with one method invocation is available
       two flavours. An array #perm# can be passed as third argument or
       not. If given, #perm# must be an array at least of size #nCols()#. It
       is used to return the permutations resulting from this removal:
       #perm[i] < 0# indicates, that the element to index #i# has been
       removed.  Otherwise, #perm[i]# is the new index of the element with
       index #i# before the removal.
    */
   void removeRows(int nums[], int n, int perm[] = 0);
   /// remove rows from #start# to #end# (including both).
   void removeRowRange(int start, int end, int perm[] = 0);

   /// remove #i#-th column.
   void removeCol(int i)
   {
      doRemoveCol(i);
   }
   /// remove #id#-th column.
   void removeCol(SPxColId id)
   {
      removeCol(number(id));
   }


   ///
   void removeCols(SPxColId id[], int n, int perm[] = 0);
   /// remove #n# #LPCol#s.
   void removeCols(int nums[], int n, int perm[] = 0);
   /// remove columns from #start# to #end# (including both).
   void removeColRange(int start, int end, int perm[] = 0);

   ///
   virtual void clear();


protected:
   //@Memo:        These  methods are use for implemementing the public remove methods
   virtual void doRemoveRow(int i);
   ///
   virtual void doRemoveCols(int perm[]);
   ///
   virtual void doRemoveRows(int perm[]);
   ///
   virtual void doRemoveCol(int i);
   //@}

public:
   /**@name IO */
   //@{
   /// read a file from #in#.
   virtual void read (std::istream& in, NameSet* rowNames = 0,
                       NameSet* colNames = 0, DIdxSet* intVars = 0);
   /// read a file in LP format from #in#.
   virtual void readLP (std::istream& in, NameSet* rowNames = 0,
                         NameSet* colNames = 0, DIdxSet* intVars = 0);
   /// read a file in MPS format from #in#.
   virtual void readMPS(std::istream& in, NameSet* rowNames = 0, NameSet* colNames = 0);
   ///
   friend std::istream& operator>>(std::istream& is, SPxLP& lp)
   {
      lp.read(is);
      return is;
   }
   ///
   friend std::ostream& operator<<(std::ostream& os, const SPxLP& lp);
   //@}


   /**@name Manipulation */
   //@{
   /// change objective vector.
   virtual void changeObj(const Vector& newObj);

   /// change #i#-th objective vector element.
   virtual void changeObj(int i, double newVal);

   /// change #id#-th objective vector element.
   virtual void changeObj(SPxColId id, double newVal)
   {
      changeObj(number(id), newVal);
   }

   /// change vector of lower bounds.
   virtual void changeLower(const Vector& newLower);

   /// change #i#-th lower bound.
   virtual void changeLower(int i, double newLower);

   /// change #id#-th lower bound.
   virtual void changeLower(SPxColId id, double newLower)
   {
      changeLower(number(id), newLower);
   }

   /// change vector of upper bounds.
   virtual void changeUpper(const Vector& newUpper);

   /// change #i#-th upper bound.
   virtual void changeUpper(int i, double newUpper);

   /// change #id#-th upper bound.
   virtual void changeUpper(SPxColId id, double newUpper)
   {
      changeUpper(number(id), newUpper);
   }

   ///
   virtual void changeBounds(const Vector& newLower, const Vector& newUpper);

   ///
   virtual void changeBounds(int i, double newLower, double newUpper);

   ///
   virtual void changeBounds(SPxColId id, double newLower, double newUpper)
   {
      changeBounds(number(id), newLower, newUpper);
   }


   /// change lhs vector for constraints.
   virtual void changeLhs(const Vector& newLhs);

   /// change #i#-th lhs value.
   virtual void changeLhs(int i, double newLhs);

   /// change #id#-th lhs value.
   virtual void changeLhs(SPxRowId id, double newLhs)
   {
      changeLhs(number(id), newLhs);
   }

   /// change rhs vector for constraints.
   virtual void changeRhs(const Vector& newRhs);

   /// change #i#-th rhs value.
   virtual void changeRhs(int i, double newRhs);

   /// change #id#-th rhs value.
   virtual void changeRhs(SPxRowId id, double newRhs)
   {
      changeRhs(number(id), newRhs);
   }

   ///
   virtual void changeRange(const Vector& newLhs, const Vector& newRhs);

   ///
   virtual void changeRange(int i, double newLhs, double newRhs);

   ///
   virtual void changeRange(SPxRowId id, double newLhs, double newRhs)
   {
      changeRange(number(id), newLhs, newRhs);
   }


   /// change #i#-th row of LP.
   virtual void changeRow(int i, const LPRow& newRow);

   /// change #id#-th row of LP.
   virtual void changeRow(SPxRowId id, const LPRow& newRow)
   {
      changeRow(number(id), newRow);
   }

   /// change #i#-th column of LP.
   virtual void changeCol(int i, const LPCol& newCol);

   /// change #id#-th column of LP.
   virtual void changeCol(SPxColId id, const LPCol& newCol)
   {
      changeCol(number(id), newCol);
   }

   /// change LP element (#i#, #j#).
   virtual void changeElement(int i, int j, double val);

   /// change LP element (#rid#, #cid#).
   virtual void changeElement(SPxRowId rid, SPxColId cid, double val)
   {
      changeElement(number(rid), number(cid), val);
   }

   /// change optimization sense to #sns#.
   virtual void changeSense(SPxSense sns)
   {
      if (sns != thesense)
         LPColSet::obj() *= -1;
      thesense = sns;
   }
   //@}


   /// consistency check.
   int isConsistent() const;

   ///
   SPxLP()
   {
      clear();
   }
   virtual ~SPxLP()
   {}

};


} // namespace soplex
#endif  // _SPXLP_H_

//-----------------------------------------------------------------------------
//Emacs Local Variables:
//Emacs mode:c++
//Emacs c-basic-offset:3
//Emacs tab-width:8
//Emacs indent-tabs-mode:nil
//Emacs End:
//-----------------------------------------------------------------------------
