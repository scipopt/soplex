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
#pragma ident "@(#) $Id: nameset.h,v 1.5 2001/11/12 17:10:02 bzfkocht Exp $"

#ifndef _NAMESET_H_
#define _NAMESET_H_


//@ ----------------------------------------------------------------------------
/*      \Section{Imports}
    Import required system include files
 */
#include <assert.h>
#include <sys/types.h>
#include <stdlib.h>


/*  and class header files
 */

#include "dataset.h"
#include "datahashtable.h"
#include "islist.h"
#include "datakey.h"

namespace soplex
{




//@ ----------------------------------------------------------------------------
/*      \Section{Class Declarations}
    The first two classes should be local classes of NameSet.
 
    Class #NameSet_Name# provide the handles (i.e. #char*#s) of names in a
    #NameSet#.
 */
class NameSet_Name
{
protected:
   friend class NameSet;
   static char deflt;                   // default zero string

public:
   const char *name;

   friend int operator==(const NameSet_Name& n1, const NameSet_Name& n2)
   {
      return (strcmp (n1.name, n2.name) == 0);
   }

   friend std::ostream& operator<<(std::ostream& out, const NameSet_Name& n)
   {
      return out << n.name;
   }

   friend int hashFunction (const NameSet_Name&);

   int isConsistent () const
   {
      return (name != 0);
   }
   NameSet_Name (const char* str)
   {
      name = str;
   }
   NameSet_Name (const NameSet_Name& str)
   {
      name = str.name;
   }
   NameSet_Name ()
   {
      name = &deflt;
   }
};

/*
    Since, #typedef IsElement<Name> CharPtr# doesn't seem to work, derive
    CharPtr explicitly:
 */
class NameSet_CharPtr : public NameSet_Name
{
protected:
   NameSet_CharPtr *the_next;
public:
   NameSet_CharPtr*& next()
   {
      return the_next;
   }
   NameSet_CharPtr*const& next() const
   {
      return the_next;
   }

   NameSet_CharPtr(const NameSet_CharPtr& org)
      : NameSet_Name(org), the_next(0)
   { }

   NameSet_CharPtr()
      : NameSet_Name(), the_next(0)
   {}

};

/** set of strings.
    Class #NameSet# implements a symbol or name table. It allows to store or
    remove names (i.e. #char*#), but does not provide means for manipulating
    stored names.
    
    Names in a #NameSet# may be accessed via numbers form 0 through #num()-1#
    and #Key#s. See \Ref{DataSet} for a description of these concepts.
    
    At a time a #NameSet# can hold a maximum of #max()# entries. This can be
    reset with method #reMax()#. If more than #max()# names are added to a
    #NameSet#, it adjusts itself automatically to the required size.  This
    implies, that references to names within a #NameSet# may become invalid if
    the #NameSet# is expanded.
 
    All names (i.e. the actual #char# strings) in a #NameSet# are stored in one
    continuous memory block of size #memMax()#. At one time #memSize()# bytes of
    it are used for actually saving names; the remaining memory is free to hold
    additional names. #memRemax()# can be used to reset #memMax()# but not lower
    than to #memSize()#. Method #memPack()# performs a garbage collection to
    gain free memory resulting from removed names.
 */
class NameSet
{
private:
   /** Identifier for set entries.
       Every name in a #NameSet# is assigned a #Key# by which it can be
       accessed (see #NameSet::operator[]#). See \Ref{DataSet::Key} for a more
       detailed description of the concept of Keys.
   */
   typedef DataKey Key; 

protected:
   /*
       Ein #NameSet# ist als #DataSet# von #CharPtr# implementiert. #CharPtr#s sind
       i.w.~(einfach verkettete) Zeiger auf Namen (#char*#), die alle in das Feld #mem#
       zeigen. Genutzte #CharPtr#s werden in der einfach verketteten Liste #list# so
       gehalten, da\ss\ ihre Zeiger aufsteigend sortiert sind. Dies erm\"oglicht eine
       effiziente Implementierung der Garbage Collection.
    */
   IsList < NameSet_CharPtr > list;           // sorted name list
   DataSet < NameSet_CharPtr > set;            // name set
   char* mem;            // string memory
   int memmax;         // #memMax()#
   int memused;        // #memSize()#
   DataHashTable < NameSet_Name, Key > hashtab;        // hashtable for names

public:
   /**@name Inquiry */
   //@{
   /// return #num#-th name of #NameSet#.
   const char* operator[](int pnum) const
   {
      return set[pnum].name;
   }
   /// return name to #key# of #NameSet#.
   const char* operator[](Key pkey) const
   {
      return set[pkey].name;
   }

   /// return nr. of names in #NameSet#.
   int num() const
   {
      return set.num();
   }
   /// return maximum nr. of names that fit into #NameSet#.
   int max() const
   {
      return set.max();
   }
   /// return maximum #Key::idx# used in #NameSet#.
   int size() const
   {
      return set.size();
   }
   /// maximum length of string memory.
   int memMax() const
   {
      return memmax;
   }
   /// used length of string memory.
   int memSize() const
   {
      return memused;
   }

   /// return #Key# to #num#-th name in #NameSet#.
   Key key(int pnum) const
   {
      return set.key(pnum);
   }
   /// return #Key# to name #str# in #NameSet#.
   Key key(const char* str) const
   {
      const NameSet_Name nam(str);
      return (*hashtab.get(nam));
   }

   /// return number of name with #key# in #NameSet#.
   int number(Key pkey) const
   {
      return set.number(pkey);
   }
   /// return number of name #str# in #NameSet#.
   int number(const char *str) const
   {
      const NameSet_Name nam(str);
      if (hashtab.has(nam))
         return number(*hashtab.get(nam));
      else
         return -1;
   }

   /// does #NameSet# have name with number #num#?.
   int has(int pnum) const
   {
      return set.has(pnum);
   }
   /// does #NameSet# have name #str#?.
   int has(const char* str) const
   {
      const NameSet_Name nam(str);
      return hashtab.has(nam);
   }
   /// does #NameSet# have name with #Key key#?.
   int has(Key pkey) const
   {
      return set.has(pkey);
   }
   //@}

   /**@name Extension */
   //@{
   ///
   void add(const char* str);
   /// add name #str# to #NameSet#.
   void add(Key& key, const char* str);

   ///
   void add(const NameSet& set);
   /// add all names in #set# to #NameSet#.
   void add(Key key[], const NameSet& set);
   //@}


   /**@name Shrinking */
   //@{
   /// remove name to #key# from #NameSet#.
   void remove(Key key);
   /// remove #i#-th name from #NameSet#.
   void remove(int pnum)
   {
      remove(key(pnum));
   }
   /// remove name #str# from #NameSet#.
   void remove(const char* str);
   /// remove #n# names with #Key#s #keys# from #NameSet#.
   void remove(Key keys[], int n);
   /// remove #n# names with numbers #nums# from #NameSet#.
   void remove(int nums[], int n);

   /// remove all names from #NameSet#.
   void clear();
   //@}


   /**@name Memory Control */
   //@{
   /// reset #max()# to #newmas#.
   void reMax(int newmax = 0);
   /// reset #memMax()# to #newmas#.
   void memRemax(int newmax = 0);
   /// garbage collection.
   void memPack();
   //@}


   /**@name Control Parameters */
   //@{
   /** memory extension factor for entries.
       When more than #max()# names are added to a #NameSet#, it is
       automatically resized to fit the additional names. Parameter
       #factor# is the factor by which the element memory is extended to do
       so.
    */
   double factor;

   /** memory extension factor for names.1
       When the names added to a #NameSet# do no longer fit into the name
       memory (i.e. the memory for saving the strings) it is automatically
       resized to fit the additional names. Parameter #memFactor# is the
       factor by which this memory is extended to do so.
    */
   double memFactor;
   //@}

   /**@name Miscellaneous */
   //@{
   /// consistency check.
   int isConsistent() const;
   ///
   NameSet& operator=(const NameSet& rhs);
   ///
   NameSet(const NameSet& old);

   /** Constructor.
       @param      max     start value for \Ref{max}()
       @param      mmax    start value for \Ref{memMax}()
       @param      fac     start value for \Ref{factor}
       @param      memFac  start value for \Ref{memFactor}
    */
   NameSet(int max = 8,
           int mmax = -1,
           double fac = 2,
           double memFac = 2);

   ///
   ~NameSet();
   //@}
}
;


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
