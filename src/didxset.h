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
#pragma ident "@(#) $Id: didxset.h,v 1.1 2001/11/06 16:18:32 bzfkocht Exp $"


#ifndef _DIDXSET_H_
#define _DIDXSET_H_

/*      \Section{DIdxSet}
*/

#include <assert.h>
#include <stdlib.h>

#include "idxset.h"

namespace soplex
{

/** dynamic index set.
    Class #DIdxSet# provids dynamic \Ref{IdxSet} in the sense, that no
    restrictions are posed on the use of methods #add()#. However, method
    #indexMem()# has been moved to the #private# members. This is because
    #DIdxSet# adds it own memory managment to class #IdxSet# and the user must
    not interfer with it.
 
    Upon construction of an #DIdxSet#, memory is allocated automatically. The
    memory consumption can be controlled with methods #max()# and #setMax()#.
    Finally, the destructor will release all allocated memory.
*/
class DIdxSet : public IdxSet
{
private:
   int*& indexMem();

public:
   /*
       Die Implementierung st\"utzt sich voll auf die Datenstrukturen der
       Basisklasse #IdxSet#. Lediglich die Methoden, bei denen Dynamizit\"at
       hinzugef\"ugt wird, m\"ussen \"uberlagert werden.
   */
   /// assignment operator.
   DIdxSet& operator=(const IdxSet& sv)
   {
      int n = sv.size();
      if (max() - size() < n)
         setMax(size() + n);
      IdxSet::operator=(sv);
      return *this;
   }

   /// assignment operator.
   DIdxSet& operator=(const DIdxSet& sv)
   {
      return operator=(IdxSet(sv));
   }

   /// add #n# uninitialized indices.
   void add(int n)
   {
      if (max() - size() < n)
         setMax(size() + n);
      IdxSet::add(n);
   }

   /// add all indices from #sv#.
   void add(const IdxSet& sv)
   {
      int n = sv.size();
      if (max() - size() < n)
         setMax(size() + n);
      IdxSet::add(sv);
   }

   /// add #n# indices from #i#.
   void add(int n, const int *i)
   {
      if (max() - size() < n)
         setMax(size() + n);
      IdxSet::add(n, i);
   }

   /** set maximum number of indices.
       This methods resets the memory consumption of the #DIdxSet# to
       #newmax#. However, if #newmax < size()#, it is reset to #size()#
       only.
    */
   void setMax(int newmax = 1);

   void addIdx(int i)
   {
      if (max() <= size())
         setMax(size() + 1);
      IdxSet::addIdx(i);
   }
   ///
   DIdxSet(const IdxSet& old);
   ///
   DIdxSet(const DIdxSet& old);
   ///
   DIdxSet(int n = 8);
   ///
   ~DIdxSet();
};

} // namespace soplex
#endif // _DIDXSET_H_

//-----------------------------------------------------------------------------
//Emacs Local Variables:
//Emacs c-basic-offset:3
//Emacs tab-width:8
//Emacs indent-tabs-mode:nil
//Emacs End:
//-----------------------------------------------------------------------------
