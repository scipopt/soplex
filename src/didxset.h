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
#pragma ident "@(#) $Id: didxset.h,v 1.6 2001/12/28 14:55:12 bzfkocht Exp $"

/**@file  didxset.h
 * @brief Dymnamic index set.
 */
#ifndef _DIDXSET_H_
#define _DIDXSET_H_

#include <assert.h>

#include "idxset.h"

namespace soplex
{

/**@brief   Dynamic index set.
   @ingroup Elementary

   Class #DIdxSet provides dynamic #IdxSet in the sense, that no
   restrictions are posed on the use of methods #add(). However, method
   #indexMem() has been moved to the private members. This is because
   #DIdxSet adds its own memory management to class #IdxSet and the user must
   not interfer with it.
   
   Upon construction of an #DIdxSet, memory is allocated automatically. The
   memory consumption can be controlled with methods #max() and #setMax().
   Finally, the destructor will release all allocated memory.
*/
class DIdxSet : public IdxSet
{
private:
   int*& indexMem();  ///< points to the allocated memory

   /// no assignment operator.
   /** The assignment operator for IdxSet needed in the public implementation
    *  below is not implemented. So it won't work anyway.
    * @todo This needs cleanup.
    */
   DIdxSet& operator=(const IdxSet& sv);
   /// no copy constructor implemented yet.
   DIdxSet(const DIdxSet& old);

public:
#if 0
   /// assignment operator
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
#endif 

   /// adds \p n uninitialized indices.
   void add(int n)
   {
      if (max() - size() < n)
         setMax(size() + n);
      IdxSet::add(n);
   }

   /// adds all indices from \p sv.
   void add(const IdxSet& sv)
   {
      int n = sv.size();
      if (max() - size() < n)
         setMax(size() + n);
      IdxSet::add(sv);
   }

   /// adds \p n indices from \p i.
   void add(int n, const int *i)
   {
      if (max() - size() < n)
         setMax(size() + n);
      IdxSet::add(n, i);
   }

   /// adds index \p i to the index set
   void addIdx(int i)
   {
      if (max() <= size())
         setMax(size() + 1);
      IdxSet::addIdx(i);
   }

   /// sets the maximum number of indices.
   /** This methods resets the memory consumption of the #DIdxSet to
    *  \p newmax. However, if \p newmax < #size(), it is reset to #size()
    *  only.
    */
   void setMax(int newmax = 1);

   /// copy constructor from #IdxSet.
   explicit DIdxSet(const IdxSet& old);

   /// default constructor. \p n gives the initial size of the index space.
   explicit DIdxSet(int n = 8);

   /// destructor.
   ~DIdxSet();
};

} // namespace soplex
#endif // _DIDXSET_H_

//-----------------------------------------------------------------------------
//Emacs Local Variables:
//Emacs mode:c++
//Emacs c-basic-offset:3
//Emacs tab-width:8
//Emacs indent-tabs-mode:nil
//Emacs End:
//-----------------------------------------------------------------------------
