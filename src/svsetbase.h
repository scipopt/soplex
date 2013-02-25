/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the class library                   */
/*       SoPlex --- the Sequential object-oriented simPlex.                  */
/*                                                                           */
/*    Copyright (C) 1996-2012 Konrad-Zuse-Zentrum                            */
/*                            fuer Informationstechnik Berlin                */
/*                                                                           */
/*  SoPlex is distributed under the terms of the ZIB Academic Licence.       */
/*                                                                           */
/*  You should have received a copy of the ZIB Academic License              */
/*  along with SoPlex; see the file COPYING. If not email to soplex@zib.de.  */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/**@file  svsetbase.h
 * @brief Set of sparse vectors.
 */

#ifndef _SVSETBASE_H_
#define _SVSETBASE_H_

#include <assert.h>

#include "spxdefines.h"
#include "svectorbase.h"
#include "classarray.h"
#include "dataset.h"
#include "datakey.h"
#include "idlist.h"

namespace soplex
{
/**@brief   Sparse vector set.
 * @ingroup Algebra
 *
 *   Class SVSetBase provides a set of sparse vectors SVectorBase. All SVectorBase%s in an SVSetBase share one big
 *   memory block for their nonzeros. This memory is reffered to as the \em nonzero \em memory. The SVectorBase%s
 *   themselves are saved in another memory block refered to as the \em vector \em memory. Both memory blocks will grow
 *   automatically if required, when adding more SVectorBase%s to the set or enlarging SVectorBase%s within the set. For
 *   controlling memory consumption, methods are provided to inquire and reset the size of the memory blocks used for
 *   vectors and nonzeros.
 *
 *   SVectorBase%s in an SVSetBase are numbered from 0 thru num()-1. They can be accessed using the index
 *   operator[](). When removing SVectorBase%s of a SVSetBase the remaining ones will be renumbered. However, all
 *   SVectorBase with a smaller number than the lowest number of the removed SVectorBase%s will remain unchanged.
 *
 *   For providing a uniform access to SVectorBase%s in a %set even if others are removed or added, SVSetBase assigns a
 *   DataKey to each SVectorBase in the %set. Such a DataKey remains unchanged as long as the corresponding SVectorBase
 *   is in the SVSetBase, no matter what other SVectorBase%s are added to or removed from the SVSetBase. Methods are
 *   provided for getting the DataKey to a SVectorBase or its number and vice versa.  Further, each add() method for
 *   enlarging an SVSetBase is provided with two signatures. One of them returns the DataKey%s assigned to the
 *   SVectorBase%s added to the SVSetBase.
 */
template < class R >
class SVSetBase : protected ClassArray < Element<R> >
{
private:

   typedef ClassArray < Element<R> > SVSetBaseArray;

   /**@class DLPSV
    * @brief SVectorBase with prev/next pointers
    * @todo  Check whether SVSetBase::DLPSV can be implemented as IdElement<SVectorBase>
    *
    *  The management of the SVectorBase%s is implemented by a DataSet<DLPSV>, the keys used externally are DataKey%s.
    *
    *  The management of nonzeros is done by a Real linked list IdList<DLPSV>, where the SVectorBase%s are kept in the
    *  order in which their indices occurr in the Array. The SVectorBase%s are kept without holes: If one is removed or
    *  moved to the end, the SVectorBase preceeding it obtains the space for all the nonzeros that previously belonged
    *  to the (re-)moved one.  However, the nonzeros in use are uneffected by this.
    */
   class DLPSV : public SVectorBase<R>
   {
   private:

      // ---------------------------------------------------------------------------------------------------------------
      /**@name Data */
      //@{

      DLPSV* thenext; ///< next SVectorBase
      DLPSV* theprev; ///< previous SVectorBase

      //@}

   public:

      // ---------------------------------------------------------------------------------------------------------------
      /**@name Construction / destruction */
      //@{

      /// Default constructor.
      DLPSV()
         : SVectorBase<R>()
      {}

      /// Copy constructor.
      DLPSV(const DLPSV& copy)
         : SVectorBase<R>(copy)
      {}

      //@}

      // ---------------------------------------------------------------------------------------------------------------
      /**@name Successor / predecessor */
      //@{

      /// Next SVectorBase.
      DLPSV*& next()
      {
         return thenext;
      }

      /// Next SVectorBase.
      DLPSV* const& next() const
      {
         return thenext;
      }

      /// Previous SVectorBase.
      DLPSV* const& prev() const
      {
         return theprev;
      }

      /// Previous SVectorBase.
      DLPSV*& prev()
      {
         return theprev;
      }

      //@}
   };

   // ------------------------------------------------------------------------------------------------------------------
   /**@name Data */
   //@{

   DataSet < DLPSV > set;  ///< %set of SVectorBase%s
   IdList < DLPSV > list;  ///< doubly linked list for non-zero management
   int possiblyUnusedMem;  ///< an estimate of the used memory due to xtends

   //@}

   // ------------------------------------------------------------------------------------------------------------------
   /**@name Control Parameters */
   //@{

   double factor;          ///< sparse vector memory enlargment factor

   //@}

   // ------------------------------------------------------------------------------------------------------------------
   /**@name Helpers */
   //@{

   /// Provides enough vector memory for \p n more SVectorBase%s.
   void ensurePSVec(int n)
   {
      if( num() + n > max() )
      {
         assert(factor > 1);

         reMax(int(factor*max()) + 8 + n);
      }
   }

   /// Provides enough nonzero memory for \p n more Elements%s.
   void ensureMem(int n)
   {
      if( memSize() + n > memMax() )
      {
         int newMax = int(SVSetBaseArray::memFactor * memMax());

         if( memSize() + n > newMax )
            newMax = memSize() + n;

         memRemax(newMax);
      }
   }

   /// Deleting a vector from the data array and the list.
   void deleteVec(DLPSV* ps)
   {
      /* delete last entries */
      if( ps == list.last() )
      {
         removeLast(ps->max());
      }
      /* merge space of predecessor with position which will be deleted, therefore we do not need to delete any memory
       * or do an expensive memory reallocation
       */
      else if( ps != list.first() )
      {
         SVectorBase<R>* prev = ps->prev();
         int sz = prev->size();

         prev->setMem(prev->max() + ps->max(), prev->mem());
         prev->set_size(sz);
      }
      /* delete the front entries of the first list entry and correct the memory pointers in the vectors; we do this by
       * merging the first both vectors, move the entries from the second vector up front, and correcting the size
       */
      else
      {
         SVectorBase<R>* next = ps->next();
         int sz = next->size();
         int bothmax = next->max() + ps->max();
         int offset = 0;

         /* the first element does not need to start at the beginning of the data array, because if the first vector is
          * extended, see xtend(), it is shifted to the end leaving unused memory behind
          */
         while( &(this->SVSetBaseArray::operator[](offset)) != ps->mem() )
         {
            ++offset;
            assert(offset < SVSetBaseArray::size());
         }

         /* move all entries of the second vector to the front */
         for( int j = 0; j <= sz; ++j )
         {
            this->SVSetBaseArray::operator[](j) = next->mem()[j];
         }

         /* correct the data memmory pointer and the maximal space */
         next->setMem(bothmax, ps->mem());

         /* correct size */
         next->set_size(sz);
      }

      list.remove(ps);
   }

   //@}

public:

   // ------------------------------------------------------------------------------------------------------------------
   /**@name Extension */
   //@{

   /// Adds \p svec to the %set.
   /** This includes copying its nonzeros to the sets nonzero memory and creating an additional SVectorBase entry in
    *  vector memory. If neccessary, the memory blocks are enlarged appropriately.
    */
   void add(const SVectorBase<R>& svec)
   {
      ensurePSVec(1);
      SVectorBase<R>* new_svec = create(svec.size());
      *new_svec = svec;
   }

   /// Adds \p svec to SVSetBase.
   /** Adds SVectorBase \p svec to the %set. This includes copying its nonzeros to the sets nonzero memory and creating
    *  an additional SVectorBase entry in vector memory. If neccessary, the memory blocks are enlarged appropriately.
    *
    *  @return \p nkey contains the DataKey, that the SVSetBase has assosicated to the new SVectorBase.
    */
   void add(DataKey& nkey, const SVectorBase<R>& svec)
   {
      ensurePSVec(1);
      SVectorBase<R>* new_svec = create(nkey, svec.size());
      *new_svec = svec;
   }

   /// Adds all \p n SVectorBase%s in the array \p svec to the %set.
   /** @pre \p svec must be not larger than \p n.
    */
   void add(const SVectorBase<R> svec[], int n)
   {
      assert(n >= 0);

      int i;
      int len;

      for( i = len = 0; i < n; ++i )
         len += svec[i].size();

      ensurePSVec(n);
      ensureMem(len);

      for( i = 0; i < n; ++i )
         *create(svec[i].size()) = svec[i];
   }

   /// Adds n SVectorBase%s to SVSetBase.
   /** Adds all \p n SVectorBase%s in the array \p svec to the %set.
    *
    *  @return \p nkey contains the DataKey%s, that the SVSetBase has assosicated to the new SVectorBase%s.
    *
    *  @pre \p nkey must be large enough to fit \p n DataKey%s.
    */
   void add(DataKey nkey[], const SVectorBase<R> svec[], int n)
   {
      add(svec, n);

      for( int i = num() - 1; --n; --i )
         nkey[n] = key(i);
   }

   /// Adds all SVectorBase%s in \p pset to an SVSetBase.
   void add(const SVSetBase<R>& pset)
   {
      int i;
      int n;
      int len;

      n = pset.num();
      for( i = len = 0; i < n; ++i )
         len += pset[i].size();

      ensurePSVec(n);
      ensureMem(len);

      for( i = 0; i < n; ++i )
         *create(pset[i].size()) = pset[i];
   }

   /// Adds all SVectorBase%s of \p pset to SVSetBase.
   /** Adds all \p n SVectorBase%s in the \p pset to an SVSetBase.
    *
    * @return \p nkey contains the DataKey%s, that the SVSetBase has assosicated to the new SVectorBase%s.
    *
    * @pre \p nkey must be large enough to fit \p pset.num() DataKey%s.
    */
   void add(DataKey nkey[], const SVSetBase<R>& pset)
   {
      add(pset);

      int i = num();
      int n = pset.num();

      while( n > 0 )
         nkey[--n] = key(--i);
   }

   /// Creates new SVectorBase in %set.
   /** The new SVectorBase will be ready to fit at least \p idxmax nonzeros.
    */
   SVectorBase<R>* create(int idxmax = -1)
   {
      DLPSV* ps;

      if( list.last() )
      {
         ps = list.last();
         removeLast(ps->max() - ps->size());
         ps->set_max( ps->size() );
      }

      if( idxmax < 0 )
      {
         ensureMem(2);
         idxmax = memMax() - memSize();
      }
      else
         ensureMem(idxmax);

      ensurePSVec(1);

      assert(idxmax > 0);
      assert(memMax() >= memSize() + idxmax);

      ps = set.create();
      list.append(ps);

      // resize the data array
      SVSetBaseArray::reSize(memSize() + idxmax);

      ps->setMem(idxmax, &SVSetBaseArray::last() - idxmax + 1);

      return ps;
   }

   /// Creates new SVectorBase in %set.
   /** The new SVectorBase will be ready to fit at least \p idxmax nonzeros.
    *
    * @return \p nkey contains the DataKey associated to the new SVectorBase.
    */
   SVectorBase<R>* create(DataKey& nkey, int idxmax = -1)
   {
      SVectorBase<R>* ps = create(idxmax);

      nkey = key(num() - 1);

      return ps;
   }

   /// Extends \p svec to fit \p newmax nonzeros.
   /** @pre \p svec must be an SVectorBase of the SVSetBase.
    */
   void xtend(SVectorBase<R>& svec, int newmax)
   {
      if( svec.max() < newmax )
      {
         if( possiblyUnusedMem * SVSetBaseArray::memFactor > memSize() )
            memPack();

         assert(has(&svec));

         DLPSV* ps = static_cast<DLPSV*>(&svec);

         if( ps == list.last() )
         {
            int sz = ps->size();
            ensureMem(newmax - ps->max());
            insert(memSize(), newmax - ps->max());
            ps->setMem(newmax, ps->mem());
            ps->set_size(sz);
         }
         else
         {
            ensureMem(newmax);
            SVectorBase<R> newps(newmax, &SVSetBaseArray::last() + 1);
            int sz = ps->size();
            insert(memSize(), newmax);
            newps = svec;

            if( ps != list.first() )
            {
               SVectorBase<R>* prev = ps->prev();
               int prevsz = prev->size();
               prev->setMem(prev->max() + ps->max(), prev->mem());
               prev->set_size(prevsz);

               possiblyUnusedMem += ps->max();
            }

            list.remove(ps);
            list.append(ps);

            ps->setMem(newmax, newps.mem());
            ps->set_size(sz);
         }
      }
   }

   /// Adds nonzero (\p idx, \p val) to \p svec of this SVSetBase.
   /** Adds one nonzero (\p idx, \p val) to SVectorBase \p svec in the SVSetBase.  If \p svec is not large enough to
    *  hold the additional nonzero, it will be automatically enlarged within the %set.
    *
    *  @pre \p svec must be an SVectorBase of the SVSetBase.
    */
   void add2(SVectorBase<R> &svec, int idx, R val)
   {
      xtend(svec, svec.size() + 1);
      svec.add(idx, val);
   }

   /// Adds \p n nonzeros to \p svec of this SVSetBase.
   /** Adds \p n nonzeros to SVectorBase \p svec in the SVSetBase. If \p svec is not large enough to hold the additional
    *  nonzeros, it will be automatically enlarged within the %set.
    *
    * @pre \p svec must be an SVectorBase of the SVSetBase.
    */
   void add2(SVectorBase<R> &svec, int n, const int idx[], const R val[])
   {
      xtend(svec, svec.size() + n);
      svec.add(n, idx, val);
   }

   //@}

   // ------------------------------------------------------------------------------------------------------------------
   /**@name Shrinking */
   //@{

   /// Removes the vector with key \p removekey from the %set.
   /** @pre \p removekey must be a key from SVSetBase
    */
   void remove(const DataKey& removekey)
   {
      deleteVec(&set[removekey]);
      set.remove(removekey);
   }

   /// Removes the vector with number \p removenum from the %set.
   /** @pre \p removenum must be a valid vector number from SVSetBase
    */
   void remove(int removenum)
   {
      remove(key(removenum));
   }

   /// Removes one SVectorBase from %set.
   /** @pre \p svec must be from SVSetBase
    */
   void remove(const SVectorBase<R> *svec)
   {
      remove(key(svec));
   }

   /// Removes multiple elements.
   /** Removes all SVectorBase%s for the SVSetBase with an index \c i such that \p perm[i] < 0. Upon completion, \p
    *  perm[i] >= 0 indicates the new index where the \c i 'th SVectorBase has been moved to due to this removal.
    *
    *  @pre \p perm must point to an array of at least num() integers.
    */
   void remove(int perm[])
   {
      int j = num();

      /* due to performance reasons we use a backwards loop to delete entries, because it could result instead of only
       * decreasing the number of elements j times in memmoving the whole array j times
       */
      for( int i = j - 1; i >= 0; --i )
      {
         if( perm[i] < 0 )
         {
            deleteVec(&set[i]);
         }
      }

      set.remove(perm);
   }

   /// Removes \p n SVectorBase%s from %set.
   /** @pre \p keys must be at least of size \p n and valid keys
    */
   void remove(const DataKey keys[], int n)
   {
      DataArray < int > perm(num());
      remove(keys, n, perm.get_ptr());
   }

   /// Removes \p n SVectorBase%s from %set.
   /** @pre \p nums must be at least of size \p n and valid vector numbers
    */
   void remove(const int nums[], int n)
   {
      DataArray < int > perm(num());
      remove(nums, n, perm.get_ptr());
   }

   ///
   void remove(const DataKey keys[], int n, int* perm)
   {
      for( int i = num() - 1; i >= 0; --i )
         perm[i] = i;

      while( n-- )
         perm[number(*keys++)] = -1;

      remove(perm);
   }

   /// Removes \p n SVectorBase%s from %set.
   /** @pre \p perm must be at least of size num()
    *  @pre \p nums must be at least of size \p n @return \p perm is the permutations resulting from this removal: \p
    *          perm[i] < 0 indicates, that the element to index \c i has been removed. Otherwise, \p perm[i] is the new
    *          index of the element with index \c i before the removal.
    */
   void remove(const int nums[], int n, int* perm)
   {
      for( int i = num() - 1; i >= 0; --i )
         perm[i] = i;

      while( n-- )
         perm[*nums++] = -1;

      remove(perm);
   }

   /// Removes all SVectorBase%s from %set.
   void clear()
   {
      SVSetBaseArray::clear();
      SVSetBaseArray::reMax(10000);
      set.clear();
      list.clear();
   }

   //@}

   // ------------------------------------------------------------------------------------------------------------------
   /**@name Access */
   //@{

   /// Gets SVectorBase by number, writeable.
   SVectorBase<R>& operator[](int n)
   {
      return set[n];
   }

   /// Gets SVectorBase by number.
   const SVectorBase<R>& operator[](int n) const
   {
      return set[n];
   }

   /// Gets SVectorBase by DataKey, writeable.
   SVectorBase<R>& operator[](const DataKey& k)
   {
      return set[k];
   }

   /// Gets SVectorBase by DataKey.
   const SVectorBase<R>& operator[](const DataKey& k) const
   {
      return set[k];
   }

   //@}

   // ------------------------------------------------------------------------------------------------------------------
   /**@name Inquiry */
   //@{

   /// Current number of SVectorBase%s.
   int num() const
   {
      return set.num();
   }

   /// Current maximum number of SVectorBase%s.
   int max() const
   {
      return set.max();
   }

   /// Gets DataKey of vector number.
   DataKey key(int n) const
   {
      return set.key(n);
   }

   /// Gets DataKey of SVectorBase.
   DataKey key(const SVectorBase<R>* svec) const
   {
      return set.key(static_cast<const DLPSV*>(svec));
   }

   /// Gets vector number of DataKey.
   int number(const DataKey& k) const
   {
      return set.number(k);
   }

   /// Gets vector number of SVectorBase.
   int number(const SVectorBase<R>* svec) const
   {
      return set.number(static_cast<const DLPSV*>(svec));
   }

   /// True iff SVSetBase contains a SVectorBase for DataKey \p k.
   bool has(const DataKey& k) const
   {
      return set.has(k);
   }

   /// True iff SVSetBase contains a SVectorBase for vector number n.
   bool has(int n) const
   {
      return set.has(n);
   }

   /// Is an SVectorBase in the %set?
   bool has(const SVectorBase<R>* svec) const
   {
      return set.has(static_cast<const DLPSV*>(svec));
   }

   //@}

   // ------------------------------------------------------------------------------------------------------------------
   /**@name Memory Management */
   //@{

   /// Used nonzero memory.
   int memSize() const
   {
      return SVSetBaseArray::size();
   }

   /// Length of nonzero memory.
   int memMax() const
   {
      return SVSetBaseArray::max();
   }

   /// Reset length of nonzero memory.
   void memRemax(int newmax)
   {
      ptrdiff_t delta = SVSetBaseArray::reMax(newmax);

      if( delta != 0 )
      {
         for( DLPSV* ps = list.first(); ps; ps = list.next(ps) )
         {
            // get new shifted nonzero memory of the SVectorBase
            Element<R>* newmem = reinterpret_cast<Element<R>*>(reinterpret_cast<char*>(ps->mem()) + delta);

            // get the size and maximum capacity of the SVectorBase
            int sz = ps->size();
            int l_max = ps->max();
            assert(l_max >= sz);

            // set new nonzero memory
            ps->setMem(l_max, newmem);
            ps->set_size(sz);
         }
      }
   }

   /// Garbage collection in nonzero memory.
   /** Pack the svectors together as tightly as possible. This removes all additional unused memory, i.e., size = max
    *  for every svector after the call.
    *
    *  Note: do *not* call isConsistent() here, because the following might happen: In SPxLP::doAddRows(const LPRowSet&
    *  p_set), when adding rows, the sizes of the vectors for the columns of the LP are increased (without yet filling
    *  in the data) to recieve the additional entries. This is done by calling xtend() above. xtend() in turn might call
    *  this method, which checks the yet unfilled positions, i.e., isConsistent() is likely to fail. In general,
    *  isConsistent() should not be called within this class, but in classes further up in the hierarchy.
    */
   void memPack()
   {
      DLPSV* ps;
      int used;
      int j;

      for( used = 0, ps = list.first(); ps; ps = list.next(ps) )
      {
         const int sz = ps->size();

         if( ps->mem() != &this->SVSetBaseArray::operator[](used) )
         {
            // cannot use memcpy, because the memory might overlap
            for( j = 0; j < sz; ++j )
               this->SVSetBaseArray::operator[](used + j) = ps->mem()[j];

            ps->setMem(sz, &this->SVSetBaseArray::operator[](used));
            ps->set_size(sz);
         }
         else
            ps->set_max(sz);

         used += sz;
      }

      SVSetBaseArray::reSize(used);

      possiblyUnusedMem = 0;
   }

   //@}

   // ------------------------------------------------------------------------------------------------------------------
   /**@name Miscellaneous */
   //@{

   /// Resets maximum number of SVectorBase%s.
   void reMax(int newmax = 0)
   {
      list.move(set.reMax(newmax));
   }

   /// Consistency check.
   bool isConsistent() const
   {
#ifdef ENABLE_CONSISTENCY_CHECKS
      DLPSV* ps;
      DLPSV* next;

      for( ps = list.first(); ps; ps = next )
      {
         if( !ps->isConsistent() )
            return MSGinconsistent("SVSetBase");

         if( ps->mem() > &last() )
            return MSGinconsistent("SVSetBase");

         next = list.next(ps);

         if( next && ps->mem() + ps->max() != next->mem() )
            return MSGinconsistent("SVSetBase");
      }

      return SVSetBaseArray::isConsistent() && set.isConsistent() && list.isConsistent();
#else
      return true;
#endif
   }

   //@}

   // ------------------------------------------------------------------------------------------------------------------
   /**@name Constructors / destructors */
   //@{

   /// Default constructor.
   explicit
   SVSetBase<R>(int pmax = -1, int pmemmax = -1, double pfac = 1.1, double pmemFac = 1.2)
      : SVSetBaseArray(0, (pmemmax > 0) ? pmemmax : 8 * ((pmax > 0) ? pmax : 8), pmemFac)
      , set((pmax > 0) ? pmax : 8)
      , possiblyUnusedMem(0)
      , factor(pfac)
   {
      assert(isConsistent());
   }

   /// Destructor
   ~SVSetBase<R>()
   {}

   /// Assignment operator.
   SVSetBase<R>& operator=(const SVSetBase<R>& rhs)
   {
      if( this != &rhs )
      {
         clear();

         if( rhs.size() > 0 )
         {
            SVSetBaseArray::operator=(rhs);
            set = rhs.set;

            DLPSV* ps;
            DLPSV* newps;

            void* delta0 = &(*(static_cast<SVSetBaseArray*>(this)))[0];
            void* delta1 = &(*(static_cast<SVSetBaseArray*>(const_cast<SVSetBase<R>*>(&rhs))))[0];
            ptrdiff_t delta = reinterpret_cast<char*>(delta0) - reinterpret_cast<char*>(delta1);

            for( ps = rhs.list.first(); ps; ps = rhs.list.next(ps) )
            {
               newps = &set[rhs.number(ps)];
               list.append(newps);
               newps->setMem(ps->max(),
                  reinterpret_cast<Element<R>*>(reinterpret_cast<char*>(ps->mem()) + delta));
               newps->set_size(ps->size());
            }
         }
      }

      assert(isConsistent());

      return *this;
   }

#if 0 ///@todo implement correctly; we need cast in SVectorBase::Element and build up SVSetBase from scratch
   /// Assignment operator.
   template < class S >
   SVSetBase<R>& operator=(const SVSetBase<S>& rhs)
   {
      if( this != &rhs )
      {
         clear();

         if( rhs.size() > 0 )
         {
            SVSetBaseArray::operator=(rhs);
            set = rhs.set;

            DLPSV* ps;
            DLPSV* newps;

            void* delta0 = &(*(static_cast<SVSetBaseArray*>(this)))[0];
            void* delta1 = &(*(static_cast<SVSetBaseArray*>(const_cast<SVSetBase<R>*>(&rhs))))[0];
            ptrdiff_t delta = reinterpret_cast<char*>(delta0) - reinterpret_cast<char*>(delta1);

            for( ps = rhs.list.first(); ps; ps = rhs.list.next(ps) )
            {
               newps = &set[rhs.number(ps)];
               list.append(newps);
               newps->setMem(ps->max(),
                  reinterpret_cast<Element<R>*>(reinterpret_cast<char*>(ps->mem()) + delta));
               newps->set_size(ps->size());
            }
         }
      }

      assert(isConsistent());

      return *this;
   }
#endif

   /// Copy constructor.
   SVSetBase<R>(const SVSetBase<R>& old)
      : SVSetBaseArray()
      , possiblyUnusedMem(old.possiblyUnusedMem)
      , factor(old.factor)
   {
      *this = old;

      assert(SVSetBase::isConsistent());
   }

#if 0
   /// Copy constructor.
   template < class S >
   SVSetBase<R>(const SVSetBase<S>& old)
      : SVSetBaseArray()
      , possiblyUnusedMem(old.possiblyUnusedMem)
      , factor(old.factor)
   {
      *this = old;

      assert(SVSetBase::isConsistent());
   }
#endif

   //@}
};

} // namespace soplex
#endif // _SVSETBASE_H_

//-----------------------------------------------------------------------------
//Emacs Local Variables:
//Emacs mode:c++
//Emacs c-basic-offset:3
//Emacs tab-width:8
//Emacs indent-tabs-mode:nil
//Emacs End:
//-----------------------------------------------------------------------------
