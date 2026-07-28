/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the class library                   */
/*       SoPlex --- the Sequential object-oriented simPlex.                  */
/*                                                                           */
/*  Copyright (c) 1996-2026 Zuse Institute Berlin (ZIB)                      */
/*                                                                           */
/*  Licensed under the Apache License, Version 2.0 (the "License");          */
/*  you may not use this file except in compliance with the License.         */
/*  You may obtain a copy of the License at                                  */
/*                                                                           */
/*      http://www.apache.org/licenses/LICENSE-2.0                           */
/*                                                                           */
/*  Unless required by applicable law or agreed to in writing, software      */
/*  distributed under the License is distributed on an "AS IS" BASIS,        */
/*  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied. */
/*  See the License for the specific language governing permissions and      */
/*  limitations under the License.                                           */
/*                                                                           */
/*  You should have received a copy of the Apache-2.0 license                */
/*  along with SoPlex; see the file LICENSE. If not email to soplex@zib.de.  */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/**@file  sorter.h
 * @brief Generic QuickSort implementation.
 */
#ifndef _SORTER_H_
#define _SORTER_H_

#include <assert.h>

namespace soplex
{
#define SOPLEX_SHELLSORTMAX 25

/** shell-sort an array of data elements; use it only for arrays smaller than 25 entries */
template < class T, class COMPARATOR >
void SPxShellsort(T* keys, int end, COMPARATOR& compare, int start = 0)
{
   static const int incs[3] = {1, 5, 19}; // sequence of increments

   assert(start >= 0);
   assert(end >= start);

   for(int k = 2; k >= 0; --k)
   {
      int h = incs[k];
      int first = start + h;

      for(int i = first; i <= end; ++i)
      {
         T tempkey = keys[i];
         int j = i;

         while(j >= first && compare(tempkey, keys[j - h]) < 0)
         {
            keys[j] = keys[j - h];
            j -= h;
         }

         keys[j] = tempkey;
      }
   }
}


/// Generic QuickSort implementation.
/** This template function sorts an array \p t holding \p n elements of
    type T using \p compare for comparisons. Class COMPARATOR must provide
    an overloaded operator()(const T& t1,const T& t2) which returns
    - < 0, if \p t1 is to appear before \p t2,
    - = 0, if \p t1 and \p t2 can appear in any order, or
    - > 0, if \p t1 is to appear after \p t2.
*/
template < class T, class COMPARATOR >
void SPxQuicksort(T* keys, int end, COMPARATOR& compare, int start = 0, bool type = true)
{
   assert(start >= 0);
   assert(end >= start);

   /* nothing to sort for at most one element */
   if(end - start <= 1)
      return;

   /* reduce end position to last element index */
   --end;

   /* use quick-sort for long lists */
   while(end - start >= SOPLEX_SHELLSORTMAX)
   {
      T pivotkey;
      T tmp;
      int lo;
      int hi;
      int mid;

      /* select pivot element */
      mid = start + (end - start) / 2; // avoid overflowing (start + end) / 2
      pivotkey = keys[mid];

      /* partition the array into elements < pivot [start,hi] and elements >= pivot [lo,end] */
      lo = start;
      hi = end;

      for(;;)
      {
         if(type)
         {
            while(lo < end && compare(keys[lo], pivotkey) < 0)
               ++lo;

            while(hi > start && compare(keys[hi], pivotkey) >= 0)
               --hi;
         }
         else
         {
            while(lo < end && compare(keys[lo], pivotkey) <= 0)
               ++lo;

            while(hi > start && compare(keys[hi], pivotkey) > 0)
               --hi;
         }

         if(lo >= hi)
            break;

         tmp = keys[lo];
         keys[lo] = keys[hi];
         keys[hi] = tmp;

         ++lo;
         --hi;
      }

      assert((hi == lo - 1) || (type && hi == start) || (!type && lo == end));

      /* skip entries which are equal to the pivot element (three partitions, <, =, > than pivot)*/
      if(type)
      {
         while(lo < end && compare(pivotkey, keys[lo]) >= 0)
            ++lo;

         /* make sure that we have at least one element in the smaller partition */
         if(lo == start)
         {
            /* everything is greater or equal than the pivot element: move pivot to the left (degenerate case) */
            assert(compare(keys[mid], pivotkey) == 0); /* the pivot element did not change its position */

            tmp = keys[lo];
            keys[lo] = keys[mid];
            keys[mid] = tmp;

            ++lo;
         }
      }
      else
      {
         while(hi > start && compare(pivotkey, keys[hi]) <= 0)
            --hi;

         /* make sure that we have at least one element in the smaller partition */
         if(hi == end)
         {
            /* everything is greater or equal than the pivot element: move pivot to the left (degenerate case) */
            assert(compare(keys[mid], pivotkey) == 0); /* the pivot element did not change its position */

            tmp = keys[hi];
            keys[hi] = keys[mid];
            keys[mid] = tmp;

            --hi;
         }
      }

      /* sort the smaller partition by a recursive call, sort the larger part without recursion */
      if(hi - start <= end - lo)
      {
         /* sort [start,hi] with a recursive call */
         if(start < hi)
         {
            SPxQuicksort(keys, hi + 1, compare, start, !type);
         }

         /* now focus on the larger part [lo,end] */
         start = lo;
      }
      else
      {
         if(lo < end)
         {
            SPxQuicksort(keys, end + 1, compare, lo, !type);
         }

         /* now focus on the larger part [start,hi] */
         end = hi;
      }

      type = !type;
   }

   /* use shell sort on the remaining small list */
   if(end - start >= 1)
   {
      SPxShellsort(keys, end, compare, start);
   }

#ifdef CHECK_SORTING

   for(int i = start; i < end; ++i)
      assert(compare(keys[i], keys[i + 1]) <= 0);

#endif
}


/**@brief  Generic implementation of Partial QuickSort.
 *
 * This template function sorts an array \p t holding \p n elements of
 * type T partially using \p compare for comparisons, i.e. ensures that
 * the \p size smallest elements are sorted to the front.
 *
 * Class COMPARATOR must provide an overloaded
 * operator()(const T& t1,const T& t2) which returns
 * - < 0, if \p t1 is to appear before \p t2,
 * - = 0, if \p t1 and \p t2 can appear in any order, or
 * - > 0, if \p t1 is to appear after \p t2.
 *
 * @param keys               array of elements to be sorted between index start and end
 * @param compare            comparator
 * @param start              index of first element in range to be sorted
 * @param end                index of last element in range to be sorted plus 1
 * @param size               guaranteed number of additionally sorted elements
 * @param start2             auxiliary start index of sub range used for recursive call (deprecated)
 * @param end2               auxiliary end index of sub range used for recursive call (disabled)
 * @param type               type of sorting, to be more flexable on degenerated cases
 * @return                   index of last element in range sorted plus 1
 */
template < class T, class COMPARATOR >
int SPxQuicksortPart(T* keys, COMPARATOR& compare, int start, int end, int size, int start2 = 0,
                     int end2 = 0, bool type = true)
{
   assert(start >= 0);
   assert(end >= start);
   assert(start2 <= end);
   assert(end2 <= end);

   /* nothing to sort for at most one element */
   if(end - start <= 1)
      return end;

   /* we assume that range {start, ..., start2-1} already contains start2-start smallest elements in sorted order */
#ifdef CHECK_SORTING

   for(int i = start; i < start2 - 1; ++i)
      assert(compare(keys[i], keys[i + 1]) <= 0);

#endif

   /* skip sorted prefix */
   if(start < start2)
      start = start2;

   /* if all remaining elements should be sorted, we simply call standard quicksort */
   if(start >= end - size - 1)
   {
      SPxQuicksort(keys, end, compare, start, type);
      return end;
   }

   T pivotkey;
   T tmp;
   int lo;
   int hi;
   int mid;

   /* reduce end position to last element index */
   --end;

   /* select pivot element */
   mid = start + (end - start) / 2; // avoid overflowing (start + end) / 2
   pivotkey = keys[mid];

   /* partition the array into elements < pivot [start,hi] and elements >= pivot [lo,end] */
   lo = start;
   hi = end;

   for(;;)
   {
      if(type)
      {
         while(lo < end && compare(keys[lo], pivotkey) < 0)
            ++lo;

         while(hi > start && compare(keys[hi], pivotkey) >= 0)
            --hi;
      }
      else
      {
         while(lo < end && compare(keys[lo], pivotkey) <= 0)
            ++lo;

         while(hi > start && compare(keys[hi], pivotkey) > 0)
            --hi;
      }

      if(lo >= hi)
         break;

      tmp = keys[lo];
      keys[lo] = keys[hi];
      keys[hi] = tmp;

      ++lo;
      --hi;
   }

   assert((hi == lo - 1) || (type && hi == start) || (!type && lo == end));

   /* skip entries which are equal to the pivot element (three partitions, <, =, > than pivot)*/
   if(type)
   {
      while(lo < end && compare(pivotkey, keys[lo]) >= 0)
         ++lo;

      /* make sure that we have at least one element in the smaller partition */
      if(lo == start)
      {
         /* everything is greater or equal than the pivot element: move pivot to the left (degenerate case) */
         assert(compare(keys[mid], pivotkey) == 0); /* the pivot element did not change its position */

         tmp = keys[lo];
         keys[lo] = keys[mid];
         keys[mid] = tmp;

         ++lo;
      }
   }
   else
   {
      while(hi > start && compare(pivotkey, keys[hi]) <= 0)
         --hi;

      /* make sure that we have at least one element in the smaller partition */
      if(hi == end)
      {
         /* everything is greater or equal than the pivot element: move pivot to the left (degenerate case) */
         assert(compare(keys[mid], pivotkey) == 0); /* the pivot element did not change its position */

         tmp = keys[hi];
         keys[hi] = keys[mid];
         keys[mid] = tmp;

         --hi;
      }
   }

#ifdef CHECK_SORTING

   for(int i = start; i < lo; ++i)
      assert(compare(keys[i], pivotkey) <= 0);

#endif

   /* if we only need to sort less than half of the "<" part, use partial sort again */
   if(start <= hi - 2 * size)
   {
      return SPxQuicksortPart(keys, compare, start, hi + 1, size, start2, end2, !type);
   }
   /* otherwise, and if we do not need to sort the ">" part, use standard quicksort on the "<" part */
   else if(start <= lo - size)
   {
      SPxQuicksort(keys, hi + 1, compare, start, !type);
      return lo;
   }
   /* otherwise we have to sort the "<" part fully (use standard quicksort) and the ">" part partially */
   else
   {
      SPxQuicksort(keys, hi + 1, compare, start, !type);
      return SPxQuicksortPart(keys, compare, lo, end + 1, start + size - lo, start2, end2, !type);
   }
}

} // namespace soplex
#endif // _SORTER_H_
