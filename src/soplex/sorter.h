/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the class library                   */
/*       SoPlex --- the Sequential object-oriented simPlex.                  */
/*                                                                           */
/*    Copyright (C) 1996-2018 Konrad-Zuse-Zentrum                            */
/*                            fuer Informationstechnik Berlin                */
/*                                                                           */
/*  SoPlex is distributed under the terms of the ZIB Academic Licence.       */
/*                                                                           */
/*  You should have received a copy of the ZIB Academic License              */
/*  along with SoPlex; see the file COPYING. If not email to soplex@zib.de.  */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/**@file  sorter.h
 * @brief Wrappers for std::sort implementation.
 */
#ifndef _SORTER_H_
#define _SORTER_H_

#include <algorithm>
#include <assert.h>

namespace soplex
{
  // A wrapper for std::sort that works like the old SPxQuicksort, except for
  // the extra variable 'type'. I don't know what 'type' does.
  template <class T, class CMPR>
  void spxSort(T* keys, int end, CMPR& compare, int start = 0)
  {
    // An auxiliary comparison function that returns a bool instead of an
    // integer. This will be used for std::sort.
    //
    // std::sort sorts according to ascending order (of the operator <). But compare(a, b)
    // behaves differently. if compare(a, b) < 0, then a should appear before b.
    // = 0, then order doesn't matter and if > 0, then a should appear after b.
    // The auxiliary function does this translation.
    auto auxCmpr = [compare](T a, T b)
                   {
                     if(compare(a, b) <= 0)
                       {
                         return true;
                       }
                     else
                       {
                         return false;
                       }
                   };
    std::sort(keys + start, keys + end, auxCmpr);
  }

  // An wrapper for std::partial sort in the style of SPxQuicksortPart. Again
  // without the variable 'type'.
  template <class T, class CMPR>
  int spxSortPart(T* keys, CMPR& compare, int start, int end, int size)
  {
    // See the comments for spxSort
    auto auxCmpr = [compare](T a, T b)
                   {
                     if(compare(a, b) <= 0)
                       {
                         return true;
                       }
                     else
                       {
                         return false;
                       }
                   };

    // All elements must be sorted
    if(size + start >= end - 1)
      {
        spxSort(keys, end, auxCmpr, start);
        return (end - 1);
      }

    std::partial_sort(keys + start, keys + size, keys + end, auxCmpr);
    // This might be a drawback of using std::partial_sort. The original
    // algorithm may have a better return value that is greater than size.
    return size;

  }

} // namespace soplex
#endif // _SORTER_H_
