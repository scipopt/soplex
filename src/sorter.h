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
#pragma ident "@(#) $Id: sorter.h,v 1.2 2001/11/06 23:31:03 bzfkocht Exp $"


#ifndef _SORTER_H_
#define _SORTER_H_

//@ ----------------------------------------------------------------------------
/*      \Section{Imports}
    Import required system include files
 */
#include <assert.h>

namespace soplex
{


//@ ----------------------------------------------------------------------------
/* \Section{Implementation}
    This is quicksort.
 */
/**@name sorting        functions */
//@{
/** Sort array.
    The template function #sort# sorts an array #t# holding #n# elements of type
    #T# using #c# for comparisions. Class #COMPARATOR# must be provide an
    overloaded #operator()(const T& t1,const T& t2)#, that returns
    \begin{description}
    \item[#<0#]         if #t1# is to appear before #t2#,
    \item[#==0#]        if #t1# and #t2# can appear in any order or
    \item[#>0#]         if #t1# is to appear after #t2#.
    \end{description}
 */

template < class T, class COMPARATOR >
void sorter_qsort(T* t, int end, COMPARATOR& compare, int start = 0)
{
   int i0, i1, j;
   double c;

   T work, mid, tmp;

   work = t[start];
   t[start] = t[(start + end) / 2];
   t[(start + end) / 2] = work;
   
   mid = t[start];
   work = t[end - 1];
   
   for (i0 = i1 = start, j = end - 1; i1 < j;)
   {
      c = compare(mid, work);
      if (c > 0)
      {
         tmp = t[i0];
         t[i0] = work;
         i0++;
         i1++;
         work = t[i1];
         t[i1] = tmp;
      }
      else if (c < 0)
      {
         t[j] = work;
         --j;
         work = t[j];
      }
      else
      {
         i1++;
         tmp = t[i1];
         t[i1] = work;
         work = tmp;
      }
   }
   
   if (start < i0 - 1)
      sorter_qsort(t, i0, compare, start);
   if (i1 + 1 < end)
      sorter_qsort(t, end, compare, i1 + 1);
}

//@}


} // namespace soplex
#endif // _SORTER_H_

//-----------------------------------------------------------------------------
//Emacs Local Variables:
//Emacs mode:c++
//Emacs c-basic-offset:3
//Emacs tab-width:8
//Emacs indent-tabs-mode:nil
//Emacs End:
//-----------------------------------------------------------------------------
