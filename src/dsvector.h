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
#pragma ident "@(#) $Id: dsvector.h,v 1.2 2001/11/06 23:31:01 bzfkocht Exp $"

/*      \Section{Imports}
 */
#ifndef _DSVECTOR_H_
#define _DSVECTOR_H_

#include <assert.h>

#include "svector.h"
#include "ssvector.h"

namespace soplex
{

/*      \Section{Class Declaration}
 */
/** Dynamic sparse vectors.
    Class #DSVector# implements dynamic sparse vectors, i.e. #SVector#s
    with an automatic memory management. This allows the user to freely #add()#
    as many nonzeros to a #DSVector# as desired, without any precautions.
    For saving memory method #setMax()# allows to reduce memory consumption to
    the amount really required.
 */
class DSVector : public SVector
{
   friend class SLinSolver;

   Element* theelem;                // here is where the memory is
   int* mem();

   void allocMem(int);
   void makeMem(int n)
   {
      if (max() - size() < ++n)
         setMax(size() + n);
   }

public:
   ///
   DSVector& operator=(const SSVector& sv)
   {
      int n = sv.size();
      clear();
      makeMem(n);
      SVector::operator=(sv);
      return *this;
   }
   ///
   DSVector& operator=(const SVector& sv)
   {
      int n = sv.size();
      clear();
      makeMem(n);
      SVector::operator=(sv);
      return *this;
   }
   ///
   DSVector& operator=(const DSVector& sv)
   {
      int n = sv.size();
      clear();
      makeMem(n);
      SVector::operator=(sv);
      return *this;
   }
   ///
   DSVector& operator=(const Vector& vec);
   ///
   DSVector& assign(const Vector& vec, double eps = 1e-16);

   ///
   void add(const SVector& sv)
   {
      int n = sv.size();
      clear();
      makeMem(n);
      SVector::add(sv);
   }

   ///
   void add(int i, double v)
   {
      makeMem(1);
      SVector::add(i, v);
   }

   ///
   void add(int n, const int *i, const double *v)
   {
      makeMem(n);
      SVector::add(n, i, v);
   }

   /** reset nonzero memory to #>= newmax#.
       This methods resets the memory consumption of the #DIdxSet# to
       #newmax#. However, if #newmax < size()#, it is reset to #size()#
       only.
    */
   void setMax(int newmax = 1);

   ///
   DSVector(const Vector& vec);

   ///
   DSVector(const SVector& old);
   ///
   DSVector(const DSVector& old);

   /** Default constructor.
       Creates a #DSVector# ready to hold #n# nonzeros. However, the memory is
       automatically enlarged, if more nonzeros are added to the #DSVector#.
    */
   DSVector(int n = 8);

   ///
   ~DSVector();

   ///
   int isConsistent() const;
};

} // namespace soplex
#endif // _DSVECTOR_H_

//-----------------------------------------------------------------------------
//Emacs Local Variables:
//Emacs mode:c++
//Emacs c-basic-offset:3
//Emacs tab-width:8
//Emacs indent-tabs-mode:nil
//Emacs End:
//-----------------------------------------------------------------------------
