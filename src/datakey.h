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
#pragma ident "@(#) $Id: datakey.h,v 1.4 2001/11/21 20:15:40 bzfkocht Exp $"

/**@file  datakey.h
 * @brief Entry identifier class for items of a #DataSet.
 */
#ifndef _DATAKEY_H_
#define _DATAKEY_H_

#include <assert.h>

namespace soplex
{
/**@brief   Entry identifier class for items of a #DataSet.
   @ingroup Elementary

   Every item in a #DataSet is assigned a #Key by which it can be
   accessed (using #DataSet::operator[]()). A #Key consists of an integer
   member #idx, which is a positive number for any valid #Key. No
   #Key::idx of an element in a #DataSet may exceed the sets #max().
   This property may be used to build arrays with additional information to
   the elements of a #DataSet.

   In addition, #Key%s provides member #info, where the programmer is
   free to store other information.
   
   Each #Key is unique for one #DataSet but different #DataSet%s may (and
   generally will) manage the same #Key%s. When an element is removed from
   a #DataSet its #Key may (and generally will) be reused for other
   elements added to the #DataSet later on.

   @todo data members should be private.
*/
class DataKey
{
public:
   signed int info: 8;                        ///< user information (8 bit)
   signed int idx : (8 * sizeof(int) - 8);    ///< (locally) unique key index

public:
   /// gets the user information (#info) attached to the #Key.
   inline int getInfo() const
   {
      return info;
   }
   /// gets the index number (#idx) of the #Key.
   inline int getIdx() const
   {
      return idx;
   }
   /// sets the user information (#info) attached to the #Key.
   inline void setInfo(int p_info) 
   {
      info = p_info;
   }
   /// sets the index number (#idx) of the #Key.
   inline void setIdx(int p_idx) 
   {
      idx = p_idx;
   }   
   /// returns TRUE, iff the #Key is valid.
   inline int isValid() const
   {
      return idx >= 0;
   }
   /// makes the #Key invalid and clears the #info field.
   inline void inValidate()
   {
      idx = -1;
      info = 0;
   }
   /// Default constructor. Constructs an invalid #Key.
   DataKey() 
      : info(0), idx(-1) 
   {}
   /// Assignment operator.
   DataKey& operator=(const DataKey& rhs)
   {
      assert(sizeof(*this) == sizeof(int));
      // dirty implementation of the week.
      *reinterpret_cast<int*>(this) = *reinterpret_cast<const int*>(&rhs);

      return *this;
   }
#if 0
   /// Copy constructor.
   explicit DataKey(const DataKey& old) 
      : info(old.info) 
      , idx(old.idx)
   {}
#endif
};

} // namespace soplex
#endif // _DATAKEY_H_

//-----------------------------------------------------------------------------
//Emacs Local Variables:
//Emacs mode:c++
//Emacs c-basic-offset:3
//Emacs tab-width:8
//Emacs indent-tabs-mode:nil
//Emacs End:
//-----------------------------------------------------------------------------
