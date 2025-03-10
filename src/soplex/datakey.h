/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the class library                   */
/*       SoPlex --- the Sequential object-oriented simPlex.                  */
/*                                                                           */
/*  Copyright (c) 1996-2025 Zuse Institute Berlin (ZIB)                      */
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

/**@file  datakey.h
 * @brief Entry identifier class for items of a DataSet.
 */
#ifndef _DATAKEY_H_
#define _DATAKEY_H_

#include <assert.h>

namespace soplex
{
/**@brief   Entry identifier class for items of a DataSet.
   @ingroup Elementary

   Every item in a DataSet is assigned a DataKey by which it can be
   accessed (using DataSet::operator[]()). A DataKey consists of an integer
   member #idx, which is a positive number for any valid DataKey. No
   #idx of an element in a DataSet may exceed the sets max().
   This property may be used to build arrays with additional information to
   the elements of a DataSet.

   In addition, #DataKey%s provide a member #info which can be used to store
   further information.

   Each DataKey is unique for one DataSet but different DataSets may (and
   generally will) manage the same #DataKey%s. When an element is removed from
   a DataSet its DataKey may (and generally will) be reused for other
   elements added to the DataSet later on.

   @todo data members should be private.
*/
class DataKey
{
public:

   //-------------------------------------
   /**@name Data */
   ///@{
   /* This was originally implemented as bitfield "signed int info: 2; signed int idx: (8 * sizeof(int) - 2);",
      however, this seems to trigger a bug with old versions of GCC/glibc on 32bit machines. */
   int info;                                  ///< user information to store values -1, 0, +1
   int idx;                                   ///< (locally) unique key index
   ///@}

public:

   //-------------------------------------
   /**@name Constructors / destructors */
   ///@{
   /// Default constructor. Constructs an invalid DataKey.
   DataKey()
      : info(0), idx(-1)
   {}
   // Full constructor
   DataKey(int p_info, int p_idx)
      : info(p_info)
      , idx(p_idx)
   {
      assert(p_info <= 1 && p_info >= -1);
   }

   ///@}

   //-------------------------------------
   /**@name Access / modification */
   ///@{
   /// gets the index number (\ref soplex::DataKey::idx "idx") of the DataKey.
   inline int getIdx() const
   {
      return idx;
   }
   /// sets the index number (\ref soplex::DataKey::idx "idx") of the DataKey.
   inline void setIdx(int p_idx)
   {
      idx = p_idx;
   }
   /// returns TRUE, iff the DataKey is valid.
   inline bool isValid() const
   {
      return idx >= 0;
   }
   /// makes the DataKey invalid and clears the \ref soplex::DataKey::info "info" field.
   inline void inValidate()
   {
      idx  = -1;
      info = 0;
   }
   ///@}

};

} // namespace soplex
#endif // _DATAKEY_H_
