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
#pragma ident "@(#) $Id: datakey.h,v 1.1 2001/11/11 20:27:29 bzfkocht Exp $"

#ifndef _DATAKEY_H_
#define _DATAKEY_H_

namespace soplex
{
   /** Entry Identifiers.
       Every item in a #DataSet# is assigned a #Key# by which it can be
       accessed (using #DataSet::operator[]#). A #Key# consists of an integer
       member #idx#, which is a positive number for any valid #Key#. No
       #Key::idx# of an element in a #DataSet# may exceed the sets #max()#.
       This property may be used to build arrays with additional information to
       the elements of a #DataSet#.

       In addition, #Key#s provides member #info#, where the programmer is
       free to store other information.

       Each #Key# is unique for one #DataSet# but different #DataSet#s may (and
       generally will) manage the same #Key#s. When an element is removed from
       a #DataSet# its #Key# may (and generally will) be reused for other
       elements added to the #DataSet# later on.
    */
class DataKey
{
public: // Should be private
   signed int info: 8;
   signed int idx : (8 * sizeof(int) - 8);

public:
   ///
   inline int getInfo() const
   {
      return info;
   }
   ///
   inline int getIdx() const
   {
      return info;
   }
   ///
   inline void setInfo(int p_info) 
   {
      info = p_info;
   }
   ///
   inline void setIdx(int p_idx) 
   {
      idx = p_idx;
   }   
   ///
   inline int isValid() const
   {
      return idx >= 0;
   }
   ///
   inline void inValidate()
   {
      idx = -1;
      info = 0;
   }
   ///    
   DataKey() : info(0), idx(-1) 
   {}
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
