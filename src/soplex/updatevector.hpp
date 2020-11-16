/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the class library                   */
/*       SoPlex --- the Sequential object-oriented simPlex.                  */
/*                                                                           */
/*    Copyright (C) 1996-2020 Konrad-Zuse-Zentrum                            */
/*                            fuer Informationstechnik Berlin                */
/*                                                                           */
/*  SoPlex is distributed under the terms of the ZIB Academic Licence.       */
/*                                                                           */
/*  You should have received a copy of the ZIB Academic License              */
/*  along with SoPlex; see the file COPYING. If not email to soplex@zib.de.  */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

namespace soplex
{

template <class R>
UpdateVector<R>& UpdateVector<R>::operator=(const UpdateVector<R>& rhs)
{
   if(this != &rhs)
   {
      theval   = rhs.theval;
      thedelta = rhs.thedelta;
      VectorBase<R>::operator=(rhs);

      assert(UpdateVector<R>::isConsistent());
   }

   return *this;
}

template <class R>
UpdateVector<R>::UpdateVector(const UpdateVector<R>& base)
   : VectorBase<R>(base)
   , theval(base.theval)
   , thedelta(base.thedelta)
{
   assert(UpdateVector<R>::isConsistent());
}

template <class R>
bool UpdateVector<R>::isConsistent() const
{
#ifdef ENABLE_CONSISTENCY_CHECKS

   if(this->dim() != thedelta.dim())
      return MSGinconsistent("UpdateVector");

   return VectorBase<R>::isConsistent() && thedelta.isConsistent();
#else
   return true;
#endif
}
} // namespace soplex
