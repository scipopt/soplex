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
#pragma ident "@(#) $Id: updatevector.cpp,v 1.3 2001/12/25 16:03:25 bzfkocht Exp $"

#include "updatevector.h"
#include "spxmessage.h"

namespace soplex
{

UpdateVector& UpdateVector::operator=(const UpdateVector& rhs)
{
   theval = rhs.theval;
   thedelta = rhs.thedelta;
   DVector::operator=(rhs);
   return *this;
}

int UpdateVector::isConsistent() const
{
   if (dim() != thedelta.dim())
      return SPXinconsistent("UpdateVector");

   return DVector::isConsistent() && thedelta.isConsistent();
}
} // namespace soplex

//-----------------------------------------------------------------------------
//Emacs Local Variables:
//Emacs mode:c++
//Emacs c-basic-offset:3
//Emacs tab-width:8
//Emacs indent-tabs-mode:nil
//Emacs End:
//-----------------------------------------------------------------------------
