/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the class library                   */
/*       SoPlex --- the Sequential object-oriented simPlex.                  */
/*                                                                           */
/*    Copyright (C) 1997-1999 Roland Wunderling                              */
/*                  1997-2002 Konrad-Zuse-Zentrum                            */
/*                            fuer Informationstechnik Berlin                */
/*                                                                           */
/*  SoPlex is distributed under the terms of the ZIB Academic Licence.       */
/*                                                                           */
/*  You should have received a copy of the ZIB Academic License              */
/*  along with SoPlex; see the file COPYING. If not email to soplex@zib.de.  */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
#pragma ident "@(#) $Id: spxdesc.cpp,v 1.13 2005/07/14 13:37:52 bzforlow Exp $"

//#define DEBUGGING 1

#include <iostream>

#include "spxdefines.h"
#include "spxbasis.h"
#include "spxsolver.h"

namespace soplex
{
 
SPxBasis::Desc::Desc(const SPxSolver& base)
{
   rowstat.reSize(base.nRows());
   colstat.reSize(base.nCols());

   if (base.rep() == SPxSolver::ROW)
   {
      stat   = &rowstat;
      costat = &colstat;
   }
   else
   {
      assert(base.rep() == SPxSolver::COLUMN);

      stat   = &colstat;
      costat = &rowstat;
   }
}

SPxBasis::Desc::Desc(const Desc& old)
   : rowstat(old.rowstat)
   , colstat(old.colstat)
{
   if (old.stat == &old.rowstat)
   {
      assert(old.costat == &old.colstat);
      
      stat   = &rowstat;
      costat = &colstat;
   }
   else
   {
      assert(old.costat == &old.rowstat);
      
      stat   = &colstat;
      costat = &rowstat;
   }
}

SPxBasis::Desc& SPxBasis::Desc::operator=(const SPxBasis::Desc& rhs)
{
   if (this != &rhs)
   {
      rowstat = rhs.rowstat;
      colstat = rhs.colstat;
      
      if (rhs.stat == &rhs.rowstat)
      {
         assert(rhs.costat == &rhs.colstat);
         
         stat   = &rowstat;
         costat = &colstat;
      }
      else
      {
         assert(rhs.costat == &rhs.rowstat);
         
         stat   = &colstat;
         costat = &rowstat;
      }
   }
   return *this;
}

void SPxBasis::Desc::reSize(int rowDim, int colDim)
{
   METHOD( "SPxBasis::Desc::reSize()" );
   rowstat.reSize(rowDim);
   colstat.reSize(colDim);
}

void SPxBasis::Desc::dump() const
{
   METHOD( "SPxBasis::Desc::dump()" );
   int i;

   // Dump regardless of the verbosity level if this method is called.
   const SPxOut::Verbosity tmp_verbosity = s_spxout.getVerbosity();
   s_spxout.setVerbosity( SPxOut::ERROR );

   s_spxout << "column status: ";
   for(i = 0; i < nCols(); i++)
      s_spxout << colStatus(i);
   s_spxout << std::endl;

   s_spxout << "row status:    ";
   for(i = 0; i < nRows(); i++)
      s_spxout << rowStatus(i); 
   s_spxout << std::endl;
   s_spxout.setVerbosity( tmp_verbosity );
}

bool SPxBasis::Desc::isConsistent() const
{
   METHOD( "SPxBasis::Desc::isConsistent()" );
   return rowstat.isConsistent() && colstat.isConsistent();
}

std::ostream& operator<<(std::ostream& os, const SPxBasis::Desc::Status& stat)
{
   char text;
   
   switch(stat)
   {
   case SPxBasis::Desc::P_ON_LOWER :
      text = 'L';
      break;
   case SPxBasis::Desc::P_ON_UPPER :
      text = 'U';
      break;
   case SPxBasis::Desc::P_FREE :
      text = 'F';
      break;
   case SPxBasis::Desc::P_FIXED :
      text = 'X';
      break;
   case SPxBasis::Desc::D_FREE :
      text = 'f';
      break;
   case SPxBasis::Desc::D_ON_UPPER :
      text = 'u';
      break;
   case SPxBasis::Desc::D_ON_LOWER :
      text = 'l';
      break;
   case SPxBasis::Desc::D_ON_BOTH :
      text = 'x';
      break;
   case SPxBasis::Desc::D_UNDEFINED :
      text = '.';
      break;
   default :
      os << std::endl << "Invalid status <" << int(stat) << ">" << std::endl;
      assert(false);
   }
   os << text;

   return os;
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
