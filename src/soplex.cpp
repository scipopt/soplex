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
#pragma ident "@(#) $Id: soplex.cpp,v 1.69 2003/03/03 08:30:07 bzfkocht Exp $"

#include <iostream>

#include "soplex.h"

namespace soplex
{
SoPlex::SoPlex(SPxSolver::Type p_type, SPxSolver::Representation p_rep)
   : m_solver(p_type, p_rep)
   , m_preScaler(0)
   , m_postScaler(0)
   , m_simplifier(0)
{
   m_solver.setSolver(&m_slu);  
   m_solver.setTester(&m_fastRT);
   m_solver.setPricer(&m_steepPR);
   m_solver.setStarter(0);

   assert(SoPlex::isConsistent());
}

void SoPlex::setPreScaler(SPxScaler* x)
{
   METHOD( "SoPlex::setScaler()" );

   m_preScaler = x;
}

void SoPlex::setPostScaler(SPxScaler* x)
{
   METHOD( "SoPlex::setScaler()" );

   m_postScaler = x;
}

void SoPlex::setSimplifier(SPxSimplifier* x)
{
   METHOD( "SoPlex::setSimplifier()" );

   m_simplifier = x;
}

Real SoPlex::objValue() const
{
   METHOD( "SoPlex::value()" );

   DVector x(nCols());

   getPrimal(x);

   return x * maxObj() * Real(spxSense());
}

SPxSolver::Status SoPlex::solve()
{
   METHOD( "SoPlex::solve()" );

   if (nRows() <= 0 && nCols() <= 0) // no problem loaded
      return SPxSolver::NO_PROBLEM;

   {  // context for working LP
      SPxLP work(*this);

      // should the LP be scaled
      if (m_preScaler != 0)
         m_preScaler->scale(work);

      // should the LP be simplified ?
      if (m_simplifier != 0)
      {
         switch(m_simplifier->simplify(work, m_solver.epsilon(), m_solver.delta()))
         {
         case SPxSimplifier::UNBOUNDED :
            m_solver.setBasisStatus(SPxBasis::UNBOUNDED);
            return SPxSolver::UNBOUNDED;
         case SPxSimplifier::INFEASIBLE :
            m_solver.setBasisStatus(SPxBasis::INFEASIBLE);
            return SPxSolver::INFEASIBLE;
         case SPxSimplifier::VANISHED :
            m_solver.setBasisStatus(SPxBasis::OPTIMAL);
            return SPxSolver::OPTIMAL;
         case SPxSimplifier::OKAY:
            break;
         default:
            abort();
         }
      }
      // should the LP be scaled after simplifing?
      if (m_postScaler != 0)
         m_postScaler->scale(work);

      m_solver.loadLP(work);
   }
   return m_solver.solve();
}

SPxSolver::Status SoPlex::getPrimal(Vector& x) const
{
   METHOD( "SoPlex::getPrimal()" );

   DVector psp_x(m_solver.nCols()); // prescaled simplified postscaled

   SPxSolver::Status stat = m_solver.getPrimal(psp_x);
   
   if (m_postScaler != 0)
      m_postScaler->unscalePrimal(psp_x);
   
   if (m_simplifier != 0)
      x = m_simplifier->unsimplifiedPrimal(psp_x);
   else
      x = psp_x;

   if (m_preScaler != 0)
      m_preScaler->unscalePrimal(x);

   return stat;
}

SPxSolver::Status SoPlex::getSlacks(Vector& s) const
{
   return m_solver.getSlacks(s);
}

SPxSolver::Status SoPlex::getDual(Vector& pi) const
{
   return m_solver.getDual(pi);
}
  
SPxSolver::Status SoPlex::getRedCost(Vector& rdcost) const
{
   return m_solver.getRedCost(rdcost);
}


void SoPlex::qualConstraintViolation(
   Real& maxviol, 
   Real& sumviol) const
{
   maxviol = 0.0;
   sumviol = 0.0;

   DVector solu( nCols() );

   getPrimal( solu );

   for( int row = 0; row < nRows(); ++row )
   {
      const SVector& rowvec = rowVector( row );

      Real val = 0.0;         

      for( int col = 0; col < rowvec.size(); ++col )
         val += rowvec.value( col ) * solu[rowvec.index( col )];

      Real viol = 0.0;

      assert(lhs( row ) <= rhs( row ));

      if (val < lhs( row )) 
         viol = fabs(val - lhs( row ));
      else
         if (val > rhs( row ))
            viol = fabs(val - rhs( row ));

      if (viol > maxviol)
         maxviol = viol;

      sumviol += viol;
   }
}

void SoPlex::qualBoundViolation(
   Real& maxviol, 
   Real& sumviol) const
{
   maxviol = 0.0;
   sumviol = 0.0;

   DVector solu( nCols() );

   getPrimal( solu );

   for( int col = 0; col < nCols(); ++col )
   {
      assert( lower( col ) <= upper( col ));

      Real viol = 0.0;

      if (solu[col] < lower( col ))
         viol = fabs( solu[col] - lower( col ));
      else
         if (solu[col] > upper( col ))
            viol = fabs( solu[col] - upper( col ));
         
      if (viol > maxviol)
         maxviol = viol;

      sumviol += viol;
   }
}

bool SoPlex::writeBasisFile(
   const char* filename, 
   const NameSet& rowNames, 
   const NameSet& colNames)
{
   std::cout << "Warning! Not fully implemented" << std::endl;
   return m_solver.writeBasisFile(filename, rowNames, colNames);
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







