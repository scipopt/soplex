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
#pragma ident "@(#) $Id: spxweightpr.h,v 1.4 2001/11/28 16:41:23 bzfpfend Exp $"


/**@file  spxweightpr.h
 * @brief Weighted pricing.
 */
#ifndef _SPXWEIGHTPR_H_
#define _SPXWEIGHTPR_H_

#include "spxpricer.h"

namespace soplex
{

/**@brief   Weighted pricing.
   @ingroup Algo
      
   Class #SPxWeightPR is an implemantation class of #SPxPricer that uses
   weights for columns and rows for selecting the Simplex pivots. The weights
   are computed by methods #computeCP() and #computeRP() which may be
   overridden by derived classes.
   
   The weights are interpreted as follows: The higher a value is, the more
   likely the corresponding row or column is set on one of its bounds.
   
   See #SPxPricer for a class documentation.
*/
class SPxWeightPR : public SPxPricer
{
protected:
   DVector cPenalty;               // column penalties
   DVector rPenalty;               // row    penalties
   DVector leavePenalty;           // penalties for leaveing alg

   const double* coPenalty;
   const double* penalty;

   double objlength;              // length of objective vector.

   SoPlex* thesolver;
   double theeps;

public:
   ///
   SoPlex* solver() const
   {
      return thesolver;
   }

   ///
   double epsilon() const
   {
      return theeps;
   }
   ///
   void setEpsilon(double eps)
   {
      theeps = eps;
   }

   ///
   void load(SoPlex* base);

   ///
   void clear()
   {
      thesolver = 0;
   }

   ///
   int selectLeave();

   ///
   SoPlex::Id selectEnter();


   ///
   virtual void addedVecs (int n);
   ///
   virtual void addedCoVecs(int n);


   ///
   virtual void removedVec(int i);
   ///
   virtual void removedCoVecs(const int perm[]);
   ///
   virtual void removedVecs(const int perm[]);
   ///
   virtual void removedCoVec(int i);


   ///
   void setType(SoPlex::Type tp);
   ///
   void setRep(SoPlex::Representation rep);
   ///
   void left4(int, SoPlex::Id)
   {}
   ///
   void entered4(SoPlex::Id, int)
   {}

   ///
   void changeObj(const Vector&)
   {}
   ///
   void changeObj(int, double)
   {}
   ///
   void changeLower(const Vector&)
   {}
   ///
   void changeLower(int, double)
   {}
   ///
   void changeUpper(const Vector&)
   {}
   ///
   void changeUpper(int, double)
   {}
   ///
   void changeLhs(const Vector&)
   {}
   ///
   void changeLhs(int, double)
   {}
   ///
   void changeRhs(const Vector&)
   {}
   ///
   void changeRhs(int, double)
   {}
   ///
   void changeRow(int, const LPRow&)
   {}
   ///
   void changeCol(int, const LPCol&)
   {}
   ///
   void changeElement(int, int, double)
   {}
   ///
   void changeSense(SoPlex::Sense)
   {}

   ///
   int isConsistent() const;


protected:
   /// compute weights for columns.
   virtual void computeCP(int start, int end);
   /// compute weights for rows.
   virtual void computeRP(int start, int end);
};


} // namespace soplex
#endif // _SPXWEIGHTPR_H_

//-----------------------------------------------------------------------------
//Emacs Local Variables:
//Emacs mode:c++
//Emacs c-basic-offset:3
//Emacs tab-width:8
//Emacs indent-tabs-mode:nil
//Emacs End:
//-----------------------------------------------------------------------------
