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
#pragma ident "@(#) $Id: spxdevexpr.h,v 1.4 2001/11/28 16:41:22 bzfpfend Exp $"


/**@file  spxdevexpr.h
 * @brief Devex pricer.
 */
#ifndef _SPXDEVEXPR_H_
#define _SPXDEVEXPR_H_


#include <assert.h>

#include "spxpricer.h"

namespace soplex
{

/**@brief   Devex pricer.
   @ingroup Algo

   The Devex Pricer for SoPlex implements an approximate steepest edge pricing,
   that does without solving an extra linear system and computing the scalar
   products.

   See #SPxPricer for a class documentation.
*/
class SPxDevexPR : public SPxPricer
{
private:
protected:
   double last;            ///< penalty, selected at last iteration.
   DVector penalty;        ///< vector of pricing penalties.
   DVector coPenalty;      ///< vector of pricing penalties.

   SoPlex* thesolver;
   double theeps;

public:
   /// returns loaded solver.
   SoPlex* solver() const
   {
      return thesolver;
   }

   /// bound violations up to #epsilon are tolerated.
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
   void setType(SoPlex::Type);

   ///
   void setRep(SoPlex::Representation rep);

   ///
   int selectLeave();
protected:
   int selectLeave(double& best, int start = 0, int incr = 1);
public:

   ///
   void left4(int n, SoPlex::Id id);
protected:
   void left4(int n, SoPlex::Id id, int start, int incr);
public:

   ///
   SoPlex::Id selectEnter();
protected:
   SoPlex::Id selectEnter(double& best, int start1 = 0, int incr1 = 1,
                          int start2 = 0, int incr2 = 1);
public:

   ///
   void entered4(SoPlex::Id id, int n);
protected:
   void entered4(SoPlex::Id id, int n, int start1, int incr1, int start2,
                 int incr2);
public:


   /// \p n vectors have been added to loaded LP.
   virtual void addedVecs (int n);
   /// \p n covectors have been added to loaded LP.
   virtual void addedCoVecs(int n);


   ///
   virtual void removedVec(int i);
   ///
   virtual void removedCoVecs(const int perm[]);
   ///
   virtual void removedVecs(const int perm[]);
   /// These  methods are use for implemementing the public remove methods.
   virtual void removedCoVec(int i);


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
};

} // namespace soplex
#endif // _SPXDEVEXPR_H_

//-----------------------------------------------------------------------------
//Emacs Local Variables:
//Emacs mode:c++
//Emacs c-basic-offset:3
//Emacs tab-width:8
//Emacs indent-tabs-mode:nil
//Emacs End:
//-----------------------------------------------------------------------------
