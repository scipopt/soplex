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
#pragma ident "@(#) $Id: spxpricer.h,v 1.2 2001/11/06 23:31:04 bzfkocht Exp $"


#ifndef _SPXPRICE_H_
#define _SPXPRICE_H_

//@ ----------------------------------------------------------------------------
/*      \Section{Imports}
    Import required system include files
 */
#include <assert.h>


/*  and class header files
 */

#include "soplex.h"

namespace soplex
{






//@ ----------------------------------------------------------------------------
/* \Section{Class Declaration}
 */

/** #SoPlex# pricer base class.
    Class #SPxPricer# is a pure virtual class defining the interface for pricer
    classes to be used by #SoPlex#. The pricers task is to select a vector to
    enter or leave the simplex basis, depending on the chosen simplex type.
 
    An #SPxPricer# is first #load#ed the #SoPlex# object for which pricing is to
    be performed for. Then depending of the #SoPlex::Type#, methods
    #selectEnter# and #entered4# (for #ENTER#ing Simplex) or #selectLeave# and
    #left4# (for #LEAVE#ing Simplex) are called by #SoPlex#. The #SPxPricer#
    object is informed of a change of the #SoPlex::Type# by calling method
    #setType#.
 */
class SPxPricer
{
public:
   /** Load LP.
       Load the solver and LP for which pricing steps are to be performed.
    */
   virtual void load(SoPlex* lp) = 0;

   ///
   virtual void clear() = 0;

   /// return loaded #SoPlex# object.
   virtual SoPlex* solver() const = 0;

   /// violations up to #epsilon# are tollerated.
   virtual double epsilon() const = 0;
   ///
   virtual void setEpsilon(double eps) = 0;

   /** Set pricing type.
       Inform pricer about (a change of) the loaded #SoPlex#'s #Type#. In
       the sequel, only the corresponding select methods may be called.
    */
   virtual void setType(SoPlex::Type) = 0;
   /** Set basis representation.
       Inform pricer about (a change of) the loaded #SoPlex#'s
       #Representation#.
    */
   virtual void setRep(SoPlex::Representation) = 0;

   /** Select index to leave basis.
       Select the index of a vector to leave the basis. The selected index
       #i#, say, must be in the range #0 <= i < solver()->dim()# and its
       tested value must fullfill #solver()->test()[i] < -epsilon()#.
    */
   virtual int selectLeave() = 0;
   /** Perform leaving pivot.
       Method #left4# is called after each simplex iteration in #LEAVE#
       mode. It informs the #SPxPricer# that the #n#-th variable has left
       the basis for #id# to come in at this position. When beeing called,
       all vectors of #SoPlex# involved in such an entering update are
       setup correctly and may be accessed via the corresponding methods
       (i.e.~#fVec()#, #pVec()# etc.). In gerneral, argument #n# will be
       the one returned by the #SPxPricer# at the previous call to
       #selectLeave()#. However, one can not rely on this.
    */
   virtual void left4(int n, SoPlex::Id id) = 0;

   /** Select Id to enter basis.
       Select the #SoPlex::Id# of a vector to enter the basis. The selected
       #SoPlex::Id id#, say, must not represent a basic index (i.e.
       #solver()->isBasic(id)# must be false). However, the corresponding
       test value needs not to be less than #-epsilon()#. If not, #SoPlex#
       will discard the pivot.\\
       Note:
       When method #selectEnter()# is called by the loaded #SoPlex#
       object, all values from #coTest()# are up to date. However, whether
       the elements of #test()# are so depends on the #SoPlex::Pricing#
       type.
    */
   virtual SoPlex::Id selectEnter() = 0;
   /** Perform entering pivot.
       Method #entered4# is called after each simplex iteration in #ENTER#
       mode. It informs the #SPxPricer# that variable #id# has entered
       at the #n#-th position. When beeing called, all vectors of #SoPlex#
       involved in such an entering update are setup correctly and may be
       accessed via the corresponding methods (i.e. #fVec()#, #pVec()#
       etc.). In gerneral, argument #id# will be the one returned by the
       #SPxPricer# at the previous call to #selectEnter()#. However, one
       can not rely on this.
    */
   virtual void entered4(SoPlex::Id id, int n) = 0;



   /**@name Extension */
   //@{
   /// #n# vectors have been added to loaded LP.
   virtual void addedVecs (int n) = 0;
   /// #n# covectors have been added to loaded LP.
   virtual void addedCoVecs(int n) = 0;
   //@}


   /**@name Shrinking */
   //@{
   /// vector #i# was removed from loaded LP.
   virtual void removedVec(int i) = 0;
   /// vectors given by #perm# have been removed from loaded LP.
   virtual void removedVecs(const int perm[]) = 0;
   /// covector #i# was removed from loaded LP.
   virtual void removedCoVec(int i) = 0;
   /// covectors given by #perm# have been removed from loaded LP.
   virtual void removedCoVecs(const int perm[]) = 0;
   //@}


   /**@name Manipulation */
   //@{
   /// change objective vector.
   virtual void changeObj(const Vector& newObj) = 0;

   /// change #i#-th objective vector element.
   virtual void changeObj(int i, double newVal) = 0;

   /// change vector of lower bounds.
   virtual void changeLower(const Vector& newLower) = 0;

   /// change #i#-th lower bound.
   virtual void changeLower(int i, double newLower) = 0;

   /// change vector of upper bounds.
   virtual void changeUpper(const Vector& newUpper) = 0;

   /// change #i#-th upper bound.
   virtual void changeUpper(int i, double newUpper) = 0;

   /// change lhs vector for constraints.
   virtual void changeLhs(const Vector& newLhs) = 0;

   /// change #i#-th lhs value.
   virtual void changeLhs(int i, double newLhs) = 0;

   /// change rhs vector for constraints.
   virtual void changeRhs(const Vector& newRhs) = 0;

   /// change #i#-th rhs value.
   virtual void changeRhs(int i, double newRhs) = 0;

   /// change #i#-th row of LP.
   virtual void changeRow(int i, const LPRow& newRow) = 0;

   /// change #i#-th column of LP.
   virtual void changeCol(int i, const LPCol& newCol) = 0;

   /// change LP element (#i#, #j#).
   virtual void changeElement(int i, int j, double val) = 0;

   /// change optimization sense to #sns#.
   virtual void changeSense(SoPlex::Sense sns) = 0;
   //@}

   ///
   virtual ~SPxPricer()
   {}

};


} // namespace soplex
#endif // _SPXPRICER_H_

//-----------------------------------------------------------------------------
//Emacs Local Variables:
//Emacs mode:c++
//Emacs c-basic-offset:3
//Emacs tab-width:8
//Emacs indent-tabs-mode:nil
//Emacs End:
//-----------------------------------------------------------------------------
