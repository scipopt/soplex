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
#pragma ident "@(#) $Id: spxpricer.h,v 1.3 2001/11/28 16:41:22 bzfpfend Exp $"


/**@file  spxpricer.h
 * @brief Abstract pricer base class.
 */
#ifndef _SPXPRICE_H_
#define _SPXPRICE_H_

#include <assert.h>

#include "soplex.h"


namespace soplex
{

/**@todo document the member variables of derived classes of SPxPricer */

/**@brief   Abstract pricer base class.
   @ingroup Algo

   Class #SPxPricer is a pure virtual class defining the interface for pricer
   classes to be used by #SoPlex. The pricers task is to select a vector to
   enter or leave the simplex basis, depending on the chosen simplex type.
   
   An #SPxPricer is first #load%ed the #SoPlex object for which pricing is to
   be performed for. Then depending of the #SoPlex::Type, methods
   #selectEnter() and #entered4() (for #ENTER%ing Simplex) or #selectLeave()
   and #left4() (for #LEAVE%ing Simplex) are called by #SoPlex. The #SPxPricer
   object is informed of a change of the #SoPlex::Type by calling method
   #setType.
*/
class SPxPricer
{
public:
   /**@name Initialization */
   //@{
   /// loads LP.
   /** Loads the solver and LP for which pricing steps are to be performed.
    */
   virtual void load(SoPlex* lp) = 0;

   /// unloads LP.
   virtual void clear() = 0;

   /// returns loaded #SoPlex object.
   virtual SoPlex* solver() const = 0;

   /// returns violation bound #epsilon.
   virtual double epsilon() const = 0;

   /// sets violation bound.
   /** Inequality violations are accepted, if their size is less than \p eps.
    */
   virtual void setEpsilon(double eps) = 0;

   /// sets pricing type.
   /** Informs pricer about (a change of) the loaded #SoPlex's #Type. In
       the sequel, only the corresponding select methods may be called.
    */
   virtual void setType(SoPlex::Type) = 0;

   /// sets basis representation.
   /** Informs pricer about (a change of) the loaded #SoPlex's
       #Representation.
   */
   virtual void setRep(SoPlex::Representation) = 0;
   //@}


   /**@name Pivoting */
   //@{
   /// returns selected index to leave basis.
   /** Selects the index of a vector to leave the basis. The selected index
       i, say, must be in the range 0 <= i < #solver()->dim() and its
       tested value must fullfill #solver()->test()[i] < -#epsilon().
    */
   virtual int selectLeave() = 0;

   /// performs leaving pivot.
   /** Method #left4() is called after each simplex iteration in #LEAVE
       mode. It informs the #SPxPricer that the \p n 'th variable has left
       the basis for \p id to come in at this position. When beeing called,
       all vectors of #SoPlex involved in such an entering update are
       setup correctly and may be accessed via the corresponding methods
       (i.e. #fVec(), #pVec(), etc.). In general, argument \p n will be
       the one returned by the #SPxPricer at the previous call to
       #selectLeave(). However, one can not rely on this.
    */
   virtual void left4(int n, SoPlex::Id id) = 0;

   /// selects Id to enter basis.
   /** Selects the #SoPlex::Id of a vector to enter the basis. The selected
       id, must not represent a basic index (i.e. #solver()->isBasic(id) must
       be false). However, the corresponding test value needs not to be less
       than #-epsilon(). If not, #SoPlex will discard the pivot.

       Note:
       When method #selectEnter() is called by the loaded #SoPlex
       object, all values from #coTest() are up to date. However, whether
       the elements of #test() are so depends on the #SoPlex::Pricing
       type.
    */
   virtual SoPlex::Id selectEnter() = 0;

   /// performs entering pivot.
   /** Method #entered4() is called after each simplex iteration in #ENTER
       mode. It informs the #SPxPricer that variable \p id has entered
       at the \p n 'th position. When beeing called, all vectors of #SoPlex
       involved in such an entering update are setup correctly and may be
       accessed via the corresponding methods (i.e. #fVec(), #pVec(),
       etc.). In general, argument \p id will be the one returned by the
       #SPxPricer at the previous call to #selectEnter(). However, one
       can not rely on this.
    */
   virtual void entered4(SoPlex::Id id, int n) = 0;
   //@}


   /**@name Extension */
   //@{
   /// \p n vectors have been added to loaded LP.
   virtual void addedVecs (int n) = 0;
   /// \p n covectors have been added to loaded LP.
   virtual void addedCoVecs(int n) = 0;
   //@}


   /**@name Shrinking */
   //@{
   /// vector \p i was removed from loaded LP.
   virtual void removedVec(int i) = 0;
   /// vectors given by \p perm have been removed from loaded LP.
   virtual void removedVecs(const int perm[]) = 0;
   /// covector \p i was removed from loaded LP.
   virtual void removedCoVec(int i) = 0;
   /// covectors given by \p perm have been removed from loaded LP.
   virtual void removedCoVecs(const int perm[]) = 0;
   //@}


   /**@name Manipulation */
   //@{
   /// objective vector has been changed to \p newObj.
   virtual void changeObj(const Vector& newObj) = 0;

   /// objective value of column \p i has been changed to \p newVal.
   virtual void changeObj(int i, double newVal) = 0;

   /// lower bound vector has been changed to \p newLower.
   virtual void changeLower(const Vector& newLower) = 0;

   /// lower bound of column \p i has been changed to \p newLower.
   virtual void changeLower(int i, double newLower) = 0;

   /// upper bound vector has been changed to \p newUpper.
   virtual void changeUpper(const Vector& newUpper) = 0;

   /// upper bound of column \p i has been changed to \p newUpper.
   virtual void changeUpper(int i, double newUpper) = 0;

   /// left hand side vector for constraints has been changed to \p newLhs.
   virtual void changeLhs(const Vector& newLhs) = 0;

   /// left hand side of row \p i has been changed to \p newLhs.
   virtual void changeLhs(int i, double newLhs) = 0;

   /// right hand side vector for constraints has been changed to \p newRhs.
   virtual void changeRhs(const Vector& newRhs) = 0;

   /// right hand side of row \p i has been changed to \p newRhs.
   virtual void changeRhs(int i, double newRhs) = 0;

   /// row \p i has been replaced by \p newRow.
   virtual void changeRow(int i, const LPRow& newRow) = 0;

   /// column \p i has been replaced by \p newCol.
   virtual void changeCol(int i, const LPCol& newCol) = 0;

   /// constraint matrix element (\p i, \p j) has been changed to \p val.
   virtual void changeElement(int i, int j, double val) = 0;

   /// optimization sense has been changed to \p sns.
   virtual void changeSense(SoPlex::Sense sns) = 0;
   //@}

   /**@name Constructors / Destructors */
   //@{
   /// destructor.
   virtual ~SPxPricer()
   {}
   //@}

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
