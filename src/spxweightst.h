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
#pragma ident "@(#) $Id: spxweightst.h,v 1.2 2001/11/06 23:31:06 bzfkocht Exp $"


#ifndef _SPXWEIGHTST_H_
#define _SPXWEIGHTST_H_

//@ ----------------------------------------------------------------------------
/*      \Section{Imports}
    Import required system include files
 */
#include <assert.h>


/*  and class header files
 */

#include "spxstarter.h"
#include "dataarray.h"

namespace soplex
{





//@ ----------------------------------------------------------------------------
/* \Section{Class Declaration}
 */

/** weighted start basis.
    Class #SPxWeightST# is an implementation of a #SPxStarter# for generating a
    Simplex starting basis. Using method #setupWeights()# it sets up arrays
    #weight# and #coWeight#, or equivalently #rowWeight# and #colWeight#.
    (#rowWeight# and #colWeight# are just pointers initialized to #weight# and
    #coWeight# according to the representation of \Ref{SoPlex} #base# passed to
    method #generate#.) 
 
    The weight values are then used to setup a starting basis for the LP:
    vectors with low values are likely to become dual (i.e. basic for a column
    basis) and such with high values are likely to become primal (i.e. nonbasic
    for a column basis).
    
    However, if a varialbe having an upper and lower bound is to become primal,
    there is still a choice for setting it either to its upper or lower bound.
    Members #rowRight# and #colUp# are used to determine where to set a primal
    variable. If #rowRight[i]# is set to a nonzero value, the right-hand side
    inequality is set tightly for the #i#-th to become primal. Analogously, If
    #colUp[j]# is nonzero, the #j#-th variable will be set to its upper bound
    if it becomes primal.
 */
class SPxWeightST : public SPxStarter
{
   DataArray < int > forbidden;

   DataArray < double > * weight;
   DataArray < double > * coWeight;
   void setPrimalStatus(SPxBasis::Desc&, const SoPlex&, const SoPlex::Id&);

protected:
   /// weight value for LP rows.
   DataArray < double > rowWeight;
   /// weight value for LP columns.
   DataArray < double > colWeight;

   /// set variable to rhs?.
   DataArray < int > rowRight;
   /// set primal variable to upper bound.
   DataArray < int > colUp;

   /**
       This method is called in order to setup the weights for all
       variables. It has been declared #virtual# in order to allow for
       derived classes to compute other weight values.
    */
   virtual void setupWeights(SoPlex& base);

public:
   ///
   void generate(SoPlex& base);

   ///
   SPxWeightST()
   {}
   ///

   int isConsistent() const;
};

} // namespace soplex
#endif // _SPXWEIGHTST_H_

//-----------------------------------------------------------------------------
//Emacs Local Variables:
//Emacs mode:c++
//Emacs c-basic-offset:3
//Emacs tab-width:8
//Emacs indent-tabs-mode:nil
//Emacs End:
//-----------------------------------------------------------------------------
