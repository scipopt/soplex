/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the class library                   */
/*       SoPlex --- the Sequential object-oriented simPlex.                  */
/*                                                                           */
/*    Copyright (C) 1996-2015 Konrad-Zuse-Zentrum                            */
/*                            fuer Informationstechnik Berlin                */
/*                                                                           */
/*  SoPlex is distributed under the terms of the ZIB Academic Licence.       */
/*                                                                           */
/*  You should have received a copy of the ZIB Academic License              */
/*  along with SoPlex; see the file COPYING. If not email to soplex@zib.de.  */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

#include <assert.h>
#include <iostream>

#include "spxdefines.h"
#include "rational.h"
#include "spxsolver.h"
#include "spxpricer.h"
#include "spxratiotester.h"
#include "spxdefaultrt.h"
#include "spxstarter.h"
#include "spxout.h"

#define MAXCYCLES 400
#define MAXSTALLS 10000
#define MAXSTALLRECOVERS 10
#define MAXREFACPIVOTS 10

namespace soplex
{

/**@todo check separately for ENTER and LEAVE algorithm */
bool SPxSolver::precisionReached(Real& newpricertol) const
{
   Real maxViolRedCost;
   Real sumViolRedCost;
   Real maxViolBounds;
   Real sumViolBounds;
   Real maxViolConst;
   Real sumViolConst;

   qualRedCostViolation(maxViolRedCost, sumViolRedCost);
   qualBoundViolation(maxViolBounds, sumViolBounds);
   qualConstraintViolation(maxViolConst, sumViolConst);

   // is the solution good enough ?
   bool reached = maxViolRedCost < opttol() && maxViolBounds < feastol() && maxViolConst < feastol();

   if (!reached)
   {
      newpricertol = thepricer->epsilon() / 10.0;

      MSG_INFO3( (*spxout), (*spxout) << "Precision not reached: Pricer tolerance = "
                           << thepricer->epsilon()
                           << " new tolerance = " << newpricertol
                           << std::endl
                           << " maxViolRedCost= " << maxViolRedCost
                           << " maxViolBounds= " << maxViolBounds
                           << " maxViolConst= " << maxViolConst
                           << std::endl
                           << " sumViolRedCost= " << sumViolRedCost
                           << " sumViolBounds= " << sumViolBounds
                           << " sumViolConst= " << sumViolConst
                           << std::endl; );
   }
   return reached;
}

SPxSolver::Status SPxSolver::solve()
{

   SPxId enterId;
   int   leaveNum;
   int   loopCount = 0;
   Real  minShift = infinity;
   int   cycleCount = 0;
   bool  priced = false;
   Real  lastDelta = 1;

   /* store the last (primal or dual) feasible objective value to recover/abort in case of stalling */
   Real  stallRefValue;
   Real  stallRefShift;
   int   stallRefIter;
   int   stallNumRecovers;

   if (dim() <= 0 && coDim() <= 0) // no problem loaded
   {
      m_status = NO_PROBLEM;
      throw SPxStatusException("XSOLVE01 No Problem loaded");
   }

   if (slinSolver() == 0) // linear system solver is required.
   {
      m_status = NO_SOLVER;
      throw SPxStatusException("XSOLVE02 No Solver loaded");
   }
   if (thepricer == 0) // pricer is required.
   {
      m_status = NO_PRICER;
      throw SPxStatusException("XSOLVE03 No Pricer loaded");
   }
   if (theratiotester == 0) // ratiotester is required.
   {
      m_status = NO_RATIOTESTER;
      throw SPxStatusException("XSOLVE04 No RatioTester loaded");
   }
   theTime->reset();
   theTime->start();

   m_numCycle = 0;
   iterCount  = 0;
   if (!isInitialized())
   {
      /*
      if(SPxBasis::status() <= NO_PROBLEM)
          SPxBasis::load(this);
       */
      /**@todo != REGULAR is not enough. Also OPTIMAL/DUAL/PRIMAL should
       * be tested and acted accordingly.
       */
      if (thestarter != 0 && status() != REGULAR)  // no basis and no starter.
         thestarter->generate(*this);              // generate start basis.

      init();

      // Inna/Tobi: init might fail, if the basis is singular
      if( !isInitialized() )
      {
         assert(SPxBasis::status() == SPxBasis::SINGULAR);
         m_status = UNKNOWN;
         return status();
      }
   }

   //setType(type());

   if (!matrixIsSetup)
      SPxBasis::load(this);

   //factorized = false;

   assert(thepricer->solver()      == this);
   assert(theratiotester->solver() == this);

   // maybe this should be done in init() ?
   thepricer->setType(type());
   theratiotester->setType(type());

   MSG_INFO3( (*spxout),
      (*spxout) << "starting value = " << value() << std::endl
             << "starting shift = " << shift() << std::endl;
   )

   if (SPxBasis::status() == SPxBasis::OPTIMAL)
      setBasisStatus(SPxBasis::REGULAR);

   m_status   = RUNNING;
   bool stop  = terminate();
   leaveCount = 0;
   enterCount = 0;
   primalCount = 0;
   boundflips = 0;
   totalboundflips = 0;

   stallNumRecovers = 0;

   /* if we run into a singular basis, we will retry from regulardesc with tighter tolerance in the ratio test */
   SPxSolver::Type tightenedtype = type();
   bool tightened = false;

   while (!stop)
   {
      const SPxBasis::Desc regulardesc = desc();

      // we need to reset these pointers to avoid unnecessary/wrong solves in leave() or enter()
      solveVector2 = 0;
      solveVector3 = 0;
      coSolveVector2 = 0;
      coSolveVector3 = 0;

      try
      {

      if (type() == ENTER)
      {
         forceRecompNonbasicValue();

         int enterCycleCount = 0;
         int enterFacPivotCount = 0;

         instableEnterVal = 0;
         instableEnterId = SPxId();
         instableEnter = false;

         stallRefIter = iteration()-1;
         stallRefShift = shift();
         stallRefValue = value();

         /* in the entering algorithm, entertol() should be maintained by the ratio test and leavetol() should be
          * reached by the pricer
          */
         Real maxpricertol = leavetol();
         Real minpricertol = 0.01 * maxpricertol;

         thepricer->setEpsilon(maxpricertol);
         priced = false;

         // to avoid shifts we restrict tolerances in the ratio test
         if( loopCount > 0 )
         {
            lastDelta = (lastDelta < entertol()) ? lastDelta : entertol();
            lastDelta *= 0.01;
            theratiotester->setDelta(lastDelta);
            assert(theratiotester->getDelta() > 0);
            MSG_DEBUG( std::cout << "decreased delta for ratiotest to: " << theratiotester->getDelta() << std::endl; )
         }
         else
         {
            lastDelta = 1;
            theratiotester->setDelta(entertol());
         }

         printDisplayLine(true);
         do
         {
            printDisplayLine();

            enterId = thepricer->selectEnter();

            if (!enterId.isValid() && instableEnterId.isValid() && lastUpdate() == 0)
            {
               /* no entering variable was found, but because of valid instableEnterId we know
                  that this is due to the scaling of the test values. Thus, we use
                  instableEnterId and SPxFastRT::selectEnter shall accept even an instable
                  leaving variable. */
               MSG_INFO3( (*spxout),
                  (*spxout) << " --- trying instable enter iteration" << std::endl;
                  )

               enterId = instableEnterId;
               instableEnter = true;
               // we also need to reset the test() or coTest() value for getEnterVals()
               assert(instableEnterVal < 0);
               if( enterId.isSPxColId() )
               {
                  int idx = number(SPxColId(enterId));
                  if( rep() == COLUMN )
                  {
                     theTest[idx] = instableEnterVal;
                     if( sparsePricingEnterCo && isInfeasibleCo[idx] == SPxPricer::NOT_VIOLATED )
                     {
                        infeasibilitiesCo.addIdx(idx);
                        isInfeasibleCo[idx] = SPxPricer::VIOLATED;
                     }
                     if( hyperPricingEnter )
                        updateViolsCo.addIdx(idx);
                  }
                  else
                  {
                     theCoTest[idx] = instableEnterVal;
                     if( sparsePricingEnter && isInfeasible[idx] == SPxPricer::NOT_VIOLATED )
                     {
                        infeasibilities.addIdx(idx);
                        isInfeasible[idx] = SPxPricer::VIOLATED;
                     }
                     if( hyperPricingEnter )
                        updateViols.addIdx(idx);
                  }
               }
               else
               {
                  int idx = number(SPxRowId(enterId));
                  if( rep() == COLUMN )
                  {
                     theCoTest[idx] = instableEnterVal;
                     if( sparsePricingEnter && isInfeasible[idx] == SPxPricer::NOT_VIOLATED )
                     {
                        infeasibilities.addIdx(idx);
                        isInfeasible[idx] = SPxPricer::VIOLATED;
                     }
                     if( hyperPricingEnter )
                        updateViols.addIdx(idx);
                  }
                  else
                  {
                     theTest[idx] = instableEnterVal;
                     if( sparsePricingEnterCo && isInfeasibleCo[idx] == SPxPricer::NOT_VIOLATED )
                     {
                        infeasibilitiesCo.addIdx(idx);
                        isInfeasibleCo[idx] = SPxPricer::VIOLATED;
                     }
                     if( hyperPricingEnter )
                        updateViolsCo.addIdx(idx);
                  }
               }
            }
            else
            {
               instableEnter = false;
            }

            if (!enterId.isValid())
            {
               // we are not infeasible and have no shift
               if (  shift() <= epsilon()
                  && ( SPxBasis::status() == SPxBasis::REGULAR 
                     || SPxBasis::status() == SPxBasis::DUAL 
                     || SPxBasis::status() == SPxBasis::PRIMAL))
               {
                  Real newpricertol = minpricertol;

                  MSG_INFO2( (*spxout), (*spxout) << " --- checking feasibility and optimality\n")
                  computeTest();
                  computeCoTest();

                  // is the solution good enough ?
                  // max three times reduced
                  if ((thepricer->epsilon() > minpricertol) && !precisionReached(newpricertol))
                  {  // no!
                     // we reduce the pricer tolerance. Note that if the pricer does not find a candiate
                     // with the reduced tolerance, we quit, regardless of the violations.
                     if (newpricertol < minpricertol)
                        newpricertol = minpricertol;

                     thepricer->setEpsilon(newpricertol);

                     MSG_INFO2( (*spxout), (*spxout) << " --- setting pricer tolerance to "
                                          << thepricer->epsilon()
                                          << std::endl; )
                  }
                  // solution seems good, no check whether we are precise enough
                  else if (lastUpdate() == 0)
                  {
                     priced = true;
                     break;
                  }
                  // We have an iterationlimit and everything looks good? Then stop!
                  // 6 is just a number picked.
                  else if (maxIters > 0 && lastUpdate() < 6)
                  {
                     priced = true;
                     break;
                  }
               }
               MSG_INFO3( (*spxout), (*spxout) << " --- solve(enter) triggers refactorization" << std::endl; )

               // if the factorization is not fresh, we better refactorize and call the pricer again; however, this can
               // create cycling, so it is performed only a limited number of times per ENTER round
               if( lastUpdate() > 0 && enterFacPivotCount < MAXREFACPIVOTS )
               {
                  factorize();

                  // if the factorization was found out to be singular, we have to quit
                  if( SPxBasis::status() < SPxBasis::REGULAR )
                  {
                     MSG_ERROR( std::cerr << "Something wrong with factorization, Basis status: " << SPxBasis::status() << std::endl; )
                     stop = true;
                     break;
                  }

                  // call pricer again
                  enterId = thepricer->selectEnter();

                  // count how often the pricer has found something only after refactorizing
                  if( enterId.isValid() )
                     enterFacPivotCount++;
               }

               if( !enterId.isValid() )
               {
                  priced = true;
                  break;
               }
            }

            /* check if we have iterations left */
            if (maxIters >= 0 && iterations() >= maxIters)
            {
               MSG_INFO2( (*spxout), (*spxout) << " --- maximum number of iterations (" << maxIters
                                 << ") reached" << std::endl; )
               m_status = ABORT_ITER;
               stop = true;
               break;
            }

            enter(enterId);
            assert((testBounds(), 1));
            thepricer->entered4(lastEntered(), lastIndex());
            stop = terminate();
            clearUpdateVecs();

            /* if a successful pivot was performed or a nonbasic variable was flipped to its other bound, we reset the
             * cycle counter
             */
            if( lastEntered().isValid() )
               enterCycleCount = 0;
            else
            {
               enterCycleCount++;
               if( enterCycleCount > MAXCYCLES )
               {
                  MSG_INFO2( (*spxout), (*spxout) << " --- abort solving due to cycling in "
                                       << "entering algorithm" << std::endl; );
                  m_status = ABORT_CYCLING;
                  stop = true;
               }
            }

            /* only if the basis has really changed, we increase the iterations counter; this is not the case when only
             * a nonbasic variable was flipped to its other bound
             */
            if( lastIndex() >= 0 )
            {
               enterCount++;
               assert(lastEntered().isValid());
            }

            /* check every MAXSTALLS iterations whether shift and objective value have not changed */
            if( (iteration() - stallRefIter) % MAXSTALLS == 0 )
            {
               if( spxAbs(value() - stallRefValue) <= epsilon() && spxAbs(shift() - stallRefShift) <= epsilon() )
               {
                  if( stallNumRecovers < MAXSTALLRECOVERS )
                  {
                     /* try to recover by unshifting/switching algorithm up to MAXSTALLRECOVERS times (just a number picked) */
                     MSG_INFO3( (*spxout), (*spxout) << " --- stalling detected - trying to recover by switching to LEAVING algorithm." << std::endl; )

                     ++stallNumRecovers;
                     break;
                  }
                  else
                  {
                     /* giving up */
                     MSG_INFO2( (*spxout), (*spxout) << " --- abort solving due to stalling in entering algorithm." << std::endl; );

                     m_status = ABORT_CYCLING;
                     stop = true;
                  }
               }
               else
               {
                  /* merely update reference values */
                  stallRefIter = iteration()-1;
                  stallRefShift = shift();
                  stallRefValue = value();
               }
            }

            //@ assert(isConsistent());
         }
         while (!stop);

         MSG_INFO3( (*spxout),
            (*spxout) << " --- enter finished. iteration: " << iteration()
                   << ", value: " << value()
                   << ", shift: " << shift()
                   << ", epsilon: " << epsilon()
                   << ", feastol: " << feastol()
                   << ", opttol: " << opttol()
                   << std::endl
                   << "ISOLVE56 stop: " << stop
                   << ", basis status: " << SPxBasis::status() << " (" << int(SPxBasis::status()) << ")"
                   << ", solver status: " << m_status << " (" << int(m_status) << ")" << std::endl;
         )

         if (!stop)
         {
            /**@todo technically it would be ok to finish already when (priced && maxinfeas + shift() <= entertol()) is
             *  satisfied; maybe at least in the case when SoPlex keeps jumping back between ENTER and LEAVE always
             *  shifting (looping), we may relax this condition here;
             *  note also that unShift may increase shift() slightly due to roundoff errors
             */
            if (shift() <= epsilon())
            {
               // factorize();
               unShift();

               Real maxinfeas = maxInfeas();

               MSG_INFO3( (*spxout),
                  (*spxout) << " --- maxInfeas: " << maxinfeas
                         << ", shift: " << shift()
                         << ", entertol: " << entertol() << std::endl;
               )

               if (priced && maxinfeas + shift() <= entertol())
               {
                  setBasisStatus(SPxBasis::OPTIMAL);
                  m_status = OPTIMAL;
                  break;
               }
            }
            setType(LEAVE);
            init();
            thepricer->setType(type());
            theratiotester->setType(type());
         }
      }
      else
      {
         assert(type() == LEAVE);

         forceRecompNonbasicValue();

         int leaveCycleCount = 0;
         int leaveFacPivotCount = 0;

         instableLeaveNum = -1;
         instableLeave = false;
         instableLeaveVal = 0;

         stallRefIter = iteration()-1;
         stallRefShift = shift();
         stallRefValue = value();

         /* in the leaving algorithm, leavetol() should be maintained by the ratio test and entertol() should be reached
          * by the pricer
          */
         Real maxpricertol = entertol();
         Real minpricertol = 0.01 * maxpricertol;

         thepricer->setEpsilon(maxpricertol);
         priced = false;

         // to avoid shifts we restrict tolerances in the ratio test
         if( loopCount > 0 )
         {
            lastDelta = (lastDelta < leavetol()) ? lastDelta : leavetol();
            lastDelta *= 0.01;
            theratiotester->setDelta(lastDelta);
            assert(theratiotester->getDelta() > 0);
            MSG_DEBUG( std::cout << "decreased delta for ratiotest to: " << theratiotester->getDelta() << std::endl; )
         }
         else
         {
            lastDelta = 1;
            theratiotester->setDelta(leavetol());
         }

         printDisplayLine(true);
         do
         {
            printDisplayLine();

            leaveNum = thepricer->selectLeave();

            if (leaveNum < 0 && instableLeaveNum >= 0 && lastUpdate() == 0)
            {
               /* no leaving variable was found, but because of instableLeaveNum >= 0 we know
                  that this is due to the scaling of theCoTest[...]. Thus, we use 
                  instableLeaveNum and SPxFastRT::selectEnter shall accept even an instable
                  entering variable. */
               MSG_INFO3( (*spxout),
                  (*spxout) << " --- trying instable leave iteration" << std::endl;
               )
            
               leaveNum = instableLeaveNum;
               instableLeave = true;
               // we also need to reset the fTest() value for getLeaveVals()
               assert(instableLeaveVal < 0);
               theCoTest[instableLeaveNum] = instableLeaveVal;

               if ( sparsePricingLeave )
               {
                  if ( isInfeasible[instableLeaveNum] == SPxPricer::NOT_VIOLATED )
                  {
                     infeasibilities.addIdx(instableLeaveNum);
                     isInfeasible[instableLeaveNum] = SPxPricer::VIOLATED;
                  }
                  if( hyperPricingLeave )
                     updateViols.addIdx(instableLeaveNum);
               }
            }
            else
            {
               instableLeave = false;
            }

            if (leaveNum < 0)
            {
               // we are not infeasible and have no shift
               if (  shift() <= epsilon()
                  && (  SPxBasis::status() == SPxBasis::REGULAR 
                     || SPxBasis::status() == SPxBasis::DUAL 
                     || SPxBasis::status() == SPxBasis::PRIMAL))
               {
                  Real newpricertol = minpricertol;

                  MSG_INFO2( (*spxout), (*spxout) << " --- checking feasibility and optimality\n")
                  computeFtest();

                  // is the solution good enough ?
                  // max three times reduced
                  if ((thepricer->epsilon() > minpricertol) && !precisionReached(newpricertol))
                  {  // no
                     // we reduce the pricer tolerance. Note that if the pricer does not find a candiate
                     // with the reduced pricer tolerance, we quit, regardless of the violations.
                     if (newpricertol < minpricertol)
                        newpricertol = minpricertol;

                     thepricer->setEpsilon(newpricertol);

                     MSG_INFO2( (*spxout), (*spxout) << " --- setting pricer tolerance to "
                                          << thepricer->epsilon()
                                          << std::endl; );
                  }
                  // solution seems good, no check whether we are precise enough
                  else if (lastUpdate() == 0)
                  {
                     priced = true;
                     break;
                  }
                  // We have an iteration limit and everything looks good? Then stop!
                  // 6 is just a number picked.
                  else if (maxIters > 0 && lastUpdate() < 6)
                  {
                     priced = true;
                     break;
                  }
               }
               MSG_INFO3( (*spxout), (*spxout) << " --- solve(leave) triggers refactorization" << std::endl; )

               // if the factorization is not fresh, we better refactorize and call the pricer again; however, this can
               // create cycling, so it is performed only a limited number of times per LEAVE round
               if( lastUpdate() > 0 && leaveFacPivotCount < MAXREFACPIVOTS )
               {
                  factorize();

                  // Inna/Tobi: if the factorization was found out to be singular, we have to quit
                  if (SPxBasis::status() < SPxBasis::REGULAR)
                  {
                     MSG_ERROR( std::cerr << "Something wrong with factorization, Basis status: " << SPxBasis::status() << std::endl; )
                     stop = true;
                     break;
                  }

                  // call pricer again
                  leaveNum = thepricer->selectLeave();

                  // count how often the pricer has found something only after refactorizing
                  if( leaveNum >= 0 )
                     leaveFacPivotCount++;
               }

               if (leaveNum < 0)
               {
                  priced = true;
                  break;
               }
            }

            /* check if we have iterations left */
            if (maxIters >= 0 && iterations() >= maxIters)
            {
               MSG_INFO2( (*spxout), (*spxout) << " --- maximum number of iterations (" << maxIters
                                 << ") reached" << std::endl; )
               m_status = ABORT_ITER;
               stop = true;
               break;
            }

            leave(leaveNum);
            assert((testBounds(), 1));
            thepricer->left4(lastIndex(), lastLeft());
            stop = terminate();
            clearUpdateVecs();

            /* if a successful pivot was performed or a nonbasic variable was flipped to its other bound, we reset the
             * cycle counter
             */
            if( lastIndex() >= 0 )
               leaveCycleCount = 0;
            else
            {
               leaveCycleCount++;
               if( leaveCycleCount > MAXCYCLES )
               {
                  MSG_INFO2( (*spxout), (*spxout) << " --- abort solving due to cycling in leaving algorithm" << std::endl; );
                  m_status = ABORT_CYCLING;
                  stop = true;
               }
            }

            /* only if the basis has really changed, we increase the iterations counter; this is not the case when only
             * a nonbasic variable was flipped to its other bound
             */
            if( lastEntered().isValid() )
            {
               leaveCount++;
               assert(lastIndex() >= 0);
            }

            /* check every MAXSTALLS iterations whether shift and objective value have not changed */
            if( (iteration() - stallRefIter) % MAXSTALLS == 0 )
            {
               if( spxAbs(value() - stallRefValue) <= epsilon() && spxAbs(shift() - stallRefShift) <= epsilon() )
               {
                  if( stallNumRecovers < MAXSTALLRECOVERS )
                  {
                     /* try to recover by switching algorithm up to MAXSTALLRECOVERS times */
                     MSG_INFO3( (*spxout), (*spxout) << " --- stalling detected - trying to recover by switching to ENTERING algorithm." << std::endl; )

                     ++stallNumRecovers;
                     break;
                  }
                  else
                  {
                     /* giving up */
                     MSG_INFO2( (*spxout), (*spxout) << " --- abort solving due to stalling in leaving algorithm" << std::endl; );

                     m_status = ABORT_CYCLING;
                     stop = true;
                  }
               }
               else
               {
                  /* merely update reference values */
                  stallRefIter = iteration()-1;
                  stallRefShift = shift();
                  stallRefValue = value();
               }
            }

            //@ assert(isConsistent());
         }
         while (!stop);

         MSG_INFO3( (*spxout),
            (*spxout) << " --- leave finished. iteration: " << iteration()
                   << ", value: " << value()
                   << ", shift: " << shift()
                   << ", epsilon: " << epsilon()
                   << ", feastol: " << feastol()
                   << ", opttol: " << opttol()
                   << std::endl
                   << "ISOLVE57 stop: " << stop
                   << ", basis status: " << SPxBasis::status() << " (" << int(SPxBasis::status()) << ")"
                   << ", solver status: " << m_status << " (" << int(m_status) << ")" << std::endl;
         )

         if (!stop)
         {
            if( shift() < minShift )
            {
               minShift = shift();
               cycleCount = 0;
            }
            else
            {
               cycleCount++;
               if( cycleCount > MAXCYCLES )
               {
                  m_status = ABORT_CYCLING;
                  throw SPxStatusException("XSOLVE13 Abort solving due to cycling");
               }
               MSG_INFO3( (*spxout),
                  (*spxout) << " --- maxInfeas: " << maxInfeas()
                         << ", shift: " << shift()
                         << ", leavetol: " << leavetol()
                         << ", cycle count: " << cycleCount << std::endl;
               )
            }

            /**@todo technically it would be ok to finish already when (priced && maxinfeas + shift() <= entertol()) is
             *  satisfied; maybe at least in the case when SoPlex keeps jumping back between ENTER and LEAVE always
             *  shifting (looping), we may relax this condition here;
             *  note also that unShift may increase shift() slightly due to roundoff errors
             */
            if (shift() <= epsilon())
            {
               cycleCount = 0;
               // factorize();
               unShift();

               Real maxinfeas = maxInfeas();

               MSG_INFO3( (*spxout),
                  (*spxout) << " --- maxInfeas: " << maxinfeas
                         << ", shift: " << shift()
                         << ", leavetol: " << leavetol() << std::endl;
               )

               // We stop if we are indeed optimal, or if we have already been
               // two times at this place. In this case it seems futile to
               // continue.
               if (loopCount > 2)
               {
                  m_status = ABORT_CYCLING;
                  throw SPxStatusException("XSOLVE14 Abort solving due to looping");
               }
               else if (priced && maxinfeas + shift() <= leavetol())
               {
                  setBasisStatus(SPxBasis::OPTIMAL);
                  m_status = OPTIMAL;
                  break;
               }
               loopCount++;
            }
            setType(ENTER);
            init();
            thepricer->setType(type());
            theratiotester->setType(type());
         }
      }
      assert(m_status != SINGULAR);

      }
      catch( const SPxException& E )
      {
         // if we stopped due to a singular basis, we reload the original basis and try again with tighter
         // tolerance (only once)
         if (m_status == SINGULAR && !tightened)
         {
            tightenedtype = type();

            if( tightenedtype == ENTER )
            {
               m_entertol = 0.01 * m_entertol;

               MSG_INFO2( (*spxout), (*spxout) << " --- basis singular: reloading basis and solving with tighter ratio test tolerance " << m_entertol << std::endl; )
            }
            else
            {
               m_leavetol = 0.01 * m_leavetol;

               MSG_INFO2( (*spxout), (*spxout) << " --- basis singular: reloading basis and solving with tighter ratio test tolerance " << m_leavetol << std::endl; )
            }

            // load original basis
            int niters = iterations();
            loadBasis(regulardesc);

            // remember iteration count
            iterCount = niters;

            // try initializing basis (might fail if starting basis was already singular)
            try
            {
               init();
               theratiotester->setType(type());
            }
            catch( const SPxException& Ex )
            {
               MSG_INFO2( (*spxout), (*spxout) << " --- reloaded basis singular, resetting original tolerances" << std::endl; )

               if( tightenedtype == ENTER )
                  m_entertol = 100.0 * m_entertol;
               else
                  m_leavetol = 100.0 * m_leavetol;

               theratiotester->setType(type());

               throw Ex;
            }

            // reset status and counters
            m_status = RUNNING;
            m_numCycle = 0;
            leaveCount = 0;
            enterCount = 0;
            stallNumRecovers = 0;

            // continue
            stop = false;
            tightened = true;
         }
         // reset tolerance to its original value and pass on the exception
         else if (tightened)
         {
            if( tightenedtype == ENTER )
               m_entertol = 100.0 * m_entertol;
            else
               m_leavetol = 100.0 * m_leavetol;

            theratiotester->setType(type());

            throw E;
         }
         // pass on the exception
         else
            throw E;
      }
   }

   // reset tolerance to its original value
   if (tightened)
   {
      if( tightenedtype == ENTER )
         m_entertol = 100.0 * m_entertol;
      else
         m_leavetol = 100.0 * m_leavetol;

      theratiotester->setType(type());
   }

   theTime->stop();
   theCumulativeTime += time();

   if (m_status == RUNNING)
   {
      m_status = ERROR;
      throw SPxStatusException("XSOLVE05 Status is still RUNNING when it shouldn't be");
   }

   MSG_INFO3( (*spxout),
      (*spxout) << "Finished solving (status=" << status()
             << ", iters=" << iterCount
             << ", leave=" << leaveCount
             << ", enter=" << enterCount
             << ", flips=" << totalboundflips;
      if( status() == OPTIMAL )
         (*spxout) << ", objValue=" << value();
      (*spxout) << ")" << std::endl;
   )

#ifdef ENABLE_ADDITIONAL_CHECKS
   /* check if solution is really feasible */
   if( status() == OPTIMAL )
   {
      int     c;
      Real    val;
      DVector sol( nCols() );

      getPrimal( sol );

      for(int row = 0; row < nRows(); ++row )
      {
         const SVector& rowvec = rowVector( row );
         val = 0.0;         
         for( c = 0; c < rowvec.size(); ++c )
            val += rowvec.value( c ) * sol[rowvec.index( c )];

         if( LT( val, lhs( row ), feastol() ) ||
             GT( val, rhs( row ), feastol() ) )
         {
            // Minor rhs violations happen frequently, so print these
            // warnings only with verbose level INFO2 and higher.
            MSG_INFO2( (*spxout), (*spxout) << "WSOLVE88 Warning! Constraint " << row
                              << " is violated by solution" << std::endl
                              << "   lhs:" << lhs( row )
                              << " <= val:" << val
                              << " <= rhs:" << rhs( row ) << std::endl; )

            if( type() == LEAVE && isRowBasic( row ) )
            {
               // find basis variable
               for( c = 0; c < nRows(); ++c )
                  if (basis().baseId(c).isSPxRowId()     
                     && (number(basis().baseId(c)) == row))
                     break;

               assert( c < nRows() );

               MSG_WARNING( (*spxout), (*spxout) << "WSOLVE90 basis idx:" << c
                                   << " fVec:" << fVec()[c]
                                   << " fRhs:" << fRhs()[c]
                                   << " fTest:" << fTest()[c] << std::endl; )
            }
         }
      }
      for(int col = 0; col < nCols(); ++col )
      {
         if( LT( sol[col], lower( col ), feastol() ) ||
             GT( sol[col], upper( col ), feastol() ) )
         {
            // Minor bound violations happen frequently, so print these
            // warnings only with verbose level INFO2 and higher.
            MSG_INFO2( (*spxout), (*spxout) << "WSOLVE91 Warning! Bound for column " << col
                                 << " is violated by solution" << std::endl
                                 << "   lower:" << lower( col )
                                 << " <= val:" << sol[col]
                                 << " <= upper:" << upper( col ) << std::endl; )

            if( type() == LEAVE && isColBasic( col ) )
            {
               for( c = 0; c < nRows() ; ++c)
                  if ( basis().baseId( c ).isSPxColId()    
                     && ( number( basis().baseId( c ) ) == col ))
                     break;

               assert( c < nRows() );
               MSG_WARNING( (*spxout), (*spxout) << "WSOLVE92 basis idx:" << c
                                   << " fVec:" << fVec()[c]
                                   << " fRhs:" << fRhs()[c]
                                   << " fTest:" << fTest()[c] << std::endl; )
            }
         }
      }
   }
#endif  // ENABLE_ADDITIONAL_CHECKS

   primalCount = ( rep() == SPxSolver::COLUMN )
     ? enterCount
     : leaveCount;

   printDisplayLine(true);
   return status();
}

void SPxSolver::testVecs()
{

   assert(SPxBasis::status() > SPxBasis::SINGULAR);

   DVector tmp(dim());

   tmp = *theCoPvec;
   multWithBase(tmp);
   tmp -= *theCoPrhs;
   if (tmp.length() > leavetol())
   {
      MSG_INFO3( (*spxout), (*spxout) << "ISOLVE93 " << iteration() << ":\tcoP error = \t"
                        << tmp.length() << std::endl; )

      tmp.clear();
      SPxBasis::coSolve(tmp, *theCoPrhs);
      multWithBase(tmp);
      tmp -= *theCoPrhs;
      MSG_INFO3( (*spxout), (*spxout) << "ISOLVE94\t\t" << tmp.length() << std::endl; )

      tmp.clear();
      SPxBasis::coSolve(tmp, *theCoPrhs);
      tmp -= *theCoPvec;
      MSG_INFO3( (*spxout), (*spxout) << "ISOLVE95\t\t" << tmp.length() << std::endl; )
   }

   tmp = *theFvec;
   multBaseWith(tmp);
   tmp -= *theFrhs;
   if (tmp.length() > entertol())
   {
      MSG_INFO3( (*spxout), (*spxout) << "ISOLVE96 " << iteration() << ":\t  F error = \t"
                           << tmp.length() << std::endl; )

      tmp.clear();
      SPxBasis::solve(tmp, *theFrhs);
      tmp -= *theFvec;
      MSG_INFO3( (*spxout), (*spxout) << "ISOLVE97\t\t" << tmp.length() << std::endl; )
   }

   if (type() == ENTER)
   {
      for (int i = 0; i < dim(); ++i)
      {
         if (theCoTest[i] < -leavetol() && isCoBasic(i))
         {
            /// @todo Error message "this shalt not be": shalt this be an assert (also below)?
            MSG_ERROR( std::cerr << "ESOLVE98 testVecs: theCoTest: this shalt not be!"
                              << std::endl
                              << "  i=" << i
                              << ", theCoTest[i]=" << theCoTest[i]
                              << ", leavetol()=" << leavetol() << std::endl; )
         }
      }

      for (int i = 0; i < coDim(); ++i)
      {
         if (theTest[i] < -leavetol() && isBasic(i))
         {
            MSG_ERROR( std::cerr << "ESOLVE99 testVecs: theTest: this shalt not be!"
                              << std::endl
                              << "  i=" << i
                              << ", theTest[i]=" << theTest[i]
                              << ", leavetol()=" << leavetol() << std::endl; )
         }
      }
   }
}


/// print display line of flying table
void SPxSolver::printDisplayLine(const bool force)
{
   MSG_INFO1( (*spxout),
      if( displayLine % (displayFreq*30) == 0 )
      {
         (*spxout) << "type |   time |   iters | facts |  shift   |    value\n";
      }
      if( force || (displayLine % displayFreq == 0) )
      {
         (type() == LEAVE) ? (*spxout) << "  L  |" : (*spxout) << "  E  |";
         (*spxout) << std::fixed << std::setw(7) << std::setprecision(1) << time() << " |";
         (*spxout) << std::scientific << std::setprecision(2);
         (*spxout) << std::setw(8) << iteration() << " | "
         << std::setw(5) << slinSolver()->getFactorCount() << " | "
         << shift() << " | "
         << std::setprecision(8) << value() + objOffset()
         << std::endl;
      }
      displayLine++;
   );
}


bool SPxSolver::terminate()
{
#ifdef ENABLE_ADDITIONAL_CHECKS
   if (SPxBasis::status() > SPxBasis::SINGULAR)
      testVecs();
#endif

   int redo = dim();
   if (redo < 1000)
      redo = 1000;

   if (iteration() > 10 && iteration() % redo == 0)
   {
#ifdef ENABLE_ADDITIONAL_CHECKS
      DVector cr(*theCoPrhs);
      DVector fr(*theFrhs);
#endif 

      if (type() == ENTER)
         computeEnterCoPrhs();
      else
         computeLeaveCoPrhs();

      computeFrhs();

#ifdef ENABLE_ADDITIONAL_CHECKS
      cr -= *theCoPrhs;
      fr -= *theFrhs;
      if (cr.length() > leavetol())
         MSG_WARNING( (*spxout), (*spxout) << "WSOLVE50 unexpected change of coPrhs "
                             << cr.length() << std::endl; )
      if (fr.length() > entertol())
         MSG_WARNING( (*spxout), (*spxout) << "WSOLVE51 unexpected change of   Frhs "
                             << fr.length() << std::endl; )
#endif

      if (updateCount > 1)
      {
         MSG_INFO3( (*spxout), (*spxout) << " --- terminate triggers refactorization"
                           << std::endl; )
         factorize();
      }
      SPxBasis::coSolve(*theCoPvec, *theCoPrhs);
      SPxBasis::solve (*theFvec, *theFrhs);

      if (pricing() == FULL)
      {
         computePvec();
         if (type() == ENTER)
            computeTest();
      }

      if (shift() > 0.0)
         unShift();
   }

   // check time limit and objective limit only for non-terminal bases
   if( SPxBasis::status() >= SPxBasis::OPTIMAL  ||
       SPxBasis::status() <= SPxBasis::SINGULAR )
   {
      m_status = UNKNOWN;
      return true;
   }

   if ( isTimeLimitReached() )
   {
      MSG_INFO2( (*spxout), (*spxout) << " --- timelimit (" << maxTime
                        << ") reached" << std::endl; )
      m_status = ABORT_TIME;
      return true;
   }

   // objLimit is set and we are running DUAL:
   // - objLimit is set if objLimit < infinity
   // - DUAL is running if rep() * type() > 0 == DUAL (-1 == PRIMAL)
   //
   // In this case we have given a objective value limit, e.g, through a
   // MIP solver, and we want stop solving the LP if we figure out that the
   // optimal value of the current LP can not be better then this objective
   // limit. More precisely:
   // - MINIMIZATION Problem
   //   We want stop the solving process if
   //   objLimit <= current objective value of the DUAL LP
   // - MAXIMIZATION Problem
   //   We want stop the solving process if 
   //   objLimit >= current objective value of the DUAL LP
   if (objLimit < infinity && type() * rep() > 0)
   {
      // We have no bound shifts; therefore, we can trust the current
      // objective value.
      // It might be even possible to use this termination value in case of
      // bound violations (shifting) but in this case it is quite difficult
      // to determine if we already reached the limit.
      if( shift() < epsilon() && maxInfeas() + shift() <= opttol() )
      {
         // SPxSense::MINIMIZE == -1, so we have sign = 1 on minimizing
         if( spxSense() * value() <= spxSense() * objLimit ) 
         {
            MSG_INFO2( (*spxout), (*spxout) << " --- objective value limit (" << objLimit
               << ") reached" << std::endl; )
            MSG_DEBUG(
               (*spxout) << " --- objective value limit reached" << std::endl
                      << " (value: " << value()
                      << ", limit: " << objLimit << ")" << std::endl
                      << " (spxSense: " << int(spxSense())
                      << ", rep: " << int(rep())
                      << ", type: " << int(type()) << ")" << std::endl;
            )
            
            m_status = ABORT_VALUE;
            return true;
         }
      }
   }

   return false;
}

SPxSolver::Status SPxSolver::getPrimal (Vector& p_vector) const
{

   if (!isInitialized())
   {
      /* exit if presolving/simplifier cleared the problem */
      if (status() == NO_PROBLEM)
         return status();
      throw SPxStatusException("XSOLVE06 Not Initialized");
   }
   if (rep() == ROW)
      p_vector = coPvec();
   else
   {
      const SPxBasis::Desc& ds = desc();

      for (int i = 0; i < nCols(); ++i)
      {
         switch (ds.colStatus(i))
         {
         case SPxBasis::Desc::P_ON_LOWER :
            p_vector[i] = SPxLP::lower(i);
            break;
         case SPxBasis::Desc::P_ON_UPPER :
         case SPxBasis::Desc::P_FIXED :
            p_vector[i] = SPxLP::upper(i);
            break;
         case SPxBasis::Desc::P_FREE :
            p_vector[i] = 0;
            break;
         case SPxBasis::Desc::D_FREE :
         case SPxBasis::Desc::D_ON_UPPER :
         case SPxBasis::Desc::D_ON_LOWER :
         case SPxBasis::Desc::D_ON_BOTH :
         case SPxBasis::Desc::D_UNDEFINED :
            break;
         default:
            throw SPxInternalCodeException("XSOLVE07 This should never happen.");
         }
      }
      for (int j = 0; j < dim(); ++j)
      {
         if (baseId(j).isSPxColId())
            p_vector[ number(SPxColId(baseId(j))) ] = fVec()[j];
      }
   }
   return status();
}

SPxSolver::Status SPxSolver::getDual (Vector& p_vector) const
{

   assert(isInitialized());

   if (!isInitialized()) 
   {
      /* exit if presolving/simplifier cleared the problem */
      if (status() == NO_PROBLEM)
         return status();
      throw SPxStatusException("XSOLVE08 No Problem loaded");
   }

   if (rep() == ROW)
   {
      int i;
      p_vector = maxRowObj();
      for (i = nCols() - 1; i >= 0; --i)
      {
         if (baseId(i).isSPxRowId())
            p_vector[ number(SPxRowId(baseId(i))) ] = fVec()[i];
      }
   }
   else
      p_vector = coPvec();

   p_vector *= Real(spxSense());

   return status();
}

SPxSolver::Status SPxSolver::getRedCost (Vector& p_vector) const
{

   assert(isInitialized());

   if (!isInitialized())
   {
      throw SPxStatusException("XSOLVE09 No Problem loaded");    
      // return NOT_INIT;
   }

   if (rep() == ROW)
   {
      int i;
      p_vector.clear();
      if (spxSense() == SPxLP::MINIMIZE)
      {
         for (i = dim() - 1; i >= 0; --i)
         {
            if (baseId(i).isSPxColId())
               p_vector[ number(SPxColId(baseId(i))) ] = -fVec()[i];
         }
      }
      else
      {
         for (i = dim() - 1; i >= 0; --i)
         {
            if (baseId(i).isSPxColId())
               p_vector[ number(SPxColId(baseId(i))) ] = fVec()[i];
         }
      }
   }
   else
   {
      p_vector = maxObj();
      p_vector -= pVec();
      if (spxSense() == SPxLP::MINIMIZE)
         p_vector *= -1.0;
   }

   return status();
}

SPxSolver::Status SPxSolver::getPrimalray (Vector& p_vector) const
{

   assert(isInitialized());

   if (!isInitialized())
   {
      throw SPxStatusException("XSOLVE10 No Problem loaded");
      // return NOT_INIT;
   }

   assert(SPxBasis::status() == SPxBasis::UNBOUNDED);
   p_vector.clear();
   p_vector = primalRay;

   return status();
}

SPxSolver::Status SPxSolver::getDualfarkas (Vector& p_vector) const
{

   assert(isInitialized());

   if (!isInitialized())
   {
      throw SPxStatusException("XSOLVE10 No Problem loaded");
      // return NOT_INIT;
   }

   assert(SPxBasis::status() == SPxBasis::INFEASIBLE);
   p_vector.clear();
   p_vector = dualFarkas;

   return status();
}

SPxSolver::Status SPxSolver::getSlacks (Vector& p_vector) const
{

   assert(isInitialized());

   if (!isInitialized())
   {
      throw SPxStatusException("XSOLVE11 No Problem loaded");
      // return NOT_INIT;
   }

   if (rep() == COLUMN)
   {
      int i;
      const SPxBasis::Desc& ds = desc();
      for (i = nRows() - 1; i >= 0; --i)
      {
         switch (ds.rowStatus(i))
         {
         case SPxBasis::Desc::P_ON_LOWER :
            p_vector[i] = lhs(i);
            break;
         case SPxBasis::Desc::P_ON_UPPER :
         case SPxBasis::Desc::P_FIXED :
            p_vector[i] = rhs(i);
            break;
         case SPxBasis::Desc::P_FREE :
            p_vector[i] = 0;
            break;
         case SPxBasis::Desc::D_FREE :
         case SPxBasis::Desc::D_ON_UPPER :
         case SPxBasis::Desc::D_ON_LOWER :
         case SPxBasis::Desc::D_ON_BOTH :
         case SPxBasis::Desc::D_UNDEFINED :
            break;
         default:
            throw SPxInternalCodeException("XSOLVE12 This should never happen.");
         }
      }
      for (i = dim() - 1; i >= 0; --i)
      {
         if (baseId(i).isSPxRowId())
            p_vector[ number(SPxRowId(baseId(i))) ] = -(*theFvec)[i];
      }
   }
   else
      p_vector = pVec();

   return status();
}

void SPxSolver::setPrimal(Vector& p_vector)
{

   if (!isInitialized())
   {
      throw SPxStatusException("XSOLVE20 Not Initialized");
   }

   if (rep() == ROW)
      coPvec() = p_vector;
   else
   {
      for (int j = 0; j < dim(); ++j)
      {
         if (baseId(j).isSPxColId())
            fVec()[j] = p_vector[ number(SPxColId(baseId(j))) ];
      }
   }
}

void SPxSolver::setDual(Vector& p_vector)
{

   assert(isInitialized());

   if (!isInitialized())
   {
      throw SPxStatusException("XSOLVE21 Not Initialized");
   }

   if (rep() == ROW)
   {
      for (int i = nCols() - 1; i >= 0; --i)
      {
         if (baseId(i).isSPxRowId())
         {
            if (spxSense() == SPxLP::MAXIMIZE)
               fVec()[i] = p_vector[ number(SPxRowId(baseId(i))) ];
            else
               fVec()[i] = -p_vector[ number(SPxRowId(baseId(i))) ];
         }
      }
   }
   else
   {
      coPvec() = p_vector;
      if (spxSense() == SPxLP::MINIMIZE)
         coPvec() *= -1.0;
   }
}

void SPxSolver::setRedCost(Vector& p_vector)
{

   assert(isInitialized());

   if (!isInitialized())
   {
      throw SPxStatusException("XSOLVE22 Not Initialized");
   }

   if (rep() == ROW)
   {
      for (int i = dim() - 1; i >= 0; --i)
      {
         if (baseId(i).isSPxColId())
         {
            if (spxSense() == SPxLP::MINIMIZE)
               fVec()[i] = -p_vector[ number(SPxColId(baseId(i))) ];
            else
               fVec()[i] = p_vector[ number(SPxColId(baseId(i))) ];
         }
      }
   }
   else
   {
      pVec() = maxObj();

      if (spxSense() == SPxLP::MINIMIZE)
         pVec() += p_vector;
      else
         pVec() -= p_vector;
   }
}

void SPxSolver::setSlacks(Vector& p_vector)
{

   assert(isInitialized());

   if (!isInitialized())
   {
      throw SPxStatusException("XSOLVE23 Not Initialized");
   }

   if (rep() == COLUMN)
   {
      for (int i = dim() - 1; i >= 0; --i)
      {
         if (baseId(i).isSPxRowId())
            (*theFvec)[i] = -p_vector[ number(SPxRowId(baseId(i))) ];
      }
   }
   else
      pVec() = p_vector;
}

SPxSolver::Status SPxSolver::status() const
{
   switch( m_status )
   {
   case UNKNOWN :      
      switch (SPxBasis::status())
      {
      case SPxBasis::NO_PROBLEM :
         return NO_PROBLEM;
      case SPxBasis::SINGULAR :
         return SINGULAR;
      case SPxBasis::REGULAR :
      case SPxBasis::DUAL :
      case SPxBasis::PRIMAL :
         return UNKNOWN;
      case SPxBasis::OPTIMAL :
         return OPTIMAL;
      case SPxBasis::UNBOUNDED :
         return UNBOUNDED;
      case SPxBasis::INFEASIBLE :
         return INFEASIBLE;
      default:
         return ERROR;
      }
   case SINGULAR : 
      return m_status;
   case OPTIMAL :
      assert( SPxBasis::status() == SPxBasis::OPTIMAL );
      /*lint -fallthrough*/
   case ABORT_CYCLING :
   case ABORT_TIME :
   case ABORT_ITER :
   case ABORT_VALUE :
   case RUNNING :
   case REGULAR :
   case NOT_INIT :
   case NO_SOLVER :
   case NO_PRICER :
   case NO_RATIOTESTER :
   case ERROR:
      return m_status;
   default:
      return ERROR;
   }
}

SPxSolver::Status SPxSolver::getResult(
   Real* p_value,
   Vector* p_primal,
   Vector* p_slacks,
   Vector* p_dual,
   Vector* reduCosts)
{
   if (p_value)
      *p_value = this->value();
   if (p_primal)
      getPrimal(*p_primal);
   if (p_slacks)
      getSlacks(*p_slacks);
   if (p_dual)
      getDual(*p_dual);
   if (reduCosts)
      getRedCost(*reduCosts);
   return status();
}
} // namespace soplex
