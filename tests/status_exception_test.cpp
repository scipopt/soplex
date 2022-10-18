/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the class library                   */
/*       SoPlex --- the Sequential object-oriented simPlex.                  */
/*                                                                           */
/*    Copyright (C) 1996-2022 Konrad-Zuse-Zentrum                            */
/*                            fuer Informationstechnik Berlin                */
/*                                                                           */
/*  SoPlex is distributed under the terms of the ZIB Academic Licence.       */
/*                                                                           */
/*  You should have received a copy of the ZIB Academic License              */
/*  along with SoPlex; see the file COPYING. If not email to soplex@zib.de.  */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

#include <iostream>
#include <string>

#include "spxlp.h"
#include "spxdefines.h"
#include "soplex.h"
#include "spxsolver.h"
#include "array.h"
#include "timer.h"
#include "spxpricer.h"
#include "spxdantzigpr.h"
#include "spxparmultpr.h"
#include "spxdevexpr.h"
#include "spxhybridpr.h"
#include "spxsteeppr.h"
#include "spxweightpr.h"
#include "spxratiotester.h"
#include "spxharrisrt.h"
#include "spxdefaultrt.h"
#include "spxfastrt.h"
#include "spxsimplifier.h"
#include "spxmainsm.h"
#include "spxscaler.h"
#include "spxequilisc.h"
#include "spxgeometsc.h"
#include "spxsumst.h"
#include "spxweightst.h"
#include "spxvectorst.h"
#include "slufactor.h"
#include "spxout.h"

using namespace soplex;

class StatusExceptionCheck
{
   bool checkXSOLVR01();
   bool checkXSOLVE01();
   bool checkXSOLVE08();
   bool checkXSOLVE09();
   bool checkXSOLVE10();
   bool checkXSOLVE11();
public:
   bool noProblemCheck();
   bool noSolverCheck();
   bool noPricerCheck();
   bool noRatioTesterCheck();
   bool notInitialized();
};

bool StatusExceptionCheck::noProblemCheck()
{
   std::cout<<std::endl<<"Test for No Problem:"<<std::endl;
   bool success = true;
   std::cout<<" XSOLVR01: ";
   if(!checkXSOLVR01())
   {
      std::cout<<"failed"<<std::endl;
      success =  false;
   }else
      std::cout<<"ok"<<std::endl;
   std::cout<<" XSOLVE01: ";
   if(!checkXSOLVE01())
   {
      std::cout<<"failed"<<std::endl;
      success =  false;
   }else
      std::cout<<"ok"<<std::endl;
   std::cout<<" XSOLVE08: ";
   if(!checkXSOLVE08())
   {
      std::cout<<"failed"<<std::endl;
      success =  false;
   }else
      std::cout<<"ok"<<std::endl;
   std::cout<<" XSOLVE09: ";
   if(!checkXSOLVE09())
   {
      std::cout<<"failed"<<std::endl;
      success =  false;
   }else
      std::cout<<"ok"<<std::endl;
   std::cout<<" XSOLVE10: ";
   if(!checkXSOLVE10())
   {
      std::cout<<"failed"<<std::endl;
      success =  false;
   }else
      std::cout<<"ok"<<std::endl;
   std::cout<<" XSOLVE11: ";
   if(!checkXSOLVE11())
   {
      std::cout<<"failed"<<std::endl;
      success =  false;
   }else
      std::cout<<"ok"<<std::endl;
   return success;
}

bool StatusExceptionCheck::checkXSOLVR01()
{
   SoPlex* solver = new SoPlex();
   try{
      solver->optimize();
   }catch(SPxStatusException& x)
   {
      delete solver;
      if(x.what().find("XSOLVR01") == 0)
         return true;
      return false;
   }
   delete solver;
   return false;
}

bool StatusExceptionCheck::checkXSOLVE01()
{
   SPxSolver* solver = new SPxSolver();
   try{
      solver->solve();
   }catch(SPxStatusException& x)
   {
      delete solver;
      if(x.what().find("XSOLVE01") == 0)
         return true;
      return false;
   }
   delete solver;
   return false;
}

bool StatusExceptionCheck::checkXSOLVE08()
{
   SPxSolver* solver = new SPxSolver();
   Vector* p_primal = 0;
   try{
      solver->getDualSol(*p_primal);
   }catch(SPxStatusException& x)
   {
      delete solver;
      if(x.what().find("XSOLVE08") == 0)
         return true;
      return false;
   }
   delete solver;
   return false;
}

bool StatusExceptionCheck::checkXSOLVE09()
{
   SPxSolver* solver = new SPxSolver();
   Vector* p_primal = 0;
   try{
      solver->getRedCostSol(*p_primal);
   }catch(SPxStatusException& x)
   {
      delete solver;
      if(x.what().find("XSOLVE09") == 0)
         return true;
      return false;
   }
   delete solver;
   return false;
}

bool StatusExceptionCheck::checkXSOLVE10()
{
   SPxSolver* solver = new SPxSolver();
   Vector* p_primal = 0;
   try{
      solver->getDualFarkas(*p_primal);
   }catch(SPxStatusException& x)
   {
      delete solver;
      if(x.what().find("XSOLVE10") == 0)
         return true;
      return false;
   }
   delete solver;
   return false;
}

bool StatusExceptionCheck::checkXSOLVE11()
{
   SPxSolver* solver = new SPxSolver();
   Vector* p_primal = 0;
   try{
      solver->getSlacks(*p_primal);
   }catch(SPxStatusException& x)
   {
      delete solver;
      if(x.what().find("XSOLVE11") == 0)
         return true;
      return false;
   }
   delete solver;
   return false;
}

bool StatusExceptionCheck::noSolverCheck()
{
   std::cout<<std::endl<<"Test for No Solver:"<<std::endl;
   std::string testString="XSOLVE02";
   std::cout<<" "<<testString<<": ";
   SPxSolver* solver = new SPxSolver();
   NameSet rownames, colnames;
   const char* filename  = "../check/LP/netlib/adlittle.mps.gz";
   solver->readFile(filename, &rownames, &colnames);
   try{
      solver->solve();
   }catch(SPxStatusException& x)
   {
      delete solver;
      if(x.what().find(testString) == 0)
      {
         std::cout<<"ok"<<std::endl;
         return true;
      }
      std::cout<<"failed"<<std::endl;
      return false;
   }
   std::cout<<"failed"<<std::endl;
   delete solver;
   return false;
}

bool StatusExceptionCheck::noPricerCheck()
{
   std::cout<<std::endl<<"Test for No Pricer:"<<std::endl;
   std::string testString="XSOLVE03";
   std::cout<<" "<<testString<<": ";
   SPxSolver* solver = new SPxSolver();
   NameSet rownames, colnames;
   const char* filename  = "../check/LP/netlib/adlittle.mps.gz";
   SLUFactor<Real> slu;
   solver->readFile(filename, &rownames, &colnames);
   solver->setBasisSolver(&slu);
   try{
      solver->solve();
   }catch(SPxStatusException& x)
   {
      delete solver;
      if(x.what().find(testString) == 0)
      {
         std::cout<<"ok"<<std::endl;
         return true;
      }
      std::cout<<"failed"<<std::endl;
      return false;
   }
   std::cout<<"failed"<<std::endl;
   delete solver;
   return false;
}

bool StatusExceptionCheck::noRatioTesterCheck()
{
   std::cout<<std::endl<<"Test for No RatioTester:"<<std::endl;
   std::string testString="XSOLVE04";
   std::cout<<" "<<testString<<": ";
   SPxSolver* solver = new SPxSolver();
   NameSet rownames, colnames;
   const char* filename  = "../check/LP/netlib/adlittle.mps.gz";
   SLUFactor<Real> slu;
   SPxSteepPR pricer;
   solver->readFile(filename, &rownames, &colnames);
   solver->setBasisSolver(&slu);
   solver->setPricer(&pricer);
   try{
      solver->solve();
   }catch(SPxStatusException& x)
   {
      delete solver;
      if(x.what().find(testString) == 0)
      {
         std::cout<<"ok"<<std::endl;
         return true;
      }
      std::cout<<"failed"<<std::endl;
      return false;
   }
   std::cout<<"failed"<<std::endl;
   delete solver;
   return false;
}

bool StatusExceptionCheck::notInitialized()
{
   std::cout<<std::endl<<"Test for NotInitialized:"<<std::endl;
   std::string testString="XSOLVE06";
   std::cout<<" "<<testString<<": ";
   SPxSolver* solver = new SPxSolver();
   NameSet rownames, colnames;
   const char* filename  = "../check/LP/netlib/adlittle.mps.gz";
   SLUFactor<Real> slu;
   SPxSteepPR pricer;
   SPxDefaultRT tester;
   Vector* v;
   solver->readFile(filename, &rownames, &colnames);
   solver->setBasisSolver(&slu);
   solver->setPricer(&pricer);
   solver->setTester(&tester);
   try{
      solver->getPrimalSol(*v);
   }catch(SPxStatusException& x)
   {
      delete solver;
      if(x.what().find(testString) == 0)
      {
         std::cout<<"ok"<<std::endl;
         return true;
      }
      std::cout<<"failed"<<std::endl;
      return false;
   }
   std::cout<<"failed"<<std::endl;
   delete solver;
   return false;
}

int main(int argc, char* argv[])
{
   StatusExceptionCheck checks;

   if(checks.noProblemCheck())
      std::cout<<"ok"<<std::endl;
   else
      std::cout<<"failed"<<std::endl;

   if(checks.noSolverCheck())
      std::cout<<"ok"<<std::endl;
   else
      std::cout<<"failed"<<std::endl;

   if(checks.noPricerCheck())
      std::cout<<"ok"<<std::endl;
   else
      std::cout<<"failed"<<std::endl;

   if(checks.noRatioTesterCheck())
      std::cout<<"ok"<<std::endl;
   else
      std::cout<<"failed"<<std::endl;

   if(checks.notInitialized())
      std::cout<<"ok"<<std::endl;
   else
      std::cout<<"failed"<<std::endl;

}
