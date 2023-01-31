/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the class library                   */
/*       SoPlex --- the Sequential object-oriented simPlex.                  */
/*                                                                           */
/*  Copyright (c) 1996-2023 Zuse Institute Berlin (ZIB)                      */
/*                                                                           */
/*  Licensed under the Apache License, Version 2.0 (the "License");          */
/*  you may not use this file except in compliance with the License.         */
/*  You may obtain a copy of the License at                                  */
/*                                                                           */
/*      http://www.apache.org/licenses/LICENSE-2.0                           */
/*                                                                           */
/*  Unless required by applicable law or agreed to in writing, software      */
/*  distributed under the License is distributed on an "AS IS" BASIS,        */
/*  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied. */
/*  See the License for the specific language governing permissions and      */
/*  limitations under the License.                                           */
/*                                                                           */
/*  You should have received a copy of the Apache-2.0 license                */
/*  along with SoPlex; see the file LICENSE. If not email to soplex@zib.de.  */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */


#include <sstream>

#include <math.h>
#include <string.h>
#include <iostream>
#include <iomanip>
#include <fstream>

#include "spxdefines.h"
#include "soplex.h"
#include "spxsolver.h"
#include <stdio.h>
#include <sys/types.h>
#include <sys/stat.h>

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

///see "soplex -help", s means Simplifier, g means Scaler, c means Starter, p for Pricer and t for Ratiotester
void set_s_g_c_p_t    (SoPlex& work,
		       const int simplifying,            // 0 for none and 1 for Main, please see "soplex -help"
		       const int scaling,                // see above
		       const int starting,
		       const int pricing,
		       const int ratiotest);
void set_simplifier   (SoPlex& work, const int simplifying);  // 0 for none and 1 for Main, please see "soplex -help"
void set_scalers      (SoPlex& work, const int scaling    );  // see above
void set_starter      (SoPlex& work, const int starting   );  // see above
void set_pricer       (SoPlex& work, const int pricing    );  // see above
void set_ratio_tester (SoPlex& work, const int ratiotest  );  // see above
double abs(double a, double b);

/** Unit test for copy constructors and assignment operators with random combinations of pricer,
 *  ratiotester, simplifier, scaler, and starter.
 *
 *  NOTE: Fails from singular bases usually come from using "bad" pricers or ratiotesters etc.,
 *  not from incorrect copying.   */

int main(int argc, const char* const argv[])
{
  if(argc < 2)
  {
    std::cerr << "You must input a test file name, e.g. quick.test." << std::endl;
    exit(0);
  }

  std::ifstream in;
  in.open(argv[argc-1]);

  if(!in)
  {
     std::cerr << "No such file "<< argv[argc-1] << " can be found," << std::endl;
    exit(0);
  }

  Param::setVerbose ( 0 );

  int times = 0;
  int succ_times = 0;
  const int num_rows_output = 6;
  const int num_type = 2;
  const int num_timepoint = 3;
  const int num_copies = 3;
  int cnt[num_rows_output];
  int pass[num_rows_output];
  int fail[num_rows_output];
  int singularBasis[num_rows_output];
  int bigErr[num_rows_output];
  std::string type[] = {"assign", "copy"};
  std::string timepoint[] = {"before_load", "after_load", "after_solve"};

  for(int num = 0; num < num_rows_output; ++num)
  {
     cnt[num]           = 0;
     pass[num]          = 0;
     fail[num]          = 0;
     singularBasis[num] = 0;
     bigErr[num]        = 0;
  }

  while(!in.eof())
  {
     std::string ffname;
     std::getline(in, ffname);
     if(ffname.length() > 1)
     {
	if(ffname[0] == '/' && ffname[0] == '/')
        {
           std::cout << ffname << "skipped\n";
           continue;
	}
	else
           times++;
	std::cout << "\n\n\nReading file " << ffname << "..." <<std::endl;
	std::string fname = "../check/LP/" + ffname;



        SoPlex lp_ori(SPxSolver::LEAVE, SPxSolver::COLUMN);
        set_s_g_c_p_t(lp_ori, 1, 1, 0, 4, 2);


        SoPlex lp[num_timepoint][num_type][num_copies] ;
        lp[0][1][0] = SoPlex(lp_ori);
        lp[0][1][1] = SoPlex(lp_ori);
        lp[0][1][2] = SoPlex(lp_ori);
        set_s_g_c_p_t(lp[0][1][0], 1, 1, 0, 0, 0);

        lp[0][0][0] = SoPlex(SPxSolver::LEAVE, SPxSolver::ROW);
        lp[0][0][1] = SoPlex(SPxSolver::ENTER, SPxSolver::COLUMN);
        lp[0][0][2] = SoPlex(SPxSolver::ENTER, SPxSolver::ROW);
        lp[0][0][0] = lp_ori;
        lp[0][0][1] = lp_ori;
        lp[0][0][2] = lp_ori;
        set_s_g_c_p_t(lp[0][0][0], 0, 4, 3, 3, 0);
        set_s_g_c_p_t(lp[0][0][2], 0, 0, 2, 4/*5*/, 2);

        lp[0][1][1].setType(SPxSolver::ENTER);
        lp[0][0][0].setRep(SPxSolver::ROW);

        if(!lp_ori.readFile(fname.c_str(), 0, 0, 0))
        {
           std::cout << "error while reading file '"  << fname.c_str() << "' -- terminating." << std::endl;
           exit(1);
        }

        // get the optimal objValue of the lp
        double optObjValue = 0.0;
        SoPlex* lp_cp_pr = new SoPlex(lp_ori);
        try{

           lp_cp_pr->optimize();
           if(lp_cp_pr->status() == SPxSolver::OPTIMAL)
           {
              std::cout << "Lp has been solved to optimality with default settings, objValue is " << lp_cp_pr->objValue() << std::endl;
           }
           else
           {
              std::cout << "Lp '"  << fname.c_str()
                        << "' has not been solved optimally with default settings -- terminating. (Choose a simpler test set!)"
                        << std::endl;
              exit(1);
           }
           optObjValue = lp_cp_pr->objValue();

           delete lp_cp_pr;

        }catch(SPxException& x) {
           std::cout << "exception caught : " << x.what() << std::endl;
        }

        lp[0][1][0].readFile(fname.c_str(), 0, 0, 0);
        lp[0][1][1].readFile(fname.c_str(), 0, 0, 0);
        lp[0][1][2].readFile(fname.c_str(), 0, 0, 0);
        set_s_g_c_p_t(lp[0][1][1], 0, 2, 2, 1, 1);
        set_s_g_c_p_t(lp[0][1][2], 1, 3, 1, 2, 2);

        lp[0][0][0].readFile(fname.c_str(), 0, 0, 0);
        lp[0][0][1].readFile(fname.c_str(), 0, 0, 0);
        lp[0][0][2].readFile(fname.c_str(), 0, 0, 0);
        set_s_g_c_p_t(lp[0][0][2], 0, 0, 2, 4/*5*/, 2);
        set_s_g_c_p_t(lp[0][0][1], 0, 0, 2, 4/*5*/, 2);
        //set_s_g_c_p_t(lp[0][0][1], 1, 2, 0, 4, 1);

        lp[0][1][2].setType(SPxSolver::ENTER);
        lp[0][0][1].setRep(SPxSolver::ROW);
        lp[0][0][1].setType(SPxSolver::LEAVE);

        /** test copy constructor and assignment function before soplex reads a lp    */


        int current_timepoint = 0;
        for(int n_type = 0; n_type < num_type; ++n_type)
        {
           int stat_idx = 2*current_timepoint + n_type;
           for(int n_time = 0; n_time < num_copies; ++n_time)
           {
              try{
                 (cnt[stat_idx])++;
                 lp[current_timepoint][n_type][n_time].optimize();
                 if(lp[current_timepoint][n_type][n_time].status() == SPxSolver::OPTIMAL)
                 {
                    double absError = abs(lp[current_timepoint][n_type][n_time].objValue(), optObjValue);
                    if(absError < lp[current_timepoint][n_type][n_time].delta() )
                    {
                       (pass[stat_idx])++;
                       std::cout << "lp_" << type[n_type] << "_" << timepoint[current_timepoint]
                                 << "_" << n_time <<  " solved, num of iterations is "
                                 << lp[current_timepoint][n_type][n_time].iteration() << std::endl;
                    }
                    else
                    {
                       if( absError / optObjValue <  lp[current_timepoint][n_type][n_time].delta() )
                       {
                          (pass[stat_idx])++;
                       }
                       else
                       {
                          (fail[stat_idx])++;
                          (bigErr[stat_idx])++;
                          std::cout << "lp_" << type[n_type] << "_" << timepoint[current_timepoint] << "_" << n_time <<  " solved "
                                    << "but does not have the same objValue with original, num of iterations is "
                                    << lp[current_timepoint][n_type][n_time].iteration() <<  ", the difference is "
                                    << abs(lp[current_timepoint][n_type][n_time].objValue(), optObjValue)  << std::endl;
                       }
                    }
                 }
              }catch(SPxException& x) {
                 (fail[stat_idx])++;
                 std::cout << "exception caught by" << " lp_" << type[n_type] << "_" << timepoint[current_timepoint] << "_" << n_time << " :"  << x.what() << std::endl;
                 std::string exceptionContent(x.what());
                 if(exceptionContent.find("XSOLVE21") != std::string::npos || exceptionContent.find("singular") != std::string::npos)
                 {
                    (singularBasis[stat_idx])++;
                 }
                 std::cout << "\n\n\n";
              }

           }
        }

        /** test copy constructor and assignment function after soplex reads a lp    */

        lp[1][1][0] = SoPlex(lp_ori);
        lp[1][1][1] = SoPlex(lp_ori);
        lp[1][1][2] = SoPlex(lp_ori);
        set_s_g_c_p_t(lp[1][1][0], 0, 2, 2, 1, 1);
        set_s_g_c_p_t(lp[1][1][1], 0, 4, 3, 3, 0);
        set_s_g_c_p_t(lp[1][1][2], 0, 0, 2, 4/*5*/, 2);

        lp[1][0][0] = SoPlex(SPxSolver::LEAVE, SPxSolver::ROW);
        lp[1][0][1] = SoPlex(SPxSolver::ENTER, SPxSolver::COLUMN);
        lp[1][0][2] = SoPlex(SPxSolver::ENTER, SPxSolver::ROW);
        set_s_g_c_p_t(lp[1][0][0], 1, 1, 0, 0, 0);
        lp[1][0][0] = lp_ori;
        lp[1][0][1] = lp_ori;
        lp[1][0][2] = lp_ori;
        set_s_g_c_p_t(lp[1][0][1], 1, 3, 1, 2, 2);
        set_s_g_c_p_t(lp[1][0][2], 1, 2, 0, 4, 1);


        lp[1][0][2].setType(SPxSolver::ENTER); // see above
        lp[1][1][2].setType(SPxSolver::ENTER);
        lp[1][0][0].setRep(SPxSolver::ROW);
        lp[1][0][1].setRep(SPxSolver::ROW);


        current_timepoint = 1;
        for(int n_type = 0; n_type < num_type; ++n_type)
        {
           int stat_idx = 2*current_timepoint + n_type;
           for(int n_time = 0; n_time < num_copies; ++n_time)
           {
              try{
                 (cnt[stat_idx])++;
                 lp[current_timepoint][n_type][n_time].optimize();
                 if(lp[current_timepoint][n_type][n_time].status() == SPxSolver::OPTIMAL)
                 {
                    double absError = abs(lp[current_timepoint][n_type][n_time].objValue(), optObjValue);
                    if(absError < lp[current_timepoint][n_type][n_time].delta() )
                    {
                       (pass[stat_idx])++;
                       std::cout << "lp_" << type[n_type] << "_"
                                 << timepoint[current_timepoint] << "_" << n_time
                                 <<  " solved, num of iterations is "
                                 << lp[current_timepoint][n_type][n_time].iteration() << std::endl;
                    }
                    else
                    {
                       if( absError / optObjValue <  lp[current_timepoint][n_type][n_time].delta() )
                       {
                          (pass[stat_idx])++;
                       }
                       else
                       {
                          (fail[stat_idx])++;
                          (bigErr[stat_idx])++;
                          std::cout << "lp_" << type[n_type] << "_" << timepoint[current_timepoint] << "_" << n_time <<  " solved "
                                    << "but does not have the same objValue with original, num of iterations is "
                                    << lp[current_timepoint][n_type][n_time].iteration() <<  ", the difference is "
                                    << abs(lp[current_timepoint][n_type][n_time].objValue(), optObjValue)  << std::endl;
                       }
                    }
                 }
              }catch(SPxException& x) {
                 (fail[stat_idx])++;
                 std::cout << "exception caught by" << " lp_" << type[n_type] << "_" << timepoint[current_timepoint] << "_" << n_time << " :"  << x.what() << std::endl;
                 std::string exceptionContent(x.what());
                 if(exceptionContent.find("XSOLVE21") != std::string::npos || exceptionContent.find("singular") != std::string::npos)
                 {
                    (singularBasis[stat_idx])++;
                 }
                 std::cout << "\n\n\n";
              }

           }
        }

        // solve the original problem

        try{

           lp_ori.optimize();
           if(lp_ori.status() == SPxSolver::OPTIMAL)
           {
              std::cout << "The original lp has been solved to optimality with default settings, num of iterations is "
                        << lp_ori.iteration() << std::endl;
           }
           else
           {
              std::cout << "The original lp '"  << fname.c_str()
                        << "' has not been solved optimally with default settings -- terminating. (Choose a simpler test set!)"
                        << std::endl;
              exit(1);
           }

        }catch(SPxException& x) {
           std::cout << "exception caught : " << x.what() << std::endl;
           std::cout << "\n\n\n";
        }


        /** test copy constructor and assignment function before soplex solved an lp    */


        lp[2][1][0] = SoPlex(lp_ori);
        lp[2][1][1] = SoPlex(lp_ori);
        lp[2][1][2] = SoPlex(lp_ori);
        set_s_g_c_p_t(lp[2][1][0], 1, 3, 1, 2, 2);
        set_s_g_c_p_t(lp[2][1][1], 0, 4, 3, 3, 0);
        set_s_g_c_p_t(lp[2][1][2], 1, 2, 0, 4, 1);


        lp[2][0][0] = SoPlex(SPxSolver::LEAVE, SPxSolver::ROW);
        lp[2][0][1] = SoPlex(SPxSolver::ENTER, SPxSolver::COLUMN);
        lp[2][0][2] = SoPlex(SPxSolver::ENTER, SPxSolver::ROW);
        set_s_g_c_p_t(lp[2][0][0], 1, 1, 0, 0, 0);

        lp[2][0][0] = lp_ori;
        lp[2][0][1] = lp_ori;
        lp[2][0][2] = lp_ori;
        set_s_g_c_p_t(lp[2][0][1], 0, 2, 2, 1, 1);
        set_s_g_c_p_t(lp[2][0][2], 0, 0, 2, 4/*5*/, 2);

        lp[2][1][0].setType(SPxSolver::ENTER);
        lp[2][0][2].setRep(SPxSolver::ROW);
        lp[2][0][1].setRep(SPxSolver::ROW);


        current_timepoint = 2;
        for(int n_type = 0; n_type < num_type; ++n_type)
        {
           int stat_idx = 2*current_timepoint + n_type;
           for(int n_time = 0; n_time < num_copies; ++n_time)
           {
              try{
                 (cnt[stat_idx])++;
                 lp[current_timepoint][n_type][n_time].optimize();
                 if(lp[current_timepoint][n_type][n_time].status() == SPxSolver::OPTIMAL)
                 {
                    double absError = abs(lp[current_timepoint][n_type][n_time].objValue(), optObjValue);
                    if(absError < lp[current_timepoint][n_type][n_time].delta() )
                    {
                       (pass[stat_idx])++;
                       std::cout << "lp_" << type[n_type] << "_" << timepoint[current_timepoint]
                                 << "_" << n_time <<  " solved, num of iterations is "
                                 << lp[current_timepoint][n_type][n_time].iteration() << std::endl;
                    }
                    else
                    {
                       if( absError / optObjValue <  lp[current_timepoint][n_type][n_time].delta() )
                       {
                          (pass[stat_idx])++;
                       }
                       else
                       {
                          (fail[stat_idx])++;
                          (bigErr[stat_idx])++;
                          std::cout << "lp_" << type[n_type] << "_" << timepoint[current_timepoint] << "_" << n_time <<  " solved "
                                    << "but does not have the same objValue with original, num of iterations is "
                                    << lp[current_timepoint][n_type][n_time].iteration() <<  ", the difference is "
                                    << abs(lp[current_timepoint][n_type][n_time].objValue(), optObjValue)  << std::endl;
                       }
                    }
                 }
              }catch(SPxException& x) {
                 (fail[stat_idx])++;
                 std::cout << "exception caught by" << " lp_"
                           << type[n_type] << "_" << timepoint[current_timepoint] << "_"
                           << n_time << " :"  << x.what() << std::endl;
                 std::string exceptionContent(x.what());
                 if(exceptionContent.find("XSOLVE21") != std::string::npos || exceptionContent.find("singular") != std::string::npos)
                 {
                    (singularBasis[stat_idx])++;
                 }
                 std::cout << "\n\n\n";
              }

           }
        }

     }
  }
  // test end

  // output the result
  std::cout << "\n\n";
  std::cout << "---------------------------------------------------------------------------------------\n";
  std::cout << "Type                                Cnt      Pass     Fail      SingularBasis    BigErr\n";
  std::cout << "---------------------------------------------------------------------------------------\n";

  for(int n_type = 0; n_type < num_type; ++n_type)
  {
     for(int n_timepoint = 0; n_timepoint < num_timepoint; ++n_timepoint)
     {
        int idx = 2*n_timepoint + n_type;
        std::cout << std::setw(14) << std::left <<  type[n_type]
                  << std::setw(14) << std::right <<  timepoint[n_timepoint]
                  << std::setw(11) << std::right <<  cnt[idx]
                  << std::setw(10) << std::right <<  pass[idx]
                  << std::setw(9) << std::right <<  fail[idx]
                  << std::setw(19) << std::right <<  singularBasis[idx]
                  << std::setw(10) << std::right <<  bigErr[idx]
                  << std::endl;
     }
  }

  int sum_cnt = 0;
  int sum_pass = 0;
  int sum_fail = 0;
  int sum_singularBasis = 0;
  int sum_bigErr = 0;

  for(int num = 0; num < num_rows_output; ++num)
  {
     sum_cnt += cnt[num];
     sum_pass += pass[num];
     sum_fail += fail[num];
     sum_singularBasis += singularBasis[num];
     sum_bigErr += bigErr[num];
  }

  std::cout << "---------------------------------------------------------------------------------------\n";
  std::cout << std::setw(14) << std::left <<  "Sum"
            << std::setw(14) << std::right <<  ""
            << std::setw(11) << std::right <<  sum_cnt
            << std::setw(10) << std::right <<  sum_pass
            << std::setw(9) << std::right <<  sum_fail
            << std::setw(19) << std::right <<  sum_singularBasis
            << std::setw(10) << std::right <<  sum_bigErr
            << std::endl;
  std::cout << "---------------------------------------------------------------------------------------\n";

}

void set_simplifier(SoPlex& work, const int simplifying)
{
   switch(simplifying)
   {
   case 1 :
     work.setSimplifier(new SPxMainSM(), true);
     break;
   case 0  :
     /*FALLTHROUGH*/
   default :
     break;
   }

}

void set_scalers(SoPlex& work, const int scaling)
{
   switch(scaling)
   {
   case 4:
     work.setPreScaler(new SPxEquiliSC(true), true);
     work.setPostScaler(new SPxGeometSC(8), true);
      break;
   case 3:
     work.setPreScaler(new SPxEquiliSC(true), true);
     work.setPostScaler(new SPxGeometSC(1), true);
      break;
   case 2 :
     work.setPreScaler(new SPxEquiliSC(true), true);
     work.setPostScaler(0, false);
      break;
   case 1 :
     work.setPreScaler(new SPxEquiliSC(false), true);
     work.setPostScaler(0, false);
      break;
   case 0 :
      /*FALLTHROUGH*/
   default :
     work.setPreScaler(0, false);
     work.setPostScaler(0, false);
      break;
   }
}

void set_starter(SoPlex& work, const int starting)
{
   switch(starting)
   {
   case 3 :
     work.setStarter(new SPxVectorST, true);
     break;
   case 2 :
     work.setStarter(new SPxSumST, true);
     break;
   case 1 :
     work.setStarter(new SPxWeightST, true);
     break;
   case 0 :
      /*FALLTHROUGH*/
   default :
     work.setStarter(0, false);
     break;
   }
}

void set_pricer(SoPlex& work, const int pricing)
{
   switch(pricing)
   {
   case 5 :
     work.setPricer(new SPxWeightPR, true);
     break;
   case 4 :
     work.setPricer(new SPxSteepPR, true);
     break;
   case 3 :
     work.setPricer(new SPxHybridPR, true);
      break;
   case 2 :
     work.setPricer(new SPxDevexPR, true);
     break;
   case 1 :
     work.setPricer(new SPxParMultPR, true);
     break;
   case 0 :
      /*FALLTHROUGH*/
   default :
     work.setPricer(new SPxDantzigPR, true);
     break;
   }
}

void set_ratio_tester(SoPlex& work, const int ratiotest)
{
   switch(ratiotest)
   {
   case 2 :
     work.setTester(new SPxFastRT, true);
     break;
   case 1 :
     work.setTester(new SPxHarrisRT, true);
     break;
   case 0 :
     /*FALLTHROUGH*/
   default:
     work.setTester(new SPxDefaultRT, true);
     break;
   }
}

void set_s_g_c_p_t    (SoPlex& work,
		       const int simplifying,
		       const int scaling,
		       const int starting,
		       const int pricing,
		       const int ratiotest)
{
  set_simplifier   (work,  simplifying);
  set_scalers      (work,  scaling    );
  set_starter      (work,  starting   );
  set_pricer       (work,  pricing    );
  set_ratio_tester (work,  ratiotest  );
}

double abs(double a, double b)
{
  if(a >= b)
    return a-b;
  else
    return b-a;
}

//-----------------------------------------------------------------------------
//Emacs Local Variables:
//Emacs mode:c++
//Emacs c-basic-offset:3
//Emacs tab-width:8
//Emacs indent-tabs-mode:nil
//Emacs End:
//-----------------------------------------------------------------------------
