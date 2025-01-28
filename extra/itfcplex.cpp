/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the class library                   */
/*       SoPlex --- the Sequential object-oriented simPlex.                  */
/*                                                                           */
/*  Copyright (c) 1996-2025 Zuse Institute Berlin (ZIB)                      */
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

#include <iostream>
#include <fstream>
#include <string.h>
#include <assert.h>

#include <cplex.h>

#include "spxalloc.h"
#include "soplex.h"
#include "slufactor.h"
#include "spxsteeppr.h"
#include "spxfastrt.h"
#include "spxweightst.h"
#include "nameset.h"
#include "didxset.h"

using namespace soplex;

class SPxCPlex : public SoPlex
{
   SLUFactor<Real>    m_slu;
   SPxSteepPR   m_price;
   SPxFastRT    m_ratio;
   bool         m_verbose;
   char*        m_probname;
   NameSet      m_colnames;
   NameSet      m_rownames;
   DIdxSet      m_intvars;

public:
    void factorize( )
    {
       SoPlex::factorize();

       if (m_verbose)
       {
          std::cout.precision(16) ;
          std::cout << (type() == LEAVE ? "L " : "E ") ;
          std::cout << basis().iteration() << ":\t" ;
          std::cout << value() ;
          std::cout << "\t(" << shift() << ")" << std::endl;;
       }
    }

   /* Was macht das? Brauchen wir das?
    * void splitLP()
    * {
    *   subcovectors.reSize( 1 ) ;
    * }
    */

   SPxCPlex()
      : SoPlex(LEAVE, COLUMN)
      , m_verbose(true)
      , m_probname(0)
   {
      loadSolver(&m_slu );
      loadTester(&m_ratio);
      loadPricer(&m_price);
      loadStarter(0);
   }
   ~SPxCPlex()
   {
      if (m_probname != 0)
         spx_free(m_probname);
   }
   // This is public in SPxBasis, but protected inherted from SoPlex.
   SPxBasis::Desc::Status dualColStatus(int i) const
   {
      return SPxBasis::dualColStatus(i);
   }
   // This is public in SPxBasis, but protected inherted from SoPlex.
   SPxBasis::Desc::Status dualRowStatus(int i) const
   {
      return SPxBasis::dualRowStatus(i);
   }
   void setProbname(const char* p_probname)
   {
      assert(p_probname != 0);
      if (m_probname != 0)
         spx_free(m_probname);
      spx_alloc(m_probname, strlen(p_probname + 1));
      strcpy(m_probname, p_probname);
   }
   void setVerbose(bool p_verbose)
   {
      m_verbose = p_verbose;
   }
};

extern "C" CPXENVptr CPXopenCPLEX(int* status_p)
{
   SPxCPlex* spx = new SPxCPlex;
   *status_p     = 0;

   return reinterpret_cast<CPXENVptr>(new SPxCPlex);
}

extern "C" int CPXcloseCPLEX(CPXENVptr *env)
{
   if (env == 0)
      return CPXERR_NULL_POINTER;
   if (*env == 0)
      return CPXERR_NO_ENVIRONMENT;

   delete reinterpret_cast<SPxCPlex*>(*env);
   *env = 0;
   return 0;
}

extern "C" char* CPXversion(CPXENVptr env)
{
   assert(env != 0);

   static char* version = "SoPlex 1.3";

   return version;
}

extern "C" CPXLPptr CPXcreateprob(CPXENVptr env, int* status_p, char* probname)
{
   if (env == 0)
   {
      *status_p = CPXERR_NO_ENVIRONMENT;
      return 0;
   }
   SPxCPlex* spx = reinterpret_cast<SPxCPlex*>(env);

   if (probname != 0)
      spx->setProbname(probname);

   *status_p = 0;

   return reinterpret_cast<CPXLPptr>(spx);
}

extern "C" int CPXfreeprob(CPXENVptr env, CPXLPptr* lp)
{
   if (env == 0)
      return CPXERR_NO_ENVIRONMENT;
   if (lp == 0)
      return CPXERR_NULL_POINTER;

   *lp = 0;

   return 0;
}

/**@todo The names are ignored. And ranges seem not to be supported.
 */
extern "C" int CPXcopylpwnames(
   CPXENVptr env,
   CPXLPptr  lp,
   int       numcols, // mac
   int       numrows, // mar
   int       objsen,
   double*   objx,
   double*   rhsx,
   char*     senx,
   int*      matbeg,
   int*      matcnt,
   int*      matind,
   double*   matval,
   double*   bdl,
   double*   bdu,
   double*   rngval,
   char**    colname,
   char**    rowname)
{
   if (env == 0)
      return CPXERR_NO_ENVIRONMENT;
   if (lp == 0)
      return CPXERR_NULL_POINTER;
   if ((numcols < 1) || (numrows < 1) || (objsen == 0))
      return CPXERR_BAD_ARGUMENT;
   if ((objx == 0) || (rhsx == 0) || (senx == 0)
      || (matbeg == 0) || (matcnt == 0) || (matind == 0) || (matval == 0)
      || (bdl == 0) || (bdu == 0) || (colname == 0) || (rowname == 0))
      return CPXERR_NULL_POINTER;

   SPxCPlex* spx = reinterpret_cast<SPxCPlex*>(env);
   LPColSet  cols(numcols);
   LPRowSet  rows(numrows);
   DSVector  colVector(numrows);
   DSVector  emptyVector(0);
   LPRow     objRow(-SPxLP::infinity, emptyVector, SPxLP::infinity) ;

   spx->addRow(objRow);

   for(int i = 0; i < numrows; i++)
   {
      switch(senx[i])
      {
      case 'L':
         rows.add(-SPxLP::infinity, emptyVector, rhsx[i]);
         break;
      case 'G':
         rows.add(rhsx[i], emptyVector, SPxLP::infinity);
         break;
      case 'E':
         rows.add(rhsx[i], emptyVector, rhsx[i]);
         break;
      case 'R':
         rows.add(rhsx[i], emptyVector, rhsx[i]);
         if (rngval == 0)
            return CPXERR_NULL_POINTER;
         break ;
      }
   }
   spx->addRows(rows);

   for(int i = 0; i < numcols; i++)
   {
      colVector.clear() ;

      for(int j = 0; j < matcnt[i]; j++)
         colVector.add(matind[matbeg[i] + j] + 1, matval[matbeg[i] + j]);

      colVector.add(0, (objsen == 1) ? -objx[i] : objx[i]);

      cols.add(objx[i], bdl[i], colVector, bdu[i]);
   }
   spx->changeSense(objsen == 1 ? SPxLP::MINIMIZE : SPxLP::MAXIMIZE);
   spx->addCols(cols);

   return 0;
}

extern "C" int CPXcopybase(
   CPXENVptr env, CPXLPptr lp, int* cstat, int* rstat)
{
   if (env == 0)
      return CPXERR_NO_ENVIRONMENT;
   if ((lp == 0) || (cstat == 0) || (rstat == 0))
      return CPXERR_NULL_POINTER;

   SPxCPlex* spx = reinterpret_cast<SPxCPlex*>(env);

   SPxBasis::Desc  desc;

   desc.reSize(spx->nRows(), spx->nCols()) ;

   if (spx->basis().solver() != spx)
      spx->basis().load(spx);
   //((SPxBasis*)&spx->basis())->load( solver ) ;

   //for(int i = spx->nRows() -1 ; i-- ; )      // obj limit in row 0!!!
   for(int i = 1; i < spx->nRows(); i++)
   {
      switch(rstat[i])
      {
      case CPX_AT_LOWER :
      case CPX_AT_UPPER :
         if (spx->rhs(i + 1) == spx->lhs(i + 1))
            desc.rowStatus(i + 1) = SPxBasis::Desc::P_FIXED;
         else if (spx->rhs(i + 1) >= SPxLP::infinity)
         {
            assert(spx->lhs(i + 1) > -SPxLP::infinity);
            desc.rowStatus(i + 1) = SPxBasis::Desc::P_ON_LOWER;
         }
         else if (spx->lhs(i + 1) <= -SPxLP::infinity)
         {
            assert(spx->rhs(i + 1) < SPxLP::infinity) ;
            desc.rowStatus(i + 1) = SPxBasis::Desc::P_ON_UPPER;
         }
         else
            abort();
         break ;
      case CPX_BASIC :
         desc.rowStatus(i + 1) = spx->dualRowStatus(i + 1);
         break;
      default:
         return CPXERR_BAD_ARGUMENT;
      }
   }
   desc.rowStatus(0) = spx->dualRowStatus(0);

   for(int i = 0; i < spx->nCols(); i++)
   {
      switch(cstat[i])
      {
      case 0:
         if (spx->upper(i) == spx->lower(i))
            desc.colStatus(i) = SPxBasis::Desc::P_FIXED;
         else
            desc.colStatus(i) = SPxBasis::Desc::P_ON_LOWER;
         break;
      case 2:
         if (spx->upper(i) == spx->lower(i))
            desc.colStatus(i) = SPxBasis::Desc::P_FIXED;
         else
            desc.colStatus(i) = SPxBasis::Desc::P_ON_UPPER;
         break;
      case 1:
         desc.colStatus(i) = spx->dualColStatus(i);
         break;
      case 3:
         desc.colStatus(i) = SPxBasis::Desc::P_FREE;
         break;
      default:
         return CPXERR_BAD_ARGUMENT;
      }
   }
   spx->load(desc);

   return 0;
}

extern "C" int CPXgetbase(
   CPXENVptr env, CPXLPptr lp, int *cstat, int *rstat)
{
   if (env == 0)
      return CPXERR_NO_ENVIRONMENT;
   if (lp == 0)
      return CPXERR_NULL_POINTER;

   SPxCPlex*             spx  = reinterpret_cast<SPxCPlex*>(env);
   const SPxBasis::Desc& desc = spx->basis().desc();

   if (cstat != 0)
   {
      for(int i = 0; i < spx->nCols(); i++)
      {
         switch(desc.colStatus(i))
         {
         case SPxBasis::Desc::P_ON_LOWER:
            cstat[i] = CPX_AT_LOWER;
            break ;
         case SPxBasis::Desc::P_ON_UPPER:
         case SPxBasis::Desc::P_FIXED:
            cstat[i] = CPX_AT_UPPER;
            break ;
         case SPxBasis::Desc::P_FREE:
            cstat[i] = CPX_FREE_SUPER;
            break ;
         case SPxBasis::Desc::D_ON_UPPER:
         case SPxBasis::Desc::D_ON_LOWER:
         case SPxBasis::Desc::D_ON_BOTH:
         case SPxBasis::Desc::D_UNDEFINED:
         case SPxBasis::Desc::D_FREE:
            cstat[i] = CPX_BASIC;
            break ;
         default:
            abort();
         }
      }
   }
   if (rstat != 0)
   {
      for(int i = 0; i < spx->nRows() - 1; i++)
      {
         switch(desc.rowStatus(i + 1))
         {
         case SPxBasis::Desc::P_ON_LOWER:
         case SPxBasis::Desc::P_ON_UPPER:
         case SPxBasis::Desc::P_FIXED:
         case SPxBasis::Desc::P_FREE:
            rstat[i] = CPX_AT_LOWER;
            break ;
         case SPxBasis::Desc::D_ON_UPPER:
         case SPxBasis::Desc::D_ON_LOWER:
         case SPxBasis::Desc::D_ON_BOTH:
         case SPxBasis::Desc::D_UNDEFINED:
         case SPxBasis::Desc::D_FREE:
            rstat[i] = CPX_BASIC;
            break ;
         default:
            abort();
         }
      }
   }
   return 0;
}

extern "C" int CPXdualopt(CPXENVptr env, CPXLPptr lp)
{
   if (env == 0)
      return CPXERR_NO_ENVIRONMENT;
   if (lp == 0)
      return CPXERR_NULL_POINTER;

   SPxCPlex* spx  = reinterpret_cast<SPxCPlex*>(env);

   if (spx->basis().status() <= 0)
      spx->setType(SoPlex::ENTER);

   spx->optimize();

   return 0;
}

/**@todo The Status codes returned are partly bullshit.
 * @todo The return value of CPXsolution is not allways ok.
 */
extern "C" int CPXsolution(
   CPXENVptr env,
   CPXLPptr  lp,
   int*      status,
   double*   obj,
   double*   primal,
   double*   dual,
   double*   slack,
   double*   redcost)
{
   if (env == 0)
      return CPXERR_NO_ENVIRONMENT;
   if (lp == 0)
      return CPXERR_NULL_POINTER;

   SPxCPlex* spx = reinterpret_cast<SPxCPlex*>(env);

   if (status != 0)
   {
      switch(spx->status())
      {
      case LPSolver::INFEASIBLE:
         *status = CPX_INFEASIBLE;
         break ;
      case LPSolver::UNBOUNDED:
         *status = CPX_UNBOUNDED;
         break;
      case LPSolver::SOLVED:
         *status = CPX_OPTIMAL;
         break;
      case LPSolver::PRIMAL:
         *status = CPX_IT_LIM_FEAS;
         break;
      case LPSolver::DUAL:
         *status = CPX_TIME_LIM_FEAS;
         break;
      default:
         *status = CPX_ABORT_INFEAS;
         break;
      }
   }

   if (obj != 0)
      *obj = spx->value();

   if (primal != 0)
   {
      Vector tmp(spx->nCols(), primal);
      spx->getPrimalSol(tmp);
   }
   if (redcost != 0)
   {
      Vector tmp(spx->nCols(), redcost);
      spx->getRdCost(tmp);
   }
   if ((slack != 0) || (dual != 0))
   {
      int     rows  = spx->nRows() ;
      Vector tmp(rows);

      if (slack != 0)
      {
         spx->getSlacks(tmp);

         for(int i = 1; i < rows; i++)
         {
            double         x = spx->rhs(i);
            double         y = -tmp[i];
            slack[i - 1] = y + ((x < SPxLP::infinity) ? x : spx->lhs(i));
         }
      }
      if (dual != 0)
      {
         spx->getDualSol(tmp);
         for(int i = 1; i < rows; i++)
            dual[i - 1] = tmp[i];
      }
   }
   return 0;
}

extern "C" int CPXgetnumcols(CPXENVptr env, CPXLPptr lp)
{
   if (env == 0)
      return CPXERR_NO_ENVIRONMENT;
   if (lp == 0)
      return CPXERR_NULL_POINTER;

   return reinterpret_cast<SPxCPlex*>(env)->nCols();
}

extern "C" int CPXgetnumrows(CPXENVptr env, CPXLPptr lp)
{
   if (env == 0)
      return CPXERR_NO_ENVIRONMENT;
   if (lp == 0)
      return CPXERR_NULL_POINTER;

   return reinterpret_cast<SPxCPlex*>(env)->nRows() - 1;
}

extern "C" int CPXgetnumnz(CPXENVptr env, CPXLPptr lp)
{
   if (env == 0)
      return CPXERR_NO_ENVIRONMENT;
   if (lp == 0)
      return CPXERR_NULL_POINTER;

   SPxCPlex* spx = reinterpret_cast<SPxCPlex*>(env);
   int       cnt = 0;

   for(int i = 1; i < spx->nRows(); i++)
      cnt += spx->rowVector(i).size();

   return cnt;
}

extern "C" int CPXgetpahse1cnt(CPXENVptr env, CPXLPptr lp)
{
   if (env == 0)
      return CPXERR_NO_ENVIRONMENT;
   if (lp == 0)
      return CPXERR_NULL_POINTER;

   // SoPlex does not have a "real" Phase I. So we declare
   // one fifth of the iterations as Phase I.
   return reinterpret_cast<SPxCPlex*>(env)->basis().iteration() / 5;
}

extern "C" int CPXgetitcnt(CPXENVptr env, CPXLPptr lp)
{
   if (env == 0)
      return CPXERR_NO_ENVIRONMENT;
   if (lp == 0)
      return CPXERR_NULL_POINTER;

   return reinterpret_cast<SPxCPlex*>(env)->basis().iteration();
}

extern "C" int CPXgetlb(
   CPXENVptr env, CPXLPptr lp, double* low, int start, int end)
{
   if (env == 0)
      return CPXERR_NO_ENVIRONMENT;
   if ((lp == 0) || (low == 0))
      return CPXERR_NULL_POINTER;

   SPxCPlex* spx = reinterpret_cast<SPxCPlex*>(env);

   if ((start < 0) || (end < 0) || (start > end) || (end >= spx->nCols()))
      return CPXERR_BAD_ARGUMENT;

   for(int i = start; i <= end; i++)
      low[i - start] = spx->lower(i);

   return 0;
}

extern "C" int CPXgetub(
   CPXENVptr env, CPXLPptr lp, double* up, int start, int end)
{
   if (env == 0)
      return CPXERR_NO_ENVIRONMENT;
   if ((lp == 0) || (up == 0))
      return CPXERR_NULL_POINTER;

   SPxCPlex* spx = reinterpret_cast<SPxCPlex*>(env);

   if ((start < 0) || (end < 0) || (start > end) || (end >= spx->nCols()))
      return CPXERR_BAD_ARGUMENT;

   for(int i = start; i <= end; i++)
      up[i - start] = spx->upper(i);

   return 0;
}

extern "C" int CPXchgbds(
   CPXENVptr env, CPXLPptr lp, int cnt, int* idx, char* lu, double* bd)
{
   if (env == 0)
      return CPXERR_NO_ENVIRONMENT;
   if ((lp == 0) || (idx == 0) || (lu == 0) || (bd == 0))
      return CPXERR_NULL_POINTER;

   SPxCPlex* spx = reinterpret_cast<SPxCPlex*>(env);

   for(int i = 0; i < cnt; i++)
   {
      if ((idx[i] < 0) || (idx[i] >= spx->nCols()))
         return CPXERR_BAD_ARGUMENT;

      switch(lu[i] )
      {
      case 'U':
         spx->changeUpper(idx[i], bd[i]);
         break ;
      case 'L':
         spx->changeLower(idx[i], bd[i]);
         break ;
      case 'B':
         spx->changeUpper(idx[i], bd[i]);
         spx->changeLower(idx[i], bd[i]);
         break ;
      default:
         return CPXERR_BAD_ARGUMENT;
      }
   }
   return 0 ;
}

/**@todo I suspect, that the case ccnt > 0 is not handled well.
 */
extern "C" int CPXaddrows(
   CPXENVptr env,
   CPXLPptr  lp,
   int       ccnt,
   int       rcnt,
   int       nzcnt,
   double*   rhs,
   char*     sns,
   int*      beg,
   int*      ind,
   double*   val,
   char**    cnames,
   char**    rnames)
{
   if (env == 0)
      return CPXERR_NO_ENVIRONMENT;
   if ((lp == 0) || (rhs == 0) || (sns == 0)
      || (beg == 0) || (ind == 0) || (val == 0))
      return CPXERR_NULL_POINTER;
   if ((ccnt < 0) || (rcnt < 1) || (nzcnt < 0))
      return CPXERR_BAD_ARGUMENT;

   assert(ccnt == 0 && "ccnt > 0 not implemented yet.");

   SPxCPlex* spx = reinterpret_cast<SPxCPlex*>(env);

   LPRowSet rset(rcnt, nzcnt);
   DSVector row(spx->nCols());

   for(int i = 0; i < rcnt; i++)
   {
      row.clear();

      int end = (i < rcnt - 1) ? beg[i + 1] : nzcnt;

      for(int k = beg[i]; k < end; k++)
         row.add(ind[k], val[k]);

      switch(sns[i])
      {
      case 'L':
         rset.add(-SPxLP::infinity, row, rhs[i]);
         break;
      case 'G':
         rset.add(rhs[i], row, SPxLP::infinity);
         break;
      case 'E':
      case 'R':
         rset.add(rhs[i], row, rhs[i]);
         break;
      default:
         return CPXERR_BAD_ARGUMENT;
      }
   }
   spx->addRows(rset);

   return 0;
}

extern "C" int CPXgetrows(
   CPXENVptr env,
   CPXLPptr  lp,
   int*      nzcnt,
   int*      beg,
   int*      ind,
   double*   val,
   int       size,
   int*      surplus,
   int       begin,
   int       end)
{
   if (env == 0)
      return CPXERR_NO_ENVIRONMENT;
   if (lp == 0)
      return CPXERR_NULL_POINTER;
   if ((size > 0) && ((beg == 0) || (ind == 0) || (val == 0)))
      return CPXERR_NULL_POINTER;

   SPxCPlex* spx = reinterpret_cast<SPxCPlex*>(env);

   if ((size < 0) || (begin < 0) || (end < begin) || (end >= spx->nRows()))
      return CPXERR_BAD_ARGUMENT;

   int cnt = 0;
   int i;

   for(i = begin; i <= end; i++)
   {
      const SVector& row = spx->rowVector(i + 1); // row 0 contains obj

      if (size - cnt < row.size())
         break;

      beg[i - begin] = cnt;

      for(int j = 0; j < row.size(); j++)
      {
         ind[cnt] = row.index(j);
         val[cnt] = row.value(j);
         cnt++;
      }
   }
   if (nzcnt)
      *nzcnt = cnt;

    if (i > end)
    {
       if (surplus)
          *surplus = size - cnt;
       return 0;
    }

    if (surplus)
       for(*surplus = 0; i <= end; i++)
          *surplus -= spx->rowVector(i + 1).size();

    return CPXERR_NEGATIVE_SURPLUS;
}

extern "C" int CPXdelsetrows(CPXENVptr env, CPXLPptr lp, int* delstat)
{
   if (env == 0)
      return CPXERR_NO_ENVIRONMENT;
   if ((lp == 0) || (delstat == 0))
      return CPXERR_NULL_POINTER;

   SPxCPlex* spx = reinterpret_cast<SPxCPlex*>(env);

   int            rnum = spx->nRows();
   DataArray<int> rem (rnum);
   DataArray<int> perm(rnum);

   for(int i = 0; i < rnum - 1; i++)
      rem[i + 1] = (delstat[i] == 1) ? -1 : 0;
   rem[0] = 0;

   spx->removeRows(rem.get_ptr(), rnum, perm.get_ptr());

   for(int i = 0; i < rnum - 1; i++)
      delstat[i] = perm[i + 1];

   return 0;
}

extern "C" int CPXsetintparam(CPXENVptr env, int whichparam, int newvalue)
{
   if (env == 0)
      return CPXERR_NO_ENVIRONMENT;

   SPxCPlex* spx = reinterpret_cast<SPxCPlex*>(env);

   switch(whichparam)
   {
   case CPX_PARAM_ADVIND :
      break;
   case CPX_PARAM_ITLIM :
      break;
   case CPX_PARAM_FASTMIP :
      break;
   case CPX_DPRIIND_FULL :
      break;
   case CPX_DPRIIND_STEEP :
      break;
   case CPX_DPRIIND_STEEPQSTART :
      break;
   case CPX_PARAM_SIMDISPLAY :
      break;
   case CPX_PARAM_SCRIND :
      spx->setVerbose(newvalue == CPX_ON);
      break;
   default :
      break;
   }
   return 0;
}

extern "C" int CPXsetdblparam(CPXENVptr env, int whichparam, double newvalue)
{
   if (env == 0)
      return CPXERR_NO_ENVIRONMENT;

   SPxCPlex* spx = reinterpret_cast<SPxCPlex*>(env);

   switch(whichparam)
   {
   case CPX_PARAM_EPRHS :
      break;
   case CPX_PARAM_OBJLLIM :
      break;
   case CPX_PARAM_OBJULIM :
      break;
   default :
      break;
   }
   return 0;
}

//============================================================================

/*ARGSUSED*/
extern "C" int CPXsetlogfile(CPXENVptr env, CPXFILEptr fp)
{
   if (env == 0)
      return CPXERR_NO_ENVIRONMENT;

   return 0;
}

extern "C" int CPXgetsense(
   CPXENVptr env, CPXLPptr lp, char *sns, int start, int end)
{
   if (env == 0)
      return CPXERR_NO_ENVIRONMENT;
   if ((lp == 0) || (sns == 0))
      return CPXERR_NULL_POINTER;

   SPxCPlex* spx = reinterpret_cast<SPxCPlex*>(env);

   if ((start < 0) || (start > end) || (end > spx->nRows()))
      return CPXERR_BAD_ARGUMENT;

   for(int i = start; i <= end; i++)
   {
      if (spx->rhs(i + 1) >= SPxLP::infinity)
         sns[i - start] = 'G' ;
      else if (spx->lhs(i + 1) <= -SPxLP::infinity)
         sns[i - start] = 'L';
      else if (spx->lhs(i + 1) == spx->rhs(i + 1))
         sns[i - start] = 'E';
      else
         sns[i - start] = 'R';
   }
   return 0;
}

extern "C" int CPXgetobjsen(CPXENVptr env, CPXLPptr lp)
{
   if (env == 0)
      return CPXERR_NO_ENVIRONMENT;
   if (lp == 0)
      return CPXERR_NULL_POINTER;

   SPxCPlex* spx = reinterpret_cast<SPxCPlex*>(env);

   return (spx->sense() == LPSolver::MAXIMIZE) ? -1 : 1;
}

/**@todo No check is made if a solution is available.
 */
extern "C" int CPXgetobjval(CPXENVptr env, CPXLPptr lp, double* objval)
{
   if (env == 0)
      return CPXERR_NO_ENVIRONMENT;
   if (lp == 0)
      return CPXERR_NULL_POINTER;

   *objval = reinterpret_cast<SPxCPlex*>(env)->value();

   return 0;
}

extern "C" int CPXgetx(
   CPXENVptr env, CPXLPptr lp, double* x, int start, int end)
{
   if (env == 0)
      return CPXERR_NO_ENVIRONMENT;
   if ((lp == 0) || (x == 0))
      return CPXERR_NULL_POINTER;

   SPxCPlex* spx = reinterpret_cast<SPxCPlex*>(env);

   if ((start < 0) || (start > end) || (end > spx->nCols()))
      return CPXERR_BAD_ARGUMENT;

   DVector tmp(spx->nCols());
   spx->getPrimalSol(tmp);

   for(int i = start; i <= end; i++)
      x[i - start] = tmp[i];

   return 0;
}

extern "C" int CPXgetpi(
   CPXENVptr env, CPXLPptr lp, double* pi, int start, int end)
{
   if (env == 0)
      return CPXERR_NO_ENVIRONMENT;
   if ((lp == 0) || (pi == 0))
      return CPXERR_NULL_POINTER;

   SPxCPlex* spx = reinterpret_cast<SPxCPlex*>(env);

   if ((start < 0) || (start > end) || (end > spx->nRows()))
      return CPXERR_BAD_ARGUMENT;

   DVector tmp(spx->nRows());
   spx->getDualSol(tmp);

   for(int i = start; i <= end; i++)
      pi[i - start] = tmp[i + 1];

   return 0;
}

extern "C" int CPXgetslack(
   CPXENVptr env, CPXLPptr lp, double* slack, int start, int end)
{
   if (env == 0)
      return CPXERR_NO_ENVIRONMENT;
   if ((lp == 0) || (slack == 0))
      return CPXERR_NULL_POINTER;

   SPxCPlex* spx = reinterpret_cast<SPxCPlex*>(env);

   if ((start < 0) || (start > end) || (end > spx->nRows()))
      return CPXERR_BAD_ARGUMENT;

   Vector tmp(spx->nRows());
   spx->getSlacks(tmp);

   for(int i = start; i <= end; i++)
   {
      double         x = spx->rhs(i + 1);
      double         y = -tmp[i + 1];
      slack[i - start] = y + ((x < SPxLP::infinity) ? x : spx->lhs(i + 1));
   }
   return 0;
}

extern "C" int CPXgetdj(
   CPXENVptr env, CPXLPptr lp, double* dj, int start, int end)
{
   if (env == 0)
      return CPXERR_NO_ENVIRONMENT;
   if ((lp == 0) || (dj == 0))
      return CPXERR_NULL_POINTER;

   SPxCPlex* spx = reinterpret_cast<SPxCPlex*>(env);

   if ((start < 0) || (start > end) || (end > spx->nCols()))
      return CPXERR_BAD_ARGUMENT;

   Vector tmp(spx->nCols());
   spx->getRdCost(tmp);

   for(int i = start; i <= end; i++)
      dj[i - start] = tmp[i];

   return 0;
}

extern "C" int CPXgetrhs(
   CPXENVptr env, CPXLPptr lp, double* rhs, int start, int end)
{
   if (env == 0)
      return CPXERR_NO_ENVIRONMENT;
   if ((lp == 0) || (rhs == 0))
      return CPXERR_NULL_POINTER;

   SPxCPlex* spx = reinterpret_cast<SPxCPlex*>(env);

   if ((start < 0) || (start > end) || (end > spx->nCols()))
      return CPXERR_BAD_ARGUMENT;

   for(int i = start; i <= end; i++)
      rhs[i-start] = (spx->rhs(i + 1) >= SPxLP::infinity)
         ? spx->lhs(i + 1) : spx->rhs(i + 1);

   return 0;
}

extern "C" int CPXgetobj(
   CPXENVptr env, CPXLPptr lp, double* obj, int start, int end)
{
   if (env == 0)
      return CPXERR_NO_ENVIRONMENT;
   if ((lp == 0) || (obj == 0))
      return CPXERR_NULL_POINTER;

   SPxCPlex* spx = reinterpret_cast<SPxCPlex*>(env);

   if ((start < 0) || (start > end) || (end > spx->nCols()))
      return CPXERR_BAD_ARGUMENT;

   for(int i = start; i <= end; i++)
      obj[i - start] = spx->obj(i);

    return 0 ;
}

extern "C" int CPXchgobj(
   CPXENVptr env, CPXLPptr lp, int cnt, int* idx, double* obj)
{
   if (env == 0)
      return CPXERR_NO_ENVIRONMENT;
   if ((lp == 0) || (idx == 0) || (obj == 0))
      return CPXERR_NULL_POINTER;

   SPxCPlex* spx = reinterpret_cast<SPxCPlex*>(env);

   while(--cnt >= 0)
      spx->changeObj(idx[cnt], obj[cnt]);

   return 0;
}

extern "C" int CPXdelrows(CPXENVptr env, CPXLPptr lp, int start, int end)
{
   if (env == 0)
      return CPXERR_NO_ENVIRONMENT;
   if (lp == 0)
      return CPXERR_NULL_POINTER;

   SPxCPlex* spx = reinterpret_cast<SPxCPlex*>(env);

   if ((start < 0) || (start > end) || (end > spx->nRows()))
      return CPXERR_BAD_ARGUMENT;

    DataArray<int> del(spx->nRows());

    // obj is 0-th row!!
    start++;
    end++;

    int i = 0;

    while(i < start)
       del[i++] = 0;
    while(i <= end)
       del[i++] = -1;
    while(i < del.size())
       del[i++] = 0;

    spx->removeRows(del.get_ptr());

    return 0;
}

#if 0

extern "C" int CPXaddcols(
   CPXENVptr env,
   CPXLPptr cplex,
   int ccnt,
   int nzcnt,
   double* objx,
   int *cmatbeg,
   int *cmatind,
   double *cmatval,
   double *bdl,
   double *bdu,
   char **cname)
{
    int                i ;
    LPColSet        cset(ccnt, nzcnt) ;
    SPxCPLEX*        optimize = (SPxCPLEX*)cplex ;

    DSVector        col( optimize->nRows() ) ;

    for( i = 0 ; i < ccnt ; ++i )
    {
        int        len ;
        if( i < ccnt-1 )
            len = cmatbeg[i+1] - cmatbeg[i] ;
        else
            len = nzcnt - cmatbeg[i] ;

        col.clear() ;
        double*        _val = &cmatval[cmatbeg[i]] ;
        int*        _idx = &cmatind[cmatbeg[i]] ;
        while( len-- )
            col.add( 1+(*_idx++), *_val++ ) ;

        if( objx[i] != 0 )
            col.add( 0, objx[i] ) ;

        cset.add(objx[i], bdl[i], col, bdu[i]) ;
    }

    optimize->addCols( cset ) ;
    return 0 ;
}

int delcols (CPXLPptr cplex, int begin, int end )
{
    int        i ;
    SPxCPLEX*                optimize = (SPxCPLEX*)cplex ;
    DataArray<int>        del( optimize->nRows() ) ;

    for( i = del.size()-1 ; i >= 0 ; --i )
        del[i] = 0 ;
    for( i = begin ; i <= end ; ++i )
        del[i] = -1 ;
    optimize->removeCols( (int*)del ) ;

    return 0 ;
}

int optimize (CPXLPptr p)
{
    return dualopt(p) ;
}


int delsetcols (CPXLPptr cplex, int *del)
{
#ifdef        DEBUG
    cout << "SPxCPLEX:        delsetcols\n" ;
#endif

    int                i ;
    SPxCPLEX*        optimize = (SPxCPLEX*)cplex ;

    for( i = optimize->nCols()-1 ; i >= 0 ; --i )
        if( del[i] )
            del[i] = -1 ;

    optimize->removeCols( del ) ;

    for( i = optimize->nCols()-1 ; i >= 0 ; --i )
        del[i] = ( del[i] < 0 ) ? 1 : 0 ;

    return 0 ;
}

int lpwrite   (CPXLPptr cplex, char *name)
{
#ifdef        DEBUG
    cout << "SPxCPLEX:        lpwrite\n" ;
#endif

    SPxCPLEX*        optimize = (SPxCPLEX*)cplex ;
    ofstream        out( name ) ;
    if( out.good() )
    {
        out << *optimize ;
        return 0 ;
    }

    cout << endl << "ERROR SPxCPLEX: " << __FILE__ << ":" << __LINE__
         << endl << endl;
    return 1 ;
}


int settilim( double time, double *ptoosmall, double *ptoobig )
{
    SPxCPLEX::maxTime = time ;
    return 0 ;
}

void   gettilim    (double *t)
{
    *t = SPxCPLEX::maxTime ;
}


int        getmethod   (CPXLPptr lp)
{
    assert(lp);
    // provisional: return CPXALG_PRIMAL unless solution() doesn't swap
    // CPX_INFEASIBLE <-> CPX_UNBOUNDED after calling dualopt():
    return CPXALG_PRIMAL;
    //return CPXALG_DUAL;
}


int        setitlim    (int lim, int *x, int *y)
{
    SPxCPLEX::maxIter = lim ;
    return 0 ;
}

void        getitlim    (int *lim )
{
    *lim = SPxCPLEX::maxIter ;
}


int getcols (CPXLPptr lp, int *nzcnt, int *cmatbeg, int *cmatind,
             double *cmatval, int cmatsz, int *surplus, int begin,
             int end)
{
    int                        i, j, n ;
    const SPxCPLEX*        solver = (const SPxCPLEX*)lp ;

#ifdef        DEBUG
    cout << "SPxCPLEX:        getcols\n" ;
#endif

    for( i = begin, n = 0 ; i <= end ; ++i )
    {
        const SVector&        col = solver->colVector(i) ;
        if( cmatsz - n < col.size() )
        {
            if( surplus )
            {
                *surplus = 0 ;
                for( ; i <= end ; ++i )
                    *surplus -= solver->colVector(i).size() ;
            }
            return 0 ;
        }
        cmatbeg[i-begin] = n ;
        for( j = 0 ; j < col.size() ; ++j )
        {
            if( col.index(j) > 0 )                        // row 0 is objective
            {
                cmatind[n] = col.index(j) - 1 ;
                cmatval[n] = col.value(j) ;
                n++ ;
            }
        }
    }

    if( nzcnt )
        *nzcnt = n ;
    if( surplus )
        *surplus = cmatsz - n ;

    return 1 ;
}


#endif
