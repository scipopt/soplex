/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the class library                   */
/*       SoPlex --- the Sequential object-oriented simPlex.                  */
/*                                                                           */
/*    Copyright (C) 2001-2002 Thorsten Koch                                  */
/*                  2001-2002 Konrad-Zuse-Zentrum                            */
/*                            fuer Informationstechnik Berlin                */
/*                                                                           */
/*  SoPlex is distributed under the terms of the ZIB Academic Licence.       */
/*                                                                           */
/*  You should have received a copy of the ZIB Academic License              */
/*  along with SoPlex; see the file COPYING. If not email to soplex@zib.de.  */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
#pragma ident "@(#) $Id: siplp.cpp,v 1.2 2002/01/27 09:32:53 bzfkocht Exp $"

#include <stdio.h>
#include <string.h>
#include <assert.h>

#include "soplex.h"
#include "spxlp.h"
#include "slufactor.h"
#include "spxsteeppr.h"
#include "spxfastrt.h"
#include "nameset.h"
#include "didxset.h"

#define DEBUG  1

#ifdef DEBUG
#define TRACE(x) fprintf(stderr, "%s\n", x)
#else
#define TRACE(x) /**/
#endif

extern "C" 
{
#include "s_lp.h"
}

using namespace soplex;

class SPxSIP : public SoPlex
{
   SLUFactor    m_slu;
   SPxSteepPR   m_price;
   SPxFastRT    m_ratio;
   bool         m_verbose;
   char*        m_probname;
   bool         m_fromscratch;  ///< use old basis indicator
   // NameSet   m_colnames;
   // NameSet   m_rownames;
   // DIdxSet   m_intvars;

protected:
   virtual void factorize()
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
public:
   bool getFromScratch() const
   {
      return m_fromscratch;
   }
   /**@todo If m_fromscratch is true, a new (slack/starter) basis
    *       should be generated.
    */
   void setFromScratch(bool fs)
   {
      m_fromscratch = fs;
   }

   /* Was macht das? Brauchen wir das?
    * void splitLP()
    * {
    *   subcovectors.reSize( 1 ) ;
    * }
    */

   SPxSIP() 
      : SoPlex(LEAVE, COLUMN)
      , m_verbose(true)
      , m_probname(0)
      , m_fromscratch(false)
   {
      setSolver(&m_slu );
      setTester(&m_ratio);
      setPricer(&m_price);
      // no starter, no simplifier
   }
   ~SPxSIP()
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
   void setProbname(const char* probname)
   {
      assert(probname != 0);
      if (m_probname != 0)
         spx_free(m_probname);
      spx_alloc(m_probname, strlen(probname + 1));
      strcpy(m_probname, probname);
   }
   void setVerbose(bool verbose)
   {
      m_verbose = verbose;
   }
};

/*ARGSUSED*/
extern "C" int SIPopenInfa(
   FILE*      /*ferr*/, 
   SIPInfaLP* infaLP, 
   SIPInfaIO* infaIO,
   SIPInfaLP  /*parinfaLP*/)
{
   TRACE("SIPopenInfa");

   assert(infaLP != 0);

   *infaLP = reinterpret_cast<SIPInfaLP>(new SPxSIP);
   *infaIO = 0;
   
   return SIP_OKAY;
      
}

/*ARGSUSED*/
extern "C" int SIPfreeInfa(
   FILE*      /*ferr*/, 
   SIPInfaLP* infaLP, 
   SIPInfaIO* /*infaIO*/,
   SIPInfaLP  /*parinfaLP*/)
{
   TRACE("SIPfreeInfa");

   assert(infaLP  != 0);
   assert(*infaLP != 0);

   delete reinterpret_cast<SPxSIP*>(*infaLP);

   *infaLP = 0;

   return SIP_OKAY;
   
}

/*ARGSUSED*/
extern "C" int SIPcloneinface(
   FILE*      /*ferr*/, 
   SIPInfaLP  /*old*/, 
   SIPInfaLP* newinface)
{
   TRACE("SIPcloneinface");

   assert(newinface != 0);

   *newinface = reinterpret_cast<SIPInfaLP>(new SPxSIP);

   return SIP_OKAY;

}

/*ARGSUSED*/
extern "C" int SIPopenLP(
   FILE*     /*ferr*/, 
   SIPInfaLP infaLP, 
   SIPLP*    lptr, 
   char*     name)
{
   TRACE("SIPopenLP");

   assert(infaLP != 0);
   assert(lptr   != 0);

   SPxSIP* spx = reinterpret_cast<SPxSIP*>(infaLP);
   
   if (name != 0)
      spx->setProbname(name);

   *lptr = reinterpret_cast<SIPLP>(spx);

   return SIP_OKAY;

}

/*ARGSUSED*/
extern "C" int SIPfreeLP(
   FILE*     /*ferr*/, 
   SIPInfaLP infaLP, 
   SIPLP*    lptr)
{
   TRACE("SIPfreeLP");

   assert(infaLP != 0);
   assert(lptr   != 0);
   
   if (*lptr != 0)
      *lptr = 0;

   return SIP_OKAY;
   
}

/*ARGSUSED*/
extern "C" int SIPcopyLP(
   FILE*     ferr, 
   SIPInfaLP infaLP, 
   SIPLP     lptr, 
   int       ncol, 
   int       nrow,
   int       objsen, 
   double*   obj, 
   double*   rhs, 
   char*     sen, 
   int*      beg,
   int*      cnt, 
   int*      ind, 
   double*   val, 
   double*   lb, 
   double*   ub,
   char**    cname, 
   char**    rname)
{
   TRACE("SIPcopyLP");

   assert(ferr   != 0);
   assert(infaLP != 0);
   assert(lptr   != 0);
   assert(ncol   >  0);
   assert(nrow   >  0);
   assert(objsen != 0);
   assert(obj    != 0);
   assert(rhs    != 0);
   assert(sen    != 0);
   assert(beg    != 0);
   assert(cnt    != 0);
   assert(ind    != 0);
   assert(val    != 0);
   assert(lb     != 0);
   assert(ub     != 0);
   assert(cname  != 0);
   assert(rname  != 0);
 
   SPxSIP*   spx = reinterpret_cast<SPxSIP*>(lptr);
   LPColSet  cols(ncol);
   LPRowSet  rows(nrow);
   DSVector  colVector(nrow);
   DSVector  emptyVector(0);
   LPRow     objRow(-SPxLP::infinity, emptyVector, SPxLP::infinity);
   int       i;

   spx->addRow(objRow);

   for(i = 0; i < nrow; i++)
   {
      switch(sen[i])
      {
      case 'L':
         rows.add(-SPxLP::infinity, emptyVector, rhs[i]);
         break;
      case 'G':
         rows.add(rhs[i], emptyVector, SPxLP::infinity);
         break;
      case 'E':
         rows.add(rhs[i], emptyVector, rhs[i]);
         break;
      case 'R':
         fprintf(ferr, "Ranges are not supported\n");
         /* assert(rng != 0);
          * if (rng[i] > 0)
          *    rows.add(rhs[i], emptyVector, rhs[i] + rng[i]);
          * else
          *    rows.add(rhs[i] + rng[i], emptyVector, rhs[i]);
          */
         break ;
      }
   }
   spx->addRows(rows);

   for(i = 0; i < ncol; i++)
   {
      colVector.clear();

      for(int j = 0; j < cnt[i]; j++)
         colVector.add(ind[beg[i] + j] + 1, val[beg[i] + j]);

      colVector.add(0, (objsen == 1) ? -obj[i] : obj[i]);

      cols.add(obj[i], lb[i], colVector, ub[i]);
   }
   spx->changeSense(objsen == 1 ? SPxLP::MINIMIZE : SPxLP::MAXIMIZE);
   spx->addCols(cols);
   
   // we do not use the names for now. They where only needed for the
   // output routines, and those should be part of SIP anyway.

   return SIP_OKAY;
   
}

/*ARGSUSED*/
extern "C" int SIPsetbase(
   FILE*     /*ferr*/, 
   SIPInfaLP infaLP, 
   SIPLP     lptr, 
   double*   /*dnorm*/,
   int*      cstat, 
   int*      rstat, 
   int       /*pricing*/)
{
   TRACE("SIPsetbase");

   assert(infaLP != 0);
   assert(lptr   != 0);
   assert(cstat  != 0);
   assert(rstat  != 0);
    
   SPxSIP* spx = reinterpret_cast<SPxSIP*>(lptr);

   SPxBasis::Desc desc;
   int i;

   desc.reSize(spx->nRows(), spx->nCols());

   if (spx->basis().solver() != spx)
      spx->basis().load(spx);

   // obj limit in row 0!!!
   for(i = 1; i < spx->nRows(); i++)
   {
      switch(rstat[i])
      {
      case 0 : // AT_LOWER
      case 2 : // AT_UPPER
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
      case 1 : // IS_BASIC
         desc.rowStatus(i + 1) = spx->dualRowStatus(i + 1);
         break;
      default:
         abort();
      }
   }
   desc.rowStatus(0) = spx->dualRowStatus(0);

   for(i = 0; i < spx->nCols(); i++)
   {
      switch(cstat[i])
      {
      case 0: // AT_LOWER
         if (spx->upper(i) == spx->lower(i))
            desc.colStatus(i) = SPxBasis::Desc::P_FIXED;
         else
            desc.colStatus(i) = SPxBasis::Desc::P_ON_LOWER;
         break;
      case 2: // AT_UPPER
         if (spx->upper(i) == spx->lower(i))
            desc.colStatus(i) = SPxBasis::Desc::P_FIXED;
         else
            desc.colStatus(i) = SPxBasis::Desc::P_ON_UPPER;
         break;
      case 1:  // IS_BASIC
         desc.colStatus(i) = spx->dualColStatus(i);
         break;
      case 3: // FREE 
         desc.colStatus(i) = SPxBasis::Desc::P_FREE;
         break;
      default:
         abort();
      }
   }  
   spx->loadBasis(desc);

   return SIP_OKAY;
   
}

/*ARGSUSED*/
extern "C" int SIPgetsol(
   FILE*     /*ferr*/, 
   SIPInfaLP infaLP, 
   SIPLP     lptr, 
   int*      solstat, 
   double*   objval,
   double*   objx, 
   double*   pi, 
   double*   slck, 
   double*   redcost)
{
   TRACE("SIPgetsol");

   assert(infaLP != 0);
   assert(lptr   != 0);
    
   SPxSIP* spx = reinterpret_cast<SPxSIP*>(lptr);

   if (solstat != 0)
      *solstat = spx->status();

   if (objval != 0)
      *objval = spx->value();

   if (objx != 0)
   {
      Vector tmp(spx->nCols(), objx);
      spx->getPrimal(tmp);
   }
   if (redcost != 0)
   {
      Vector tmp(spx->nCols(), redcost);
      spx->getRdCost(tmp);
   }
   if ((slck != 0) || (pi != 0))
   {
      int     rows  = spx->nRows() ;
      DVector tmp(rows);

      if (slck != 0)
      {
         spx->getSlacks(tmp);

         for(int i = 1; i < rows; i++)
         {
            double    x = spx->rhs(i);
            double    y = -tmp[i];
            slck[i - 1] = y + ((x < SPxLP::infinity) ? x : spx->lhs(i)); 
         }
      }
      if (pi != 0)
      {
         spx->getDual(tmp);
         for(int i = 1; i < rows; i++)
            pi[i - 1] = tmp[i];
      }
   }
   return SIP_OKAY;

}

/*ARGSUSED*/
extern "C" int SIPgetbase(
   FILE*     /*ferr*/, 
   SIPInfaLP infaLP, 
   SIPLP     lptr,
   double*   /*dnorm*/, 
   int*      cstat, 
   int*      rstat)
{
   TRACE("SIPgetbase");

   assert(infaLP != 0);
   assert(lptr   != 0);
    
   SPxSIP* spx = reinterpret_cast<SPxSIP*>(lptr);

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
   return SIP_OKAY;

}

/*ARGSUSED*/
extern "C" int SIPsetintparLP(
   SIPInfaLP infaLP, 
   int       type, 
   int       ival)
{
   TRACE("SIPsetintparLP");

   assert(infaLP != 0);
    
   SPxSIP* spx = reinterpret_cast<SPxSIP*>(infaLP);
   int restat  = SIP_OKAY;

   switch(type) 
   {
   case SIP_FROMSCRATCH:
      spx->setFromScratch(ival == SIP_ON);
      break;
   case SIP_LPITLIM:
      spx->setTerminationIter(ival);
      break;
   case SIP_FASTMIP:
      break;
   case SIP_PRICING:
      break;
   case SIP_LPINFO:
      break;
   default:
      restat = SIP_LPERROR;
      break;
   }
   return restat;

}/* END SETPARLP */

/*ARGSUSED*/
extern "C" int SIPsetdblparLP(
   SIPInfaLP infaLP, 
   int       type, 
   double    dval)
{
   TRACE("SIPsetdblparLP");

   assert(infaLP != 0);
    
   SPxSIP* spx = reinterpret_cast<SPxSIP*>(infaLP);
   int restat  = SIP_OKAY;

   switch(type) 
   {
   case SIP_FEASTOL:
      spx->setDelta(dval);
      break;
   case SIP_LOBJLIM:
      spx->setTerminationValue(dval);
      break;
   case SIP_UOBJLIM:
      spx->setTerminationValue(dval);
      break;
   case SIP_LPTILIM:
      spx->setTerminationTime(dval);
      break;
   default:
      restat = SIP_LPERROR;
      break;
   }
   return restat;
}

/*ARGSUSED*/
extern "C" int SIPgetintparLP(
   SIPInfaLP infaLP, 
   int       type, 
   SIPLP     /*lptr*/,
   int*      ival)
{
   TRACE("SIPgetintparLP");

   assert(infaLP != 0);
   assert(ival   != 0);

   SPxSIP* spx = reinterpret_cast<SPxSIP*>(infaLP);
   int restat  = SIP_OKAY;

   switch(type) 
   {
   case SIP_FROMSCRATCH:
      *ival = spx->getFromScratch() ? SIP_ON : SIP_OFF;
      break;
   case SIP_LPNROW:      
      *ival = spx->nRows();
      break;
   case SIP_LPIT:
      *ival = spx->basis().iteration();
      break;
   case SIP_LPIT1:
      // There is no phase I in SoPlex, so we guess a value.
      *ival = spx->basis().iteration() / 5;
      break;
   default:
      restat = SIP_LPERROR;
      break;
   }
   return restat;
}

/*ARGSUSED*/
extern "C" int SIPgetdblparLP(
   SIPInfaLP infaLP, 
   int       type, 
   SIPLP     /*lptr*/,
   double*   dval)
{
   TRACE("SIPgetdblparLP");

   assert(infaLP != 0);
   assert(dval   != 0);

   SPxSIP* spx = reinterpret_cast<SPxSIP*>(infaLP);
   int restat  = SIP_OKAY;

   switch(type)
   {
   case SIP_FEASTOL:
      *dval = spx->delta();
      break;
   default:
      restat = SIP_LPERROR;
      break;
   }
   return restat;
}

/* returns TRUE iff LP solution is stable, ie
   no numerical troubles occured */
/*ARGSUSED*/
extern "C" int SIPisStable(int solstat)
{
   TRACE("SIPisStable");

   return (solstat != SoPlex::SINGULAR) && (solstat != SoPlex::ERROR);
}

/* returns TRUE iff LP has been solved by presolve */ 
/*ARGSUSED*/
extern "C" int SIPisPresolved(int /*status*/)
{
   TRACE("SIPisPresolved");
   // There is no presolve in SoPlex
   return 0; 
}

/* returns TRUE iff LP is infeasible */
/*ARGSUSED*/
extern "C" int SIPisInfeas(int status)
{
   TRACE("SIPisInfeas");

   return (status == SoPlex::INFEASIBLE) || (status == SoPlex::UNBOUNDED); 
}

/* returns TRUE iff LP is solved to optimality */
/*ARGSUSED*/
extern "C" int SIPisOptimal(int status)
{
   TRACE("SIPisOptimal");

   return status == SoPlex::OPTIMAL;
}


/* returns TRUE iff a feasible solution has been found */
/*ARGSUSED*/
extern "C" int SIPisValidBound(int status)
{
   TRACE("SIPisValidBound");

   return (status == SoPlex::OPTIMAL)
      ||  (status == SoPlex::ABORT_TIME)
      ||  (status == SoPlex::ABORT_ITER)
      ||  (status == SoPlex::ABORT_VALUE);
}

/* returns TRUE iff objective limit is exceeded */
/*ARGSUSED*/
extern "C" int SIPexObjlim(int status)
{
   TRACE("SIPexObjlim");

   return (status == SoPlex::UNBOUNDED)
      ||  (status == SoPlex::ABORT_VALUE);
}

/* returns TRUE iff error occured during LP solve */ 
/*ARGSUSED*/
extern "C" int SIPerrorLP(int status)
{
   TRACE("SIPerrorLP");

   return ( status < 0 );
}

/* TRUE iff iteration limit has been exceeded */
/*ARGSUSED*/
extern "C" int SIPiterlim(int status)
{
   TRACE("SIPiterlim");

   return (status == SoPlex::ABORT_VALUE);
}

/*ARGSUSED*/
extern "C" int SIPwriteLP(
   FILE*     ferr, 
   SIPInfaLP infaLP, 
   SIPLP     lptr, 
   char*     fname)
{
   TRACE("SIPwriteLP");

   assert(ferr   != 0);
   assert(infaLP != 0);
   assert(lptr   != 0);
   assert(fname  != 0);

   fprintf(ferr, "SIPwriteLP not implemented\n");

   return SIP_LPERROR;
#if 0
   int restat;
   
   restat = CPXlpwrite ( (CPXENVptr) infaLP, (CPXLPptr) lptr, fname);
   if ( restat != 0 ) {
      fprintf (ferr, "Error in CPXlpwrite, %d returned\n", restat);
      return SIP_LPERROR;
   }
   return SIP_OKAY;
#endif
}

/*ARGSUSED*/
extern "C" int SIPwriteB(
   FILE*     ferr, 
   SIPInfaLP infaLP, 
   SIPLP     lptr, 
   char*     fname)
{
   TRACE("SIPwriteB");

   assert(ferr   != 0);
   assert(infaLP != 0);
   assert(lptr   != 0);
   assert(fname  != 0);

   fprintf(ferr, "SIPwriteB not implemented\n");

   return SIP_LPERROR;
#if 0
   int restat;
   
   restat = CPXmbasewrite ((CPXENVptr) infaLP, (CPXLPptr) lptr, fname);
   if ( restat != 0 ) {
      fprintf (ferr, "Error in CPXmbasewrite (), %d returned\n", restat);
      return SIP_LPERROR;
   }

   return SIP_OKAY;
#endif
}

/*ARGSUSED*/
extern "C" int SIPgetlb(
   FILE*     /*ferr*/, 
   SIPInfaLP infaLP, 
   SIPLP     lptr,
   double*   lb, 
   int       beg, 
   int       end)
{
   TRACE("SIPgetlb");

   assert(infaLP != 0);
   assert(lptr   != 0);
   assert(lb     != 0);

   SPxSIP* spx = reinterpret_cast<SPxSIP*>(lptr);

   assert(beg    >= 0);
   assert(end    >= beg);
   assert(end    <  spx->nCols());

   for(int i = beg; i <= end; i++)
      lb[i - beg] = spx->lower(i);

   return SIP_OKAY; 
}

/*ARGSUSED*/
extern "C" int SIPgetub(
   FILE*     /*ferr*/, 
   SIPInfaLP infaLP, 
   SIPLP     lptr,
   double*   ub, 
   int       beg, 
   int       end)
{
   TRACE("SIPgetub");

   assert(infaLP != 0);
   assert(lptr   != 0);
   assert(ub     != 0);

   SPxSIP* spx = reinterpret_cast<SPxSIP*>(lptr);

   assert(beg    >= 0);
   assert(end    >= beg);
   assert(end    <  spx->nCols());
 
   for(int i = beg; i <= end; i++)
      ub[i - beg] = spx->upper(i);

   return SIP_OKAY;
}

/*ARGSUSED*/
extern "C" int SIPchgbds(
   FILE*     /*ferr*/, 
   SIPInfaLP infaLP, 
   SIPLP     lptr,
   int       cnt, 
   int*      ind, 
   char*     lu, 
   double*   bd)
{
   TRACE("SIPchgbds");

   assert(infaLP != 0);
   assert(lptr   != 0);
   assert(cnt    >  0);
   assert(ind    != 0);
   assert(lu     != 0);
   assert(bd     != 0);

   SPxSIP* spx = reinterpret_cast<SPxSIP*>(lptr);

   for(int i = 0; i < cnt; i++)
   {
      assert(ind[i] >= 0);
      assert(ind[i] <  spx->nCols());

      switch(lu[i] )
      {
      case 'U':
         spx->changeUpper(ind[i], bd[i]);
         break ;
      case 'L':
         spx->changeLower(ind[i], bd[i]);
         break ;
      case 'B':
         spx->changeUpper(ind[i], bd[i]);
         spx->changeLower(ind[i], bd[i]);
         break ;
      default:
         abort();
      }
   }
   return SIP_OKAY;
}

/*ARGSUSED*/
extern "C" void SIPchgobjsen(
   FILE*     /*ferr*/, 
   SIPInfaLP infaLP, 
   SIPLP     lptr, 
   int       objsen)
{
   TRACE("SIPchgobjsen");

   assert(infaLP != 0);
   assert(lptr   != 0);
   assert(objsen != 0);

   SPxSIP* spx = reinterpret_cast<SPxSIP*>(lptr);

   spx->changeSense(objsen == 1 ? SoPlex::MINIMIZE : SoPlex::MAXIMIZE);
}

/*ARGSUSED*/
extern "C" int SIPdelrows(
   FILE*     /*ferr*/, 
   SIPInfaLP infaLP, 
   SIPLP     lptr, 
   int*      dstat)
{
   TRACE("SIPdelrows");

   assert(infaLP != 0);
   assert(lptr   != 0);
   assert(dstat  != 0);

   SPxSIP* spx = reinterpret_cast<SPxSIP*>(lptr);

   int            rnum = spx->nRows();
   DataArray<int> rem (rnum);
   DataArray<int> perm(rnum);
   int            i;

   for(i = 0; i < rnum - 1; i++)
      rem[i + 1] = (dstat[i] == 1) ? -1 : 0;
   rem[0] = 0;

   spx->removeRows(rem.get_ptr(), rnum, perm.get_ptr());

   for(i = 0; i < rnum - 1; i++)
      dstat[i] = perm[i + 1];

   return SIP_OKAY;
}

/*ARGSUSED*/         
extern "C" int SIPgetBind(
   FILE*     ferr, 
   SIPInfaLP /*infaLP*/, 
   SIPLP     /*lptr*/,
   int*      /*head*/, 
   double*   /*x*/)
{
   TRACE("SIPgetBind");

   fprintf(ferr, "SIPgetBind not implemented\n");

#if 0
   restat = CPXgetbhead ((CPXENVptr) infaLP, (CPXLPptr) lptr, head, x);
#endif

   return SIP_LPERROR;
   
}/* END GETBIND */

/* returns a row of some matrix:
   SIP_ABm1:  *val contains the ith row of A_B^-1 (non-sparse)
              *ind, *nnonz are not used
   SIP_ABm1A: *val contains the ith row of A_B^-1 A (non-sparse)
              *ind, *nnonz are not used
   SIP_A:     *val and *ind contain the ith row of A (sparse)
              *nnonz contains on input the least size of *val and *ind
              and on output the actual size of *val and *ind 
 */
/*ARGSUSED*/         
extern "C" int SIPgetrow(
   FILE*     ferr, 
   SIPInfaLP infaLP, 
   SIPLP     lptr,
   int       of, 
   int       i, 
   double*   val, 
   int*      ind, 
   int*      nnonz)
{
   TRACE("SIPgetrow");

   assert(ferr   != 0);
   assert(infaLP != 0);
   assert(lptr   != 0);

   SPxSIP* spx = reinterpret_cast<SPxSIP*>(lptr);

   switch(of) 
   {
   case SIP_ABm1:
      fprintf(ferr, "SIPgetrow ABm1 not implemented\n");
      return SIP_LPERROR;
   case SIP_ABm1A:
      fprintf(ferr, "SIPgetrow ABm1A not implemented\n");
      return SIP_LPERROR;
   case SIP_A:
      assert(val   != 0);
      assert(ind   != 0);
      assert(nnonz != 0);
      assert(i     >= 0);
      assert(i     <  spx->nRows() - 1);
      // needed because of creation of Svector
      {
         int cnt            = 0;     
         const SVector& row = spx->rowVector(i + 1); // row 0 contains obj

         assert(*nnonz < row.size());

         for(int j = 0; j < row.size(); j++)
         {
            ind[cnt] = row.index(j);
            val[cnt] = row.value(j);
            cnt++;
         }
         *nnonz = cnt;
      }
      break;
   default:
      fprintf(ferr, "Unknown matrix %d.\n", of); 
      return SIP_LPERROR;
   }
   return SIP_OKAY;

#if 0
   switch ( of ) {
      case SIP_ABm1:
         restat = CPXbinvrow ((CPXENVptr) infaLP, (CPXLPptr) lptr, i, val);
         if ( restat != 0 ) {
            fprintf (ferr, "Error in CPXbinvrow(), %d returned.\n", restat);
            return SIP_LPERROR;
         }
         break;
      case SIP_ABm1A:
         restat = CPXbinvarow ((CPXENVptr) infaLP, (CPXLPptr) lptr, i, val);
         if ( restat != 0 ) {
            fprintf (ferr, "Error in CPXbinarow(), %d returned.\n", restat); 
            return SIP_LPERROR;
         }
         break;
#endif   
}

/*ARGSUSED*/
extern "C" int SIPoptLP(
   SIPInfaLP infaLP, 
   SIPLP     lptr)
{
   TRACE("SIPoptLP");

   assert(infaLP != 0);
   assert(lptr   != 0);

   SPxSIP*        spx  = reinterpret_cast<SPxSIP*>(lptr);
   SoPlex::Status stat = spx->solve();

   return (stat == SoPlex::ERROR) 
      ||  (stat == SoPlex::SINGULAR) 
      ||  (stat == SoPlex::NO_PROBLEM);  
}

/*ARGSUSED*/
extern "C" int SIPstrongbranch(
   FILE*     ferr, 
   SIPInfaLP infaLP, 
   SIPLP     lptr, 
   int*      cand, 
   int       ncand,
   double*   down, 
   double*   up, 
   int       itlim)
{
   TRACE("SIPstrongbranch");

   assert(ferr   != 0);
   assert(infaLP != 0);
   assert(lptr   != 0);
   assert(cand   != 0);
   assert(ncand  >  0);
   assert(down   != 0);
   assert(up     != 0);
   assert(itlim  >  0);

   // SPxSIP* spx = reinterpret_cast<SPxSIP*>(lptr);

   fprintf(ferr, "SIPstrongbranch not yet implemented\n");

   return SIP_LPERROR;
#if 0
   int restat;
   
   restat = CPXstrongbranch ((CPXENVptr) infaLP, (CPXLPptr) lptr, 
                             cand, ncand, down, up, itlim);
   if ( restat != 0 ) {
      fprintf (ferr, "Error in CPXstrongbranch (), %d returned\n", restat);
      return SIP_LPERROR;
   }
   return SIP_OKAY;
#endif   
}

/*ARGSUSED*/
extern "C" int SIPaddrow(
   FILE*     ferr, 
   SIPInfaLP infaLP, 
   SIPLP     lptr, 
   int       rcnt, 
   int       nzcnt,
   double*   rhs, 
   char*     sns, 
   int*      beg, 
   int*      ind, 
   double*   val,
   char**    /*name*/)    // we don't use the names yet.
{
   TRACE("SIPaddrow");

   assert(ferr   != 0);
   assert(infaLP != 0);
   assert(lptr   != 0);
   assert(rhs    != 0);
   assert(sns    != 0);
   assert(beg    != 0);
   assert(ind    != 0);
   assert(val    != 0);
   assert(rcnt   >= 0);
   assert(nzcnt  >= 0);

   SPxSIP* spx = reinterpret_cast<SPxSIP*>(lptr);

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
         abort();
      }
   }
   spx->addRows(rset);

   return SIP_OKAY;
}

//-----------------------------------------------------------------------------
//Emacs Local Variables:
//Emacs mode:c++
//Emacs c-basic-offset:3
//Emacs tab-width:8
//Emacs indent-tabs-mode:nil
//Emacs End:
//-----------------------------------------------------------------------------
