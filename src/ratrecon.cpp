/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the class library                   */
/*       SoPlex --- the Sequential object-oriented simPlex.                  */
/*                                                                           */
/*    Copyright (C) 1996-2013 Konrad-Zuse-Zentrum                            */
/*                            fuer Informationstechnik Berlin                */
/*                                                                           */
/*  SoPlex is distributed under the terms of the ZIB Academic Licence.       */
/*                                                                           */
/*  You should have received a copy of the ZIB Academic License              */
/*  along with SoPlex; see the file COPYING. If not email to soplex@zib.de.  */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

#ifdef SOPLEX_WITH_GMP

#include <iostream>
#include <assert.h>

#include "ratrecon.h"

#define BITMAX 20
#define EPS 0.000001

namespace soplex
{
   void GMPqv_init(mpq_t * vect,int dim);
   void GMPqv_clear(mpq_t *vect,int dim);
   void GMPv_init(mpz_t * vect,int dim);
   void GMPv_clear(mpz_t *vect,int dim);
   int Reconstruct(mpq_t *x,mpz_t * xnum, mpz_t denom,int dim);

   bool reconstructSol(SolRational& solution)
   {
      DVectorRational buffer;

      if( solution.hasPrimal() )
      {
         buffer.reDim((solution._primal).dim());
         solution.getPrimal(buffer);
         reconstructVector(buffer);
         solution._primal = buffer;

         buffer.reDim((solution._slacks).dim());
         solution.getSlacks(buffer);
         reconstructVector(buffer);
         solution._slacks = buffer;
      }
      if( solution.hasPrimalray() )
      {
         buffer.reDim((solution._primalray).dim());
         solution.getPrimalray(buffer);
         reconstructVector(buffer);
         solution._primalray = buffer;
      }
      if( solution.hasDual() )
      {
         buffer.reDim((solution._dual).dim());
         solution.getDual(buffer);
         reconstructVector(buffer);
         solution._dual = buffer;

         buffer.reDim((solution._redcost).dim());
         solution.getRedcost(buffer);
         reconstructVector(buffer);
         solution._redcost = buffer;
      }
      if( solution.hasDualfarkas() )
      {
         buffer.reDim((solution._dualfarkas).dim());
         solution.getDualfarkas(buffer);
         reconstructVector(buffer);
         solution._dualfarkas = buffer;
      }
      return true;
   }


   bool reconstructVector(VectorBase<Rational>& input)
   {
      int i, dim, rval = true;
      mpq_t * resvec; //reconstructed vector storage
      mpz_t * xnum; //numerator of input vector
      mpz_t denom; //common denominator of input vector

      dim = input.dim();

      // convert vector to mpz format
      resvec = (mpq_t*) malloc( dim * sizeof(mpq_t));
      xnum = (mpz_t*) malloc( dim * sizeof(mpz_t));
      GMPqv_init(resvec,dim);
      GMPv_init(xnum,dim);
      mpz_init_set_ui(denom, 1);

      // find common denominator
      for( i = 0; i < dim; i++ )
         mpz_lcm(denom,denom,mpq_class(input[i]).get_den_mpz_t());

      for( i = 0; i < dim; i++ )
      {
         mpz_mul(xnum[i], denom, mpq_class(input[i]).get_num_mpz_t());
         mpz_divexact(xnum[i], xnum[i], mpq_class(input[i]).get_den_mpz_t());
      }

      //reconstruct
      rval = Reconstruct(resvec,xnum,denom,dim);

      Rational tmprat;
      for( i = 0; i < dim; i++ )
      {
         tmprat = mpq_class(resvec[i]);
         input[i]=tmprat;
      }

      mpz_clear(denom);
      GMPv_clear(xnum,dim);
      GMPqv_clear(resvec,dim);
      return rval;
   }


   int Reconstruct(mpq_t *x,mpz_t * xnum, mpz_t denom,int dim)
   {
      /* This reconstruction routine will set x equal to the mpq vector where */
      /* each component is the best rational approximation of xnum / denom */
      /* with denominator no larger than bound */

      mpz_t delta;
      mpz_init_set_ui(delta,1);
      mpz_t xi;
      mpz_init(xi);

      mpz_t vmod;	/* this is the variable modulous that changes with delta */
      mpz_init(vmod);

      int i;
      int done,skip,rrclose;
      int bitcnt,step;
      int verbose = 0;

      mpz_t dbound;
      mpz_init(dbound);

      mpz_t bound;
      mpz_init(bound);
      mpz_div_ui(bound,bound,2);
      mpz_sqrt(bound,denom);

      mpz_t a,b,c,d,q,r0,r1,anext,bnext,cnext,dnext,r0next,r1next,temp;
      mpz_init(a);
      mpz_init(b);
      mpz_init(c);
      mpz_init(d);
      mpz_init(q);
      mpz_init(r0);
      mpz_init(r1);
      mpz_init(anext);
      mpz_init(bnext);
      mpz_init(cnext);
      mpz_init(dnext);
      mpz_init(r0next);
      mpz_init(r1next);
      mpz_init(temp);
      mpq_t qtemp;
      mpq_init(qtemp);

      double Df,Dq,Dtemp;
      int Da,Db,Dc,Dd,Danext,Dcnext;

      for(i =0; i<dim;i++)
      {
         mpz_set_ui(delta,1);
         /* if x is in range, assign it */

         if(mpz_sgn(xnum[i])==0)
         {
            mpq_set_num(x[i],xnum[i]);
            mpq_set_den(x[i],delta);
            mpq_canonicalize(x[i]);
         }
         else if(mpz_cmp(delta,bound)>0)
         {
            mpz_mul(xi,xnum[i],delta);
            mpz_fdiv_qr(xi,temp,xi,denom);
            if(mpz_cmp_ui(temp,0.5)>0)
               mpz_add_ui(xi,xi,1);
            /* get xi to be rrclosest integer to xnum[i]*delta/denom */

            mpq_set_num(x[i],xi);
            mpq_set_den(x[i],delta);
            mpq_canonicalize(x[i]);
         }
         else
         {
            mpz_set(xi,xnum[i]);
            mpz_mul(xi,xi,delta);

            mpz_cdiv_q(dbound,bound,delta);

            mpz_set(r0,xi);
            mpz_set(r1,denom);

            mpz_gcd(temp,r0,r1);
            mpz_divexact(r0,r0,temp);
            mpz_divexact(r1,r1,temp);

            mpz_fdiv_qr(q,temp,r0,r1);
            mpz_set(r0,r1);
            mpz_set(r1,temp);

            mpz_set_ui(a,1);
            mpz_set_ui(b,0);
            mpz_set(c,q);
            mpz_set_ui(d,1);

            done =0;
            skip = 0;
            rrclose = 0;

            if(mpz_sgn(r1)==0)
            {
               done =1;
            }
            while(!done)
            {

               if( !skip && !rrclose)
               {
                  /* Do a part with doubles */
                  /* Convert r0/r1 to double */
                  mpq_set_num(qtemp,r0);
                  mpq_set_den(qtemp,r1);

                  Df=mpq_get_d(qtemp);
                  bitcnt = 1.44*(log(floor(Df)));

                  /* if bitcnt is small enough do a dbl round */
                  /* otherwise this will send us to an exact step */
                  if(bitcnt >= BITMAX )
                  {
                     /* in this case, get out and do a full precision step */
                     /* printf("bitcnt large, skipping next loop \n"); */
                     skip =1;
                  }
                  else
                  {
                     /*  printf("bitcnt good, doing dbl EEA \n"); */
                     Da=1;Db=0;Dc=0;Dd=1;
                     step =0;

                     /* while bitcount low, do EEA with dbls */
                     while(bitcnt<BITMAX && Df > EPS && step <15)
                     {

                        /*  printf("Doing EEA step \n"); */

                        Dq=floor(Df);
                        /*  printf("Df = %f, Dq = %f\n",Df,Dq); */

                        Dtemp=(Df-Dq);
                        if(Dtemp>EPS)
                        {
                           Df=1/Dtemp;
                           bitcnt += 1.44*(log(Dq))+1;
                           if(bitcnt<BITMAX)
                           {
                              Danext = (int)(Da*Dq+Db+ EPS );
                              Dcnext = (int)(Dc*Dq+Dd+ EPS );
                              Db=Da;
                              Dd=Dc;
                              Da=Danext;
                              Dc=Dcnext;
                              step++;
                           }
                        }
                        else
                           bitcnt= BITMAX +1;
                     }

                     /* update Q and ri using DQ*/
                     /* make Qtemp = Q x DQ */

                     if(Da==1)
                        skip = 1;
                     else
                     {
                        mpz_mul_ui(anext,a,Da);
                        mpz_addmul_ui(anext,b,Dc);

                        if(mpz_cmp(anext,dbound)>0)
                        {
                           /*  printf("next dbound large, moving to all exact\n"); */
                           rrclose = 1; /* rrclose to dbound do rest in full prec */
                        }
                        else
                        {
                           /* printf("updating exact Q\n");	 */
                           /* printf("abcd = %d %d %d %d\n",Da,Db,Dc,Dd);  */
                           /* if it is in correct range assign */

                           mpz_mul_ui(bnext,a,Db);
                           mpz_addmul_ui(bnext,b,Dd);

                           mpz_mul_ui(cnext,c,Da);
                           mpz_addmul_ui(cnext,d,Dc);

                           mpz_mul_ui(dnext,c,Db);
                           mpz_addmul_ui(dnext,d,Dd);

                           mpz_set(a,anext);
                           mpz_set(b,bnext);
                           mpz_set(c,cnext);
                           mpz_set(d,dnext);
                           /* and update remainders as well */

                           mpz_mul_ui(r0next,r0,Dd);
                           mpz_submul_ui(r0next,r1,Db);

                           mpz_mul_ui(r1next,r1,Da);
                           mpz_submul_ui(r1next,r0,Dc);

                           mpz_set(r0,r0next);
                           mpz_set(r1,r1next);

                           if(step%2==1) /* if step number is odd */
                           {
                              mpz_neg(r1,r1);
                              mpz_neg(r0,r0);
                           }
                        }
                     }
                  }

               }
               else
               {
                  /* Do a step in full precision */

                  skip=0; /* tell it to try dbl step again if it can */

                  mpz_fdiv_qr(q,temp,r0,r1);
                  mpz_set(r0,r1);
                  mpz_set(r1,temp);

                  mpz_mul(anext,a,q);
                  mpz_mul(cnext,c,q);
                  mpz_add(anext,anext,b);
                  mpz_add(cnext,cnext,d);
                  step++;

                  if(mpz_cmp(anext,dbound)>0)
                     done =1;
                  else
                  {
                     mpz_set(b,a);
                     mpz_set(d,c);
                     mpz_set(a,anext);
                     mpz_set(c,cnext);
                  }
                  if(mpz_sgn(r1)==0)
                     done = 1;
               }
            }

            /* Assign Values */
            mpq_set_num(x[i],c);
            mpz_mul(delta,delta,a);
            mpq_set_den(x[i],delta);
            mpq_canonicalize(x[i]);
            if(verbose)
            {
               printf("x[%d] set to ",i);
               mpq_out_str(stdout,10,x[i]);
               printf("\n");
            }
         }
      }

      mpz_clear(xi);
      mpz_clear(delta);
      mpz_clear(dbound);
      mpz_clear(bound);
      mpz_clear(a);
      mpz_clear(b);
      mpz_clear(c);
      mpz_clear(d);
      mpz_clear(q);
      mpz_clear(r0);
      mpz_clear(r1);
      mpz_clear(anext);
      mpz_clear(bnext);
      mpz_clear(cnext);
      mpz_clear(dnext);
      mpz_clear(r0next);
      mpz_clear(r1next);
      mpz_clear(temp);
      mpq_clear(qtemp);
      return true;
   }

   void GMPqv_init(mpq_t * vect,int dim)
   {
      int i;
      for(i=0;i<dim;i++)
      {
         mpq_init(vect[i]);
      }
   }

   void GMPqv_clear(mpq_t *vect,int dim)
   {
      int i;
      for(i=0;i<dim ;i++)
      {
         mpq_clear(vect[i]);
      }
   }

   void GMPv_init(mpz_t * vect,int dim)
   {
      int i;
      for(i=0;i<dim;i++)
      {
         mpz_init(vect[i]);
      }
   }

   void GMPv_clear(mpz_t *vect,int dim)
   {
      int i;
      for(i=0;i<dim ;i++)
      {
         mpz_clear(vect[i]);
      }
   }

} // namespace soplex

#endif // SOPLEX_WITH_GMP

//-----------------------------------------------------------------------------
//Emacs Local Variables:
//Emacs mode:c++
//Emacs c-basic-offset:3
//Emacs tab-width:8
//Emacs indent-tabs-mode:nil
//Emacs End:
//-----------------------------------------------------------------------------
