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
#pragma ident "@(#) $Id: sipinput.cpp,v 1.2 2002/01/23 17:47:01 bzfkocht Exp $"

#include <fstream>
#include <stdio.h>
#include <string.h>
#include <assert.h>

#include "spxlp.h"
#include "nameset.h"
#include "didxset.h"

extern "C"
{
#include "s_def.h"
#include "s_io.h"
}

using namespace soplex;

/*ARGSUSED*/
extern "C" int SIPinput(
   FILE*         ferr, 
   SIPInfaIO     /*infaIO*/, 
   const char*   filename,
   int*          ncol, 
   int*          nrow, 
   int*          objsen,
   double**      obj, 
   double**      rhs, 
   char**        sen, 
   int**         beg, 
   int**         cnt, 
   int**         ind, 
   double**      val, 
   double**      lb, 
   double**      ub, 
   double**      rng,
   char***       cname, 
   char**        cstore, 
   char***       rname, 
   char**        rstore,
   unsigned int* cstoresz, 
   unsigned int* rstoresz, 
   char**        ctype)
{
   assert(ferr     != 0);
   assert(filename != 0);
   assert(ncol     != 0);
   assert(nrow     != 0);
   assert(objsen   != 0);
   assert(obj      != 0);
   assert(rhs      != 0);
   assert(sen      != 0);
   assert(beg      != 0);
   assert(cnt      != 0);
   assert(ind      != 0);
   assert(val      != 0);
   assert(lb       != 0);
   assert(ub       != 0);
   assert(rng      != 0);
   assert(cname    != 0);
   assert(cstore   != 0);
   assert(rname    != 0);
   assert(rstore   != 0);
   assert(cstoresz != 0);
   assert(rstoresz != 0);
   assert(ctype    != 0);

   std::ifstream file(filename);
   SPxLP         lp;
   NameSet       colnames;
   NameSet       rownames;
   DIdxSet       intvars;
   int           i;

   if (!file || !lp.read(file, &rownames, &colnames, &intvars))
   {
      fprintf(ferr, "Could not read file %s\n", filename);
      return SIP_NOFILE;
   }
   rownames.memPack();
   colnames.memPack();

   *ncol   = lp.nCols();
   *nrow   = lp.nRows();
   *objsen = lp.spxSense() == SPxLP::MAXIMIZE ? -1 : 1;  

   if (0 == (*obj = static_cast<double*>(malloc(*ncol * sizeof(**obj)))))
      return SIP_NOMEMORY;
   else
   {
      Vector tmpobj(*ncol, *obj);
      tmpobj = lp.maxObj() * lp.spxSense();
   }

   if (0 == (*rhs = static_cast<double*>(malloc(*nrow * sizeof(**rhs)))))
      return SIP_NOMEMORY;

   if (0 == (*rng = static_cast<double*>(malloc(*nrow * sizeof(**rng)))))
      return SIP_NOMEMORY;

   if (0 == (*sen = static_cast<char*>(malloc(*nrow * sizeof(**sen)))))
      return SIP_NOMEMORY;

   for(i = 0; i < *nrow; i++)
   {
      switch(lp.rowType(i))
      {
      case LPRow::LESS_EQUAL :
         (*sen)[i] = 'L';
         (*rhs)[i] = lp.rhs(i);
         (*rng)[i] = 0.0;
         break; 
      case LPRow::EQUAL :
         (*sen)[i] = 'E';
         (*rhs)[i] = lp.rhs(i);
         (*rng)[i] = 0.0;
         break; 
      case LPRow::GREATER_EQUAL :
         (*sen)[i] = 'G';
         (*rhs)[i] = lp.lhs(i);
         (*rng)[i] = 0.0;
         break; 
      case LPRow::RANGE :
         (*sen)[i] = 'R';
         (*rhs)[i] = lp.rhs(i);
         (*rng)[i] = lp.rhs(i) - lp.lhs(i);
         break; 
      default :
         abort();
      }
   }   
   
   if (  (0 != (*beg = static_cast<int*>(malloc(*ncol * sizeof(**beg)))))
      && (0 != (*cnt = static_cast<int*>(malloc(*ncol * sizeof(**cnt))))))
      return SIP_NOMEMORY;

   int nzo = 0;
   int k   = 0;

   for(i = 0; i < *ncol; i++)
      nzo += lp.colVector(i).size();

   if (  (0 != (*ind = static_cast<int*>(malloc(nzo * sizeof(**ind))))
      && (0 != (*val = static_cast<double*>(malloc(nzo * sizeof(**val)))))))
      return SIP_NOMEMORY;

   for(i = 0; i < *ncol; i++)
   {
      SVector col = lp.colVector(i);
      
      (*beg)[i] = k;
      (*cnt)[i] = col.size();            
      (*ind)[k] = col.index(i);
      (*val)[k] = col.value(i);
      k++;
   }
   assert(k == nzo);

   if (0 == (*lb = static_cast<double*>(malloc(*ncol * sizeof(**lb)))))
      return SIP_NOMEMORY;

   memcpy(lb, lp.lower().get_const_ptr(), *ncol);

   if (0 == (*ub = static_cast<double*>(malloc(*ncol * sizeof(**ub)))))
      return SIP_NOMEMORY;

   memcpy(ub, lp.upper().get_const_ptr(), *ncol);

   assert(colnames.num() == *ncol);

   *cstoresz = colnames.memSize();

   if (0 == (*cstore = 
      static_cast<char*>(malloc(*cstoresz * sizeof(**cstore)))))
      return SIP_NOMEMORY;

   if (0 == (*cname = static_cast<char**>(malloc(*ncol * sizeof(**cname)))))
      return SIP_NOMEMORY;
   
   k = 0;

   for(i = 0; i < *ncol; i++)
   {
      (*cname)[i] = &((*cstore)[k]);
      strcpy((*cname)[i], colnames[i]);
      k += strlen(colnames[i]) + 1;
      assert(static_cast<unsigned int>(k) < *cstoresz);
   }

   assert(rownames.num() == *nrow);

   *rstoresz = rownames.memSize();

   if (0 == (*rstore = 
      static_cast<char*>(malloc(*rstoresz * sizeof(**rstore)))))
      return SIP_NOMEMORY;

   if (0 == (*rname = static_cast<char**>(malloc(*nrow * sizeof(**rname)))))
      return SIP_NOMEMORY;
   
   k = 0;

   for(i = 0; i < *nrow; i++)
   {
      (*rname)[i] = &((*rstore)[k]);
      strcpy((*rname)[i], rownames[i]);
      k += strlen(rownames[i]) + 1;
      assert(static_cast<unsigned int>(k) < *rstoresz);
   }

   if (0 == (*ctype = static_cast<char*>(malloc(*ncol * sizeof(**ctype)))))
      return SIP_NOMEMORY;

   for(i = 0; i < *ncol; i++)
      (*ctype)[i] = 'C';

   for(i = 0; i < intvars.size(); i++)
   {
      int idx = intvars.index(i);

      (*ctype)[i] = 
         (lp.lower(idx) == 0.0 && lp.upper(idx) == 1.0) ? 'B' : 'I';
   }
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










