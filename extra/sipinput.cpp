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
#pragma ident "@(#) $Id: sipinput.cpp,v 1.1 2002/01/22 16:48:31 bzfkocht Exp $"

#include <iostream>
#include <fstream>
#include <string.h>
#include <assert.h>

#include "spxlp.h"
#include "nameset.h"
#include "didxset.h"

using namespace soplex;

extern "C" int SIPinput(
   FILE*         ferr, 
   SIPInfaIO     infaIO, 
   char*         filename,
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
   std::ifstream file(filename);
   SPxLP         lp;
   NameSet       colnames;
   NameSet       rownames;
   DIdxSet       intvars;
   int           i;

   if (!file || !lp.read(file, rownames, colnames, intvars))
   {
      fprintf (ferr, "Could not read file %s\n", filename);
      return SIP_NOFILE;
   }
   rownames.memPack();
   colnames.memPack();

   *ncol   = lp.nCols();
   *nrow   = lp.nRows();
   *objsen = lp.spxSense() == MAXIMIZE ? -1 : 1;  

   if (0 == (*obj = malloc(*ncol * sizeof(**obj))))
      return SIP_NOMEMORY;
   else
   {
      Vector tmpobj(*ncol, *obj) = lp.maxObj() * lp.spxSense();
   }

   if (0 == (*rhs = malloc(*nrow * sizeof(**rhs))))
      return SIP_NOMEMORY;

   if (0 == (*rng = malloc(*nrow * sizeof(**rng))))
      return SIP_NOMEMORY;

   if (0 == (*sen = malloc(*nrow * sizeof(**sen))))
      return SIP_NOMEMORY;

   for(i = 0; i < *nrow; i++)
   {
      switch(lp.rowType(i))
      {
      case LESS_EQUAL :
         (*sen)[i] = 'L';
         (*rhs)[i] = lp.rhs(i);
         (*rng)[i] = 0.0;
         break; 
      case EQUAL :
         (*sen)[i] = 'E';
         (*rhs)[i] = lp.rhs(i);
         (*rng)[i] = 0.0;
         break; 
      case GREATER_EQUAL :
         (*sen)[i] = 'G';
         (*rhs)[i] = lp.lhs(i);
         (*rng)[i] = 0.0;
         break; 
      case RANGE :
         (*sen)[i] = 'R';
         (*rhs)[i] = lp.rhs(i);
         (*rng)[i] = lp.rhs(i) - lp.lhs(i);
         break; 
      default :
         abort();
      }
   }   
   
   if (  (0 != (*beg = malloc(*ncol * sizeof(**beg))))
      && (0 != (*cnt = malloc(*ncol * sizeof(**cnt)))))
      return SIP_NOMEMORY;

   int nzo = 0;
   int cnt = 0;
   int i;

   for(i = 0; i < *ncol; i++)
      nzo += lp.colVector(i).size();

   if (  (0 != (*ind = malloc(nzo * sizeof(**ind))))
      && (0 != (*val = malloc(nzo * sizeof(**val)))))
      return SIP_NOMEMORY;

   for(i = 0; i < *ncol; i++)
   {
      SVector col = lp.colVector(i);
      
      (*beg)[i]   = cnt;
      (*cnt)[i]   = col.size();            
      (*ind)[cnt] = col.index(i);
      (*val)[cnt] = col.value(i);
      cnt++;
   }
   assert(cnt == nzo);

   if (0 == (*lb = malloc(*ncol * sizeof(**lb))))
      return SIP_NOMEMORY;
   else
   {
      Vector tmplb(*ncol, *lb) = lp.lower();
   }

   if (0 == (*ub = malloc(*ncol * sizeof(**ub))))
      return SIP_NOMEMORY;
   else
   {
      Vector tmpub(*ncol, *ub) = lp.upper();
   }

   assert(colnames.num() == *ncol);

   *cstoresz = colnames.memSize();

   if (0 == (*cstore = malloc(*cstoresz * sizeof(**cstore))))
      return SIP_NOMEMORY;

   if (0 == (*cname = malloc(*ncol * sizeof(**cname))))
      return SIP_NOMEMORY;
   
   cnt = 0;

   for(i = 0; i < *ncol; i++)
   {
      (*cname)[i] = &((*cstore)[cnt]);
      strcpy((*cname)[i], colnames[i]);
      cnt += strlen(colnames[i]) + 1;
      assert(cnt < *cstoresz);
   }

   assert(rownames.num() == *nrow);

   *rstoresz = rownames.memSize();

   if (0 == (*rstore = malloc(*rstoresz * sizeof(**rstore))))
      return SIP_NOMEMORY;

   if (0 == (*rname = malloc(*nrow * sizeof(**rname))))
      return SIP_NOMEMORY;
   
   cnt = 0;

   for(i = 0; i < *nrow; i++)
   {
      (*rname)[i] = &((*rstore)[cnt]);
      strcpy((*rname)[i], rownames[i]);
      cnt += strlen(rownames[i]) + 1;
      assert(cnt < *rstoresz);
   }

   if (0 == (*ctype = malloc(*ncol * sizeof(**ctype))))
      return SIP_NOMEMORY;

   for(i = 0; i < *ncol; i++)
      (*ctype)[i] = 'C';

   for(i = 0; i < intvars.size(); i++)
   {
      int idx = intvars.index[i];

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










