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
#pragma ident "@(#) $Id: cmdline.cpp,v 1.1 2001/11/06 16:18:31 bzfkocht Exp $"

/*      \Section{Complex Methods}
 */

/*  Import system include files
 */
#include <stdlib.h>
#include <string.h>
#include <stdarg.h>
#include <iostream>

/*  ... and class header files
 */
#include "cmdline.h"

namespace soplex
{


//@ ----------------------------------------------------------------------------
/*  \Section{The Parser}
 */
void CmdLine::usageLine(CmdLine::A* arg)
{
   switch (arg->type)
   {
   case A::LIST :
      std::cerr << "       \t<list>\t";
      break;
   case A::SWITCH :
      std::cerr << arg->copt << "  --" << arg->sopt << "\t\t";
      break;
   case A::INT :
      std::cerr << arg->copt << "  --" << arg->sopt
      << "\t<int>\t";
      break;
   case A::DOUBLE :
      std::cerr << arg->copt << "  --" << arg->sopt
      << "\t<real>\t";
      break;
   case A::STRING :
      std::cerr << arg->copt << "  --" << arg->sopt
      << "\t<str>\t";
      break;
   case A::OPT_LIST :
      std::cerr << arg->copt << "  --" << arg->sopt
      << "\t<list>\t";
      break;
   default:
      std::cerr << arg->copt << "  --" << arg->sopt;
      std::cerr << "\t\t";
      break;
   }
}

void CmdLine::usage(
   const char* prog,
   const char* help,
   DataArray < A* > & alist
)
{
   std::cerr << prog << ": " << help << std::endl;
   std::cerr << "USAGE: \n" << std::endl;
   std::cerr << '\t' << "[h  --help" << "\t\tdisplay this help message]" << std::endl;
   for (int i = 1; i < alist.size(); ++i)
   {
      A* arg = (A*)alist[i];
      std::cerr << "\t";
      if (arg->optional)
         std::cerr << '[';
      else
         std::cerr << ' ';
      usageLine(arg);
      std::cerr << arg->hlp;
      if (arg->optional)
         std::cerr << ']' << std::endl;
      else
         std::cerr << std::endl;
   }
}

int CmdLine::doOption(
   int argc,
   char** argv,
   const char help[],
   DataArray < A* > & alist,
   int& n,
   A* alisti,
   const char* start
)
{
   switch (alisti->type)
   {
   case A::SWITCH :
      *alisti->data.onoff = 1;
      break;
   case A::INT :
      if (*start)
         *alisti->data.integer = atoi(start);
      else if (++n < argc)
         *alisti->data.integer = atoi(argv[n]);
      else
      {
         std::cerr << "ERROR: in " << n << "-th option '" << argv[n] << "'\n";
         usage(argv[0], help, alist);
         return n;
      }
      break;
   case A::DOUBLE :
      if (*start)
         *alisti->data.real = atof(start);
      else if (++n < argc)
         *alisti->data.real = atof(argv[n]);
      else
      {
         std::cerr << "ERROR: in " << n << "-th option '" << argv[n] << "'\n";
         usage(argv[0], help, alist);
         return n;
      }
      break;
   case A::STRING :
      if (*start)
         *alisti->data.string = start;
      else if (++n < argc)
         *alisti->data.string = argv[n];
      else
      {
         std::cerr << "ERROR: in " << n << "-th option '" << argv[n] << "'\n";
         usage(argv[0], help, alist);
         return n;
      }
      break;
   case A::OPT_LIST :
      if (++n < argc)
      {
         *alisti->data.list.strings = (const char**) & argv[n];
         for (; n < argc && argv[n][0] != '-'; ++n)
            (*alisti->data.list.num)++;
         if (n < argc)
            --n;
      }
      break;
   default:
      std::cerr << "CmdLine::parser ERROR: unknown option type found for option '"
      << argv[n] << "'\n";
      exit(-1);
   }
   return 0;
}

int CmdLine::parse(int argc, char** argv, const char help[], DataArray < A* > & alist)
{
   /*  Parse command line
    */
   int n, i;
   for (n = 1; n < argc; ++n)
   {
      if (argv[n][0] == '-')
      {
         const char* start = NULL;
         if (argv[n][1] == '-')
         {
            const char* opt = &argv[n][2];
            if (strcmp(opt, "help") == 0)
            {
               usage(argv[0], help, alist);
               return -1;
            }
            for (i = alist.size() - 1; i > 0; --i)
            {
               if (strcmp(((A*)alist[i])->sopt, opt) == 0)
               {
                  start = opt;
                  while (*start) start++;
                  break;
               }
            }
         }
         else
         {
            char opt = argv[n][1];
            if (opt == 'h')
            {
               usage(argv[0], help, alist);
               return -1;
            }
            for (i = alist.size() - 1; i > 0; --i)
            {
               if (((A*)alist[i])->copt == opt)
               {
                  start = &argv[n][2];
                  break;
               }
            }
         }
         if (i <= 0)
         {
            std::cerr << argv[0] << " ERROR: unknown option '" << argv[n] << "'\n";
            usage(argv[0], help, alist);
            return n;
         }
         A* alisti = (A*)alist[i];
         alisti->done = 1;
         if ((i = CmdLine::doOption(argc, argv, help, alist, n, alisti, start)))
            return i;
      }
      else
      {
         for (i = alist.size() - 1; i >= 0; --i)
         {
            A* alisti = (A*)alist[i];
            if (alisti->type == A::LIST)
            {
               *alisti->data.list.strings = (const char**) & argv[n];
               *alisti->data.list.num = argc - n;
               n = argc;
               alisti->done = 1;
               break;
            }
         }
         if (i < 0)
         {
            std::cerr << "ERROR: in " << n << "-th argument '" << argv[n] << "'\n";
            usage(argv[0], help, alist);
            return n;
         }
      }
   }

   /*  Test all required arguments
    */
   int miss = 0;
   for (i = alist.size() - 1; i > 0; --i)
   {
      A* alisti = (A*)alist[i];
      if (alisti->optional == 0 && alisti->done == 0)
      {
         std::cerr << "ERROR: missing argument -" << alisti->copt << std::endl;
         miss = 1;
      }
   }
   if (miss)
   {
      usage(argv[0], help, alist);
      return n;
   }

   return 0;
}


int CmdLine::parse(int argc, char** argv, const char help[], Parameter* args)
{
   DataArray < A* > alist(0, 20, 2);
   int i;

   A help_arg;
   help_arg.copt = 'h';
   help_arg.type = A::SWITCH;
   alist.append(&help_arg);

   if (args)
      for (i = 0; args[i]; ++i)
         alist.append(args[i].arg);

   return parse(argc, argv, help, alist);
}

int CmdLine::parse(int argc, char** argv, const char help[], ...)
{
   va_list ap;
   va_start (ap, help);

   /*  Setup option space
    */
   DataArray < A* > alist(0, 20, 2);
   A* arg;

   A help_arg;
   help_arg.copt = 'h';
   help_arg.type = A::SWITCH;
   alist.append(&help_arg);

   for (;;)
   {
      arg = (A*) va_arg(ap, void*);
      if (arg == 0 || arg->type == A::TERMINATE)
         break;
      for (int i = alist.size() - 1; i >= 0; --i)
      {
         if (((A*)alist[i])->copt == arg->copt)
         {
            std::cerr << "CmdLine::parser ERROR: doubly specified "
            << alist.size() - 1 << "-th option character '"
            << arg->copt << "'\n";
            exit(-1);
         }
      }
      alist.append(arg);
   }

   return parse(argc, argv, help, alist);
}

//@ ----------------------------------------------------------------------------
/*  \SubSection{Parameter construction Methods}
 */
CmdLine::Parameter::Parameter(
   char opt_char,
   const char opt_string[],
   const char help[],
   char& onoff,
   int opt
)
{
   arg = new Arg;
   arg->optional = opt;
   arg->done = 0;
   arg->copt = opt_char;
   arg->sopt = opt_string;
   arg->hlp = help;
   arg->type = Arg::SWITCH;
   arg->data.onoff = &onoff;
   onoff = 0;
}

CmdLine::Parameter::Parameter(
   char opt_char,
   const char* opt_string,
   const char* help,
   int& val,
   int opt
)
{
   arg = new Arg;
   arg->optional = opt;
   arg->done = 0;
   arg->copt = opt_char;
   arg->sopt = opt_string;
   arg->hlp = help;
   arg->type = Arg::INT;
   arg->data.integer = &val;
}

CmdLine::Parameter::Parameter(
   char opt_char,
   const char* opt_string,
   const char* help,
   double& val,
   int opt
)
{
   arg = new Arg;
   arg->optional = opt;
   arg->done = 0;
   arg->copt = opt_char;
   arg->sopt = opt_string;
   arg->hlp = help;
   arg->type = Arg::DOUBLE;
   arg->data.real = &val;
}

CmdLine::Parameter::Parameter(
   char opt_char,
   const char* opt_string,
   const char* help,
   const char*& string,
   int opt
)
{
   arg = new Arg;
   arg->optional = opt;
   arg->done = 0;
   arg->copt = opt_char;
   arg->sopt = opt_string;
   arg->hlp = help;
   arg->type = Arg::STRING;
   arg->data.string = &string;
   string = 0;
}

CmdLine::Parameter::Parameter(
   char opt_char,
   const char* opt_string,
   const char* help,
   int& num,
   const char**& strings,
   int opt
)
{
   arg = new Arg;
   arg->optional = opt;
   arg->done = 0;
   arg->copt = opt_char;
   arg->sopt = opt_string;
   arg->hlp = help;
   arg->type = Arg::OPT_LIST;
   arg->data.list.num = &num;
   arg->data.list.strings = &strings;
   num = 0;
}

CmdLine::Parameter::Parameter(
   const char* help,
   int& num,
   const char**& strings,
   int opt
)
{
   static char none[] = "";
   arg = new Arg;
   arg->optional = opt;
   arg->done = 1;
   arg->copt = 0;
   arg->sopt = none;
   arg->hlp = help;
   arg->type = Arg::LIST;
   arg->data.list.num = &num;
   arg->data.list.strings = &strings;
   num = 0;
}

CmdLine::Parameter::Parameter()
{
   arg = 0;
   /*
       arg = new Arg;
       arg->type = Arg::TERMINATE;
   */
}

CmdLine::Parameter::~Parameter()
{
   delete arg;
}
} // namespace soplex

//-----------------------------------------------------------------------------
//Emacs Local Variables:
//Emacs c-basic-offset:3
//Emacs tab-width:8
//Emacs indent-tabs-mode:nil
//Emacs End:
//-----------------------------------------------------------------------------
