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
#pragma ident "@(#) $Id: cmdline.h,v 1.2 2001/11/06 23:31:00 bzfkocht Exp $"


#ifndef _CMDLINE_H_
#define _CMDLINE_H_

//@ ----------------------------------------------------------------------------
/*  \Section{Imports}
    Import required system include files ...
 */
#include <assert.h>


/*  ... and class header files
 */

#include "dataarray.h"

namespace soplex
{

//@ ----------------------------------------------------------------------------
/* \Section{Class Declaration}
 */

/** command line argument parser.
    This class provides a command line argument parser for arguments list given
    in the #argc#, #argv# manner. It provides optional and required parameters,
    for the following parameter types: #int#, #double#, #string#, list of
    #strings# and simple boolean switches. Further, it provides a full help
    system for command line options.
 
    Parsing a command line is a matter of calling on of the static methods
    #parse()# with the command line, a help message and a full descriptor of all
    possible parameters.
 */
class CmdLine
{
public:
   /** command line argument desriptor.
       Instances of class #Parameter# are used to describe (optional or
       required) parameters to be parsed. The type of parameter is
       specified at construction time. Generally, the constructor of
       #Parameter# provides
       \begin{itemize}
       \item option character      that will be used for defining shortcut
                                   parameters (such as #-h#). This may be
                                   ommitted by passing 0.
       \item option strings        that will be used for verbose parameter
                                   selection (such as #--help#).
       \item help strings          will be output by #--help# in order to
                                   describe the parameter
       \item reference(s)          to the variables that will contain the
                                   parsed values.
       \item optional switch       to indicate optional parameters.
       \end{itemize}
       The options #-h# and #--help# are allways provided. No two
       parameters may use the same character. The type of the parameter is
       determined from the type of the passed variable. Lists will be
       returned using a pair #int len#, #char** strings# where #len# returns
       the number of strings found in the list, and #strings# points to the
       first argument of this list in #argv#.
    */
   class Parameter
   {
      friend class CmdLine;
      struct Arg
      {
         enum Type {
            SWITCH,
            INT,
            DOUBLE,
            STRING,
            OPT_LIST,
            LIST,
            TERMINATE
         } type;
         union Data {
            char* onoff;
            int* integer;
            double* real;
            const char** string;
            struct
            {
               int* num;
               const char*** strings;
            }
            list;
         } data;
         char copt;
         const char* sopt;
         const char* hlp;
         int optional;
         int done;
      };
      Arg* arg;

public:
      /// construct switch argument.
      Parameter(char opt_char,
                 const char opt_string[],
                 const char help[],
                 char& onoff,
                 int optional = 0
              );

      /// construct integer valued option argument.
      Parameter(char opt_char,
                 const char* opt_string,
                 const char* help,
                 int& val,
                 int optional = 0
              );

      /// construct double valued option argument.
      Parameter(char opt_char,
                 const char* opt_string,
                 const char* help,
                 double& val,
                 int optional = 0
              );

      /// construct string valued option argument.
      Parameter(char opt_char,
                 const char* opt_string,
                 const char* help,
                 const char*& string,
                 int optional = 0
              );

      /// construct string list valued option argument.
      Parameter(char opt_char,
                 const char* opt_string,
                 const char* help,
                 int& num,
                 const char**& strings,
                 int optional = 0
              );

      /// construct string list valued argument.
      Parameter(const char* help,
                 int& num,
                 const char**& strings,
                 int optional = 0
              );

      /// terminating argument.
      Parameter();

      ~Parameter();
      operator void*()
      {
         return arg;
      }
   };

   /** parses arguments in #argv#.
       Function #parse()# is called with the argument list to be passed in
       #argc# and #argv#, a help text in #help# and a list of any number of
       \Ref{Parameter} that are to be understood by the command line parser.
       This list must be terminated by an #Parameter# constructructed with
       its default constructor.  On success, the parser return 0, otherwise
       the number of passed argument that failed to match the possible
       arguments. After successfull termination, the values pointed to by
       the #Parameter#list are setup according to what has been passed in
       #argv#.
    */
   static int parse(int argc, char** argv, const char help[], ...);

   /** parses arguments in #argv#.
    *  Parses arguments in #argv# with 0-terminated list of
    *  parameters in #args#.
    */
   static int parse(int argc, char** argv, const char help[], Parameter* args);

   ///
   int isConsistent() const
   {
      return 1;
   }

private:
   typedef Parameter::Arg A;
   static int doOption(
      int argc,
      char** argv,
      const char help[],
      DataArray < A* > & alist,
      int& n,
      CmdLine::A* alisti,
      const char* start
  );
   static void usageLine(A* arg);
   static void usage(
      const char* prog,
      const char* help,
      DataArray < A* > & alist
  );
   static int parse(int argc, char** argv, const char help[], DataArray < A* > & args);
};

} // namespace soplex
#endif // _CMDLINE_H_

//-----------------------------------------------------------------------------
//Emacs Local Variables:
//Emacs mode:c++
//Emacs c-basic-offset:3
//Emacs tab-width:8
//Emacs indent-tabs-mode:nil
//Emacs End:
//-----------------------------------------------------------------------------
