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

/**@file  spxlpbase_real.cpp
 * @brief Saving LPs with Real values in a form suitable for SoPlex.
 */

#include <assert.h>
#include <stdio.h>
#include <ctype.h>
#include <iostream>

#include "spxdefines.h"
#include "spxlpbase.h"
#include "spxout.h"
#include "mpsinput.h"
#include "exceptions.h"

namespace soplex
{
// ---------------------------------------------------------------------------------------------------------------------
//  Specialization for reading LP format
// ---------------------------------------------------------------------------------------------------------------------

#define LPF_MAX_LINE_LEN  8192     ///< maximum length of a line (8190 + \\n + \\0)

/// Is \p c a \c space, \c tab, \c nl or \c cr ?
static inline bool LPFisSpace(int c)
{
   return (c == ' ') || (c == '\t') || (c == '\n') || (c == '\r');
}



/// Is there a number at the beginning of \p s ?
static bool LPFisValue(const char* s)
{
   return ((*s >= '0') && (*s <= '9')) || (*s == '+') || (*s == '-') || (*s == '.');
}



/// Is there a possible column name at the beginning of \p s ?
static bool LPFisColName(const char* s)
{
   // strchr() gives a true for the null char.
   if( *s == '\0' )
      return false;

   return ((*s >= 'A') && (*s <= 'Z'))
      || ((*s >= 'a') && (*s <= 'z'))
      || (strchr("!\"#$%&()/,;?@_'`{}|~", *s) != 0);
}



/// Is there a comparison operator at the beginning of \p s ?
static bool LPFisSense(const char* s)
{
   return (*s == '<') || (*s == '>') || (*s == '=');
}



static bool LPFisInfinity(const char* s)
{
   return ((s[0] == '-') || (s[0] == '+'))
      && (tolower(s[1]) == 'i')
      && (tolower(s[2]) == 'n')
      && (tolower(s[3]) == 'f');
}



static bool LPFisFree(const char* s)
{
   return (tolower(s[0]) == 'f')
      && ( tolower(s[1]) == 'r')
      && ( tolower(s[2]) == 'e')
      && ( tolower(s[3]) == 'e');
}



/// Read the next number and advance \p pos.
/** If only a sign is encountered, the number is assumed to be \c sign * 1.0.  This routine will not catch malformatted
 *  numbers like .e10 !
 */
static Real LPFreadValue(char*& pos, SPxOut* spxout)
{
   assert(LPFisValue(pos));

   char        tmp[LPF_MAX_LINE_LEN];
   const char* s = pos;
   char*       t;
   Real        value = 1.0;
   bool        has_digits = false;
   bool        has_emptyexponent = false;

   // 1. sign
   if( (*s == '+') || (*s == '-') )
      s++;

   // 2. Digits before the decimal dot
   while( (*s >= '0') && (*s <= '9') )
   {
      has_digits = true;
      s++;
   }

   // 3. Decimal dot
   if( *s == '.' )
   {
      s++;

      // 4. If there was a dot, possible digit behind it
      while( (*s >= '0') && (*s <= '9') )
      {
         has_digits = true;
         s++;
      }
   }

   // 5. Exponent
   if( tolower(*s) == 'e' )
   {
      has_emptyexponent = true;
      s++;

      // 6. Exponent sign
      if( (*s == '+') || (*s == '-') )
         s++;

      // 7. Exponent digits
      while( (*s >= '0') && (*s <= '9') )
      {
         has_emptyexponent = false;
         s++;
      }
   }
   assert(s != pos);

   if( has_emptyexponent )
   {
      MSG_WARNING( (*spxout), (*spxout) << "WLPFRD01 Warning: found empty exponent in LP file - check for forbidden variable names with initial 'e' or 'E'\n"; )
   }

   if( !has_digits )
      value = (*pos == '-') ? -1.0 : 1.0;
   else
   {
      for( t = tmp; pos != s; pos++ )
         *t++ = *pos;
      *t = '\0';
      value = atof(tmp);
   }

   pos += s - pos;

   assert(pos == s);

   MSG_DEBUG( std::cout << "DLPFRD01 LPFreadValue = " << value << std::endl; )

   if( LPFisSpace(*pos) )
      pos++;

   return value;
}



/// Read the next column name from the input.
/** The name read is looked up and if not found \p emptycol
 *  is added to \p colset. \p pos is advanced behind the name.
 *  @return The Index of the named column.
 */
static int LPFreadColName(char*& pos, NameSet* colnames, LPColSetBase<Real>& colset, const LPColBase<Real>* emptycol, SPxOut* spxout)
{
   assert(LPFisColName(pos));
   assert(colnames != 0);

   char        name[LPF_MAX_LINE_LEN];
   const char* s = pos;
   int         i;
   int         colidx;

   // These are the characters that are not allowed in a column name.
   while( (strchr("+-.<>= ", *s) == 0) && (*s != '\0') )
      s++;

   for( i = 0; pos != s; i++, pos++ )
      name[i] = *pos;

   name[i] = '\0';

   if( (colidx = colnames->number(name)) < 0 )
   {
      // We only add the name if we got an empty column.
      if( emptycol == 0 )
         MSG_WARNING( (*spxout), (*spxout) << "WLPFRD02 Unknown variable \"" << name << "\" "; )
      else
      {
         colidx = colnames->num();
         colnames->add(name);
         colset.add(*emptycol);
      }
   }

   MSG_DEBUG( std::cout << "DLPFRD03 LPFreadColName [" << name << "] = " << colidx << std::endl; )

   if( LPFisSpace(*pos) )
      pos++;

   return colidx;
}



/// Read the next <,>,=,==,<=,=<,>=,=> and advance \p pos.
static int LPFreadSense(char*& pos)
{
   assert(LPFisSense(pos));

   int sense = *pos++;

   if( (*pos == '<') || (*pos == '>') )
      sense = *pos++;
   else if( *pos == '=' )
      pos++;

   MSG_DEBUG( std::cout << "DLPFRD04 LPFreadSense = " << static_cast<char>(sense) << std::endl; )

   if( LPFisSpace(*pos) )
      pos++;

   return sense;
}



/// Is the \p keyword present in \p buf ? If yes, advance \p pos.
/** \p keyword should be lower case. It can contain optional sections which are enclosed in '[' ']' like "min[imize]".
 */
static bool LPFhasKeyword(char*& pos, const char* keyword)
{
   int i;
   int k;

   assert(keyword != 0);

   for( i = 0, k = 0; keyword[i] != '\0'; i++, k++ )
   {
      if( keyword[i] == '[' )
      {
         i++;

         // Here we assumed that we have a ']' for the '['.
         while( (tolower(pos[k]) == keyword[i]) && (pos[k] != '\0') )
         {
           k++;
           i++;
         }
         while( keyword[i] != ']' )
            i++;
         --k;
      }
      else
      {
         if( keyword[i] != tolower(pos[k]) )
            break;
      }
   }

   // we have to be at the end of the keyword and the word found on the line also has to end here.  Attention: The
   // LPFisSense is a kludge to allow LPFhasKeyword also to process Inf[inity] keywords in the bounds section.
   if( keyword[i] == '\0' && (pos[k] == '\0' || LPFisSpace(pos[k]) || LPFisSense(&pos[k])) )
   {
      pos += k;

      MSG_DEBUG( std::cout << "DLPFRD05 LPFhasKeyword: " << keyword << std::endl; )

      return true;
   }

   return false;
}



/// If \p buf start with "name:" store the name in \p rownames and advance \p pos.
static bool LPFhasRowName(char*& pos, NameSet* rownames)
{
   const char* s = strchr(pos, ':');

   if( s == 0 )
      return false;

   int dcolpos = int(s - pos);

   int end;
   int srt;

   // skip spaces between name and ":"
   for( end = dcolpos-1; end >= 0; end-- )
      if( pos[end] != ' ')
         break;

   // are there only spaces in front of the ":" ?
   if( end < 0 )
   {
      pos = &(pos[dcolpos+1]);
      return false;
   }

   // skip spaces in front of name
   for( srt = end-1; srt >= 0; srt-- )
      if( pos[srt] == ' ' )
         break;

   // go back to the non-space character
   srt++;

   assert(srt <= end && pos[srt] != ' ');

   char name[LPF_MAX_LINE_LEN];
   int i;
   int k = 0;

   for( i = srt; i <= end; i++ )
      name[k++] = pos[i];

   name[k] = '\0';

   if( rownames != 0 )
      rownames->add(name);

   pos = &(pos[dcolpos+1]);

   return true;
}



static Real LPFreadInfinity(char*& pos)
{
   assert(LPFisInfinity(pos));

   Real sense = (*pos == '-') ? -1.0 : 1.0;

   (void) LPFhasKeyword(++pos, "inf[inity]");

   return sense * infinity;
}



/// Read LP in "CPLEX LP File Format".
   /** The specification is taken from the ILOG CPLEX 7.0 Reference Manual, Appendix E, Page 527.
    *
    *  This routine should read (most?) valid LP format files.  What it will not do, is find all cases where a file is ill
    *  formed.  If this happens it may complain and read nothing or read "something".
    *
    *  Problem: A line ending in '+' or '-' followed by a line starting with a number, will be regarded as an error.
    *
    *  The reader will accept the keyword INT[egers] as a synonym for GEN[erals] which is an undocumented feature in CPLEX.
    *
    *  A difference to the CPLEX reader, is that no name for the objective row is required.
    *
    * The manual says the maximum allowed line length is 255 characters, but CPLEX does not complain if the lines are
    * longer.
    *
    *  @return true if the file was read correctly
    */
template <>
bool SPxLPBase<Real>::readLPF(
   std::istream& p_input,                ///< input stream.
   NameSet*      p_rnames,               ///< row names.
   NameSet*      p_cnames,               ///< column names.
   DIdxSet*      p_intvars)              ///< integer variables.
{
   enum
   {
      START, OBJECTIVE, CONSTRAINTS, BOUNDS, INTEGERS, BINARIES
   } section = START;

   NameSet* rnames;                      ///< row names.
   NameSet* cnames;                      ///< column names.

   LPColSetBase<Real> cset;              ///< the set of columns read.
   LPRowSetBase<Real> rset;              ///< the set of rows read.
   LPColBase<Real> emptycol;             ///< reusable empty column.
   LPRowBase<Real> row;                  ///< last assembled row.
   DSVectorBase<Real> vec;               ///< last assembled vector (from row).

   Real val = 1.0;
   int colidx;
   int sense = 0;

   char buf[LPF_MAX_LINE_LEN];
   char tmp[LPF_MAX_LINE_LEN];
   char line[LPF_MAX_LINE_LEN];
   int lineno = 0;
   bool unnamed = true;
   bool finished = false;
   bool other;
   bool have_value = true;
   int i;
   int k;
   char* s;
   char* pos;
   char* pos_old = 0;

   if( p_cnames )
      cnames = p_cnames;
   else
   {
      cnames = 0;
      spx_alloc(cnames);
      cnames = new (cnames) NameSet();
   }

   cnames->clear();

   if( p_rnames )
      rnames = p_rnames;
   else
   {
      try
      {
         rnames = 0;
         spx_alloc(rnames);
         rnames = new (rnames) NameSet();
      }
      catch( const SPxMemoryException& x )
      {
         if( !p_cnames )
         {
            cnames->~NameSet();
            spx_free(cnames);
         }
         throw x;
      }
   }

   rnames->clear();

   SPxLPBase<Real>::clear(); // clear the LP.

   //--------------------------------------------------------------------------
   //--- Main Loop
   //--------------------------------------------------------------------------
   for(;;)
   {
      // 0. Read a line from the file.
      if( !p_input.getline(buf, sizeof(buf)) )
      {
         if( strlen(buf) == LPF_MAX_LINE_LEN - 1 )
         {
            MSG_ERROR( std::cerr << "ELPFRD06 Line exceeds " << LPF_MAX_LINE_LEN - 2
                            << " characters" << std::endl; )
         }
         else
         {
            MSG_ERROR( std::cerr << "ELPFRD07 No 'End' marker found" << std::endl; )
            finished = true;
         }
         break;
      }
      lineno++;
      i   = 0;
      pos = buf;

      MSG_DEBUG( std::cout << "DLPFRD08 Reading line " << lineno
                        << " (pos=" << pos << ")" << std::endl; )

      // 1. Remove comments.
      if( 0 != (s = strchr(buf, '\\')) )
         *s = '\0';

      // 2. Look for keywords.
      if( section == START )
      {
         if( LPFhasKeyword(pos, "max[imize]") )
         {
            changeSense(SPxLPBase<Real>::MAXIMIZE);
            section = OBJECTIVE;
         }
         else if( LPFhasKeyword(pos, "min[imize]") )
         {
            changeSense(SPxLPBase<Real>::MINIMIZE);
            section = OBJECTIVE;
         }
      }
      else if( section == OBJECTIVE )
      {
         if( LPFhasKeyword(pos, "s[ubject][   ]t[o]")
            || LPFhasKeyword(pos, "s[uch][    ]t[hat]")
            || LPFhasKeyword(pos, "s[.][    ]t[.]")
            || LPFhasKeyword(pos, "lazy con[straints]") )
         {
            // store objective vector
            for( int j = vec.size() - 1; j >= 0; --j )
               cset.maxObj_w(vec.index(j)) = vec.value(j);
            // multiplication with -1 for minimization is done below
            vec.clear();
            have_value = true;
            val = 1.0;
            section = CONSTRAINTS;
         }
      }
      else if( section == CONSTRAINTS &&
              (LPFhasKeyword(pos, "s[ubject][   ]t[o]")
              || LPFhasKeyword(pos, "s[uch][    ]t[hat]")
              || LPFhasKeyword(pos, "s[.][    ]t[.]")) )
      {
         have_value = true;
         val = 1.0;
      }
      else
      {
         if( LPFhasKeyword(pos, "lazy con[straints]") )
            ;
         else if( LPFhasKeyword(pos, "bound[s]") )
            section = BOUNDS;
         else if( LPFhasKeyword(pos, "bin[ary]") )
            section = BINARIES;
         else if( LPFhasKeyword(pos, "bin[aries]") )
            section = BINARIES;
         else if( LPFhasKeyword(pos, "gen[erals]") )
            section = INTEGERS;
         else if( LPFhasKeyword(pos, "int[egers]") ) // this is undocumented
            section = INTEGERS;
         else if( LPFhasKeyword(pos, "end") )
         {
            finished = true;
            break;
         }
         else if( LPFhasKeyword(pos, "s[ubject][   ]t[o]") // second time
            || LPFhasKeyword(pos, "s[uch][    ]t[hat]")
            || LPFhasKeyword(pos, "s[.][    ]t[.]")
            || LPFhasKeyword(pos, "lazy con[straints]") )
         {
            // In principle this has to checked for all keywords above,
            // otherwise we just ignore any half finished constraint
            if( have_value )
               goto syntax_error;

            have_value = true;
            val = 1.0;
         }
      }

      // 3a. Look for row names in objective and drop it.
      if( section == OBJECTIVE )
         LPFhasRowName(pos, 0);

      // 3b. Look for row name in constraint and store it.
      if( section == CONSTRAINTS )
         if( LPFhasRowName(pos, rnames) )
            unnamed = false;

      // 4a. Remove initial spaces.
      while( LPFisSpace(pos[i]) )
         i++;

      // 4b. remove spaces if they do not appear before the name of a vaiable.
      for( k = 0; pos[i] != '\0'; i++ )
         if( !LPFisSpace(pos[i]) || LPFisColName(&pos[i + 1]) )
            tmp[k++] = pos[i];

      tmp[k] = '\0';

      // 5. Is this an empty line ?
      if( tmp[0] == '\0' )
         continue;

      // 6. Collapse sequences of '+' and '-'. e.g ++---+ => -
      for( i = 0, k = 0; tmp[i] != '\0'; i++ )
      {
         while( ((tmp[i] == '+') || (tmp[i] == '-')) && ((tmp[i + 1] == '+') || (tmp[i + 1] == '-')) )
         {
            if( tmp[i++] == '-' )
               tmp[i] = (tmp[i] == '-') ? '+' : '-';
         }
         line[k++] = tmp[i];
      }
      line[k] = '\0';

      //-----------------------------------------------------------------------
      //--- Line processing loop
      //-----------------------------------------------------------------------
      pos = line;

      MSG_DEBUG( std::cout << "DLPFRD09 pos=" << pos << std::endl; )

      // 7. We have something left to process.
      while( (pos != 0) && (*pos != '\0') )
      {
         // remember our position, so we are sure we make progress.
         pos_old = pos;

         // now process the sections
         switch( section )
         {
         case OBJECTIVE:
            if( LPFisValue(pos) )
            {
               Real pre_sign = 1.0;

               /* Already having here a value could only result from being the first number in a constraint, or a sign
                * '+' or '-' as last token on the previous line.
                */
               if( have_value )
               {
                  if( NE(spxAbs(val), 1.0) )
                     goto syntax_error;

                  if( EQ(val, -1.0) )
                     pre_sign = val;
               }
               have_value = true;
               val = LPFreadValue(pos, spxout) * pre_sign;
            }
            if( *pos == '\0' )
               continue;

            if( !have_value || !LPFisColName(pos) )
               goto syntax_error;

            have_value = false;
            colidx = LPFreadColName(pos, cnames, cset, &emptycol, spxout);
            vec.add(colidx, val);
            break;
         case CONSTRAINTS:
            if( LPFisValue(pos) )
            {
               Real pre_sign = 1.0;

               /* Already having here a value could only result from being the first number in a constraint, or a sign
                * '+' or '-' as last token on the previous line.
                */
               if( have_value )
               {
                  if( NE(spxAbs(val), 1.0) )
                     goto syntax_error;

                  if( EQ(val, -1.0) )
                     pre_sign = val;
               }

               have_value = true;
               val = LPFreadValue(pos, spxout) * pre_sign;

               if( sense != 0 )
               {
                  if( sense == '<' )
                  {
                     row.setLhs(-infinity);
                     row.setRhs(val);
                  }
                  else if( sense == '>' )
                  {
                     row.setLhs(val);
                     row.setRhs(infinity);
                  }
                  else
                  {
                     assert(sense == '=');

                     row.setLhs(val);
                     row.setRhs(val);
                  }
                  row.setRowVector(vec);
                  rset.add(row);
                  vec.clear();

                  if( !unnamed )
                     unnamed = true;
                  else
                  {
                     char name[16];
                     sprintf(name, "C%d", rset.num());
                     rnames->add(name);
                  }
                  have_value = true;
                  val = 1.0;
                  sense = 0;
                  pos = 0;
                  // next line
                  continue;
               }
            }
            if( *pos == '\0' )
               continue;

            if( have_value )
            {
               if( LPFisColName(pos) )
               {
                  colidx = LPFreadColName(pos, cnames, cset, &emptycol, spxout);

                  if( val != 0.0 )
                  {
                     // Do we have this index already in the row?
                     int n = vec.number(colidx);

                     // if not, add it
                     if( n < 0 )
                        vec.add(colidx, val);
                     // if yes, add them up and remove the element if it amounts to zero
                     else
                     {
                        assert(vec.index(n) == colidx);

                        val += vec.value(n);

                        if( val == 0.0 )
                           vec.remove(n);
                        else
                           vec.value(n) = val;

                        assert(cnames->has(colidx));

                        MSG_WARNING( (*spxout), (*spxout) << "WLPFRD10 Duplicate index "
                                            << (*cnames)[colidx]
                                            << " in line " << lineno
                                            << std::endl; )
                     }
                  }
                  have_value = false;
               }
               else
               {
                  // We have a row like c1: <= 5 with no variables. We can not handle 10 <= 5; issue a syntax error.
                  if( val != 1.0 )
                     goto syntax_error;

                  // If the next thing is not the sense we give up also.
                  if( !LPFisSense(pos) )
                     goto syntax_error;

                  have_value = false;
               }
            }
            assert(!have_value);

            if( LPFisSense(pos) )
               sense = LPFreadSense(pos);
            break;
         case BOUNDS:
            other = false;
            sense = 0;

            if( LPFisValue(pos) )
            {
               val = LPFisInfinity(pos) ? LPFreadInfinity(pos) : LPFreadValue(pos, spxout);

               if( !LPFisSense(pos) )
                  goto syntax_error;

               sense = LPFreadSense(pos);
               other = true;
            }

            if( !LPFisColName(pos) )
               goto syntax_error;

            if( (colidx = LPFreadColName(pos, cnames, cset, 0, spxout)) < 0 )
            {
               MSG_WARNING( (*spxout), (*spxout) << "WLPFRD11 in Bounds section line "
                                   << lineno << " ignored" << std::endl; )
               *pos = '\0';
               continue;
            }

            if( sense )
            {
               if( sense == '<' )
                  cset.lower_w(colidx) = val;
               else if( sense == '>' )
                  cset.upper_w(colidx) = val;
               else
               {
                  assert(sense == '=');
                  cset.lower_w(colidx) = val;
                  cset.upper_w(colidx) = val;
               }
            }

            if( LPFisFree(pos) )
            {
               cset.lower_w(colidx) = -infinity;
               cset.upper_w(colidx) =  infinity;
               other = true;
               pos += 4;  // set position after the word "free"
            }
            else if( LPFisSense(pos) )
            {
               sense = LPFreadSense(pos);
               other = true;

               if( !LPFisValue(pos) )
                  goto syntax_error;

               val = LPFisInfinity(pos) ? LPFreadInfinity(pos) : LPFreadValue(pos, spxout);

               if( sense == '<' )
                  cset.upper_w(colidx) = val;
               else if( sense == '>' )
                  cset.lower_w(colidx) = val;
               else
               {
                  assert(sense == '=');
                  cset.lower_w(colidx) = val;
                  cset.upper_w(colidx) = val;
               }
            }

            /* Do we have only a single column name in the input line?  We could ignore this savely, but it is probably
             * a sign of some other error.
             */
            if( !other )
               goto syntax_error;
            break;
         case BINARIES:
         case INTEGERS:
            if( (colidx = LPFreadColName(pos, cnames, cset, 0, spxout)) < 0 )
            {
               MSG_WARNING( (*spxout), (*spxout) << "WLPFRD12 in Binary/General section line " << lineno << " ignored" << std::endl; )
            }
            else
            {
               if( section == BINARIES )
               {
		  if( cset.lower(colidx) < 0.0 )
		  {
		     cset.lower_w(colidx) = 0.0;
		  }
		  if( cset.upper(colidx) > 1.0 )
		  {
		     cset.upper_w(colidx) = 1.0;
		  }
               }

               if( p_intvars != 0 )
                  p_intvars->addIdx(colidx);
            }
            break;
         case START:
            MSG_ERROR( std::cerr << "ELPFRD13 This seems to be no LP format file" << std::endl; )
            goto syntax_error;
         default:
            throw SPxInternalCodeException("XLPFRD01 This should never happen.");
         }

         if( pos == pos_old )
            goto syntax_error;
      }
   }

   assert(isConsistent());

   addCols(cset);
   assert(isConsistent());

   addRows(rset);
   assert(isConsistent());

syntax_error:
   if( finished )
   {
      MSG_INFO2( (*spxout), (*spxout) << "Finished reading " << lineno << " lines" << std::endl; )
   }
   else
      MSG_ERROR( std::cerr << "ELPFRD15 Syntax error in line " << lineno << std::endl; )

   if( p_cnames == 0 )
      spx_free(cnames);
   if( p_rnames == 0 )
      spx_free(rnames);

   return finished;
}



// ---------------------------------------------------------------------------------------------------------------------
// Specialization for reading MPS format
// ---------------------------------------------------------------------------------------------------------------------

/// Process NAME section.
static void MPSreadName(MPSInput& mps, SPxOut* spxout)
{
   do
   {
      // This has to be the Line with the NAME section.
      if( !mps.readLine() || (mps.field0() == 0) || strcmp(mps.field0(), "NAME") )
         break;

      // Sometimes the name is omitted.
      mps.setProbName((mps.field1() == 0) ? "_MPS_" : mps.field1());

      MSG_INFO2( (*spxout), (*spxout) << "IMPSRD01 Problem name   : " << mps.probName() << std::endl; )

      // This has to be a new section
      if( !mps.readLine() || (mps.field0() == 0) )
         break;

      if( !strcmp(mps.field0(), "ROWS") )
         mps.setSection(MPSInput::ROWS);
      else if( !strncmp(mps.field0(), "OBJSEN", 6) )
         mps.setSection(MPSInput::OBJSEN);
      else if( !strcmp(mps.field0(), "OBJNAME") )
         mps.setSection(MPSInput::OBJNAME);
      else
         break;

      return;
   }
   while(false);

   mps.syntaxError();
}



/// Process OBJSEN section. This Section is an ILOG extension.
static void MPSreadObjsen(MPSInput& mps)
{
   do
   {
      // This has to be the Line with MIN or MAX.
      if( !mps.readLine() || (mps.field1() == 0) )
         break;

      if( !strcmp(mps.field1(), "MIN") )
         mps.setObjSense(MPSInput::MINIMIZE);
      else if( !strcmp(mps.field1(), "MAX") )
         mps.setObjSense(MPSInput::MAXIMIZE);
      else
         break;

      // Look for ROWS or OBJNAME Section
      if( !mps.readLine() || (mps.field0() == 0) )
         break;

      if( !strcmp(mps.field0(), "ROWS") )
         mps.setSection(MPSInput::ROWS);
      else if( !strcmp(mps.field0(), "OBJNAME") )
         mps.setSection(MPSInput::OBJNAME);
      else
         break;

      return;
   }
   while(false);

   mps.syntaxError();
}



/// Process OBJNAME section. This Section is an ILOG extension.
static void MPSreadObjname(MPSInput& mps)
{
   do
   {
      // This has to be the Line with the name.
      if( !mps.readLine() || (mps.field1() == 0) )
         break;

      mps.setObjName(mps.field1());

      // Look for ROWS Section
      if( !mps.readLine() || (mps.field0() == 0) )
         break;

      if( strcmp(mps.field0(), "ROWS") )
         break;

      mps.setSection(MPSInput::ROWS);

      return;
   }
   while(false);

   mps.syntaxError();
}



/// Process ROWS section.
static void MPSreadRows(MPSInput& mps, LPRowSetBase<Real>& rset, NameSet& rnames, SPxOut* spxout)
{
   LPRowBase<Real> row;

   while( mps.readLine() )
   {
      if( mps.field0() != 0 )
      {
         MSG_INFO2( (*spxout), (*spxout) << "IMPSRD02 Objective name : " << mps.objName() << std::endl; )

         if( strcmp(mps.field0(), "COLUMNS") )
            break;

         mps.setSection(MPSInput::COLUMNS);

         return;
      }

      if( (mps.field1() == 0) || (mps.field2() == 0) )
         break;

      if( *mps.field1() == 'N' )
      {
         if( *mps.objName() == '\0' )
            mps.setObjName(mps.field2());
      }
      else
      {
         if( rnames.has(mps.field2()) )
            break;

         rnames.add(mps.field2());

         switch( *mps.field1() )
         {
         case 'G':
            row.setLhs(0.0);
            row.setRhs(infinity);
            break;
         case 'E':
            row.setLhs(0.0);
            row.setRhs(0.0);
            break;
         case 'L':
            row.setLhs(-infinity);
            row.setRhs(0.0);
            break;
         default:
            mps.syntaxError();
            return;
         }

         rset.add(row);
      }

      assert((*mps.field1() == 'N') || (rnames.number(mps.field2()) == rset.num() - 1));
   }

   mps.syntaxError();
}



/// Process COLUMNS section.
static void MPSreadCols(MPSInput& mps, const LPRowSetBase<Real>& rset, const NameSet&  rnames, LPColSetBase<Real>& cset, NameSet& cnames, DIdxSet* intvars)
{
   Real val;
   int idx;
   char colname[MPSInput::MAX_LINE_LEN] = { '\0' };
   LPColBase<Real> col(rset.num());
   DSVectorBase<Real> vec;

   col.setObj(0.0);
   vec.clear();

   while( mps.readLine() )
   {
      if( mps.field0() != 0 )
      {
         if( strcmp(mps.field0(), "RHS") )
            break;

         if( colname[0] != '\0' )
         {
            col.setColVector(vec);
            cset.add(col);
         }

         mps.setSection(MPSInput::RHS);

         return;
      }

      if( (mps.field1() == 0) || (mps.field2() == 0) || (mps.field3() == 0) )
         break;

      // new column?
      if( strcmp(colname, mps.field1()) )
      {
         // first column?
         if( colname[0] != '\0' )
         {
            col.setColVector(vec);
            cset.add(col);
         }

         // save copy of string (make sure string ends with \0)
         strncpy(colname, mps.field1(), MPSInput::MAX_LINE_LEN-1);
         colname[MPSInput::MAX_LINE_LEN-1] = '\0';
         cnames.add(colname);
         vec.clear();
         col.setObj(0.0);
         col.setLower(0.0);
         col.setUpper(infinity);

         if( mps.isInteger() )
         {
            assert(cnames.number(colname) == cset.num());

            if( intvars != 0 )
               intvars->addIdx(cnames.number(colname));

            // for Integer variable the default bounds are 0/1
            col.setUpper(1.0);
         }
      }

      val = atof(mps.field3());

      if( !strcmp(mps.field2(), mps.objName()) )
         col.setObj(val);
      else
      {
         if( (idx = rnames.number(mps.field2())) < 0 )
            mps.entryIgnored("Column", mps.field1(), "row", mps.field2());
         else if( val != 0.0 )
            vec.add(idx, val);
      }

      if( mps.field5() != 0 )
      {
         assert(mps.field4() != 0);

         val = atof(mps.field5());

         if( !strcmp(mps.field4(), mps.objName()) )
            col.setObj(val);
         else
         {
            if( (idx = rnames.number(mps.field4())) < 0 )
               mps.entryIgnored("Column", mps.field1(), "row", mps.field4());
            else if( val != 0.0 )
               vec.add(idx, val);
         }
      }
   }

   mps.syntaxError();
}



/// Process RHS section.
static void MPSreadRhs(MPSInput& mps, LPRowSetBase<Real>& rset, const NameSet& rnames, SPxOut* spxout)
{
   char rhsname[MPSInput::MAX_LINE_LEN] = { '\0' };
   char addname[MPSInput::MAX_LINE_LEN] = { '\0' };
   int idx;
   Real val;

   while( mps.readLine() )
   {
      if( mps.field0() != 0 )
      {
         MSG_INFO2( (*spxout), (*spxout) << "IMPSRD03 RHS name       : " << rhsname  << std::endl; );

         if( !strcmp(mps.field0(), "RANGES") )
            mps.setSection(MPSInput::RANGES);
         else if( !strcmp(mps.field0(), "BOUNDS") )
            mps.setSection(MPSInput::BOUNDS);
         else if( !strcmp(mps.field0(), "ENDATA") )
            mps.setSection(MPSInput::ENDATA);
         else
            break;

         return;
      }

      if( ((mps.field2() != 0) && (mps.field3() == 0)) || ((mps.field4() != 0) && (mps.field5() == 0)) )
         mps.insertName("_RHS_");

      if( (mps.field1() == 0) || (mps.field2() == 0) || (mps.field3() == 0) )
         break;

      if( *rhsname == '\0' )
         strcpy(rhsname, mps.field1());

      if( strcmp(rhsname, mps.field1()) )
      {
         if( strcmp(addname, mps.field1()) )
         {
            assert(strlen(mps.field1()) < MPSInput::MAX_LINE_LEN);
            strcpy(addname, mps.field1());
            MSG_INFO3( (*spxout), (*spxout) << "IMPSRD07 RHS ignored    : " << addname << std::endl );
         }
      }
      else
      {
         if( (idx = rnames.number(mps.field2())) < 0 )
            mps.entryIgnored("RHS", mps.field1(), "row", mps.field2());
         else
         {
            val = atof(mps.field3());

            // LE or EQ
            if( rset.rhs(idx) < infinity )
               rset.rhs_w(idx) = val;
            // GE or EQ
            if( rset.lhs(idx) > -infinity )
               rset.lhs_w(idx) = val;
         }

         if( mps.field5() != 0 )
         {
            if( (idx = rnames.number(mps.field4())) < 0 )
               mps.entryIgnored("RHS", mps.field1(), "row", mps.field4());
            else
            {
               val = atof(mps.field5());

               // LE or EQ
               if( rset.rhs(idx) < infinity )
                  rset.rhs_w(idx) = val;
               // GE or EQ
               if( rset.lhs(idx) > -infinity )
                  rset.lhs_w(idx) = val;
            }
         }
      }
   }

   mps.syntaxError();
}



/// Process RANGES section.
static void MPSreadRanges(MPSInput& mps,  LPRowSetBase<Real>& rset, const NameSet& rnames, SPxOut* spxout)
{
   char rngname[MPSInput::MAX_LINE_LEN] = { '\0' };
   int idx;
   Real val;

   while( mps.readLine() )
   {
      if( mps.field0() != 0 )
      {
         MSG_INFO2( (*spxout), (*spxout) << "IMPSRD04 Range name     : " << rngname << std::endl; );

         if( !strcmp(mps.field0(), "BOUNDS") )
            mps.setSection(MPSInput::BOUNDS);
         else if( !strcmp(mps.field0(), "ENDATA") )
            mps.setSection(MPSInput::ENDATA);
         else
            break;

         return;
      }

      if( ((mps.field2() != 0) && (mps.field3() == 0)) || ((mps.field4() != 0) && (mps.field5() == 0)) )
         mps.insertName("_RNG_");

      if( (mps.field1() == 0) || (mps.field2() == 0) || (mps.field3() == 0) )
         break;

      if( *rngname == '\0' )
      {
         assert(strlen(mps.field2()) < MPSInput::MAX_LINE_LEN);
         strcpy(rngname, mps.field1());
      }

      /* The rules are:
       * Row Sign   LHS             RHS
       * ----------------------------------------
       *  G   +/-   rhs             rhs + |range|
       *  L   +/-   rhs - |range|   rhs
       *  E   +     rhs             rhs + range
       *  E   -     rhs + range     rhs
       * ----------------------------------------
       */
      if( !strcmp(rngname, mps.field1()) )
      {
         if( (idx = rnames.number(mps.field2())) < 0 )
            mps.entryIgnored("Range", mps.field1(), "row", mps.field2());
         else
         {
            val = atof(mps.field3());

            // EQ
            if( (rset.lhs(idx) > -infinity) && (rset.rhs_w(idx) <  infinity) )
            {
               assert(rset.lhs(idx) == rset.rhs(idx));

               if( val >= 0 )
                  rset.rhs_w(idx) += val;
               else
                  rset.lhs_w(idx) += val;
            }
            else
            {
               // GE
               if( rset.lhs(idx) > -infinity )
                  rset.rhs_w(idx)  = rset.lhs(idx) + spxAbs(val);
               // LE
               else
                  rset.lhs_w(idx)  = rset.rhs(idx) - spxAbs(val);
            }
         }

         if( mps.field5() != 0 )
         {
            if( (idx = rnames.number(mps.field4())) < 0 )
               mps.entryIgnored("Range", mps.field1(), "row", mps.field4());
            else
            {
               val = atof(mps.field5());

               // EQ
               if( (rset.lhs(idx) > -infinity) && (rset.rhs(idx) <  infinity) )
               {
                  assert(rset.lhs(idx) == rset.rhs(idx));

                  if( val >= 0 )
                     rset.rhs_w(idx) += val;
                  else
                     rset.lhs_w(idx) += val;
               }
               else
               {
                  // GE
                  if( rset.lhs(idx) > -infinity )
                     rset.rhs_w(idx)  = rset.lhs(idx) + spxAbs(val);
                  // LE
                  else
                     rset.lhs_w(idx)  = rset.rhs(idx) - spxAbs(val);
               }
            }
         }
      }
   }

   mps.syntaxError();
}



/// Process BOUNDS section.
static void MPSreadBounds(MPSInput& mps, LPColSetBase<Real>& cset, const NameSet& cnames, DIdxSet* intvars, SPxOut* spxout)
{
   DIdxSet oldbinvars;
   char bndname[MPSInput::MAX_LINE_LEN] = { '\0' };
   int  idx;
   Real val;

   while( mps.readLine() )
   {
      if( mps.field0() != 0 )
      {
         MSG_INFO2( (*spxout), (*spxout) << "IMPSRD05 Bound name     : " << bndname << std::endl; )

         if( strcmp(mps.field0(), "ENDATA") )
            break;

         mps.setSection(MPSInput::ENDATA);

         return;
      }

      // Is the value field used ?
      if(  (!strcmp(mps.field1(), "LO"))
         || (!strcmp(mps.field1(), "UP"))
         || (!strcmp(mps.field1(), "FX"))
         || (!strcmp(mps.field1(), "LI"))
         || (!strcmp(mps.field1(), "UI")) )
      {
         if( (mps.field3() != 0) && (mps.field4() == 0) )
            mps.insertName("_BND_", true);
      }
      else
      {
         if( (mps.field2() != 0) && (mps.field3() == 0) )
            mps.insertName("_BND_", true);
      }

      if( (mps.field1() == 0) || (mps.field2() == 0) || (mps.field3() == 0) )
         break;

      if( *bndname == '\0' )
      {
         assert(strlen(mps.field2()) < MPSInput::MAX_LINE_LEN);
         strcpy(bndname, mps.field2());
      }

      // Only read the first Bound in section
      if( !strcmp(bndname, mps.field2()) )
      {
         if( (idx = cnames.number(mps.field3())) < 0 )
            mps.entryIgnored("column", mps.field3(), "bound", bndname);
         else
         {
            val = (mps.field4() == 0) ? 0.0 : atof(mps.field4());

            // ILOG extension (Integer Bound)
            if( mps.field1()[1] == 'I' )
            {
               if( intvars != 0 )
                  intvars->addIdx(idx);

               // if the variable has appeared in the MARKER section of the COLUMNS section then its default bounds were
               // set to 0,1; the first time it is declared integer we need to change to default bounds 0,infinity
               if( oldbinvars.number(idx) < 0 )
               {
                  cset.upper_w(idx) = infinity;
                  oldbinvars.addIdx(idx);
               }
            }

            switch( *mps.field1() )
            {
            case 'L':
               cset.lower_w(idx) = val;
               break;
            case 'U':
               cset.upper_w(idx) = val;
               break;
            case 'F':
               if( mps.field1()[1] == 'X' )
               {
                  cset.lower_w(idx) = val;
                  cset.upper_w(idx) = val;
               }
               else
               {
                  cset.lower_w(idx) = -infinity;
                  cset.upper_w(idx) = infinity;
               }
               break;
            case 'M':
               cset.lower_w(idx) = -infinity;
               break;
            case 'P':
               cset.upper_w(idx) = infinity;
               break;
               // Ilog extension (Binary)
            case 'B':
               cset.lower_w(idx) = 0.0;
               cset.upper_w(idx) = 1.0;

               if( intvars != 0 )
                  intvars->addIdx(idx);
               break;
            default:
               mps.syntaxError();
               return;
            }
         }
      }
   }

   mps.syntaxError();
}



/// Read LP in MPS File Format.
/**
 *  The specification is taken from the IBM Optimization Library Guide and Reference, online available at
 *  http://www.software.ibm.com/sos/features/libuser.htm and from the ILOG CPLEX 7.0 Reference Manual, Appendix E, Page
 *  531.
 *
 *  This routine should read all valid MPS format files.  What it will not do, is find all cases where a file is ill
 *  formed.  If this happens it may complain and read nothing or read "something".
 *
 *  @return true if the file was read correctly.
 */
#define INIT_COLS 10000 ///< initialy allocated columns.
#define INIT_NZOS 100000 ///< initialy allocated non zeros.
template <>
bool SPxLPBase<Real>::readMPS(
   std::istream& p_input,           ///< input stream.
   NameSet*      p_rnames,          ///< row names.
   NameSet*      p_cnames,          ///< column names.
   DIdxSet*      p_intvars)         ///< integer variables.
{
   LPRowSetBase<Real>& rset = *this;
   LPColSetBase<Real>& cset = *this;
   NameSet* rnames;
   NameSet* cnames;

   if( p_cnames )
      cnames = p_cnames;
   else
   {
      cnames = 0;
      spx_alloc(cnames);
      cnames = new (cnames) NameSet();
   }

   cnames->clear();

   if( p_rnames )
      rnames = p_rnames;
   else
   {
      try
      {
         rnames = 0;
         spx_alloc(rnames);
         rnames = new (rnames) NameSet();
      }
      catch( const SPxMemoryException& x)
      {
         if( !p_cnames )
         {
            cnames->~NameSet();
            spx_free(cnames);
         }
         throw x;
      }
   }

   rnames->clear();

   SPxLPBase<Real>::clear(); // clear the LP.

   cset.memRemax(INIT_NZOS);
   cset.reMax(INIT_COLS);

   MPSInput mps(p_input);

   MPSreadName(mps, spxout);

   if( mps.section() == MPSInput::OBJSEN )
      MPSreadObjsen(mps);

   if( mps.section() == MPSInput::OBJNAME )
      MPSreadObjname(mps);

   if( mps.section() == MPSInput::ROWS )
      MPSreadRows(mps, rset, *rnames, spxout);

   addedRows(rset.num());

   if( mps.section() == MPSInput::COLUMNS )
      MPSreadCols(mps, rset, *rnames, cset, *cnames, p_intvars);

   if( mps.section() == MPSInput::RHS )
      MPSreadRhs(mps, rset, *rnames, spxout);

   if( mps.section() == MPSInput::RANGES )
      MPSreadRanges(mps, rset, *rnames, spxout);

   if( mps.section() == MPSInput::BOUNDS )
      MPSreadBounds(mps, cset, *cnames, p_intvars, spxout);

   if( mps.section() != MPSInput::ENDATA )
      mps.syntaxError();

   if( mps.hasError() )
      clear();
   else
   {
      changeSense(mps.objSense() == MPSInput::MINIMIZE ? SPxLPBase<Real>::MINIMIZE : SPxLPBase<Real>::MAXIMIZE);

      MSG_INFO2( (*spxout), (*spxout) << "IMPSRD06 Objective sense: " << ((mps.objSense() == MPSInput::MINIMIZE) ? "Minimize\n" : "Maximize\n") );

      added2Set(
         *(reinterpret_cast<SVSetBase<Real>*>(static_cast<LPRowSetBase<Real>*>(this))),
         *(reinterpret_cast<SVSetBase<Real>*>(static_cast<LPColSetBase<Real>*>(this))),
         cset.num());
      addedCols(cset.num());

      assert(isConsistent());
   }

   if( p_cnames == 0 )
   {
      cnames->~NameSet();
      spx_free(cnames);
   }

   if( p_rnames == 0 )
   {
      rnames->~NameSet();
      spx_free(rnames);
   }

   return !mps.hasError();
}



// ---------------------------------------------------------------------------------------------------------------------
// Specialization for writing LP format
// ---------------------------------------------------------------------------------------------------------------------

// get the name of a row or construct one
static const char* LPFgetRowName(
   const SPxLPBase<Real>& p_lp,
   int                    p_idx,
   const NameSet*         p_rnames,
   char*                  p_buf,
   int                    p_num_written_rows
   )
{
   assert(p_buf != 0);
   assert(p_idx >= 0);
   assert(p_idx <  p_lp.nRows());

   if( p_rnames != 0 )
   {
      DataKey key = p_lp.rId(p_idx);

      if( p_rnames->has(key) )
         return (*p_rnames)[key];
   }

   sprintf(p_buf, "C%d", p_num_written_rows);

   return p_buf;
}



// get the name of a column or construct one
static const char* getColName(
   const SPxLPBase<Real>& p_lp,
   int                    p_idx,
   const NameSet*         p_cnames,
   char*                  p_buf
   )
{
   assert(p_buf != 0);
   assert(p_idx >= 0);
   assert(p_idx <  p_lp.nCols());

   if( p_cnames != 0 )
   {
      DataKey key = p_lp.cId(p_idx);

      if( p_cnames->has(key) )
         return (*p_cnames)[key];
   }

   sprintf(p_buf, "x%d", p_idx);

   return p_buf;
}



// write an SVector
#define NUM_ENTRIES_PER_LINE 5
static void LPFwriteSVector(
   const SPxLPBase<Real>&   p_lp,       ///< the LP
   std::ostream&            p_output,   ///< output stream
   const NameSet*           p_cnames,   ///< column names
   const SVectorBase<Real>& p_svec )    ///< vector to write
{

   char name[16];
   int num_coeffs = 0;

   for( int j = 0; j < p_lp.nCols(); ++j )
   {
      const Real coeff = p_svec[j];

      if( coeff == 0 )
         continue;

      if( num_coeffs == 0 )
         p_output << coeff << " " << getColName(p_lp, j, p_cnames, name);
      else
      {
         // insert a line break every NUM_ENTRIES_PER_LINE columns
         if( num_coeffs % NUM_ENTRIES_PER_LINE == 0 )
            p_output << "\n\t";

         if( coeff < 0 )
            p_output << " - " << -coeff;
         else
            p_output << " + " << coeff;

         p_output << " " << getColName(p_lp, j, p_cnames, name);
      }

      ++num_coeffs;
   }
}



// write the objective
static void LPFwriteObjective(
   const SPxLPBase<Real>& p_lp,       ///< the LP
   std::ostream&          p_output,   ///< output stream
   const NameSet*         p_cnames    ///< column names
   )
{

   const int sense = p_lp.spxSense();

   p_output << ((sense == SPxLPBase<Real>::MINIMIZE) ? "Minimize\n" : "Maximize\n");
   p_output << "  obj: ";

   const VectorBase<Real>& obj = p_lp.maxObj();
   DSVectorBase<Real> svec(obj.dim());
   svec.operator=(obj);
   svec *= Real(sense);
   LPFwriteSVector(p_lp, p_output, p_cnames, svec);
   p_output << "\n";
}



// write non-ranged rows
static void LPFwriteRow(
   const SPxLPBase<Real>&   p_lp,       ///< the LP
   std::ostream&            p_output,   ///< output stream
   const NameSet*           p_cnames,   ///< column names
   const SVectorBase<Real>& p_svec,     ///< vector of the row
   const Real&              p_lhs,      ///< lhs of the row
   const Real&              p_rhs       ///< rhs of the row
   )
{

   LPFwriteSVector(p_lp, p_output, p_cnames, p_svec);

   if( p_lhs == p_rhs )
      p_output << " = " << p_rhs;
   else if( p_lhs <= -infinity )
      p_output << " <= " << p_rhs;
   else
   {
      assert(p_rhs >= infinity);
      p_output << " >= " << p_lhs;
   }

   p_output << "\n";
}



// write all rows
static void LPFwriteRows(
   const SPxLPBase<Real>& p_lp,       ///< the LP
   std::ostream&          p_output,   ///< output stream
   const NameSet*         p_rnames,   ///< row names
   const NameSet*         p_cnames   ///< column names
   )
{

   char name[16];

   p_output << "Subject To\n";

   int num_written_rows = 0;  // num_written_rows > nRows with ranged rows
   for( int i = 0; i < p_lp.nRows(); ++i )
   {
      const Real lhs = p_lp.lhs(i);
      const Real rhs = p_lp.rhs(i);

      if( lhs > -infinity && rhs < infinity && lhs != rhs )
      {
         // ranged row -> write two non-ranged rows
         p_output << " " << LPFgetRowName(p_lp, i, p_rnames, name, ++num_written_rows) << "_1 : ";
         LPFwriteRow(p_lp, p_output, p_cnames, p_lp.rowVector(i), lhs, infinity);

         p_output << " " << LPFgetRowName(p_lp, i, p_rnames, name, ++num_written_rows) << "_2 : ";
         LPFwriteRow(p_lp, p_output, p_cnames, p_lp.rowVector(i), -infinity, rhs);
      }
      else
      {
         p_output << " " << LPFgetRowName(p_lp, i, p_rnames, name, ++num_written_rows) << " : ";
         LPFwriteRow(p_lp, p_output, p_cnames, p_lp.rowVector(i), lhs, rhs);
      }
   }
}



// write the variable bounds
// (the default bounds 0 <= x <= infinity are not written)
static void LPFwriteBounds(
   const SPxLPBase<Real>&   p_lp,       ///< the LP to write
   std::ostream&            p_output,   ///< output stream
   const NameSet*           p_cnames    ///< column names
   )
{

   char name[16];

   p_output << "Bounds\n";

   for( int j = 0; j < p_lp.nCols(); ++j )
   {
      const Real lower = p_lp.lower(j);
      const Real upper = p_lp.upper(j);

      if( lower == upper )
      {
         p_output << "  "   << getColName(p_lp, j, p_cnames, name) << " = "  << upper << '\n';
      }
      else if( lower > -infinity )
      {
         if( upper < infinity )
         {
            // range bound
            if( lower != 0 )
               p_output << "  "   << lower << " <= "
                        << getColName(p_lp, j, p_cnames, name)
                        << " <= " << upper << '\n';
            else
               p_output << "  "   << getColName(p_lp, j, p_cnames, name)
                        << " <= " << upper << '\n';
         }
         else if( lower != 0 )
            p_output << "  " << lower << " <= "
                     << getColName(p_lp, j, p_cnames, name)
                     << '\n';
      }
      else if( upper < infinity )
         p_output << "   -Inf <= "
                  << getColName(p_lp, j, p_cnames, name)
                  << " <= " << upper << '\n';
      else
         p_output << "  "   << getColName(p_lp, j, p_cnames, name)
                  << " free\n";
   }
}



// write the generals section
static void LPFwriteGenerals(
   const SPxLPBase<Real>&   p_lp,         ///< the LP to write
   std::ostream&            p_output,     ///< output stream
   const NameSet*           p_cnames,     ///< column names
   const DIdxSet*           p_intvars     ///< integer variables
   )
{

   char name[16];

   if( p_intvars == NULL || p_intvars->size() <= 0 )
      return;  // no integer variables

   p_output << "Generals\n";

   for( int j = 0; j < p_lp.nCols(); ++j )
      if( p_intvars->number(j) >= 0 )
         p_output << "  " << getColName(p_lp, j, p_cnames, name) << "\n";
}



/// Write LP in LP Format.
template <>
void SPxLPBase<Real>::writeLPF(
   std::ostream&  p_output,          ///< output stream
   const NameSet* p_rnames,          ///< row names
   const NameSet* p_cnames,          ///< column names
   const DIdxSet* p_intvars          ///< integer variables
   ) const
{

   p_output << std::setprecision(15);
   LPFwriteObjective(*this, p_output, p_cnames);
   LPFwriteRows(*this, p_output, p_rnames, p_cnames);
   LPFwriteBounds(*this, p_output, p_cnames);
   LPFwriteGenerals(*this, p_output, p_cnames, p_intvars);

   p_output << "End" << std::endl;
}



// ---------------------------------------------------------------------------------------------------------------------
// Specialization for writing MPS format
// ---------------------------------------------------------------------------------------------------------------------

static void MPSwriteRecord(
   std::ostream&  os,
   const char*    indicator,
   const char*    name,
   const char*    name1  = 0,
   const Real     value1 = 0.0,
   const char*    name2  = 0,
   const Real     value2 = 0.0
   )
{
   char buf[81];

   sprintf(buf, " %-2.2s %-8.8s", (indicator == 0) ? "" : indicator, (name == 0)      ? "" : name);
   os << buf;

   if( name1 != 0 )
   {
      sprintf(buf, "  %-8.8s  %.15" REAL_FORMAT, name1, value1);
      os << buf;

      if( name2 != 0 )
      {
         sprintf(buf, "   %-8.8s  %.15" REAL_FORMAT, name2, value2);
         os << buf;
      }
   }

   os << std::endl;
}



static Real MPSgetRHS(Real left, Real right)
{
   Real rhsval;

   if( left > -infinity ) /// This includes ranges
      rhsval = left;
   else if( right <  infinity )
      rhsval = right;
   else
      throw SPxInternalCodeException("XMPSWR01 This should never happen.");

   return rhsval;
}



static const char* MPSgetRowName(
   const SPxLPBase<Real>& lp,
   int                   idx,
   const NameSet*        rnames,
   char*                 buf
   )
{
   assert(buf != 0);
   assert(idx >= 0);
   assert(idx <  lp.nRows());

   if( rnames != 0 )
   {
      DataKey key = lp.rId(idx);

      if( rnames->has(key) )
         return (*rnames)[key];
   }

   sprintf(buf, "C%d", idx);

   return buf;
}



/// Write LP in MPS format.
/** @note There will always be a BOUNDS section, even if there are no bounds.
 */
template <>
void SPxLPBase<Real>::writeMPS(
   std::ostream&  p_output,          ///< output stream.
   const NameSet* p_rnames,          ///< row names.
   const NameSet* p_cnames,          ///< column names.
   const DIdxSet* p_intvars          ///< integer variables.
   ) const
{

   const char*    indicator;
   char           name [16];
   char           name1[16];
   char           name2[16];
   bool           has_ranges = false;
   int            i;
   int            k;

   // --- NAME Section ---
   p_output << "NAME          MPSDATA" << std::endl;

   // --- ROWS Section ---
   p_output << "ROWS" << std::endl;

   for( i = 0; i < nRows(); i++ )
   {
      if( lhs(i) == rhs(i) )
         indicator = "E";
      else if( (lhs(i) > -infinity) && (rhs(i) < infinity) )
      {
         indicator = "E";
         has_ranges = true;
      }
      else if( lhs(i) > -infinity )
         indicator = "G";
      else if( rhs(i) <  infinity )
         indicator = "L";
      else
         throw SPxInternalCodeException("XMPSWR02 This should never happen.");

      MPSwriteRecord(p_output, indicator, MPSgetRowName(*this, i, p_rnames, name));
   }

   MPSwriteRecord(p_output, "N", "MINIMIZE");

   // --- COLUMNS Section ---
   p_output << "COLUMNS" << std::endl;

   bool has_intvars = (p_intvars != 0) && (p_intvars->size() > 0);

   for( int j = 0; j < (has_intvars ? 2 : 1); j++ )
   {
      bool is_intrun = has_intvars && (j == 1);

      if( is_intrun )
         p_output << "    MARK0001  'MARKER'                 'INTORG'" << std::endl;

      for( i = 0; i < nCols(); i++ )
      {
         bool is_intvar = has_intvars && (p_intvars->number(i) >= 0);

         if( ( is_intrun && !is_intvar) || (!is_intrun &&  is_intvar) )
             continue;

         const SVectorBase<Real>& col = colVector(i);
         int colsize2 = (col.size() / 2) * 2;

         assert(colsize2 % 2 == 0);

         for( k = 0; k < colsize2; k += 2 )
            MPSwriteRecord(p_output, 0, getColName(*this, i, p_cnames, name),
               MPSgetRowName(*this, col.index(k), p_rnames, name1), col.value(k),
               MPSgetRowName(*this, col.index(k + 1), p_rnames, name2), col.value(k + 1));

         if( colsize2 != col.size() )
            MPSwriteRecord(p_output, 0, getColName(*this, i, p_cnames, name),
               MPSgetRowName(*this, col.index(k), p_rnames, name1), col.value(k));

         if( isNotZero(maxObj(i)) )
            MPSwriteRecord(p_output, 0, getColName(*this, i, p_cnames, name), "MINIMIZE", -maxObj(i));
      }

      if( is_intrun )
         p_output << "    MARK0001  'MARKER'                 'INTEND'" << std::endl;
   }

   // --- RHS Section ---
   p_output << "RHS" << std::endl;

   i = 0;
   while( i < nRows() )
   {
      Real rhsval1 = 0.0;
      Real rhsval2 = 0.0;

      for( ; i < nRows(); i++ )
         if( (rhsval1 = MPSgetRHS(lhs(i), rhs(i))) != 0.0 )
            break;

      if( i < nRows() )
      {
         for( k = i + 1; k < nRows(); k++ )
         {
            if( (rhsval2 = MPSgetRHS(lhs(k), rhs(k))) != 0.0 )
               break;
         }

         if( k < nRows() )
         {
            MPSwriteRecord(p_output, 0, "RHS", MPSgetRowName(*this, i, p_rnames, name1), rhsval1,
               MPSgetRowName(*this, k, p_rnames, name2), rhsval2);
         }
         else
            MPSwriteRecord(p_output, 0, "RHS", MPSgetRowName(*this, i, p_rnames, name1), rhsval1);

         i = k + 1;
      }
   }

   // --- RANGES Section ---
   if( has_ranges )
   {
      p_output << "RANGES" << std::endl;

      for( i = 0; i < nRows(); i++ )
      {
         if( (lhs(i) > -infinity) && (rhs(i) < infinity) )
            MPSwriteRecord(p_output, "", "RANGE", MPSgetRowName(*this, i, p_rnames, name1), rhs(i) - lhs(i));
      }
   }

   // --- BOUNDS Section ---
   p_output << "BOUNDS" << std::endl;

   for( i = 0; i < nCols(); i++ )
   {
      // skip variables that do not appear in the objective function or any constraint
      const SVectorBase<Real>& col = colVector(i);

      if( col.size() == 0 && isZero(maxObj(i)) )
         continue;

      if( lower(i) == upper(i) )
      {
         MPSwriteRecord(p_output, "FX", "BOUND", getColName(*this, i, p_cnames, name1), lower(i));
         continue;
      }

      if( (lower(i) <= -infinity) && (upper(i) >= infinity) )
      {
         MPSwriteRecord(p_output, "FR", "BOUND", getColName(*this, i, p_cnames, name1));
         continue;
      }

      if( lower(i) != 0.0 )
      {
         if( lower(i) > -infinity )
            MPSwriteRecord(p_output, "LO", "BOUND", getColName(*this, i, p_cnames, name1), lower(i));
         else
            MPSwriteRecord(p_output, "MI", "BOUND", getColName(*this, i, p_cnames, name1));
      }

      if( has_intvars && (p_intvars->number(i) >= 0) )
      {
         // Integer variables have default upper bound 1.0, but we should write
         // it nevertheless since CPLEX seems to assume infinity otherwise.
         MPSwriteRecord(p_output, "UP", "BOUND", getColName(*this, i, p_cnames, name1), upper(i));
      }
      else
      {
         // Continous variables have default upper bound infinity
         if( upper(i) < infinity )
            MPSwriteRecord(p_output, "UP", "BOUND", getColName(*this, i, p_cnames, name1), upper(i));
      }
   }

   // --- ENDATA Section ---
   p_output << "ENDATA" << std::endl;

   // Output warning when writing a maximisation problem
   if( spxSense() == SPxLPBase<Real>::MAXIMIZE )
   {
      MSG_WARNING( (*spxout), (*spxout) << "XMPSWR03 Warning: objective function inverted when writing maximization problem in MPS file format\n" );
   }
}



// ---------------------------------------------------------------------------------------------------------------------
//  Explicit instantiation
// ---------------------------------------------------------------------------------------------------------------------

template class SPxLPBase < Real >;

} // namespace soplex
