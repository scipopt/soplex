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
#pragma ident "@(#) $Id: spxlpfread.cpp,v 1.5 2001/11/21 09:30:15 bzfkocht Exp $"

/**@file  spxlpfread.cpp
 * @brief Read LP format files.
 */
#include <assert.h>
#include <ctype.h>
#include <iostream>

#include "spxlp.h"

#define MAX_LINE_LEN  256         ///< maximum length of a line

#define INIT_COLS     10000       ///< initialy allocated columns.
#define INIT_ROWS     10000       ///< initialy allocated rows.
#define INIT_NZOS     100000      ///< initialy allocated non zeros.
#define INIT_NAME_MEM 100000      ///< initialy memory for names.

namespace soplex
{
/// Is \p c a \c space, \c tab, \c nl or \c cr ?
static inline bool isSpace(int c)
{
   return (c == ' ') || (c == '\t') || (c == '\n') || (c == '\r'); 
}

/// Is there a number at the beginning of \p s ?
static bool isValue(const char* s)
{
   return ((*s >= '0') && (*s <= '9'))  
      || (*s == '+') || (*s == '-') || (*s == '.');
}

/// Is there a posiible column name at the beginning of \p s ?
static bool isColName(const char* s)
{
   return ((*s >= 'A') && (*s <= 'Z'))
      || ((*s >= 'a') && (*s <= 'z'))
      || (strchr("!\"#$%&()/,;?@_'`{}|~", *s) != 0);
}

/// Is there a comparison operator at the beginning of \p s ?
static bool isSense(const char* s)
{
   return (*s == '<') || (*s == '>') || (*s == '=');
}

static bool isInfinity(const char* s)
{
   if ((*s == '-') || (*s == '+'))
      s++;

   return (tolower(s[0]) == 'i') 
      && ( tolower(s[1]) == 'n') 
      && ( tolower(s[2]) == 'f');
}

/// Read the next number and advance \p pos.
/** If only a sign is encountered, the number is assumed to 
 *  be \c sign * 1.0. 
 * This routine will not catch malformatted numbers like .e10 !
 */
static double readValue(char*& pos)
{
   assert(isValue(pos));

   char        tmp[MAX_LINE_LEN];
   const char* s = pos;
   char*       t;
   double      value = 1.0;
   bool        has_digits = false;

   // 1. sign 
   if ((*s == '+') || (*s == '-'))
      s++;

   // 2. Digits before the decimal dot
   while((*s >= '0') && (*s <= '9'))
   {
      has_digits = true;
      s++;
   }
   // 3. Decimal dot
   if (*s == '.')
   {
      s++;

      // 4. If there was a dot, posible digit behind it
      while((*s >= '0') && (*s <= '9'))
      {
         has_digits = true;
         s++;
      }
   }
   // 5. Exponent
   if (tolower(*s) == 'e')
   {
      s++;

      // 6. Exponent sign 
      if ((*s == '+') || (*s == '-'))
         s++;

      // 7. Exponent digits
      while((*s >= '0') && (*s <= '9'))
         s++;      
   }
   assert(s != pos);
   
   if (!has_digits)
      value = (*pos == '-') ? -1.0 : 1.0;
   else
   {
      for(t = tmp; pos != s; pos++)
         *t++ = *pos;   
      *t = '\0';
      value = atof(tmp);
   }
   pos += s - pos;

   assert(pos == s);

   // std::cout << "readValue = " << value << std::endl;

   return value;
}

/// Read the next column name from the input.
/** The name read is looked up and if not found \p emptycol
 *  is added to \p colset. \p pos is advanced behind the name.
 *  @return The Index of the named column.
 */
static int readColName(
   char*& pos, NameSet* colnames, LPColSet& colset, LPCol& emptycol)
{
   assert(isColName(pos));
   assert(colnames != 0);

   char        name[MAX_LINE_LEN];
   const char* s = pos;
   int         i;
   int         colidx;

   // This are the characters that are not allowed in a column name.
   while((strchr("+-.<>=", *s) == 0) && (*s != '\0'))
      s++;

   for(i = 0; pos != s; i++, pos++)
      name[i] = *pos;

   name[i] = '\0';

   if ((colidx = colnames->number(name)) < 0)
   {
      colidx = colnames->num();
      colnames->add(name);
      colset.add(emptycol);
   }
   // std::cout << "readColName [" << name << "] = " << colidx << std::endl;

   return colidx;
}

/// Read the next <,>,=,==,<=,=<,>=,=> and advance \p pos.
static int readSense(char*& pos)
{
   assert(isSense(pos));

   int sense = *pos++;

   if ((*pos == '<') || (*pos == '>'))
      sense = *pos++;
   else if (*pos == '=')
      pos++;

   // std::cout << "readSense = " << static_cast<char>(sense) << std::endl;

   return sense;
}

/// Is the \p keyword present in \p buf ? If yes, advance \p pos.
/** \p keyword should be lower case. It can contain optional sections
 *  which are enclosed in '[' ']' like "min[imize]". 
 */
static bool hasKeyword(char*& pos, const char* keyword)
{
   int i;
   int k;

   assert(keyword != 0);

   for(i = 0, k = 0; keyword[i] != '\0'; i++, k++)
   {
      if (keyword[i] == '[')
      {
         i++;

         // Here we assumed that we have a ']' for the '['.
         while((tolower(pos[k]) == keyword[i]) && (pos[k] != '\0'))
         {
           k++;
           i++;
         }
         while(keyword[i] != ']')
            i++;         
         --k;
      }
      else
      {
         if (keyword[i] != tolower(pos[k]))
            break;
      }
   }
   if (keyword[i] == '\0')
   {
      pos += k;

      // std::cout << "hasKeyowrd: " << keyword << std::endl;
      return true;
   }
   return false;
}

/// If \p buf start with "name:" store the name in \p rownames 
/// and advance \p pos. 
static int hasRowName(char*& pos, NameSet* rownames)
{
   assert(rownames != 0);

   const char* s = strchr(pos, ':');

   if (s == 0)
      return false;

   int end = s - pos;
   int srt = end - 1;
      
   for(; srt >= 0; srt--)
      if (pos[srt] == ' ')
         break;

   srt++;

   char name[MAX_LINE_LEN]; 
   int  i = srt;
   int  k = 0;

   for(i = srt; i < end; i++)
      name[k++] = pos[i];

   name[k] = '\0';

   rownames->add(name);

   pos += end + 1;

   return true;
}

static double readInfinity(char*& pos)
{
   assert(isInfinity(pos));

   double sense;

   if ((*pos == '-') || (*pos == '+'))
   {
      sense = (*pos == '-') ? -1.0 : 1.0;
      pos++;
   }
   hasKeyword(pos, "inf[inity]");

   return sense;
}

/// Read LP in "CPLEX LP File Format".
/** 
 *  The specification is taken from the
 *  ILOG CPLEX 7.0 Reference Manual, Appendix E, Page 527.
 *
 *  This routine should read (most?) valid LP format files. 
 *  What it will not do, is find all cases where a file is ill formed. 
 *  If this happens it may complain and read nothing or read "something".
 *
 *  @todo Keyword \c free in BOUNDS section is not implemented.
 *  @todo Keyword \c Infinity in BOUNDS section is not implemented.
 *        Problem is that inf is a valid variable name in the format.
 *        readColName should be changed with a parameter to select, if 
 *        new names should be atomatically added to the column set.
 */  
void SPxLP::readLP(
   istream& p_input, 
   NameSet* p_rnames,               ///< row names.
   NameSet* p_cnames,               ///< column names.
   DIdxSet* p_intvars)              ///< integer variables.
{
   enum 
   { 
      START, OBJECTIVE, CONSTRAINTS, BOUNDS, INTEGERS, BINARYS 
   } section = START;

   NameSet*  rnames;                ///< row names.
   NameSet*  cnames;                ///< column names.

   LPCol     emptycol;              ///< reusable empty column.
   LPColSet  cset;                  ///< the set of columns read.
   LPRow     row;                   ///< last assembled row.
   LPRowSet  rset;                  ///< the set of rows read.
   DSVector& vec = row.rowVector(); ///< last assembled vector (from row).
   double    val = 1.0;
   int       colidx;
   int       sense = 0;

   char      buf [MAX_LINE_LEN];
   char      tmp [MAX_LINE_LEN];
   char      line[MAX_LINE_LEN];
   int       lineno = 0;
   bool      finished = false;
   int       i;
   int       k;
   char*     s;
   char*     pos;

   cnames = (p_cnames != 0) 
      ? p_cnames : new NameSet(INIT_COLS, INIT_NAME_MEM);

   cnames->clear();

   rnames = (p_rnames != 0)
      ? p_rnames : new NameSet(INIT_ROWS, INIT_NAME_MEM);

   rnames->clear();

   clear(); // clear the LP.

   //--------------------------------------------------------------------------
   //--- Main Loop
   //--------------------------------------------------------------------------
   for(;;)
   {      
      // 0. Read a line from the file.
      if (p_input.getline(buf, sizeof(buf)) == 0)
      {
         std::cerr << "No 'End' marker found" << std::endl;
         finished = true;
         break;
      }
      lineno++;
      i = 0;
      pos = buf;
      val = 1.0;

      // cout << "Reading line " << lineno << std::endl;
      // cout << pos << std::endl;

      // 1. Remove comments.
      if (0 != (s = strchr(buf, '\\')))
         *s = '\0';

      // 2. Look for keywords. 
      if (section == START)
      {
         if (hasKeyword(pos, "max[imize]"))
         {
            changeSense(SPxLP::MAXIMIZE);
            section = OBJECTIVE;
         } 
         else if (hasKeyword(pos, "min[imize]"))
         {
            changeSense(SPxLP::MINIMIZE);
            section = OBJECTIVE;
         } 
      }
      else if (section == OBJECTIVE)
      {
         if (hasKeyword(pos, "s[ubject][   ]t[o]")
            || hasKeyword(pos, "s[uch][    ]t[hat]")
            || hasKeyword(pos, "s[.][    ]t[.]"))
         {
            // store objective vector            
            for(int j = vec.size() - 1; j >= 0; --j)
               cset.obj(vec.index(j)) = vec.value(j);
            vec.clear();
            section = CONSTRAINTS;
         }
      }
      else
      {
         if (hasKeyword(pos, "bound[s]"))
            section = BOUNDS;
         else if (hasKeyword(pos, "bin[arys]"))
            section = BINARYS;
         else if (hasKeyword(pos, "gen[erals]"))
            section = INTEGERS;
         else if (hasKeyword(pos, "end"))
         {
            finished = true;
            break;
         }
      }

      // 3. Look for row names.
      if ((section == OBJECTIVE) || (section == CONSTRAINTS))
         hasRowName(pos, rnames);

      // 4. Remove spaces.
      if ((section == BINARYS) || (section == INTEGERS))
      {
         // only inital spaces
         while(isSpace(*pos))
            pos++;

         strcpy(tmp, pos);
      }
      else
      {
         // all spaces
         for(k = 0; pos[i] != '\0'; i++)
            if (!isSpace(pos[i]))
               tmp[k++] = pos[i];

         tmp[k] = '\0';
      }

      // 5. Is this a empty line ?
      if (tmp[0] == '\0')
         continue;

      // 6. Collapse sequences of '+' and '-'. e.g ++---+ => -
      for(i = 0, k = 0; tmp[i] != '\0'; i++)
      {
         while(((tmp[i] == '+') || (tmp[i] == '-')) 
            && ((tmp[i + 1] == '+') || (tmp[i + 1] == '-')))
         {
            if (tmp[i++] == '-')
               tmp[i] = (tmp[i] == '-') ? '+' : '-';
         }
         line[k++] = tmp[i];
      }
      line[k] = '\0';

      //-----------------------------------------------------------------------
      //--- Line processing loop
      //-----------------------------------------------------------------------
      pos = line;
      
      // cout << "we have [" << pos << "]" << std::endl;

      // 7. We have something left to process. 
      while((pos != 0) && (*pos != '\0'))
      {
         // Now process the sections 
         switch(section)
         {
         case OBJECTIVE :
            if (isValue(pos))
               val = readValue(pos);
            
            if (isColName(pos))
            {
               colidx = readColName(pos, cnames, cset, emptycol);
               vec.add(colidx, val);
            }
            break;
         case CONSTRAINTS :
            if (isValue(pos))
            {
               val = readValue(pos);
               
               if (sense != 0)
               {
                  if (sense == '<')
                  { 
                     row.lhs() = -SPxLP::infinity; 
                     row.rhs() = val;
                  }
                  else if (sense == '>')
                  {
                     row.lhs() = val;
                     row.rhs() = SPxLP::infinity;
                  }
                  else 
                  {
                     assert(sense == '=');
                  
                     row.lhs() = val;
                     row.rhs() = val;
                  }
                  rset.add(row);
                  vec.clear();
                  sense = 0;
                  pos   = 0;
                  // next line
                  continue;
               }         
            }
            if (isColName(pos))
            {
               colidx = readColName(pos, cnames, cset, emptycol);

               if (val != 0.0)
                  vec.add(colidx, val);
            }
            if (isSense(pos))
               sense = readSense(pos);
            break;
         case BOUNDS :
            sense = 0;

            if (isValue(pos))
            {
               val = readValue(pos);
               
               if (!isSense(pos))
                  goto syntax_error;

               sense = readSense(pos);
            }
            if (!isColName(pos))
               goto syntax_error;
            
            colidx = readColName(pos, cnames, cset, emptycol);
            
            if (sense)
            {
               if (sense == '<') 
                  cset.lower(colidx) = val;
               else if (sense == '>')
                  cset.upper(colidx) = val;
               else
               {
                  assert(sense == '=');
                  cset.lower(colidx) = val;
                  cset.upper(colidx) = val;
               }
            }
            if (isSense(pos))
            {
               sense = readSense(pos);
               
               if (!isValue(pos))
                  goto syntax_error;
               
               val = readValue(pos);
               
               if (sense == '<') 
                  cset.upper(colidx) = val;
               else if (sense == '>')
                  cset.lower(colidx) = val;
               else
               {
                  assert(sense == '=');
                  cset.lower(colidx) = val;
                  cset.upper(colidx) = val;
               }
            }
            break;
         case BINARYS  :
         case INTEGERS :
            colidx = readColName(pos, cnames, cset, emptycol);

            if (section == BINARYS)
            {
               cset.lower(colidx) = 0.0;
               cset.upper(colidx) = 1.0;
            }
            if (p_intvars != 0)
               p_intvars->addIdx(colidx);

            while(isSpace(*pos))
               pos++;
            break;
         default :
            abort();
         }
      }
   }
   //--------------------------------------------------------------------------

   assert(isConsistent());

   addCols(cset);
   assert(isConsistent());

   addRows(rset); 
   assert(isConsistent());

 syntax_error:
   if (finished)
      std::cout << "Finished reading " << lineno << " lines" << std::endl;
   else
      std::cerr << "Syntax error in line " << lineno << std::endl;

   if (p_cnames == 0)
      delete cnames;
   if (p_rnames == 0)
      delete rnames;

   // std::cout << *this;
}
} // namespace soplex







