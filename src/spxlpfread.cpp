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
#pragma ident "@(#) $Id: spxlpfread.cpp,v 1.1 2001/11/19 09:30:47 bzfkocht Exp $"

#include <ctype.h>
#include <assert.h>
#include <iostream>

#include "spxlp.h"

namespace soplex
{
class LPFReader
{
private:
   enum Section { START, OBJECTIVE, LHS, RHS, BOUNDS, INTEGERS, BINARYS };

   static const int MaxLineLen  = 256;    ///< maximum input line length.
   static const int InitCols    = 10000;  ///< initialy allocated columns.
   static const int InitRows    = 10000;  ///< initialy allocated rows.
   static const int InitNzos    = 100000; ///< initialy allocated non zeros.
   static const int InitNameMem = 100000; ///< initialy memory for names.

   istream&        m_input;               ///< where to read input from.
   char            m_line[MaxLineLen];    ///< last read line.
   int             m_lineno;              ///< current line number.
   char*           m_pos;                 ///< current position in line.
   double          m_coeff;               ///< last read coefficient. 
   int             m_colIdx;              ///< index of last read column name.
   int             m_sense;               ///< last read sense.

   LPCol           m_emptyCol;            ///< reusable empty column.
   LPColSet        m_colSet;              ///< the set of columns read.

   LPRow           m_row;                 ///< last assembled row.
   LPRowSet        m_rowSet;              ///< the set of rows read.
   DSVector&       m_vec;                 ///< last assembled vector (from row).

   Section         m_section;

   SPxLP::SPxSense m_optSense;            ///< objective sense (min or max?).
   NameSet*        m_colNames;            ///< column names.
   NameSet*        m_rowNames;            ///< row names.
   DIdxSet*        m_intVars;             ///< integer variables.
   bool            m_freeColNames;        ///< should we free m_colNames.
   bool            m_freeRowNames;        ///< should we free m_rowNames.
   bool            m_error;               ///< was an error encountered.

   /// read a line from the file.
   bool readLine();


   bool isSense();
   bool isCoeff();
   bool isColName();

   bool readSense();

   /// read the next coefficent (if any) from the input.
   bool readCoeff();

   /// read the next identifier from the input.
   bool readColName();

   /// does the line start with a keyoword ?
   int hasKeyword(const char* keyword, const char* buf);

   int hasRowName(const char* buf);

   bool syntaxError();

public:
   /// constructor.
   LPFReader(istream& p_input, NameSet* p_rowNames, 
      NameSet* p_colNames, DIdxSet* p_intVars);

   /// destructor.
   ~LPFReader();

   bool readLP();
};

LPFReader::LPFReader(
   istream& p_input, 
   NameSet* p_rowNames, 
   NameSet* p_colNames, 
   DIdxSet* p_intVars)
   : m_input(p_input)
   , m_lineno(0)
   , m_pos(0)
   , m_emptyCol()
   , m_colSet(InitCols, 1)
   , m_row(100)
   , m_rowSet(InitRows, InitNzos)
   , m_vec(m_row.rowVector())
   , m_section(START)
   , m_optSense(SPxLP::MAXIMIZE)
   , m_colNames(p_colNames)
   , m_rowNames(p_rowNames)
   , m_intVars(p_intVars)
   , m_freeColNames(false)
   , m_freeRowNames(false)
   , m_error(false)
{
   assert(m_input != 0);

   if (m_colNames == 0)
   {
      m_colNames     = new NameSet(InitCols, InitNameMem);
      m_freeColNames = true;
   }
   assert(m_colNames != 0);

   if (m_rowNames == 0)
   {
      m_rowNames     = new NameSet(InitRows, InitNameMem);
      m_freeRowNames = true;
   }
   assert(m_rowNames != 0);
}

LPFReader::~LPFReader()
{
   if (m_freeColNames)
      delete m_colNames;
   if (m_freeRowNames)
      delete m_rowNames;
}

bool LPFReader::readLP()
{
   for(;;)
   {
      if ((m_pos == 0) || (*m_pos == '\0'))
         if (readLine())
            break;

      // Now process the sections 
      switch(m_section)
      {
      case OBJECTIVE :
         if (isCoeff())
            readCoeff();

         if (isColName())
         {
            readColName();
            m_vec.add(m_colIdx, m_coeff);
         }
         break;
      case LHS :
         if (isCoeff())
            readCoeff();

         if (isColName())
         {
            readColName();
            m_vec.add(m_colIdx, m_coeff);
         }
         if (isSense())
         {
            // Hier muesste man wissen of gerade ein ColName gekommen ist
            // andernfalls hat man ein range.
            readSense();
            m_section = RHS;
         }
         break;
      case RHS :
         if (isCoeff())
         {
            readCoeff();

            if (m_sense == '<')
            { 
               m_row.lhs() = -SPxLP::infinity; 
               m_row.rhs() = m_coeff;
            }
            else if (m_sense == '>')
            {
               m_row.lhs() = m_coeff;
               m_row.rhs() = SPxLP::infinity;
            }
            else if (m_sense == '=')
            {
               m_row.lhs() = m_coeff;
               m_row.rhs() = m_coeff;
            }
            else if (m_sense == 'R')
            {
               std::cerr << "Ranges not yet implemented" << std::endl;
               // m_row.lhs() = ???;
               // m_row.rhs() = m_coeff;
            }
            else
               abort();

            // add inequality 
            m_rowSet.add(m_row);
            m_vec.clear();

            m_section = LHS;
            continue;
         }
         return syntaxError();
      case BOUNDS :
         // gucken ob wir schon einen namen haben bzw sense haben und
         // dann entsprechend das richtige tun.
         if (isCoeff())
            readCoeff();

         if (isColName())
         {
            readColName();
         }
         if (isSense())
         {
            readSense();
         }
         break;
      case BINARYS :
         if (isColName())
         {
            readColName();

            if ((m_intVars != 0) && (m_colIdx >= 0))
               m_intVars->addIdx(m_colIdx);
         }
         break;
      case INTEGERS :
         if (isColName())
         {
            readColName();

            if ((m_intVars != 0) && (m_colIdx >= 0))
               m_intVars->addIdx(m_colIdx);
         }
         break;
      case START :
      default :
         abort();
      }
   }
   if (!m_error)
   {
      assert(m_colSet.isConsistent());
      assert(m_rowSet.isConsistent());
   }

   return m_error;
}

bool LPFReader::readLine()
{
   char buf[MaxLineLen];
   char tmp[MaxLineLen];
   int  i;
   int  k;

   do
   {      
      if (m_input.getline(buf, sizeof(buf)) == 0)
         return true;

      m_lineno++;
      i = 0;

      // 1. Remove comments.
      char* s = strchr(buf, '\\');

      if (s != 0)
         *s = '\0';

      // 2. look for keywords. 
      switch(m_section)
      {
      case START :
         if ((i = hasKeyword("max[imize]", buf)) > 0)
         {
            m_optSense = SPxLP::MAXIMIZE;
            m_section  = OBJECTIVE;
            m_vec.clear();
         } 
         else if ((i = hasKeyword("min[imize]", buf)) > 0)
         {
            m_optSense = SPxLP::MINIMIZE;
            m_section  = OBJECTIVE;
            m_vec.clear();
         } 
         break;
      case OBJECTIVE :
         i = hasKeyword("s[ubject][   ]t[o]", buf);
         if (i == 0)
            i = hasKeyword("s[uch][    ]t[hat]", buf);
         if (i == 0)
            i = hasKeyword("s[.][    ]t[.]", buf);
         if (i > 0)
         {
            // store objective vector            
            for(int j = m_vec.size() - 1; j >= 0; --j)
               m_colSet.obj(m_vec.index(j)) = m_vec.value(j);
            m_vec.clear();
            m_section = LHS;
         }
      case LHS :
      case BOUNDS :
      case BINARYS :
      case INTEGERS :
         if (hasKeyword("bounds", buf))
         {
            m_section = BOUNDS;
         } 
         if (hasKeyword("bin[arys]", buf))
         {
            m_section = BINARYS;
         } 
         if (hasKeyword("int[egers]", buf))
         {
            m_section = INTEGERS;
         } 
         if (hasKeyword("end", buf))
         {
            return false;
         } 
         break;
      case RHS :
      default :
         syntaxError();
         return false;
      }
      // 3. look for row names.
      i += hasRowName(&buf[i]);

      // 4. remove spaces.
      for(k = 0; buf[i] != '\0'; i++)
         if (!isspace(buf[i]))
            tmp[k++] = buf[i];

      tmp[k] = '\0';

   } // 5. Is this a empty line ?
   while(tmp[0] == '\0');

   // 6. collapse sequences of '+' and '-'. e.g ++---+ => -
   for(i = 0, k = 0; tmp[i] != '\0'; i++)
   {
      while(((tmp[i] == '+') || (tmp[i] == '-')) 
         && ((tmp[i + 1] == '+') || (tmp[i + 1] == '-')))
      {
         if (tmp[i++] == '-')
            tmp[i] = (tmp[i] == '-') ? '+' : '-';
      }
      m_line[k++] = tmp[i];
   }
   m_line[k] = '\0';

   // 7. We have something to process. 
   m_pos = m_line;

   return false;
}

bool LPFReader::isCoeff()
{
   return isdigit(*m_pos) 
      || (*m_pos == '+') || (*m_pos == '-') || (*m_pos == '.');
}

bool LPFReader::isColName()
{
   return isalpha(*m_pos) || (*m_pos == '_');
}

bool LPFReader::isSense()
{
   return (*m_pos == '<') || (*m_pos == '>') || (*m_pos == '=');
}

bool LPFReader::readCoeff()
{
   if (!isCoeff())
      return syntaxError();

   char        tmp[MaxLineLen];
   const char* s = m_pos;
   char*       t;

   // 1. sign 
   if ((*s == '+') || (*s == '-'))
      s++;

   // 2. Digits before the decimal dot
   while(isdigit(*s))
      s++;

   // 3. Decimal dot
   if (*s == '.')
   {
      s++;

      // 4. If there was a dot, posible digit behind it
      while(isdigit(*s))
         s++;
   }
   // 5. Exponent
   if (tolower(*s) == 'e')
   {
      s++;

      // 6. Exponent sign 
      if ((*s == '+') || (*s == '-'))
         s++;

      // 7. Exponent digits
      while(isdigit(*s))
         s++;      
   }
   if (s == m_pos)
      m_coeff = 1.0;
   else
   {
      for(t = tmp; m_pos != s; m_pos++)
         *t++ = *m_pos;
   
      *t = '\0';
      
      m_coeff = atof(tmp);
   }
   assert(m_pos == s);

   return false;
}

bool LPFReader::readColName()
{
   if (!isColName())
      return syntaxError();

   char        name[MaxLineLen];
   const char* s = m_pos;
   int         i;

   while((*s != '+') && (*s != '-') && (*s != ':') && (*s != '\0'))
      s++;

   for(i = 0; m_pos != s; i++, m_pos++)
      name[i] = *m_pos;

   name[i] = '\0';

   if ((m_colIdx = m_colNames->number(name)) < 0)
   {
      m_colIdx = m_colNames->num();
      m_colNames->add(name);
      m_colSet.add(m_emptyCol);
   }
   return false;
}

bool LPFReader::readSense()
{
   if (m_pos[0] == '<') 
      m_sense = m_pos[0];
   else if (m_pos[0] == '>') 
      m_sense = m_pos[0];
   else if (m_pos[0] != '=')
      return syntaxError();
   else 
   {
      assert(m_pos[0] == '=');

      if (m_pos[1] == '<') 
         m_sense = m_pos[1];
      else if (m_pos[1] == '>') 
         m_sense = m_pos[1];
      else 
         m_sense = m_pos[0];
   }
   m_pos++;

   if ((m_pos[0] == '<') || (m_pos[0] == '>') || (m_pos[0] == '=')) 
      m_pos++;

   return false;
}

#if 0
bool LPFReader::hasKeyword(const char* keyword)
{
   int i;
   int k;

   for(i = 0, k = 0; keyword[i] != '\0'; i++, k++)
   {
      if (keyword[i] == '*')
      {
         i++;
         while((tolower(m_line[k]) != keyword[i]) && (m_line[k] != '\0'))
            k++;
      }
      if (keyword[i] != tolower(m_line[k]))
         break;
   }
   return keyword[i] == '\0';
}
#endif

int LPFReader::hasKeyword(const char* keyword, const char* buf)
{
   int i;
   int k;

   assert(keyword != 0);
   assert(buf     != 0);

   for(i = 0, k = 0; keyword[i] != '\0'; i++, k++)
   {
      if (keyword[i] == '[')
      {
         i++;
         // Here we assumed that we have a ']' for the '['.
         while((tolower(buf[k]) == keyword[i]) && (buf[k] != '\0'))
            k++;
         while(keyword[i] != ']')
            i++;         
         k--;
         i--;
      }
      else
         if (keyword[i] != tolower(buf[k]))
            break;
   }
   return keyword[i] == '\0' ? k : 0;
}

int LPFReader::hasRowName(const char* buf)
{
   assert(buf != 0);

   const char* s = strchr(buf, ':');

   if (s == 0)
      return 0;

   int end = s - buf;
   int srt = end - 1;
      
   for(; srt >= 0; srt--)
      if (isspace(buf[srt]))
         break;

   srt++;

   char name[MaxLineLen]; 
   int  i = srt;
   int  k = 0;

   for(i = srt; i < end; i++)
      name[k++] = buf[i];

   name[k] = '\0';

   m_rowNames->add(name);

   return end + 1;
}

bool LPFReader::syntaxError()
{
   std::cerr << "Syntax error in line " << m_lineno << std::endl;

   m_error = true;

   return true;
}

void SPxLP::readLP(istream& is, NameSet* rn, NameSet* cn, DIdxSet* iv)
{
#if 0
    if( cn )
 colnames = cn ;
    else
 colnames = new NameSet(10000, 100000) ;
    colnames->clear()  ;

    if( rn )
 rownames = rn ;
    else
 rownames = new NameSet(10000, 100000) ;
    rownames->clear()  ;


    clear() ;
    changeSense( optSense ) ;
    assert( isConsistent() ) ;
    addCols ( _colset ) ;
    assert( isConsistent() ) ;
    addRows ( _rowset ) ;

    assert( isConsistent() ) ;

    if( cn == 0 )
 delete colnames ;
    if( rn == 0 )
 delete rownames ;
#endif

    LPFReader reader(is, rn, cn, iv);

    if (reader.readLP())
       std::cerr << "An Error occurred" << std::endl;
    else
       std::cout << "Yipppiiiieeeeeeee" << std::endl;

    abort();
}



} // namespace soplex
