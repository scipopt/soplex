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
#pragma ident "@(#) $Id: spxmpsread.cpp,v 1.2 2001/12/10 19:05:54 bzfkocht Exp $"

/**@file  spxmpsread.cpp
 * @brief Read MPS format files.
 */
#include <assert.h>
#include <string.h>
#include <iostream>

#include "spxlp.h"

#define MAX_LINE_LEN  256

#define INIT_COLS     10000       ///< initialy allocated columns.
#define INIT_ROWS     10000       ///< initialy allocated rows.
#define INIT_NZOS     100000      ///< initialy allocated non zeros.
#define INIT_NAME_MEM 100000      ///< initialy memory for names.

#define PATCH_CHAR    '_'
#define BLANK         ' '

#define DEBUG 0                   // Setting this generates a lot of output

namespace soplex
{ 

class MPSInput
{
public:
   enum Section
   {
      NAME, OBJSEN, OBJNAME, ROWS, COLUMNS, RHS, RANGES, BOUNDS, ENDATA
   };

private:
   Section         m_section;
   std::istream&   m_input;
   int             m_lineno;
   SPxLP::SPxSense m_objsense;
   bool            m_has_error;
   char            m_buf[MAX_LINE_LEN];
   const char*     m_f0;
   const char*     m_f1;
   const char*     m_f2;
   const char*     m_f3;
   const char*     m_f4;
   const char*     m_f5;
   char            m_probname[MAX_LINE_LEN];
   char            m_objname [MAX_LINE_LEN];

public:
   MPSInput(std::istream& p_input)
      : m_section(NAME)
      , m_input(p_input)
      , m_lineno(0)
      , m_objsense(SPxLP::MINIMIZE)
      , m_has_error(false)
   {
      m_probname[0] = '\0';
      m_objname [0] = '\0';
   }
   Section         section()  const { return m_section; }
   int             lineno()   const { return m_lineno; }
   const char*     field0()   const { return m_f0; }
   const char*     field1()   const { return m_f1; }
   const char*     field2()   const { return m_f2; }
   const char*     field3()   const { return m_f3; }
   const char*     field4()   const { return m_f4; }
   const char*     field5()   const { return m_f5; }
   const char*     probName() const { return m_probname; }
   const char*     objName()  const { return m_objname; }
   SPxLP::SPxSense objSense() const { return m_objsense; }
   bool            hasError() const { return m_has_error; }

   void setSection(Section p_section)
   {
      m_section = p_section;
   }
   void setProbName(const char* p_probname)
   {
      strcpy(m_probname, p_probname);
   }
   void setObjName(const char* p_objname)
   {
      strcpy(m_objname, p_objname);
   }
   void setObjSense(SPxLP::SPxSense sense)
   {
      m_objsense = sense;
   }
   void syntaxError() 
   {
      std::cerr << "Syntax error in line " << m_lineno << std::endl;
      m_section = ENDATA;
      m_has_error = true;
   }
   void entryIgnored(
      const char* what, const char* what_name, 
      const char* row_name)
   {
      std::cerr << "Warning line " << m_lineno << ": "
                << what << " \"" << what_name << "\"" 
                << " for row \"" 
                << row_name << "\" ignored" << std::endl;
   }
   bool readLine();
   void insertName(const char* name);
};

/// fill the line from \p pos up to column 80 with blanks.
static void clear_from(char* buf, int pos)
{
   for(int i = pos; i < 80; i++)
      buf[i] = BLANK;
   buf[80] = '\0';
}

/// change all blanks inside a field to #PATCH_CHAR.
static void patch_field(char* buf, int beg, int end)
{
   while((beg <= end) && (buf[end] == BLANK))
      end--;

   while((beg <= end) && (buf[beg] == BLANK))
      beg++;

   for(int i = beg; i <= end; i++)
      if (buf[i] == BLANK)
         buf[i] = PATCH_CHAR;
}

/// read a MPS format data line and parse the fields.
bool MPSInput::readLine()
{
   int   len;
   int   space;
   char* s;
   bool  is_marker;

   do
   {
      m_f0 = m_f1 = m_f2 = m_f3 = m_f4 = m_f5 = 0;
      is_marker = false;
   
      // Read until we have a not comment line.
      do
      {
         if (m_input.getline(m_buf, sizeof(m_buf)) == 0)
            return false;
        m_lineno++;

#if DEBUG
        std::cout << "Line " << m_lineno << " " << m_buf << std::endl;
#endif // DEBUG

      } 
      while(*m_buf == '*');

      /* Normalize line
       */
      len = strlen(m_buf);

      for(int i = 0; i < len; i++)
         if ((m_buf[i] == '\t') || (m_buf[i] == '\n') || (m_buf[i] == '\r'))
            m_buf[i] = BLANK;
      
      if (len < 80)
         clear_from(m_buf, len);

      assert(strlen(m_buf) >= 80);

      /* Look for new section
       */
      if (*m_buf != BLANK)
      {
         m_f0 = strtok(&m_buf[0], " ");

         assert(m_f0 != 0);

         m_f1 = strtok(0, " ");

         return true;
      }

      /* Test for fixed format comments
       */
      if ((m_buf[14] == '$') && (m_buf[13] == ' '))
         clear_from(m_buf, 14);
      else if ((m_buf[39] == '$') && (m_buf[38] == ' '))
         clear_from(m_buf, 39);

      /* Test for fixed format
       */
      space = m_buf[12] | m_buf[13] | m_buf[14] 
         | m_buf[22] | m_buf[23] | m_buf[24] 
         | m_buf[36] | m_buf[37] | m_buf[38] | m_buf[39]
         | m_buf[61] | m_buf[62] | m_buf[63] | m_buf[64];
      
      if (space == BLANK)
      {
         /* We assume fixed format, so we patch possible embedded spaces.
          */
         patch_field(m_buf,  4, 12);
         patch_field(m_buf, 14, 22);
         patch_field(m_buf, 39, 47);
      }
      s = &m_buf[1];
      
      /* At this point it is not clear if we have a indicator field.
       * If there is none (e.g. empty) f1 will be the first name field.
       * If there is one, f2 will be the first name field.
       * 
       * Initially comment marks '$' ar only allowed in the beginning
       * of the 2nd and 3rd name field. We test all fields but the first.
       * This makes no difference, since if the $ is at the start of a value
       * field, the line will be errornous anyway.
       */
      do
      {
         if (0 == (m_f1 = strtok(s, " ")))
            break;
         
         if ((0 == (m_f2 = strtok(0, " "))) || (*m_f2 == '$'))
         {
            m_f2 = 0;
            break;      
         }
         if (!strcmp(m_f2, "'MARKER'"))
            is_marker = true;
            
         if ((0 == (m_f3 = strtok(0, " "))) || (*m_f3 == '$'))
         {
            m_f3 = 0;
            break;      
         }
         if (!strcmp(m_f3, "'MARKER'"))
            is_marker = true;

         if ((0 == (m_f4 = strtok(0, " "))) || (*m_f4 == '$'))
         {
            m_f4 = 0;
            break;      
         }
         if ((0 == (m_f5 = strtok(0, " "))) || (*m_f5 == '$'))
            m_f5 = 0;
      }
      while(false);
   }
   while(is_marker);

#if DEBUG
   std::cout << "-----------------------------------------------" << std::endl;
   std::cout << "f0=" << ((m_f0 == 0) ? "nil" : m_f0) << std::endl; 
   std::cout << "f1=" << ((m_f1 == 0) ? "nil" : m_f1) << std::endl; 
   std::cout << "f2=" << ((m_f2 == 0) ? "nil" : m_f2) << std::endl; 
   std::cout << "f3=" << ((m_f3 == 0) ? "nil" : m_f3) << std::endl; 
   std::cout << "f4=" << ((m_f4 == 0) ? "nil" : m_f4) << std::endl; 
   std::cout << "f5=" << ((m_f5 == 0) ? "nil" : m_f5) << std::endl; 
   std::cout << "-----------------------------------------------" << std::endl;
#endif
   return true;
}

/// Insert \p name as field 1 and shift all other fields up.
void MPSInput::insertName(const char* name)
{
   m_f5 = m_f4;
   m_f4 = m_f3;
   m_f3 = m_f2;
   m_f2 = m_f1;
   m_f1 = name;
}

/// Process NAME section.
static void readName(MPSInput& mps)
{
   do
   {
      // This has to be the Line with the NAME section.
      if (!mps.readLine() 
         || (mps.field0() == 0) || (strcmp(mps.field0(), "NAME")))
         break;

      mps.setProbName(mps.field1());

      std::cout << "Problem name   : " << mps.field1() << std::endl;
 
      // This hat to be a new section
      if (!mps.readLine() || (mps.field0() == 0))
         break;

      if (!strcmp(mps.field0(), "ROWS"))
         mps.setSection(MPSInput::ROWS);
      else if (!strcmp(mps.field0(), "OBJSEN"))
         mps.setSection(MPSInput::OBJSEN);
      else if (!strcmp(mps.field0(), "OBJNAME"))
         mps.setSection(MPSInput::OBJNAME);
      else
         break;

      return;
   }
   while(false);

   mps.syntaxError();
}

/// Process OBJSEN section. This Section is an ILOG extension.
static void readObjsen(MPSInput& mps)
{
   do
   {
      // This has to be the Line with MIN or MAX.
      if (!mps.readLine() || (mps.field1() == 0))
         break;

      if (strcmp(mps.field1(), "MIN"))
         mps.setObjSense(SPxLP::MINIMIZE);
      else if (strcmp(mps.field1(), "MAX"))
         mps.setObjSense(SPxLP::MAXIMIZE);
      else
         break;

      // Look for ROWS or OBJNAME Section
      if (!mps.readLine() || (mps.field0() == 0))
         break;

      if (!strcmp(mps.field0(), "ROWS"))
         mps.setSection(MPSInput::ROWS);
      else if (!strcmp(mps.field0(), "OBJNAME"))
         mps.setSection(MPSInput::OBJNAME);
      else
         break;

      return;
   }
   while(false);

   mps.syntaxError();
}

/// Process OBJNAME section. This Section is an ILOG extension.
static void readObjname(MPSInput& mps)
{
   do
   {
      // This has to be the Line with the name.
      if (!mps.readLine() || (mps.field1() == 0))
         break;

      mps.setObjName(mps.field1());

      // Look for ROWS Section
      if (!mps.readLine() || (mps.field0() == 0))
         break;

      if (strcmp(mps.field0(), "ROWS"))
         break;

      mps.setSection(MPSInput::ROWS);
      return;
   }
   while(false);

   mps.syntaxError();
}

/// Process ROWS section. 
static void readRows(
   MPSInput& mps,
   LPRowSet& rset,
   NameSet&  rnames)
{
   LPRow row;

   while(mps.readLine())
   {
      if (mps.field0() != 0)
      {
         std::cout << "Objective name : " << mps.objName() << std::endl;

         if (strcmp(mps.field0(), "COLUMNS"))
            break;

         mps.setSection(MPSInput::COLUMNS);
         return;
      }
      if (*mps.field1() == 'N')
      {
         if (*mps.objName() == '\0')
            mps.setObjName(mps.field2());
      }
      else
      {
         if (rnames.has(mps.field2()))
            break;

         rnames.add(mps.field2());
            
         switch(*mps.field1())
         {
         case 'G' :
            row.lhs() = 0.0;
            row.rhs() = SPxLP::infinity;
            break;
         case 'E' :
            row.lhs() = 0.0;
            row.rhs() = 0.0;
            break;
         case 'L' :
            row.lhs() = -SPxLP::infinity;
            row.rhs() = 0.0;
            break;
         default :
            mps.syntaxError();
            return;
         }
         rset.add(row);
      }
      assert((*mps.field1() == 'N') 
         || (rnames.number(mps.field2()) == rset.num() - 1));
   }
   mps.syntaxError();
}

/// Process COLUMNS section. 
static void readCols(
   MPSInput& mps,
   LPRowSet& rset,
   NameSet&  rnames,
   LPColSet& cset,
   NameSet&  cnames)
{
   double val;
   int    idx;
   char   colname[MAX_LINE_LEN] = { '\0' };
   LPCol  col(rset.num());

   col.obj()   = 0.0;
   col.lower() = 0.0;
   col.upper() = SPxLP::infinity;
   col.colVector().clear();

   while(mps.readLine())
   {
      if (mps.field0() != 0)
      {
         if (strcmp(mps.field0(), "RHS"))
            break;

         if (colname[0] != '\0')
            cset.add(col);

         mps.setSection(MPSInput::RHS);
         return;
      }
      if ((mps.field1() == 0) || (mps.field2() == 0) || (mps.field3() == 0))
         break;

      // new column?
      if (strcmp(colname, mps.field1()))
      {
         // first column?
         if (colname[0] != '\0')
            cset.add(col);

         strcpy(colname, mps.field1());
         cnames.add(colname);
         col.colVector().clear();
         col.obj() = 0.0;
      }
      val = atof(mps.field3());

      if (!strcmp(mps.field2(), mps.objName()))
         col.obj() = val;
      else 
      {
         if ((idx = rnames.number(mps.field2())) < 0)
            mps.entryIgnored("Column", mps.field1(), mps.field2());
         else
            if (val != 0.0)
               col.colVector().add(idx, val);
      }
      if (mps.field5() != 0)
      {
         assert(mps.field4() != 0);

         val = atof(mps.field5());

         if (!strcmp(mps.field4(), mps.objName()))
            col.obj() = val;
         else 
         {
            if ((idx = rnames.number(mps.field4())) < 0)
               mps.entryIgnored("Column", mps.field1(), mps.field4());
            else
               if (val != 0.0)
                  col.colVector().add(idx, val);
         }
      }
   }
   mps.syntaxError();
}

/// Process RHS section. 
static void readRhs(
   MPSInput& mps,
   LPRowSet& rset,
   NameSet&  rnames)
{
   char   rhsname[MAX_LINE_LEN] = { '\0' };
   int    idx;
   double val;

   while(mps.readLine())
   {
      if (mps.field0() != 0)
      {
         std::cout << "RHS name       : " << rhsname << std::endl;

         if (!strcmp(mps.field0(), "RANGES"))
            mps.setSection(MPSInput::RANGES);
         else if (!strcmp(mps.field0(), "BOUNDS"))
            mps.setSection(MPSInput::BOUNDS);
         else if (!strcmp(mps.field0(), "ENDATA"))
            mps.setSection(MPSInput::ENDATA);
         else
            break;

         return;
      }
      if (((mps.field2() != 0) && (mps.field3() == 0))
         || ((mps.field4() != 0) && (mps.field5() == 0)))
         mps.insertName("_RHS_");
      
      if ((mps.field1() == 0) || (mps.field2() == 0) || (mps.field3() == 0))
         break;

      if (*rhsname == '\0')
         strcpy(rhsname, mps.field1());
      
      if (!strcmp(rhsname, mps.field1()))
      {
         if ((idx = rnames.number(mps.field2())) < 0)
            mps.entryIgnored("RHS", mps.field1(), mps.field2());
         else
         {
            val = atof(mps.field3());

            if (rset.rhs(idx) < SPxLP::infinity)
               rset.rhs(idx) = val;
            if (rset.lhs(idx) > -SPxLP::infinity)
               rset.lhs(idx) = val;
         }
      }
      if (mps.field5() != 0)
      {
         if ((idx = rnames.number(mps.field4())) < 0)
            mps.entryIgnored("RHS", mps.field1(), mps.field4());
         else
         {
            val = atof(mps.field5());
         
            if (rset.rhs(idx) < SPxLP::infinity)
               rset.rhs(idx) = val;
            if (rset.lhs(idx) > -SPxLP::infinity)
               rset.lhs(idx) = val;
         }
      }
   }
   mps.syntaxError();
}

/// Process RANGES section. 
static void readRanges(
   MPSInput& mps,
   LPRowSet& rset,
   NameSet&  rnames)
{
   char   rngname[MAX_LINE_LEN] = { '\0' };
   int    idx;
   double val;

   while(mps.readLine())
   {
      if (mps.field0() != 0)
      {
         std::cout << "Range name     : " << rngname << std::endl;

         if (!strcmp(mps.field0(), "BOUNDS"))
            mps.setSection(MPSInput::BOUNDS);
         else if (!strcmp(mps.field0(), "ENDATA"))
            mps.setSection(MPSInput::ENDATA);
         else
            break;

         return;
      }
      if (((mps.field2() != 0) && (mps.field3() == 0))
         || ((mps.field4() != 0) && (mps.field5() == 0)))
         mps.insertName("_RNG_");

      if ((mps.field1() == 0) || (mps.field2() == 0) || (mps.field3() == 0))
         break;

      if (*rngname == '\0')
         strcpy(rngname, mps.field1());
      
      if (!strcmp(rngname, mps.field1()))
      {
         if ((idx = rnames.number(mps.field2())) < 0)
            mps.entryIgnored("Range", mps.field1(), mps.field2());
         else
         {
            val = atof(mps.field3());

            if (val >= 0)
            {
               if (rset.lhs(idx) > -SPxLP::infinity)
                  rset.rhs(idx) = rset.lhs(idx) + val;
               else
                  rset.lhs(idx) = rset.rhs(idx) - val;
            }
            else
            {
               assert(rset.rhs(idx) == rset.lhs(idx));
               rset.lhs(idx) += val;
            }
         }
         if (mps.field5() != 0)
         {
            if ((idx = rnames.number(mps.field4())) < 0)
               mps.entryIgnored("Range", mps.field1(), mps.field4());
            else
            {
               val = atof(mps.field5());

               if (val >= 0)
               {
                  if (rset.lhs(idx) > -SPxLP::infinity)
                     rset.rhs(idx) = rset.lhs(idx) + val;
                  else
                     rset.lhs(idx) = rset.rhs(idx) - val;
               }
               else
               {
                  assert(rset.rhs(idx) == rset.lhs(idx));
                  rset.lhs(idx) += val;
               }
            }
         }
      }
   }
   mps.syntaxError();
}

/// Process BOUNDS section. 
static void readBounds(
   MPSInput& mps,
   LPColSet& cset,
   NameSet&  cnames)
{
   char   bndname[MAX_LINE_LEN] = { '\0' };
   int    idx;
   double val;

   while(mps.readLine())
   {
      if (mps.field0() != 0)
      {
         std::cout << "Bound name     : " << bndname << std::endl;         

         if (strcmp(mps.field0(), "ENDATA"))
            break;

         mps.setSection(MPSInput::ENDATA);
         return;
      }
      if ((mps.field2() != 0) && (mps.field3() == 0))
         mps.insertName("_BND_");

      if ((mps.field1() == 0) || (mps.field2() == 0) || (mps.field3() == 0))
         break;

      if (*bndname == '\0')
         strcpy(bndname, mps.field2());
      
      if (!strcmp(bndname, mps.field2()))
      {
         idx = cnames.number(mps.field3());         
         val = (mps.field4() == 0) ? 0.0 : atof(mps.field4());

         switch(*mps.field1())
         {
         case 'L':
            cset.lower(idx) = val;
            break;
         case 'U':
            cset.upper(idx) = val;
            break;
         case 'F':
            if (mps.field1()[1] == 'X')
            {
               cset.lower(idx) = val;
               cset.upper(idx) = val;
            }
            else
            {
               cset.lower(idx) = -SPxLP::infinity;
               cset.upper(idx) = SPxLP::infinity;
            }
            break;
         case 'M':
            cset.lower(idx) = -SPxLP::infinity;
            break;
         case 'P':
            cset.upper(idx) = SPxLP::infinity;
            break;
         default:
            mps.syntaxError();
            return;
         }
      }
   }
   mps.syntaxError();
}

/// Read LP in "MPS File Format".
/** 
 *  The specification is taken from the
 *
 *  IBM Optimization Library Guide and Reference
 *
 *  Online available at http://www.software.ibm.com/sos/features/libuser.htm
 *
 *  and from the 
 *
 *  ILOG CPLEX 7.0 Reference Manual, Appendix E, Page 531.
 *
 *  This routine should read all valid MPS format files. 
 *  What it will not do, is find all cases where a file is ill formed. 
 *  If this happens it may complain and read nothing or read "something".
 */  
void SPxLP::readMPS(
   std::istream& p_input, 
   NameSet*      p_rnames,          ///< row names.
   NameSet*      p_cnames)          ///< column names.
{
   LPRowSet& rset = *this;
   LPColSet& cset = *this;
   NameSet*  rnames;                ///< row names.
   NameSet*  cnames;                ///< column names.

   cnames = (p_cnames != 0) 
      ? p_cnames : new NameSet(INIT_COLS, INIT_NAME_MEM);

   cnames->clear();

   rnames = (p_rnames != 0)
      ? p_rnames : new NameSet(INIT_ROWS, INIT_NAME_MEM);

   rnames->clear();

   clear(); // clear the LP.

   cset.memRemax(INIT_NZOS);
   cset.reMax(INIT_COLS);

   MPSInput mps(p_input);

   readName(mps);

   if (mps.section() == MPSInput::OBJSEN)
      readObjsen(mps);

   if (mps.section() == MPSInput::OBJNAME)
      readObjname(mps);

   if (mps.section() == MPSInput::ROWS)
      readRows(mps, rset, *rnames);

   addedRows(rset.num());

   if (mps.section() == MPSInput::COLUMNS)
      readCols(mps, rset, *rnames, cset, *cnames);

   if (mps.section() == MPSInput::RHS)
      readRhs(mps, rset, *rnames);

   if (mps.section() == MPSInput::RANGES)
      readRanges(mps, rset, *rnames);

   if (mps.section() == MPSInput::BOUNDS)
      readBounds(mps, cset, *cnames);

   if (mps.section() != MPSInput::ENDATA)
      mps.syntaxError();

   if (mps.hasError())
      clear();
   else
   {
      changeSense(mps.objSense());

      std::cout << "Objective sense: " 
                << ((mps.objSense() == MINIMIZE) ? "Minimize" : "Maximize") 
                << std::endl;         

      added2Set(
         *(reinterpret_cast<SVSet*>(static_cast<LPRowSet*>(this))), 
         *(reinterpret_cast<SVSet*>(static_cast<LPColSet*>(this))), 
         cset.num());
      addedCols(cset.num());
      assert(isConsistent());
   }
}

} // namespace soplex

//-----------------------------------------------------------------------------
//Emacs Local Variables:
//Emacs mode:c++
//Emacs c-basic-offset:3
//Emacs tab-width:8
//Emacs indent-tabs-mode:nil
//Emacs End:
//-----------------------------------------------------------------------------








