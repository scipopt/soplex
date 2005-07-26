/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the class library                   */
/*       SoPlex --- the Sequential object-oriented simPlex.                  */
/*                                                                           */
/*    Copyright (C) 1997-1999 Roland Wunderling                              */
/*                  1997-2002 Konrad-Zuse-Zentrum                            */
/*                            fuer Informationstechnik Berlin                */
/*                                                                           */
/*  SoPlex is distributed under the terms of the ZIB Academic Licence.       */
/*                                                                           */
/*  You should have received a copy of the ZIB Academic License              */
/*  along with SoPlex; see the file COPYING. If not email to soplex@zib.de.  */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
#pragma ident "@(#) $Id: mpsinput.h,v 1.9 2005/07/26 17:02:21 bzforlow Exp $"

/**@file  mpsinput.h
 * @brief Read MPS format files.
 */
#include <iostream>

#include "spxlp.h"


namespace soplex
{ 

/**@class MPSInput
   
   Reads MPS input files. A right-hand side for the objective function is
   allowed but ignored.
 */
class MPSInput
{
public:

   //-----------------------------------
   /**@name Types */
   //@{
   ///
   enum Section
   {
      NAME, OBJSEN, OBJNAME, ROWS, COLUMNS, RHS, RANGES, BOUNDS, ENDATA
   };
   ///
   enum { MAX_LINE_LEN = 256 };
   //@}
   
private:

   //-----------------------------------
   /**@name Private data */
   //@{
   ///
   Section         m_section;
   ///
   std::istream&   m_input;
   ///
   int             m_lineno;
   ///
   SPxLP::SPxSense m_objsense;
   ///
   bool            m_has_error;
   ///
   char            m_buf[MAX_LINE_LEN];
   ///
   const char*     m_f0;
   ///
   const char*     m_f1;
   ///
   const char*     m_f2;
   ///
   const char*     m_f3;
   ///
   const char*     m_f4;
   ///
   const char*     m_f5;
   ///
   char            m_probname[MAX_LINE_LEN];
   ///
   char            m_objname [MAX_LINE_LEN];
   ///
   bool            m_is_integer;
   ///
   bool            m_is_new_format;
   //@}

public:

   //-----------------------------------
   /**@name Construction / destruction */
   //@{
   ///
   explicit
   MPSInput( std::istream& p_input )
      : m_section       ( NAME )
      , m_input         ( p_input )
      , m_lineno        ( 0 )
      , m_objsense      ( SPxLP::MINIMIZE )
      , m_has_error     ( false )
      , m_is_integer    ( false )
      , m_is_new_format ( false )
   {
      m_f0 = m_f1 = m_f2 = m_f3 = m_f4 = m_f5 = 0;

      m_buf     [0] = '\0';
      m_probname[0] = '\0';
      m_objname [0] = '\0';
   }
   ///
   virtual 
   ~MPSInput()
   {}
   //@}
   
   //-----------------------------------
   /**@name Access */
   //@{
   ///
   Section         section()   const { return m_section; }
   ///
   int             lineno()    const { return m_lineno; }
   ///
   const char*     field0()    const { return m_f0; }
   ///
   const char*     field1()    const { return m_f1; }
   ///
   const char*     field2()    const { return m_f2; }
   ///
   const char*     field3()    const { return m_f3; }
   ///
   const char*     field4()    const { return m_f4; }
   ///
   const char*     field5()    const { return m_f5; }
   ///
   const char*     probName()  const { return m_probname; }
   ///
   const char*     objName()   const { return m_objname; }
   ///
   SPxLP::SPxSense objSense()  const { return m_objsense; }
   ///
   bool            hasError()  const { return m_has_error; }
   ///
   bool            isInteger() const { return m_is_integer; }
   //@}

   //-----------------------------------
   /**@name Modification */
   //@{
   ///
   void setSection(Section p_section)
   {
      m_section = p_section;
   }
   ///
   void setProbName(const char* p_probname)
   {
      strcpy(m_probname, p_probname);
   }
   ///
   void setObjName(const char* p_objname)
   {
      strcpy(m_objname, p_objname);
   }
   ///
   void setObjSense(SPxLP::SPxSense sense)
   {
      m_objsense = sense;
   }
   //@}

   //-----------------------------------
   /**@name Warnings and Errors */
   //@{
   ///
   void syntaxError() 
   {
      MSG_ERROR( spxout << "Syntax error in line " << m_lineno << std::endl; )
      m_section = ENDATA;
      m_has_error = true;
   }
   ///
   void entryIgnored(
      const char* what, const char* what_name, 
      const char* entity, const char* entity_name)
   {
      MSG_WARNING( spxout << "Warning: line " << m_lineno << ": "
                          << what << " \"" << what_name << "\"" 
                          << " for " << entity << " \"" 
                          << entity_name << "\" ignored" << std::endl; )
   }
   //@}

   //-----------------------------------
   /**@name Helpers */
   //@{
   /// reads an MPS format data line and parse the fields.
   bool readLine();
   /// Inserts \p name as field 1 and shifts all other fields up.
   void insertName( const char* name, 
                    bool second = false );
   //@}
};
} // namespace soplex

//-----------------------------------------------------------------------------
//Emacs Local Variables:
//Emacs mode:c++
//Emacs c-basic-offset:3
//Emacs tab-width:8
//Emacs indent-tabs-mode:nil
//Emacs End:
//-----------------------------------------------------------------------------


