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
#pragma ident "@(#) $Id: spxout.h,v 1.6 2005/11/01 21:27:04 bzforlow Exp $"

/**@file  spxout.h
 * @brief Wrapper for different output streams and verbosity levels.
 */
#ifndef _SPXOUT_H_
#define _SPXOUT_H_

#include <iostream>
#include <iomanip>
#include "spxdefines.h"

// ----------------------------------------------------------------------
//    class SPxOut
// ----------------------------------------------------------------------

namespace soplex 
{

/**@class SPxOut
   @ingroup Elementary

   @brief Wrapper for several output streams. 
   A verbosity level is used to decide which stream to use and whether to
   really print a given message. Regardless of whether the verbosity level
   is set via a manipulator or via the member function, it is persistent
   until a new value is set.

   Most ostream member functions (e.g., @c precision()) are not provided here;
   use the corresponding stream manipulators (e.g., @c setprecision())
   instead. These are passed on to the <em>current</em> ostream, which is
   chosen according to the verbosity level. In particular, this means that the
   first element in an output stream should always be the verbosity. For
   instance, use
   @code
      spxout << verb( SPxOut::WARNING ) << std::setw( 15 ) << 42 << std::endl;
   @endcode
   or 
   @code
      spxout.setVerbosity( SPxOut::WARNING );
      spxout << std::setw( 15 ) << 42 << std::endl;
   @endcode
   instead of
   @code
      spxout << std::setw( 15 ) << verb( SPxOut::WARNING ) << 42 << std::endl;
   @endcode
   in order to make sure that @c std::setw( 15 ) is applied to the warning stream.
*/
class SPxOut
{
public:


   //-----------------------------------
   /**@name Output control types */
   //@{
   /// Verbosity level
   typedef enum 
   {
      // Note: the implementation uses the fact that ERROR == 0 
      // and that the verbosity levels are subsequent numbers.
      // If you change this, change the implementation as well.
      ERROR    = 0, 
      WARNING  = 1,
      VERBOSE1 = 2,
      VERBOSE2 = 3,
      VERBOSE3 = 4,
      DEBUG    = 5
   } Verbosity;

   /// helper struct for the output operator
   struct struct_Verbosity 
   { 
      /// verbosity level
      Verbosity v_; 
   };
   //@}

   //-----------------------------------
   /**@name Construction / destruction */
   //@{
   /// constructor
   SPxOut();
   /// destructor
   virtual ~SPxOut();
   //@}

   //-----------------------------------
   /**@name Verbosity */
   //@{
   ///
   virtual void 
   setVerbosity( const Verbosity& v )
   {
      m_verbosity = v;
   }
   ///
   inline Verbosity
   getVerbosity()
      const
   {
      return m_verbosity;
   }
   //@}

   //----------------------------------------
   /**@name Wrappers for the current stream */
   //@{
   ///
   inline bool good() const
   {
      return getCurrentStream().good();
   }
   ///
   inline bool operator ! () const
   {
      return ! getCurrentStream();
   }
   //@}

   //-----------------------------------
   /**@name Getting / setting streams */
   //@{
   /// Sets the stream for the specified verbosity level.
   virtual void
   setStream( const Verbosity& verbosity,
               std::ostream&   stream )
   {
      m_streams[ verbosity ] = &stream;
   }
   /// Returns the stream for the specified verbosity level.
   inline std::ostream&
   getStream( const Verbosity& verbosity )
      const
   {
      return *(m_streams[ verbosity ]);
   }
   /// Returns the stream for the current verbosity.
   inline std::ostream&
   getCurrentStream()
      const
   {
      return getStream( getVerbosity() );
   }
   //@}

private:

   //-----------------------------------
   /**@name Private data */
   //@{
   /// verbosity level
   Verbosity               m_verbosity;
   /// array of pointers to internal streams, indexed by verbosity level
   std::ostream**          m_streams;
   //@}

   //-----------------------------------
   /**@name Blocked */
   //@{
   /// copy constructor
   SPxOut( const SPxOut& );
   /// assignment operator
   SPxOut& operator=( const SPxOut& );
   //@}
};

   // ---------------------------------------------------------
   //    manipulators
   // ---------------------------------------------------------


   //-------------------------------------------
   /**@name Verbosity manipulator
       This implementation is done similar to the one for setw(), setprecision(),
       etc. in the standard file iomanip. For instance, the non-menber function
       #verb(v) returns a struct struct_Severity which contains only the 
       verbosity level. Calling 
       @code
            SPxOut spxout;
            spxout << verb( SPxOut::ERROR ) << "This is an error!" << std::endl;
       @endcode
       passes such a struct to the output operator defined below, which
       extracts the verbosity level from the struct and passes it to the 
       member function SPxOut::setVerbosity(). 
   */
   //@{
   /// manipulator to be used in an output statement
   inline SPxOut::struct_Verbosity
   verb( const SPxOut::Verbosity&  v )
   {
      SPxOut::struct_Verbosity verbosity;
      verbosity.v_ = v;
      return verbosity;
   }

   /// output operator with verbosity level struct
   inline SPxOut& 
   operator<< ( SPxOut& stream, 
                const SPxOut::struct_Verbosity&  verbosity )
   {
      stream.setVerbosity( verbosity.v_ );
      return stream;
   }
   //@}

   //--------------------------------------------------------
   /**@name Standard manipulators and output of other types */
   //@{
   /// Passes standard manipulators without arguments, like @c std::endl, 
   /// @c std::flush, or @c std::ios::right, to the current stream.
   inline SPxOut&
   operator<< ( SPxOut&       _spxout, 
                std::ostream& (*manip)( std::ostream& ) )
   {
      if ( _spxout.getVerbosity() <= Param::verbose() )
         _spxout.getCurrentStream() << manip;
      return _spxout;
   }

   /// Passes everything else to the current stream. This includes 
   /// basic types (char, char*, long, int,...) as well as structs which
   /// correspond to manipulators with arguments, such as the struct _Setw 
   /// for the setw() manipulator.
   template< typename T >
   inline SPxOut&
   operator<< ( SPxOut& _spxout, const T&  v )
   {
      if ( _spxout.getVerbosity() <= Param::verbose() )
         _spxout.getCurrentStream() << v;
      return _spxout;
   }
   //@}

   // ---------------------------------------------------------
   //    declare global SPxOut instance
   // ---------------------------------------------------------

   //-----------------------------------
   /**@name Global instance */
   //@{
   ///
   extern SPxOut spxout;
   //@}

}      // namespace soplex


#endif // _SPXOUT_H_

//-----------------------------------------------------------------------------
//Emacs Local Variables:
//Emacs mode:c++
//Emacs c-basic-offset:3
//Emacs tab-width:8
//Emacs indent-tabs-mode:nil
//Emacs End:
//-----------------------------------------------------------------------------
   
