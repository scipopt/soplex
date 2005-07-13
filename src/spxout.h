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
#pragma ident "@(#) $Id: spxout.h,v 1.1 2005/07/13 14:18:50 bzforlow Exp $"

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

   This is a wrapper for several output streams, where a verbosity level is
   used to decide which stream to use and whether to really print a given
   message. Regardless of whether the verbosity level is set via a manipulator
   or via the member function, it is persistent until a new value is set.

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
      spxout.set_verbosity( SPxOut::WARNING );
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
   /** Output control types */
   //@{
   /// Verbosity level
   typedef enum 
   {
      // Note: the implementation uses the fact that ERROR == 0 
      // and that the verbosity levels are subsequent numbers.
      // If you change this, change the implementation as well.
      ERROR    = 0, 
      WARNING  = 1,
      INFO1    = 2,
      INFO2    = 3,
      INFO3    = 4,
      DEBUG    = 5
   } 
   Verbosity;

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
   SPxOut()
      : m_verbosity( ERROR )
      , m_streams( new std::ostream*[ DEBUG+1 ] )
   {
      m_streams[ ERROR ] = m_streams[ WARNING ] = &std::cerr;
      for ( int i = INFO1; i <= DEBUG; ++i )
         m_streams[ i ] = &std::cout;
   }
   /// destructor
   virtual ~SPxOut()
   {
      delete [] m_streams;
   }
   //@}

   //-----------------------------------
   /**@name Verbosity */
   //@{
   ///
   virtual void 
   set_verbosity( const Verbosity v )
   {
      m_verbosity = v;
   }
   ///
   inline Verbosity
   get_verbosity()
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
      return get_current_stream().good();
   }
   ///
   inline bool operator ! () 
   {
      return ! get_current_stream();
   }
   //@}

   //-----------------------------------
   /**@name Getting / setting streams */
   //@{
   /// Sets the stream for the specified verbosity level.
   virtual void
   set_stream( const Verbosity verbosity,
               std::ostream&   stream )
   {
      m_streams[ verbosity ] = &stream;
   }
   /// Returns the stream for the specified verbosity level.
   inline std::ostream&
   get_stream( const Verbosity verbosity )
      const
   {
      return *(m_streams[ verbosity ]);
   }
   /// Returns the stream for the current verbosity.
   inline std::ostream&
   get_current_stream()
      const
   {
      return get_stream( get_verbosity() );
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
   /**@name Verbosity manipulator */
   //@{
   /// This implementation is done similar to the one for setw(), setprecision(),
   /// etc. in the standard file iomanip. For instance, the non-menber function
   /// verb( s ) returns a struct struct_Severity which contains only the 
   /// verbosity level. Calling 
   /// @code
   ///      SPxOut spxout;
   ///      spxout << verb( SPxOut::ERROR ) << "This is an error!" << std::endl;
   /// @endcode
   /// passes such a struct to the output operator defined below, which
   /// extracts the verbosity level from the struct and passes it to the 
   /// member function SPxOut::set_verbosity(). 

   /// manipulator to be used in an output statement
   inline SPxOut::struct_Verbosity
   verb( SPxOut::Verbosity v )
   {
      SPxOut::struct_Verbosity verb;
      verb.v_ = v;
      return verb;
   }

   /// output operator with verbosity level struct
   inline SPxOut& 
   operator<< ( SPxOut& stream, 
                const SPxOut::struct_Verbosity verb )
   {
      stream.set_verbosity( verb.v_ );
      return stream;
   }
   //@}

   //--------------------------------------------------------
   /**@name Standard manipulators and output of other types */
   //@{
#define SEND_TO_CURRENT_STREAM( spxout, v ) \
   if ( spxout.get_verbosity() <= Param::verbose() )\
      spxout.get_current_stream() << v; \
   return spxout;

   /// Passes standard manipulators without arguments, like @c std::endl, 
   /// @c std::flush, or @c std::ios::right, to the current stream.
   inline SPxOut&
   operator<< ( SPxOut&       spxout, 
                std::ostream& (*manip)( std::ostream& ) )
   {
      SEND_TO_CURRENT_STREAM( spxout, manip );
   }

   /// Passes everything else to the current stream. This includes 
   /// basic types (char, char*, long, int,...) as well as structs which
   /// correspond to manipulators with arguments, such as the struct _Setw 
   /// for the setw() manipulator.
   template< typename T >
   inline SPxOut&
   operator<< ( SPxOut& spxout, T v )
   {
      SEND_TO_CURRENT_STREAM( spxout, v );
   }
   //@}

   // ---------------------------------------------------------
   //    static SPxOut instance
   // ---------------------------------------------------------

   static soplex::SPxOut s_spxout;

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
   
