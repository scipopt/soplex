#ifndef _TRACEMETHOD_H_
#define _TRACEMETHOD_H_


#include <iostream>


#define FILE_NAME_COL  60

namespace soplex
{

/**@class TraceMethod 
   @brief Helper class to trace the calling of methods.

   Prints out a message with the current filename and line to indicate
   that some method has been called.
*/
class TraceMethod
{
private:

   //------------------------------------
   /**@name Data */
   //@{
   /// the current indentation
   static int s_indent;
   //@}

public:

   //------------------------------------
   /**@name Construction / destruction */
   //@{
   /// constructor
   TraceMethod(const char* s, const char* file, int line );
 
  /// destructor
   virtual ~TraceMethod()
   {
      s_indent--;
   }
   //@}
};


};


#endif
