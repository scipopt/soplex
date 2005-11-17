

#include "tracemethod.h"
#include "spxout.h"


namespace soplex
{
   /// constructor
   TraceMethod::TraceMethod(const char* s, const char* file, int line )
   {
      int i;
 
      spxout << "\t";
      
      for(i = 0; i < s_indent; i++)
         spxout << ".";      
      
      spxout << s;
      
      for(i = strlen(s) + s_indent; i < FILE_NAME_COL - 8; i++)
         spxout << "_";             
      spxout << "[" << file << ":" << line << "]" << std::endl; 
      s_indent++;
   }
}
