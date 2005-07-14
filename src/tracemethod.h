#ifndef _TRACEMETHOD_H_
#define _TRACEMETHOD_H_


#include <iostream>


#define FILE_NAME_COL  60

namespace soplex
{

class TraceMethod
{
private:
   static int s_indent;

public:
   TraceMethod(const char* s, const char* file, int line )
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
   virtual ~TraceMethod()
   {
      s_indent--;
   }
};


};


#endif
