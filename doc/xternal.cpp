//-----------------------------------------------------------------------------
/**@mainpage SoPlex 
   @version  1.2
   @author   Roland Wunderling
   @author   Andreas Bley
   @author   Tobias Pfender
   @author   Thorsten Koch

   @section The Sequential object-oriented simplex class library.

   The SoPlex class library comprises classes that may be categorized into
   three different types:

   - Elementary classes are provided for general purpose use in
     projects way bejond the scope of numerical software or linear
     programming.
   - Linear algebra classes provide basic data types for (sparse)
     linear algebra computations. However, their functionality is
     restricted to simple operations such as addition and scaling.
     For complex tasks, such as solving linear systems of equations,
     algorithmic classes are provided instead.
   - Algorithmic classes serve for implementing maybe a variaty of
     algorithms for solving numerical (sub-)problems.

  This software has been implemented as a part of Roland Wunderlings 
  Ph.D. thesis 
  "Paralleler und Objektorientierter Simplex-Algorithmus"
  http://www.zib.de/PaperWeb/abstracts/TR-96-09

  SoPlex is implememted in C++. The code should be compliend with the 
  ANSI standard. Exceptions, RTTI and STL are not 
  used in the code. Everything is in one namespace \em soplex.

  Note, that you are allowed to retreive SoPlex only for research purpose as a
  member of a \em noncommercial and \em academic institution.
*/
//-----------------------------------------------------------------------------
/**@namespace soplex
   @brief     Everything should be within this namespace.

   We have put the whole class library in the namespace soplex.
   If anything here is defined outside, this is a mistake and 
   should be reported. 
*/
//-----------------------------------------------------------------------------
/**@defgroup Elementary Elementary Classes
   @brief    General purpose classes.
   
   Elementary classes are provided for general purpose use in
   projects way bejond the scope of numerical software or linear
   programming.
*/
//-----------------------------------------------------------------------------
/**@defgroup Algebra Linear Algebra Classes
   @brief Basic data types for linear algebra computations.
   
   Linear algebra classes provide basic data types for (sparse)
   linear algebra computations. However, their functionality is
   restricted to simple operations such as addition and scaling.
   For complex tasks, such as solving linear systems of equations,
   algorithmic classes are provided instead.
*/
//-----------------------------------------------------------------------------
/**@defgroup Algo Algorithmic Classes
   @brief Implementation of numerical algorithms.   
   
   Algorithmic classes serve for implementing a variaty of
   algorithms for solving numerical (sub-)problems.
*/
//-----------------------------------------------------------------------------
/**@page DataObjects Data Objects 

    \em Data \em objects refer to C++ objects that do not allocate any
    resources, particularly that do not allocate any memory.  This
    makes them behave just like ordinary C structures, in that both,
    the copy constructor and assignment operator are equivalent to a
    memcopy of the memory taken by the object. Examples for data
    objects are all builtin types such as #int or #double or
    \e simple classes such as #complex.
 
    We distinguish \em data \em objects from general C++ objects that
    may include some allocation of resources. (Note, that for general
    C++ objects that do allocate resources, this must be respected by
    providing apropriate copy constructor and assignment operators.)
    An example for a general C++ class is class DataArray.
 
    The distingtion between data and general C++ objects becomes
    relevant when using such objects in container classes such as
    DataArray or Array.  
*/
//-----------------------------------------------------------------------------
/**@page Where is SoPlex running?
   We have tested SoPlex at least to compile with the following 
   Compilers:
   <TABLE>
   <TR><TD>Vendor</TD><TD>Version  </TD><TD>OS                </TD></TR>
   <TR><TD>Gnu   </TD><TD>2.95.3   </TD><TD>SuSE 7.3/x86 Linux</TD></TR>
   <TR><TD>Gnu   </TD><TD>2.96     </TD><TD>SuSE 7.1/AXP Linux</TD></TR>
   <TR><TD>Intel </TD><TD>5.0.1    </TD><TD>SuSE 7.3/x86 Linux</TD></TR>
   <TR><TD>Compaq</TD><TD>6.2-024  </TD><TD>Tru64 5.0         </TD></TR>
   <TR><TD>Compaq</TD><TD>6.3-010  </TD><TD>SuSE 7.1/AXP Linux</TD></TR>
   <TR><TD>Sun   </TD><TD>WS6U2 5.3</TD><TD>Solaris 7         </TD></TR>
   <TR><TD>SGI   </TD><TD>7.3.1.1m </TD><TD>IRIX 6.5          </TD></TR>
   <TR><TD></TD><TD></TD><TD></TD></TR>
   <TR><TD></TD><TD></TD><TD></TD></TR>
   </TABLE>
 */
//-----------------------------------------------------------------------------
/**@page FAQ Frequently Asked Questions

   Here are some answers that can not be answered from the code alone:

   -# Why is <iostream> used but <assert.h> and not <cassert>
     
      The reason is twofold. From the theoretical point we were not
      able to exactly find out in TC++PL in what namespace cassert 
      should load it's declarations. Shurely in std. But since this are
      normally functions with C linkage this won't work.
      The some, like assert are macros, which have no namespace.
      The practical point is, that the compiler vendors seem to be 
      unsure also. Most put everything in both namespaces std and global.
      So there is no advantage in using <cassert>. Compaq even left them
      off because it seemed unclear to them.
      So our reasoning was: If it is a C++ header we use the correct form
      without the .h and in std. If it is a C header, according to the
      standard the .h header has to be there and uses global namespace.
      That seems acceptable to us, especially for C functions.
      
   -# Why is malloc/free sometimes used and not new/delete.
      
      Because there is no realloc with new/delete. Because malloc
      is faster. And we only use it for builtin types.
      If you do not like this descision, it's easy to change spxalloc.h
      to use new/delete.

   -# Can SoPlex solve Integer Programs (IP's).
      
      No. You need an IP-Solver for this. Most IP-Solver use LP-Solvers
      as a subroutine and do some kind of Branch-and-Bound.

   -# Is there a Windows version ?

      The code is tested to compile under some version of Visual C.
      We do \b not provide any make or project files for VC.

   -# What is the academic license ?

      Essentially, you can do what you want, except: Give SoPlex to 
      anybody else. Do anything with SoPlex that makes money (this 
      includes saving your money on any commercial activity).
      See the license for details http://pengpuffblub.html

   -# What is the commercial license ?
  
      If you want to make with SoPlex you need one. We give the license,
      you pay for it. The amount is completely negotiable depending on
      what you want to do with SoPlex, which rights you want and what
      you are willing to tell us.

   -# I want I primal and a dual simplex, where are they ?

      That is quite easy. You can set ENTERing and LEAVEing algorithm and
      COLUMN and ROW representation.

      <TABLE>
      <TR><TD>&nbsp;</TD><TD>ENTER </TD><TD>LEAVE </TD></TR>
      <TR><TD>COLUMN</TD><TD>Primal</TD><TD>Primal</TD></TR>
      <TR><TD>ROW   </TD><TD>Dual  </TD><TD>Dual  </TD></TR>
      </TABLE>

      COLUMN oriented is the "usual" representation.
      then   Entering is the Primal and Leaving is the Dual algorithm.
      In ROW oriented representation, we have in principle the
      explicit dual and then the algorithms reverse.

*/      
      


