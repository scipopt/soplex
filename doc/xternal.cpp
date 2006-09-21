//-----------------------------------------------------------------------------
/**@mainpage SoPlex 
   @version  1.3.1
   @author   Roland Wunderling
   @author   Andreas Bley
   @author   Tobias Achterberg
   @author   Thorsten Koch
   @author   Benjamin Hiller
   @author   Sebastian Orlowski
   @author   Andreas Tuchscherer

   @section The Sequential object-oriented simplex class library.

   This software has been implemented as a part of Roland Wunderlings 
   Ph.D. thesis "Paralleler und Objektorientierter Simplex-Algorithmus"
   which can be found at http://www.zib.de/PaperWeb/abstracts/TR-96-09
   (in German).

   SoPlex is implemented in C++. The code should be compliant with the 
   current ANSI standard. Exceptions, RTTI and STL (other then iostream) 
   are not used. Everything is in one namespace \em soplex.

   Note, that you are allowed to retrieve SoPlex only for research purpose 
   as a member of a \em non-commercial and \em academic institution.

   - \ref RUN      "Where does it run"

   - \ref INST     "Installation"

   - \ref FLAGS    "How to use the flags in the example program"

   - \ref FAQ      "Frequently asked questions"

   - \ref PROG     "Programming with SoPlex"
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
   projects way beyond the scope of numerical software or linear
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
   
   Algorithmic classes serve for implementing a variety of
   algorithms for solving numerical (sub-)problems.
*/
//-----------------------------------------------------------------------------
/**@page DataObjects Data Objects 

    \em Data \em objects refer to C++ objects that do not allocate any
    resources, particularly that do not allocate any memory.  This
    makes them behave just like ordinary C structures, in that both,
    the copy constructor and assignment operator are equivalent to a
    memcopy of the memory taken by the object. Examples for data
    objects are all builtin types such as \c int or \c double or
    \e simple classes such as \c complex.
 
    We distinguish \em data \em objects from general C++ objects that
    may include some allocation of resources. (Note that for general
    C++ objects that do allocate resources, this must be respected by
    providing appropriate copy constructor and assignment operators.)
    An example for a general C++ class is class DataArray.
 
    The distinction between data and general C++ objects becomes
    relevant when using such objects in container classes such as
    DataArray or Array.  
*/
//-----------------------------------------------------------------------------
/**@page RUN On which Platforms is SoPlex running

   The current release of SoPlex was tested to compile with the following
   compilers / architectures:
   <TABLE>
   <TR><TD>Architecture</TD><TD>CPU              </TD><TD>OS                </TD><TD>Compiler             </TD></TR>
   <TR><TD>x86         </TD><TD>Pentium 4 3.2 GHz</TD><TD>Suse Linux 10.0   </TD><TD>GCC 4.0.2            </TD></TR>
   <TR><TD>x86         </TD><TD>Pentium 4 3.8 GHz</TD><TD>Suse Linux 10.0   </TD><TD>Intel 9.0.026        </TD></TR>
   <TR><TD>x86         </TD><TD>Pentium 4 3.2 GHz</TD><TD>Windows XP SP2    </TD><TD>MS Visual Studio 2003</TD></TR>
   <TR><TD>x86_64      </TD><TD>AMD Opteron 875  </TD><TD>SLES 9            </TD><TD>GCC 4.0.2            </TD></TR>
   <TR><TD>PowerPC     </TD><TD>G5 2.3           </TD><TD>MacOS-X 10.4.4    </TD><TD>GCC 4.0.0            </TD></TR>
   <TR><TD>Power4      </TD><TD>Power4 1.7 GHz   </TD><TD>AIX 5.1           </TD><TD>VisualAge 6          </TD></TR>
   <TR><TD>Alpha       </TD><TD>21264a 750 MHz   </TD><TD>OSF1 V5.1         </TD><TD>Compaq 6.5-39        </TD></TR>
   <TR><TD>Sparc       </TD><TD>Ultra 60         </TD><TD>Solaris 9         </TD><TD>Sun Studio 11        </TD></TR>
   </TABLE>

   We expect SoPlex to be compilable using any sufficiently recent compiler.
   Remember, your mileage may vary.
 */
//-----------------------------------------------------------------------------
/**@page INST Installation
 
 \section Prequisites Prequisites 
 You need the following programs to compile SoPlex:

  - C++ compiler (e.g. http://www.gnu.org/software/gcc)
  - gnu make http://www.gnu.org/software/make
  - gnu awk http://www.gnu.org/software/gawk (if you want the \c check target)
  - doxygen http://www.doxygen.org (if you want to generate the documentation)

 After receiving SoPlex you have to uncompress it with \c gzip 
 (http://www.gnu.org/software/gzip) and 
 unpack it with \c tar (http://www.gnu.org/software/tar)
 into a directory. Then change to that directory.

 \section Tested Linux/x86, Linux/AXP, Tru64, Solaris, IRIX, HP-UX
 If you are working under one those OSes you should try the following:
 
 \c gmake \c COMP=xxx \c OPT=yyy

 with \c xxx being one of \c gnu, \c intel, \c sun, \c hp, \c sgi, \c compaq
 or \c ibm
 and \c yyy one of \c dbg (if you want the debug version) or \c opt 
 (if you want the optimized version).

 This should generate a binary in the \c bin subdirectory. 

 \section Others If the previous section was not for you
 First the Makefile tries to find out the OS and the 
 architecture by using the \c uname command. 
 If your OS or architecture is missing, please update the 
 Makefile and submit the change to me.

 Then a submakefile from \c make/make.OS.ARCH.COMP.OPT is included in
 the main Makefile. You should adapt the compiler flags as needed.
 Be especially careful with the \c AR setting since some C++ compilers do
 not like using the standard \c \c ar program together with code that uses
 templates.

 If this all does not work, change to the \c src directory and type

 \c CC \c *.cpp 

 This should do the trick. Adding \c -DNDEBUG gives you a non debugging 
 version. Add flags as needed. 

 \section Testing Testing the Binary
 After you compiled the binary you should get the Netlib LP files at
 http://www.zib.de/Optimization/Software/Soplex/netlib.tar.gz and
 unpack them in the \c check directory. Then you can try

 \c gmake \c COMP=xxx \c OPT=dbg \c quick

 \c gmake \c COMP=xxx \c OPT=opt \c quick

 \c gmake \c COMP=xxx \c OPT=opt \c check

 Use the \c check target together with \c OPT=dbg only if you have really 
 a lot of time. \c quick should run in a few minutes and \c check will 
 need between less than one hour and a day, depending on your machine.

 \c quick should report no fails at all. \c check should report no fails in the
 \c LC and \c EC columns and in the \c LR and \c ER columns only with 
 the instances greenbea, greenbeb, pilot-ja or pilot87. 
 One or two fails is normal, above four is probably a problem with the 
 compiler. Look how big the error is. Try again with less optimization. 
 Our results are available at 
 http://www.zib.de/Optimization/Software/Soplex/results
 to compare with.

 \section Documentation  Generating the documentation
 If you have \c doxygen (and \c dot) installed, you just can say 

 \c gmake \c doc

 After that the documentation should be in doc/html.

 \section Installation Installation
 The binary is in the \c bin directory, the library in \c lib and all
 headers are in \c src. Feel free to install them at a suitable place.

 \section Naming Naming of the OPT Variable

  - dbg (debugging) -DNDEBUG is \b not set. Optimization is mostly off
    and debugging info is generated.

  - std (standard) -DNDEBUG is set. All optimizations may be switched on
    that do \b not \b alter the floating point behaviour or may 
    otherwise by result in wrong code.

  - opt (optimized) -DNDEBUG is set. All optimizations may be switched on,
    as long as the code seems to run correctly. The code should run
    on the relevant architectures. Best is something like -fast that
    uses the right optimizations for the architecture that is used.

  - opt-XXX (optimized for XXX) like opt-p4. This includes optimization
    for a specific processor. 

  - prf (profile) like opt, but generate profile data.

  - XXX-ld (long-double) same as XXX but with long doubles.
 */
//-----------------------------------------------------------------------------
/**@page FAQ Frequently Asked Questions

   Here are some answers that cannot be answered from the code alone:

   <ol>
   <li> Why is \<iostream\> used but \<assert.h\> and not \<cassert\> ?
     
      The reason is twofold. From the theoretical point we were not
      able to exactly find out in TC++PL in what namespace cassert 
      should load its declarations. Surely in std. But since these are
      normally functions with C linkage, this won't work.
      Then some of them, like assert, are macros which have no namespace.
      The practical point is that the compiler vendors seem to be 
      unsure as well. Most of them put everything in both namespaces std 
      and global. So there is no advantage in using \<cassert\>. Compaq 
      even left them off because it seemed unclear to them.
      So our reasoning was: If it is a C++ header we use the correct form
      without the .h and in std. If it is a C header, according to the
      standard the .h header has to be there and uses the global namespace.
      That seems acceptable to us, especially for C functions.
   </li>
      
   <li> Why is malloc/free sometimes used and not new/delete ?
      
      Because there is no realloc with new/delete. Because malloc
      is faster. And we only use it for builtin types or so called 
      \ref DataObjects "Data Objects" .
      If you do not like this decision, it is quite easy to change 
      spxalloc.h such as to use new/delete.
   </li>

   <li> Can SoPlex solve Integer Programs (IPs) ?
      
      No. You need an IP-Solver for this. Most IP-Solver use LP-Solvers
      as a subroutine and do some kind of Branch-and-Bound. For
      instance, you can use \c SCIP (Solving Constraint Integer Programs) 
      together with SoPlex to solve IPs. \c SCIP can be obtained at 
      http://scip.zib.de/ and is distributed under the ZIB academic
      license, like SoPlex.
   </li>

   <li> Is there a Windows version ?

      The code is tested to compile under Microsoft Visual Studio 2003.
      However, we do \b not provide any Makefiles or project files for VC++.
      It should also be possible to create a Windows version by compiling it
      with the GCC under <a href="http://www.cygwin.com/">Cygwin</a>.
   </li>

   <li> I want a primal and a dual simplex, where are they ?

      That is nearly easy. You can set SoPlex to use the  
      ENTERing and LEAVEing algorithm and
      COLUMN and ROW representation.

      <TABLE>
      <TR><TD>&nbsp;</TD><TD>ENTER </TD><TD>LEAVE </TD></TR>
      <TR><TD>COLUMN</TD><TD>Primal</TD><TD>Dual  </TD></TR>
      <TR><TD>ROW   </TD><TD>Dual  </TD><TD>Primal</TD></TR>
      </TABLE>

      COLUMN oriented is the "usual" representation.
      Then Entering is the Primal and Leaving is the Dual algorithm.
      In ROW oriented representation, we have in principle the
      explicit dual and then the algorithms are reversed.

      The only problem is that SoPlex is a composite simplex algorithm.
      That means it switches between entering and leaving algorithm
      as it needs. So all you can choose is which algorithm is used first,
      but then an arbitrary number of switches may occur. (Even so, often
      no switch at all happens.)
   </li>

   <li> I got an segment violation or a signal 11.

      If all of the test instances from Netlib work, but your LP gives this
      problem, mail your LP in as an gzip'ed MPS of LPF file and we will
      check. If you have this problem also with the test instances, check your
      stack space: \c ulimit \c -s will report the current size in kilobytes.
      Try a higher value. If this doesn't help, maybe your compiler is broken.
      Try compiling without optimization. 
   </li>
       
   <li> SoPlex means \em Sequential Simplex. Is there a parallel version
        available? 
 
      No. There was done research in this direction. You can find most of 
      the results in http://www.zib.de/PaperWeb/abstracts/TR-96-09 and 
      http://www.zib.de/PaperWeb/abstracts/SC-95-45 .
   </li>

   <li> Is there a wrapper class/library to use SoPlex instead of CPLEX ? 

      No. 
   </li>

   <li> Is there an interface for COIN ?

      Yes, have a look at 
      http://www.coin-or.org/documentation.html#OSI
   </li>

   <li> How can I make LP generation easier?
   
      You can use \c ZIMPL, available at http://www.zib.de/koch/zimpl/
      under the ZIB academic license. It takes a (human readable) file
      describing the linear program together with a data file as input and 
      generates LPs or MIPs in LPF- or MPS-format.

   </li>

   <li> So, how do I start working with SoPlex?

       To just read LPs from a file and solve them without modifying the LPs, 
       compile SoPlex and use it on the command line. Calling the binary
       without parameters gives you some usage information. To use SoPlex as a
       callable library, see src/example.cpp (from which the SoPleX binary is
       compiled) on how to construct an LP, solve it, and extract the
       solution. For further details, you have to consult the doxygen
       documentation of the code (or even the code itself).
   </li>

   <li> What kind of preprocessing is done?

      In the default setting SoPlex presolves the given LP in order to
      simplify the problem as far as possible by removing rows,
      columns, and bounds. Moreover, infeasibility or unboundedness may
      be detected. The simplification cannot be undone, but given a
      primal/dual solution for the simplified LP, the simplifier can
      reconstruct the primal/dual solution and basis of the
      unsimplified LP.
      
      Handled are:
      - empty rows / columns
      - unconstrained rows
      - row singletons
      - forcing rows
      - zero objective column singletons
      - (implied) free column singletons
      - doubleton equations combined with a column singleton
      - (implicitly) fixed columns
      - redundant lhs / rhs
      - redundant variable bounds
      - variables that are free in one direction
      - (weakly) dominated columns
      - duplicate rows / columns

      The preprocessing can be switched off by using the option \c -s0.
   </li>

   </ol>
*/           
//-----------------------------------------------------------------------------
/**@page FLAGS How to use the flags in the example program

   Here are some tips on which flags to use with the example program:

   If you have more constraints (rows) than variables (cols) it is 
   a good idea to try the \c -r flag to choose a row-wise representation
   of the basis.

   Setting the different epsilons to a smaller value like 1e-18 or 1e-20
   (using \c -zz or \c -zu) might improve the quality of the
   solution, but may also slow down the program. Setting the epsilons to
   bigger values may speed up the algorithm, but values greater than 1e-12 
   are definitely a bad idea.

   Setting \c -d to smaller values like 1e-7 or 1e-8 will improve the
   quality of the solution, but it will take longer. Values smaller
   then 1e-9 are not recommended. The \c -d value should be
   substantial bigger then the \c -z values.

   If the default settings are too slow, using \c -e eventually together 
   with \c -p1 might improve the running time.  
*/
//-----------------------------------------------------------------------------
/**@page PROG Programming with SoPlex 
   
   The SoPlex class library comprises classes that may be categorized into
   three different types:

   - Elementary classes are provided for general purpose use in
     projects way beyond the scope of numerical software or linear
     programming.
   - Linear algebra classes provide basic data types for (sparse)
     linear algebra computations. However, their functionality is
     restricted to simple operations such as addition and scaling.
     For complex tasks, such as solving linear systems of equations,
     algorithmic classes are provided instead.
   - Algorithmic classes serve for implementing maybe a variety of
     algorithms for solving numerical (sub-)problems.

   The following sections are dedicated to users who want to
   provide own pricers, ratio test, start basis generation codes or
   LP simplifiers to use with SoPlex or who want to derive own
   implementations (e.g. parallel versions) using SoPlex.

   @section Virtualizing the Representation
   The primal Simplex on the columnwise representation is
   structurally equivalent to the dual Simplex on the rowwise
   representation and vice versa (see above). Hence, it is
   desirable to treat both cases in a very similar manner. This
   is supported by the programmer's interface of SoPlex which
   provides access methods for all internal data in two ways: one
   is relative to the "physical" representation of the LP in
   rows and columns, while the other is relative to the chosen
   basis representation. If e.g. a soplex::SPxPricer is
   written using the second type of methods only (which will
   generally be the case), the same code can be used for running
   SoPlex's simplex algorithm for both representations. 
   We will now give two examples for this
   abstraction from the chosen representation.

   Methods \c vector() will return a column or a row vector,
   corresponding to the chosen basis representation. 
   The other "vectors" will be referred to as \em covectors:
     
   <TABLE>
   <TR><TD>&nbsp;  </TD><TD>ROW      </TD><TD>COLUMN   </TD></TR>
   <TR><TD>vector  </TD><TD>rowVector</TD><TD>colVector</TD></TR>
   <TR><TD>coVector</TD><TD>colVector</TD><TD>rowVector</TD></TR>
   </TABLE>
    
   Weather the soplex::SPxBasis::Desc::Status of a variable indicates that the
   corresponding vector is in the basis matrix or not also depends on the
   chosen representation. Hence, methods \c isBasic() are provided to get the
   correct answer for both representations.  
   
   @section Simplex Vectors and Bounds
   The Simplex algorithms keeps three vectors which are associated to each basis.
   Two of them are required for the pricing, while the third one is needed for
   detecting feasibility of the basis. For all three vectors, bounds are
   defined. The Simplex algorithm changes the basis until all three vectors
   satisfy their bounds, which means that the optimal solution has been found.
    
   With each update of the basis, also the three vectors need to be
   updated. This is best supported by the use of UpdateVectors.
    
   @subsection Variables
   The Simplex algorithm works with two types of variables, primals and
   duals.  The primal variables are associated with each column of
   an LP, whereas the dual variables are associated with each row.
   However, for each row a slack variable must be added to the set of
   primals (to represent inequalities), and a reduced cost variable must be
   added for each column (to represent upper or lower bounds). Note, that
   mathematically, one dual variable for each bound (upper and lower) should
   be added. However, this variable would always yield the same value and
   can, hence, be implemented as one.
    
   To summarize, we have a primal variable for each LP column and row
   (i.e., its slack) as well as a dual variable for each LP row and column
   (i.e., its bounds). However, not all these values need to be stored and
   computed, since the structure of the Simplex algorithms allow to
   keep them implicitly.
      
   If the SPxBasis's Status of a row or column is one of \c P_ON_LOWER,
   \c P_ON_UPPER, \c P_FIXED or \c P_FREE, the value of the corresponding
   primal variable is the lower, upper or both bound(s) or 0, respectively.
   The corresponding dual variable needs to be computed. Equivalently, for
   a Status of \c D_FREE, \c D_ON_UPPER, \c D_ON_LOWER, \c D_ON_BOTH or
   \c D_UNDEFINED, the corresponding dual variable is 0, whereas the primal 
   one needs to be computed.

   The following vectors are declared for holding the values to be computed:
   \c primRhs, \c primVec (with dimension \c nCols()) for the primal
   variables, and \c dualRhs, \c dualVec (with dimension \c nRows()) for the 
   dual variables. The additional variable \c addvec (with dimension \c coDim() 
   depends on the representation.

   @subsection Bounds 
   Primal and dual variables are bounded (including \f$\pm\infty\f$ as
   bounds).  If all primal variables are within their bounds, the
   Simplex basis is said to be primal feasible. Analogously, if all
   dual variables are within their bounds, its is called dual
   feasible.  If a basis is both, primal and dual feasible, the
   optimal solution has been found.

   In the dual Simplex, the basis is maintained dual feasible, while
   primal feasibility is improved via basis updates. However, for
   numerical reasons dual feasibility must be relaxed from time to time.
   Equivalently, primal feasibility will be relaxed to
   retain numerical stability in the primal Simplex algorithm.

   Relaxation of (dual or primal) feasibility is achieved by
   relaxing the bounds of primal or dual variables. However, for each
   type of Simplex only the corresponding bounds need to be
   relaxed. Hence, we define only one vector of upper and lower bound
   for each row and column and initialize it with primal or dual
   bound, depending on the Simplex type (see \c theURbound,
   \c theLRbound, \c theUCbound, \c theLCbound). 
*/
//-----------------------------------------------------------------------------


