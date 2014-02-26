//-----------------------------------------------------------------------------
/**@mainpage Overview
   @version  2.0.0

   @section Main The Sequential object-oriented simplex class library.

   This software has been implemented as a part of Roland Wunderling's 
   Ph.D. thesis "Paralleler und Objektorientierter Simplex-Algorithmus"
   which can be found at http://www.zib.de/PaperWeb/abstracts/TR-96-09
   (in German).

   SoPlex is part of the SCIP Optimization Suite.  A tutorial article for
   getting started with the SCIP Optimization Suite is available as ZIB
   technical report 12-27 <a href="http://scip.zib.de/doc/ZR-12-27.pdf">here</a>.

   SoPlex is implemented in C++. The code should be compliant with the 
   current ANSI standard. RTTI and STL (other then iostream) are not used. 
   Everything is in one namespace \em soplex.

   - \ref INST     "Installation"

   - \ref FAQ      "Frequently asked questions"

   - \ref SHELL    "How to use the SoPlex command line"

   - \ref LIB      "How to use SoPlex as a callable library"

   - \ref PROG     "Programming with SoPlex"

   @author   Roland Wunderling
   @author   Tobias Achterberg
   @author   Timo Berthold
   @author   Andreas Bley
   @author   Dennis Elbrächter
   @author   Ambros Gleixner
   @author   Wei Huang
   @author   Benjamin Hiller
   @author   Thorsten Koch
   @author   Matthias Miltenberger
   @author   Sebastian Orlowski
   @author   Marc Pfetsch
   @author   Eva Ramlow
   @author   Daniel Steffy
   @author   Andreas Tuchscherer

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
/**@page INST Installation

   See the INSTALL file in the distribution.
 */
//-----------------------------------------------------------------------------
/**@page FAQ Frequently Asked Questions
 * \htmlinclude faq.inc
 */           
//-----------------------------------------------------------------------------
/**@page SHELL How to use the SoPlex command line

   Running the command line version of SoPlex without any arguments displays
   a list of options.  You can write a parameter file with default parameters
   by using option --saveset=FILENAME.set.  After changing parameter values
   in this file you can use it by with --loadset=FILENAME.set.
   The most frequently used parameters have abbreviations as explained in the
   initial help.

   The old command line interface of version 1.x is available when you
   compile with LEGACY=true (using the provided Makefile) or compile with the
   preprocessor define SOPLEX_LEGACY.
*/
//-----------------------------------------------------------------------------
/**@page LIB How to use SoPlex as a callable library

   The main interface is given by the class \ref soplex::SoPlex, which handles the
   construction and modification of an LP, the solving process, allows to
   access and change parameters, and retrive solution information.

   With version 2.0, the SoPlex class has been updated significantly compared
   to the 1.x version.  Although this is deprecated, it is still possible to use
   the old interface class by compiling with LEGACY=true (using the provided
   Makefile) or defining the preprocessor flag SOPLEX_LEGACY.
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

   @section Representation Virtualizing the Representation
   The primal Simplex on the columnwise representation is
   structurally equivalent to the dual Simplex on the rowwise
   representation and vice versa (see below). Hence, it is
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
    
   Whether the soplex::SPxBasis::Desc::Status of a variable indicates that the
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
   updated. This is best supported by the use of \c UpdateVectors.
    
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
   dual variables. The additional variable \c addvec (with dimension \c coDim())
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
/**@page IR Iterative Refinement

   Since version 1.7, SoPlex provides the new feature \em iterative \em refinement that
   allows for computing extended-precision solutions beyond the limits of
   standard floating-point arithmetic.  It may be particularly helpful for
   numerically troublesome LPs and applications that require solutions within
   tight feasibility tolerances.

   For a detailed explanation of the iterative refinement algorithm see the ZIB
   technical report 12-19 available online at

   http://opus4.kobv.de/opus4-zib/frontdoor/index/index/docId/1545

   Vaguely speaking, at each round of the iterative refinement process, SoPlex
   computes the primal and dual violation of the current solution in exact
   arithmetic and uses this to set up an LP with modified bounds, sides, and
   objective function.  The solution to this LP is then used to correct the
   solution and reduce primal and dual violation.  This is repeated until the
   requested feasibility and optimality tolerance is reached.

   For this, SoPlex needs to be compiled with GMP support, see the INSTALL file
   in the distribution.

 */
//-----------------------------------------------------------------------------


