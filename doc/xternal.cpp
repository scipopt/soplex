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
/**@page Platforms Where is SoPlex running?
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
   <TR><TD>HP    </TD><TD>A.03.27  </TD><TD>HP-UX 11.00       </TD></TR>
   <TR><TD></TD><TD></TD><TD></TD></TR>
   <TR><TD></TD><TD></TD><TD></TD></TR>
   </TABLE>
 */
//-----------------------------------------------------------------------------
/**@page FAQ Frequently Asked Questions

   Here are some answers that can not be answered from the code alone:

   <ol>
   <li> Why is <iostream> used but <assert.h> and not <cassert>
     
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
      
   <li> Why is malloc/free sometimes used and not new/delete.
      
      Because there is no realloc with new/delete. Because malloc
      is faster. And we only use it for builtin types.
      If you do not like this descision, it's easy to change spxalloc.h
      to use new/delete.

   <li> Can SoPlex solve Integer Programs (IP's).
      
      No. You need an IP-Solver for this. Most IP-Solver use LP-Solvers
      as a subroutine and do some kind of Branch-and-Bound.

   <li> Is there a Windows version ?

      The code is tested to compile under some version of Visual C.
      We do \b not provide any make or project files for VC.

   <li> What is the academic license ?

      Essentially, you can do what you want, except: Give SoPlex to 
      anybody else. Do anything with SoPlex that makes money (this 
      includes saving your money on any commercial activity).
      See the license for details http://pengpuffblub.html

   <li> What is the commercial license ?
  
      If you want to make with SoPlex you need one. We give the license,
      you pay for it. The amount is completely negotiable depending on
      what you want to do with SoPlex, which rights you want and what
      you are willing to tell us.

   <li> I want I primal and a dual simplex, where are they ?

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
   </ol>
*/           
//-----------------------------------------------------------------------------
/**@page Programming Programming with SoPlex 
   
   The following sections are dedicated to users, that want to
   provide own pricers, ratio test, start basis generation codes or
   LP simplifiers to use with #SoPlex or that want to derive own
   implementations (e.g. parallel versions) using SoPlex.

   @section Virtualizing the Representation
   The primal Simplex on the columnwise representation is
   structurely equivalent to the dual Simplex on the rowwise
   representation and vice versa (see above). Hence, it is
   desirable to treat both cases in a very similar manner. This
   is supported by the programmers interface of #Soplex, that
   provides access methods for all internal data in two ways: one
   is relative to the "physical" representation of the LP in
   rows and columns, while the other is relative to the chosen
   basis #Representation. If e.g. a #SPxPricer is
   written using the second type of methods only (which will
   generally be the case), the same code can be used for running
   SoPlex's the Simplex algorithm for both
   #Representation%s. We will now give two examples for this
   virtualization from the chosen representation.

   Methods #vector will return a column or a row vector,
   corresponding to the chosen basis representation. 
   The other "vectors" will be referred to as \em covectors:
     
   <TABLE>
   <TR><TD>&nbsp;  </TD><TD>ROW      </TD><TD>COLUMN   </TD></TR>
   <TR><TD>vector  </TD><TD>rowVector</TD><TD>colVector</TD></TR>
   <TR><TD>coVector</TD><TD>colVector</TD><TD>rowVector</TD></TR>
   </TABLE>
    
   Weather the #SPxBasis::Desc::Status of a variable indicates that the
   corresponding #vector is in the basis matrix or not also depends on the
   chosen representation. Hence, methods #isBasic() are provided to get the
   correct answer for both representations.  
   
   @section Simplex Vectors and Bounds
   The Simplex algorithms keeps three vectors, that are defined to a basis.
   Two of them are required for the pricing, while the other is needed for
   detecting feasibility of the basis. For all three vectors, bounds are
   defined. The Simplex alogrithm changes the basis until all three vectors
   satisfy their bounds. In this case the optimal solution is found.
    
   Whith each update of the basis, also the three vectors need to be
   updated. This is best supported by the use of #UpdateVector%s.
    
   @subsection Variables
   The Simplex algorithm works with two types of variables, primals and
   duals.  The primal variables are the ones associated to each column of
   an LP, whereas the dual variables are the ones associated to each row.
   However, to each row a slack variable must be added to the set of
   primals (to represent inequalities), and a reduced cost variable must be
   added for each column (to represent upper or lower bounds). Note, that
   mathematically, on dual variable for each bound (upper and lower) should
   be added. However, this variable would always yield the same value and
   can, hence, be implemented as one.
    
   To summarize, we have a primal variable for each LP column and row
   (i.e. its slack) as well as a dual variable for each LP row and column
   (i.e. its bounds). However, not all these values need to be stored and
   computed, since the structure of the Simplex algorithms allow to
   implicitely keep them.
      
   If the #SPxBasis's #Status of a row or column is one of #P_ON_LOWER,
   #P_ON_UPPER, #P_FIXED or #P_FREE the value of the corresponding
   primal variable is the lower, upper or both bound(s) or 0, respectively.
   The corresponding dual variable needs to be computed. Equivalently, for
   a #Status of #D_FREE, #D_ON_UPPER, #D_ON_LOWER, #D_ON_BOTH or
   #D_UNDEFINED the corresponding dual variable is 0, whereas the primal 
   one needs to be computed.

   We declare the following vectors for holding the values to be computed.
   Primal variables (dimension #nCols()): #primRhs, #primVec.
   Dual variables (dimension #nRows()): #dualRhs, #dualVec.
   Additional variables depending on representation (dimension coDim()):
   #addvec.

   @subsection Bounds 
   Dual and primal variables are bounded (including \f$\pm\infty\f$ as
   bounds).  If all primal variables are within their bounds, the
   Simplex basis is said to be primal feasible. Analogously, if all
   dual variables are within their bounds, its is called dual
   feasible.  If a basis is both, primal and dual feasible, the
   optimal solution is been found.

   In the dual Simplex, the basis is maintained to be dual, while
   primal feasibility is improved via basis updates. However, for
   numerical reasons dual feasibility must from time to time be
   relaxed.  Equivalently, primal feasibility will be relaxed to
   retain numerical stability in the primal Simplex algorithm.

   Relaxation of (dual or primal) feasibility is acchieved by
   enlarging the bounds to primal or dual variables. However, for each
   type of Simplex only the corresponding bounds need to be
   enlarged. Hence, we define only one vector of upper and lower bound
   for each row and column and initialize it with primal or dual
   bound, depending on the Simplex type.  (see #theURbound,
   #theLRbound, #theUCbound, #theLCbound.) 
*/
//-----------------------------------------------------------------------------
/**@page License ZIB ACADEMIC LICENSE

This license for ZIB software is designed to guarantee freedom to share and 
change software for academic use, but restricting commercial firms
from exploiting your knowhow for their benefit. The precise terms and 
conditions for using, copying, distribution, and modification follow. 

Terms and Conditions for Using, Copying, Distribution, and Modification

The "Program" below refers to source, object and executable code, and a 
"work based on the Program" means either the Program or any
derivative work under copyright law: that is a work containing the Program 
or a portion of it, either verbatim or with modifications and/or translated
into another language. Each licensee is addressed as "you". 

<ol>
<li> This license applies to you only if you are a member of a noncommercial 
     and academic institution, e.g., a university. The license expires as
     soon as you are no longer a member of this institution. 
<li> Every publication and presentation for which work based on the Program 
     or its output has been used must contain an appropriate citation
     and aknowledgement of the author(s) of the Program. 
<li> You may copy and distribute the Program or work based on the Program 
     in source, object, or executable form provided that you also meet
     all of the following conditions: 

 <ol>
 <li> You must cause any work that you distribute or publish, that in whole 
      or in part contains or is derived from the Program or any
      part thereof, to be licensed as a whole at no charge under the terms of
      this License. You must accompany it with this unmodified license text. 

      These requirements apply to the Program or work based on the Program 
      as a whole. If identifiable sections of that work are
      not derived from the Program, and can be reasonably considered 
      independent and separate works in themselves, this License
      does not apply to those sections when you distribute them as 
      separate works. But when you distribute the same sections as
      part of a whole which is a work based on the Program, the distribution 
      of the whole must be on the terms of this License,
      whose permissions for other licensees extend to the entire whole and, 
      thus, to each and every part regardless of who wrote it. 

 <li> You must cause the modified files to carry prominent notices stating 
      that you changed the files and the date of any change. 
 <li> You must keep track of access to the Program (e.g., similar to the 
      registration procedure at ZIB). 
 <li> You must accompany it with the complete corresponding machine-readable 
      source code. 

      The source code for a work means the preferred form of the work for 
      making modifications to it. For an executable work,
      complete source code means all the source code for all modules it 
      contains, plus any associated interface definition files, plus
      the scripts used to control compilation and installation of the 
      executable. However, as a special exception, the source code
      distributed need not include anything that is normally distributed 
      (in either source or binary form) with the major components
      (compiler, kernel, and so on) of the operating system on which the 
      executable runs, unless that component itself accompanies 
      the executable. 
 </ol>

<li> You may not copy, modify, sublicense, or distribute the Program except as 
     expressly provided under this License. Any attempt otherwise to
     copy, modify, sublicense, or distribute the Program is void and will 
     automatically terminate your rights under this License. However, parties
     who have received copies or rights from you under this License will 
     not have their licenses terminated so long as such parties remain in full
     compliance. 
<li> You are not required to accept this License, since you have not signed 
     it. However, nothing else grants you permission to use, modify, or
     distribute the Program or its derivative works. These actions are 
     prohibited by law if you do not accept this License. Therefore, by using,
     modifying or distributing the Program (or any work based on the Program), 
     you indicate your acceptance of this License to do so and all its
     terms and conditions for copying, distributing or modifying the Program 
     or works based on it. 
<li> Each time you redistribute the Program (or any work based on the Program),
     the recipient automatically receives a license from the original
     licensor to copy, distribute or modify the Program subject to these terms 
     and conditions. You may not impose any further restrictions on the
     recipient's exercise of the rights granted herein. You are not 
     responsible for enforcing compliance by third parties to this License. 
<li> If, as a consequence of a court judgment or allegation of patent 
     infringement or for any other reason (not limited to patent issues), 
     conditions are imposed on you (whether by court order, agreement, or 
     otherwise) that contradict the conditions of this License, they do not 
     excuse you from the conditions of this License. 
<li> If you wish to incorporate parts of the Program into other programs whose 
     distribution conditions are different, write to ZIB to ask for 
     permission. 

<li> Because the program is licensed free of charge, there is no warranty for 
     the program to the extent permitted by applicable law. The
     copyright holders provide the program "as is" without warranty of any 
     kind, either expressed or implied, including, but not limited to, the
     implied warranties of merchantability and fitness for a particular 
     purpose. The entire risk as to the quality and performance of the 
     program is with you. Should the program prove defective, you assume the 
     cost of all necessary servicing, repair, or correction. 
<li> In no event will any copyright holder, or any other party who may modify 
     and/or redistribute the program as permitted above, be liable to
     you for damages, including any general, special, incidental or 
     consequential damages arising out of the use or inability to use the 
     program (including but not limited to loss of data or data being rendered 
     inaccurate or losses sustained by you or third parties or a failure of the
     program to operate with any other programs), even if such holder or other 
     party has been advised of the possibility of such damages. 
</ol>
*/




