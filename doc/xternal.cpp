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

  This software has been implemented as a part of my Ph.D. thesis 
  "Paralleler und Objektorientierter Simplex-Algorithmus"
  http://www.zib.de/PaperWeb/abstracts/TR-96-09

  Note, that you are allowed to retreive SoPlex only for research purpose as a
  member of a {\em noncommercial and academic} institution.
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
/**@defgroup Algebra Linea Algebra Classes
   @brief Basic data types for linear algebra computations.
   
   Linear algebra classes provide basic data types for (sparse)
   linear algebra computations. However, their functionality is
   restricted to simple operations such as addition and scaling.
   For complex tasks, such as solving linear systems of equations,
   algorithmic classes are provided instead.
*/
//-----------------------------------------------------------------------------
/**@defgroup Algorithmic Algorithmic Classes
   @brief Implementation of numerical algorithms.   
   
   Algorithmic classes serve for implementing a variaty of
   algorithms for solving numerical (sub-)problems.
*/
//-----------------------------------------------------------------------------

/*
{\em Data Objects} refer to C++ objects that do not allocate any resources,
particularly that do not allocate any memory.  This makes them behave just
like ordinary C structures, in that both, the copy constructor and
assignment operator are equivalent to a memcopy of the memory taken by the
object. Examples for data objects are all builtin types such as {\tt\strut int} or
{\tt\strut double} or ``simple'' classes such as {\tt\strut complex}.

We distinguish {\em data objects} from general C++ objects that may include
some allocation of resources. (Note, that for general C++ objects that do
allocate resources, this must be respected by providing apropriate copy
constructor and assignment operators.) An example for a general C++ class is
class DataArray ($\rightarrow$1.4, {\em page \pageref{cxx.1.4}}).

The distingtion between data and general C++ objects becomes relevant when
using such objects in container classes such as DataArray ($\rightarrow$1.4, {\em page \pageref{cxx.1.4}}) or
Array ($\rightarrow$1.1, {\em page \pageref{cxx.1.1}}).
*/
