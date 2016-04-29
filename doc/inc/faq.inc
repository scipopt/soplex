<!-- this file is used to generate the local doxygen documentation -->
<!-- using make doc from within scip or soplex -->

<style>
.reveal:hover {
    text-shadow: 1px 1px 1px #777;
}
.reveal {
    color:#06c;
    padding-top: 1em;
    cursor: pointer;
    margin: 5px 0;
}
.answer {
    background-color: #fff;
    padding-left:   1em;
}
</style>

<h4>Frequently asked questions about SoPlex</h4><ol>  <li>
    <div id="installationproblems" class="reveal">
      I have installation problems. What can I do?    </div>
    <div id="installationproblems_ans" class="answer">
       <p>
Please have a look at the <code>INSTALL</code> file. If you are using SoPlex
version 2.0.0 or higher you should try to recompile with <code>LEGACY=true</code>.
This will restore the previous command line interface and disable some new
features like iterative refinement with rational arithmetic.
</p>
    </div>
  </li>
    <li>
    <div id="howtocontribute" class="reveal">
      How can I contribute?    </div>
    <div id="howtocontribute_ans" class="answer">
       <p>
Please report
<a href="http://soplex.zib.de/#bugs">bugs and problems</a>.<br/>

In fact, a couple of improvements of SoPlex are
possible. If you are interested in helping please contact
us.
</p>
    </div>
  </li>
    <li>
    <div id="howtostart" class="reveal">
      So, how do I start working with SoPlex?    </div>
    <div id="howtostart_ans" class="answer">
       <p>
To just read LPs from a file and solve them without modifying the LPs,
compile SoPlex and use it on the command line. Calling the binary
without parameters gives you some usage information. To use SoPlex as a
callable library, see src/example.cpp or src/soplexmain.cpp (from which
the SoPlex binary is compiled) on how to construct an LP, solve it, and
extract the solution. For further details, you have to consult the doxygen
documentation of the code (or even the code itself).
Another good starting point might also be
the SoPlex interface of <a href="http://scip.zib.de">SCIP</a>.
</p>
    </div>
  </li>
    <li>
    <div id="usesoplexinsidemipsolvers" class="reveal">
      Can I use SoPlex with Branch-and-cut or MIP solvers?    </div>
    <div id="usesoplexinsidemipsolvers_ans" class="answer">
       <p>
Yes, but the corresponding software needs an interface.

Currently, the following frameworks have an
interface to SoPlex.
</p>
<ul>
<li>
<a  href="http://scip.zib.de">SCIP</a> - Solving
Constraint Integer Programs.
</li>
<li>
<a href="http://www.informatik.uni-koeln.de/abacus/">ABACUS</a>
has an older interface to 1.2.1 only.
</li>
<li>
There is also an interface to <a
href="https://projects.coin-or.org/Osi">COIN
OSI</a>.
</li>
</ul>
    </div>
  </li>
    <li>
    <div id="mallocvsnew" class="reveal">
      Why is malloc/free sometimes used and not new/delete?    </div>
    <div id="mallocvsnew_ans" class="answer">
       <p>
Because there is no realloc with new/delete. Because malloc
is faster. And we only use it for builtin types or so called
"Data Objects" .
If you do not like this decision, it is quite easy to change
spxalloc.h such as to use new/delete.
</p>
    </div>
  </li>
    <li>
    <div id="solvingips" class="reveal">
      Can SoPlex solve Integer Programs (IPs)?    </div>
    <div id="solvingips_ans" class="answer">
       <p>
No. You need an IP-Solver for this. Most IP-Solver use LP-Solvers
as a subroutine and do some kind of Branch-and-Bound. For
instance, you can use  SCIP (Solving Constraint Integer Programs)
together with SoPlex to solve IPs.  SCIP can be obtained at
<a href="http://scip.zib.de/">here</a> and is
distributed under the ZIB academic license, like SoPlex.
</p>

    </div>
  </li>
    <li>
    <div id="windoof" class="reveal">
      Is there a Windows version?    </div>
    <div id="windoof_ans" class="answer">
       <p>
The code is tested to compile under Microsoft Visual Studio 2008 and 2010.
It should also be possible to create a Windows version by compiling it
with the GCC under <a href="http://www.cygwin.com/">Cygwin</a> or minGW
</p>
    </div>
  </li>
    <li>
    <div id="primaldualsimplex" class="reveal">
      I want a primal and a dual simplex, where are they?    </div>
    <div id="primaldualsimplex_ans" class="answer">
       <p>
SoPlex implements a composite simplex, i.e., it switches between primal and dual
simplex and the user can only specify the starting algorithm via the integer
parameter "int:algorithm = 0 (primal) / 1 (dual)".
</p>
<p>
When using the low-level class SPxSolver, you can specify to use the
ENTERing and LEAVEing algorithm and COLUMN and ROW representation.
</p>
<center>
<table border="1" rules="all">
<tr><td>&nbsp;</td><td>ENTER </td><td>LEAVE </td></tr>
<tr><td>COLUMN</td><td>Primal</td><td>Dual  </td></tr>
<tr><td>ROW   </td><td>Dual  </td><td>Primal</td></tr>
</table>
</center>
<p>
COLUMN oriented is the "usual" representation.
Then Entering is the Primal and Leaving is the Dual algorithm.
In ROW oriented representation, we have in principle the
explicit dual and then the algorithms are reversed.
</p>
    </div>
  </li>
    <li>
    <div id="segmentationfault" class="reveal">
      I got a segment violation or a signal 11.    </div>
    <div id="segmentationfault_ans" class="answer">
       <p>
If all of the test instances from Netlib work, but your LP gives this
problem, mail your LP in as an gzip'ed MPS of LP file and we will
check. If you have this problem also with the test instances, check your
stack space:  ulimit  -s will report the current size in kilobytes.
Try a higher value. If this doesn't help, maybe your compiler is broken.
Try compiling without optimization.
</p>
    </div>
  </li>
    <li>
    <div id="numericaltroubles" class="reveal">
      I got messages or exceptions stating &quot;This shall not be.&quot; or similar.    </div>
    <div id="numericaltroubles_ans" class="answer">
       <p>
Such messages are an indication of numerical trouble, which may happen
from time to time especially on LPs with weird numbers. If the problem
occurs with the binary, try using another scaler (call soplex without
parameters to see how to specify this). If the problem occurs while using
the callable library of SoPlex, you can try to add a scaler and a
preprocessor yourself (see the SoPlex class) unless you solve the same LP
repeatedly with small changes (as in a Branch-and-cut process, for
instance). If you send us an LP or MPS file with the offending linear
program together with a description of your settings, we may have a
closer look at it to see if we can do something about it.
</p>

    </div>
  </li>
    <li>
    <div id="basisfiles" class="reveal">
      Can SoPlex handle basis files?    </div>
    <div id="basisfiles_ans" class="answer">
       <p>
Yes, the binary version allows to read and write basis files using the
options <code>--readbas=</code> and <code>writebas=</code>, respectively. The library
version allows to read or write a basis file from within the code as well.
This can be useful for debugging. If you cannot reconstruct an error occurring
in the library code with the binary version, try to save a basis and try again.
Note that the LP needs to be in MPS format in order to ensure compatability of the
basis and the LP (the LP format cannot store ranged rows).
</p>
    </div>
  </li>
    <li>
    <div id="parallelversion" class="reveal">
      SoPlex means &quot;Sequential Simplex&quot;. Is there a parallel version
available?    </div>
    <div id="parallelversion_ans" class="answer">
       <p>
No. Some research was done in this direction. You can find most of
the results in
<a href="http://www.zib.de/PaperWeb/abstracts/TR-96-09">
http://www.zib.de/PaperWeb/abstracts/TR-96-09</a> and
<a href="http://www.zib.de/PaperWeb/abstracts/SC-95-45">
http://www.zib.de/PaperWeb/abstracts/SC-95-45</a>.
</p>
    </div>
  </li>
    <li>
    <div id="soplexwrapperclass" class="reveal">
      Is there a wrapper class/library to use SoPlex instead of CPLEX ?    </div>
    <div id="soplexwrapperclass_ans" class="answer">
       <p>
No.
</p>
    </div>
  </li>
    <li>
    <div id="lpgeneration" class="reveal">
      How can I make LP generation easier?    </div>
    <div id="lpgeneration_ans" class="answer">
       <p>
You can use <a href="http://zimpl.zib.de/">ZIMPL</a>, available at
under the ZIB academic license. It takes a (human readable) file
describing the linear program together with a data file as input and
generates LPs or MIPs in LP- or MPS-format.
</p>
    </div>
  </li>
    <li>
    <div id="preprocessing" class="reveal">
      What kind of preprocessing is done?    </div>
    <div id="preprocessing_ans" class="answer">
       <p>
In the default setting SoPlex presolves the given LP in order to
simplify the problem as far as possible by removing rows,
columns, and bounds. Moreover, infeasibility or unboundedness may
be detected. The simplification cannot be undone, but given a
primal/dual solution for the simplified LP, the simplifier can
reconstruct the primal/dual solution and basis of the
unsimplified LP.
</p>
<p>
Handled are:
</p>
<ul>
<li> empty rows / columns                                  </li>
<li> unconstrained rows                                    </li>
<li> row singletons                                        </li>
<li> forcing rows                                          </li>
<li> zero objective column singletons                      </li>
<li> (implied) free column singletons                      </li>
<li> doubleton equations combined with a column singleton  </li>
<li> (implicitly) fixed columns                            </li>
<li> redundant lhs / rhs                                   </li>
<li> redundant variable bounds                             </li>
<li> variables that are free in one direction              </li>
<li> (weakly) dominated columns                            </li>
<li> duplicate rows / columns                              </li>
</ul>

<p>
The preprocessing can be switched off by using the option  <code>-s0</code>.
</p>    </div>
  </li>
  </ol><br/>
<script>
$(".answer").hide();
$(".reveal").click(function() {
    var questionId = "#"+this.id+"_ans";
    $(questionId).toggle("fast");
});
</script>