<!DOCTYPE HTML PUBLIC "-//W3C//DTD HTML 4.0 Transitional//EN"
"http://www.w3.org/TR/REC-html40/loose.dtd">
<!-- $Id: soplex.php,v 1.8 2002/04/10 06:48:56 bzfkocht Exp $ -->
<?php
$counter = "cgi/counter.dat";  # the data storage file

if (file_exists($counter)) {
	$fp = fopen($counter, 'r');
	$buffer = fread($fp, filesize($counter));
	fclose($fp);
	$visitorCount = intval($buffer);
}
else
{
	$visitorCount = 0;
}
$fp = fopen($counter, 'w');
fwrite($fp, ++$visitorCount);
fclose($fp);
?>
<html lang="en">
<head>
<meta name="author" content="Thorsten Koch">
<meta name="description" content="SoPlex Homepage">
<meta name="keywords" content="SoPlex, Simplex, Solver, LP, Linear Programming, Mathematical Programming, Optimization">
<meta name="date" content="$Date: 2002/04/10 06:48:56 $">
<meta name="robots" content="index|follow">
<meta http-equiv="content-type" content="text/html;charset=iso-8859-1">
<title>SoPlex</title> 
</head>
<body bgcolor="#FFFFFF" text="#000000" 
      link="#0000FF" alink="#FF0000" vlink="#800080">
<p>
<!-- <img src="http://www.zib.de/global/images/zib_logo1.gif" -->

<a href="http://www.zib.de/index.en.html">
   <img src="images/zib_logo1.gif" 
   alt="Goto ZIB" border=0 align=right width=86 height=118>
</a>
</p>  
<h2>Konrad-Zuse-Zentrum für Informationstechnik Berlin</h2>
<h3>Division Scientific Computing</h3>
<h3>Department Optimization</h3>
<hr noshade>
<h1>SoPlex</h1>
<p>
<big>The <b>S</b>equential <b>o</b>bject-oriented sim<b>plex</b> 
class library</big>
</p>
<hr noshade>
<h3>News</h3>
<table>
<tr>
 <td><img src="images/newest.gif" alt="New" border=0></td>
 <td>09. Apr 2002</td>
 <td>Version 1.2.1 is released. <a href="notes-121.txt">Release Notes</a></td>
</tr>
<tr>
 <td>&nbsp;</td>
 <td>14. Jan 2002</td>
 <td>Version 1.2.0 is released. <a href="notes-120.txt">Release Notes</a></td>
</tr>
</table>
<hr noshade>
<h3>What is SoPlex?</h3>
<p>
SoPlex is an implementation of the revised simplex algorithm.
It features primal and dual solving routines for linear programs and 
is implemented as a C++ class library that can be used with other 
programs. An example program to solve standalone linear programs 
given in 
<a href="http://www6.software.ibm.com/sos/features/featur11.htm">MPS</a>
or LP-Format files is also included.</p>
<p>
SoPlex has been implemented as a part of Roland Wunderling's 
Ph.D. thesis 
<a href="http://www.zib.de/PaperWeb/abstracts/TR-96-09">
<em>Paralleler und Objektorientierter Simplex-Algorithmus</em></a>
(in German).
</p>
<h3>Where does is run?</h3>
<p>
SoPlex is now completely implemented in C++. 
The code should be compliant with the current ANSI standard. 
Exceptions, RTTI, and STL (other then iostream) 
are not used. With a decent modern C++ compiler you should have a chance.
<br>
We have tested SoPlex with compilers from 
<a href="http://www.gnu.org/software/gcc">GNU</a>,
<a href="http://www.compaq.com/products/software/compilers/candcxx.html">Compaq</a>, 
<a 
href="http://support.intel.com/support/performancetools/c/v5/linux/index.htm">Intel</a>,
<a href="http://www.sun.com/forte/cplusplus">SUN</a>, 
<a href="http://h21007.www2.hp.com/dspp/tech/tech_TechSoftwareDetailPage_IDX/1,1703,1740,00.html">HP</a>,
<a href="http://www.sgi.com/developers/devtools/languages/mipspro.html">SGI</a>,
<a href="http://www-4.ibm.com/software/ad/vacpp">IBM</a>, 
and even M$.
</p>
<h3>What are the license terms?</h3>
<p>
SoPlex is distributed under the 
<a href="academic.txt">ZIB Academic License</a>.
You are allowed to retrieve SoPlex only for research purpose 
as a member of a <em>non-commercial</em> and <em>academic</em> institution.
</p>
<p>
<b>Any publication for which SoPlex is used must include an
acknowledgment and a reference to the Ph.D. thesis:
Roland Wunderling, 
<em>Paralleler und Objektorientierter Simplex-Algorithmus</em>,
<a href="http://www.zib.de/PaperWeb/abstracts/TR-96-09">
ZIB technical report TR 96-09</a>, Berlin 1996</b>
</p>
<p>
If you are not applicable for the academic license, 
here are some notes on 
<a href="commercial.html">commercial licensing</a>.
</p>
<p>
<hr noshade>
<h3>Download</h3>
<p>
The latest Version is 1.2.1. Register and 
<a href="register.html">download</a> the complete source 
code and documentation.
</p>
<p>
Here are the 
<a href="http://elib.zib.de/pub/Packages/mp-testdata/lp/netlib-lp">Netlib</a>
LP files, already decompressed and assembled as an archive 
<a href="netlib.tar.gz">netlib.tar.gz</a>. These instances are used by the
<em>check</em> target in the Makefile. For comparison here is a directory 
with the log files of our <a href="results">results</a>.
</p>
<h3>Bugs</h3>
<p>
If you find one, it would be nice if you 
send a description together with a data file that shows the
problem or even better a working fix to 
<a href="mailto://koch@zib.de">koch@zib.de</a>.
</p>
<h3>Mailing list</h3>
<p>
<a href="mailto:soplex@zib.de">soplex@zib.de</a>.
To subscribe send "subscribe soplex" in the body to 
<a href="mailto:majordomo@zib.de">majordomo@zib.de</a>.
</p>
<h3>Documentation</h3>
<p>
You can browse the complete (that means this is all we have) 
<a href="html/index.html">documentation</a>. There you will find
information on how to compile, install, use, and modify SoPlex.
</p>
<hr noshade>
<h3>Links</h3>
<dl>
<dt><a href="http://plato.la.asu.edu/guide.html">
    Decision Tree for Optimization Software</a></dt>
<dt><a href="http://www-neos.mcs.anl.gov">
   NEOS Server for Optimization</a></dt>
<dt><a href="http://elib.zib.de/">
    The electronical library for mathematical software (eLib)</a></dt>
<dt><a href="http://optnet.itwm.uni-kl.de/opt-net/">
    Opt-Net</a></dt>
<dt><a href="http://www.cudenver.edu/~hgreenbe/glossary/glossary.html">
    Mathematical Programming Glossary© by Harvey J. Greenberg</a></dt>
<dt><a href="http://www.informs.org/Resources">
   INFORMS OR/MS Resource Collection</a></dt>
<dt><a href="http://www.cise.ufl.edu/~davis/sparse">
   University of Florida Sparse Matrix Collection</a></dt>
</dl>
<p>
<hr noshade>
<a href="http://www.anybrowser.org/campaign/">
<img src="images/anybrowser3.gif" alt="Best viewed with any browser" 
border="0" height="31" width="88" align="bottom"></a>
<a href="http://www.gimp.org">
<img src="images/gfx_by_gimp.gif" alt="Graphics by GIMP" 
border="0" height="36" width="90" align="bottom"></a> 
<img src="images/msfree.gif" alt="100% Microsoft Free" 
border="0" width="100" height="31" align="bottom">    
<a href="http://www.php.net">
<img src="images/phpsmall.gif" alt="Powered by PHP" 
border="0" height="31" width="88" align="bottom"></a>
<a href="http://petition.eurolinux.org">
<img src="images/patent_button.gif" alt="No ePatents" 
border="0" width="88" height="36" align="bottom"></a>    
<a href="http://www.gnu.org">
<img src="images/gnuspirit.png" alt="The GNU Project" 
 border="0" width="100" height="36" align="bottom"></a>
<a href="http://validator.w3.org/check/referer">
<img src="images/valid-html40.png" alt="Valid HTML 4.0!"
border="0" height="31" width="88" align="bottom"></a>
<hr noshade><address><font size=-1>
Last Update $Date: 2002/04/10 06:48:56 $ by
<a href="/personal/personal.pl?name=koch">Thorsten Koch</a>
<br>&copy; 2002 by Konrad-Zuse-Zentrum für Informationstechnik Berlin (ZIB)<br>
http://www.zib.de/Optimization/Software/Soplex/soplex.php
</font>
</address>
</body>
</html>


