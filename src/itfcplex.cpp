/*@ ----------------------------------------------------------------------------

   CPLEX api using SoPlex

   Identification:
   $Id: itfcplex.cpp,v 1.1 2001/11/30 15:49:20 bzfkocht Exp $

   Program history:
   $Log: itfcplex.cpp,v $
   Revision 1.1  2001/11/30 15:49:20  bzfkocht
   Der Rest des alten cplex interfaces.

// Revision 1.2  1995/07/05  09:36:32  bzfwunde
// minor changes
//
// Revision 1.1.1.1  1995/03/31  14:53:53  bzfwunde
// tested Version running with set packing
//
// Revision 1.1.1.1  1995/03/09  15:47:02  bzfwunde
// Initial version:
//     Tested for rowwise simplex
//     Error in columnwise part --- probably in doplex
//

    ----------------------------------------------------------------------------
 */
/*	\Section{Complex Methods}
 */
#ifndef PRICER
// #define	PRICER	SPxSteepPR
#define	PRICER	SPxHybridPR
#endif

#ifndef	COLUMN_SIMPLEX
#define	COLUMN_SIMPLEX
#endif


/*  Import system include files
 */
#include <assert.h>
#include <iostream.h>
#include <fstream.h>


/*  and class header files
 */
#ifndef	SUBDIR_INCLUDE

#include "spxhybridpr.hh"
#include "spxsteeppr.hh"
#include "spxfastrt.hh"
#include "slufactor.hh"
#include "spxweightst.hh"

#else 	// #SUBDIR_INCLUDE#

#include "spxhybridpr/spxhybridpr.hh"
#include "spxsteeppr/spxsteeppr.hh"
#include "spxfastrt/spxfastrt.hh"
#include "slufactor/slufactor.hh"
#include "spxweightst/spxweightst.hh"

#endif	// #SUBDIR_INCLUDE#

extern "C"
{
#include "spxcplex.h"
}

int	SPX_verbose = 0 ;

//@ ----------------------------------------------------------------------------

static void*	Malloc	( int size )
{
    void	*x = malloc(size + (size<=0)) ;
    if( x == 0 )
    {
	cerr << "ERROR: SPxCPLEX could not allocate memory\n" ;
	exit(-1) ;
    }
    return x ;
}

//@ ----------------------------------------------------------------------------
/*	\SubSection{SPxCPLEX class}
 */
class SPxCPLEX : public SoPlex
{
    SLUFactor		slu ;
    PRICER		price ;
    SPxFastRT		ratio ;
    SPxWeightST		start ;
public:
    static int		objSet ;
    static double	objULim ;
    static double	objLLim ;
    static double	maxTime ;
    static int		maxIter ;
    static int		readMPS ;

    void	factorize( )
		{
		    SoPlex::factorize( ) ;
		    if( SPX_verbose )
		    {
			cout.precision(16) ;
			cout << (type() == LEAVE ? "L " : "E ") ;
			cout << basis().iteration() << ":\t" ;
			cout << value() ;
			cout << "\t(" << shift() << ")\n" << flush;
		    }
		}

    int		terminate()
    		{
		    if( maxTime > 0  &&  time() > maxTime )
			return 1 ;
		    if( maxIter > 0  &&  iteration() > maxIter )
			return 1 ;
		    return SoPlex::terminate() ;
		}

    void	splitLP()
    		{
		    subcovectors.reSize( 1 ) ;
		}

    SPxCPLEX()
#ifdef	COLUMN_SIMPLEX
    : SoPlex( SoPlex::LEAVE, SoPlex::COLUMN )
#else
    : SoPlex( SoPlex::ENTER, SoPlex::ROW )
#endif
    {
	load( (SLinSolver*)&slu ) ;
	load( &ratio ) ;
	load( &price ) ;
	load( &start ) ;
    }
} ;

int	SPxCPLEX::objSet  =  0 ;
double	SPxCPLEX::objULim =  SPxLP::infinity ;
double	SPxCPLEX::objLLim = -SPxLP::infinity ;
double	SPxCPLEX::maxTime = -1 ;
int	SPxCPLEX::maxIter = -1 ;
int	SPxCPLEX::readMPS = 0 ;

//@ ----------------------------------------------------------------------------
/*	\SubSection{CPLEX Function Implementations}
 */

int openCPLEX( void ) 
{
  return(0);
}

int closeCPLEX( void )
{
  return 0 ;
}

int setlogfile(FILE *fp)
{
  return(0);
}

int lpmread (char *name, int *pmac, int *pmar, int *pobjsen, double **pobjx,
	     double **prhsx, char **psenx, int **pmatbeg, int **pmatcnt,
	     int **pmatind, double **pmatval, double **pbdl, double **pbdu,
	     char **pobjname, char **prhsname, char ***pcname, char **pcstore,
	     char ***prname, char **prstore, int *pmacsz, int *pmarsz,
	     int *pmatsz, unsigned *pcstorsz, unsigned *prstorsz, char** pctype)
{
#ifdef	DEBUG
    cout << "SPxCPLEX:	lpread\n" ;
#endif

    ifstream	file( name ) ;
    if( file.good() )
    {
	int		cnt, i, j ;
	SPxLP		lp ;
	const SPxLP&	clp = lp ;
	DIdxSet		intVars ;

	if( SPxCPLEX::readMPS )
	    lp.readMPS( file ) ;
	else
	    lp.readLP ( file, 0, 0, &intVars ) ;

	for( i = lp.nRows()-1, cnt = 0 ; i >= 0 ; --i )
	    cnt += clp.rowVector(i).size() ;
	
	*pmar     = lp.nRows() ;
	*pmac     = lp.nCols() ;
	*pmarsz   = lp.nRows() ;
	*pmacsz   = lp.nCols() ;
	*pmatsz   = cnt ;
	*pcstorsz = 0 ;
	*prstorsz = 0 ;

	*pobjx     = (double*)Malloc( (*pmacsz+2) * sizeof(double) ) ;
	*pbdl      = (double*)Malloc( (*pmacsz+2) * sizeof(double) ) ;
	*pbdu      = (double*)Malloc( (*pmacsz+2) * sizeof(double) ) ;
	*pmatbeg   = (int*)   Malloc( (*pmacsz+2) * sizeof(int) ) ;
	*pmatcnt   = (int*)   Malloc( (*pmacsz+2) * sizeof(int) ) ;
	*prhsx     = (double*)Malloc( (*pmarsz+2) * sizeof(double) ) ;
	*psenx     = (char*)  Malloc( (*pmarsz+2) * sizeof(char) ) ;
	*pmatind   = (int*)   Malloc( (*pmatsz+2) * sizeof(int) ) ;
	*pmatval   = (double*)Malloc( (*pmatsz+2) * sizeof(double) ) ;

	if( pctype )
	    *pctype = (char*)Malloc( (*pmacsz) * sizeof(char) ) ;

	*pcname    = 0 ;
	*prname    = 0 ;
	*pcstore   = 0 ;
	*prstore   = 0 ;

	*pobjsen   = (lp.spxSense() == SPxLP::MINIMIZE) ? 1 : -1 ;

	for( i = lp.nRows()-1 ; i >= 0 ; --i )
	{
	    (*prhsx)[i] = (lp.rhs(i) >= SPxLP::infinity) ? lp.lhs(i) : lp.rhs(i) ;

	    if( lp.rhs(i) >= SPxLP::infinity )
		(*psenx)[i] = 'G' ;
	    else if( lp.lhs(i) <= -SPxLP::infinity )
		(*psenx)[i] = 'L' ;
	    else if( lp.lhs(i) == lp.rhs(i) )
		(*psenx)[i] = 'E' ;
	    else
		(*psenx)[i] = 'R' ;
	}

	for( i = cnt = 0 ; i < lp.nCols() ; ++i )
	{
	    const SVector&	colVector = clp.colVector(i) ;

	    (*pbdu)[i]    = lp.upper(i) ;
	    (*pbdl)[i]    = lp.lower(i) ;
	    (*pobjx)[i]   = lp.obj(i) ;
	    (*pmatbeg)[i] = cnt ;
	    (*pmatcnt)[i] = colVector.size() ;

	    for( j = 0 ; j < colVector.size() ; ++j )
	    {
		(*pmatind)[cnt] = colVector.index(j) ;
		(*pmatval)[cnt] = colVector.value(j) ;
		++cnt ;
	    }
	}

	if( pctype )
	{
	    cerr << endl << "ERROR SPxCPLEX: " << __FILE__ << ":" << __LINE__ ;
	    cerr << "processing integer variables not yet implemented" << endl ;
	}

	return 0 ;
    }

    cout << endl << "ERROR SPxCPLEX: " << __FILE__ << ":" << __LINE__
	 << endl << endl;
    return 1 ;
}

int lpread  (char *name, int *pmac, int *pmar, int *pobjsen, double **pobjx,
	     double **prhsx, char **psenx, int **pmatbeg, int **pmatcnt,
	     int **pmatind, double **pmatval, double **pbdl, double **pbdu,
	     char **pobjname, char **prhsname, char ***pcname, char **pcstore,
	     char ***prname, char **prstore, int *pmacsz, int *pmarsz,
	     int *pmatsz, unsigned *pcstorsz, unsigned *prstorsz)
{
    return lpmread ( name, pmac, pmar, pobjsen, pobjx,
		     prhsx, psenx, pmatbeg, pmatcnt,
		     pmatind, pmatval, pbdl, pbdu,
		     pobjname, prhsname, pcname, pcstore,
		     prname, prstore, pmacsz, pmarsz,
		     pmatsz, pcstorsz, prstorsz, 0 ) ;
}

CPXLPptr loadprob (char *probname, int mac, int mar, int mae, int objsen,
		   double *objx, double *rhsx, char *senx, int *matbeg,
		   int *matcnt, int *matind, double *matval, double *bdl,
		   double *bdu, double *rngval, int *nrowind, int *etype,
		   int *enzbeg, int *enzcnt, int *enzind, double *enzval,
		   char *dataname, char *objname, char *rhsname, char *rngname,
		   char *bndname, char **cname, char *cstore, char **rname, char
		   *rstore, char **ename, char *estore, int macsz, int marsz,
		   int matsz, int maesz, int enzsz, unsigned cstorsz,
		   unsigned rstorsz, unsigned estorsz)
{
#ifdef	DEBUG
    cout << "SPxCPLEX:	loadprob\n" ;
#endif

    int		i, j ;
    LPColSet	cols(mac, matsz) ;
    LPRowSet	rows(mar) ;
    SPxCPLEX*	spx = new SPxCPLEX ;

    if
    (
	objx	== 0	||
	bdl	== 0	||
	bdu	== 0	||
	rhsx	== 0	||
	senx	== 0	||
	matbeg	== 0	||
	matcnt	== 0	||
	matind	== 0	||
	matval	== 0	||
	mae	!= 0 
    )
    {
	cout << endl << "ERROR SPxCPLEX: " << __FILE__ << ":" << __LINE__
	     << endl << endl;
	return 0 ;
    }


    DSVector	emptyVector(0) ;
    LPRow objRow( -SPxLP::infinity, emptyVector, SPxLP::infinity ) ;
    spx->addRow( objRow ) ;
    for( i = 0 ; i < mar ; ++i )
    {
	switch( senx[i] )
	{
	case 'L':
	    rows.add(-SPxLP::infinity, emptyVector, rhsx[i]) ;
	    break ;
	case 'G':
	    rows.add(rhsx[i], emptyVector, SPxLP::infinity) ;
	    break ;
	case 'E':
	    rows.add(rhsx[i], emptyVector, rhsx[i]) ;
	    break ;
	case 'R':
	    rows.add(rhsx[i], emptyVector, rhsx[i]) ;
	    assert( rngval ) ;
	    cout << endl << "ERROR SPxCPLEX: " << __FILE__ << ":" << __LINE__
		 << endl << endl;
	    break ;
	}
    }

    spx->addRows ( rows ) ;

    DSVector	colVector(mar) ;
    for( i = 0 ; i < mac ; ++i )
    {
	double*	val = &(matval[matbeg[i]]) ;
	int*	idx = &(matind[matbeg[i]]) ;
	colVector.clear() ;
	for( j = matcnt[i] - 1 ; j >= 0 ; --j )
	    colVector.add( (*idx++)+1, *val++ ) ;
	if( objx[i] )
	{
	    if( objsen == 1 )
		colVector.add( 0, -objx[i] ) ;
	    else
		colVector.add( 0,  objx[i] ) ;
	}
	cols.add( objx[i], bdl[i], colVector, bdu[i] ) ;
    }
    spx->changeSense( objsen == 1 ? SPxLP::MINIMIZE : SPxLP::MAXIMIZE ) ;
    spx->addCols( cols ) ;

    return	CPXLPptr(spx) ;
}

CPXLPptr loadlp(char   *probname,
	       int    mac,
	       int    mar,
	       int    objsen,
	       double *objx,
	       double *rhsx,
	       char   *senx,
	       int    *matbeg,
	       int    *matcnt,
	       int    *matind,
	       double *matval,
	       double *bdl,
	       double *bdu,
	       double *rngval,
	       int    macsz,
	       int    marsz,
	       int    matsz)
{
    return loadprob( probname, mac, mar, 0, objsen, objx, rhsx, senx, matbeg,
		   matcnt, matind, matval, bdl, bdu, rngval, 0, 0, 0, 0, 0, 0,
		   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, macsz, marsz, matsz, 0, 0, 0,
		   0, 0) ;
}


/*
void unloadprob (CPXLPptr* cplex)
{
#ifdef	DEBUG
    cout << "SPxCPLEX:	unloadprob\n" ;
#endif
    freeprob(cplex) ;
}
 */

void freeprob   (CPXLPptr* cplex)
{
#ifdef	DEBUG
    cout << "SPxCPLEX:	freeprob\n" ;
#endif
    SPxCPLEX*	solve = (SPxCPLEX*)*cplex ;
    *cplex = 0 ;
    delete solve ;
}

int setmlim  (int, int *, int *)
{
#ifdef	DEBUG
    cout << "Ah, good software does not need such a call ... ;-)\n" ;
#endif
    return 0 ;
}

int setnlim  (int, int *, int *)
{
#ifdef	DEBUG
    cout << "Ah, good software does not need such a call ... ;-)\n" ;
#endif
    return 0 ;
}

int setnzlim (int, int *, int *)
{
#ifdef	DEBUG
    cout << "Ah, good software does not need such a call ... ;-)\n" ;
#endif
    return 0 ;
}

void getmlim  (int *i)
{
#ifdef	DEBUG
    cout << "Ah, good software does not need such a call ... ;-)\n" ;
#endif
    *i = 9999999 ;
}

void getnlim  (int *i)
{
#ifdef	DEBUG
    cout << "Ah, good software does not need such a call ... ;-)\n" ;
#endif
    *i = 9999999 ;
}

void getnzlim (int *i)
{
#ifdef	DEBUG
    cout << "Ah, good software does not need such a call ... ;-)\n" ;
#endif
    *i = 9999999 ;
}


int setdpriind  (int, int *, int *)
{
#ifdef	DEBUG
    cout << "Ah, good software does not need such a call ... ;-)\n" ;
#endif
    return 0 ;
}

int getmac   (CPXLPptr cplex)
{
#ifdef	DEBUG
    cout << "SPxCPLEX:	getmac\n" ;
#endif

    SPxCPLEX*	solve = (SPxCPLEX*)cplex ;
    return solve->nCols() ;
}

int getmar   (CPXLPptr cplex)
{
#ifdef	DEBUG
    cout << "SPxCPLEX:	getmar\n" ;
#endif

    SPxCPLEX*	solve = (SPxCPLEX*)cplex ;
    return solve->nRows() - 1 ;
}

int getmat   (CPXLPptr cplex)
{
    int	i, cnt ;
    const SPxCPLEX*	solve = (const SPxCPLEX*)cplex ;

#ifdef	DEBUG
    cout << "SPxCPLEX:	getmat\n" ;
#endif

    i = solve->nRows()-1 ;
    for( cnt = 0 ; i > 0 ; --i )
	cnt += solve->rowVector(i).size() ;
    return cnt ;
}

int getbase  (CPXLPptr cplex, int *col, int *row)
{
    SPxCPLEX*			solve = (SPxCPLEX*)cplex ;
    const SPxBasis::Desc&	desc  = solve->basis().desc() ;

#ifdef	DEBUG
    cout << "SPxCPLEX:	getbase\n" ;
#endif

    int	i ;
    if(col)
    {
	for( i = solve->nCols()-1 ; i >= 0 ; --i )
	    switch( desc.colStatus(i) )
	    {
	    case SPxBasis::Desc::P_ON_LOWER:
		col[i] = 0 ;
		break ;
	    case SPxBasis::Desc::P_ON_UPPER:
	    case SPxBasis::Desc::P_FIXED:
		col[i] = 2 ;
		break ;
	    case SPxBasis::Desc::P_FREE:
		col[i] = 3 ;
		break ;
	    case SPxBasis::Desc::D_ON_UPPER:
	    case SPxBasis::Desc::D_ON_LOWER:
	    case SPxBasis::Desc::D_ON_BOTH:
	    case SPxBasis::Desc::D_UNDEFINED:
	    case SPxBasis::Desc::D_FREE:
		col[i] = 1 ;
		break ;
	    default:
		cout << endl << "ERROR SPxCPLEX: "
		     << __FILE__ << ":" << __LINE__ << endl << endl;
		break ;
	    }
    }

    if(row)
    {
	for( i = solve->nRows()-2 ; i >= 0 ; --i )
	    switch( desc.rowStatus(i+1) )
	    {
	    case SPxBasis::Desc::P_ON_LOWER:
	    case SPxBasis::Desc::P_ON_UPPER:
	    case SPxBasis::Desc::P_FIXED:
	    case SPxBasis::Desc::P_FREE:
		row[i] = 0 ;
		break ;
	    case SPxBasis::Desc::D_ON_UPPER:
	    case SPxBasis::Desc::D_ON_LOWER:
	    case SPxBasis::Desc::D_ON_BOTH:
	    case SPxBasis::Desc::D_UNDEFINED:
	    case SPxBasis::Desc::D_FREE:
		row[i] = 1 ;
		break ;
	    default:
		cout << endl << "ERROR SPxCPLEX: "
		     << __FILE__ << ":" << __LINE__ << endl << endl;
		break ;
	    }
    }

#ifdef	REAL_DEBUG
{
    int	i ;
    cout << "SPxCPLEX:	getbase\n" ;
    if( row )
    {
	cout << "rows:   " ;
	for( i = 1 ; i < nRows()-1 ; ++i )
	    cout << row[i] << ' ' ;
	cout << row[i] << '\n' ;
    }
    if( col )
    {
	cout << "cols:   " ;
	for( i = 0 ; i < solve->nCols()-1 ; ++i )
	    cout << col[i] << ' ' ;
	cout << col[i] << '\n' ;
    }
}
#endif

    return 0 ;
}

int getsense (CPXLPptr cplex, char *sns, int start, int end)
{
    int		i ;
    SPxCPLEX*	solve = (SPxCPLEX*)cplex ;

#ifdef	DEBUG
    cout << "SPxCPLEX:	getsense\n" ;
#endif

    if( start < 0  ||  end >= solve->nRows()-1 )
    {
	cout << endl << "ERROR SPxCPLEX: " << __FILE__ << ":" << __LINE__
	     << endl << endl;
	return 1 ;
    }

    for( i = start ; i <= end ; ++i )
    {
	if( solve->rhs(i+1) >= SPxLP::infinity )
	    sns[i-start] = 'G' ;
	else if( solve->lhs(i+1) <= -SPxLP::infinity )
	    sns[i-start] = 'L' ;
	else if( solve->lhs(i+1) == solve->rhs(i+1) )
	    sns[i-start] = 'E' ;
	else
	    sns[i-start] = 'R' ;
    }

    return 0 ;
}

int getrhs   (CPXLPptr cplex, double *rhs, int start, int end)
{
    int		i ;
    SPxCPLEX*	solve = (SPxCPLEX*)cplex ;

#ifdef	DEBUG
    cout << "SPxCPLEX:	getrhs\n" ;
#endif

    if( start < 0  ||  end >= solve->nRows()-1 )
    {
	cout << endl << "ERROR SPxCPLEX: " << __FILE__ << ":" << __LINE__
	     << endl << endl;
	return 1 ;
    }

    for( i = start ; i <= end ; ++i )
    {
	if( solve->rhs(i+1) >= SPxLP::infinity )
	    rhs[i-start] = solve->lhs(i+1) ;
	else 
	    rhs[i-start] = solve->rhs(i+1) ;
    }

    return 0 ;
}

int getobj   (CPXLPptr cplex, double *obj, int start, int end)
{
    int		i ;
    SPxCPLEX*	solve = (SPxCPLEX*)cplex ;

#ifdef	DEBUG
    cout << "SPxCPLEX:	getobj\n" ;
#endif

    if( start < 0  ||  end >= solve->nCols() )
    {
	cout << endl << "ERROR SPxCPLEX: " << __FILE__ << ":" << __LINE__
	     << endl << endl;
	return 1 ;
    }

    for( i = start ; i <= end ; ++i )
	obj[i-start] = solve->obj(i) ;

    return 0 ;
}

int getbdl   (CPXLPptr cplex, double *low, int start, int end)
{
    int		i ;
    SPxCPLEX*	solve = (SPxCPLEX*)cplex ;

#ifdef	DEBUG
    // cout << "SPxCPLEX:	getbdl\n" ;
#endif

    if( start < 0  ||  end >= solve->nCols() )
    {
	cout << endl << "ERROR SPxCPLEX: " << __FILE__ << ":" << __LINE__
	     << endl << endl;
	return 1 ;
    }

    for( i = start ; i <= end ; ++i )
	low[i-start] = solve->lower(i) ;

    return 0 ;
}

int getbdu   (CPXLPptr cplex, double *up, int start, int end)
{
    int		i ;
    SPxCPLEX*	solve = (SPxCPLEX*)cplex ;

#ifdef	DEBUG
    // cout << "SPxCPLEX:	getbdu\n" ;
#endif

    if( start < 0  ||  end >= solve->nCols() )
    {
	cout << endl << "ERROR SPxCPLEX: " << __FILE__ << ":" << __LINE__
	     << endl << endl;
	return 1 ;
    }

    for( i = start ; i <= end ; ++i )
	up[i-start] = solve->upper(i) ;

    return 0 ;
}

int chgobj (CPXLPptr cplex, int cnt, int *idx, double *obj)
{
    SPxCPLEX*	solve = (SPxCPLEX*)cplex ;

#ifdef	DEBUG
    cout << "SPxCPLEX:	chgbds\n" ;
#endif

    while( --cnt >= 0 )
	solve->changeObj( idx[cnt], obj[cnt] ) ;

    return 0 ;
}

int chgbds (CPXLPptr cplex, int cnt, int *index, char *lu, double *bd)
{
    SPxCPLEX*	solve = (SPxCPLEX*)cplex ;

#ifdef	DEBUG
    cout << "SPxCPLEX:	chgbds\n" ;
#endif

    while( --cnt >= 0 )
    {
	int	idx = index[cnt] ;
	if( idx < 0  ||  idx >= solve->nCols() )
	{
	    cout << endl << "ERROR SPxCPLEX: " << __FILE__ << ":" << __LINE__
		 << endl << endl;
	    return 1 ;
	}

	switch( lu[cnt] )
	{
	case 'U':
	    solve->changeUpper( idx, bd[cnt] ) ;
	    break ;
	case 'L':
	    solve->changeLower( idx, bd[cnt] ) ;
	    break ;
	case 'B':
	    solve->changeUpper( idx, bd[cnt] ) ;
	    solve->changeLower( idx, bd[cnt] ) ;
	    break ;
	default:
	    cout << endl << "ERROR SPxCPLEX: " << __FILE__ << ":" << __LINE__
		 << endl << endl;
	    return 1 ;
	}
    }

    return 0 ;
}

int getrows  (CPXLPptr cplex, int *nzcnt, int *beg, int *ind, double *val,
	      int size, int *surplus, int begin, int end)
{
    int	i, j, cnt ;
    const SPxCPLEX*	solve = (const SPxCPLEX*)cplex ;

#ifdef	DEBUG
    cout << "SPxCPLEX:	getrows\n" ;
#endif

    for( i = begin, cnt = 0 ; i <= end ; ++i )
    {
	const SVector&	row = solve->rowVector(i+1) ;	// row 0 contains obj
	if( size - cnt < row.size() )
	    break ;
	beg[i-begin] = cnt ;
	for( j = 0 ; j < row.size() ; ++j )
	{
	    ind[cnt] = row.index(j) ;
	    val[cnt] = row.value(j) ;
	    cnt++ ;
	}
    }

    if(nzcnt)
	*nzcnt = cnt ;

    if( i > end )
    {
	if( surplus )
	*surplus = size - cnt ;
	return 0 ;
    }

    if( surplus )
    {
	*surplus = 0 ;
	for( ; i <= end ; ++i )
	    *surplus -= solve->rowVector(i+1).size() ;
    }

    return 1 ;
}

int delsetrows (CPXLPptr cplex, int *del)
{
#ifdef	DEBUG
    cout << "SPxCPLEX:	delsetrows\n" ;
#endif

    SPxCPLEX*	solve = (SPxCPLEX*)cplex ;
    int		rnum  = solve->nRows()-1 ;
    DataArray<int>	rem( rnum+1 ) ;
    int	i ;

    for( i = rnum-1 ; i >= 0 ; --i )
	if( del[i] )
	    rem[i+1] = -1 ;
	else
	    rem[i+1] = 0 ;
    rem[0] = 0 ;

    solve->removeRows( rem ) ;
    return 0 ;
}

int addrows (CPXLPptr cplex, int ccnt, int rcnt, int nzcnt,
	     double *rhs, char *sns,
             int *beg, int *ind, double *val,
	     char **cnames, char **rnames)
{
#ifdef	DEBUG
    cout << "SPxCPLEX:	addrows\n" ;
#endif

    int		i ;
    LPRowSet	rset(rcnt, nzcnt) ;
    SPxCPLEX*	solve = (SPxCPLEX*)cplex ;

    DSVector	row( solve->nCols() ) ;

    for( i = 0 ; i < rcnt ; ++i )
    {
	int	len ;
	if( i < rcnt-1 )
	    len = beg[i+1] - beg[i] ;
	else
	    len = nzcnt - beg[i] ;

	row.clear() ;
	double*	_val = &val[beg[i]] ;
	int*	_idx = &ind[beg[i]] ;
	while( len-- )
	    row.add( *_idx++, *_val++ ) ;

	switch( sns[i] )
	{
	case 'L':
	    rset.add(-SPxLP::infinity, row, rhs[i]) ;
	    break ;
	case 'G':
	    rset.add(rhs[i], row, SPxLP::infinity) ;
	    break ;
	case 'E':
	case 'R':
	    rset.add(rhs[i], row, rhs[i]) ;
	    break ;
	}
    }

    solve->addRows( rset ) ;
    return 0 ;
}

int delrows (CPXLPptr cplex, int begin, int end )
{
    int	i ;
    SPxCPLEX*		solve = (SPxCPLEX*)cplex ;
    DataArray<int>	del( solve->nRows() ) ;

    for( i = del.size()-1 ; i >= 0 ; --i )
	del[i] = 0 ;
    for( i = begin ; i <= end ; ++i )
	del[i+1] = -1 ;			// obj is 0-th row!!
    solve->removeRows( (int*)del ) ;

    return 0 ;
}

int addcols (CPXLPptr cplex, int ccnt, int nzcnt, double* objx,
	     int *cmatbeg, int *cmatind, double *cmatval,
	     double *bdl, double *bdu, char **cname)
{
#ifdef	DEBUG
    cout << "SPxCPLEX:	addcols\n" ;
#endif

    int		i ;
    LPColSet	cset(ccnt, nzcnt) ;
    SPxCPLEX*	solve = (SPxCPLEX*)cplex ;

    DSVector	col( solve->nRows() ) ;

    for( i = 0 ; i < ccnt ; ++i )
    {
	int	len ;
	if( i < ccnt-1 )
	    len = cmatbeg[i+1] - cmatbeg[i] ;
	else
	    len = nzcnt - cmatbeg[i] ;

	col.clear() ;
	double*	_val = &cmatval[cmatbeg[i]] ;
	int*	_idx = &cmatind[cmatbeg[i]] ;
	while( len-- )
	    col.add( 1+(*_idx++), *_val++ ) ;

	if( objx[i] != 0 )
	    col.add( 0, objx[i] ) ;

	cset.add(objx[i], bdl[i], col, bdu[i]) ;
    }

    solve->addCols( cset ) ;
    return 0 ;
}

int delcols (CPXLPptr cplex, int begin, int end )
{
    int	i ;
    SPxCPLEX*		solve = (SPxCPLEX*)cplex ;
    DataArray<int>	del( solve->nRows() ) ;

    for( i = del.size()-1 ; i >= 0 ; --i )
	del[i] = 0 ;
    for( i = begin ; i <= end ; ++i )
	del[i] = -1 ;
    solve->removeCols( (int*)del ) ;

    return 0 ;
}

int optimize (CPXLPptr p)
{
    return dualopt(p) ;
}

int dualopt (CPXLPptr cplex)
{
#ifdef	DEBUG
    cout << "SPxCPLEX:	dualopt\n" ;
#endif
    SPxCPLEX*	solve = (SPxCPLEX*)cplex ;

    if( solve->basis().status() <= 0 )
	solve->setType( SoPlex::ENTER );

    if( SPxCPLEX::objSet )
    {
	if( SPxCPLEX::objULim != solve->rhs(0) )
	    solve->changeRhs( 0, SPxCPLEX::objULim ) ;
	if( SPxCPLEX::objLLim != solve->lhs(0) )
	    solve->changeLhs( 0, SPxCPLEX::objLLim ) ;
    }

#ifdef	DEBUG
    ofstream	lpfile( "tmp.lp" ) ;
    lpfile << *solve << endl ;
#endif

    solve->solve() ;
    return 0 ;
}

int solution  (CPXLPptr cplex, int *status, double *obj, double *primal,
	       double *dual, double *slack, double *redcost)
{
    double	x, y ;
    int		i ;
    SPxCPLEX*	solve = (SPxCPLEX*)cplex ;

#ifdef	DEBUG
    cout << "SPxCPLEX:	solution\n" ;
#endif

    if(status)
    {
	switch( solve->status() )
	{
	case LPSolver::INFEASIBLE:
	    *status = CPX_INFEASIBLE ;
	    break ;
	case LPSolver::UNBOUNDED:
	    *status = CPX_UNBOUNDED ;
	    break ;
	case LPSolver::SOLVED:
	    *status = CPX_OPTIMAL ;
	    break ;
	case LPSolver::PRIMAL:
	    *status = 5 ;
	    break ;
	case LPSolver::DUAL:
	    *status = 7 ;
	    break ;
	default:
	    *status = -1 ;
	    break ;
	}
    }

    if(obj)	*obj = solve->value() ;
    if(primal)
    {
	Vector	tmp(solve->nCols(), primal) ;
	solve->getPrimal ( tmp ) ;
    }
    if(redcost)
    {
	Vector	tmp(solve->nCols(), redcost) ;
	solve->getRdCost ( tmp ) ;
    }

    if( slack || dual )
    {
	int	rnum  = solve->nRows() ;
	DVector	tmp( rnum ) ;
	if(slack)
	{
	    solve->getSlacks( tmp ) ;
	    for( i = rnum-1 ; i > 0 ; --i )
	    {
		x = solve->rhs(i) ;
		y = -tmp[i] ;
		slack[i-1] = y + ((x < SPxLP::infinity) ? x : solve->lhs(i)) ; 
	    }
	}
	if(dual)
	{
	    solve->getDual( tmp ) ;
	    for( i = rnum-1 ; i > 0 ; --i )
		dual[i-1] = tmp[i] ;
	}
    }

#ifdef	REAL_DEBUG
    {
	int	i ;
	if(primal)
	{
	    cout << "primal        = (" ;
	    for (i = 0 ; i < solve->nCols() ; ++i)
		cout << primal[i] << ", " ;
	    cout << primal[i] << ")\n" ;
	}
	if(redcost)
	{
	    cout << "reduced costs = (" ;
	    for (i = 0 ; i < solve->nCols() ; ++i)
		cout << redcost[i] << ", " ;
	    cout << redcost[i] << ")\n" ;
	}
	if(slack)
	{
	    cout << "slacks        = (" ;
	    for (i = 0 ; i < solve->nRows()-1 ; ++i)
		cout << slack[i] << ", " ;
	    cout << slack[i] << ")\n" ;
	}
	if(dual)
	{
	    cout << "dual          = (" ;
	    for (i = 0 ; i < solve->nRows()-1 ; ++i)
		cout << dual[i] << ", " ;
	    cout << dual[i] << ")\n" ;
	}
    }
#endif

    return 0 ;
}

int delsetcols (CPXLPptr cplex, int *del)
{
#ifdef	DEBUG
    cout << "SPxCPLEX:	delsetcols\n" ;
#endif

    int		i ;
    SPxCPLEX*	solve = (SPxCPLEX*)cplex ;

    for( i = solve->nCols()-1 ; i >= 0 ; --i )
	if( del[i] )
	    del[i] = -1 ;

    solve->removeCols( del ) ;

    for( i = solve->nCols()-1 ; i >= 0 ; --i )
	del[i] = ( del[i] < 0 ) ? 1 : 0 ;

    return 0 ;
}

int lpwrite   (CPXLPptr cplex, char *name)
{
#ifdef	DEBUG
    cout << "SPxCPLEX:	lpwrite\n" ;
#endif

    SPxCPLEX*	solve = (SPxCPLEX*)cplex ;
    ofstream	out( name ) ;
    if( out.good() )
    {
	out << *solve ;
	return 0 ;
    }

    cout << endl << "ERROR SPxCPLEX: " << __FILE__ << ":" << __LINE__
	 << endl << endl;
    return 1 ;
}

int mpsread (char *name, int *pmac, int *pmar, int *pmae, int *pobjsen,
	     double **pobjx, double **prhsx, char **psenx, int **pmatbeg,
	     int **pmatcnt, int **pmatind, double **pmatval, double **pbdl,
	     double **pbdu, double **prngval, int **pnrowind, int **petype,
	     int **penzbeg, int **penzcnt, int **penzind, double **penzval,
	     char **pdataname, char **pobjname, char **prhsname, char **prngname,
	     char **pbndname, char ***pcname, char **pcstore, char ***prname,
	     char **prstore, char ***pename, char **pestore, int *pmacsz, int
	     *pmarsz, int *pmatsz, int *pmaesz, int *penzsz, unsigned *pcstorsz,
	     unsigned *prstorsz, unsigned *pestorsz)
{
#ifdef	DEBUG
    cout << "SPxCPLEX:	mpsread\n" ;
#endif
    int	r ;
    SPxCPLEX::readMPS = 1 ;
    r = lpread(name, pmac, pmar, pobjsen, pobjx, prhsx, psenx, pmatbeg, pmatcnt,
	     pmatind, pmatval, pbdl, pbdu, pobjname, prhsname, pcname, pcstore,
	     prname, prstore, pmacsz, pmarsz, pmatsz, pcstorsz, prstorsz) ;
    SPxCPLEX::readMPS = 0 ;
    return r ;
}

int setscr_ind  (int val)
{
#ifdef	DEBUG
    cerr << "SPxCPLEX (setscr_ind): not implemented ... \n" ;
#endif	// DEBUG
    return 0 ;
}

int setobjulim (double x, double *y, double *z)
{
#ifdef	DEBUG
    cout << "SPxCPLEX:	setobjulim\n" ;
#endif	// DEBUG
    SPxCPLEX::objULim = x ;
    SPxCPLEX::objSet  = 1 ;
    return 0 ;
}

int	setobjllim  (double x, double *a, double *b)
{
#ifdef	DEBUG
    cout << "SPxCPLEX:	setobjllim\n" ;
#endif	// DEBUG
    SPxCPLEX::objLLim = x ;
    SPxCPLEX::objSet  = 1 ;
    return 0 ;
}

void	getobjulim( double *x )
{
    *x = SPxCPLEX::objULim ;
}

void	getobjllim( double *x )
{
    *x = SPxCPLEX::objLLim ;
}


int setcraind  ( int, int*, int* )
{
#ifdef	DEBUG
    cerr << "SPxCPLEX (setcraind): crash procedures are not implemented ... \n" ;
#endif	// DEBUG
    return 0 ;
}

int setadvind  ( int )
{
#ifdef	DEBUG
    cerr << "SPxCPLEX (setadvind): basis advance is default in SPxCPLEX ... \n" ;
#endif	// DEBUG
    return 0 ;
}

int savread (char *x01, int *x02, int *x03, int *x04, int *x05, double **x06,
             double **x07, char **x08, int **x09, int **x10, int **x11,
             double **x12, double **x13, double **x14, double **x15, int **x16,
             int **x17, int **x18, int **x19, int **x20, double **x21,
             char **x22, char **x23, char **x24, char **x25, char **x26,
             char ***x27, char **x28, char ***x29, char **x30, char ***x31,
             char **x32, int *x33, int *x34, int *x35, int *x36, int *x37,
             unsigned *x38, unsigned *x39, unsigned *x40, int **x41, int **x42)
{
    cerr << "Sorry, cannot read CPLEX sav format, aborting ...\n" ;
    return -1 ;
}

int loadbase(CPXLPptr cplex, int *cstat, int *rstat)
{
    int			i ;
    SPxBasis::Desc	desc ;
    SPxCPLEX*		solver = (SPxCPLEX*)cplex ;
    int			rnum  = solver->nRows() ;

    desc.reSize( solver->nRows(), solver->nCols() ) ;

    if( solver->basis().solver() != solver )
	((SPxBasis*)&solver->basis())->load( solver ) ;

#ifdef	DEBUG
    cout << "SPxCPLEX:	loadbase\n" ;
#endif

    for( i = rnum-1 ; i-- ; )		// obj limit in row 0!!!
    {
	switch(rstat[i])
	{
	case 0:
	case 2:
	    if( solver->rhs(i+1) == solver->lhs(i+1) )
		desc.rowStatus(i+1) = SPxBasis::Desc::P_FIXED ;
	    else if( solver->rhs(i+1) >= SPxLP::infinity )
	    {
		assert( solver->lhs(i+1) > -SPxLP::infinity ) ;
		desc.rowStatus(i+1) = SPxBasis::Desc::P_ON_LOWER ;
	    }
	    else if( solver->lhs(i+1) <= -SPxLP::infinity )
	    {
		assert( solver->rhs(i+1) < SPxLP::infinity ) ;
		desc.rowStatus(i+1) = SPxBasis::Desc::P_ON_UPPER ;
	    }
	    else
	    {
		cout << endl << "ERROR SPxCPLEX: " << __FILE__ << ":" << __LINE__
		     << endl << endl;
	    }
	    break ;
	case 1:
	    desc.rowStatus(i+1) = solver->dualRowStatus(i+1) ;
	    break ;
	default:
	    cout << endl << "ERROR SPxCPLEX: " << __FILE__ << ":" << __LINE__
		 << endl << endl;
	    return 1 ;
	}
    }
    desc.rowStatus(0) = solver->dualRowStatus(0) ;

    for( i = solver->nCols()-1 ; i >= 0 ; i-- )
    {
	switch(cstat[i])
	{
	case 0:
	    if( solver->upper(i) == solver->lower(i) )
		desc.colStatus(i) = SPxBasis::Desc::P_FIXED ;
	    else
		desc.colStatus(i) = SPxBasis::Desc::P_ON_LOWER ;
	    break ;
	case 2:
	    if( solver->upper(i) == solver->lower(i) )
		desc.colStatus(i) = SPxBasis::Desc::P_FIXED ;
	    else
		desc.colStatus(i) = SPxBasis::Desc::P_ON_UPPER ;
	    break ;
	case 1:
	    desc.colStatus(i) = solver->dualColStatus(i) ;
	    break ;
	case 3:
	    desc.colStatus(i) = SPxBasis::Desc::P_FREE ;
	    break ;
	default:
	    cout << endl << "ERROR SPxCPLEX: " << __FILE__ << ":" << __LINE__
		 << endl << endl;
	    return 1 ;
	}
    }

    solver->load( desc ) ;

    return 0 ;
}

int settilim( double time, double *ptoosmall, double *ptoobig )
{
    SPxCPLEX::maxTime = time ;
    return 0 ;
}

void   gettilim    (double *t)
{
    *t = SPxCPLEX::maxTime ;
}

int	getobjsen   (CPXLPptr cplex)
{
    SPxCPLEX*	solver = (SPxCPLEX*)cplex ;
    if( solver->sense() == LPSolver::MAXIMIZE )
	return -1 ;
    else
	return 1 ;
}

int	getobjval   (CPXLPptr lp, double *z)
{
    assert(lp);
    assert(z);
    // this suffices at the moment:
    return solution(lp, 0, z, 0, 0, 0, 0);
}

int	getx        (CPXLPptr lp, double *x, int from, int upto)
{
    assert(lp);
    assert(x);
    assert(0 <= from && from <= upto && upto < getmac(lp));
    // this suffices at the moment:
    assert(from == 0 && upto == getmac(lp)-1);
    return solution(lp, 0, 0, x, 0, 0, 0);
}

int	getpi       (CPXLPptr lp, double *pi, int from, int upto)
{
    assert(lp);
    assert(pi);
    assert(0 <= from && from <= upto && upto < getmac(lp));
    // this suffices at the moment:
    assert(from == 0 && upto == getmar(lp)-1);
    return solution(lp, 0, 0, 0, pi, 0, 0);
}

int	getslack    (CPXLPptr lp, double *slack, int from, int upto)
{
    assert(lp);
    assert(slack);
    assert(0 <= from && from <= upto && upto < getmac(lp));
    // this suffices at the moment:
    assert(from == 0 && upto == getmar(lp)-1);
    return solution(lp, 0, 0, 0, 0, slack, 0);
}

int	getdj       (CPXLPptr lp, double *dj, int from, int upto)
{
    assert(lp);
    assert(dj);
    assert(0 <= from && from <= upto && upto < getmac(lp));
    // this suffices at the moment:
    assert(from == 0 && upto == getmac(lp)-1);
    return solution(lp, 0, 0, 0, 0, 0, dj);
}

int	getmethod   (CPXLPptr lp)
{
    assert(lp);
    // provisional: return CPXALG_PRIMAL unless solution() doesn't swap
    // CPX_INFEASIBLE <-> CPX_UNBOUNDED after calling dualopt():
    return CPXALG_PRIMAL;
    //return CPXALG_DUAL;
}


void	setitfoind(int, int*, int*)
{}

void	getitfoind(int*)
{}

void	getdpriind(int*)
{}

void	getadvind(int*)
{}

int	setitlim    (int lim, int *x, int *y)
{
    SPxCPLEX::maxIter = lim ;
    return 0 ;
}

void	getitlim    (int *lim )
{
    *lim = SPxCPLEX::maxIter ;
}

int getitc (CPXLPptr lp)
{
    SPxCPLEX*	solver = (SPxCPLEX*)lp ;
    return solver->basis().iteration() ;
}

int getcols (CPXLPptr lp, int *nzcnt, int *cmatbeg, int *cmatind,
	     double *cmatval, int cmatsz, int *surplus, int begin,
	     int end)
{
    int			i, j, n ;
    const SPxCPLEX*	solver = (const SPxCPLEX*)lp ;

#ifdef	DEBUG
    cout << "SPxCPLEX:	getcols\n" ;
#endif

    for( i = begin, n = 0 ; i <= end ; ++i )
    {
	const SVector&	col = solver->colVector(i) ;
	if( cmatsz - n < col.size() )
	{
	    if( surplus )
	    {
		*surplus = 0 ;
		for( ; i <= end ; ++i )
		    *surplus -= solver->colVector(i).size() ;
	    }
	    return 0 ;
	}
	cmatbeg[i-begin] = n ;
	for( j = 0 ; j < col.size() ; ++j )
	{
	    if( col.index(j) > 0 )			// row 0 is objective
	    {
		cmatind[n] = col.index(j) - 1 ;
		cmatval[n] = col.value(j) ;
		n++ ;
	    }
	}
    }

    if( nzcnt )
	*nzcnt = n ;
    if( surplus )
	*surplus = cmatsz - n ;

    return 1 ;
}
