#include "soplex.h"

using namespace soplex;

extern "C" void *SoPlex_create()
{
   SoPlex *so = new SoPlex();
   return so;
}

extern "C" void SoPlex_free(void *soplex)
{
   SoPlex *so = (SoPlex *)(soplex);
   delete so;
}

extern "C" void SoPlex_clearLPReal(void *soplex)
{
   SoPlex *so = (SoPlex *)(soplex);
   so->clearLPReal();
}

extern "C" int SoPlex_numRows(void *soplex)
{
   SoPlex *so = (SoPlex *)(soplex);
   return so->numRows();
}

extern "C" int SoPlex_numCols(void *soplex)
{
   SoPlex *so = (SoPlex *)(soplex);
   return so->numCols();
}

extern "C" void SoPlex_setRational(void *soplex)
{
   SoPlex *so = (SoPlex *)(soplex);
   so->setIntParam(SoPlex::READMODE, SoPlex::READMODE_RATIONAL);
   so->setIntParam(SoPlex::SOLVEMODE, SoPlex::SOLVEMODE_RATIONAL);
   so->setIntParam(SoPlex::CHECKMODE, SoPlex::CHECKMODE_RATIONAL);
   so->setIntParam(SoPlex::SYNCMODE, SoPlex::SYNCMODE_AUTO);
   so->setRealParam(SoPlex::FEASTOL, 0.0);
   so->setRealParam(SoPlex::OPTTOL, 0.0);
}

extern "C" void SoPlex_setIntParam(void *soplex, int paramcode, int paramvalue)
{
   SoPlex *so = (SoPlex *)(soplex);
   so->setIntParam((SoPlex::IntParam)paramcode, paramvalue);
}

extern "C" void SoPlex_addColReal(void *soplex, double* colentries, int colsize, int nnonzeros, double objval, double lb, double ub)
{
	SoPlex* so = (SoPlex*)(soplex);
	DSVector col(nnonzeros);

	/* add nonzero entries to column vector */
	for( int i = 0; i < colsize; ++i )
	{
		if( colentries[i] != 0.0 )
			col.add(i, colentries[i]);
	}
	
	so->addColReal(LPCol(objval, col, ub, lb));
}

extern "C" void SoPlex_addColRational(void *soplex, int* colnums, int* coldenoms, int colsize, int nnonzeros, int objvalnum, int objvaldenom, int lbnum, int lbdenom, int ubnum, int ubdenom)
{
	SoPlex* so = (SoPlex*)(soplex);
	DSVectorRational col(nnonzeros);
	
	/* get rational lower bound */
	Rational lb;
	lb = lbnum;
	lb /= lbdenom;
	
	/* get rational upper bound */
	Rational ub;
	ub = ubnum;
	ub /= ubdenom;

	/* get rational objective value */
	Rational objval;
	objval = objvalnum;
	objval /= objvaldenom;

	/* add nonzero entries to column vector */
	for( int i = 0; i < colsize; ++i )
	{
		if( colnums[i] != 0 )
		{
			/* get rational nonzero entry */
			Rational r;
			r = colnums[i];
			r /= coldenoms[i];

			col.add(i, r);
		}
	}
	
	so->addColRational(LPColRational(objval, col, ub, lb));
}

extern "C" void SoPlex_addRowReal(void *soplex, double* rowentries, int rowsize, int nnonzeros, double lb, double ub)
{
	SoPlex* so = (SoPlex*)(soplex);
	DSVector row(nnonzeros);

	/* add nonzero entries to row vector */
	for( int i = 0; i < rowsize; ++i )
	{
		if( rowentries[i] != 0.0 )
			row.add(i, rowentries[i]);
	}
	
	so->addRowReal(LPRow(lb, row, ub));
}

extern "C" void SoPlex_addRowRational(void *soplex, int* rownums, int* rowdenoms, int rowsize, int nnonzeros, int lbnum, int lbdenom, int ubnum, int ubdenom)
{
	SoPlex* so = (SoPlex*)(soplex);
	DSVectorRational row(nnonzeros);

	/* get rational lower bound */
	Rational lb;
	lb = lbnum;
	lb /= lbdenom;
	
	/* get rational upper bound */
	Rational ub;
	ub = ubnum;
	ub /= ubdenom;

	/* add nonzero entries to row vector */
	for( int i = 0; i < rowsize; ++i )
	{
		if( rownums[i] != 0 )
		{
			/* get rational nonzero entry */
			Rational r;
			r = rownums[i];
			r /= rowdenoms[i];

			row.add(i, r);
		}
	}
	
	so->addRowRational(LPRowRational(lb, row, ub));
}

extern "C" void SoPlex_getPrimalReal(void *soplex, double* primal, int dim)
{
	SoPlex* so = (SoPlex*)(soplex);
	so->getPrimalReal(primal, dim);
}

extern "C" char* SoPlex_getPrimalRationalString(void *soplex, int dim)
{
    SoPlex* so = (SoPlex*)(soplex);
    VectorRational primal(dim);
    std::string primalstring;
    so->getPrimalRational(primal);

    for( int i = 0; i < dim; ++i )
    {
        primalstring.append(rationalToString(primal[i], 0));
        primalstring.append(" ");
    }
    return const_cast<char*>( primalstring.c_str());
}

extern "C" void SoPlex_getDualReal(void *soplex, double* dual, int dim)
{
	SoPlex* so = (SoPlex*)(soplex);
    so->getDualReal(dual, dim);
}

extern "C" int SoPlex_optimize(void *soplex)
{
   SoPlex* so = (SoPlex*)(soplex);
   return so->optimize();
}

extern "C" void SoPlex_changeObjReal(void *soplex, double* obj, int dim)
{
   SoPlex* so = (SoPlex*)(soplex);
   Vector objective(dim, obj);
   return so->changeObjReal(objective);
}

extern "C" void SoPlex_changeObjRational(void *soplex, int* objnums, int* objdenoms, int dim)
{
   SoPlex* so = (SoPlex*)(soplex);
   Rational* objrational = new Rational [dim];

   /* create rational objective vector */
   for( int i = 0; i < dim; ++i )
   {
        Rational r;
        r = objnums[i];
        r /= objdenoms[i];
        objrational[i] = r;
    }

   VectorRational objective(dim, objrational);
   return so->changeObjRational(objective);
}

extern "C" void SoPlex_changeLhsReal(void *soplex, double* lhs, int dim)
{
    SoPlex* so = (SoPlex*)(soplex);
    Vector lhsvec(dim, lhs);
    return so->changeLhsReal(lhsvec);
}

extern "C" void SoPlex_changeLhsRational(void *soplex, int* lhsnums, int* lhsdenoms, int dim)
{
    SoPlex* so = (SoPlex*)(soplex);
    Rational* lhsrational = new Rational [dim];

    /* create rational lhs vector */
    for( int i = 0; i < dim; ++i)
    {
        Rational r;
        r = lhsnums[i];
        r /= lhsdenoms[i];
        lhsrational[i] = r;
    }

    VectorRational lhs(dim, lhsrational);
    return so->changeLhsRational(lhs);
}

extern "C" void SoPlex_changeRhsReal(void *soplex, double* rhs, int dim)
{
    SoPlex* so = (SoPlex*)(soplex);
    Vector rhsvec(dim, rhs);
    return so->changeRhsReal(rhsvec);
}

extern "C" void SoPlex_changeRhsRational(void *soplex, int* rhsnums, int* rhsdenoms, int dim)
{
    SoPlex* so = (SoPlex*)(soplex);
    Rational* rhsrational = new Rational [dim];

    /* create rational rhs vector */
    for( int i = 0; i < dim; ++i)
    {
        Rational r;
        r = rhsnums[i];
        r /= rhsdenoms[i];
        rhsrational[i] = r;
    }

    VectorRational rhs(dim, rhsrational);
    return so->changeRhsRational(rhs);
}

extern "C" void SoPlex_writeFileReal(void *soplex, char* filename)
{
    SoPlex* so = (SoPlex*)(soplex);
    so->writeFile(filename);
}

extern "C" double SoPlex_objValueReal(void *soplex)
{
    SoPlex* so = (SoPlex*)(soplex);
    return so->objValueReal();
}

extern "C" char* SoPlex_objValueRationalString(void *soplex)
{
    SoPlex* so = (SoPlex*)(soplex);
    return const_cast<char*>(rationalToString(so->objValueRational(), 0).c_str());
}

extern "C" void SoPlex_changeBoundsReal(void *soplex, double* lb, double* ub, int dim)
{
    SoPlex* so = (SoPlex*)(soplex);
    Vector lbvec(dim, lb);
    Vector ubvec(dim, ub);
    return so->changeBoundsReal(lbvec, ubvec);
}

extern "C" void SoPlex_changeVarBoundsReal(void *soplex, int colidx, double lb, double ub)
{
    SoPlex* so = (SoPlex*)(soplex);
    std::cout << lb;
    std::cout << '\n';
    std::cout << colidx;
    std::cout << '\n';
    std::cout << ub;
    std::cout << '\n';
    return so->changeBoundsReal(colidx, lb, ub);
}

extern "C" void SoPlex_changeVarBoundsRational(void *soplex, int colidx, int lbnum, int lbdenom, int ubnum, int ubdenom)
{
    SoPlex* so = (SoPlex*)(soplex);

    Rational lb;
    lb = lbnum;
    lb /= lbdenom;

    Rational ub;
    ub = ubnum;
    lb /= ubdenom;

    return so->changeBoundsRational(colidx, lb, ub);
}

extern "C" void SoPlex_changeVarUpperReal(void *soplex, int colidx, double ub)
{
    SoPlex* so = (SoPlex*)(soplex);
    return so->changeLowerReal(colidx, ub);
}

extern "C" void SoPlex_getUpperReal(void *soplex, double* ub, int dim)
{
    SoPlex* so = (SoPlex*)(soplex);
    Vector ubvec(dim, ub);

    so->getLowerReal(ubvec);

    for( int i = 0; i < dim; ++i )
        ub[i] = ubvec[i];
}
