#include "soplex.h"
#include <iostream>

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

extern "C" void SoPlex_addColRational(void *soplex, long* colnums, long* coldenoms, int colsize, int nnonzeros, long objvalnum, long objvaldenom, long lbnum, long lbdenom, long ubnum, long ubdenom)
{
	SoPlex* so = (SoPlex*)(soplex);
	DSVectorRational col(nnonzeros);
	
	/* get rational lower bound */
	mpq_t lb;
    mpq_init(lb);
    mpq_set_si(lb, lbnum, lbdenom);
    Rational lower = lb;
    mpq_clear(lb);
	
	/* get rational upper bound */
	mpq_t ub;
    mpq_init(ub);
    mpq_set_si(ub, ubnum, ubdenom);
    Rational upper = ub;
    mpq_clear(ub);
	
    /* get rational objective value */	
	mpq_t obj;
    mpq_init(obj);
    mpq_set_si(obj, objvalnum, objvaldenom);
    Rational objval = obj;
    mpq_clear(obj);

	/* add nonzero entries to column vector */
	for( int i = 0; i < colsize; ++i )
	{
		if( colnums[i] != 0 )
		{
			/* get rational nonzero entry */
	        mpq_t r;
            mpq_init(r);
            mpq_set_si(r, colnums[i], coldenoms[i]);
            Rational colentry = r;
            mpq_clear(r);

			col.add(i, colentry);
		}
	}
	
	so->addColRational(LPColRational(objval, col, upper, lower));
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

extern "C" void SoPlex_addRowRational(void *soplex, long* rownums, long* rowdenoms, int rowsize, int nnonzeros, long lbnum, long lbdenom, long ubnum, long ubdenom)
{
	SoPlex* so = (SoPlex*)(soplex);
	DSVectorRational row(nnonzeros);

	/* get rational lower bound */
	mpq_t lb;
    mpq_init(lb);
    mpq_set_si(lb, lbnum, lbdenom);
    Rational lower = lb;
    mpq_clear(lb);

	/* get rational upper bound */	
	mpq_t ub;
    mpq_init(ub);
    mpq_set_si(ub, ubnum, ubdenom);
    Rational upper = ub;
    mpq_clear(ub);

	/* add nonzero entries to row vector */
	for( int i = 0; i < rowsize; ++i )
	{
		if( rownums[i] != 0 )
		{
			/* get rational nonzero entry */
	        mpq_t r;
            mpq_init(r);
            mpq_set_si(r, rownums[i], rowdenoms[i]);
			Rational rowentry = r;
            mpq_clear(r);

			row.add(i, rowentry);
		}
	}
	
	so->addRowRational(LPRowRational(lower, row, upper));
}

extern "C" void SoPlex_getPrimalReal(void *soplex, double* primal, int dim)
{
	SoPlex* so = (SoPlex*)(soplex);
	so->getPrimalReal(primal, dim);

    //for( int i = 0; i < dim; ++i)
    //    std::cout << primal[i];
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

extern "C" void SoPlex_changeObjRational(void *soplex, long* objnums, long* objdenoms, int dim)
{
   SoPlex* so = (SoPlex*)(soplex);
   Rational* objrational = new Rational [dim];

   /* create rational objective vector */
   for( int i = 0; i < dim; ++i )
   {
	    mpq_t r;
        mpq_init(r);
        mpq_set_si(r, objnums[i], objdenoms[i]);
	    Rational objentry = r;
        mpq_clear(r);
        objrational[i] = objentry;
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

extern "C" void SoPlex_changeLhsRational(void *soplex, long* lhsnums, long* lhsdenoms, int dim)
{
    SoPlex* so = (SoPlex*)(soplex);
    Rational* lhsrational = new Rational [dim];

    /* create rational lhs vector */
    for( int i = 0; i < dim; ++i)
    {
	    mpq_t r;
        mpq_init(r);
        mpq_set_si(r, lhsnums[i], lhsdenoms[i]);
	    Rational lhsentry = r;
        mpq_clear(r);
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

extern "C" void SoPlex_changeRhsRational(void *soplex, long* rhsnums, long* rhsdenoms, int dim)
{
    SoPlex* so = (SoPlex*)(soplex);
    Rational* rhsrational = new Rational [dim];

    /* create rational rhs vector */
    for( int i = 0; i < dim; ++i)
    {

	    mpq_t r;
        mpq_init(r);
        mpq_set_si(r, rhsnums[i], rhsdenoms[i]);
	    Rational rhsentry = r;
        mpq_clear(r);
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
    //std::cout << lb;
    //std::cout << '\n';
    //std::cout << colidx;
    //std::cout << '\n';
    //std::cout << ub;
    //std::cout << '\n';
    return so->changeBoundsReal(colidx, lb, ub);
}

extern "C" void SoPlex_changeVarBoundsRational(void *soplex, int colidx, long lbnum, long lbdenom, long ubnum, long ubdenom)
{
    SoPlex* so = (SoPlex*)(soplex);

    mpq_t lb;
    mpq_init(lb);
    mpq_set_si(lb, lbnum, lbdenom);
	Rational lower = lb;
    mpq_clear(lb);
    
    mpq_t ub;
    mpq_init(ub);
    mpq_set_si(ub, ubnum, ubdenom);
	Rational upper = ub;
    mpq_clear(ub);

    return so->changeBoundsRational(colidx, lower, upper);
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
