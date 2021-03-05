
#include "soplex.h"

using namespace soplex;

void* SoPlex_create()
{
   SoPlex* so = new SoPlex();
   return so;
}

void SoPlex_free(void *soplex)
{
   SoPlex* so = (SoPlex*)(soplex);
   delete so;
}

int SoPlex_numRows(void *soplex)
{
   SoPlex* so = (SoPlex*)(soplex);
   return so->numRows();
}

int SoPlex_numCols(void *soplex)
{
   SoPlex* so = (SoPlex*)(soplex);
   return so->numCols();
}

void SoPlex_setRational(void *soplex)
{
   SoPlex* so = (SoPlex*)(soplex);
   so->setIntParam(SoPlex::READMODE, SoPlex::READMODE_RATIONAL);
   so->setIntParam(SoPlex::SOLVEMODE, SoPlex::SOLVEMODE_RATIONAL);
   so->setIntParam(SoPlex::CHECKMODE, SoPlex::CHECKMODE_RATIONAL);
   so->setIntParam(SoPlex::SYNCMODE, SoPlex::SYNCMODE_AUTO);
   so->setRealParam(SoPlex::FEASTOL, 0.0);
   so->setRealParam(SoPlex::OPTTOL, 0.0);
}

void SoPlex_setIntParam(void *soplex, int paramcode, int paramvalue)
{
   SoPlex* so = (SoPlex*)(soplex);
   so->setIntParam((SoPlex::IntParam)paramcode, paramvalue);
}
