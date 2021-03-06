
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
