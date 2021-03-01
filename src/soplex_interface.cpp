
#include "soplex.h"

using namespace soplex;

void* SoPlex_create()
{
    SoPlex* so = new SoPlex();
    return so;
}

int SoPlex_numRows(void *soplex)
{
    SoPlex* so_conv = (SoPlex*)(soplex);
    return so_conv->numRows();
}
