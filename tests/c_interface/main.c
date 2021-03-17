#include <stdio.h>
#include <soplex_interface.h>
#include <assert.h>

int main(void)
{
   /* create LP via columns */
	void *soplex = SoPlex_create();
   double infty = 10e+20;
   double colentries1[] = {-1.0};
   double colentries2[] = {1.0};
   double lhs[] = {-10.0};
   double primal[] = {0.0,0.0};

   /* minimize */
   SoPlex_setIntParam(soplex, 0, -1);

   /* add columns */
   SoPlex_addColReal(soplex, colentries1, 1, 1, 1.0, 0.0, infty);
   SoPlex_addColReal(soplex, colentries2, 1, 1, 1.0, -infty, infty);
   assert(SoPlex_numRows(soplex) == 1);
   assert(SoPlex_numCols(soplex) == 2);

   /* set lhs of constraint */
   SoPlex_changeLhsReal(soplex, lhs, 1);

   /* optimize and get solution and objective value */
   int result = SoPlex_optimize(soplex);
   assert(result == 1);
   SoPlex_getPrimalReal(soplex, primal, 2);
   assert(primal[0] == 0.0 && primal[1] == -10.0);

	SoPlex_free(soplex);

   /* create LP via rows */
   void *soplex2 = SoPlex_create();
   double rowentries1[] = {-1.0, 1.0};
   double rowentries2[] = {1.0, 0.0};
   double rowentries3[] = {0.0, 1.0};
   double obj[] = {1.0, 1.0};

   /* minimize */
   SoPlex_setIntParam(soplex2, 0, -1);

   /* add row */
   SoPlex_addRowReal(soplex2, rowentries1, 2, 2, -10, infty);
   SoPlex_addRowReal(soplex2, rowentries2, 2, 1, 0.0, infty);
   SoPlex_addRowReal(soplex2, rowentries3, 2, 1, -infty, infty);
   assert(SoPlex_numRows(soplex2) == 3);
   assert(SoPlex_numCols(soplex2) == 2);

   /* add objective */
   SoPlex_changeObjReal(soplex2, obj, 2);

   /* optimize and get solution and objective value */
   result = SoPlex_optimize(soplex2);
   assert(result == 1);
   SoPlex_getPrimalReal(soplex2, primal, 2);
   assert(primal[0] == 0.0 && primal[1] == -10.0);

	SoPlex_free(soplex2);
}
