#include <stdio.h>
#include <soplex_interface.h>
#include <assert.h>

void test_real(void)
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
   assert(SoPlex_objValueReal(soplex) == -10.0);

	SoPlex_free(soplex);

   /* create LP via rows */

   void *soplex2 = SoPlex_create();
   double rowentries1[] = {-1.0, 1.0};
   double lb[] = {0.0, -infty};
   double ub[] = {infty, infty};
   double obj[] = {1.0, 1.0};

   /* minimize */
   SoPlex_setIntParam(soplex2, 0, -1);

   /* add row */
   SoPlex_addRowReal(soplex2, rowentries1, 2, 2, -10.0, infty);

   /* add variable bounds */
   SoPlex_changeBoundsReal(soplex2, lb, ub, 2);
   assert(SoPlex_numRows(soplex2) == 1);
   assert(SoPlex_numCols(soplex2) == 2);

   /* add objective */
   SoPlex_changeObjReal(soplex2, obj, 2);

   /* optimize and get solution and objective value */
   result = SoPlex_optimize(soplex2);
   assert(result == 1);
   SoPlex_getPrimalReal(soplex2, primal, 2);
   assert(primal[0] == 0.0 && primal[1] == -10.0);
   assert(SoPlex_objValueReal(soplex2) == -10.0);

	SoPlex_free(soplex2);
}

int main(void)
{
   printf("testing real... \n");
   test_real();

   printf("\n");
   printf("testing rational... \n");
}
