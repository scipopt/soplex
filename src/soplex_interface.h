#include "soplex/config.h"

#ifdef __cplusplus
extern "C" {
#endif

/** creates new SoPlex struct **/
SOPLEX_EXPORT
void* SoPlex_create();

/** frees SoPlex struct **/
SOPLEX_EXPORT
void SoPlex_free(void* soplex);

/** reads LP file in LP or MPS format according to READMODE parameter; returns true on success **/
SOPLEX_EXPORT
int SoPlex_readInstanceFile(void* soplex, const char* filename);

/** reads basis information from filename and returns true on success **/
SOPLEX_EXPORT
int SoPlex_readBasisFile(void* soplex, const char* filename);

/** reads settings from filename and returns true on success **/
SOPLEX_EXPORT
int SoPlex_readSettingsFile(void* soplex, const char* filename);

/** clears the (floating point) LP **/
SOPLEX_EXPORT
void SoPlex_clearLPReal(void* soplex);

/** returns number of rows **/
SOPLEX_EXPORT
int SoPlex_numRows(void* soplex);

/** returns number of columns **/
SOPLEX_EXPORT
int SoPlex_numCols(void* soplex);

/** enables rational solving mode  **/
SOPLEX_EXPORT
void SoPlex_setRational(void* soplex);

/** sets boolean parameter value **/
SOPLEX_EXPORT
void SoPlex_setBoolParam(void* soplex, int paramcode, int paramvalue);

/** sets integer parameter value **/
SOPLEX_EXPORT
void SoPlex_setIntParam(void* soplex, int paramcode, int paramvalue);

/** sets real parameter value **/
SOPLEX_EXPORT
void SoPlex_setRealParam(void* soplex, int paramcode, double paramvalue);

/** returns value of integer parameter **/
SOPLEX_EXPORT
int SoPlex_getIntParam(void* soplex, int paramcode);

/** adds a single (floating point) column **/
SOPLEX_EXPORT
void SoPlex_addColReal(void* soplex, double* colentries, int colsize, int nnonzeros, double objval,
                       double lb, double ub);

/** removes a single (floating point) column **/
SOPLEX_EXPORT
void SoPlex_removeColReal(void* soplex, int colidx);

/** adds a single rational column **/
SOPLEX_EXPORT
void SoPlex_addColRational(void* soplex, long* colnums, long* coldenoms, int colsize, int nnonzeros,
                           long objvalnum, long objvaldenom, long lbnum, long lbdenom, long ubnum, long ubdenom);

/** adds a single (floating point) row **/
SOPLEX_EXPORT
void SoPlex_addRowReal(void* soplex, double* rowentries, int rowsize, int nnonzeros, double lb,
                       double ub);

/** removes a single (floating point) row **/
SOPLEX_EXPORT
void SoPlex_removeRowReal(void* soplex, int rowidx);

/** adds a single rational row **/
SOPLEX_EXPORT
void SoPlex_addRowRational(void* soplex, long* rownums, long* rowdenoms, int rowsize, int nnonzeros,
                           long lbnum, long lbdenom, long ubnum, long ubdenom);

/** gets primal solution **/
SOPLEX_EXPORT
void SoPlex_getPrimalReal(void* soplex, double* primal, int dim);

/** Returns rational primal solution in a char pointer.
*   The caller needs to ensure the char array is freed.
**/
SOPLEX_EXPORT
char* SoPlex_getPrimalRationalString(void* soplex, int dim);

/** gets dual solution **/
SOPLEX_EXPORT
void SoPlex_getDualReal(void* soplex, double* dual, int dim);

/** gets reduced cost vector **/
SOPLEX_EXPORT
void SoPlex_getRedCostReal(void* soplex, double* rc, int dim);

/** optimizes the given LP and returns solver status **/
SOPLEX_EXPORT
int SoPlex_optimize(void* soplex);

/** returns the current solver status **/
SOPLEX_EXPORT
int SoPlex_getStatus(void* soplex);

/** returns the time spent in last call to solve **/
SOPLEX_EXPORT
double SoPlex_getSolvingTime(void* soplex);

/** returns the number of iteration in last call to solve **/
SOPLEX_EXPORT
int SoPlex_getNumIterations(void* soplex);

/** changes objective function vector to obj **/
SOPLEX_EXPORT
void SoPlex_changeObjReal(void* soplex, double* obj, int dim);

/** changes rational objective function vector to obj **/
SOPLEX_EXPORT
void SoPlex_changeObjRational(void* soplex, long* objnums, long* objdenoms, int dim);

/** changes left-hand side vector for constraints to lhs **/
SOPLEX_EXPORT
void SoPlex_changeLhsReal(void* soplex, double* lhs, int dim);

/** changes left-hand side of a row to lhs **/
SOPLEX_EXPORT
void SoPlex_changeRowLhsReal(void* soplex, int rowidx, double lhs);

/** changes rational left-hand side vector for constraints to lhs **/
SOPLEX_EXPORT
void SoPlex_changeLhsRational(void* soplex, long* lhsnums, long* lhsdenoms, int dim);

/** changes right-hand side vector for constraints to rhs **/
SOPLEX_EXPORT
void SoPlex_changeRhsReal(void* soplex, double* rhs, int dim);

/** changes right-hand side of a row to rhs **/
SOPLEX_EXPORT
void SoPlex_changeRowRhsReal(void* soplex, int rowidx, double rhs);

/** changes rational right-hand side vector for constraints to rhs **/
SOPLEX_EXPORT
void SoPlex_changeRhsRational(void* soplex, long* rhsnums, long* rhsdenoms, int dim);

/** changes both sides for constraints to given lhs and rhs **/
SOPLEX_EXPORT
void SoPlex_changeRangeReal(void* soplex, double* lhs, double* rhs, int dim);

/** changes both sides of a row to given lhs and rhs **/
SOPLEX_EXPORT
void SoPlex_changeRowRangeReal(void* soplex, int rowidx, double lhs, double rhs);

/** write LP to file; LP or MPS format is chosen from the extension in filename **/
SOPLEX_EXPORT
void SoPlex_writeFileReal(void* soplex, char* filename);

/** returns the objective value if a primal solution is available **/
SOPLEX_EXPORT
double SoPlex_objValueReal(void* soplex);

/** Returns the rational objective value (as a string) if a primal solution is available.
*   The caller needs to ensure the char array is freed.
**/
SOPLEX_EXPORT
char* SoPlex_objValueRationalString(void* soplex);

/** changes vectors of column bounds to lb and ub **/
SOPLEX_EXPORT
void SoPlex_changeBoundsReal(void* soplex, double* lb, double* ub, int dim);

/** changes bounds of a column to lb and ub **/
SOPLEX_EXPORT
void SoPlex_changeVarBoundsReal(void* soplex, int colidx, double lb, double ub);

/** changes rational bounds of a column to lbnum/lbdenom and ubnum/ubdenom **/
SOPLEX_EXPORT
void SoPlex_changeVarBoundsRational(void* soplex, int colidx, long lbnum, long lbdenom, long ubnum,
                                    long ubdenom);

/** changes vector of lower bounds to lb **/
SOPLEX_EXPORT
void SoPlex_changeLowerReal(void* soplex, double* lb, int dim);

/** changes lower bound of column to lb **/
SOPLEX_EXPORT
void SoPlex_changeVarLowerReal(void* soplex, int colidx, double lb);

/** gets lower bound vector of columns into lb **/
SOPLEX_EXPORT
void SoPlex_getLowerReal(void* soplex, double* lb, int dim);

/** gets objective vector into obj **/
SOPLEX_EXPORT
void SoPlex_getObjReal(void* soplex, double* obj, int dim);

/** changes vector of upper bounds to ub **/
SOPLEX_EXPORT
void SoPlex_changeUpperReal(void* soplex, double* ub, int dim);

/** changes upper bound of column to ub **/
SOPLEX_EXPORT
void SoPlex_changeVarUpperReal(void* soplex, int colidx, double ub);

/** gets upper bound vector of columns into ub **/
SOPLEX_EXPORT
void SoPlex_getUpperReal(void* soplex, double* ub, int dim);

/** returns status of row
 *  0 -> row is set to its upper bound
 *  1 -> row is set to its lower bound
 *  2 -> row is fixed to its identical bounds
 *  4 -> row is basic
 *  5 -> nothing known about basis status
 **/
SOPLEX_EXPORT
int SoPlex_basisRowStatus(void* soplex, int rowidx);

/** returns status of column
 *  0 -> column is set to its upper bound
 *  1 -> column is set to its lower bound
 *  2 -> column is fixed to its identical bounds
 *  3 -> column is free and fixed to zero
 *  4 -> column is basic
 *  5 -> nothing known about basis status
 **/
SOPLEX_EXPORT
int SoPlex_basisColStatus(void* soplex, int colidx);

/** get non-zero entries and indices of row i **/
SOPLEX_EXPORT
void SoPlex_getRowVectorReal(void* soplex, int i, int* nnonzeros, long* indices, double* coefs);

/** get non-zero entries and indices of rational row i **/
SOPLEX_EXPORT
void SoPlex_getRowVectorRational(void* soplex, int i, int* nnonzeros, long* indices, long* coefsnum,
                                 long* coefsdenom);

/** get lower and upper bounds of row i **/
SOPLEX_EXPORT
void SoPlex_getRowBoundsReal(void* soplex, int i, double* lb, double* ub);

/** get rational lower and upper bounds of row i **/
SOPLEX_EXPORT
void SoPlex_getRowBoundsRational(void* soplex, int i, long* lbnum, long* lbdenom, long* ubnum,
                                 long* ubdenom);

#ifdef __cplusplus
}
#endif
