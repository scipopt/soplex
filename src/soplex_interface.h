
#ifdef __cplusplus
extern "C" {
#endif

/** creates new SoPlex struct **/
void* SoPlex_create();

/** frees SoPlex struct **/
void SoPlex_free(void* soplex);

/** reads LP file in LP or MPS format according to READMODE parameter; returns true on success **/
int SoPlex_readInstanceFile(void* soplex, const char* filename);

/** reads basis information from filename and returns true on success **/
int SoPlex_readBasisFile(void* soplex, const char* filename);

/** reads settings from filename and returns true on success **/
int SoPlex_readSettingsFile(void* soplex, const char* filename);

/** clears the (floating point) LP **/
void SoPlex_clearLPReal(void* soplex);

/** returns number of rows **/
int SoPlex_numRows(void* soplex);

/** returns number of columns **/
int SoPlex_numCols(void* soplex);

/** enables rational solving mode  **/
void SoPlex_setRational(void* soplex);

/** sets boolean parameter value **/
void SoPlex_setBoolParam(void* soplex, int paramcode, int paramvalue);

/** sets integer parameter value **/
void SoPlex_setIntParam(void* soplex, int paramcode, int paramvalue);

/** sets real parameter value **/
void SoPlex_setRealParam(void* soplex, int paramcode, double paramvalue);

/** returns value of integer parameter **/
int SoPlex_getIntParam(void* soplex, int paramcode);

/** adds a single (floating point) column **/
void SoPlex_addColReal(void* soplex, double* colentries, int colsize, int nnonzeros, double objval,
                       double lb, double ub);

/** removes a single (floating point) column **/
void SoPlex_removeColReal(void* soplex, int colidx);

/** adds a single rational column **/
void SoPlex_addColRational(void* soplex, long* colnums, long* coldenoms, int colsize, int nnonzeros,
                           long objvalnum, long objvaldenom, long lbnum, long lbdenom, long ubnum, long ubdenom);

/** adds a single (floating point) row **/
void SoPlex_addRowReal(void* soplex, double* rowentries, int rowsize, int nnonzeros, double lb,
                       double ub);

/** removes a single (floating point) row **/
void SoPlex_removeRowReal(void* soplex, int rowidx);

/** adds a single rational row **/
void SoPlex_addRowRational(void* soplex, long* rownums, long* rowdenoms, int rowsize, int nnonzeros,
                           long lbnum, long lbdenom, long ubnum, long ubdenom);

/** gets primal solution **/
void SoPlex_getPrimalReal(void* soplex, double* primal, int dim);

/** Returns rational primal solution in a char pointer.
*   The caller needs to ensure the char array is freed.
**/
char* SoPlex_getPrimalRationalString(void* soplex, int dim);

/** gets dual solution **/
void SoPlex_getDualReal(void* soplex, double* dual, int dim);

/** gets reduced cost vector **/
void SoPlex_getRedCostReal(void* soplex, double* rc, int dim);

/** optimizes the given LP and returns solver status **/
int SoPlex_optimize(void* soplex);

/** returns the current solver status **/
int SoPlex_getStatus(void* soplex);

/** returns the time spent in last call to solve **/
double SoPlex_getSolvingTime(void* soplex);

/** returns the number of iteration in last call to solve **/
int SoPlex_getNumIterations(void* soplex);

/** changes objective function vector to obj **/
void SoPlex_changeObjReal(void* soplex, double* obj, int dim);

/** changes rational objective function vector to obj **/
void SoPlex_changeObjRational(void* soplex, long* objnums, long* objdenoms, int dim);

/** changes left-hand side vector for constraints to lhs **/
void SoPlex_changeLhsReal(void* soplex, double* lhs, int dim);

/** changes left-hand side of a row to lhs **/
void SoPlex_changeRowLhsReal(void* soplex, int rowidx, double lhs);

/** changes rational left-hand side vector for constraints to lhs **/
void SoPlex_changeLhsRational(void* soplex, long* lhsnums, long* lhsdenoms, int dim);

/** changes right-hand side vector for constraints to rhs **/
void SoPlex_changeRhsReal(void* soplex, double* rhs, int dim);

/** changes right-hand side of a row to rhs **/
void SoPlex_changeRowRhsReal(void* soplex, int rowidx, double rhs);

/** changes rational right-hand side vector for constraints to rhs **/
void SoPlex_changeRhsRational(void* soplex, long* rhsnums, long* rhsdenoms, int dim);

/** changes both sides for constraints to given lhs and rhs **/
void SoPlex_changeRangeReal(void* soplex, double* lhs, double* rhs, int dim);

/** changes both sides of a row to given lhs and rhs **/
void SoPlex_changeRowRangeReal(void* soplex, int rowidx, double lhs, double rhs);

/** write LP to file; LP or MPS format is chosen from the extension in filename **/
void SoPlex_writeFileReal(void* soplex, char* filename);

/** returns the objective value if a primal solution is available **/
double SoPlex_objValueReal(void* soplex);

/** Returns the rational objective value (as a string) if a primal solution is available.
*   The caller needs to ensure the char array is freed.
**/
char* SoPlex_objValueRationalString(void* soplex);

/** changes vectors of column bounds to lb and ub **/
void SoPlex_changeBoundsReal(void* soplex, double* lb, double* ub, int dim);

/** changes bounds of a column to lb and ub **/
void SoPlex_changeVarBoundsReal(void* soplex, int colidx, double lb, double ub);

/** changes rational bounds of a column to lbnum/lbdenom and ubnum/ubdenom **/
void SoPlex_changeVarBoundsRational(void* soplex, int colidx, long lbnum, long lbdenom, long ubnum,
                                    long ubdenom);

/** changes vector of lower bounds to lb **/
void SoPlex_changeLowerReal(void* soplex, double* lb, int dim);

/** changes lower bound of column to lb **/
void SoPlex_changeVarLowerReal(void* soplex, int colidx, double lb);

/** gets lower bound vector of columns into lb **/
void SoPlex_getLowerReal(void* soplex, double* lb, int dim);

/** gets objective vector into obj **/
void SoPlex_getObjReal(void* soplex, double* obj, int dim);

/** changes vector of upper bounds to ub **/
void SoPlex_changeUpperReal(void* soplex, double* ub, int dim);

/** changes upper bound of column to ub **/
void SoPlex_changeVarUpperReal(void* soplex, int colidx, double ub);

/** gets upper bound vector of columns into ub **/
void SoPlex_getUpperReal(void* soplex, double* ub, int dim);

/** returns status of row
 *  0 -> row is set to its upper bound
 *  1 -> row is set to its lower bound
 *  2 -> row is fixed to its identical bounds
 *  4 -> row is basic
 *  5 -> nothing known about basis status
 **/
int SoPlex_basisRowStatus(void* soplex, int rowidx);

/** returns status of column
 *  0 -> column is set to its upper bound
 *  1 -> column is set to its lower bound
 *  2 -> column is fixed to its identical bounds
 *  3 -> column is free and fixed to zero
 *  4 -> column is basic
 *  5 -> nothing known about basis status
 **/
int SoPlex_basisColStatus(void* soplex, int colidx);

/** get non-zero entries and indices of row i **/
void SoPlex_getRowVectorReal(void* soplex, int i, int* nnonzeros, long* indices, double* coefs);

/** get non-zero entries and indices of rational row i **/
void SoPlex_getRowVectorRational(void* soplex, int i, int* nnonzeros, long* indices, long* coefsnum,
                                 long* coefsdenom);

/** get lower and upper bounds of row i **/
void SoPlex_getRowBoundsReal(void* soplex, int i, double* lb, double* ub);

/** get rational lower and upper bounds of row i **/
void SoPlex_getRowBoundsRational(void* soplex, int i, long* lbnum, long* lbdenom, long* ubnum,
                                 long* ubdenom);

#ifdef __cplusplus
}
#endif
