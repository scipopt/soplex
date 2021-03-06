
#ifdef __cplusplus
extern "C" {
#endif

void* SoPlex_create();

void SoPlex_free(void *soplex);

int SoPlex_numRows(void *soplex);

int SoPlex_numCols(void *soplex);

void SoPlex_setRational(void *soplex);

void SoPlex_setIntParam(void *soplex, int paramcode, int paramvalue);

#ifdef __cplusplus
}
#endif
