#include <soplex_interface.h>

int main(void)
{
	void *soplex = SoPlex_create();
	int n = SoPlex_numRows(soplex);
	SoPlex_free(soplex);
	printf("Testing %d\n", n);
}
