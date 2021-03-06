#include <soplex_interface.h>

int main(void)
{
	void *soplex = SoPlex_create();
	SoPlex_free(soplex);
	printf("Testing\n");
}