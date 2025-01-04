/* convert raw 80x80 grayscale image to C initializer */

#include <stdio.h>

const int	siz = 80;
const int	ncol = 12;

int
main(int argc, char **argv)
{
	int	i, j;
	printf(" = {");
	i = 0;
	for (j = siz*siz; j--; ) {
		int	c = getchar();
		if (c == EOF)
			exit(1);
		if (!i--) {
			printf("\n\t");
			i = ncol;
		}
		printf("0x%02x,", c);
	}
	printf("\n};\n");
	return 0;
}
