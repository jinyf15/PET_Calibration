#include "headFile.h"

void A_multiply_x(float sa[], unsigned int ija[], double x[], double b[],unsigned int n)
{
// b = A*x; // index of x and b start from 0; sa ijat start from 1;
	unsigned int i,k;

    if (ija[1] != n+2)
	{
        printf("mismatched vector and matrix in sprstx\n");
        getchar();
        exit(-1);
	}

	for (i=1;i<=n;i++)
	{
		b[i-1]=sa[i]*x[i-1];  // diagonal element;

		for (k=ija[i];k<=ija[i+1]-1;k++)
		{
                b[i-1] += sa[k]*x[ija[k]-1];
		}
	}
}
