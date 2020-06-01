#include"headFile.h"

void TranA_multiply_x(float sa[], unsigned int ija[], double x[], double b[],unsigned int n){
// b =A' * X; // attention index of x and b are start from zeros;
	unsigned int i,j,k;

	if (ija[1] != n+2)
	{
        printf("mismatched vector and matrix in sprstx\n");
        getchar();
        exit(-1);
	}
	for (i=1;i<=n;i++){
		b[i-1] = sa[i]*x[i-1];
	}

	for (i=1;i<=n;i++) {
		for (k=ija[i]; k<=(ija[i+1]-1); ++k) {
			j = ija[k];
			b[j-1] += sa[k]*x[i-1];
		}
	}
}
