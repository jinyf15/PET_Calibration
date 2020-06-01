#include <math.h>
#define NRANSI
#include "nrutil.h"
#define ITMAX 200
#define EPS 1.0e-10
#define FREEALL free_vector(xi,1,n);free_vector(h,1,n);free_vector(g,1,n);

void frprmn(float p[], unsigned long n, float ftol, unsigned long *iter, float *fret,
	float (*func)(float []), void (*dfunc)(float [], float []))
{
	void linmin(float p[], float xi[], unsigned long n, float *fret, float (*func)(float []));
	unsigned long j,its;
	float gg,gam,fp,dgg;
	float *g,*h,*xi;

	//added by MLJ
	FILE *fp1; unsigned long i;

	g=vector(1,n);
	h=vector(1,n);
	xi=vector(1,n);
	fp=(*func)(p);

	(*dfunc)(p,xi);
	for (j=1;j<=n;j++) {
		g[j] = -xi[j];
		xi[j]=h[j]=g[j];
	}
			

	//getchar();

	for (its=1;its<=ITMAX;its++) {
		*iter=its;
		
	//added by MLJ for debugging
	printf("return 1: fret, fp, ftol, EPS: %e, %e, %e, %e\n", *fret, fp, ftol, EPS);
	fp1=fopen("grad_test","wb");
	for(i=1; i<=NIMG; i++) { fwrite(&(xi[i]), sizeof(float), 1, fp1); }	
	fclose(fp1);

	fp1=fopen("p_test","wb");
	for(i=1; i<=NIMG; i++) { fwrite(&(p[i]), sizeof(float), 1, fp1); }	
	fclose(fp1);
	
		linmin(p,xi,n,fret,func);

	//added by MLJ for debugging
	printf("return 2: fret, fp, ftol, EPS: %e, %e, %e, %e\n", *fret, fp, ftol, EPS);

		if (2.0*fabs(*fret-fp) <= ftol*(fabs(*fret)+fabs(fp)+EPS)) {
			FREEALL
			//added by MLJ for debugging
			printf("return 1: fret, fp, ftol, EPS: %e, %e, %e, %e\n", *fret, fp, ftol, EPS);
			return;
		}
		fp= *fret;
		(*dfunc)(p,xi);
		dgg=gg=0.0;
		for (j=1;j<=n;j++) {
			gg += g[j]*g[j];
			dgg += (xi[j]+g[j])*xi[j];
		}
		if (gg == 0.0) {
			FREEALL
			return;
		}
		gam=dgg/gg;
		for (j=1;j<=n;j++) {
			g[j] = -xi[j];
			xi[j]=h[j]=g[j]+gam*h[j];
		}
	}
	nrerror("Too many iterations in frprmn");
}
#undef ITMAX
#undef EPS
#undef FREEALL
#undef NRANSI
