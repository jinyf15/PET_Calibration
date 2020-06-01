#include <stdio.h>
#include <math.h>
#define NRANSI
#include "nrutil.h"
#define EPS 1.0e-14
#include "./snrm1_multipoints.c"
void linbcg1_multipoints(unsigned long n, float b[], float x[], int itol, float tol,
	int itmax, int *iter, float *err)
{
	void asolve_multipoints(unsigned long n, float b[], float x[], int itrnsp);
	void atimes_multipoints(unsigned long n, float x[], float r[], int itrnsp);
	float snrm1_multipoints(unsigned long n, float sx[], int itol);
	unsigned long j;
	float ak,akden,bk,bkden,bknum,bnrm,dxnrm,xnrm,zm1nrm,znrm;
	float *p,*pp,*r,*rr,*z,*zz;
	char filename[500];
	FILE *fp1;

	p=vector(1,n);
	pp=vector(1,n);
	r=vector(1,n);
	rr=vector(1,n);
	z=vector(1,n);
	zz=vector(1,n);


	*iter=0;
	atimes_multipoints(n,x,r,1);

	for (j=1;j<=n;j++) {
		r[j]=b[j]-r[j];
		rr[j]=r[j];
	}
	if (itol == 1) {
		bnrm=snrm1_multipoints(n,b,itol);
		asolve_multipoints(n,r,z,0);
	}
	else if (itol == 2) {
		asolve_multipoints(n,b,z,0);
		bnrm=snrm1_multipoints(n,z,itol);
		asolve_multipoints(n,r,z,0);
	}
	else if (itol == 3 || itol == 4) {
		asolve_multipoints(n,b,z,0);
		bnrm=snrm1_multipoints(n,z,itol);
		asolve_multipoints(n,r,z,0);
		znrm=snrm1_multipoints(n,z,itol);
	} else nrerror("illegal itol in linbcg");
	while (*iter <= itmax) {
		++(*iter);
		asolve_multipoints(n,rr,zz,1);
		for (bknum=0.0,j=1;j<=n;j++) bknum += z[j]*rr[j];
		if (*iter == 1) {
			for (j=1;j<=n;j++) {
				p[j]=z[j];
				pp[j]=zz[j];
			}
		}
		else {
			bk=bknum/bkden;
			for (j=1;j<=n;j++) {
				p[j]=bk*p[j]+z[j];
				pp[j]=bk*pp[j]+zz[j];
			}
		}
		bkden=bknum;
		atimes_multipoints(n,p,z,1);
		for (akden=0.0,j=1;j<=n;j++) akden += z[j]*pp[j];
		ak=bknum/akden;
		atimes_multipoints(n,pp,zz,1);
		for (j=1;j<=n;j++) {
			x[j] += ak*p[j];
			r[j] -= ak*z[j];
			rr[j] -= ak*zz[j];
		}
		asolve_multipoints(n,r,z,0);
		if (itol == 1)
			*err=snrm1_multipoints(n,r,itol)/bnrm;
 		else if (itol == 2)
			*err=snrm1_multipoints(n,z,itol)/bnrm;
		else if (itol == 3 || itol == 4) {
			zm1nrm=znrm;
			znrm=snrm1_multipoints(n,z,itol);
			if (fabs(zm1nrm-znrm) > EPS*znrm) {
				dxnrm=fabs(ak)*snrm1_multipoints(n,p,itol);
				*err=znrm/fabs(zm1nrm-znrm)*dxnrm;
			} else {
				*err=znrm/bnrm;
				continue;
			}
			xnrm=snrm1_multipoints(n,x,itol);
			if (*err <= 0.5*xnrm) *err /= xnrm;
			else {
				*err=znrm/bnrm;
				continue;
			}
		}
		printf("iter=%4d err=%12.6f\n",*iter,*err);
		
		sprintf(filename, ".\\results\\x_%d", *iter);
		fp1=fopen(filename, "wb"); 
		for (j=1; j<=n; j++) fwrite(&(x[j]), sizeof(float), 1, fp1);			
		fclose(fp1);



	if (*err <= tol) break;
	}

	free_vector(p,1,n);
	free_vector(pp,1,n);
	free_vector(r,1,n);
	free_vector(rr,1,n);
	free_vector(z,1,n);
	free_vector(zz,1,n);
}
#undef EPS
#undef NRANSI
