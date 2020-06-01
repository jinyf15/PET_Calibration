#include <stdio.h>
#include <math.h>
#define NRANSI
#include "nrutil.h"
#define EPS 1.0e-14


void linbcg1_multipoints_eff_rebin2_filter(unsigned long n, unsigned long m, int NROW, int NCOL, double b[], double x[], int itol, float tol,
	int itmax, int *iter, float *err, int fwhm, int pt)
{
	void atimes_multipoints_diag_rebin(unsigned long mm, int NROW, int NCOL, float x[], float rr[], float diag[], int itrnsp);
	void atimes_multipoints_2_rebin(unsigned long mm, int NROW, int NCOL,  float x1[], float x2[], float rr1[], float rr2[], int itrnsp);
	float snrm1_multipoints(unsigned long n, float sx[], int itol);
	unsigned long j, i;
	float ak,akden,bk,bkden,bknum,bnrm,dxnrm,xnrm,zm1nrm,znrm;
	double *p,*pp,*r,*rr,*z,*zz, *diag, *diag_tmp, *g;
    char filename[1500];
	FILE *fp1;

	p=dvector(1,n);
	pp=dvector(1,n);
	r=dvector(1,n);
	rr=dvector(1,n);
	g=dvector(1,n);
	z=dvector(1,n);
	zz=dvector(1,n);
	diag=dvector(1,n);
	diag_tmp=dvector(1, m);

	*iter=0;
	atimes_multipoints_diag_rebin2(n,NROW, NCOL, x,r, diag_tmp, 1);

	sprintf(filename, "%s\\results\\r_ucrb", ROOT_DIR);
	fp1=fopen(filename, "wb"); 
	for (i=1; i<=n; i++) fwrite(&(r[i]), sizeof(double), 1, fp1);
	for (i=1; i<=m; i++) fwrite(&(diag_tmp[i]), sizeof(double), 1, fp1);
	fclose(fp1);

	for (j=1; j<=n/m; j++){
		for (i=1; i<=m; i++){
			diag[i+(j-1)*m]=diag_tmp[i];
		}
	}

	for (j=1;j<=n;j++) {
		r[j]=b[j]-r[j];
		rr[j]=r[j];
	}
	if (itol == 1) {
		bnrm=snrm1_multipoints(n,b,itol);
		for (i=1; i<=n; i++) z[i]=(diag[i] != 0.0 ? r[i]/diag[i] : r[i]);
	}
	else if (itol == 2) {
		for (i=1; i<=n; i++) z[i]=(diag[i] != 0.0 ? b[i]/diag[i] : b[i]);
		bnrm=snrm1_multipoints(n,z,itol);
		for (i=1; i<=n; i++) z[i]=(diag[i] != 0.0 ? r[i]/diag[i] : r[i]);
	}
	else if (itol == 3 || itol == 4) {
		for (i=1; i<=n; i++) z[i]=(diag[i] != 0.0 ? b[i]/diag[i] : b[i]);
		bnrm=snrm1_multipoints(n,z,itol);
		for (i=1; i<=n; i++) z[i]=(diag[i] != 0.0 ? r[i]/diag[i] : r[i]);
		znrm=snrm1_multipoints(n,z,itol);
	} else nrerror("illegal itol in linbcg");
	while (*iter <= itmax) {
		++(*iter);
		for (i=1; i<=n; i++) zz[i]=(diag[i] != 0.0 ? rr[i]/diag[i] : rr[i]);
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
		atimes_multipoints_2_rebin2(n,NROW, NCOL, p,pp,z,zz,1);

		for (akden=0.0,j=1;j<=n;j++) akden += z[j]*pp[j];
		ak=bknum/akden;
		
		for (j=1;j<=n;j++) {
			x[j] += ak*p[j];
			r[j] -= ak*z[j];
			rr[j] -= ak*zz[j];
		}
		for (i=1; i<=n; i++) z[i]=(diag[i] != 0.0 ? r[i]/diag[i] : r[i]);
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
		printf("************************************************\n");
		printf("pt=%d fwhm=%d iter=%4d err=%12.6f\n", pt, fwhm, *iter,*err);
		printf("************************************************\n");

		if ((*iter/100*100)==*iter){
			for (j=1; j<=n; j++) g[j]=0;
			atimes_multipoints_2_rebin2(n, NROW, NCOL, x, x, g, g, 0);

			sprintf(filename, "%s\\results\\xxr%d_%d_%d", ROOT_DIR, pt, fwhm, *iter);
			fp1=fopen(filename, "wb"); 
			for (j=1; j<=n; j++) fwrite(&(x[j]), sizeof(double), 1, fp1);			
			fclose(fp1);

			sprintf(filename, "%s\\results\\ggr%d_%d_%d", ROOT_DIR, pt, fwhm, *iter);
			fp1=fopen(filename, "wb"); 
			for (j=1; j<=n; j++) fwrite(&(g[j]), sizeof(double), 1, fp1);			
			fclose(fp1);
		}

		sprintf(filename, "%s\\results\\error\\err%d_%d_%d", ROOT_DIR, pt, fwhm, *iter);
		fp1=fopen(filename, "wb"); 
		fwrite(err, sizeof(float), 1, fp1);			
		fclose(fp1);

	if (*err <= tol) break;
	}

	free_dvector(p,1,n);
	free_dvector(pp,1,n);
	free_dvector(r,1,n);
	free_dvector(rr,1,n);
	free_dvector(g,1,n);
	free_dvector(z,1,n);
	free_dvector(zz,1,n);
	free_dvector(diag, 1, n);
	free_dvector(diag_tmp, 1, m);
}
#undef EPS
#undef NRANSI
