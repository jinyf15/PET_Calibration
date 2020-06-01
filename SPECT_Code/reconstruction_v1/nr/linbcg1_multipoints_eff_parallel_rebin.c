#include <stdio.h>
#include <math.h>
#define NRANSI
#include "nrutil.h"
#define EPS 1.0e-14
#include <sys/timeb.h>
#include <time.h>
//#include "./snrm1_multipoints.c"

void linbcg1_multipoints_eff_parallel_rebin(unsigned long n, unsigned long m, int NROW, int NCOL, float b[], float x[], int itol, float tol,
	int itmax, int *iter, float *err)
{
	void atimes_multipoints_diag_rebin(unsigned long mm, int NROW, int NCOL, float x[], float rr[], float diag[], int itrnsp);
	void atimes_multipoints_2_rebin(unsigned long mm, int NROW, int NCOL,float x1[], float x2[], float rr1[], float rr2[], int itrnsp);
	void atimes_multipoints_2_parallel_rebin(unsigned long mm,int NROW, int NCOL, float x1[], float x2[], float rr1[], float rr2[], int itrnsp);
	float snrm1_multipoints(unsigned long n, float sx[], int itol);
	unsigned long j, i;
	float ak,akden,bk,bkden,bknum,bnrm,dxnrm,xnrm,zm1nrm,znrm;
	float *p,*pp,*r,*rr,*z,*zz, *diag, *diag_tmp;
    char filename[500], filename2[500], filename1[500];
	FILE *fp1, *fp2;
	int na, nos, nsubdet, nf_current, nf, if_master, duration, time_period, if_period;
	int nt1, nt2, nt3, nt4, nt5, nt6, nt7, nn, ntt6, ntt7, ww_ready, rr_ready;
	float nt8;
	extern int niter;
	struct _timeb tstruct;

	duration=600; // (second)
	time_period=10; // (milisecond) 
	if_period=1;

	p=vector(1,n);
	pp=vector(1,n);
	r=vector(1,n);
	rr=vector(1,n);
	z=vector(1,n);
	zz=vector(1,n);
	diag=vector(1,n);
	diag_tmp=vector(1, m);


////////////////////////////////////////////////////////////////////////////////
	/** calculate file_read_order (start) **/
	sprintf(filename2, "%s\\data\\file_read_order_UCRB", ROOT_DIR);
	nf=0;
	if((fp2=fopen(filename2, "rb"))==NULL) {
		fp2=fopen(filename2, "wb");
		for (nos=1; nos<=NOS; nos++) {
			for(na=(nos-1); na<NANG; na=na+NOS) { //no. of subset
				for (nsubdet=1; nsubdet<=NSUBDET; nsubdet++) {
					nf++;
					fwrite(&(nf), sizeof(int), 1, fp2);	
					fwrite(&(nos), sizeof(int), 1, fp2);
					fwrite(&(na), sizeof(int), 1, fp2);
					fwrite(&(nsubdet), sizeof(int), 1, fp2);
				}
			}
		}
		fclose(fp2);
		printf("file_read_order is ready!\n");
	}
	else {fclose(fp2);}
	/** calculate file_read_order (end) **/

/////////////////////////////////////////////////////////////////////////*/
	*iter=0;
	
	sprintf(filename2, "%s\\data\\r_UCRB", ROOT_DIR);
	if((fp2=fopen(filename2, "rb"))==NULL) {
		atimes_multipoints_diag_rebin(n,NROW, NCOL, x,r, diag_tmp, 1);
		fp2=fopen(filename2, "wb");
		for (i=1; i<=n; i++) fwrite(&(r[i]), sizeof(float), 1, fp2);
		for (i=1; i<=m; i++) fwrite(&(diag_tmp[i]), sizeof(float), 1, fp2);
		fclose (fp2);
	}
	else {
		for (i=1; i<=n; i++) fread(&(r[i]), sizeof(float), 1, fp2);
		for (i=1; i<=m; i++) fread(&(diag_tmp[i]), sizeof(float), 1, fp2);
		fclose (fp2);
	}


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
		//printf("*********** before atimes \n n=%d, bk=%f, bknum=%f\n", *iter, bk, bknum); getchar();
//		atimes_multipoints_2(n,p,pp,z,zz,1);

/**************************************************/
		/*** save bknum for current iteration (start) ***/
		sprintf(filename2, "%s\\data\\bknum", ROOT_DIR);
		if (if_period) {
			while (1){
				_ftime( &tstruct );
				if (tstruct.millitm/time_period*time_period==tstruct.millitm-time_period/2){
					break;
				}
			}
		}
		fp2=fopen(filename2,"rb");
		if (fp2!=NULL){
			fread(&(nt7), sizeof(int), 1, fp2);	
			fread(&(nt8), sizeof(float), 1, fp2);	
			fclose(fp2);
			if ((nt7)==(*iter-1)) {
				//printf("start to save latest bknum\n");
				if (if_period) {
					while (1){
						ftime( &tstruct );
						if (tstruct.millitm/time_period*time_period==tstruct.millitm){
							break;
						}
					}
				}

				fp2=fopen(filename2, "wb");
				if(fp2!=NULL ) {
					fwrite(&(*iter), sizeof(int), 1, fp2);	
					fwrite(&(bknum), sizeof(float), 1, fp2);
				}
				fclose(fp2);

				/*** save vectors x, r, rr (start) ***/
				sprintf(filename2, "%s\\data\\x_r_rr_%d", ROOT_DIR, *iter);
				if (if_period) {
					while (1){
						_ftime( &tstruct );
						if (tstruct.millitm/time_period*time_period==tstruct.millitm){
							break;
						}
					}
				}
				fp2=fopen(filename2,"wb");
				for(i=1; i<=n; i++) fwrite(&(x[i]), sizeof(float), 1, fp2);
				for(i=1; i<=n; i++) fwrite(&(r[i]), sizeof(float), 1, fp2);
				for(i=1; i<=n; i++) fwrite(&(rr[i]), sizeof(float), 1, fp2);
				fclose(fp2);
				/*** save vectors x, r, rr (end) ***/
			}
		}
		else {
			//printf("start to save bknum\n");
			if (if_period) {
				while (1){
					ftime( &tstruct );
					if (tstruct.millitm/time_period*time_period==tstruct.millitm){
						break;
					}
				}
			}

			fp2=fopen(filename2, "wb");
			if(fp2!=NULL ) {
				fwrite(&(*iter), sizeof(int), 1, fp2); 
				fwrite(&(bknum), sizeof(float), 1, fp2); 
			}
			fclose(fp2);

			/*** save vectors x, r, rr (start) ***/
			sprintf(filename2, "%s\\data\\x_r_rr_%d", ROOT_DIR, *iter);
			if (if_period) {
				while (1){
					_ftime( &tstruct );
					if (tstruct.millitm/time_period*time_period==tstruct.millitm){
						break;
					}
				}
			}
			fp2=fopen(filename2,"wb");
			for(i=1; i<=n; i++) fwrite(&(x[i]), sizeof(float), 1, fp2);
			for(i=1; i<=n; i++) fwrite(&(r[i]), sizeof(float), 1, fp2);
			for(i=1; i<=n; i++) fwrite(&(rr[i]), sizeof(float), 1, fp2);
			fclose(fp2);
			/*** save vectors x, r, rr (end) ***/
		}
		/*** save bknum for current iteration (end) ***/

		//printf("in!!!!!!!!!!!1\n");
		atimes_multipoints_2_parallel_rebin(n,NROW, NCOL, p,pp,z,zz,1, *iter);
		*iter=niter;
		//printf("out!!!!!!!1\n");getchar();

		/*** recall vector bknum (start) ***/
		sprintf(filename2, "%s\\data\\bknum", ROOT_DIR);
		if (if_period) {
			while (1){
				_ftime( &tstruct );
				if (tstruct.millitm/time_period*time_period==tstruct.millitm-time_period/2){
					break;
				}
			}
		}
		fp2=fopen(filename2,"rb");
		if (fp2!=NULL){
			fread(&(nt7), sizeof(int), 1, fp2);	
			fread(&(nt8), sizeof(float), 1, fp2);
			fclose(fp2);
			if (nt7==(*iter)) {
				bknum=nt8;
			}
		}
		bkden=bknum;
		/*** recall vector bknum (end) ***/

		/*** recall vector x, r, rr (start) ***/
		while (1) {
			sprintf(filename2, "%s\\data\\x_r_rr_%d", ROOT_DIR, niter);
			if (if_period) {
				while (1){
					_ftime( &tstruct );
					if (tstruct.millitm/time_period*time_period==tstruct.millitm-time_period/2){
						break;
					}
				}
			}
			fp2=fopen(filename2,"rb");
			if (fp2!=NULL){
				for(i=1; i<=n; i++) fread(&(x[i]), sizeof(float), 1, fp2);
				for(i=1; i<=n; i++) fread(&(r[i]), sizeof(float), 1, fp2);	
				for(i=1; i<=n; i++) fread(&(rr[i]), sizeof(float), 1, fp2);	
				fclose(fp2);
				break;
			}
			else {
				for(nn=1;nn<3e7; nn++)
					;
				printf("."); 
			}
		}

		sprintf(filename2, "%s\\data\\x_r_rr_%d", ROOT_DIR, niter);
		nn=0; while ((remove(filename2) == -1 ) && ((nn++)<100)) printf("delete error: %s\n", filename2);
		/*** recall vector x, r, rr (end) ***/

		/*** recall vector p and pp (start) ***/
		while (1) {
			sprintf(filename2, "%s\\data\\x12_ni%d", ROOT_DIR, niter);
			if (if_period) {
				while (1){
					_ftime( &tstruct );
					if (tstruct.millitm/time_period*time_period==tstruct.millitm-time_period/2){
						break;
					}
				}
			}
			fp2=fopen(filename2,"rb");
			if (fp2!=NULL){
				for(i=1; i<=n; i++) fread(&(p[i]), sizeof(float), 1, fp2);
				for(i=1; i<=n; i++) fread(&(pp[i]), sizeof(float), 1, fp2);	
				fclose(fp2);
				break;
			}
			else {
				for(nn=1;nn<3e7; nn++)
					;
				printf("."); 
			}
		}

		sprintf(filename2, "%s\\data\\x12_ni%d", ROOT_DIR, niter);
		nn=0; while ((remove(filename2) == -1 ) && ((nn++)<100)) printf("delete error: %s\n", filename2);

		sprintf(filename1, "%s\\data\\x12_file_ready_UCRB", ROOT_DIR);
		if (if_period) {
			while (1){
				_ftime( &tstruct );
				if (tstruct.millitm/time_period*time_period==tstruct.millitm-time_period/2){
					break;
				}
			}
		}

		fp2=fopen(filename1, "rb");
		if (fp2!=NULL) {
			fread(&(ntt6), sizeof(int), 1, fp2);
			fread(&(rr_ready), sizeof(int), 1, fp2);
			fread(&(ww_ready), sizeof(int), 1, fp2);
			fread(&(ntt7), sizeof(int), 1, fp2);
			fclose(fp2);
		}

		if ((ntt6==(niter))&&(rr_ready==1)&&(ww_ready==0)&&(ntt7==(niter))) {
			if (if_period) {
				while (1){
					_ftime( &tstruct );
					if (tstruct.millitm/time_period*time_period==tstruct.millitm){
						break;
					}
				}
			}
			fp2=fopen(filename1, "wb");
			if (fp2!=NULL) {
				ntt6=niter; fwrite(&(ntt6), sizeof(int), 1, fp2);
				rr_ready=0; fwrite(&(rr_ready), sizeof(int), 1, fp2);
				ww_ready=1; fwrite(&(ww_ready), sizeof(int), 1, fp2);
				ntt7=niter; fwrite(&(ntt6), sizeof(int), 1, fp2);
				fclose(fp2);
			}
		}
		else {printf("error!!");}
		
		/*** recall vector p and pp (end) ***/

/***************************************************/

		for (akden=0.0,j=1;j<=n;j++) akden += z[j]*pp[j];
		ak=bknum/akden;
		//printf("*********** after atimes \n n=%d, ak=%f, bknum=%f\n", *iter, ak, bknum); getchar();

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
		printf("iter=%4d err=%12.6f\n",*iter,*err);
		printf("************************************************\n");

		if ((*iter/1*1)==*iter){
		sprintf(filename, "%s\\results\\x_%d", ROOT_DIR, *iter);
		fp1=fopen(filename, "wb"); 
		for (j=1; j<=n; j++) fwrite(&(x[j]), sizeof(float), 1, fp1);			
		fclose(fp1);
		}

		sprintf(filename, "%s\\results\\err_%d", ROOT_DIR, *iter);
		fp1=fopen(filename, "wb"); 
		fwrite(err, sizeof(float), 1, fp1);			
		fclose(fp1);

	if (*err <= tol) break;
	}

	free_vector(p,1,n);
	free_vector(pp,1,n);
	free_vector(r,1,n);
	free_vector(rr,1,n);
	free_vector(z,1,n);
	free_vector(zz,1,n);
	free_vector(diag, 1, n);
	free_vector(diag_tmp, 1, m);
}
#undef EPS
#undef NRANSI
