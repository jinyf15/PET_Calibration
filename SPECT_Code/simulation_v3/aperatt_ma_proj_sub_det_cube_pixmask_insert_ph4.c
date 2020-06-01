/* code description (start) */
/*
The code is used to simulate response functions in 3D geometry, which is suitable for the condition that aperture plane and detector plane are not parallel. 
In the code the penetration factor is considered. N1, N2, N3, N4 is the step number for penetration calculation.
new int value isave is used to determine save pattern of response function. isave=1, rfs is saved in current disk, isave=2, rfs is saved in all disks in local pc, isave=3, rfs is saved in the disk determined by the file: path.txt in COMMON_DIR\\data.
The necessary parameters is transfered by file of dat_a_MPH.txt. The number of pinhole will be calculated automatically and show on the screen.
In the meantime, the two intermediate variables (i and ii) will also show. i==ii means calculation is not wrong.
Take care: 	hscMAX and hscMIN should be changed for different geometry and don't forget to change.
*/
/* code description (end) */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <malloc.h>
#include <time.h>

// general 
#include "./nr/nr.h"
#include "./nr/nrutil.h"
//#include "./setting1.h"
#include "./setting2.h"
#include "./setting3.h"

// settings for pinhole geometry
#define PI 3.14159265458979
#define NPH_DET1 500
#define MAX_DM_DET 10

//global variables and sub-routines
float *sen, *sp, *dp, *source, *img, *proj,*det_mask, *det_pix_mask; // det_mask: mark all useful pixels on the detector; det_pix_mask: mark detector pixels satisfy certain creterion	
float pixelprob_det1(float *buf, float *dpos);
float pixelprob_det1NL(float *buf, float *dpos, float *spos, float *map, float *dirt, int *index );
float pixelprob_det2ML_sum(float *buf, float *dpos, float *spos, float *map, float *dirt, float *cos_D, int DM_DET, float mu_det );
float pixelprob_det2ML_sum_ph_insert(float *buf, float *dpos, float *spos, float *map, float *dirt, float *cos_D, int DM_DET, float mu_det );
float penetration_map_layer(float r1, float r2, float r3, float ct, float cht,float mu, float *map,float *dirt , float tanAlpha);
float penetration_map_layer_insert_ph(float r1, float r2, float r3, float ct, float cht, float pct, float mu1, float mu2, float map[],float *dirt, float tanAlpha, float drill_angle);
float penetration_map_layer_symetric(float r1, float r2, float r3, float ct, float cht, float pct, float mu1, float mu2, float map[],float *dirt, float tanAlpha, float drill_angle);
void save_path();
float pixelprob_det2ML_sum_ph_insert3(float *buf, float *dpos, float *subdpos, float *spos, float *map, float *dirt, float *cos_D, int DM_DET, float mu_det, int fn);
void penetration_map_compound_eye(float *,float *, float *,float *,float *, float *,float *,float *, float *,float *,int,int,float*,float);
										
int aperatt_ma_proj_sub_det_cube_pixmask_insert_ph3(int *iistart_out, int *nph_out)
{
	extern long idum;
	extern int  DET_MASK_CENTERX1, DET_MASK_CENTERX2, DET_MASK_CENTERY1, DET_MASK_CENTERY2, DET_MASK_RADIUS1, DET_MASK_RADIUS2,NMAX1, NMAX2, NMAXr1, NMAXr2, DM_DET1, DM_DET2, NMAX, NMAXr;
	extern float PH_RADIUS1, PH_RADIUS2, PH_CHT1, PH_CHT2, PH_CT1, PH_CT2, LEAD_CHT1, LEAD_CHT2, MU_APER1, MU_APER2, PB_APER1, PB_APER2, MU_DET1, MU_DET2, DDZ_DET1, DDZ_DET2, APERALPHA1, APERALPHA2, APERALPHA_LIMIT1, APERALPHA_LIMIT2;
	extern int NANG, NSX, NSY, NSZ, NIMG, DN_DET1, NDET_DET1, NSUBDET, NSIZE, NCPU, NOS, NOS, NCOPY, IDETPOS,IMAP;
	extern float DANG, SDX, SDY, SDZ, DDX_DET1, DDY_DET1;
	extern char CALIB_FILE_DIR[256];

	long int nn, na, nt;
	// variables for defining the aperture 
	float sx, sy, sz, scx, scy, scz, sdx, sdy, sdz, thresh, czz; //source
	float da; //detector pixel
	float cx, cy, ct,cht; //collimator slab
	float *aaa;
	float cz[NPH_DET1]={0.},  hx[NPH_DET1]={0.},  hy[NPH_DET1]={0.}; //center and radius // === CHANGE ===//
	float hz1, r1, 
		  hz2, r2,r4,r5,
		  hz3, r3;
	float dx_fg[NPH_DET1]={0.},  dy_fg[NPH_DET1]={0.}; // geng fu 
	float bhx[NPH_DET1][2], bhy[NPH_DET1][2];

	// other variables 
	int p, q, nx, ny, nsx, nsy, nsz, ndx, ndy, nos, istart, ifinished,  iistart, nf_current,nt3,nt4,nt5;
	int sn_fg;
	long int i, j, k, sn, dn, dncx, dncy, *dn_buf, l, n, ii, in, nf, iii;
	float mu, mu_pb, prob, deltaOmega, tanTheta, *proj_buf, tmp, tmp1, tmp2, tmp3, tmp4,xx,yy;
	float r, rmax, t1, t2, t3;
	float d, dcx, dcy, dcz, dx, dx0,dy,dy0, dz,dz0, ddx, ddy, ddz,hdist, alpha; //detector pixel
	float *compx,*compy; // center position for compound eye;
	float *radius_cp;// radius of compound eye;
	short int *buf_hn, nph;
	int	dpi;
	char filename[1000], filename1[1000], filename2[1000], filename3[1000], srfs_path[1000];
	FILE *fp1, *fp2, *fp3, *fp4, *fp5,*fp_fg;
	float f, source_max, *sat;
	extern float *sen, *sp, *dp, *source, *img, *proj;	
	float deltan, dpcx, dpcy, rp, dsq, pp, pr, prt, sen_max, pr_max, pr_tmp[2000], dn_tmp[2000], tx, ty; 

	float cosTheta;
	float hs, hsp,hsa, hsc[NPH_DET1]={0.};			// distance from source to pinhole and detector plane
	float esaL,elaL,esap,elap;
	float buf[15], dpos[3*MAX_DM_DET],spos[3],apos[3];		// parameter definition for probability calculation subroutine
	float proj_cx,proj_cy,dpy1,dpy2,dpx1,dpx2;
	
	/** variable definition for detector plane orientation  **/
	float cosAlpha,sinAlpha; // rotation angle of detector plate between axis and y&z surface
	float cosPhi,sinPhi;	// rotation angle of detector plate between  axis and x&z surface
	float cosBeta,sinBeta;	// rotation angle of pixel on detector plate
	float cosThetaD, sinThetaD;		//incident angle to detector surface

	/** variable definition for aperture plane rotation  **/
    float AperSinAlpha[NPH_DET1]={0.};			// definition of aperture rotation angle
	float AperSinPhi[NPH_DET1]={0.};
	float AperCosAlpha[NPH_DET1]={1.};
	float AperCosPhi[NPH_DET1]={1.};

	/** vector definition  **/
	float detX[3]={0.},detY[3]={0.};		// X and Y pixel direction in detector
	float detAXx[NPH_DET1]={1.},detAXy[NPH_DET1]={0.},detAXz[NPH_DET1]={0.};
	float detAYx[NPH_DET1]={0.},detAYy[NPH_DET1]={1.},detAYz[NPH_DET1]={0.},unit_detAY;	// X and Y pixel direction of detector on Aperture surface
	float detAZx[NPH_DET1]={0.},detAZy[NPH_DET1]={0.},detAZz[NPH_DET1]={1.};
	float detZ[3]={0.};						// n vector of detector plane
	float nsrx,nsry,nsrz;							// n vector of source&pinhole
	float esa[3]={0.},ela[3]={0.};			// n vector of long and short axis of ellipse
	float DDX_ph[NPH_DET1]={0.},DDY_ph[NPH_DET1]={0.};

	/** variables for penetration factor  **/
	float rho, theta, h, mm1, mm2, mm3, mm4, mm5, tt1, tt2, tt3, tt4, dirt_L, *map[NPH_DET1], Rlimit[NPH_DET1],rlimit, Beta[NPH_DET1], treshR, tanAlpha;
	float alphaA, betaA, gammaA, alphaD, betaD, gammaD, hscMAX, hscMIN, AperAlpha, tanBeta[NPH_DET1];
	int  m,index[2], ph_order[NPH_DET1]={0}, tt, nm;
	float cos_D[MAX_DM_DET],dirt[4],*det_prob, ll, nxtmp[2], nytmp[2], nztmp[2], nxx, nyy, nzz, nx_pb, ny_pb, nz_pb;

	/** system response function **/
	long int npix, spri, sprj, sprk, *ijat, isave, totalcounts, msize, *sequence, isequence, dtmpx, dtmpy;
    time_t time1;
	float *dsubp;

	/** temp variables **/
	long int sn_evt, dn_evt, dsubn, nevt, ntmp, ntmp1, ntmp2, ntmp3, ntmp4, ntmp5, fn, ntmp6,ntmp7,ntmp8,ntmp9, nsubdet, nt1, nt2, nna, nnsubdet;
	float *sou_mask, minhsc, totalpixel, pixstep, dp1, dp2, *dxx0, *dyy0, dpx, dpy;

	int DM_DET,DET_MASK_CENTERX,DET_MASK_CENTERY,DET_MASK_RADIUS;
	float PH_RADIUS,PH_CHT,PH_CT,LEAD_CHT,MU_APER,PB_APER,MU_DET,DDZ_DET,APERALPHA,APERALPHA_LIMIT;
	float *subp, *subdpos, tmpdx, tmpdy, tmpdz; 

	/** =1: random srf locations, single copy only; =2: complete copy on each HD; isave=2, rfs is saved in all disks in local pc; isave=3, rfs is saved in the disk determined by the path.txt **/
	isave=3;isequence=1;

	/** initializatons **/
	det_prob=vector(0,NDET_DET1-1);
	det_mask=vector(0, NDET_DET1-1);
	det_pix_mask=vector(0, NDET_DET1-1);
	proj=vector(0, NDET_DET1-1);
	proj_buf=vector(0, NDET_DET1-1);
	dn_buf=lvector(0, NDET_DET1-1);
	sou_mask=vector(0, NIMG-1);
	sen=vector(0, NIMG-1);
	source=vector(0, NIMG-1);
	img=vector(0, NIMG-1);
	sp=vector(0, NIMG*3-1);
	dxx0=vector(0, NDET_DET1-1);
	dyy0=vector(0, NDET_DET1-1);
	compx=vector(0,NPH_DET1); // center position for each compound eye;
	compy=vector(0,NPH_DET1);
	radius_cp=vector(0,NPH_DET1);
	aaa=vector(0,NPH_DET1*10+16);
	fn=4;
	sequence=lvector(1, (NANG*NSUBDET));

	for(i=0; i<=NDET_DET1-1; i++) { det_mask[i]=0;det_pix_mask[i]=0;det_prob[i]=0.;}
	for(i=0; i<=NIMG-1; i++) { sou_mask[i]=0;}
	totalcounts=0;
	/** read in aperture and detector information (start) **/
	if(0) {
		nph=31;
		sprintf(filename1, "%s\\data\\aaa", COMMON_DIR);
		fp1=fopen(filename1, "rb"); 
		printf("using %s as the system geometry\n", filename1);
		for(i=1; i<=14+6*nph; i++) {fread(&(aaa[i]), sizeof(float), 1, fp1); printf("i=%d %e\n", i, aaa[i]);}
		printf("a =%f\n", aaa[197]); getchar();
		fclose(fp1);
	}
	if(1) {
		sprintf(filename1, "%s", CALIB_FILE_DIR);
		fp1=fopen(filename1, "r"); 
		printf("using %s as the system geometry with insert ph\n", filename1);
		i=0;
		while (!(feof(fp1))) {
			i++;
			fscanf(fp1, "%f", &(aaa[i]));printf("%e %d\n", aaa[i], feof(fp1));
			//if((i-11)%10==0&&(i-11)/10==150) {printf("%d nph \n",(i-11)/10);getchar();}
		}
		nph=(i-16)/10;
		printf("i =%d nph=%d ii=%d\n", i, nph, nph*10+16); //getchar();
		fclose(fp1);
		*nph_out=nph;
	}
	/** read in aperture and detector information (end) **/
	if(1) {
		sprintf(filename1, "%s\\data\\compound_eye.txt", COMMON_DIR );
		fp1=fopen(filename1, "r"); 
		printf("using %s as the compound eye parameter\n", filename1);
		i=0;
		if(fp1==NULL) 
		{
			printf(" error opening file");
			getchar();
		}
		else
		{
			while (!(feof(fp1))) {
				i++;
				if(i%3==1) {fscanf(fp1, "%f", &(compx[i/3]));printf("%e %d\n", compx[i/3], feof(fp1));}
				if(i%3==2) {fscanf(fp1, "%f", &(compy[i/3]));printf("%e %d\n", compy[i/3], feof(fp1));}
				if(i%3==0) {fscanf(fp1, "%f", &(radius_cp[i/3-1]));printf("%e %d\n", radius_cp[i/3], feof(fp1));}

			}
		}
		}
		printf("i=%d,nph=%d\n",i,nph);


	if (nph==192) {
		DM_DET=DM_DET1;
		NMAX=NMAX1;
		NMAXr=NMAXr1;
		PH_RADIUS=PH_RADIUS1;
		DET_MASK_CENTERX=DET_MASK_CENTERX1;
		DET_MASK_CENTERY=DET_MASK_CENTERY1;
		DET_MASK_RADIUS=DET_MASK_RADIUS1;
		PH_CHT=PH_CHT1;
		PH_CT=PH_CT1;
		LEAD_CHT=LEAD_CHT1;
		MU_APER=MU_APER1;
		PB_APER=PB_APER1;
		MU_DET=MU_DET1;
		DDZ_DET=DDZ_DET1;
		APERALPHA=APERALPHA1;
		APERALPHA_LIMIT=APERALPHA_LIMIT1;
	}
	else if (nph==7){
		DM_DET=DM_DET2;
		DET_MASK_CENTERX=DET_MASK_CENTERX2;
		DET_MASK_CENTERY=DET_MASK_CENTERY2;
		DET_MASK_RADIUS=DET_MASK_RADIUS2;
		NMAX=NMAX2;
		NMAXr=NMAXr2;
		PH_RADIUS=PH_RADIUS2;
		PH_CHT=PH_CHT2;
		PH_CT=PH_CT2;
		LEAD_CHT=LEAD_CHT2;
		MU_APER=MU_APER2;
		PB_APER=PB_APER2;
		MU_DET=MU_DET2;
		DDZ_DET=DDZ_DET2;
		APERALPHA=APERALPHA2;
		APERALPHA_LIMIT=APERALPHA_LIMIT2;
	}
	else {
		printf("there is a mistake to read ph information!");getchar();getchar();
	}
	printf("DM_DET=%d\n", DM_DET);
	printf("NMAX=%d\n", NMAX);
	printf("NMAXr=%d\n", NMAXr);
	printf("DET_MASK_CENTERX=%d\n", DET_MASK_CENTERX);
	printf("DET_MASK_CENTERY=%d\n", DET_MASK_CENTERY);
	printf("DET_MASK_RADIUS=%d\n", DET_MASK_RADIUS);
	printf("PH_RADIUS=%f\n", PH_RADIUS);
	printf("PH_CHT=%f\n", PH_CHT);
	printf("MU_APER=%f\n", MU_APER);
	printf("MU_DET=%f\n", MU_DET);
	printf("DDZ_DET=%f\n", DDZ_DET);
	printf("APERALPHA=%f\n", APERALPHA);
	printf("APERALPHA_LIMIT=%f\n", APERALPHA_LIMIT);
	printf("NSX=%d, NSY=%d, NSZ=%d\n", NSX, NSY, NSZ);
	printf("DN_DET1=%d\n", DN_DET1);
	printf("NANG=%d, NOS=%d, NSUBDET=%d\n", NANG, NOS, NSUBDET);

	subdpos=vector(0, DM_DET*3*fn*fn-1);
	dsubp=vector(0, DM_DET*NDET_DET1*3*fn*fn-1);
	dp=vector(0, NDET_DET1*DM_DET*3-1);
	sat=vector(1, NMAX);
	ijat=lvector(1, NMAX);
	for(i=1; i<=NMAX; i++) { sat[i]=0.; ijat[i]=0; }

	if (isequence) {
		sprintf(filename2, "%s\\data\\file_read_order",BUF_DIR1);
		nf=0;
		fp2=fopen(filename2, "wb");
		for(n=1; n<=1; n++) {		
			for (nos=1; nos<=NOS; nos++) {
				for(na=(nos-1); na<NANG; na=na+NOS) { //no. of subset
					for (nsubdet=1; nsubdet<=NSUBDET; nsubdet++) {
						nf++;
						fwrite(&(nf), sizeof(int), 1, fp2);	
						fwrite(&(nos), sizeof(int), 1, fp2);
						fwrite(&(na), sizeof(int), 1, fp2);
						fwrite(&(nsubdet), sizeof(int), 1, fp2);
						fwrite(&(n), sizeof(int), 1, fp2);	
					}
				}
			}
		}
		fclose(fp2);
		printf("file_read_order is ready!");
	}
    
	/**** defining basic geometry (start) ** **/
	/**  detector geometry settings **/ 
	ddx=DDX_DET1; ddy=DDY_DET1;		ddz=DDZ_DET;//detector pixel dimension 
	dcx=aaa[1]; dcy=aaa[2]; dcz=aaa[3];
	da=ddx*ddy;						//detector pixel size
	dz=dcz;							//detector z position 

	/**  detector angle setting  **/
	/* sinAlpha=sin(0.);	cosAlpha=cos(0.);	// the angle between norm of detector plane and global y&z surface (Rotate angle between axis and y&z surface)
	sinPhi=sin(0.);		cosPhi=cos(0.);			// the angle between the projection of the norm on the global y&z surface and global z axis (Rotate angle between axis1 (projection of axis on ) and x&z surface)
	sinBeta=sin(0.);	cosBeta=cos(0.);	*/	// the angle detector plane rotates around the norm (Rotate angle between axis1 (projection of axis on ) and x&z surface)
    sinAlpha=sin(aaa[5]);	cosAlpha=cos(aaa[5]);	// the angle between norm of detector plane and global y&z surface (Rotate angle between axis and y&z surface)
	sinPhi=sin(aaa[6]);		cosPhi=cos(aaa[6]);			// the angle between the projection of the norm on the global y&z surface and global z axis (Rotate angle between axis1 (projection of axis on ) and x&z surface)
	sinBeta=sin(aaa[4]);	cosBeta=cos(aaa[4]);		// the angle detector plane rotates around the norm (Rotate angle between axis1 (projection of axis on ) and x&z surface)
	
	/**  source geometry **/
    sdx=SDX; sdy=SDY;  sdz=SDZ; //source pixel dimension			// === CHANGE ===//
	//scx=3.2717; scy= -2.9958 ; scz=0.;	//centre of the source object		//**** CHANGE FOR RECON (2) ****// 	
	scx=aaa[10*nph+12+2]; scy= aaa[10*nph+13+2]; scz=aaa[10*nph+14+2];
	//scx=-1.9999005; scy=-4.8940496 ; scz=0.41262797;

	/** threshold for acceptable SRF elements **/
	thresh=1.e-12;
	/**** defining basic geometry (end) ** **/

	/**** calculating detector positions (start) ****/
	dn=0;
	if (1) {
		if (1) {
			if (nph==192) sprintf(filename1, "%s\\data\\det_pix_mask1.txt", COMMON_DIR); 
			else if (nph==7) sprintf(filename1, "%s\\data\\det_pix_mask2.txt", COMMON_DIR); 
			fp1=fopen(filename1, "r"); 
			while (fp1==NULL) {
				Sleep(1000);
				printf("."); 
				fp1=fopen(filename, "r"); 
			}
			i=0;
			while (!(feof(fp1))) {
				fscanf(fp1, "%f", &(det_mask[i]));
				i++;
			}
			if (!(i==NDET_DET1+1)) {printf("error in data reading! i=%d\n", i);getchar();}
			fclose(fp1);
		}
		else {
			dn=0;
			for(ndy=0; ndy<DN_DET1; ndy++) {
				for(ndx=0; ndx<DN_DET1; ndx++) {
					det_mask[dn]=1;	
					dn++;
				}
			}
		}
	}
	else {
		sprintf(filename, "%s\\data\\det_pix_mask", COMMON_DIR); fp1=fopen(filename, "wb"); 
		printf("DETMASKX=%d, DETMASKY=%d, DETMASKRADIUS=%d",DET_MASK_CENTERX, DET_MASK_CENTERY, DET_MASK_RADIUS);
		for(ndy=0; ndy<DN_DET1; ndy++) {
			for(ndx=0; ndx<DN_DET1; ndx++) {
				/** detector pixel det_mask. =1: valid pixel; =0: invalid **/
				// since the detection area is round
				//if (sqrt((ndx-130)*(ndx-130)+(ndy-130)*(ndy-130))<130) {	// for det1
				if (sqrt((ndx-DET_MASK_CENTERX)*(ndx-DET_MASK_CENTERX)+(ndy-DET_MASK_CENTERY)*(ndy-DET_MASK_CENTERY))<DET_MASK_RADIUS) {	// for det2
				//if (sqrt((ndx-256)*(ndx-256)+(ndy-256)*(ndy-256))<256) {	// for det1
					det_mask[dn]=1;	
				}
				dn++;
			}
		}
		fwrite(dp, sizeof(float), NDET_DET1*3, fp1); fclose(fp1);
	}
	

	if (IDETPOS) {
		sprintf(filename1, "%s\\data\\det_position.txt", COMMON_DIR);
		printf("using %s as the detector position\n", filename1);
		fp1=fopen(filename1, "r"); 
		while (fp1==NULL) {
			Sleep(1000);
			printf("."); 
			fp1=fopen(filename, "r"); 
		}
		i=0;
		while (!(feof(fp1))) {
			if (i/2*2==i){fscanf(fp1, "%f", &(dxx0[i/2]));}
			else {fscanf(fp1, "%f", &(dyy0[(i-1)/2]));}
			i++;
		}
		if (!(i==NDET_DET1*2+1)) {printf("error in data reading! i=%d\n", i);getchar();}
		fclose(fp1);
	}
	/*for (i=0; i<NDET_DET1; i++)
		printf("i=%d, x=%f, y=%f\n", i, dxx0[i],dyy0[i]);
	getchar();*/

	dn=0;dsubn=0;
	for (dpi=0;dpi<DM_DET;dpi++){
		for(ndy=0; ndy<DN_DET1; ndy++) {
			for(ndx=0; ndx<DN_DET1; ndx++) {	
				if(dn!=(dpi*NDET_DET1+DN_DET1*ndy+ndx)) {	printf("\n %ld\t%ld\n", dn, ndx*ndy+ndx); getchar();}

				if (IDETPOS) {
					dx0=dxx0[DN_DET1*ndy+ndx];
					dy0=dyy0[DN_DET1*ndy+ndx];
				}
				else {
					dx0=(ndx-(DN_DET1-1)*.5)*ddx;
					dy0=(ndy-(DN_DET1-1)*.5)*ddy;
				}
				if(nph==192){	dz0=ddz*(dpi-(DM_DET-1)*.5);	} // crystal from near to far from the source
				if(nph==7){	dz0=ddz*(-dpi+(DM_DET-1)*.5);	}

				dx=dcx+cosAlpha*cosBeta*dx0+(-sinBeta*cosPhi+sinPhi*sinAlpha*cosBeta)*dy0+(sinBeta*sinPhi+cosBeta*sinAlpha*cosPhi)*dz0;							//x position of the pixel
				dy=dcy+cosAlpha*sinBeta*dx0+(cosPhi*cosBeta+sinPhi*sinAlpha*sinBeta)*dy0+(-cosBeta*sinPhi+sinBeta*sinAlpha*cosPhi)*dz0;			//y position of the pixel
				dz=dcz-sinAlpha*dx0+sinPhi*cosAlpha*dy0+cosAlpha*cosPhi*dz0;			//z position of the pixel
				dp[dn*3]=dx; dp[dn*3+1]=dy; dp[dn*3+2]=dz;
				//printf("dn %d dp %f %f %f\n",dn,dp[dn*3],dp[dn*3+1],dp[dn*3+2]);//getchar();
				dn++;
			
				for (j=0;j<fn;j++){
					for (i=0;i<fn;i++){
						dx=dcx+cosAlpha*cosBeta*(dx0+(i-(fn-1)/2)*ddx/fn)+(-sinBeta*cosPhi+sinPhi*sinAlpha*cosBeta)*(dy0+(j-(fn-1)/2)*ddy/fn)+(sinBeta*sinPhi+cosBeta*sinAlpha*cosPhi)*dz0;							//x position of the pixel
						dy=dcy+cosAlpha*sinBeta*(dx0+(i-(fn-1)/2)*ddx/fn)+(cosPhi*cosBeta+sinPhi*sinAlpha*sinBeta)*(dy0+(j-(fn-1)/2)*ddy/fn)+(-cosBeta*sinPhi+sinBeta*sinAlpha*cosPhi)*dz0;			//y position of the pixel
						dz=dcz-sinAlpha*(dx0+(i-(fn-1)/2)*ddx/fn)+sinPhi*cosAlpha*(dy0+(j-(fn-1)/2)*ddy/fn)+cosAlpha*cosPhi*dz0;			//z position of the pixel
						dsubp[dsubn*3]=dx; dsubp[dsubn*3+1]=dy; dsubp[dsubn*3+2]=dz;
						dsubn++;
					}
				}
			}
		}
	}
	     
	/**** calculating detector positions (end) ****/

	/**** directional vectors for the X, Y, and Z axes of the detector (vector n of detector plane (0,0,1)@cosAlpha=sinPhi=0) (start) **/
	detX[0]=cosAlpha*cosBeta;
	detX[1]=cosAlpha*sinBeta;
	detX[2]=-sinAlpha;
	detY[0]=-sinBeta*cosPhi+sinPhi*sinAlpha*cosBeta;
	detY[1]=cosPhi*cosBeta+sinPhi*sinAlpha*sinBeta;
	detY[2]=sinPhi*cosAlpha;
	detZ[0]=sinBeta*sinPhi+cosBeta*sinAlpha*cosPhi;
	detZ[1]=sinBeta*sinAlpha*cosPhi-cosBeta*sinPhi;
	detZ[2]=cosAlpha*cosPhi;
	printf("detX %f %f %f \n",detX[0],detX[1],detX[2]);
	printf("detY %f %f %f \n",detY[0],detY[1],detY[2]);
	printf("detZ %f %f %f \n",detZ[0],detZ[1],detZ[2]);
	//getchar();
	/**** directional vectors for the X, Y, and Z axes of the detector (end) **/

	/**** pinhole geometry settings (start)  ****/
	ct=PH_CT,cht=PH_CHT;						//thichness									// === CHANGE ===//
	hz1=cz[4]+ct*0.5; r1=2.02;	//top circle radius							// === CHANGE ===//
	hz2=cz[4];        r2=PH_RADIUS;	//pinhole radius							// === CHANGE ===//
	hz3=cz[4]-ct*0.5; r3=2.02;	//bottom circle radius						// === CHANGE ===//
	hdist=6.;					//distance between pinholes. seven in total	// === CHANGE ===//
	mu=MU_APER;					//attenuation coef. of Tungsten, in (mm-1)	// === CHANGE ===//
	mu_pb=PB_APER;
	/**  pinhole basic settings (end)  **/	

	fp1=fopen("ph_loca_det1", "w");	

	minhsc=NSX*SDX*10;

	/*cx=0; cy=0; czz=0;
	for (k=0; k<nph; k++){
		cx=cx+aaa[11+k*6+1];
		cy=cy+aaa[11+k*6+2];
		czz=czz+aaa[11+k*6+3];
	}
	cx=cx/nph;cy=cy/nph;czz=czz/nph;
	printf("cx=%f cy=%f cz=%f\n", cx, cy, czz);

	nx_pb=0;ny_pb=0;nz_pb=0;
	printf("k=0 x=%f y=%f z=%f\n",aaa[12], aaa[13], aaa[14]);
	for (k=1; k<nph; k=k+2) {
		for (tt=0; tt<=1; tt++) {
			hx[k+tt]=aaa[11+(k+tt)*6+1];
			hy[k+tt]=aaa[11+(k+tt)*6+2];
			cz[k+tt]=aaa[11+(k+tt)*6+3];

			nxtmp[tt]=hx[k+tt]-cx;
			nytmp[tt]=hy[k+tt]-cy;
			nztmp[tt]=cz[k+tt]-czz;
			//printf("k=%d x=%f y=%f z=%f dist=%f\n", k+tt, hx[k+tt], hy[k+tt], cz[k+tt], sqrt((hx[k+tt]-cx)*(hx[k+tt]-cx)+(hy[k+tt]-cy)*(hy[k+tt]-cy)+(cz[k+tt]-czz)*(cz[k+tt]-czz)));
		}
		nxx=nytmp[0]*nztmp[1]-nytmp[1]*nztmp[0];
		nyy=nztmp[0]*nxtmp[1]-nztmp[1]*nxtmp[0];
		nzz=nxtmp[0]*nytmp[1]-nxtmp[1]*nytmp[0];

		if (k>=7 && k<=11) {
			nxx=-nxx;
			nyy=-nyy;
			nzz=-nzz;
		}
		tmp=sqrt(nxx*nxx+nyy*nyy+nzz*nzz);
		nx_pb=nx_pb+nxx/tmp;
		ny_pb=ny_pb+nyy/tmp;
		nz_pb=nz_pb+nzz/tmp;
		//printf("k=%d %f %f %f\n", k, nxx/sqrt(nxx*nxx+nyy*nyy+nzz*nzz),nyy/sqrt(nxx*nxx+nyy*nyy+nzz*nzz),nzz/sqrt(nxx*nxx+nyy*nyy+nzz*nzz));
	}//getchar();
	nx_pb/=(nph-1)/2;
	ny_pb/=(nph-1)/2;
	nz_pb/=(nph-1)/2;
	tmp=sqrt(nx_pb*nx_pb+ny_pb*ny_pb+nz_pb*nz_pb);
	nx_pb/=tmp;
	ny_pb/=tmp;
	nz_pb/=tmp;*/
	//printf("nx=%f ny=%f nz=%f %f\n", nx_pb, ny_pb, nz_pb, nx_pb*nx_pb+ny_pb*ny_pb+nz_pb*nz_pb);getchar();
		
	for (k=0; k<nph; k++){

		hx[k]=aaa[11+k*10+1];
		hy[k]=aaa[11+k*10+2];
		cz[k]=aaa[11+k*10+3];
	
		hsc[k]=sqrt(hx[k]*hx[k]+hy[k]*hy[k]+cz[k]*cz[k]);
		minhsc=min(minhsc, hsc[k]);

		//printf("k=%d dist=%f\n", k, sqrt((hx[k]-cx)*(hx[k]-cx)+(hy[k]-cy)*(hy[k]-cy)+(cz[k]-czz)*(cz[k]-czz)));
		//printf("%d %f %f %f %f\n",k,hx[k],hy[k],cz[k],hsc[k]);
		detAZx[k]=aaa[11+k*10+4];
		detAZy[k]=aaa[11+k*10+5];
		detAZz[k]=aaa[11+k*10+6]; 

		//printf("%f %f %f\n", detAZx[k], detAZy[k], detAZz[k]);

		/*nxx=hx[k]-cx;
		nyy=hy[k]-cy;
		nzz=cz[k]-czz;

		tmp=sqrt(nxx*nxx+nyy*nyy+nzz*nzz);
		nxx/=tmp;
		nyy/=tmp;
		nzz/=tmp;
		
		//printf("%f %f\n", nxx*nxx+nyy*nyy+nzz*nzz,nx_pb*nx_pb+ny_pb*ny_pb+nz_pb*nz_pb);getchar();

		detAXx[k]=nyy*nz_pb-nzz*ny_pb;
		detAXy[k]=nzz*nx_pb-nxx*nz_pb;
		detAXz[k]=nxx*ny_pb-nyy*nx_pb;
		tmp=sqrt(detAXx[k]*detAXx[k]+detAXy[k]*detAXy[k]+detAXz[k]*detAXz[k]);
		detAXx[k]/=tmp;
		detAXy[k]/=tmp;
		detAXz[k]/=tmp;

		detAYx[k]=detAZy[k]*detAXz[k]-detAZz[k]*detAXy[k];
		detAYy[k]=detAZz[k]*detAXx[k]-detAZx[k]*detAXz[k];
		detAYz[k]=detAZx[k]*detAXy[k]-detAZy[k]*detAXx[k];
		tmp=sqrt(detAYx[k]*detAYx[k]+detAYy[k]*detAYy[k]+detAYz[k]*detAYz[k]);
		detAYx[k]/=tmp;
		detAYy[k]/=tmp;
		detAYz[k]/=tmp;*/
		detAYx[k]=aaa[11+k*10+7];
		detAYy[k]=aaa[11+k*10+8];
		detAYz[k]=aaa[11+k*10+9];

		Beta[k]=aaa[11+k*10+10];

		printf("k=%d, cx=%f cy=%f cz=%f beta=%f\n", k,hx[k],hy[k],cz[k], Beta[k]);
		if(k==0) printf("detAZx=%f detAZy=%f detAZz=%f detAYx=%f detAYy=%f detAYz=%f\n",detAZx[k],detAZy[k],detAZz[k], detAYx[k],detAYy[k],detAYz[k]);

		//printf("%f %f\n", detAXx[k]*detAXx[k]+detAXy[k]*detAXy[k]+detAXz[k]*detAXz[k], detAYx[k]*detAYx[k]+detAYy[k]*detAYy[k]+detAYz[k]*detAYz[k]);getchar();
		/*detAYx[k]=0.;
		detAYy[k]=detAZz[k]/sqrt(detAZz[k]*detAZz[k]+detAZy[k]*detAZy[k]);
		detAYz[k]=-detAZy[k]/sqrt(detAZz[k]*detAZz[k]+detAZy[k]*detAZy[k]);*/
		//======= detAX = detAY * detAZ
		detAXx[k]=detAYy[k]*detAZz[k]-detAYz[k]*detAZy[k];
		detAXy[k]=detAYz[k]*detAZx[k]-detAYx[k]*detAZz[k];
		detAXz[k]=detAYx[k]*detAZy[k]-detAYy[k]*detAZx[k];

		tmp=detAZx[k]*detX[0]+detAZy[k]*detX[1]+detAZz[k]*detX[2];
		DDX_ph[k]=DDX_DET1*sqrt(1-tmp*tmp);
		tmp=detAZx[k]*detY[0]+detAZy[k]*detY[1]+detAZz[k]*detY[2];
		DDY_ph[k]=DDY_DET1*sqrt(1-tmp*tmp);
		//printf("%d %e %e\n", k, DDX_ph[k], DDY_ph[k]);
	}
	fclose(fp1);
	printf("minimum hsc = %f\n", minhsc);
	printf("<aperatt_ma_proj_sub_det>: Total no. of pinholes open: %d\n", nph);
	/**** pinhole geometry settings (end)  ****/
 
	/****  calculating map for penetration factor (start)  ****/
	printf("start to calculate map for penetration factor!\n");
	/** intermediate variables definition **/

	//tanAlpha=(r1-r2)/(ct/2);				// Tan value of half of acceptance angle
	AperAlpha=PI*APERALPHA/180.;					// half of acceptance angle
	tanAlpha=tan(PI*APERALPHA/180.);	// open angle of pinhole
	
	/*for (k=1; k<nph; k=k+2) {
		tmp=detAXx[k]*detAXx[k+1]+detAXy[k]*detAXy[k+1]+detAXz[k]*detAXz[k+1];
		printf("coscos=%f\n", tmp);
	}*/



	/*Beta[0]=0;
	if (nph==1) {
		for (m=1; m<7; m++) Beta[m]=10.9807*PI/180.;
		for (m=7; m<12; m++) Beta[m]=17.9079*PI/180.;
		for (m=13; m<nph; m++) Beta[m]=20.4623*PI/180.;
	}
	else if (nph==7) {
		for (m=1; m<7; m++) Beta[m]=17.4472*PI/180.;
	}*/

	for (m=0; m<nph; m++) {
		ll=LEAD_CHT/cos(Beta[m])/2.-((ct-cht)/2.*tanAlpha+r2)*tan(Beta[m]);
		tanBeta[m]=(ct-cht)/2*tanAlpha/(ll-cht/2);
	}

	//r4=r2-cht/2.*tan(AperAlpha);		//smaller circle radius	
	r5=r2+(ct-cht)/2*tanAlpha;
	printf("r5=%f\n",r5);
	/** calculation steps settings  **/
	rlimit=2*PH_RADIUS;
	dirt[0]=1.1*rlimit/N1;//(2.)/N1;				// rho:	distance between ph center and projection on the pinhole plane
	dirt[1]=1.1*rlimit/N2;//(2.)/N2;			// theta:	angle between norm of ph plane and the line from pinhole to source
	dirt[2]=PI*APERALPHA_LIMIT/180./N3;					// cos(beta-phi)	: cos value of difference between beta and phi. beta is the azimuthal angle of the projection of line from source to ph center on the ph plane. phi is the azimuthal angle of the projection of the ray on the ph plane. 
	dirt[3]=PI/N4;		
	
	//printf("%e\n %e\n", dirt[2], atan(tan(PI*45./180.))/N3);getchar();
	printf("\n M1 %d M2 %d r1 %e r2 %e r3 %e ",M1,M2,r1,r2,r3);
	printf("\n ct %e cht %e mu %e ",ct,cht,mu);
	printf("\n dirt[3] %e dirt[2] %e dirt[1] %e dirt[0] %e \n",dirt[3], dirt[2],dirt[1],dirt[0]);//getchar();
	/** calculation steps settings  **/

	/** read calculation_order file and decide running sequence (start) **/
	if (isequence) {
		if(1) {
			sprintf(filename2, "%s\\data\\srfs_order", COMMON_DIR);
			if((fp2=fopen(filename2, "rb"))!=NULL) {
				fread(&(iistart), sizeof(int), 1, fp2);
				fclose(fp2);

				fp2=fopen(filename2, "wb");
				iistart++; fwrite(&(iistart), sizeof(int), 1, fp2);	
				fclose(fp2);
				printf("start number is %d\n", iistart); // getchar();
			}
			else {
				fp2=fopen(filename2, "wb");
				iistart=1; fwrite(&(iistart), sizeof(int), 1, fp2);	
				fclose(fp2);
				printf("start number is %d\n", iistart); // getchar();
			}
		}

		for (i=1; i<=(NANG*NSUBDET); i++) { sequence[i]=0;}
		for (i=iistart; i<=(NANG*NSUBDET); i=i+NCOPY) { sequence[i]=1;printf("%d\n",i);}	// nf_current=1, sequence should be calculated, =0, should not be calculated.
		*iistart_out=iistart;

		if (iistart==1 || iistart==9 || iistart==17 || iistart==25) save_path();
	}
	/** read penetration map **/
	if (IMAP) {
		for (m=0; m<nph; m++) { 
			if (nph==7) sprintf(filename,"%s\\data\\map insert ph\\map_insert_ph%d", BUF_DIR1,m);
			else if (nph==192) sprintf(filename,"%s\\data\\map insert ph\\map_insert_ph%d", BUF_DIR1,m);
			//sprintf(filename,".\\data\\map_insert_ph1");
			printf("%s\n", filename);
			fp2=fopen(filename, "rb");
			for (i=0; i<N1*N2*N3*N4; i++) {
				fread(&(map[m][i]), sizeof(float), 1, fp2);
			}
			fread(&(dirt[0]), sizeof(float), 1, fp2);
			fread(&(dirt[1]), sizeof(float), 1, fp2);
			fread(&(dirt[2]), sizeof(float), 1, fp2);
			fread(&(dirt[3]), sizeof(float), 1, fp2);
			fclose(fp2);
		}
		
	}
	else {
		sprintf(filename2, "%s\\data\\map_order_start",BUF_DIR1);
		if((fp2=fopen(filename2, "rb"))!=NULL) {
			fread(&(iii), sizeof(int), 1, fp2);
			fclose(fp2);
			fp2=fopen(filename2, "wb");
			iii++; fwrite(&(iii), sizeof(int), 1, fp2);	
			fclose(fp2);
		}
		else {
			fp2=fopen(filename2, "wb");
			iii=0; fwrite(&(iii), sizeof(int), 1, fp2);	
			fclose(fp2);
		}

		for (m=iii; m<nph; m=m+NCPU) {
			m=0;
			printf("start to simulate ph# %d\n", m); // getchar();
			//penetration_map_layer_insert_ph(r1, r2, r3, ct, cht, LEAD_CHT, mu, mu_pb, map[m],dirt, tanAlpha, (Beta[m]));// creat multiple layer map
			//penetration_map_layer_symetric(r1, r2, r3, ct, cht, LEAD_CHT, mu, mu_pb, map[m],dirt, tanAlpha, (Beta[m]));// creat multiple layer map
			map[m]=vector(0, N1*N2*N3*N4-1);			// Map for penetration factor
			penetration_map_compound_eye(detAXx,detAXy, detAXz,detAYx,detAYy,detAYz,detAZx,detAZy, detAZz,dirt,m,nph,map[m],mu);
			ph_order[m]=1;
			sprintf(filename,"%s\\data\\map insert ph\\map_insert_ph%d", BUF_DIR1, m);
			printf("%s\n", filename);
			fp2=fopen(filename, "wb");
			for (i=0; i<N1*N2*N3*N4; i++) {
				fwrite(&(map[m][i]), sizeof(float), 1, fp2);
			}
			fwrite(&(dirt[0]), sizeof(float), 1, fp2);
			fwrite(&(dirt[1]), sizeof(float), 1, fp2);
			fwrite(&(dirt[2]), sizeof(float), 1, fp2);
			fwrite(&(dirt[3]), sizeof(float), 1, fp2);
			fclose(fp2);

			sprintf(filename2, "%s\\data\\map_order_finished_ph%d", BUF_DIR1, m);
			if((fp2=fopen(filename2, "wb"))!=NULL) {
				fclose(fp2);
			}
			printf("map for %d has been finished!\n", m);	
		}

		for (m=0; m<nph; m++) { 
			sprintf(filename1, "%s\\data\\map_order_finished_ph%d", BUF_DIR1, m);
			if (ph_order[m]==1) continue;
			
			while (1) {
				if((fp1=fopen(filename1, "rb"))!=NULL) {
					fclose(fp1);
					break;				
				}
				else {
					printf(".");
					Sleep(10000);
				}
			}

			sprintf(filename,"%s\\data\\map insert ph\\map_insert_ph%d", BUF_DIR1, m);
			printf("read map from %s\n", filename);
			if((fp2=fopen(filename, "rb"))!=NULL) {
				for (i=0; i<N1*N2*N3*N4; i++) {
					fread(&(map[m][i]), sizeof(float), 1, fp2);
				}
				fclose(fp2);
			}
			else printf("fail to open file %s\n", filename);
		}
	}
	printf("penetration factor map has been read!\n");



	/**  Rlimit definition **/
	for (m=0; m<nph; m++) {
		Rlimit[m]=0;
		for (i=0; i<N1; i++) {
			for (j=0; j<N2; j++) {
				k=0;l=0;
				if (exp(-map[m][i*M1+j*M2+k*N4+l])>0.1) {
					xx=(i-N1/2)*dirt[0];
					yy=(j-N2/2)*dirt[1];
					if (sqrt(xx*xx+yy*yy)>Rlimit[m]) {
						//printf("%e\n", map[i*M1+j*M2+k*N4+l]);
						Rlimit[m]=sqrt(xx*xx+yy*yy);
						//printf("%d %d %d %e %f %f %f\n", m,i,j,Rlimit[m], xx, yy, exp(-map[m][i*M1+j*M2+k*N4+l]));
					}
				}
			}
		}
		//Rlimit[m]=Rlimit[m]+0.4;
		Rlimit[m]=rlimit;
		printf("m=%d, Rlimit=%f\n", m, Rlimit[m]);//getchar();

	} 
	
	/****  calculating map for penetration factor (end)  ****/

	printf("<aperatt_ma_proj_sub_det>: start simulation ...\n");
restart:
	for (nos=1; nos<=NOS; nos++) {
		for(na=(nos-1); na<NANG; na=na+NOS) { //no. of subset
			//for(na=0; na<NANG; na++) {
			//na=14;
			printf("\n<aperatt_ma_proj_sub_det>: na: %d\n", na);

			if (0) {
				sprintf(filename1, "%s\\data\\pixel_mask", COMMON_DIR); 
				fp1=fopen(filename1, "rb"); 

				while (fp1==NULL) {
					printf("waiting for projection file %s ...\n", filename1);	
					Sleep(1000);
					printf("."); 
					fp1=fopen(filename1, "rb"); 
				}

				dn=0; ntmp=0; ntmp9=0;
				for(ndy=0; ndy<DN_DET1; ndy++) {
					for(ndx=0; ndx<DN_DET1; ndx++) {	
						fread(&(tmp), sizeof(float), 1, fp1);
						
						if (det_mask[dn]==1) {
							if (tmp==0.) {
								det_pix_mask[dn]=ntmp9%NSUBDET+101;
								ntmp9++;
							}
							else {
								//det_pix_mask[dn]=(int)((ntmp+0.5)/pixstep)+1;
								det_pix_mask[dn]=ntmp%NSUBDET+1;
								ntmp++;
							}
						}
						dn++;
					}
				}
				fclose(fp1);ntmp=0;
			}
			else {
				dn=0; ntmp=0; 
				for(ndy=0; ndy<DN_DET1; ndy++) {
					for(ndx=0; ndx<DN_DET1; ndx++) {	
						
						if (det_mask[dn]==1) {

								//det_pix_mask[dn]=(int)((ntmp+0.5)/pixstep)+1;
								det_pix_mask[dn]=ntmp%NSUBDET+1;
								ntmp++;

						}
						dn++;
					}
				}
				ntmp=0;
			}

			if (nph==192) sprintf(filename, "%s\\data\\det1_det_pix\\det_pix_mask%d", COMMON_DIR, na+1); 
			else if (nph==7) sprintf(filename, "%s\\data\\det2_det_pix\\det_pix_mask%d", COMMON_DIR, na+1); 
			printf("save det_pix_mask to file %s\n", filename);
			fp3=fopen(filename, "wb"); 
			while (fp3==NULL) {
				Sleep(1000);
				printf("."); 
				fp3=fopen(filename, "wb"); 
			}
			for (i=0; i<NDET_DET1; i++) fwrite(&(det_pix_mask[i]), sizeof(float), 1,	fp3);
			fclose(fp3);

			printf("done!\n");	
			printf("Number of pixels masked off by threshold: %ld\n", ntmp);

			/**** determine angle-dependent det_pix_mask (end) ****/

			/**** calculating source positions (start) ****/
			//now turning the source by a angle around X axis
			//Y direction: parallel to det; Z direction: papendicular to detector surface and pointing to the source.			sn=0; sen_max=0.;
			for(i=0; i<=NIMG-1; i++) { sou_mask[i]=0;}
			sn=0;
			for(nsz=0; nsz<NSZ; nsz++) {
				sz=scz+(nsz-(NSZ-1)*.5)*sdz; 
				for(nsy=0; nsy<NSY; nsy++) {
					sy=scy+(nsy-(NSY-1)*.5)*sdy; 
					for(nsx=0; nsx<NSX; nsx++) {
						sx=scx+(nsx-(NSX-1)*.5)*sdx; 
						d=sqrt((sy-scy)*(sy-scy)+(sz-scz)*(sz-scz));

						/** after rotation around the X axis **/
						if((sz-scz)>0.) alpha=acos((sy-scy)/d);  //angle between source point and y axis
						else alpha=2.*PI-acos((sy-scy)/d);
						alpha+=DANG*na; //angle after rotating
						//alpha-=DANG*na; //angle after rotating
						if(d>0){
							sp[sn*3]=sx; sp[sn*3+1]=scy+d*cos(alpha); sp[sn*3+2]=scz+d*sin(alpha);
						}
						else{ // d==0;
							sp[sn*3]=sx; sp[sn*3+1]=scy; sp[sn*3+2]=scz;
						}
						//sp[sn*3]=sx; sp[sn*3+1]=sy; sp[sn*3+2]=sz; //=============== without rotation for cubic geometry
						/** define source mask. =1: simulated with non-zero detection probability; =0: assume not detectable (do not simulate) **/
						if (d<minhsc*0.9) sou_mask[sn]=1.0;

						sn++;
					}
				}
			}	
			sprintf(filename, "%s\\data\\sou_mask%d", BUF_DIR1, na); 
			fp3=fopen(filename, "wb"); 
			while (fp3==NULL) {
				Sleep(1000);
				printf("."); 
				fp3=fopen(filename, "wb"); 
			}
			fwrite(sou_mask, sizeof(float), NIMG,	fp3);
			fclose(fp3);
			/**** calculating source positions (end) ****/

			for (nsubdet=1; nsubdet<=NSUBDET; nsubdet++) {
				if (isequence) {
					sprintf(filename2, "%s\\data\\file_read_order", BUF_DIR1);
					fp2=fopen(filename2, "rb");
					if(fp2!=NULL ) {
						while ((fread(&(nt1), sizeof(int), 1, fp2)==1) && (fread(&(nt2), sizeof(int), 1, fp2)==1) && (fread(&(nt3), sizeof(int), 1, fp2)==1) && (fread(&(nt4), sizeof(int), 1, fp2)==1) && (fread(&(nt5), sizeof(int), 1, fp2)==1)) {
							//printf("%d\t%d %d %d %d\n", nt1, nt2,nt3,nt4,nt5); 
							//printf("%d\t%d %d\n", nos, na,nsubdet); getchar();
							if(nt2==nos && nt3==na && nt4==nsubdet && nt5==1) {
								nf_current=nt1;	
								//printf("%d\n", nf_current);
								break;
							}
						}
						fclose(fp2);
					}

					if (sequence[nf_current]==0) {	
						continue;   
					}
					else			
					{
						printf("\n*** Step 1: start reading in srf for: %d\tnos: %d\tna: %d\tnsubdet: %d\tni: %d\n", nf_current, nt2, nt3, nt4, nt5); // getchar();
					}
				}

				for (i=0; i<NIMG; i++) sen[i]=0.;
				ntmp1=(nsubdet-1)*(DN_DET1/NSUBDET);
				ntmp2=(nsubdet)*(DN_DET1/NSUBDET);
				//printf("ntmp1, ntmp2: %d, %d\n", ntmp1, ntmp2);
			
				/* make sure the file is needed (start) */
				if(0) {
					// see if the file has been used in recon
					sprintf(filename2, "%s\\srfs\\file_recon_used", COMMON_DIR);
					//printf("check whether the file is already used and recorded in <%s>\n", filename2);
					fp2=fopen(filename2, "rb");
					//printf("ok1\n");
					if((fp2!=NULL) && (!feof(fp2))) {
						while ((fread(&(nt1), sizeof(int), 1, fp2)==1) && (fread(&(nt2), sizeof(int), 1, fp2)==1)) {
							//printf("%d\t%d\n", nt1, nt2); //getchar();
							if(nt1==na && nt2==nsubdet) {
								printf("match found 1 (na=%d, nsubdet=%d) in file_recon_used --> skip\n", na, nsubdet);
								fclose(fp2);
								goto next;
							}
						}
						fclose(fp2);
						printf("match not found (na=%d, nsubdet=%d) in file_recon_used\n", na, nsubdet);
					}
					else
						printf("error opening file %s\n", filename2);
				}

				if(1) {
					// see if the srf was started by other clones
					sprintf(filename2, "%s\\srfs\\file_MC_finished", COMMON_DIR);
					fp2=fopen(filename2, "rb");
					
					while (fp2==NULL) {
						Sleep(1000);
						printf("."); 
						fp2=fopen(filename2, "rb");
					}
					//printf("ok2\n");
					if((fp2!=NULL) && (!feof(fp2))) {
						while ((fread(&(nt1), sizeof(int), 1, fp2)==1) && (fread(&(nt2), sizeof(int), 1, fp2)==1)) {
							//printf("%d\t%d\n", nt1, nt2); //getchar();
							if(nt1==na && nt2==nsubdet) {
								printf("match found 2 (na=%d, nsubdet=%d) in file_MC_finished --> skip\n", na, nsubdet);
								fclose(fp2);
								goto next;
							}
						}
						fclose(fp2);
						printf("match not found (na=%d, nsubdet=%d) in file_MC_finished\n", na, nsubdet);
					}
					else {
						printf("error opening file %s ... restarting from beginning ...\n", filename2);
						close(fp2);
						goto restart;
					}
				}

				if(1) {
					// see if the srf was started by other clones
					sprintf(filename2, "%s\\srfs\\file_MC_started", COMMON_DIR);
					fp2=fopen(filename2, "rb");

					while (fp2==NULL) {
						Sleep(1000);
						printf("."); 
						fp2=fopen(filename2, "rb");
					}
					//printf("ok2\n");
					if((fp2!=NULL) && (!feof(fp2))) {
						while ((fread(&(nt1), sizeof(int), 1, fp2)==1) && (fread(&(nt2), sizeof(int), 1, fp2)==1)) {
							//printf("%d\t%d\n", nt1, nt2); //getchar();
							if(nt1==na && nt2==nsubdet) {
								printf("match found 2 (na=%d, nsubdet=%d) in file_MC_started --> skip\n", na, nsubdet);
								fclose(fp2);
								goto next;
							}
						}
						fclose(fp2);
						printf("match not found (na=%d, nsubdet=%d) in file_MC_started\n", na, nsubdet);
					}
					else
						printf("error opening file %s\n", filename2);
				}

				istart=0; ifinished=0;
				sprintf(filename2, "%s\\srfs\\file_MC_started", COMMON_DIR);
				fp2=fopen(filename2, "ab");  
				fwrite(&(na), sizeof(int), 1, fp2);	
				fwrite(&(nsubdet), sizeof(int), 1, fp2);
				fclose(fp2);
				istart=1;
				/* make sure the file is needed (end) */
				printf("na=%d Dang=%f alpha=%f\n", na, DANG, alpha);

				/** calculating projection parameters (start) **/
				for(i=1; i<=NMAX; i++) { sat[i]=0.; ijat[i]=0; }
				//ijat[1]=NIMG+2;
				//sprk=NIMG+1;
				ijat[1]=NSIZE+2;
				sprk=NSIZE+1;		
				for(dn=0; dn<NDET_DET1; dn++) proj[dn]=0.;
				for(sn=0; sn<NIMG; sn++) {
					//sn=10;
					if((sn/10000*10000)==sn) { printf("sn: %d %ld estimated NMAX=%ld\t", sn,sprk, (long int)(sprk/(sn+1)*NIMG)); time(&time1); printf("%s", asctime(localtime(&time1)));}
					npix=0;
					
					//if (sn!=764105) continue;

					if (sou_mask[sn]>0.5) {
						sx=sp[sn*3]; 
						sy=sp[sn*3+1]; 
						sz=sp[sn*3+2]; 
						//sx=0;sy=0;sz=0;
						//sx=-3.81;sy=0.32;sz=-0.19;
						//sx=-5.83;sy=-0.50;sz=0.26;             //in
						//sx=-5.41;sy=-0.50;sz=0.26;          //out
						//printf("\n sx=%f,sy=%f,sz=%f",sx,sy,sz); getchar();
						if (sqrt((sz-scz)*(sz-scz)+(sy-scy)*(sy-scy))>minhsc) {printf("source pixel is out of ph position\n"); getchar();}

						/* projection parameters (start) */			
						//f=(dz-sz)/(cz-sz);
						for(k=0; k<nph; k++) {
						//for(k=173;k<174;k++){
							spos[0]=(sx-hx[k])*detAXx[k]+(sy-hy[k])*detAXy[k]+(sz-cz[k])*detAXz[k]; //x position of source pixel in pinhole geometry 
							spos[1]=(sx-hx[k])*detAYx[k]+(sy-hy[k])*detAYy[k]+(sz-cz[k])*detAYz[k];
							spos[2]=(sx-hx[k])*detAZx[k]+(sy-hy[k])*detAZy[k]+(sz-cz[k])*detAZz[k];// === CHANGE ===//
							//printf("sx=%f, sy=%f, sz=%f\n", spos[0], spos[1], spos[2]);
							/*if (fabs(spos[0])<0.2 && fabs(spos[1])<0.2) {
								printf("%d %f %f\n", sn, spos[0], spos[1]);getchar();
							}
							else continue;*/
     						dsq=(sx-hx[k])*(sx-hx[k])+(sy-hy[k])*(sy-hy[k])+(sz-cz[k])*(sz-cz[k]); //square of source-pinhole distance
							//printf("dist=%f\n", dsq);getchar();
							/*sprintf(filename, "%s\\srfs\\dps_na%d_sn%d", BUF_DIR1,na,sn); fp3=fopen(filename, "wb"); 
							fwrite(&(dsq), sizeof(float),1 ,fp3);
							fclose(fp3);continue;*/

							hs=fabs(detAZx[k]*(sx-hx[k])+detAZy[k]*(sy-hy[k])+detAZz[k]*(sz-cz[k]));   // distance source to pinhole plane
							//================================ n vector  of source&pinhole  
							nsrx=(sx-hx[k])/sqrt(dsq);
							nsry=(sy-hy[k])/sqrt(dsq);
							nsrz=(sz-cz[k])/sqrt(dsq);
							cosTheta=fabs(hs/sqrt(dsq));
							cosThetaD=fabs(nsrx*detZ[0]+nsry*detZ[1]+nsrz*detZ[2]); //incident angle to detector surface
							sinThetaD=sqrt(1-cosThetaD*cosThetaD);

							buf[5]=r2;					// radius of ph
							//buf[15]=r4;
							
							if(cosTheta<=0.9999) buf[6]=acos(cosTheta);		// angle between norm of ph plane and the line from source to pinhole center
							if(buf[6]>PI/2.) buf[6]=PI-buf[6];  //angle between source point and y axis
							if(cosTheta>0.9999) buf[6]=0.;
							buf[7]=ddx;				// projection of x dimension of detector pixel on the pinhole plane
							buf[8]=ddy;				// projection of y dimension of detector pixel on the pinhole plane
							buf[9]=Rlimit[k];				// calculation limit in rho
							
							//buf[10]=atan(tanBeta[k]);			// half of acceptance angle of pinhole	
							//buf[15]=r2-cht/2*tanBeta[k];
							buf[10]=atan(tan(APERALPHA*PI/180));			// half of acceptance angle of pinhole	
							buf[15]=r2-cht/2*tan(APERALPHA*PI/180);
							buf[11]=hs;					// distance between source and ph plane
							buf[13]=ddz;				// detector layer thickness
							hsa=fabs((sx-hx[k])*detZ[0]+(sy-hy[k])*detZ[1]+(sz-cz[k])*detZ[2]);// distance bwtween source and ph on detector plane
							ntmp3=DN_DET1-1;ntmp4=0;ntmp5=DN_DET1-1;ntmp6=0;
							buf[3]=(r2*r2)/(4*dsq)*cosTheta; // prob of falling through a pinhole
							dp2=1000;dpx2=1000;dpy2=1000;
							for(dpi=0; dpi<DM_DET; dpi++) {
							//for(dpi=0; dpi<1; dpi++) {
							/**  Defination of buf and index for subroutine of penetration factor (start)  **/
								hsp=fabs(detZ[0]*(sx-dcx)+detZ[1]*(sy-dcy)+detZ[2]*(sz-dcz))+ddz*(dpi-(DM_DET-1)*.5); // distance between source and detector plate
								f=hsp/hsa;	
								buf[0]=sx+(hx[k]-sx)*f;		// x of the center of the projection
								buf[1]=sy+(hy[k]-sy)*f;		// y of the center of the projection
								buf[4]=sz+(cz[k]-sz)*f;		// z value of detector plane
								buf[2]=f*(Rlimit[k]);			// radius of the projection using Rlimit for penetration factor

								if (IDETPOS){
									proj_cx=(buf[0]-dcx)*detX[0]+(buf[1]-dcy)*detX[1]+(buf[4]-dcz)*detX[2];
									proj_cy=(buf[0]-dcx)*detY[0]+(buf[1]-dcy)*detY[1]+(buf[4]-dcz)*detY[2];
									//printf("proj_cx %f, y%f\n",proj_cx,proj_cy);
									dtmpx=(int)(proj_cx/ddx+DN_DET1*0.5);
									dtmpy=(int)(proj_cy/ddy+DN_DET1*0.5);
									if((proj_cy>=dyy0[0]&&proj_cy<=dyy0[NDET_DET1-1])&&(proj_cx>=dxx0[0]&&proj_cx<=dxx0[DN_DET1-1])){ // projection center are in the detector space;
										for(ndy=max(0,dtmpy-2); ndy<=min(dtmpy+2,DN_DET1); ndy++) {
											for(ndx=max(0,dtmpx-2); ndx<=min(dtmpx+2,DN_DET1); ndx++) {
												dn=dpi*NDET_DET1+DN_DET1*ndy+ndx;
												dpx=(dp[dn*3]-buf[0])*detX[0]+(dp[dn*3+1]-buf[1])*detX[1]+(dp[dn*3+2]-buf[4])*detX[2];
												dpy=(dp[dn*3]-buf[0])*detY[0]+(dp[dn*3+1]-buf[1])*detY[1]+(dp[dn*3+2]-buf[4])*detY[2];
												dp1=dpx*dpx+dpy*dpy;
												if (dp1<dp2) {dncx=ndx;dncy=ndy;dp2=dp1;}
											}
										}
										//printf("dncx %d dncy %d \n",dncx,dncy);
									}
									else if(proj_cx<dxx0[0]&&proj_cy>=dyy0[1]&&proj_cy<=dyy0[NDET_DET1-1]){
										dncx=(int)((proj_cx-dxx0[0])/ddx); // x direction projection center is outside the detector;
										dncy=dtmpy;
                                        for(ndy=max(0,dtmpy-2); ndy<=min(dtmpy+2,DN_DET1); ndy++){
											dn=dpi*NDET_DET1+DN_DET1*ndy;
											dpy=(dp[dn*3]-buf[0])*detY[0]+(dp[dn*3+1]-buf[1])*detY[1]+(dp[dn*3+2]-buf[4])*detY[2];
											dpy1=dpy*dpy;
											if(dpy1<dpy2){
											dpy2=dpy1;
											dncy=ndy;											}
										}
									//printf("1 dncx %d dncy %d \n",dncx,dncy);	
									}
									else if(proj_cx>dxx0[DN_DET1-1]&&proj_cy>=dyy0[1]&&proj_cy<=dyy0[NDET_DET1-1]){ // x direction projection center is outside the detector;
										dncx=(int)((proj_cx-dxx0[DN_DET1-1])/ddx+DN_DET1);
										dncy=dtmpy;
										for(ndy=max(0,dtmpy-2); ndy<=min(dtmpy+2,DN_DET1); ndy++){
											dn=dpi*NDET_DET1+DN_DET1*ndy+DN_DET1-1;
											dpy=(dp[dn*3]-buf[0])*detY[0]+(dp[dn*3+1]-buf[1])*detY[1]+(dp[dn*3+2]-buf[4])*detY[2];
											dpy1=dpy*dpy;
											if(dpy1<dpy2){
											dpy2=dpy1;
											dncy=ndy;											
											}
										}
										//printf("2 dncx %d dncy %d \n",dncx,dncy);
									}
									else if(proj_cy<dyy0[0]&&proj_cx>=dxx0[0]&&proj_cx<=dxx0[DN_DET1-1]){ // y direction out of detector;
										dncy=(int)((proj_cy-dyy0[0])/ddy);
										dncx=dtmpx;
										for(ndx=max(0,dtmpx-2); ndx<=min(dtmpx+2,DN_DET1); ndx++) {
											dn=dpi*NDET_DET1+ndx;
											dpx=(dp[dn*3]-buf[0])*detX[0]+(dp[dn*3+1]-buf[1])*detX[1]+(dp[dn*3+2]-buf[4])*detX[2];
											dpx1=dpx*dpx;
										    if (dpx1<dpx2){
												dncx=ndx;
										        dpx2=dpx1;
											}
										}
										//printf("3 dncx %d dncy %d \n",dncx,dncy);
									}
									else if(proj_cy>dyy0[NDET_DET1-1]&&proj_cx>=dxx0[0]&&proj_cx<=dxx0[DN_DET1-1]){
										dncy=(int)((proj_cy-dyy0[NDET_DET1-1])/ddy+DN_DET1);
										dncx=dtmpx;
										for(ndx=max(0,dtmpx-2); ndx<=min(dtmpx+2,DN_DET1); ndx++) {
											dn=dpi*NDET_DET1+NDET_DET1-1+ndx;
											dpx=(dp[dn*3]-buf[0])*detX[0]+(dp[dn*3+1]-buf[1])*detX[1]+(dp[dn*3+2]-buf[4])*detX[2];
											dpx1=dpx*dpx;
										    if (dpx1<dpx2){
												dncx=ndx;
										        dpx2=dpx1;
											}
										}
										//printf(" 4 dncx %d dncy %d \n",dncx,dncy);
									}
									else if((proj_cy<dyy0[0]||proj_cy>dyy0[NDET_DET1-1])&&(proj_cx>dxx0[0]||proj_cx<dxx0[DN_DET1-1])){// x y direction are out of detector;
										    if(proj_cx<dxx0[0]) dncx=(int)((proj_cx-dxx0[0])/ddx);
											if(proj_cx>dxx0[DN_DET1-1]) dncx=(int)((proj_cx-dxx0[DN_DET1-1])/ddx+DN_DET1);
											if(proj_cy<dyy0[0]) dncy=(int)((proj_cy-dyy0[0])/ddy);
											if(proj_cy>dyy0[NDET_DET1-1]) dncy=(int)((proj_cy-dyy0[NDET_DET1-1])/ddy+DN_DET1);
											//printf(" 5 dncx %d dncy %d \n",dncx,dncy);
											
									}
															
         //                           dtmpx=(int)(proj_cx/ddx+DN_DET1*0.5);
									//dtmpy=(int)(proj_cy/ddy+DN_DET1*0.5);
									//dp2=1000;

									//if (dtmpx>12 && dtmpy>65 && dtmpx<17 && dtmpy<70){
								   /* printf("k %d dtmpy=%d&& dtmpx=%d ok\n",k,dtmpy,dtmpx); 
									printf("buf[0] %f,buf[1] %f, buf[4] %f\n",buf[0],buf[1],buf[4]);
									printf("dcx %f dcy %f dcz %f\n",dcx,dcy,dcz);
									printf("detX[0] %f %f %f\n",detX[0],detX[1],detX[2]);
									printf("detY[0] %f %f %f\n",detY[0],detY[1],detY[2]);*/
									//printf("%d %d %d %d",max(0,dtmpy-2),min(dtmpy+2,DN_DET1),max(0,dtmpx-2),min(dtmpx+2,DN_DET1));getchar();
									/*for(ndy=max(0,dtmpy-2); ndy<=min(dtmpy+2,DN_DET1); ndy++) {
										for(ndx=max(0,dtmpx-2); ndx<=min(dtmpx+2,DN_DET1); ndx++) {
											dn=dpi*NDET_DET1+DN_DET1*ndy+ndx;
											dpx=(dp[dn*3]-buf[0])*detX[0]+(dp[dn*3+1]-buf[1])*detX[1]+(dp[dn*3+2]-buf[4])*detX[2];
											dpy=(dp[dn*3]-buf[0])*detY[0]+(dp[dn*3+1]-buf[1])*detY[1]+(dp[dn*3+2]-buf[4])*detY[2];
											dp1=dpx*dpx+dpy*dpy;
										if (dp1<dp2) {dncx=ndx;dncy=ndy;dp2=dp1;printf("dp1 %f dp2 %f dncx %d dncy %d",dp1, dp2, dncx,dncy);getchar();}
										}
									}*/

									/*if (!((dncx==dtmpx)&&(dncy==dtmpy))) {
										printf("sn=%d k=%d\n", sn, k);
									}*/
								}
								else {	//??? dcz layer
									dncx=(int)(((buf[0]-dcx)*detX[0]+(buf[1]-dcy)*detX[1]+(buf[4]-dcz)*detX[2])/ddx+DN_DET1*0.5);
									dncy=(int)(((buf[0]-dcx)*detY[0]+(buf[1]-dcy)*detY[1]+(buf[4]-dcz)*detY[2])/ddy+DN_DET1*0.5);
								}
								//printf("%d %d\n", dncx, dncy);
								deltan=fabs(buf[2]/cosTheta/cosThetaD/DDX_DET1)+2;
								ntmp3=min(dncy-deltan, ntmp3);
								ntmp4=max(dncy+deltan, ntmp4);
								ntmp5=min(dncx-deltan,ntmp5);
								ntmp6=max(dncx+deltan, ntmp6);
							//printf("\n deltan %f dncx %d y %d ",deltan, dncx, dncy);getchar();
							}	

							ntmp3=max(ntmp3, 0);
							ntmp4=min(ntmp4, DN_DET1-1);
							ntmp5=max(ntmp5,0);
							ntmp6=min(ntmp6, DN_DET1-1);

							for(ndy=ntmp3; ndy<=ntmp4; ndy++) {
								for(ndx=ntmp5; ndx<=ntmp6; ndx++) {
						    //for(ndy=(dncy-2); ndy<=(dncy+2); ndy++) {
								//for(ndx=(dncx-2); ndx<=(dncx+2); ndx++) {
									//if(((ndx>=0)&&(ndx<DN_DET1))&&((ndy>=ntmp1)&&(ndy<ntmp2))) {
									//if (ndx>12 && ndy>65 && ndx<17 && ndy<70){
										//printf("dncy=%d&& dncx=%d ok\n",dncy,dncx); getchar();}
									dn=DN_DET1*ndy+ndx;
									//if (dn!=13624) continue;

									//dn=DN_DET1*206+124;	
									//if ((det_pix_mask[dn]==nsubdet)||((nsubdet==1)&&(det_pix_mask[dn]==101))) {
									
									if (((det_pix_mask[dn]==nsubdet)||(det_pix_mask[dn]==100+nsubdet))&&fabs(dxx0[dn]-compx[k])<=radius_cp[k]&&fabs(dyy0[dn]-compy[k])<=radius_cp[k]){
									//if (((det_pix_mask[dn]==nsubdet)||(det_pix_mask[dn]==100+nsubdet))){	//// projection on ellipse short and long - axis
										////elap= ela[0]*(dpos[0]-buf[0]) +ela[1]*(dpos[1]-buf[1]) + ela[2]*(dpos[2]-buf[4]);
										////esap= esa[0]*(dpos[0]-buf[0]) +esa[1]*(dpos[1]-buf[1]) + esa[2]*(dpos[2]-buf[4]);
										//	if(abs(dxx0[dn]-compx[k])<=radius_cp[k]) printf("it is here,fucking error xx %f \n",(dxx0[dn]-compx[k]));
										//	if(abs(dyy0[dn]-compy[k])<=radius_cp[k]) printf("it is here,fucking error yy %f \n",abs(dyy0[dn]-compy[k]));
										//printf("ddx0 %f ddy0 %f  dn %d ndx %d ndy %d",dxx0[dn],dyy0[dn],dn,ndx,ndy);
										//printf("\ndp x %f y %f z % \n",dp[dn*3],dp[dn*3+1],dp[3*dn+2]);
										//printf("k=%d\n",k);
										//printf("compx %f compx %f  %f radius %f \n",compx[k],compy[k],dcz,radius_cp[k]);
										//printf("pinholex %f y %f  z%f\n",hx[k],hy[k],cz[k]);
										//printf("detAZ %f %f %f\n",detAZx[k],detAZy[k],detAZz[k]);
										//printf("detAZ %f %f %f\n",detAYx[k],detAYy[k],detAYz[k]);
										//printf("detAZ %f %f %f\n",detAXx[k],detAXy[k],detAXz[k]);
										//printf("line from pinhol to detector %f %f %f\n",compx[k]-hx[k],compy[k]-hy[k],dcz-cz[k]);
										//getchar();


										// vector detector pixel	(dpos[0]-buf[0])  (dpos[1]-buf[1])  (dpos[2]-buf[4])  
										// to record subpixel position in pinhole coordinate
										nt=0;
										for(dpi=0; dpi<DM_DET; dpi++) {
											ntmp7=(dn+dpi*NDET_DET1)*3;ntmp8=3*dpi;
											//printf("\n dn %d dp %e %e %e",dn,dp[ntmp7],dp[ntmp7+1],dp[ntmp7+2]);
											dpos[0+ntmp8]=(dp[ntmp7]-hx[k])*detAXx[k]+(dp[ntmp7+1]-hy[k])*detAXy[k]+(dp[ntmp7+2]-cz[k])*detAXz[k]; //x position of source pixel in pinhole geometry 
											dpos[1+ntmp8]=(dp[ntmp7]-hx[k])*detAYx[k]+(dp[ntmp7+1]-hy[k])*detAYy[k]+(dp[ntmp7+2]-cz[k])*detAYz[k];
											dpos[2+ntmp8]=(dp[ntmp7]-hx[k])*detAZx[k]+(dp[ntmp7+1]-hy[k])*detAZy[k]+(dp[ntmp7+2]-cz[k])*detAZz[k];
											nsrx=(sx-dp[ntmp7]);
											nsry=(sy-dp[ntmp7+1]);
											nsrz=(sz-dp[ntmp7+2]);
											cos_D[dpi]=fabs(nsrx*detZ[0]+nsry*detZ[1]+nsrz*detZ[2])/sqrt(nsrx*nsrx+nsry*nsry+nsrz*nsrz); // for each pixel
											nm=0;
											for (j=0;j<fn;j++){
												for (i=0;i<fn;i++){
													tmpdx=dsubp[((dn+dpi*NDET_DET1)*fn*fn+nm)*3];
													tmpdy=dsubp[((dn+dpi*NDET_DET1)*fn*fn+nm)*3+1];
													tmpdz=dsubp[((dn+dpi*NDET_DET1)*fn*fn+nm)*3+2];
													subdpos[nt*3]=(tmpdx-hx[k])*detAXx[k]+(tmpdy-hy[k])*detAXy[k]+(tmpdz-cz[k])*detAXz[k]; //x position of source pixel in pinhole geometry 
													subdpos[nt*3+1]=(tmpdx-hx[k])*detAYx[k]+(tmpdy-hy[k])*detAYy[k]+(tmpdz-cz[k])*detAYz[k];
													subdpos[nt*3+2]=(tmpdx-hx[k])*detAZx[k]+(tmpdy-hy[k])*detAZy[k]+(tmpdz-cz[k])*detAZz[k];
													nt++;nm++;
												}
											}
										}
											
										//pr=pixelprob_det1(buf, dpos); //calculate the prob of incident on give pixel																			
										//pr=pixelprob_det1NL(buf, dpos, spos, map, dirt, index);	
										//pr=pixelprob_det2ML(buf, dpos, spos, map, dirt, index);
										pr=pixelprob_det2ML_sum_ph_insert3(buf, dpos, subdpos, spos, map[k], dirt, cos_D, DM_DET, MU_DET, fn);	
										//printf("dn=%d, prob=%12.6e\n", dn, pr);

										/*if (sn==3 && k==18 && dn==75627) {
											printf("%d %d %d %e\n", sn, k, dn, pr);
											printf("%d %e\n", DM_DET, MU_DET);
											printf("%e %e %e\n", spos[0], spos[1], spos[2]);
											for (i=1;i<=13; i++) printf("%d %e\n", i, buf[i]);
											getchar();
										}*/
										//if((pr>1)) {printf("sn=%d k=%d dn=%d\n", sn, k, dn);getchar();getchar();}
										if(pr>thresh) {											
											sen[sn]+=pr;
											//proj_buf[npix]=pr;
											//printf("sn %d dn %d pr %e\n", sn, dn, pr);//getchar();
											if(det_pix_mask[dn]==nsubdet) {
												proj[dn]=pr;
												//det_prob[dn]=pr;
												//printf("\n dn=%d, prob=%e", dn, pr);getchar();
												dn_buf[npix]=dn; // record number of pixels with prob >0
												prt+=pr;
												npix++;
											}											
										}
										else {
											proj[dn]=0.;
										}
									}// det_mask
								}//ndy
							}//ndx
						}// k nph
					}// sn source
					//printf("xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx\n");
					//printf("sn=%d, sen=%e\n", sn, sen[sn]);getchar();
					//=============== fg debug probability file
					if(0){
					sprintf(filename, "%s\\srfs\\det_sn%d_ph%d_insert_ph", BUF_DIR1,sn, k); fp3=fopen(filename, "wb"); 
					for (i=0; i<NDET_DET1; i++) fwrite(&(proj[i]), sizeof(float),1 ,	fp3);
					fclose(fp3);
					printf("\n <aperatt_ma_proj_sub_det_cube_pixmask.c> det_prob %s", filename);
					getchar();getchar();getchar();getchar();getchar();getchar();getchar();
					}
					//=================================================

					//if(sen[sn]>sen_max) sen_max=sen[sn];
					/* projection parameters (end) */	

					if(1) {
						n=0;
						spri=sn+1;
						if(spri<=NDET_DET1) { 
							sat[spri]=proj[spri-1]; 
							proj[spri-1]=0.; 
						} 
						else sat[spri]=0.;

						//off-diagnal elements
						for(l=0; l<npix; l++) {
							sprj=dn_buf[l]+1;
							//for(sprj=1; sprj<=NDET_DET1; sprj++) {
							if((fabs(proj[sprj-1])>=thresh) && (sprj!=spri)) {
								if(++sprk>NMAX) { printf("error\n"); getchar(); }
								sat[sprk]=proj[sprj-1]; 
								ijat[sprk]=sprj;
								proj[sprj-1]=0.;
							}
						}
						ijat[spri+1]=sprk+1; 
					}
					/* optimised version (end) */

					/* original version (start) */
					//diagnal elements
					if(0) {
						spri=sn+1;
						if(spri<=NDET_DET1) { sat[spri]=proj[spri-1]; proj[spri-1]=0.;  } 

						//off-diagnal elements
						for(sprj=1; sprj<=NDET_DET1; sprj++) {
							if((fabs(proj[sprj-1])>=thresh) && (sprj!=spri)) {
								if(++sprk>NMAX) { printf("error\n"); getchar(); }
								sat[sprk]=proj[sprj-1]; proj[sprj-1]=0.;
								ijat[sprk]=sprj;
							}
						}
						ijat[spri+1]=sprk+1; 
					}
				}		

				for (spri=(NIMG+1); spri<=NSIZE; spri++) {
					ijat[spri+1]=ijat[NIMG+1];
				}
				if(1) {
					// see if the srf is in the started_not_finished list
					sprintf(filename2, "%s\\srfs\\file_MC_started_not_finished", COMMON_DIR);
					fp2=fopen(filename2, "rb");

					while (fp2==NULL) {
						Sleep(1000);
						printf("."); 
						fp2=fopen(filename2, "rb");
					}
					//printf("ok2\n");
					if((fp2!=NULL) && (!feof(fp2))) {
						while ((fread(&(nt1), sizeof(int), 1, fp2)==1) && (fread(&(nt2), sizeof(int), 1, fp2)==1)) {
							//printf("%d\t%d\n", nt1, nt2); //getchar();
							if(nt1==na && nt2==nsubdet) {
								printf("match found 3 (na=%d, nsubdet=%d) in file_MC_started_not_finished\n", na, nsubdet);
							}
						}
						fclose(fp2);
						printf("match not found (na=%d, nsubdet=%d) in file_MC_started\n", na, nsubdet);
					}
					else
						printf("error opening file %s\n", filename2);
				}

				if (isave==1) { // =1: save srfs at random locations evenly distributed over all HDs
					sprintf(filename, "%s\\srfs\\sen_na%d_det1_nsubdet%d", BUF_DIR1, na, nsubdet); fp3=fopen(filename, "wb"); 
					//sprintf(filename, "%s\\sen_na%d_det1_nsubdet%d", SRFS_DIR, na, nsubdet); fp3=fopen(filename, "wb"); 
					printf("<aperatt_ma_proj_sub_det>: writing sensitivity to file %s ... ", filename);
					fwrite(sen, sizeof(float), NIMG,	fp3);
					fclose(fp3);
					printf("finished!\n");
					/** calculating projection parameters (end) **/

					msize=ijat[ijat[1]-1]-1;
					printf("<aperatt_ma_proj_sub_det>: msize: %ld, %ld\n", msize, ijat[1]-2);
					
					sprintf(filename, "%s\\srfs\\ijat_na%d_det1_nsubdet%d", BUF_DIR1, na, nsubdet); fp1=fopen(filename, "wb"); 
					//sprintf(filename, "%s\\ijat_na%d_det1_nsubdet%d", SRFS_DIR, na, nsubdet); fp1=fopen(filename, "wb"); 
					printf("<aperatt_ma_proj_sub_det>: writing srf_%d to %s\n", na, filename);
					ntmp2=1;
					for(i=1; i<=msize; i++) {
						ntmp2=fwrite(&(ijat[i]), sizeof(unsigned long), 1, fp1);
						if(ntmp2!=1) {
							printf("<aperatt_ma_proj_sub_det>: error in writing ijat, na=%d and nsubdet=%d\n", na, nsubdet); 						
							fclose(fp1);

							ntmp2=0;
							printf("<aperatt_ma_proj_sub_det.c>: since file %s is not properly finished, it is deleted\n", filename);
							if( remove(filename) == -1 )
								printf("<aperatt_ma_proj_sub_det.c>: error 11: delete error: %s\n", filename);
							//getchar();
							goto restart;
						}
					}
					fclose(fp1);

					sprintf(filename, "%s\\srfs\\sat_na%d_det1_nsubdet%d", BUF_DIR1, na, nsubdet); fp1=fopen(filename, "wb"); 
					//sprintf(filename, "%s\\sat_na%d_det1_nsubdet%d", SRFS_DIR, na, nsubdet); fp1=fopen(filename, "wb"); 
					printf("<aperatt_ma_proj_sub_det>: writing srf_%d to %s\n", na, filename);
					ntmp1=1;
					for(i=1; i<=msize; i++) {
						ntmp1=fwrite(&(sat[i]), sizeof(float), 1, fp1);
						if(ntmp1!=1) {
							printf("<aperatt_ma_proj_sub_det>: error in writing sat, na=%d and nsubdet=%d\n", na, nsubdet); //getchar();
							ntmp1=0;
							printf("<aperatt_ma_proj_sub_det.c>: since file %s is not properly finished, it is deleted\n", filename);
							if( remove(filename) == -1 )
								printf("<aperatt_ma_proj_sub_det.c>: error 11: delete error: %s\n", filename);
							//getchar();
							goto restart;
						}
					}
					fclose(fp1);
				}
				else if (isave==2){ //isave=2: a complete copy on each HD -- save time in reconstruction
					for (m=1; m<=4; m++) {
						if (m==1) sprintf(filename, "%s\\srfs\\sen_na%d_det1_nsubdet%d", ROOT_DIR1, na, nsubdet); 
						if (m==2) sprintf(filename, "%s\\srfs\\sen_na%d_det1_nsubdet%d", ROOT_DIR2, na, nsubdet); 
						if (m==3) sprintf(filename, "%s\\srfs\\sen_na%d_det1_nsubdet%d", ROOT_DIR3, na, nsubdet); 
						if (m==4) sprintf(filename, "%s\\srfs\\sen_na%d_det1_nsubdet%d", ROOT_DIR4, na, nsubdet); 
						
						fp3=fopen(filename, "wb"); 
						//sprintf(filename, "%s\\sen_na%d_det1_nsubdet%d", SRFS_DIR, na, nsubdet); fp3=fopen(filename, "wb"); 
						printf("<aperatt_ma_proj_sub_det>: writing sensitivity to file %s ... ", filename);
						fwrite(sen, sizeof(float), NIMG,	fp3);
						fclose(fp3);
					}
					printf("finished saving sensitivity file ...\n");

					msize=ijat[ijat[1]-1]-1;
					printf("<aperatt_ma_proj_sub_det>: msize: %ld, %ld\n", msize, ijat[1]-2);
					
					for (m=1; m<=4; m++) {
						if (m==1) sprintf(filename, "%s\\srfs\\ijat_na%d_det1_nsubdet%d", ROOT_DIR1, na, nsubdet); 
						if (m==2) sprintf(filename, "%s\\srfs\\ijat_na%d_det1_nsubdet%d", ROOT_DIR2, na, nsubdet); 
						if (m==3) sprintf(filename, "%s\\srfs\\ijat_na%d_det1_nsubdet%d", ROOT_DIR3, na, nsubdet); 
						if (m==4) sprintf(filename, "%s\\srfs\\ijat_na%d_det1_nsubdet%d", ROOT_DIR4, na, nsubdet); 
						fp1=fopen(filename, "wb"); 
						//sprintf(filename, "%s\\ijat_na%d_det1_nsubdet%d", SRFS_DIR, na, nsubdet); fp1=fopen(filename, "wb"); 
						printf("<aperatt_ma_proj_sub_det>: writing srf_%d to %s\n", na, filename);
						ntmp2=1;
						for(i=1; i<=msize; i++) {
							ntmp2=fwrite(&(ijat[i]), sizeof(unsigned long), 1, fp1);
							if(ntmp2!=1) {
								printf("<aperatt_ma_proj_sub_det>: error in writing ijat, na=%d and nsubdet=%d\n", na, nsubdet); 						
								fclose(fp1);

								ntmp2=0;
								printf("<aperatt_ma_proj_sub_det.c>: since file %s is not properly finished, it is deleted\n", filename);
								if( remove(filename) == -1 )
									printf("<aperatt_ma_proj_sub_det.c>: error 11: delete error: %s\n", filename);
								//getchar();
								goto restart;
							}
						}
						fclose(fp1);
					}

					for (m=1; m<=4; m++) {
						if (m==1) sprintf(filename, "%s\\srfs\\sat_na%d_det1_nsubdet%d", ROOT_DIR1, na, nsubdet); 
						if (m==2) sprintf(filename, "%s\\srfs\\sat_na%d_det1_nsubdet%d", ROOT_DIR2, na, nsubdet); 
						if (m==3) sprintf(filename, "%s\\srfs\\sat_na%d_det1_nsubdet%d", ROOT_DIR3, na, nsubdet); 
						if (m==4) sprintf(filename, "%s\\srfs\\sat_na%d_det1_nsubdet%d", ROOT_DIR4, na, nsubdet); 

						fp1=fopen(filename, "wb"); 
						//sprintf(filename, "%s\\sat_na%d_det1_nsubdet%d", SRFS_DIR, na, nsubdet); fp1=fopen(filename, "wb"); 
						printf("<aperatt_ma_proj_sub_det>: writing srf_%d to %s\n", na, filename);
						ntmp1=1;
						for(i=1; i<=msize; i++) {
							ntmp1=fwrite(&(sat[i]), sizeof(float), 1, fp1);
							if(ntmp1!=1) {
								printf("<aperatt_ma_proj_sub_det>: error in writing sat, na=%d and nsubdet=%d\n", na, nsubdet); //getchar();
								ntmp1=0;
								printf("<aperatt_ma_proj_sub_det.c>: since file %s is not properly finished, it is deleted\n", filename);
								if( remove(filename) == -1 )
									printf("<aperatt_ma_proj_sub_det.c>: error 11: delete error: %s\n", filename);
								//getchar();
								goto restart;
							}
						}
						fclose(fp1);
					}
				}
				else if (isave==3) {	//isave=3, rfs is saved in the disk determined by the file: path.txt in COMMON_DIR\\data
					sprintf(filename2, "%s\\data\\path.txt", COMMON_DIR);
					fp1=fopen(filename2, "r"); 
					printf("read path from %s\n", filename2);
					for (i=1; i<=NANG*NSUBDET; i++) {
						fscanf(fp1, "%d", &(nna));
						fscanf(fp1, "%d", &(nnsubdet));
						fscanf(fp1, "%s", srfs_path);
						if ((nna==na)&&(nnsubdet==nsubdet)) { 
							printf("srfs should be saved in %s\n", srfs_path);
							break;
						}
					}
					fclose(fp1);

					if (0) {
						if (srfs_path[8]=='1') sprintf(filename, "%s\\srfs\\sen_na%d_det1_nsubdet%d", ROOT_DIR1, na, nsubdet); 
						if (srfs_path[8]=='2') sprintf(filename, "%s\\srfs\\sen_na%d_det1_nsubdet%d", ROOT_DIR2, na, nsubdet);
						if (srfs_path[8]=='3') sprintf(filename, "%s\\srfs\\sen_na%d_det1_nsubdet%d", ROOT_DIR3, na, nsubdet);
						if (srfs_path[8]=='4') sprintf(filename, "%s\\srfs\\sen_na%d_det1_nsubdet%d", ROOT_DIR4, na, nsubdet);
						if (srfs_path[8]=='5') sprintf(filename, "%s\\srfs\\sen_na%d_det1_nsubdet%d", ROOT_DIR5, na, nsubdet); 
						if (srfs_path[8]=='6') sprintf(filename, "%s\\srfs\\sen_na%d_det1_nsubdet%d", ROOT_DIR6, na, nsubdet);
						if (srfs_path[8]=='7') sprintf(filename, "%s\\srfs\\sen_na%d_det1_nsubdet%d", ROOT_DIR7, na, nsubdet);
						if (srfs_path[8]=='8') sprintf(filename, "%s\\srfs\\sen_na%d_det1_nsubdet%d", ROOT_DIR8, na, nsubdet);
					}
					else {
						if (nph==192) sprintf(filename, "%s\\srfs\\det1\\sen_na%d_det1_nsubdet%d", COMMON_DIR, na, nsubdet);
						else if (nph==7) sprintf(filename, "%s\\srfs\\det2\\sen_na%d_det1_nsubdet%d", COMMON_DIR, na, nsubdet);
					}
					fp3=fopen(filename, "wb"); 
					while (fp3==NULL) {
						Sleep(1000);
						printf("."); 
						fp3=fopen(filename, "wb"); 
					}
					printf("<aperatt_ma_proj_sub_det>: writing sensitivity to file %s ... ", filename);
					fwrite(sen, sizeof(float), NIMG,fp3);
					fclose(fp3);
					printf("finished!\n");
					/** calculating projection parameters (end) **/

					msize=ijat[ijat[1]-1]-1;
					printf("<aperatt_ma_proj_sub_det>: msize: %ld, %ld\n", msize, ijat[1]-2);
					
					if (nph==192) {
						if (srfs_path[8]=='1') sprintf(filename, "%s\\srfs\\det1\\ijat_na%d_det1_nsubdet%d", ROOT_DIR1, na, nsubdet); 
						if (srfs_path[8]=='2') sprintf(filename, "%s\\srfs\\det1\\ijat_na%d_det1_nsubdet%d", ROOT_DIR2, na, nsubdet);
						if (srfs_path[8]=='3') sprintf(filename, "%s\\srfs\\det1\\ijat_na%d_det1_nsubdet%d", ROOT_DIR3, na, nsubdet);
						if (srfs_path[8]=='4') sprintf(filename, "%s\\srfs\\det1\\ijat_na%d_det1_nsubdet%d", ROOT_DIR4, na, nsubdet);
						if (srfs_path[8]=='5') sprintf(filename, "%s\\srfs\\det1\\ijat_na%d_det1_nsubdet%d", ROOT_DIR5, na, nsubdet); 
						if (srfs_path[8]=='6') sprintf(filename, "%s\\srfs\\det1\\ijat_na%d_det1_nsubdet%d", ROOT_DIR6, na, nsubdet);
						if (srfs_path[8]=='7') sprintf(filename, "%s\\srfs\\det1\\ijat_na%d_det1_nsubdet%d", ROOT_DIR7, na, nsubdet);
						if (srfs_path[8]=='8') sprintf(filename, "%s\\srfs\\det1\\ijat_na%d_det1_nsubdet%d", ROOT_DIR8, na, nsubdet);
					}
					else if (nph==7) {
						if (srfs_path[8]=='1') sprintf(filename, "%s\\srfs\\det2\\ijat_na%d_det1_nsubdet%d", ROOT_DIR1, na, nsubdet); 
						if (srfs_path[8]=='2') sprintf(filename, "%s\\srfs\\det2\\ijat_na%d_det1_nsubdet%d", ROOT_DIR2, na, nsubdet);
						if (srfs_path[8]=='3') sprintf(filename, "%s\\srfs\\det2\\ijat_na%d_det1_nsubdet%d", ROOT_DIR3, na, nsubdet);
						if (srfs_path[8]=='4') sprintf(filename, "%s\\srfs\\det2\\ijat_na%d_det1_nsubdet%d", ROOT_DIR4, na, nsubdet);
						if (srfs_path[8]=='5') sprintf(filename, "%s\\srfs\\det2\\ijat_na%d_det1_nsubdet%d", ROOT_DIR5, na, nsubdet); 
						if (srfs_path[8]=='6') sprintf(filename, "%s\\srfs\\det2\\ijat_na%d_det1_nsubdet%d", ROOT_DIR6, na, nsubdet);
						if (srfs_path[8]=='7') sprintf(filename, "%s\\srfs\\det2\\ijat_na%d_det1_nsubdet%d", ROOT_DIR7, na, nsubdet);
						if (srfs_path[8]=='8') sprintf(filename, "%s\\srfs\\det2\\ijat_na%d_det1_nsubdet%d", ROOT_DIR8, na, nsubdet);
					}
					fp1=fopen(filename, "wb"); 
					
					while (fp1==NULL) {
						Sleep(1000);
						printf("."); 
						fp1=fopen(filename, "wb"); 
					}
					//sprintf(filename, "%s\\ijat_na%d_det1_nsubdet%d", SRFS_DIR, na, nsubdet); fp1=fopen(filename, "wb"); 
					printf("<aperatt_ma_proj_sub_det>: writing srf_%d to %s\n", na, filename);
					ntmp2=1;
					for(i=1; i<=msize; i++) {
						ntmp2=fwrite(&(ijat[i]), sizeof(unsigned long), 1, fp1);
						if(ntmp2!=1) {
							printf("<aperatt_ma_proj_sub_det>: error in writing ijat, na=%d and nsubdet=%d\n", na, nsubdet); 						
							fclose(fp1);

							ntmp2=0;
							printf("<aperatt_ma_proj_sub_det.c>: since file %s is not properly finished, it is deleted\n", filename);
							if( remove(filename) == -1 )
								printf("<aperatt_ma_proj_sub_det.c>: error 11: delete error: %s\n", filename);
							//getchar();
							goto restart;
						}
					}
					fclose(fp1);

					if (nph==192) {
						if (srfs_path[8]=='1') sprintf(filename, "%s\\srfs\\det1\\sat_na%d_det1_nsubdet%d", ROOT_DIR1, na, nsubdet); 
						if (srfs_path[8]=='2') sprintf(filename, "%s\\srfs\\det1\\sat_na%d_det1_nsubdet%d", ROOT_DIR2, na, nsubdet);
						if (srfs_path[8]=='3') sprintf(filename, "%s\\srfs\\det1\\sat_na%d_det1_nsubdet%d", ROOT_DIR3, na, nsubdet);
						if (srfs_path[8]=='4') sprintf(filename, "%s\\srfs\\det1\\sat_na%d_det1_nsubdet%d", ROOT_DIR4, na, nsubdet);
						if (srfs_path[8]=='5') sprintf(filename, "%s\\srfs\\det1\\sat_na%d_det1_nsubdet%d", ROOT_DIR5, na, nsubdet); 
						if (srfs_path[8]=='6') sprintf(filename, "%s\\srfs\\det1\\sat_na%d_det1_nsubdet%d", ROOT_DIR6, na, nsubdet);
						if (srfs_path[8]=='7') sprintf(filename, "%s\\srfs\\det1\\sat_na%d_det1_nsubdet%d", ROOT_DIR7, na, nsubdet);
						if (srfs_path[8]=='8') sprintf(filename, "%s\\srfs\\det1\\sat_na%d_det1_nsubdet%d", ROOT_DIR8, na, nsubdet);
					}
					else if (nph==7) {
						if (srfs_path[8]=='1') sprintf(filename, "%s\\srfs\\det2\\sat_na%d_det1_nsubdet%d", ROOT_DIR1, na, nsubdet); 
						if (srfs_path[8]=='2') sprintf(filename, "%s\\srfs\\det2\\sat_na%d_det1_nsubdet%d", ROOT_DIR2, na, nsubdet);
						if (srfs_path[8]=='3') sprintf(filename, "%s\\srfs\\det2\\sat_na%d_det1_nsubdet%d", ROOT_DIR3, na, nsubdet);
						if (srfs_path[8]=='4') sprintf(filename, "%s\\srfs\\det2\\sat_na%d_det1_nsubdet%d", ROOT_DIR4, na, nsubdet);
						if (srfs_path[8]=='5') sprintf(filename, "%s\\srfs\\det2\\sat_na%d_det1_nsubdet%d", ROOT_DIR5, na, nsubdet); 
						if (srfs_path[8]=='6') sprintf(filename, "%s\\srfs\\det2\\sat_na%d_det1_nsubdet%d", ROOT_DIR6, na, nsubdet);
						if (srfs_path[8]=='7') sprintf(filename, "%s\\srfs\\det2\\sat_na%d_det1_nsubdet%d", ROOT_DIR7, na, nsubdet);
						if (srfs_path[8]=='8') sprintf(filename, "%s\\srfs\\det2\\sat_na%d_det1_nsubdet%d", ROOT_DIR8, na, nsubdet);
					}
					fp1=fopen(filename, "wb"); 
					while (fp1==NULL) {
						Sleep(1000);
						printf("."); 
						fp1=fopen(filename, "wb"); 
					}
					//sprintf(filename, "%s\\sat_na%d_det1_nsubdet%d", SRFS_DIR, na, nsubdet); fp1=fopen(filename, "wb"); 
					printf("<aperatt_ma_proj_sub_det>: writing srf_%d to %s\n", na, filename);
					ntmp1=1;
					for(i=1; i<=msize; i++) {
						ntmp1=fwrite(&(sat[i]), sizeof(float), 1, fp1);
						if(ntmp1!=1) {
							printf("<aperatt_ma_proj_sub_det>: error in writing sat, na=%d and nsubdet=%d\n", na, nsubdet); //getchar();
							ntmp1=0;
							printf("<aperatt_ma_proj_sub_det.c>: since file %s is not properly finished, it is deleted\n", filename);
							if( remove(filename) == -1 )
								printf("<aperatt_ma_proj_sub_det.c>: error 11: delete error: %s\n", filename);
							//getchar();
							goto restart;
						}
					}
					fclose(fp1);
				}

				ifinished=0;
				if((ntmp1==1) && (ntmp2==1)) {
					//if(ntmp3==1) {
					sprintf(filename2, "%s\\srfs\\file_MC_finished", COMMON_DIR);
					fp2=fopen(filename2, "ab");  
					fwrite(&(na), sizeof(int), 1, fp2);	
					fwrite(&(nsubdet), sizeof(int), 1, fp2);
					fclose(fp2);
					ifinished=1;
				}

				if ((istart==1) && (ifinished!=1)) {
					printf("<aperatt_ma_proj_sub_det>: error 1: started but not finished ...\n");
					sprintf(filename2, "%s\\srfs\\file_MC_started_not_finished", COMMON_DIR);
					fp2=fopen(filename2, "ab");  
					fwrite(&(na), sizeof(int), 1, fp2);	
					fwrite(&(nsubdet), sizeof(int), 1, fp2);
					fclose(fp2);
					//getchar();
					goto restart;
				}
				goto restart;
next:
				continue;
			}
		}
	}
	/** now calculate the list-mode srf based on the events generated above (end)  **/
	
	free_vector(sat, 1, NMAX);
	free_lvector(ijat, 1, NMAX);	
	free_vector(det_prob, 0, NDET_DET1-1);
	free_vector(det_mask, 0, NDET_DET1-1);
	free_vector(det_pix_mask, 0, NDET_DET1-1);
	free_vector(sou_mask, 0, NIMG-1);

	free_vector(proj, 0, NDET_DET1-1);
	free_vector(proj_buf, 0, NDET_DET1-1);
	free_lvector(dn_buf, 0, NDET_DET1-1);
	free_vector(source, 0, NIMG-1);
	free_vector(sen, 0, NIMG-1);
	free_vector(img, 0, NIMG-1);
	for (i=0; i<nph; i++) 	free_vector(map[i], 0, N1*N2*N3*N4-1);
	free_vector(dp, 0, NDET_DET1*DM_DET *3-1);
	free_vector(sp, 0, NIMG*3-1);
    free_vector(compx,0,NPH_DET1); // center position for each compound eye;
	free_vector(compy,0,NPH_DET1);
	free_vector(radius_cp,0,NPH_DET1);
	free_vector(aaa,0,10*NPH_DET1+16);
	return iistart;
}

/*float pixelprob_det1(float *buf, float *dpos)
{
	float cx, cy, cr, dx, dy, dz, d, x, y, pr, a, fn, ddx, ddy, dxi, dyi;
	int i, j;	
	cx=buf[0]; cy=buf[1]; cr=buf[2]; 
	dx=dpos[0]; dy=dpos[1]; 
	
	//printf("cx=%f, cy=%f, cr=%f, dx=%f, dy=%f DDX_DET1=%f\n", cx, cy, cr, dx, dy, DDX_DET1); getchar();
	/** pre-exclude if outside a square zone **/
/*	if((dx>(cx+cr+DDX_DET1*0.5))||(dx<(cx-cr-DDX_DET1*0.5))) return 0.;
	if((dy>(cy+cr+DDY_DET1*0.5))||(dy<(cy-cr-DDY_DET1*0.5))) return 0.;
	
	/** calculating the prob for totally inculded pixels **/
/*	a=DDX_DET1*DDY_DET1; //size of a pixel
	d=sqrt((cx-dx)*(cx-dx)+(cy-dy)*(cy-dy));
	if(d<(cr-DDX_DET1*sqrt(2.))) {
		//printf("<pixelprob_det1>: total inclusion, prob: %.12f\n", a/(PI*cr*cr)*buf[3]);
		//printf("<pixelprob_det1>: %f %f %f %f\n", a, cr, (PI*cr*cr), buf[3]);
		return a/(PI*cr*cr)*buf[3];
	}

	/** calculating the prob for partially inculded pixels (start) **/	
/*	pr=0.;
	fn=4.;
	ddx=DDX_DET1/fn;
	ddy=DDY_DET1/fn;
	dxi=-DDX_DET1*0.5+DDX_DET1*0.5/fn+dx; // center of the lower left sub-pixel
	dyi=-DDY_DET1*0.5+DDY_DET1*0.5/fn+dy;
	for(i=0; i<fn; i++)
		for(j=0; j<fn; j++) {
			x=i*ddx+dxi;
			y=j*ddy+dyi;
			d=sqrt((cx-x)*(cx-x)+(cy-y)*(cy-y));
			if(d<cr) pr+=a;
		}
	pr=pr/(fn*fn)/(PI*cr*cr)*buf[3];
	return pr;
	/** calculating the prob for partially inculded pixels (end) **/	
/*}

float pixelprob_det1NL(float *buf, float *dpos, float *spos, float *map, float *dirt, int *index )
{
	float a, d, sc, x, y, pr, fn, ddx, ddy, dxi, dyi;
	float sx, sy, sz;	// source position
	float dx, dy, dz;   // detector pixel position
	float cx, cy;		// projection of center of ph in the detector surface
	float DDXph, DDYph;	// shift displacement in subpixel in ph plane coordinates
	float h, b;		
	float theta, cos_beta_phi;	//theta cos(beta-phi)
	int i, j, k, l, ii, jj;	
	float t1, t2, t3, dirt_L, cr, r2, rho, f, ll, Rlimit, alpha;

	/**  defination of parameters (start)  **/
/*	r2=buf[5];	Rlimit=buf[9]; 
	DDXph=buf[7]; DDYph=buf[8];   
	dx=dpos[0]; dy=dpos[1]; dz=dpos[2];
	sx=spos[0]; sy=spos[1]; sz=spos[2];
	j=index[0]; l=index[1]; ll=buf[12]; 
	theta=buf[6]; alpha=buf[10];
	t3=buf[13];

	h=buf[11];	// distance between the source and ph surface
	b=fabs(dz);	// distance between the detector pixel and ph surface
	wher
	
	cx=sx+(-sx)*sc;		// x of the center of the projection 
	cy=sy+(-sy)*sc;		// y of the center of the projection
	cr=r2*sc;			// radius of the projection of pinhole
	/**  defination of parameters (end)  **/

	/** pre-exclude if ray cannot pass through the pinhole **/
/*	if ((theta>alpha)||(theta==alpha)) {		
		return 0;
	}
	
	/** pre-exclude if outside the rho limit **/
/*	a=DDXph*DDYph;
	d=sqrt((dx-cx)*(dx-cx)+(dy-cy)*(dy-cy));
	if (d>Rlimit*sc){	
		return 0;
	}	

	/** calculating the prob for totally inculded pixels **/
/*	if(d<(cr-sqrt((DDXph*DDXph)+(DDYph*DDYph)))) {
		return a/(PI*cr*cr)*buf[3];
	}
	
	/** calculating the prob for partially inculded pixels (start) **/	
/*	pr=0.;
	fn=4.;
	ddx=DDXph/fn;
	ddy=DDYph/fn;
	dxi=-DDXph*0.5+DDXph*0.5/fn+dx; // center of the lower left sub-pixel
	dyi=-DDYph*0.5+DDYph*0.5/fn+dy;
	
	for(ii=0; ii<fn; ii++) {
		for(jj=0; jj<fn; jj++) {
			x=ii*ddx+dxi;
			y=jj*ddy+dyi;
			
			d=sqrt((x-cx)*(x-cx)+(y-cy)*(y-cy));
			t1=d/sc;
			t2=sqrt((x-sx)*(x-sx)+(y-sy)*(y-sy))/sc;

			if ((t1 > Rlimit)||(t1==Rlimit)) { pr=pr;}	// outside of the rho limit
			else {
				if ((d <cr)||(d==cr)) { pr= pr+a;}		// totally included subpixels
				else {									// partially included subpixels
			
					cos_beta_phi=(t1*t1+t3*t3-t2*t2)/2/t1/t3;
					
					i=(int)((t1-r2)/dirt[0]);
					k=(int)((cos_beta_phi+1)/dirt[2]);

					if ((i==N1)||(j==N2)||(k==N3)||(l==N4)||(i<0)||(j<0)||(k<0)||(l<0)||(i>N1)||(j>N2)||(k>N3)||(l>N4)){
						printf("\n !!!! error!  i=%d, j=%d, k=%d, l=%d\n", i, j, k, l);
						printf("t1=%f, t2=%f t3=%f, cos_beta_phi=%f\n", t1, t2, t3, cos_beta_phi); 
						printf("x=%f, y=%f, sx=%f, sy=%f sz %f\n", x, y, sx, sy, sz);
						printf("dirt[3] %10.10f dirt[2] %10.10f dirt[1] %10.10f dirt[0] %10.10f",dirt[3],dirt[2],dirt[1],dirt[0]); 
						printf(" sc %f d %f t1 %f           dx=%f, dy=%f dz %f\n",sc,d,t1, dx, dy, dz);
						printf(" ll %f",ll);
						getchar();

					}

					f=(ll-l)*map[i*M1+j*M2+k*N4+l]+(l+1-ll)*map[i*M1+j*M2+k*N4+l+1];

					pr=pr+a*f;
					
				}
			}
		}
	}

	pr = pr/fn/fn/(PI*cr*cr)*buf[3];
	return pr;
	/** calculating the prob for partially inculded pixels (end) **/	

//}
