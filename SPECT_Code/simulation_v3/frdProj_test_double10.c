#include "headerFile"


int frdProj_test_double10(double PSys[],double PDet[],double *PAttMap[],double Source[],double Proj[])
{
	double sx,sy,sz,spos[3]; /* source pixel postion in the global geometry*/
	/* Detector parameters setup start*/
	time_t start, finish;
	double duration;
	int DN_DET1,DN_DETX1,DN_DETY1,DM_DET,NDET_DET1,fn,IDETPOS,nsubdet;
	double ddx,DDX_DET1,ddy,DDY_DET1,ddz;
	double MU_DET;/* attenuation coeffention of detect*/
	double *dp,*dsubp,*subdpos,*det_pix_mask,*dpos;
	double detX[3],detY[3],detZ[3];/* detector norm vector in global coordination*/
	double dcx,dcy,dcz; /*/ detector center position*/
    /*/ detector parameters end*/

    /*/ pinhole parameters start*/
	int nph; /*/ numbers of pinhole;*/
	double r2,cht,ct,APERALPHA; /*/ the pinhole structure parameters*/
    double  *hx,*hy,*cz; /*/ pointer of vector for pinhole center postion in the global coordination;*/
	double Rlimit; /*/ projection limite*/
	double *detAXx, *detAXy, *detAXz,*detAYx, *detAYy, *detAYz,*detAZx, *detAZy, *detAZz;/*/ pointer of vector for Pinhole num in the global geometry*/
	double *DDX_ph,*DDY_ph;
	/*/ pinhole parameters end*/

	/*/ other variables*/
	int npix; /*/return value:the number of pixel, whose projection probabilities are larger than thresh.*/
	double thresh;/*/ the threshold of probabilities*/
	double pr; /*/ the probabilities of each pixel.*/
	int n,i,j,k,dpi,ndy,ndx,deltan,index111;/*/ loop index*/
	int ntmp3,ntmp4,ntmp5,ntmp6,ntmp7,ntmp8;
	int dn,dncx,dncy,dtmpx,dtmpy;
	long int nt,nm;
	double dsq,hs,hsa,hsp,f,nsrx,nsry,nsrz,dpx,dpy,dp1,dp2;
	double cosTheta,cosThetaD,sinThetaD, *cos_D;
	double buf[16];
	double tmp,tmpdx,tmpdy,tmpdz,*N_tmp;
	int N_step[4];/*/ the vector to store the N1,N2,N3,N4 of map;*/
	/*/*/
	time(&start);
    sx=Source[0],sy=Source[1], sz=Source[2];/*/ soure position in the global coordinat */
	/*/ detector parameters inializing start*/
	DN_DET1=(int)(PSys[0]);/*/ numbers of Detector in one demension*/
	DM_DET=(int)(PSys[1]);/*/ number of detector layers*/
	nph=(int)(PSys[17]);/*/ numbers of pinhole;*/
	DN_DETX1=DN_DET1;
	DN_DETY1=(int)(PSys[23+nph*12]);
	NDET_DET1=DN_DETX1*DN_DETY1;/*/ numbers of detector pixel;*/


	fn=(int)(PSys[2]);/*/ detector pixel further devided number*/
	IDETPOS=(int)(PSys[3]); /*/if IDETPOS=0, there is no gap between the pixel;*/
	nsubdet=(int)(PDet[0]);/*/ subdetector index;*/
	ddx=PSys[13];/* detector pix size;*/
	DDX_DET1=ddx;/*/ detector pix size;*/
	ddy=PSys[14];/*/ detector pix size;*/
	DDY_DET1=ddy;/*/ detector pix size;*/
	ddz=PSys[15];/*/ detector pix size;*/
	MU_DET=PSys[16];/*/ attenuation coeffention of detect*/
	dp=dvector(0,3*NDET_DET1*DM_DET-1);/*/ detector pixel postion in Global coordination*/
	dpos=dvector(0,3*DM_DET-1);
    dsubp=dvector(0,DM_DET*NDET_DET1*3*fn*fn-1);/*/ detector pixel are futhered divided by fn*fn; this vector shows each subpixel's position in global geometry.*/
	subdpos=dvector(0,fn*fn*DM_DET*3-1);/*/ detector subpixel in phinhole geometry.//reivised by Guanfeng*/
	det_pix_mask=dvector(0,DM_DET*NDET_DET1-1);/*/ subdetector index and  if value==0, then it is a defect detector*/
	for(i=0;i<3;i++){/*/ detector norm vector in the global coordination;*/
		detX[i]=PSys[4+i];
		detY[i]=PSys[7+i];
		detZ[i]=PSys[10+i];
	}
	for(dn=0;dn<NDET_DET1*DM_DET;dn++){
	det_pix_mask[dn]=(PDet[4*dn+1]);/*/ subdetector index or bad detector index(==0)*/
	dp[3*dn]=PDet[4*dn+2];/*/ detector position in x direction in global coordination;*/
	dp[3*dn+1]=PDet[4*dn+3];/*/detector position in y direction in global coordination;*/
	dp[3*dn+2]=PDet[4*dn+4];/*/detector position in z direction in global coordination;*/
	/*/mexPrintf("dn %d dp %f %f %f\n",dn, dp[3*dn],dp[3*dn+1],dp[3*dn+2]);*/
	/*/mexPrintf("%f %f %f\n",PDet[dn*4+2],PDet[4*dn+3],PDet[4*dn+4]);*/
	/*/mexPrintf("%f %f\n",det_pix_mask[dn],PDet[4*dn+1]);*/
	}

	dcx=(dp[0]+dp[3*NDET_DET1*DM_DET-3])/2.0;/*/ detector center position in x direction for global coordination. /*/
	dcy=(dp[1]+dp[3*NDET_DET1*DM_DET-2])/2.0;/*detector center position in y direction /*/
	dcz=(dp[2]+dp[3*NDET_DET1*DM_DET-1])/2.0;/*detector center position in z direction /*/
	/*mexPrintf("dcx %f dcy %f dcz %f\n",dcx,dcy,dcz);return;/*/
	/* Detector parameters inistialize end/*/

	/* Pinhole parameters initializeing start/*/
	nph=(int)(PSys[17]);/* numbers of pinhole/*/
	hx=dvector(0,nph-1);/* pinholes center position in the global geometry(x)/*/
	hy=dvector(0,nph-1);/* pinholes center postion in the global geometry(y)/*/
	cz=dvector(0,nph-1);/* pinhloes center postion in the global geometry(z)/*/
	Rlimit=PSys[22+12*nph];/* calculation limit in pinhole/*/
	detAXx=dvector(0,nph-1),detAXy=dvector(0,nph-1),detAXz=dvector(0,nph-1);/* Pinhole local X-axis in global coordination/*/
	detAYx=dvector(0,nph-1),detAYy=dvector(0,nph-1),detAYz=dvector(0,nph-1);/* Pihhole local Y-axis in global coordination/*/
	detAZx=dvector(0,nph-1),detAZy=dvector(0,nph-1),detAZz=dvector(0,nph-1);/* Pinhole local Z-axis in global coordination/*/
	DDX_ph=dvector(0,nph-1),DDY_ph=dvector(0,nph-1); /* the detector pixel size projection to pinhole;/*/
	r2=PSys[18+12*nph],cht=PSys[18+12*nph+1],ct=PSys[18+12*nph+2],APERALPHA=PSys[18+12*nph+3];/*the pinhole structure parameters/*/
	for(k=0;k<nph;k++){
	hx[k]=PSys[18+k*12];
	hy[k]=PSys[18+k*12+1];
	cz[k]=PSys[18+k*12+2];
	detAXx[k]=PSys[21+k*12];
	detAXy[k]=PSys[21+k*12+1];
	detAXz[k]=PSys[21+k*12+2];
	detAYx[k]=PSys[24+k*12];
	detAYy[k]=PSys[24+k*12+1];
	detAYz[k]=PSys[24+k*12+2];
    detAZx[k]=PSys[27+k*12];
	detAZy[k]=PSys[27+k*12+1];
	detAZz[k]=PSys[27+k*12+2];

	tmp=detAZx[k]*detX[0]+detAZy[k]*detX[1]+detAZz[k]*detX[2];
	DDX_ph[k]=ddx*sqrt(1-tmp*tmp);
	tmp=detAZx[k]*detY[0]+detAZy[k]*detY[1]+detAZz[k]*detY[2];
	DDY_ph[k]=ddy*sqrt(1-tmp*tmp);
	}

	/* Pinhole Parameters setup ended/*/
	/* other parameters initializeing/*/
	npix=0;/* return value: the number of pixel, whose probabilities are larger than threshold/*/
	pr=0;
	thresh=1.e-12;
	cos_D=dvector(0,DM_DET-1);/* incident angle;/*/
	for(n=0;n<15;n++){
		buf[n]=0;
	}
	/*mexPrintf("%f\d",PAttMap[i][0]);/*/
	N_tmp=PAttMap[nph];/* the pointer to the vector, which stores N1,N2,N3,N4 of each map;/*/
    index111=0;

	for(k=0; k<nph; k++) {  /*printf("%d k\n",k);/*/
							spos[0]=(sx-hx[k])*detAXx[k]+(sy-hy[k])*detAXy[k]+(sz-cz[k])*detAXz[k]; /*x position of source pixel in pinhole geometry /*/
							spos[1]=(sx-hx[k])*detAYx[k]+(sy-hy[k])*detAYy[k]+(sz-cz[k])*detAYz[k];
							spos[2]=(sx-hx[k])*detAZx[k]+(sy-hy[k])*detAZy[k]+(sz-cz[k])*detAZz[k];/*// === CHANGE ===*/
							spos[1]=(sx-hx[k])*detAYx[k]+(sy-hy[k])*detAYy[k]+(sz-cz[k])*detAYz[k];
							spos[2]=(sx-hx[k])*detAZx[k]+(sy-hy[k])*detAZy[k]+(sz-cz[k])*detAZz[k];/* === CHANGE ===//*/
                            N_step[0]=(int)N_tmp[0+4*k];N_step[1]=(int)(N_tmp[1+4*k]);N_step[2]=(int)(N_tmp[2+4*k]);N_step[3]=(int)(N_tmp[3+4*k]);/* N1,N2,N3,N4 for kth pinhole/*/
     						dsq=(sx-hx[k])*(sx-hx[k])+(sy-hy[k])*(sy-hy[k])+(sz-cz[k])*(sz-cz[k]); /*square of source-pinhole distance/*/
							hs=fabs(detAZx[k]*(sx-hx[k])+detAZy[k]*(sy-hy[k])+detAZz[k]*(sz-cz[k]));   /* distance source to pinhole plane/*/
							/*================================ n vector  of source&pinhole  /*/
							nsrx=(sx-hx[k])/sqrt(dsq);
							nsry=(sy-hy[k])/sqrt(dsq);
							nsrz=(sz-cz[k])/sqrt(dsq);
							cosTheta=fabs(hs/sqrt(dsq));
							cosThetaD=fabs(nsrx*detZ[0]+nsry*detZ[1]+nsrz*detZ[2]); /*incident angle to detector surface/*/
							sinThetaD=sqrt(1-cosThetaD*cosThetaD);

							buf[5]=r2;					/* radius of ph/*/

							if(cosTheta<=0.9999) buf[6]=acos(cosTheta);		/* angle between norm of ph plane and the line from source to pinhole center/*/
							if(buf[6]>PI/2.) buf[6]=PI-buf[6];  /*angle between source point and y axis/*/
							if(cosTheta>0.9999) buf[6]=0.;

							buf[7]=ddx;				/* projection of x dimension of detector pixel on the pinhole plane/*/
							buf[8]=ddy;				/* projection of y dimension of detector pixel on the pinhole plane/*/
							buf[9]=Rlimit;				/* calculation limit in rho/*/

							buf[10]=atan(tan(APERALPHA*PI/180));			/* half of acceptance angle of pinhole	*/
							buf[15]=r2-cht/2*tan(APERALPHA*PI/180);
							buf[11]=hs;					/* distance between source and ph plane/*/
							buf[13]=ddz;				/* detector layer thickness/*/
							hsa=fabs((sx-hx[k])*detZ[0]+(sy-hy[k])*detZ[1]+(sz-cz[k])*detZ[2]);/* distance bwtween source and ph on detector plane/*/
							ntmp3=DN_DETY1-1;ntmp4=0;ntmp5=DN_DETX1-1;ntmp6=0;
							buf[3]=(r2*r2)/(4*dsq)*cosTheta; /* prob of falling through a pinhole/*/
							dp2=1000;
							for(dpi=0; dpi<DM_DET; dpi++) {
							/**  Defination of buf and index for subroutine of penetration factor (start)  **/
								hsp=fabs(detZ[0]*(sx-dcx)+detZ[1]*(sy-dcy)+detZ[2]*(sz-dcz))+ddz*(dpi-(DM_DET-1)*.5); /* distance between source and detector plate/*/
								f=hsp/hsa;

								buf[0]=sx+(hx[k]-sx)*f;		/* x of the center of the projection/*/
								buf[1]=sy+(hy[k]-sy)*f;		/* y of the center of the projection/*/
								buf[4]=sz+(cz[k]-sz)*f;		/* z value of detector plane/*/
								buf[2]=f*(Rlimit);			/* radius of the projection using Rlimit for penetration factor/*/

								if (IDETPOS&&0){
									dtmpx=(int)(((buf[0]-dcx)*detX[0]+(buf[1]-dcy)*detX[1]+(buf[4]-dcz)*detX[2])/ddx+DN_DETX1*0.5);
									dtmpy=(int)(((buf[0]-dcx)*detY[0]+(buf[1]-dcy)*detY[1]+(buf[4]-dcz)*detY[2])/ddy+DN_DETY1*0.5);
									/*dp2=1000;/*/

									for(ndy=fmax(0,dtmpy-2); ndy<=fmin(dtmpy+2,DN_DETY1); ndy++) {
										for(ndx=fmax(0,dtmpx-2); ndx<=fmin(dtmpx+2, DN_DETX1); ndx++) {
											dn=dpi*NDET_DET1+DN_DETX1*ndy+ndx;
											dpx=(dp[dn*3]-buf[0])*detX[0]+(dp[dn*3+1]-buf[1])*detX[1]+(dp[dn*3+2]-buf[4])*detX[2];
											dpy=(dp[dn*3]-buf[0])*detY[0]+(dp[dn*3+1]-buf[1])*detY[1]+(dp[dn*3+2]-buf[4])*detY[2];
											dp1=dpx*dpx+dpy*dpy;

											if (dp1<dp2) {dncx=ndx;dncy=ndy;dp2=dp1;}
										}
									}

								}
								else {	/*??? dcz layer/*/
									dncx=(int)(((buf[0]-dcx)*detX[0]+(buf[1]-dcy)*detX[1]+(buf[4]-dcz)*detX[2])/ddx+DN_DETX1*0.5);
									dncy=(int)(((buf[0]-dcx)*detY[0]+(buf[1]-dcy)*detY[1]+(buf[4]-dcz)*detY[2])/ddy+DN_DETY1*0.5);
								}
								/*printf("%d %d\n", dncx, dncy);/*/
								deltan=(int)(fabs(buf[2]/cosTheta*cosThetaD/DDX_DET1)+3);
								ntmp3=fmin(dncy-deltan, ntmp3);
								ntmp4=fmax(dncy+deltan, ntmp4);
								ntmp5=fmin(dncx-deltan,ntmp5);
								ntmp6=fmax(dncx+deltan, ntmp6);
							}

							ntmp3=fmax(ntmp3, 0);
							ntmp3=fmin(ntmp3,DN_DETY1-1);

							ntmp4=fmin(ntmp4, DN_DETY1-1);
							ntmp4=fmax(ntmp4,0);

							ntmp5=fmax(ntmp5,0);
							ntmp5=fmin(ntmp5, DN_DETY1-1);

							ntmp6=fmin(ntmp6, DN_DETX1-1);
							ntmp6=fmax(ntmp6,0);


							for(ndy=ntmp3; ndy<=ntmp4; ndy++) {
								for(ndx=ntmp5; ndx<=ntmp6; ndx++) {

									dn=DN_DETX1*ndy+ndx;

									if ((det_pix_mask[dn]==nsubdet)||(det_pix_mask[dn]==100+nsubdet)) {

										nt=0;
										for(dpi=0; dpi<DM_DET; dpi++) {
											ntmp7=(dn+dpi*NDET_DET1)*3;ntmp8=3*dpi;
											dpos[0+ntmp8]=(dp[ntmp7]-hx[k])*detAXx[k]+(dp[ntmp7+1]-hy[k])*detAXy[k]+(dp[ntmp7+2]-cz[k])*detAXz[k]; /*x position of source pixel in pinhole geometry /*/
											dpos[1+ntmp8]=(dp[ntmp7]-hx[k])*detAYx[k]+(dp[ntmp7+1]-hy[k])*detAYy[k]+(dp[ntmp7+2]-cz[k])*detAYz[k];
											dpos[2+ntmp8]=(dp[ntmp7]-hx[k])*detAZx[k]+(dp[ntmp7+1]-hy[k])*detAZy[k]+(dp[ntmp7+2]-cz[k])*detAZz[k];
											nsrx=(sx-dp[ntmp7]);
											nsry=(sy-dp[ntmp7+1]);
											nsrz=(sz-dp[ntmp7+2]);
											cos_D[dpi]=fabs(nsrx*detZ[0]+nsry*detZ[1]+nsrz*detZ[2])/sqrt(nsrx*nsrx+nsry*nsry+nsrz*nsrz); /* for each pixel/*/
											nm=0;
											for (j=0;j<fn;j++){
												for (i=0;i<fn;i++){
													tmpdx=dp[0+ntmp7]+(i-(fn-1)/2.0)*ddx/fn*detX[0]+(j-(fn-1)/2.0)*ddy/fn*detY[0]+0*detZ[0];
													tmpdy=dp[1+ntmp7]+(i-(fn-1)/2.0)*ddx/fn*detX[1]+(j-(fn-1)/2.0)*ddy/fn*detY[1]+0*detZ[1];
													tmpdz=dp[2+ntmp7]+(i-(fn-1)/2.0)*ddx/fn*detX[2]+(j-(fn-1)/2.0)*ddy/fn*detY[2]+0*detZ[2];

													subdpos[nt*3]=(tmpdx-hx[k])*detAXx[k]+(tmpdy-hy[k])*detAXy[k]+(tmpdz-cz[k])*detAXz[k]; /*x position of source pixel in pinhole geometry /*/
													subdpos[nt*3+1]=(tmpdx-hx[k])*detAYx[k]+(tmpdy-hy[k])*detAYy[k]+(tmpdz-cz[k])*detAYz[k];
													subdpos[nt*3+2]=(tmpdx-hx[k])*detAZx[k]+(tmpdy-hy[k])*detAZy[k]+(tmpdz-cz[k])*detAZz[k];
													nt++;nm++;
												}
											}
										}
										pr=pixelprob_det2ML_sum_ph_insert2_double10(buf, dpos, subdpos, spos, PAttMap[k], N_step, cos_D, DM_DET, MU_DET, fn);
										index111=index111+1;
										if(pr>thresh) {

											 if(det_pix_mask[dn]==nsubdet) {
												Proj[2*npix]=dn+1;/*dn start from zero, the actually index of detector of dn should plus 1/*/
												Proj[2*npix+1]=pr;
												npix++;
											 }
										}
									}/* det_mask/*/
								}/*ndy/*/
							}
}/* k nph/*/
						free_dvector(dp,0,3*NDET_DET1*DM_DET-1);
						free_dvector(dpos,0,3*DM_DET-1);
						free_dvector(dsubp,0,DM_DET*NDET_DET1*3*fn*fn-1);
						free_dvector(subdpos,0,fn*fn*DM_DET*3-1);/*revised by Guanfeng/*/
						free_dvector(det_pix_mask,0,NDET_DET1-1);
						free_dvector(cos_D,0,DM_DET-1);
						free_dvector(detAXx,0,nph-1);
						free_dvector(detAXy,0,nph-1);
						free_dvector(detAXz,0,nph-1);
						free_dvector(detAYx,0,nph-1);
						free_dvector(detAYy,0,nph-1);
						free_dvector(detAYz,0,nph-1);
						free_dvector(detAZx,0,nph-1);
						free_dvector(detAZy,0,nph-1);
						free_dvector(detAZz,0,nph-1);
						free_dvector(DDX_ph,0,nph-1);
						free_dvector(DDY_ph,0,nph-1);
						time(&finish);
						duration=1.0e+20*difftime (finish,start);
						return npix;
}
