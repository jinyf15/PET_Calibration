#include "headFile.h";

int genSubSRF(struct detector * det,int iDet,int iSubDet,struct source *sou,int iSouP,struct pinMap *ph,struct sysSRF *srf,int iThread)
{
    double sx,sy,sz,spos[3]; /* source pixel postion in the global geometry*/
	/* Detector parameters setup start*/
	time_t start, finish,time1;
	double duration;
	int DN_DETX1,DN_DETY1,DM_DET,NDET_DET1,fn,IDETPOS,nsubdet;
	double ddx,DDX_DET1,ddy,DDY_DET1,ddz;
	double MU_DET;/* attenuation coeffention of detect*/
	double *dp,*dsubp,*subdpos,*dpos;
	int * det_pix_mask;
	double *detX;
	double *detY;
	double *detZ;/* detector norm vector in global coordination*/
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

	// system response function parameter;
    double *proj;
	int *pixIndex;
    unsigned int spri;
    unsigned int sprj;
    unsigned int sprk;
    double *sen;
    unsigned int NIMG;
    unsigned int NSIZE;
    unsigned long int n;
    unsigned int l;


	/*/ other variables*/
	int npix; /*/return value:the number of pixel, whose projection probabilities are larger than thresh.*/
	double thresh;/*/ the threshold of probabilities*/
	double pr; /*/ the probabilities of each pixel.*/
	int i,j,k,dpi,ndy,ndx,deltan,index111;/*/ loop index*/
	int ntmp3,ntmp4,ntmp5,ntmp6,ntmp7,ntmp8;
	int dn,dncx,dncy,dtmpx,dtmpy;
	long int nt,nm;
	double dsq,hs,hsa,hsp,f,nsrx,nsry,nsrz,dpx,dpy,dp1,dp2;
	double cosTheta,cosThetaD,sinThetaD, *cos_D;
	double buf[16];
	double tmp,tmpdx,tmpdy,tmpdz,*N_tmp;
	int *N_step;/*/ the vector to store the N1,N2,N3,N4 of map;*/
	unsigned int sn; // index of source voxel;
	int phIndex;
	char fileName[1000];
	FILE *fp1;

	/*/*/
	time(&start);

	/*/ detector parameters inializing start*/

	DM_DET = det[iDet].NZ;/*/ number of detector layers*/
	nph = det[iDet].nphPerDet;/*/ numbers of pinhole;*/
	DN_DETX1 = det[iDet].NX;
	DN_DETY1 = det[iDet].NY;
	NDET_DET1 = DN_DETX1*DN_DETY1;/*/ numbers of detector pixel;*/


	fn = det[iDet].nDivide;/*/ detector pixel further devided number*/
	IDETPOS = 1; /*/if IDETPOS=0, there is no gap between the pixel;*/
	nsubdet = iSubDet;/*/ subdetector index;*/

	ddx = det[iDet].dx;/* detector pix size;*/
	DDX_DET1 = ddx;/*/ detector pix size;*/
	ddy = det[iDet].dy;/*/ detector pix size;*/
	DDY_DET1 = ddy;/*/ detector pix size;*/
	ddz = det[iDet].dz;/*/ detector pix size;*/

	MU_DET = det[iDet].attenCoef;/*/ attenuation coeffention of detect*/

	dp = det[iDet].dp;/*/ detector pixel postion in Global coordination*/
	dpos = (double *)malloc(sizeof(double)*3*DM_DET);
    //dsubp = (double *)malloc(sizeof(double)*DM_DET*NDET_DET1*3*fn*fn);/*/ detector pixel are futhered divided by fn*fn; this vector shows each subpixel's position in global geometry.*/
	subdpos = (double *)malloc(sizeof(double)*3*DM_DET*fn*fn); /*/ detector subpixel in phinhole geometry.//reivised by Guanfeng*/

    dcx = det[iDet].detCent[0]; /*/ detector center position in x direction for global coordination. /*/
	dcy = det[iDet].detCent[1];/*detector center position in y direction /*/
	dcz = det[iDet].detCent[2];/*detector center position in z direction /*/
	/*mexPrintf("dcx %f dcy %f dcz %f\n",dcx,dcy,dcz);return;/*/


	det_pix_mask = det[iDet].detPixMask;/*/ subdetector index and  if value==0, then it is a defect detector*/


    detX = det[iDet].detX;/*/ detector norm vector in the global coordination;*/
    detY = det[iDet].detY;
    detZ = det[iDet].detZ;



	/* Detector parameters inistialize end/*/

	/* Pinhole parameters initializeing start/*/
	nph = det[iDet].nphPerDet;/* numbers of pinhole/*/
	hx = (double *)malloc(sizeof(double)*nph);/* pinholes center position in the global geometry(x)/*/
	hy = (double *)malloc(sizeof(double)*nph);/* pinholes center postion in the global geometry(y)/*/
	cz = (double *)malloc(sizeof(double)*nph);/* pinhloes center postion in the global geometry(z)/*/


	detAXx = (double *)malloc(sizeof(double)*nph);
	detAXy = (double *)malloc(sizeof(double)*nph);
	detAXz = (double *)malloc(sizeof(double)*nph);/* Pinhole local X-axis in global coordination/*/
	detAYx = (double *)malloc(sizeof(double)*nph);
	detAYy = (double *)malloc(sizeof(double)*nph);
	detAYz = (double *)malloc(sizeof(double)*nph);/* Pihhole local Y-axis in global coordination/*/
	detAZx = (double *)malloc(sizeof(double)*nph);
	detAZy = (double *)malloc(sizeof(double)*nph);
	detAZz = (double *)malloc(sizeof(double)*nph);/* Pinhole local Z-axis in global coordination/*/

	DDX_ph = (double *)malloc(sizeof(double)*nph);
	DDY_ph = (double *)malloc(sizeof(double)*nph);

	//r2=PSys[18+12*nph],cht=PSys[18+12*nph+1],ct=PSys[18+12*nph+2],APERALPHA=PSys[18+12*nph+3];/*the pinhole structure parameters/*/

	for(k=0;k<nph;k++)
	{
        hx[k] = det[iDet].phCenter[k][0]; // pinhole center position
        hy[k] = det[iDet].phCenter[k][1];
        cz[k] = det[iDet].phCenter[k][2];

        detAXx[k] = det[iDet].phX[k][0]; // local x
        detAXy[k] = det[iDet].phX[k][1];
        detAXz[k] = det[iDet].phX[k][2];

        detAYx[k] = det[iDet].phY[k][0]; // local y
        detAYy[k] = det[iDet].phY[k][1];
        detAYz[k] = det[iDet].phY[k][2];

        detAZx[k] = det[iDet].phZ[k][0]; // local z
        detAZy[k] = det[iDet].phZ[k][1];
        detAZz[k] = det[iDet].phZ[k][2];


        tmp=detAZx[k]*detX[0]+detAZy[k]*detX[1]+detAZz[k]*detX[2];
        DDX_ph[k]=ddx*sqrt(1-tmp*tmp);
        tmp=detAZx[k]*detY[0]+detAZy[k]*detY[1]+detAZz[k]*detY[2];
        DDY_ph[k]=ddy*sqrt(1-tmp*tmp);
	}

	/* Pinhole Parameters setup ended/*/
	/* other parameters initializeing/*/

    thresh = 1.e-10;
    cos_D = (double *)malloc(sizeof(double)*DM_DET* fn);/*/incident angle;/*/

    proj = (double *)malloc(sizeof(double)*NDET_DET1);
    pixIndex = (int *)malloc(sizeof(int)*NDET_DET1);

    NIMG = sou[iSouP].NX * sou[iSouP].NY * sou[iSouP].NZ;
    sen = (double *)malloc(sizeof(double)*NIMG);

    for(dn = 0; dn < NDET_DET1; dn ++)
    {
        proj[dn] = 0;
        pixIndex[dn] = -10000;
    }

	for(sn=0; sn<NIMG;sn++)
	//for(sn=0; sn<1;sn++)
	{
        if(sn%(10000)==0&&iThread==0)
        {
            time(&(time1));
            printf("fininshed %u,total %u,srf size %u, maxiumusize %ld,time %s\n",sn,NIMG,sprk,srf[iThread].nMax,asctime(localtime(&time1)));
        }
        //printf("sn %d \n",sn);
        sx = sou[iSouP].sp[3*sn+0];
        sy = sou[iSouP].sp[3*sn+1];
        sz = sou[iSouP].sp[3*sn+2];
        //sx = 0; sy = 0; sz = 0;

        npix = 0;
        pr = 0;
        sen[sn] = 0;

        for(n=0;n<15;n++)
        {
            buf[n]=0;
        }


        /*mexPrintf("%f\d",PAttMap[i][0]);/*/
        //N_tmp = PAttMap[nph];/* the pointer to the vector, which stores N1,N2,N3,N4 of each map;/*/
        index111 = 0;

        for(k=0; k<nph; k++)
        {
            phIndex = det[iDet].phIndex[k]; // which pinhole is in pinMap for kth pinhole detector
            r2 = ph[phIndex].radius; // radius of pinhole
            cht = ph[phIndex].channel;
            ct = ph[phIndex].length;
            APERALPHA = ph[phIndex].openAngle*180/PI;/*half of total open angle 60degree/180*pi pinhole ,this should be 30*pi/180*/
            N_step = ph[phIndex].NStep; // number of step in each demension;
            Rlimit = ph[phIndex].Rlimit;


            spos[0] = (sx-hx[k])*detAXx[k] + (sy-hy[k])*detAXy[k] + (sz-cz[k])*detAXz[k]; /*x position of source pixel in pinhole geometry /*/
            spos[1] = (sx-hx[k])*detAYx[k] + (sy-hy[k])*detAYy[k] + (sz-cz[k])*detAYz[k];
            spos[2] = (sx-hx[k])*detAZx[k] + (sy-hy[k])*detAZy[k] + (sz-cz[k])*detAZz[k];/*// === CHANGE ===*/

            //spos[1] = (sx-hx[k])*detAYx[k] + (sy-hy[k])*detAYy[k] + (sz-cz[k])*detAYz[k];
            //spos[2] = (sx-hx[k])*detAZx[k] + (sy-hy[k])*detAZy[k] + (sz-cz[k])*detAZz[k];/* === CHANGE ===//*/


            dsq = (sx-hx[k])*(sx-hx[k]) + (sy-hy[k])*(sy-hy[k]) + (sz-cz[k])*(sz-cz[k]); /*square of source-pinhole distance/*/
            hs = fabs(detAZx[k]*(sx-hx[k]) + detAZy[k]*(sy-hy[k]) + detAZz[k]*(sz-cz[k]));   /* distance source to pinhole plane/*/
            /*================================ n vector  of source&pinhole  /*/
            nsrx = (sx-hx[k])/sqrt(dsq);
            nsry = (sy-hy[k])/sqrt(dsq);
            nsrz = (sz-cz[k])/sqrt(dsq);

            cosTheta = fabs(hs/sqrt(dsq));
            cosThetaD = fabs(nsrx*detZ[0]+nsry*detZ[1]+nsrz*detZ[2]); /*incident angle to detector surface/*/
            sinThetaD = sqrt(1-cosThetaD*cosThetaD);

            buf[5] = r2;					/* radius of ph/*/

            if(cosTheta<=0.9999) buf[6]=acos(cosTheta);		/* angle between norm of ph plane and the line from source to pinhole center/*/
            if(buf[6]>PI/2.) buf[6]=PI-buf[6];  /*angle between source point and y axis/*/

            if(cosTheta>0.9999) buf[6]=0.;

            buf[7] = ddx;				/* projection of x dimension of detector pixel on the pinhole plane/*/
            buf[8] = ddy;				/* projection of y dimension of detector pixel on the pinhole plane/*/
            buf[9] = Rlimit;				/* calculation limit in rho/*/

            buf[10] = atan(tan(APERALPHA*PI/180));			/* half of acceptance angle of pinhole	*/
            buf[15] = r2-cht/2*tan(APERALPHA*PI/180);
            buf[11] = hs;					/* distance between source and ph plane/*/
            buf[13] = ddz;				/* detector layer thickness/*/
            hsa = fabs((sx-hx[k])*detZ[0] + (sy-hy[k])*detZ[1] + (sz-cz[k])*detZ[2]);/* distance bwtween source and ph on detector plane/*/

            ntmp3 = DN_DETY1-1;
            ntmp4 = 0;
            ntmp5 = DN_DETX1-1;
            ntmp6 = 0;
            buf[3] = (r2*r2)/(4*dsq)*cosTheta; /* prob of falling through a pinhole/*/
            dp2 = 1000;
            for(dpi=0; dpi<DM_DET; dpi++)
            {
            /**  Defination of buf and index for subroutine of penetration factor (start)  **/
                hsp = fabs(detZ[0]*(sx-dcx)+detZ[1]*(sy-dcy)+detZ[2]*(sz-dcz))+ddz*(dpi-(DM_DET-1)*.5); /* distance between source and detector plate/*/
                f = hsp/hsa;

                buf[0] = sx+(hx[k]-sx)*f;		/* x of the center of the projection/*/
                buf[1] = sy+(hy[k]-sy)*f;		/* y of the center of the projection/*/
                buf[4] = sz+(cz[k]-sz)*f;		/* z value of detector plane/*/
                buf[2] = f*(Rlimit);			/* radius of the projection using Rlimit for penetration factor/*/

                if (IDETPOS&&0)
                {
                    dtmpx=(int)(((buf[0]-dcx)*detX[0]+(buf[1]-dcy)*detX[1]+(buf[4]-dcz)*detX[2])/ddx+DN_DETX1*0.5);
                    dtmpy=(int)(((buf[0]-dcx)*detY[0]+(buf[1]-dcy)*detY[1]+(buf[4]-dcz)*detY[2])/ddy+DN_DETY1*0.5);
                    /*dp2=1000;/*/

                    for(ndy=fmax(0,dtmpy-2); ndy<=fmin(dtmpy+2,DN_DETY1); ndy++)
                    {
                        for(ndx=fmax(0,dtmpx-2); ndx<=fmin(dtmpx+2, DN_DETX1); ndx++)
                        {
                            dn=dpi*NDET_DET1+DN_DETX1*ndy+ndx;
                            dpx=(dp[dn*3]-buf[0])*detX[0]+(dp[dn*3+1]-buf[1])*detX[1]+(dp[dn*3+2]-buf[4])*detX[2];
                            dpy=(dp[dn*3]-buf[0])*detY[0]+(dp[dn*3+1]-buf[1])*detY[1]+(dp[dn*3+2]-buf[4])*detY[2];
                            dp1=dpx*dpx+dpy*dpy;

                            if (dp1<dp2) {dncx=ndx;dncy=ndy;dp2=dp1;}
                        }
                    }

                }
                else
                {	/*??? dcz layer/*/
                    dncx = (int)(((buf[0]-dcx)*detX[0] + (buf[1]-dcy)*detX[1] + (buf[4]-dcz)*detX[2])/ddx + DN_DETX1*0.5);
                    dncy = (int)(((buf[0]-dcx)*detY[0] + (buf[1]-dcy)*detY[1] + (buf[4]-dcz)*detY[2])/ddy + DN_DETY1*0.5);
                }
                /*printf("%d %d\n", dncx, dncy);/*/
                deltan = (int)(fabs(Rlimit/cosTheta*cosThetaD/DDX_DET1)+3);
                ntmp3 = fmin(dncy-deltan, ntmp3);
                ntmp4 = fmax(dncy+deltan, ntmp4);
                ntmp5 = fmin(dncx-deltan,ntmp5);
                ntmp6 = fmax(dncx+deltan, ntmp6);
            }

            ntmp3 = fmax(ntmp3, 0);
            ntmp3 = fmin(ntmp3,DN_DETY1-1);

            ntmp4 = fmin(ntmp4, DN_DETY1-1);
            ntmp4 = fmax(ntmp4,0);

            ntmp5 = fmax(ntmp5,0);
            ntmp5 = fmin(ntmp5, DN_DETY1-1);

            ntmp6 = fmin(ntmp6, DN_DETX1-1);
            ntmp6 = fmax(ntmp6,0);

            for(ndy=ntmp3; ndy<=ntmp4; ndy++)
            {
                for(ndx=ntmp5; ndx<=ntmp6; ndx++)
                {

                    dn = DN_DETX1*ndy+ndx;
                    //printf("DN_DETX1 %d ny %d nx %d dn %d det_pix_mask[dn] %d nsubdet %d\n",DN_DETX1,ndy,ndx,dn,det_pix_mask[dn],nsubdet);getchar();

                    if ((det_pix_mask[dn]==nsubdet)||(det_pix_mask[dn]==100+nsubdet))
                    {

                        nt = 0;

                        for(dpi=0; dpi<DM_DET; dpi++)
                        {
                            ntmp7 = (dn+dpi*NDET_DET1)*3; // pixel index after consider the detector thickness;
                            ntmp8 = 3*dpi;
                            dpos[0+ntmp8] = (dp[ntmp7]-hx[k])*detAXx[k]+(dp[ntmp7+1]-hy[k])*detAXy[k]+(dp[ntmp7+2]-cz[k])*detAXz[k]; /*x position of source pixel in pinhole geometry /*/
                            dpos[1+ntmp8] = (dp[ntmp7]-hx[k])*detAYx[k]+(dp[ntmp7+1]-hy[k])*detAYy[k]+(dp[ntmp7+2]-cz[k])*detAYz[k];
                            dpos[2+ntmp8] = (dp[ntmp7]-hx[k])*detAZx[k]+(dp[ntmp7+1]-hy[k])*detAZy[k]+(dp[ntmp7+2]-cz[k])*detAZz[k];

                            //printf("dn %d dpi %d dpos %f %f %f\n",dn,dpi,dpos[0+ntmp8],dpos[1+ntmp8],dpos[2+ntmp8]);getchar();

                            nsrx = (sx-dp[ntmp7]);
                            nsry = (sy-dp[ntmp7+1]);
                            nsrz = (sz-dp[ntmp7+2]);
                            cos_D[dpi] = fabs(nsrx*detZ[0]+nsry*detZ[1] + nsrz*detZ[2]) / sqrt(nsrx*nsrx + nsry*nsry + nsrz*nsrz); /* for each pixel/*/
                            nm=0;
                            for (j=0;j<fn;j++)
                            {
                                for (i=0;i<fn;i++)
                                {
                                    tmpdx = dp[0 + ntmp7]+(i-(fn-1)/2.0)*ddx/fn*detX[0]+(j-(fn-1)/2.0)*ddy/fn*detY[0];
                                    tmpdy = dp[1 + ntmp7]+(i-(fn-1)/2.0)*ddx/fn*detX[1]+(j-(fn-1)/2.0)*ddy/fn*detY[1];
                                    tmpdz = dp[2 + ntmp7]+(i-(fn-1)/2.0)*ddx/fn*detX[2]+(j-(fn-1)/2.0)*ddy/fn*detY[2];

                                    subdpos[nt*3+0] = (tmpdx-hx[k])*detAXx[k] + (tmpdy-hy[k])*detAXy[k] + (tmpdz-cz[k])*detAXz[k]; /*x position of source pixel in pinhole geometry /*/
                                    subdpos[nt*3+1] = (tmpdx-hx[k])*detAYx[k] + (tmpdy-hy[k])*detAYy[k] + (tmpdz-cz[k])*detAYz[k];
                                    subdpos[nt*3+2] = (tmpdx-hx[k])*detAZx[k] + (tmpdy-hy[k])*detAZy[k] + (tmpdz-cz[k])*detAZz[k];
                                    nt++;nm++;
                                }
                            }
                        }

                        //pr = 1.e-7;
                        pr=pixelprob_det2ML_sum_ph_insert2_double11(buf, dpos, subdpos, spos, ph[phIndex].map, N_step, cos_D, DM_DET, MU_DET, fn);
                        index111=index111+1;
                        if(pr>thresh)
                        {

                             if(det_pix_mask[dn]==nsubdet)
                             {
                                proj[dn] = pr;
                                sen[sn] += pr;
                                pixIndex[npix] = dn;
                                //Proj[2*npix]=dn+1;/*dn start from zero, the actually index of detector of dn should plus 1/*/
                                //Proj[2*npix+1]=pr;
                                npix++;
                             }
                        }
                    }/* det_mask/*/
                }/*ndy/*/
            }


        }/* k nph/*/

        if(0)
        {
            sprintf(fileName,"%s//proj_sn%d_inph%d_iDet%d",DATA_DIR,sn,k,iDet);
            printf("fileName,%s",fileName);
            fp1 = fopen(fileName,"wb");
            fwrite(proj,sizeof(double),DN_DETX1*DN_DETY1,fp1);
            fclose(fp1);
            printf("HERE!save psf to %s \n",fileName);
            getchar();

        }
        //printf("sens %f sn %ld sx %f sy %f sz %f\n",sen[sn],sn,sx,sy,sz);



        // save the response function;

        if(sn==0)
        {
            for(n =0; n<=srf[iThread].nMax; n++)
            {
                srf[iThread].index[n] = 0;
                srf[iThread].value[n] = 0;
            }

            NSIZE = fmax(NIMG,NDET_DET1);
            srf[iThread].index[1] = NSIZE+2; // index vector start from 1; not 0;
            sprk = NSIZE+1;
            //printf("NSIZE %ld NIMG %ld NDET_DET1 %ld in SRF %ld ",NSIZE,NIMG,NDET_DET1,srf[iThread].index[1]);getchar();

        }

        n=0;
        spri=sn+1;

        if(spri<=NDET_DET1)  // diagnal elements;
        {
            srf[iThread].value[spri]=(float) proj[spri-1];
            proj[spri-1]=0.;
        }
        else srf[iThread].value[spri]=0.;

        //off-diagnal elements
        for(l=0; l<npix; l++)
        {
            sprj=pixIndex[l]+1;
            pixIndex[l] = 0;

            if((fabs(proj[sprj-1])>=thresh) && (sprj!=spri))
            {

                if(++sprk>srf[iThread].nMax) { printf("error, the maxtri is too small to store the response function!!! \n"); getchar(); }
                srf[iThread].value[sprk]=(float)proj[sprj-1];
                srf[iThread].index[sprk]=sprj;
                proj[sprj-1]=0.;

            }
        }
        srf[iThread].index[spri+1]=sprk+1;


    } // end of sn


    for (spri=(NIMG+1); spri<=NSIZE; spri++) // the size of maxtrix is NSIZExNSIZE;NSIZE is maxium of (NIMG,NDET_DET1);
    {
        srf[iThread].index[spri+1]=srf[iThread].index[NIMG+1];
    }

    //printf("NSIZE %ld mSIZE %ld \n",srf[iThread].index[1]-1,srf[iThread].index[srf[iThread].index[1]-1]-1);
    //printf("before TranSRF.c\n");getchar();getchar();getchar();




   tranSRF(sen,iThread);

    //printf("gensubSRF.c readyTran %d \n",srf[iThread+1].readyTran);
    //printf("gensubSRF.c readySave %d \n",srf[iThread+1].readySave);

    //srf[iThread+1].readyTran = 1;
   // srf[iThread+1].readySave = 0;






    free(dpos);
    //free(dsubp);
    free(subdpos);
    free(hx);
    free(hy);
    free(cz);
    free(cos_D);
    free(detAXx);
    free(detAXy);
    free(detAXz);
    free(detAYx);
    free(detAYy);
    free(detAYz);
    free(detAZx);
    free(detAZy);
    free(detAZz);
    free(DDX_ph);
    free(DDY_ph);
    free(proj);
    free(pixIndex);
    free(sen);

    if(iThread ==0){
        time(&finish);
        duration=difftime (finish,start);
        printf("the time needs to calculate one srf is %f\n",duration);
    }

    return(1);
}
