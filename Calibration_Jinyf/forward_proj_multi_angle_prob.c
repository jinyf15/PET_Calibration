#include <math.h>
#include "mex.h"
#include "matrix.h"
#include "stdio.h"
#include <string.h>

double abs_t(double a)
{
    double b;
    if(a >=0)
        b=a;
    else
        b=0-a;
    return b;
}

int max(int a[])
{
    int sss = 0;
    int i;
    for(i = 0;i < 4;i++)
    {

        //mexPrintf("%d ",a[i]);
        if (a[i]>sss)
        {
            sss = a[i];
        }
    }

    return sss;
}


int check_if_inside_quadrange(double *pt1, double *pt2, double *pt3, double *pt4, double *pt,int fn)
{

	double k1, k2, k3, k4, a1, a2, a3, a4, b1, b2, b3, b4, c1, c2, c3, c4;
    int n=0;
    int sn;

//     mexPrintf("pt1:(%f,%f), pt2(%f,%f), pt3(%f,%f), pt4(%f,%f)\n",pt1[0],pt1[1],pt2[0],pt2[1],pt3[0],pt3[1],pt4[0],pt4[1]);
	k1=(pt2[1]-pt1[1])/(pt2[0]-pt1[0]);
	a1=k1*(pt3[0]-pt1[0])+pt1[1]-pt3[1];
	b1=k1*(pt4[0]-pt1[0])+pt1[1]-pt4[1];

	k2=(pt3[1]-pt2[1])/(pt3[0]-pt2[0]);
	a2=k2*(pt1[0]-pt2[0])+pt2[1]-pt1[1];
	b2=k2*(pt4[0]-pt2[0])+pt2[1]-pt4[1];

	k3=(pt4[1]-pt3[1])/(pt4[0]-pt3[0]);
	a3=k3*(pt1[0]-pt3[0])+pt3[1]-pt1[1];
	b3=k3*(pt2[0]-pt3[0])+pt3[1]-pt2[1];

	k4=(pt1[1]-pt4[1])/(pt1[0]-pt4[0]);
	a4=k4*(pt3[0]-pt4[0])+pt4[1]-pt3[1];
	b4=k4*(pt2[0]-pt4[0])+pt4[1]-pt2[1];



	for(sn=0;sn<fn*fn;sn++)
    {
//         mexPrintf("pt index:%d, (%f,%f)\n",sn+1,pt[sn*2],pt[sn*2+1]);
		c1=k1*(pt[sn*2]-pt1[0])+pt1[1]-pt[sn*2+1];
		if (a1*c1<0 || b1*c1<0 ) continue;

		c2=k2*(pt[sn*2]-pt2[0])+pt2[1]-pt[sn*2+1];
		if (a2*c2<0 || b2*c2<0 ) continue;


		c3=k3*(pt[sn*2]-pt3[0])+pt3[1]-pt[sn*2+1];
		if (a3*c3<0 || b3*c3<0 ) continue;

		c4=k4*(pt[sn*2]-pt4[0])+pt4[1]-pt[sn*2+1];

		if (a4*c4<0 || b4*c4<0) continue;

		n=n+1;
//         mexPrintf("n:%d\n",n);
	}

	return n;
}

double forward_proj_pet_parallel(double dx1, double dy1, double dz1, double dx2,double dy2, double dz2,
                                 double *dpc,int dn1,double sx, double sy, double sz,  double *det1Z, double *det2X, double *det2Y, double *det2Z,
                                 double *pt, double *reqpara,int fn,int layer1,int layer2)
{
    int DM_DET1, DM_DET2 ,NDET_DET1,NDET_DET2,cnt;
    double DDZ_DET1,DDZ_DET2,mu_det1,mu_det2,s1,s2;
    double pr;
    int m,nn;
    double difx,dify,difz;
    double dsq1,hs1,cosTheta1,cosTheta2;

    double p1,ltmp,w1,w2;
    double sxt,syt,szt;
    double x_tmp,y_tmp,z_tmp;
    double pt1[2],pt2[2],pt3[2],pt4[2];
    double A2,p2;
    double dxct[4],dyct[4],dzct[4];



    DM_DET1=(int)reqpara[0];
    DM_DET2=(int)reqpara[1];
    DDZ_DET1=reqpara[2];
    DDZ_DET2=reqpara[3];
    NDET_DET1=(int)reqpara[4];
    NDET_DET2=(int)reqpara[5];
    mu_det1=reqpara[6];
    mu_det2=reqpara[7];
    s1=reqpara[8];
    s2=reqpara[9];

    pr=0;

      difx=sx-dx1; dify=sy-dy1; difz=sz-dz1;
      dsq1=difx*difx+dify*dify+difz*difz;
      hs1=abs_t(det1Z[0]*difx+det1Z[1]*dify+det1Z[2]*difz);
      cosTheta1=hs1/sqrt(dsq1);
      cosTheta2=abs_t(det2Z[0]*difx+det2Z[1]*dify+det2Z[2]*difz)/sqrt(dsq1);

      p1=s1*0.07957747154594766788444188168626/dsq1*cosTheta1*1000;
      ltmp=mu_det1*DDZ_DET1/cosTheta1;
      w1=(1-exp(0-ltmp))*exp((0-layer1)*ltmp);
      p1=p1*w1;

       //mexPrintf("dx1:%f, dy1: %f, dz1:%f, difx: %f, dify: %f, difz: %f\n", dx1,dy1,dz1,difx,dify,difz);
 //      mexPrintf("det1Z:(%f %f %f), det2Z:(%f %f %f)\n", det1Z[0],det1Z[1],det1Z[2],det2Z[0],det2Z[1],det2Z[2]);
       //mexPrintf("dsq1: %f, hs1: %f, cosTheta1: %f, p1: %f, ltmp: %f, w1: %f, p1:%f\n",dsq1,hs1,cosTheta1,p1,ltmp,w1,p1);
      ltmp=mu_det2*DDZ_DET2/cosTheta2;

          sxt=sx*det2X[0]+sy*det2X[1]+sz*det2X[2]-dx2*det2X[0]-dy2*det2X[1]-dz2*det2X[2];
          syt=sx*det2Y[0]+sy*det2Y[1]+sz*det2Y[2]-dx2*det2Y[0]-dy2*det2Y[1]-dz2*det2Y[2];
          szt=sx*det2Z[0]+sy*det2Z[1]+sz*det2Z[2]-dx2*det2Z[0]-dy2*det2Z[1]-dz2*det2Z[2];

          w2=(1-exp(0-ltmp))*exp((0-layer2)*ltmp);
//           mexPrintf("i: %d, j: %d\n", i,j);

          for(m=0;m<4;m++)
          {
              x_tmp=dpc[(dn1-1)*12+m*3];
              y_tmp=dpc[(dn1-1)*12+m*3+1];
              z_tmp=dpc[(dn1-1)*12+m*3+2];

              dxct[m]=(x_tmp-dx2)*det2X[0]+(y_tmp-dy2)*det2X[1]+(z_tmp-dz2)*det2X[2];
              dyct[m]=(x_tmp-dx2)*det2Y[0]+(y_tmp-dy2)*det2Y[1]+(z_tmp-dz2)*det2Y[2];
              dzct[m]=(x_tmp-dx2)*det2Z[0]+(y_tmp-dy2)*det2Z[1]+(z_tmp-dz2)*det2Z[2];
              //mexPrintf("x, y, z:%f %f %f\ndx, dy, dz:%f %f %f\n",x_tmp,y_tmp,z_tmp,dxct[m],dyct[m],dzct[m]);
          }

          pt1[0]=(sxt-dxct[0])*(0-dzct[0])/(szt-dzct[0])+dxct[0];

          pt2[0]=(sxt-dxct[1])*(0-dzct[1])/(szt-dzct[1])+dxct[1];

          pt3[0]=(sxt-dxct[2])*(0-dzct[2])/(szt-dzct[2])+dxct[2];
          pt4[0]=(sxt-dxct[3])*(0-dzct[3])/(szt-dzct[3])+dxct[3];


          pt1[1]=(syt-dyct[0])*(0-dzct[0])/(szt-dzct[0])+dyct[0];

          pt2[1]=(syt-dyct[1])*(0-dzct[1])/(szt-dzct[1])+dyct[1];
          pt3[1]=(syt-dyct[2])*(0-dzct[2])/(szt-dzct[2])+dyct[2];
          pt4[1]=(syt-dyct[3])*(0-dzct[3])/(szt-dzct[3])+dyct[3];

          A2=abs_t((pt3[0]-pt1[0])*(pt4[1]-pt2[1])*1000-(pt4[0]-pt2[0])*(pt3[1]-pt1[1])*1000)*0.5;
          //mexPrintf("pt1: %f, %f\n pt2: %f, %f\n pt3: %f, %f\n pt4: %f, %f\n",pt1[0],pt1[1],pt2[0],pt2[1],pt3[0],pt3[1],pt4[0],pt4[1]);
          //mexPrintf("abs: %f",abs_t(pt3[0]-pt1[0]));
          nn=check_if_inside_quadrange(pt1,pt2,pt3,pt4,pt,fn);

          p2=nn*s2*100/fn/fn/A2;

          if(p2>1) p2=1.0;

          p2=p2*w2;

          pr=p1*p2;
          //mexPrintf("w2: %f, nn: %d, s2: %f, A2: %f, fn %d, p2: %f, pr: %f\n",w2,nn,s2,A2,fn,p2,pr);



    return pr;

}

void forward_proj_multi_angle_prob(double *mdata,int nlen_data, int nlen_sr,
        double *dp, double *dpc0,double *spx,
        double *spy,double *spz, double *Mtr, double *pt, double *reqpara,int fn,double *prob,double *spc)//FILE *fp,FILE *fp1)
{
    int i,j;
    int dn1,dn2;
    double dx1,dx2,dy1,dy2,dz1,dz2;
    double det1Z[3],det2X[3],det2Y[3],det2Z[3];
    int PET1, PET2, flag, layer1, layer2;
    double sx,sy,sz;
    double sprob;
    double pr,logpr;
    double inputpara[10];
    //double dpc[16*16*4*3*4];
    int NDET_DET[4], DM_DET[4];
    int len_data,len_sr;
    double *dpc;
    len_data=nlen_data;
    len_sr=nlen_sr;

    pr=0;
     logpr=0; //we wanna change the log operation outside the C code.
//       logpr=1;

//     FILE *fp,*fp1;
//     fp=fopen("prob2_perevent.txt","a");
//     fp1=fopen("prob2_perpt.txt","a");
//
//     fprintf(fp,"prob. per event\n");
//     fprintf(fp1,"prob. per sr.pixle\n");
    for(i = 0;i < 4;i++)
    {
        NDET_DET[i] = reqpara[2+5*i];
        DM_DET[i] = reqpara[5*i];
    }
   int maxsize = max(NDET_DET)*max(DM_DET);
   dpc = (double *)malloc(maxsize*12*sizeof(double));
   for(i=0;i<len_data;i++)
    {
        pr=0;
        flag = 0;
        for(j = 0;j < 4;j++)
            if ((int)mdata[4*i+j])
                flag = flag + (1<<j);
        PET1 = (flag-(flag&(flag-1)))/2;
        PET2 = ((flag - (1<<PET1))>>2) + 1;
        //mexPrintf("PET1=%d PET2=%d\n",PET1,PET2);
        dn1=(int)mdata[4*i+PET1];
        dn2=(int)mdata[4*i+PET2];
        layer1 = DM_DET[PET1]-(dn1%(NDET_DET[PET1]/4))/(NDET_DET[PET1]/4)-1;
        layer2 = DM_DET[PET2]-(dn2%(NDET_DET[PET1]/4))/(NDET_DET[PET2]/4)-1;
        dx1=dp[(dn1-1)*3+PET1*maxsize*3];
        dy1=dp[(dn1-1)*3+1+PET1*maxsize*3];
        dz1=dp[(dn1-1)*3+2+PET1*maxsize*3];
        dx2=dp[(dn2-1)*3+PET2*maxsize*3];
        dy2=dp[(dn2-1)*3+1+PET2*maxsize*3];
        dz2=dp[(dn2-1)*3+2+PET2*maxsize*3];
        //mexPrintf("5\n");
        for(j = 0;j < 5;j++)
        {
            inputpara[j*2] = reqpara[PET1*5+j];
            inputpara[j*2+1] = reqpara[PET2*5+j];
        }
        for(j = 0;j < 3;j++)
        {
            det1Z[j] = Mtr[PET1*9+6+j];
            det2X[j] = Mtr[PET2*9+j];
            det2Y[j] = Mtr[PET2*9+3+j];
            det2Z[j] = Mtr[PET2*9+6+j];
        }

        //mexPrintf("PET1=%d PET2=%d\n",PET1,PET2);
        //mexPrintf("%f %f %f %f %f %f\n",dx1,dy1,dz1,dx2,dy2,dz2);

        memcpy(dpc,dpc0+PET1*maxsize*12,maxsize*12*sizeof(double));
        //for (j = 0;j < maxsize*12;j++)
          //  mexPrintf("%d\n",dpc[j]-dpc0[j+PET1*maxsize*12]);
        for(j=0;j<len_sr;j++)
        {
            sx=spx[j];sy=spy[j];sz=spz[j];
            sprob=spc[j];
            //mexPrintf("i=%d j=%d",i,j);
            pr=pr+sprob*forward_proj_pet_parallel(dx1,dy1,dz1,dx2,dy2,dz2,dpc,dn1,sx,sy,sz,
                                                  det1Z,det2X,det2Y,det2Z,pt,inputpara,fn,layer1,layer2);
//             fprintf(fp1,"%f\n",pr);
         }

         //mexPrintf("i=%d prob=%f\n",i,pr);
         if(pr>0)
             logpr=logpr+log(pr);
	     //logpr = logpr*pr;
         else
             logpr=logpr-100;
	     //logpr=logpr-60;
             //logpr = logpr;
//         fprintf(fp,"%f\n",logpr);

    }
    prob[0]=logpr;
    free(dpc);dpc = NULL;
//     fprintf(fp,"%f\n",logpr);
//     fclose(fp);
//     fclose(fp1);
    return;
}

void mexFunction( int nlhs, mxArray *plhs[],
		  int nrhs, const mxArray*prhs[] )
{
    double *mdata;
    int len_data,len_sr,fn;

    double *dp,*dpc,*spc;
    double *spx,*spy,*spz;
    double *Mtr,*pt,*reqpara;
    double *prob;

    int i,j;

    //FILE *fp,*fp1;


    mdata=mxGetPr(prhs[0]);
    len_data=(int)(*mxGetPr(prhs[1]));
    len_sr=(int)(*mxGetPr(prhs[2]));




    dp=mxGetPr(prhs[3]);

    dpc=mxGetPr(prhs[4]);

    spx=mxGetPr(prhs[5]);
    spy=mxGetPr(prhs[6]);
    spz=mxGetPr(prhs[7]);

    Mtr = mxGetPr(prhs[8]);

    pt=mxGetPr(prhs[9]);

    reqpara=mxGetPr(prhs[10]);
    fn = (int)(*(mxGetPr(prhs[11])));
    spc=mxGetPr(prhs[12]);

    plhs[0]=mxCreateDoubleMatrix( 1, 1,mxREAL);

    prob=mxGetPr(plhs[0]);

    forward_proj_multi_angle_prob(mdata,len_data,len_sr,dp,dpc,spx,spy,spz,Mtr,pt,reqpara,fn,prob,spc);//fp,fp1);

    return;
}





