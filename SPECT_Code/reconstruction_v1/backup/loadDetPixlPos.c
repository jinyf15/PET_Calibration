#include"headFile.h"

void loadDetPixPos(struct detector *det,int nDet,int iDet)
{
    int NX;
    int NY;
    int NZ;
    int iPix;
    int iZ;
    double bufPx;
    double bufPy;
    double dz;
    char fileName[1000];
    FILE *fp1;

    NX = det[iDet].NX;
    NY = det[iDet].NY;
    NZ = det[iDet].NZ;
    dz = det[iDet].dz;


    det[iDet].dp = (double *)malloc(sizeof(double)*NZ*NX*NY*3);

    sprintf(fileName,"%s//det_position%d.txt",DATA_DIR,iDet);
    fp1=fopen(fileName,"r");
    checkFile(fp1,fileName);

    //printf("NX,NY, %d %d\n",NX,NY);getchar();

    for(iPix=0;iPix<NX*NY;iPix++)
    {
        fscanf(fp1,"%le\n",&bufPx);
        fscanf(fp1,"%le\n",&bufPy);

        for(iZ = 0; iZ < NZ; iZ++)
        {
            //printf("z direction dz %f NZ %d %f\n",dz,NZ,iZ*dz-(NZ/2.0-0.5)*dz);
            det[iDet].dp[iZ*NX*NY*3+3*iPix+0] = det[iDet].detCent[0]+bufPx*det[iDet].detX[0]+bufPy*det[iDet].detY[0] + (iZ*dz-(NZ/2.0-0.5)*dz)*det[iDet].detZ[0];
            det[iDet].dp[iZ*NX*NY*3+3*iPix+1] = det[iDet].detCent[1]+bufPx*det[iDet].detX[1]+bufPy*det[iDet].detY[1] + (iZ*dz-(NZ/2.0-0.5)*dz)*det[iDet].detZ[1];
            det[iDet].dp[iZ*NX*NY*3+3*iPix+2] = det[iDet].detCent[2]+bufPx*det[iDet].detX[2]+bufPy*det[iDet].detY[2] + (iZ*dz-(NZ/2.0-0.5)*dz)*det[iDet].detZ[2];
        }

    //    if(iPix==0||iPix==NX*NY-1)
    //    {
    //        printf("iPix %d %le %le %le\n",iPix,det[iDet].dp[3*iPix+0],det[iDet].dp[3*iPix+1],det[iDet].dp[3*iPix+2]);
    //        getchar();
    //    }


    }

}
