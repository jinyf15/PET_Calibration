#include"headFile.h"

void subPhInit(struct pinMap *ph,int iph)
{
    char fileName[1000];
    FILE *fp1;
    int nSize;
    double buffer;
    int ii;


    sprintf(fileName,"%s//phMap//map%d",DATA_DIR,iph);
    fp1 = fopen(fileName,"rb");
    checkFile(fp1,fileName);


    fread(&buffer,sizeof(double),1,fp1); // number of pinhole in detector
    ph[iph].nph = (int) buffer;

    fread(&buffer,sizeof(double),1,fp1); // which ph it is.
    ph[iph].iph = (int) buffer;

    fread(&buffer,sizeof(double),1,fp1);
    ph[iph].radius = (double) buffer;

    fread(&buffer,sizeof(double),1,fp1);
    ph[iph].Rlimit = (double) buffer;

    fread(&buffer,sizeof(double),1,fp1);
    ph[iph].openAngle = (double) buffer;

    fread(&buffer,sizeof(double),1,fp1);
    ph[iph].angleLimit = (double) buffer;

    fread(&buffer,sizeof(double),1,fp1);
    ph[iph].incidenceAngle = (double) buffer;

    fread(&buffer,sizeof(double),1,fp1);
    ph[iph].length = (double) buffer;

    fread(&buffer,sizeof(double),1,fp1);
    ph[iph].channel = (double) buffer;

    fread(&buffer,sizeof(double),1,fp1);
    ph[iph].leadThickness = (double) buffer;


    nSize = 1;
    for(ii=0;ii<4;ii++)  // read the number of step in each demension of map;
    {
        fread(&buffer,sizeof(double),1,fp1); // which ph it is.
        ph[iph].NStep[ii] = (int) buffer;
        nSize = nSize*((int) buffer);

    }


    ph[iph].map = (double *)malloc(sizeof(double)*(nSize+4));

    for(ii=0;ii<nSize+4;ii++)  // read the number of step in each demension of map;
    {
        fread(&buffer,sizeof(double),1,fp1); // which ph it is.
        ph[iph].map[ii] = (double) buffer;
        if(ii==0) printf(" lading map:first element of map%d %f \n",iph,ph[iph].map[ii]);
    }

    fclose(fp1);

    // attached the step size to end of map
    ph[iph].dirt[0] = ph[iph].map[nSize+0];
    ph[iph].dirt[1] = ph[iph].map[nSize+1];
    ph[iph].dirt[2] = ph[iph].map[nSize+2];
    ph[iph].dirt[3] = ph[iph].map[nSize+3];



    if(ph[iph].iph == 0)
    {
        printf("nph %d inph %d r %lf Rlimit %lf open %lf angllimt %lf\n",ph[iph].nph,ph[iph].iph,ph[iph].radius,ph[iph].Rlimit,ph[iph].openAngle/PI*180,ph[iph].angleLimit/PI*180);
        printf("N1 %d N2 %d N3 %d N4 %d dirt %lf %lf %lf %lf\n",ph[iph].NStep[0],ph[iph].NStep[1],ph[iph].NStep[2],ph[iph].NStep[3],ph[iph].dirt[0],ph[iph].dirt[1],ph[iph].dirt[2],ph[iph].dirt[3]);
        printf("incidentangle %f length %lf channel %lf leadThickness %lf\n\n",ph[iph].incidenceAngle/PI*180.0,ph[iph].length,ph[iph].channel,ph[iph].leadThickness);
        //printf("press anykey to continue\n");
        //printf("here");
        getchar();

    }


}
