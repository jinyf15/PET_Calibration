#include"headFile.h"

void phGlobalPosi(struct detector *det,int iDet,double *sysPa)
{
    int iph;

    for(iph=0;iph<det[iDet].nphPerDet;iph++)
    {

    det[iDet].phCenter[iph] = (double *)malloc(sizeof(double)*3);

    det[iDet].phX[iph] = (double *)malloc(sizeof(double)*3);
    det[iDet].phY[iph] = (double *)malloc(sizeof(double)*3);
    det[iDet].phZ[iph] = (double *)malloc(sizeof(double)*3);

    det[iDet].phCenter[iph][0] = sysPa[11+iph*10];
    det[iDet].phCenter[iph][1] = sysPa[12+iph*10];
    det[iDet].phCenter[iph][2] = sysPa[13+iph*10];

    det[iDet].phZ[iph][0] = sysPa[14+iph*10];
    det[iDet].phZ[iph][1] = sysPa[15+iph*10];
    det[iDet].phZ[iph][2] = sysPa[16+iph*10];

    det[iDet].phY[iph][0] = sysPa[17+iph*10];
    det[iDet].phY[iph][1] = sysPa[18+iph*10];
    det[iDet].phY[iph][2] = sysPa[19+iph*10];

    // phX=phY X PHZ;
    det[iDet].phX[iph][0] = det[iDet].phY[iph][1]*det[iDet].phZ[iph][2]-det[iDet].phY[iph][2]*det[iDet].phZ[iph][1];
    det[iDet].phX[iph][1] = det[iDet].phY[iph][2]*det[iDet].phZ[iph][0]-det[iDet].phY[iph][0]*det[iDet].phZ[iph][2];
    det[iDet].phX[iph][2] = det[iDet].phY[iph][0]*det[iDet].phZ[iph][1]-det[iDet].phY[iph][1]*det[iDet].phZ[iph][0];

    printf("nph_Det %d iph_Det %d iph %d phCenter %lf %lf %lf \n",det[iDet].nphPerDet,iph,det[iDet].phIndex[iph], det[iDet].phCenter[iph][0],det[iDet].phCenter[iph][1],det[iDet].phCenter[iph][2]);
    printf("nph_Det %d iph_Det %d iph %d phX %lf %lf %lf \n",det[iDet].nphPerDet,iph,det[iDet].phIndex[iph],det[iDet].phX[iph][0],det[iDet].phX[iph][1],det[iDet].phX[iph][2]);
    printf("nph_Det %d iph_Det %d iph %d phY %lf %lf %lf \n",det[iDet].nphPerDet,iph,det[iDet].phIndex[iph],det[iDet].phY[iph][0],det[iDet].phY[iph][1],det[iDet].phY[iph][2]);
    printf("nph_Det %d iph_Det %d iph %d phZ %lf %lf %lf \n",det[iDet].nphPerDet,iph,det[iDet].phIndex[iph],det[iDet].phZ[iph][0],det[iDet].phZ[iph][1],det[iDet].phZ[iph][2]);



    }

    printf("press anykey to continue\n");
   // getchar();


}

