#include"headFile.h"

void calDetGlobalPos(struct detector *det,double *sysPar,int iDet)
{
    double sinAlpha,cosAlpha;
    double sinPhi,cosPhi;
    double sinBeta,cosBeta;
    double detCent[3];

    // euler angle of detector plane;
    sinAlpha=sin(sysPar[4]);
    cosAlpha=cos(sysPar[4]);
    sinPhi=sin(sysPar[5]+sysPar[10]);
    cosPhi=cos(sysPar[5]+sysPar[10]);
    sinBeta=sin(sysPar[3]);
    cosBeta=cos(sysPar[3]);

    // global axis of det
    det[iDet].detX[0]=cosAlpha*cosBeta;
    det[iDet].detX[1]=cosAlpha*sinBeta;
    det[iDet].detX[2]=-sinAlpha;
    det[iDet].detY[0]=-sinBeta*cosPhi+sinPhi*sinAlpha*cosBeta;
    det[iDet].detY[1]=cosPhi*cosBeta+sinPhi*sinAlpha*sinBeta;
    det[iDet].detY[2]=sinPhi*cosAlpha;
    det[iDet].detZ[0]=sinBeta*sinPhi+cosBeta*sinAlpha*cosPhi;
    det[iDet].detZ[1]=sinBeta*sinAlpha*cosPhi-cosBeta*sinPhi;
    det[iDet].detZ[2]=cosAlpha*cosPhi;

    // global center of detetor;
    detCent[0]=sysPar[0];
    detCent[1]=sysPar[1];
    detCent[2]=sysPar[2];


    det[iDet].detCent[0]=det[iDet].detX[0]*detCent[0]+det[iDet].detY[0]*detCent[1]+det[iDet].detZ[0]*detCent[2];
    det[iDet].detCent[1]=det[iDet].detX[1]*detCent[0]+det[iDet].detY[1]*detCent[1]+det[iDet].detZ[1]*detCent[2];
    det[iDet].detCent[2]=det[iDet].detX[2]*detCent[0]+det[iDet].detY[2]*detCent[1]+det[iDet].detZ[2]*detCent[2];

}
