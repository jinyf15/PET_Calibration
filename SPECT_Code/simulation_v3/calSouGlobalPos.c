#include"headFile.h"

void calSouGlobalPos(struct source *sou,int nPos,int iPos,int NX, int NY, int NZ, double dx, double dy,double dz)
{
    double sinAlpha,cosAlpha;
    double sinPhi,cosPhi;
    double sinBeta,cosBeta;
    unsigned long int nSou;
    unsigned long int iSou;
    int iX;
    int iY;
    int iZ;
    double sx;
    double sy;
    double sz;


    double souX[3];
    double souY[3];
    double souZ[3];

    sou[iPos].nPosi = nPos;
    sou[iPos].iPosi = iPos;

    sou[iPos].dx = dx;
    sou[iPos].dy = dy;
    sou[iPos].dz = dz;

    sou[iPos].NX = NX;
    sou[iPos].NY = NY;
    sou[iPos].NZ = NZ;



    sou[iPos].sp = (double *)malloc(sizeof(double)*NX*NY*NZ*3);

    sinAlpha = sin(sou[iPos].eulerAg[1]*PI/180.0); // rotation sequence z(beta)->y(alpha)->x(phi);
    cosAlpha = cos(sou[iPos].eulerAg[1]*PI/180.0);
    sinPhi   = sin(sou[iPos].eulerAg[2]*PI/180.0);
    cosPhi   = cos(sou[iPos].eulerAg[2]*PI/180.0);
    sinBeta  = sin(sou[iPos].eulerAg[0]*PI/180.0);
    cosBeta  = cos(sou[iPos].eulerAg[0]*PI/180.0);
    //printf("%f %f \n%f %f \n%f %f \n",sinAlpha,cosAlpha,sinPhi,cosPhi,sinBeta,cosBeta);
    //getchar();
    // global axis of det
    souX[0] = cosAlpha*cosBeta;
    souX[1] = cosAlpha*sinBeta;
    souX[2] = -sinAlpha;

    souY[0] = -sinBeta*cosPhi+sinPhi*sinAlpha*cosBeta;
    souY[1] = cosPhi*cosBeta+sinPhi*sinAlpha*sinBeta;
    souY[2] = sinPhi*cosAlpha;

    souZ[0] = sinBeta*sinPhi+cosBeta*sinAlpha*cosPhi;
    souZ[1] = sinBeta*sinAlpha*cosPhi-cosBeta*sinPhi;
    souZ[2] = cosAlpha*cosPhi;


    iSou = 0;
    for(iZ = 0;iZ < NZ; iZ ++)
    {
        for(iY = 0;iY < NY; iY ++ )
        {
            for(iX = 0; iX < NX; iX ++)
            {
                sx = dx*(iX-(NX-1)/2.0);
                sy = dy*(iY-(NY-1)/2.0);
                sz = dz*(iZ-(NZ-1)/2.0);

                sou[iPos].sp[3*iSou+0] = sx*souX[0]+sy*souY[0]+sz*souZ[0]+sou[iPos].souCenter[0]; // the position of source pixel in global geometry;
                sou[iPos].sp[3*iSou+1] = sx*souX[1]+sy*souY[1]+sz*souZ[1]+sou[iPos].souCenter[1];
                sou[iPos].sp[3*iSou+2] = sx*souX[2]+sy*souY[2]+sz*souZ[2]+sou[iPos].souCenter[2];
                iSou++;
                //printf("sou %f %f %f\n",sou[iPos].sp[3*iSou+0],sou[iPos].sp[3*iSou+1],sou[iPos].sp[3*iSou+2]);getchar();
            }
        }
    }



}
