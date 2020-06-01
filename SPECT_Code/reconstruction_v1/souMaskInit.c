#include"headFile.h"

void souMaskInit(struct source *sou,int nPos, int iPos)
{
    unsigned long int iSou;
    int NX;
    int NY;
    int NZ;

    NX = sou[iPos].NX;
    NY = sou[iPos].NY;
    NZ = sou[iPos].NZ;

    sou[iPos].souMask = (int *)malloc(sizeof(int)*NX*NY*NZ);

    for(iSou = 0; iSou < NX*NY*NZ; iSou ++)
    {
        sou[iPos].souMask[iSou] = 1;

    }

}
