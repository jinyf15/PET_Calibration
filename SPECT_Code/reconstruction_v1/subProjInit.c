#include"headFile.h"

struct projection *subProjInit(int nDet,int NX,int NY,int nSubDet,int nSou)
{
    int iDet;
    int iSou;
    int projIndex;
    struct projection *proj;

    proj = (struct projection *)malloc(sizeof(struct projection)*nDet*nSou);

    projIndex = 0;
    for(iSou = 0; iSou < nSou; iSou++)
    {

        for(iDet=0; iDet<nDet; iDet++)
        {
            proj[projIndex].nDet = nDet;
            proj[projIndex].iDet = iDet;
            proj[projIndex].nPosi = nSou;
            proj[projIndex].iPosi = iSou;
            proj[projIndex].nSubDet = nSubDet;

            proj[projIndex].NX = NX;
            proj[projIndex].NY = NY;

            loadProjection(proj,projIndex);
            loadDetPixMaskForRecon(proj,projIndex);


            projIndex++;
        }

    }


    return(proj);

}
