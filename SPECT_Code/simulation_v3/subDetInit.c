#include"headFile.h"

struct detector *subDetInit(int nDet,int NX,int NY,int NZ,int nSubDet,int nDivide,double dx,double dy,double dz,double attenCoef)
{
    int iDet;
    struct detector *det;

    det = (struct detector *)malloc(sizeof(struct detector)*nDet);

    for(iDet=0;iDet<nDet;iDet++)
    {
        det[iDet].iDet = iDet;
        det[iDet].nDet = nDet;
        det[iDet].nDivide = nDivide;
        det[iDet].nSubDet = nSubDet;

        det[iDet].NX = NX;
        det[iDet].NY = NY;
        det[iDet].NZ = NZ;

        det[iDet].dx = dx;
        det[iDet].dy = dy;
        det[iDet].dz = dz;

        det[iDet].attenCoef = attenCoef;


        loadDetFile(det,nDet,iDet); // calculated detector center position and global coordinate;

        loadDetPixPos(det,nDet,iDet); // calulated detector pixel global position;

        loadDetPixMask(det,nDet,iDet);

        //loadDetFile(det,nDet,NX,NY,NZ,nSubDet,dx,dy,dz);

    }

    return(det);

}
