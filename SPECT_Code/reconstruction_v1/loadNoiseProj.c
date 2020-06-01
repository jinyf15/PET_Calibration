#include"headFile.h"
#include"globalVariable.h"

void loadNoiseProj(struct projection * proj,int index)
{
    int iSouP;
    int nSouP;
    int iDet;
    int nDet;
    int ii;
    int dN;
    char fileName[1000];
    FILE *fp1;
    int iPix;
    int nRead;
    int projIndex;

    dN = proj[0].NX * Proj[0].NY;
    nSouP = proj[0].nPosi;
    nDet = proj[0].nDet;

    for(iSouP = 0;iSouP < nSouP; iSouP++)
    {
        for(iDet = 0; iDet < nDet; iDet++)
        {

            projIndex = iSouP * nDet + iDet;


            sprintf(fileName,"%s//proj_iDet%d_iPosi%d_iSet%d", &(proj[projIndex].projFolder[0]), iDet, iSouP,index);


            printf("loading projection %s\n",fileName);

            fp1 = fopen(fileName,"rb");
            checkFile(fp1,fileName);

            nRead = fread(&(proj[projIndex].detImage[0]),sizeof(double),dN,fp1);
            if(nRead != dN)
            {
                printf("error % reading file from %s,number of data request %d,actualy read %d\n",fileName,dN,nRead);
                getchar();getchar();
                exit(-1);
            }


            fclose(fp1);


        }

    }
}
