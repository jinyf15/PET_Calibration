#include"headFile.h"
#include"globalVariable.h"

void genNoiseProj(struct projection * Proj,int index)
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

    dN = Proj[0].NX * Proj[0].NY;
    nSouP = Proj[0].nPosi;
    nDet = Proj[0].nDet;

    for(iSouP = 0;iSouP < nSouP; iSouP++)
    {
        for(iDet = 0; iDet < nDet; iDet++)
        {

            ii = iSouP * nDet + iDet;
            sprintf(fileName,"%s//mproj_iDet%d_iPosi%d",&(Proj[0].projFolder[0]),iDet,iSouP);  // mean projection
            fp1 = fopen(fileName,"rb");
            nRead = fread( &( Proj[ii].detImage[0] ), sizeof(double), dN,fp1);
            fclose(fp1);

            if(nRead =! dN)
            {
                printf("error in reading file from %s\n",fileName);
                getchar();
                getchar();
                exit(-1);
            }



            for(iPix = 0; iPix < dN; iPix++) // generate poisson number
            {
                Proj[ii].detImage[iPix] = (double) poidev( (float) Proj[ii].detImage[iPix], &(randSeed));
            }

            sprintf(fileName,"%s//proj_iDet%d_iPosi%d_iSet%d",&(Proj[0].projFolder[0]),iDet,iSouP,index);  // projection
            fp1 = fopen(fileName,"wb");
            fwrite( &( Proj[ii].detImage[0] ), sizeof(double), dN,fp1);
            fclose(fp1);


        }

    }
}
