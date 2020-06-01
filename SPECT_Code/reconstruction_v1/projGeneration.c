#include "headFile.h"
#include "globalVariable.h"

void projGeneration(void * threadID)
{
    int nIter;
    int iIter;
    int iLoad;
    int nLoad;
    int iThread;
    int iDet;
    int nDet;
    int nSouP;
    int iSouP;
    int nSubDet;
    int iSubDet;
    int iSRF;
    int iProj;
    unsigned int nRow;
    double *imgTMP;
    double *projTMP;
    double *projRatio;
    int dNx;
    int dNy;
    int dN;
    int iPix;
    unsigned int nIMG;
    unsigned int sn;
    int nSubSet;
    int ii,jj,kk;




    char fileName[1000];
    FILE *fp1;

    iThread = (int) threadID;

    printf("thread ID %d\n",iThread);

    nIter = PSeq[0][0].nIter;
    nDet = PSeq[0][0].nDet;
    nSouP = PSeq[0][0].nSouP;
    nSubDet = PSeq[0][0].nSubDet;
    nSubSet = PSeq[0][0].nSubSet;

    dNx = Proj[0].NX;
    dNy = Proj[0].NY;
    dN = dNx * dNy;
    nIMG = Img.NX * Img.NY * Img.NZ;
    nRow = fmax(nIMG,dN);

    projTMP = (double*) malloc(sizeof(double)*nRow);
    projRatio = (double*) malloc(sizeof(double)*nRow);
    imgTMP = (double*) malloc(sizeof(double) *nRow);

    for(ii = 0; ii<nRow; ii++)
    {
        projTMP[ii] = 0;
        projRatio[ii] = 0;
        imgTMP[ii] = 0;
    }

    sprintf(fileName,"%s//source",DATA_DIR);

    fp1 = fopen(fileName,"rb");
    checkFile(fp1,fileName);
    fread(imgTMP,sizeof(double),nIMG,fp1);

    fclose(fp1);

    for(iSouP = 0; iSouP < nSouP; iSouP++)
    {
        for(iDet = 0; iDet < nDet; iDet++)
        {
           iProj = iSouP * nDet + iDet;

            for(iPix = 0; iPix < dN; iPix++)
            {
                Proj[iProj].detImage[iPix] = 0;
            }
        }
    }




    for(iIter = 0; iIter < nSubSet; iIter++)
    {
        nLoad = PSeq[iIter][0].nLoad;  // acculation load for each iter; each subset of OSEM is one iter;

        for(iLoad = 0; iLoad < nLoad; iLoad++)
        {
            if(1)
            {
                iSouP = PSeq[iIter][iLoad].iSouP;
                iDet = PSeq[iIter][iLoad].iDet;
                iSubDet = PSeq[iIter][iLoad].iSubDet;
                printf("proj: iSouP %d iDet %d iSubDet %d\n",iSouP,iDet,iSubDet);

                iSRF = iSouP * nDet * nSubDet + iDet * nSubDet + iSubDet;
                iProj = iSouP * nDet + iDet;

                nRow = fmax(Img.NX * Img.NY * Img.NZ, Proj[iProj].NX * Proj[iProj].NY);


                // forward Projection; A' * X = y; attention system repsonse function is saved as colom as detector index;
                //for(sn = 0; sn < nIMG; sn ++) imgTMP[sn] = Img.image[sn];

                TranA_multiply_x( &(Srf[iSRF].value[0]), &(Srf[iSRF].index[0]), imgTMP, projTMP, nRow );


                for(iPix = 0; iPix < dN; iPix++)
                {
                    Proj[iProj].detImage[iPix] += projTMP[iPix];
                }
            }
        }


    }


    for(iSouP = 0;iSouP < nSouP; iSouP++)
    {
        for(iDet = 0; iDet < nDet; iDet++)
        {

            ii = iSouP * nDet + iDet;
            sprintf(fileName,"%s//mproj_iDet%d_iPosi%d",&(Proj[0].projFolder[0]),iDet,iSouP);  // mean projection
            fp1 = fopen(fileName,"wb");
            fwrite( &( Proj[ii].detImage[0] ), sizeof(double), dN,fp1);
            fclose(fp1);


//            for(iPix = 0; iPix < dN; iPix++) // generate non-stationary poisson noise
//            {
//                jj = ii*dN + iPix;
//                Proj[ii].detImage[iPix] = (double) poidev( (float) Proj[ii].detImage[iPix], &(randSeed));
//            }
//
//            sprintf(fileName,"%s//proj_iDet%d_iPosi%d",&(Proj[0].projFolder[0]),iDet,iSouP);
//            fp1 = fopen(fileName,"wb");
//            fwrite( &( Proj[ii].detImage[0] ), sizeof(double), dN,fp1);
//            fclose(fp1);


        }

    }

    free(projTMP);
    free(projRatio);
    free(imgTMP);
    getchar();

}
