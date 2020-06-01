#include "headFile.h"
#include "globalVariable.h"

void Recon_Cal(void * threadID)
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
    int ii;




    char fileName[1000];
    FILE *fp1;

    iThread = (int) threadID;

    printf("thread ID %d\n",iThread);

    nIter = PSeq[0][0].nIter;
    nDet = PSeq[0][0].nDet;
    nSouP = PSeq[0][0].nSouP;
    nSubDet =PSeq[0][0].nSubDet;

    dNx = Proj[0].NX;
    dNy = Proj[0].NY;
    dN = dNx * dNy;
    nIMG = Img.NX * Img.NY * Img.NZ;;
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

    for(iIter = 0; iIter < nIter; iIter++)
    {
        nLoad = PSeq[iIter][0].nLoad;  //  load for each iter; each subset of OSEM is one iter;
        for(iLoad = 0; iLoad < nLoad; iLoad++)
        {
            if(iThread == PSeq[iIter][iLoad].iThread)
            {
                iSouP = PSeq[iIter][iLoad].iSouP;
                iDet = PSeq[iIter][iLoad].iDet;
                iSubDet = PSeq[iIter][iLoad].iSubDet;
                //printf("recon: iSouP %d iDet %d iSubDet %d, iSubDet %d\n",iSouP,iDet,iSubDet);

                iSRF = iSouP * nDet * nSubDet + iDet * nSubDet + iSubDet;
                iProj = iSouP * nDet + iDet;


                nRow = fmax(Img.NX * Img.NY * Img.NZ, Proj[iProj].NX * Proj[iProj].NY);


                // forward Projection; A' * X = y; attention system repsonse function is saved as colom as detector index;
                for(sn = 0; sn < nIMG; sn ++) imgTMP[sn] = Img.image[sn];

                TranA_multiply_x( &(Srf[iSRF].value[0]), &(Srf[iSRF].index[0]), imgTMP, projTMP, nRow );

                // calucation of ratio
                for(iPix = 0; iPix < dN; iPix++)
                {
                    if(iSubDet = Proj[iProj].detPixMask[iPix])
                    {
                        if(projTMP[iPix] > 0)
                        {
                            projRatio[iPix] = Proj[iProj].detImage[iPix] / projTMP[iPix];

                        }

                    }

                    else
                    {
                        projRatio[iPix] = 1; // there is problem here to handle bad pixels;
                    }
                }

                // backProjection of the ratio;

                A_multiply_x( &(Srf[iSRF].value[0]), &(Srf[iSRF].index[0]), projRatio, imgTMP, nRow );

                pthread_mutex_lock(&comuMutex); // update Imgbuf; needs to lock it;

                for(sn = 0; sn <nIMG; sn++)
                {
                    ImgBuf.image[sn] += imgTMP[sn];
                }

                PSeq[iIter][0].nFinished++; // to record how many load has been finished;

                if(PSeq[iIter][0].nFinished == PSeq[iIter][0].nLoad)
                {
                    Img.imgReady = 0;                        // after finished all the work load, reset imgReady for update;
                    ImgBuf.imgReady = 1;
                }

                pthread_mutex_unlock(&comuMutex);
//
            }
        }


        while(iIter<nIter && (PSeq[iIter][0].nFinished < nLoad||Img.imgReady != 1) ) // waiting for other threads finishing working;
        {


            usleep(10000);
            //printf("wait for other thread finished cal iter %d ithread %d\n",iIter,threadID);

        }

    }


    free(projTMP);
    free(projRatio);
    free(imgTMP);


    pthread_exit(NULL);


}
