#include "headFile.h"
#include "globalVariable.h"

void srfsForThread(void * ID)
{
    int iLoad;
    int iDet;
    int iSouP;
    int iIter;
    int iSubDet;
    int nLoad;
    int threadId;


    threadId = (int)ID;
    nLoad = PSeq[threadId][0].nLoad;

    printf("number of calculation load %d for thread %d \n",nLoad,threadId);
    for(iLoad =0;iLoad<nLoad;iLoad++)
    {
        iDet = PSeq[threadId][iLoad].iDet;
        iSouP = PSeq[threadId][iLoad].iSouP;
        iSubDet = PSeq[threadId][iLoad].iSubDet;
        printf("iThread %d iLoad %d iSoup %d iDet %d iSubDet %d\n",threadId,iLoad,iSouP,iDet,iSubDet);

        Srf[threadId].iDet = iDet;
        Srf[threadId].iSou = iSouP;
        Srf[threadId].iSubDet = iSubDet;
        Srf[threadId].iSRF =  PSeq[threadId][iLoad].iSeq;
        Srf[threadId].nSRF =  PSeq[threadId][iLoad].totalSeq; // in simulation nIter =1;
        //sleep(1);

        genSubSRF(Det,iDet,iSubDet,Sou,iSouP,Ph,Srf,threadId);
        //genSubSRF(Det,iDet,iSubDet,Sou,iSouP,Ph,Srf);
    }

    pthread_exit(NULL);
}
