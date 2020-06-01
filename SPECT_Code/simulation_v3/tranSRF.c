#include "headFile.h"
#include "globalVariable.h"

int tranSRF(double *bufSen,int iThread){

    int nThread;
    unsigned  int mSize;
    unsigned int tmp;
    unsigned int kk;
    unsigned int sn;
    unsigned int NIMG;



    nThread = PSeq[iThread][0].nThread;

    while(Srf[iThread+nThread].readyTran <= 0.5){

        printf("thread %d is waiting for transfer\n",iThread);
        sleep(10);
    }

    tmp = Srf[iThread].index[1];
    mSize = Srf[iThread].index[tmp-1] - 1; // size of the matrix;
    printf("ithread %d, transfer response function, size mSize %d\n",iThread,mSize);
    for(kk=1; kk<=mSize;kk++){

        Srf[iThread +nThread].value[kk] =  Srf[iThread].value[kk];
        Srf[iThread +nThread].index[kk] =  Srf[iThread].index[kk];
    }

    NIMG = Sen[iThread].NX * Sen[iThread].NY *Sen[iThread].NZ;


    for(sn = 0; sn< NIMG; sn++){
        Sen[iThread].image[sn] = bufSen[sn];
    }



    Srf[iThread+nThread].iDet = Srf[iThread].iDet;
    Srf[iThread+nThread].iSou = Srf[iThread].iSou;
    Srf[iThread+nThread].iSubDet = Srf[iThread].iSubDet;
    Srf[iThread+nThread].iSRF = Srf[iThread].iSRF;
    Srf[iThread+nThread].nSRF = Srf[iThread].nSRF;
    Srf[iThread+nThread].mSize = mSize;


    pthread_mutex_lock(&mutexSRF); // lock variable;

    Srf[iThread+nThread].readyTran = 0;
    Srf[iThread+nThread].readySave = 1;

    pthread_mutex_unlock(&mutexSRF);


    return(1);

}
