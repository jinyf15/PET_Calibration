#include "headFile.h"

struct source *senInit(struct parellSequence **PSeq,struct source *Sou)
{

    int iThread;
    int nThread;
    struct source *sen;
    unsigned long int sn;
    unsigned long int NIMG;

    nThread = PSeq[0][0].nThread;

    printf("senInit.c nThread %d \n",nThread);

    sen = (struct source*) malloc(sizeof(struct source)*nThread);

    NIMG = Sou[0].NX*Sou[0].NY*Sou[0].NZ;


    for(iThread = 0; iThread < nThread; iThread++)
    {
       sen[iThread].NX = Sou[0].NX;
       sen[iThread].NY = Sou[0].NY;
       sen[iThread].NZ = Sou[0].NZ;


       sen[iThread].image = (double*)malloc(sizeof(double)*NIMG);
       printf("iThread %d nThread %d NX NY NZ %d %d %d NIMG %d\n",iThread,nThread,sen[iThread].NX,sen[iThread].NY,sen[iThread].NZ,NIMG);
    }


    return(sen);
}
