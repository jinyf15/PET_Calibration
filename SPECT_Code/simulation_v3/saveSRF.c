#include "headFile.h"
#include "globalVariable.h"

int saveSRF(){

    int nThread;
    unsigned int mSize;
    unsigned int tmp;
    unsigned int kk;
    unsigned int NIMG;
    int *finishLoad;
    int ii;
    int iSouce;
    int Ifinish;
    int iThread;
    char fileName[1000];
    FILE *fp1;
    int iSRF,iSou,iDet,iSubDet;
    int status;

    double *senBuf;


    nThread = PSeq[iThread][0].nThread;
    printf("nThread %d\n",nThread);


    finishLoad = (int *) malloc(sizeof(int)*nThread);

    NIMG = Sou[0].NX * Sou[0].NY * Sou[0].NZ;
    //printf("source NIMg %d\n",NIMG);getchar();
    senBuf = (double *) malloc( sizeof(double) * NIMG);
    for(ii = 0 ;ii < NIMG; ii ++) senBuf[ii] = 0;

    for(ii=0;ii<nThread;ii++){
        finishLoad[ii] = 0;
    }

    while(1)
    {         // whether work finished or not?
        Ifinish = 1;

        for(ii=0;ii<nThread;ii++)
        {

            if(finishLoad[ii]<PSeq[ii][0].nLoad)
            {
                Ifinish = 0;
            }
        }

        if(Ifinish ==1) break; // finished, break loop;

        for(iThread=0; iThread < nThread; iThread++)
        {

            if(finishLoad[iThread] == PSeq[iThread][0].nLoad)
            {
                printf("inside finished??\n");
                continue; // the work for iThread is finished;

            }

            printf("thread %d is waiting for save\n",iThread);

             while(Srf[iThread+nThread].readySave <= 0.5)
             {

                sleep(10);

             }

            mSize = Srf[iThread+nThread].mSize; // size of the matrix;
            iSRF = Srf[iThread+nThread].iSRF;
            iSou = Srf[iThread+nThread].iSou;
            iDet = Srf[iThread+nThread].iDet;
            iSubDet = Srf[iThread+nThread].iSubDet;


            sprintf(fileName,"%s//srf//sat_iSou%d_iDet%d_iSubDet%d",DATA_DIR,iSou,iDet,iSubDet);
            printf("save srf fileName sat  %s\n",fileName);

            fp1 = fopen(fileName,"wb");
            checkFile(fp1,fileName);

            status = fwrite(&(Srf[iThread +nThread].value[1]),sizeof(float),mSize,fp1); // index for value and index vector are start from 1;

            if(status!=mSize)
            {
                printf("error in saveSRF for %s",fileName);
                getchar();
            }

            fclose(fp1);

            sprintf(fileName,"%s//srf//ijat_iSou%d_iDet%d_iSubDet%d",DATA_DIR,iSou,iDet,iSubDet);
            fp1 = fopen(fileName,"wb");
            checkFile(fp1,fileName);

            status = fwrite(&(Srf[iThread +nThread].index[1]),sizeof(unsigned int),mSize,fp1);

            if(status!=mSize)
            {
                printf("error in saveSRF for %s",fileName);
                getchar();getchar();
                exit(-1);
            }

            fclose(fp1);
        //
            NIMG = Sen[iThread].NX * Sen[iThread].NY *Sen[iThread].NZ;

            sprintf(fileName,"%s//srf//sen_iSou%d_iDet%d_iSubDet%d",DATA_DIR,iSou,iDet,iSubDet);
            fp1 = fopen(fileName,"wb");
            checkFile(fp1,fileName);

            status = fwrite(&(Sen[iThread].image[0]),sizeof(double),NIMG,fp1);


            if(status!=NIMG)
            {
                printf("error in saveSRF for %s",fileName);
                getchar();
            }

            fclose(fp1);

            for(ii = 0; ii < NIMG ; ii++) senBuf[ii] += Sen[iThread].image[ii];



            ++finishLoad[iThread];


            pthread_mutex_lock(&mutexSRF); // lock variable;

            Srf[iThread+nThread].readyTran = 1;
            Srf[iThread+nThread].readySave = 0;

            pthread_mutex_unlock(&mutexSRF);


        } // end of for iThread loop;

    }


    sprintf(fileName,"%s//srf//sen",DATA_DIR);
    fp1 = fopen(fileName,"wb");
    checkFile(fp1,fileName);

    status = fwrite(&(senBuf[0]),sizeof(double),NIMG,fp1);


    if(status!=NIMG)
    {
        printf("error in saveSRF for %s",fileName);
        getchar();
    }

    fclose(fp1);



}
