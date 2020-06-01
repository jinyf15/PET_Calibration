#include"headFile.h"

struct sysSRF *sysSRFInit(struct parellSequence **pSeq)
{
    int nSrf;
    int iDet;
    int nDet;
    int iSouP;
    int nSouP;
    int iSubDet;
    int nSubDet;
    int iSrf;
    int nThread;
    unsigned long int  nMax;
    int checkMemory;
    struct sysSRF *srf;


    printf("nDet %d nSoup %d pSeq %d\n",pSeq[0][0].nDet,pSeq[0][0].nSouP,pSeq[0][0].nThread);getchar();
    nSrf = pSeq[0][0].nDet*pSeq[0][0].nSouP*pSeq[0][0].nSubDet;

    printf("attention!!! the maxium memory of system that could be used is  %d GB\n ",PC_MAX_MEMORY);
    printf("if smaller than this value, by press number: 1 \n else press: 0\n ");
    scanf("%d",&checkMemory);

    if(checkMemory==1)
    {
        printf("the the loading memory of program is larger than system memory.\n program will exist\n");
        printf("please change the maxium system memory in head file <system.h>");
        getchar();
        exit(-1);
    }


    nMax = PC_MAX_MEMORY*pow(10,9)/nSrf/4/2; // wihtou rebin process;
    printf("maxium number of non zero elecment srf could be stored nMax %ld %d %ld\n",nMax,nSrf);




    nThread = pSeq[0][0].nThread;
    printf("number of thread %d\n",nThread);

    if(nThread*2>nSrf) // why need 2 times smaller than nSrf; check the line of allocate memomory for srf;
    {
        printf("error! the nthread is not two time smaller than nSrf");
        getchar();
    }

    srf = (struct sysSRF *)malloc(sizeof(struct sysSRF)*nThread*2); // vector size of srfs is two times of nthread; index nthread to 2*nthread is buffer srf for save srf;

    for(iSrf=0;iSrf< nThread;iSrf++)
    {
        srf[iSrf].nMax = nMax;
        srf[iSrf].mSize = nMax;
        srf[iSrf].value = (float *)malloc(sizeof(float)*nMax);
        srf[iSrf].index = (unsigned int *)malloc(sizeof(unsigned int)*nMax);

        srf[iSrf+nThread].nMax = nMax;
        srf[iSrf+nThread].mSize = nMax;
        srf[iSrf+nThread].value = (float *)malloc(sizeof(float)*nMax);
        srf[iSrf+nThread].index = (unsigned int *)malloc(sizeof(unsigned int)*nMax);
        srf[iSrf+nThread].readyTran = 1;
        srf[iSrf+nThread].readySave = 0;
    }


return(srf);

}
