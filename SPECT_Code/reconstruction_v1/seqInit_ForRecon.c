#include "headFile.h"


struct parellSequence **seqInit_ForRecon(struct projection *proj)
{
    FILE *fp1;
    FILE *fp;
    char fileName[1000];
    int maxlength, CharRetCd, arglen;
	char filename1[1000], oneline[256], TagValue[256], TagName[256];
	char* ptr2;
	char* ptr1;
	int nDet;
	int nSouP;
	int nSubDet;
	int nThread;
	int nIter;
	int iThread,ii;
	int nSubSet_OSEM;
	int iSubSet;
	int nSRF;
	struct parellSequence **pSeq;
	int flag;
	int **osemPara;

    // read the parallel sequecne setting files
    maxlength = 256;


    nDet  = proj[0].nDet;
    nSubDet = proj[0].nSubDet;
    nSouP = proj[0].nPosi;
    nSRF = nDet*nSubDet*nSouP;


    osemPara = loadReconFile(nDet,nSubDet,nSouP,nSRF,&nIter,&nSubSet_OSEM,&nThread);

    pSeq = subSeqInit_ForRecon(nDet,nSouP,nSubDet,nIter,nThread,nSubSet_OSEM,osemPara);

    //printf("here4 %d %d\n",1,osemPara[0][0]);getchar();
    //printf("loading OSEM setting done, press enter key to continue\n");
    //getchar();




    return(pSeq);




}
