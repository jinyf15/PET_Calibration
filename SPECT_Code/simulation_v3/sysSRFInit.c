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

    FILE *fp1;
    FILE *fp;
    char fileName[1000];
    int maxlength, CharRetCd, arglen;
	char  oneline[256], TagValue[256], TagName[256];
	char* ptr2;
	char* ptr1;
	int loadMemory;

    maxlength=256;

    // read file setting files to get nMax value;
    sprintf(fileName,"%s/setting.txt",DATA_DIR);
    fp1=fopen(fileName,"r");

    checkFile(fp1,fileName);

    // import the value of iMap;
    while (1)
    {
		CharRetCd = fgets (oneline, maxlength, fp1);
		if(!CharRetCd) break;

		ptr2 = oneline;
		ptr1 = strstr(oneline,"=");

		if(ptr1)
		{
			ptr2 = ptr1+1;
			arglen = strlen (ptr2);
			if(arglen<=0 ||  arglen>=256) continue;
			*ptr1 = 0;
		}
		strcpy(TagName, oneline);
		strcpy(TagValue, ptr2);

		if(strstr(TagName, "nMax"))	nMax=(unsigned long int)atof(TagValue);

    }

    fclose(fp1);
    printf("maxium number of non zero elecment srf could be stored nMax %ld \n",nMax);


    nSrf = pSeq[0][0].nDet*pSeq[0][0].nSouP*pSeq[0][0].nSubDet;
    nThread = pSeq[0][0].nThread;

    loadMemory = 4*nThread*nMax*sizeof(float)/(1.e9);


    printf("attention!!! the code is going to load memory  %d GB\n system maxium memory %d \n ",loadMemory,PC_MAX_MEMORY);
    if(loadMemory>=PC_MAX_MEMORY){
        printf("the the loading memory of program is larger than system memory.\n program will exist\n");
        printf("please change the nMax in the file setting.txt");
        getchar();
        exit(-1);

    }

    printf("press anykey to continue\n");
    getchar();


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
