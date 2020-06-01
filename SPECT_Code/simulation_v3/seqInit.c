#include "headFile.h"


struct parellSequence **seqInit(struct detector *det,struct source *sou)
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
	struct parellSequence **pSeq;
	int flag;

    // read the parallel sequecne setting files
    maxlength = 256;


    nDet  = det[0].nDet;
    nSubDet = det[0].nSubDet;
    nSouP = sou[0].nPosi;
    nIter = 1; // for simulation it set to be 1;





    printf("parallel sequence setting\n");

    sprintf(fileName,"%s//setting.txt",DATA_DIR);
    fp1=fopen(fileName,"r");


    if(fp1==NULL)
    {
        printf("%s could NOT be found\n",fileName);
        getchar();
    }

    flag = 0;
    while (1) // find the section of source setting in setting file
    {
		CharRetCd = fgets (oneline, maxlength, fp1);
        ptr1 = strstr(oneline,"parallel setting");

        if(ptr1)
        {
            flag = 1;
            break;
        }
    }

    if(flag == 0)
    {
        printf("error could find the parallel setting in %s\n",fileName);
        getchar();
        exit(-1);
    }

    while (1)
    {
		CharRetCd = fgets (oneline, maxlength, fp1);
		if(!CharRetCd) break;

        ptr1 = strstr(oneline,"parallel setting end");

		if(ptr1)
		{
            break;
		}

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

//		if(strstr(TagName, "nDet"))	 nDet  = (int)atof(TagValue);
//		if(strstr(TagName, "nSouP")) nSouP = (int)atof(TagValue);
//		if(strstr(TagName, "nSubDet"))	nSubDet = (int)atof(TagValue);

		if(strstr(TagName, "nThread"))	nThread = (int)atof(TagValue);

//		if(strstr(TagName, "nIter"))	nIter = (int)atof(TagValue);
    }
   // read file ends

    printf("nDet%d nSoup%d nSubDet%d nThread%d nIter%d\n",nDet,nSouP,nSubDet,nThread,nIter);
    pSeq = subSeqInit(nDet,nSouP,nSubDet,nIter,nThread); // initiliaze the sequence vector;

    printf("press anykey to continue\n");
    getchar();

    return(pSeq);






}
