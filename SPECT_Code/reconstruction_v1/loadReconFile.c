#include"headFile.h"

int  **loadReconFile(int nDet,int nSubDet,int nSouP,int nSRF,int *nIter, int *nSubSet_OSEM,int *nThread)
{
    FILE *fp1;
    FILE *fp;
    char fileName[1000];
    int maxlength, CharRetCd, arglen;
	char filename1[1000], oneline[256], TagValue[256], TagName[256];
	char* ptr2;
	char* ptr1;

	int iSubSet;
	int iDet;
	int iSouP;
	int iSubDet;
	int iSRF;
	int flag;
	int ii;

	int **osemPara;

    maxlength = 256;
     sprintf(fileName,"%s//setting.txt",DATA_DIR);
    fp1=fopen(fileName,"r");

    checkFile(fp1,fileName);

    while (1)
    {
		CharRetCd = fgets (oneline, maxlength, fp1);
		if(!CharRetCd) break;

        ptr1 = strstr(oneline,"reconstruction setting end");

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

		if(strstr(TagName, "nIter"))	*nIter = (int)atof(TagValue);
		if(strstr(TagName, "nSubSet_OSEM"))	 *nSubSet_OSEM= (int)atof(TagValue);
		if(strstr(TagName, "nThread"))	 *nThread= (int)atof(TagValue);


		ptr1 = strstr(oneline,"OSEM parameter start");
        if(ptr1)
		{
            printf("nSubSet_OSEM %d nSRF %d\n", *nSubSet_OSEM,nSRF);
            osemPara = (int **) malloc(sizeof(int*) * ( nSRF ));
            for(iSRF = 0; iSRF < nSRF; iSRF++)
            {
                osemPara[iSRF] = (int*) malloc(sizeof(int) * 5);
            }

            CharRetCd = fgets (oneline, maxlength, fp1);


            for(iSRF = 0 ;iSRF < nSRF; iSRF++)
            {
                fscanf(fp1,"%d\t %d\t %d\t %d\t %d\n",&(iSubSet),&(ii),&(iSouP),&(iDet),&(iSubDet));
                printf("iSubset %d ii %d iSou %d idet %d isubdet %d\n",iSubSet,ii,iSouP,iDet,iSubDet);

                osemPara[iSRF][0] = iSubSet;
                osemPara[iSRF][1] = ii;
                osemPara[iSRF][2] = iSouP;
                osemPara[iSRF][3] = iDet;
                osemPara[iSRF][4] = iSubDet;
            }
		}

    }

    fclose(fp1);

    return(osemPara);


}
