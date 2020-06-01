#include"headFile.h"

void loadPhAxis(int phIndex,double *phX,double *phY,double *phZ)
{
    FILE *fp1;
    FILE *fp;
    char fileName[1000];
    int maxlength, CharRetCd, arglen;
	char  oneline[256], TagValue[256], TagName[256];
	char* ptr2;
	char* ptr1;
	int flag;
	char CharTmp[256];
	char *pEnd;



    maxlength=256;

    sprintf(fileName,"%s/phMap/%d_profile.txt",DATA_DIR,phIndex);
    fp1=fopen(fileName,"rt");
    printf("fileName %s\n",fileName);

    checkFile(fp1,fileName);

    flag = 0;
    while (1)
    {
		CharRetCd = fgets (oneline, maxlength, fp1);
        ptr1 = strstr(oneline,"asix parameter");

        if(ptr1)
        {
            flag = 1;
            break;
        }
    }

    if(flag == 0)
    {
        printf("error could find in axis parameter %s\n",fileName);
        getchar();
        exit(-1);
    }

    CharRetCd = fgets (oneline, maxlength, fp1); // x axis
    phX[0] = strtod(&(oneline[1]),&pEnd);
    phX[1] = strtod(pEnd,&pEnd);
    phX[2] = strtod(pEnd,NULL);

    CharRetCd = fgets (oneline, maxlength, fp1); // y axis
    phY[0] = strtod(&(oneline[1]),&pEnd);
    phY[1] = strtod(pEnd,&pEnd);
    phY[2] = strtod(pEnd,NULL);

    CharRetCd = fgets (oneline, maxlength, fp1); // y axis
    phZ[0] = strtod(&(oneline[1]),&pEnd);
    phZ[1] = strtod(pEnd,&pEnd);
    phZ[2] = strtod(pEnd,NULL);




    fclose(fp1);
}
