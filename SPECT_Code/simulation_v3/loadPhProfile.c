#include"headFile.h"

void loadPhProfile(int phIndex,int nPoints,double *phProfCenter,double *phProfDown,double *phProfUp)
{
    FILE *fp1;
    FILE *fp;
    char fileName[1000];
    int maxlength, CharRetCd, arglen;
	char  oneline[256], TagValue[256], TagName[256];
	char* ptr2;
	char* ptr1;
	int flag;
	char *pEnd;
	int ii,pointIndex;



    maxlength=256;

    sprintf(fileName,"%s/phMap/%d_profile.txt",DATA_DIR,phIndex);
    fp1=fopen(fileName,"rt");
    printf("fileName %s\n",fileName);

    checkFile(fp1,fileName);


    // file line profile close to object;
    flag = 0;
    while (1)
    {
		CharRetCd = fgets (oneline, maxlength, fp1);
        ptr1 = strstr(oneline,"profile1");

        if(ptr1)
        {
            flag = 1;
            break;
        }
    }

    if(flag == 0)
    {
        printf("error could find in ph profile1 %s\n",fileName);
        getchar();
        exit(-1);
    }

    for(ii = 0;ii<nPoints;ii++)
    {
        CharRetCd = fgets(oneline, maxlength, fp1);
        pointIndex = (int)strtod(oneline,&pEnd);
        if (ii != pointIndex)
        {
            printf("error in read ph profile ii %d pointInde %d, those two value should be equal",ii,pointIndex);
            getchar();
            exit(-1);
        }

        phProfUp[3*ii] = strtod(pEnd, &pEnd);
        phProfUp[3*ii + 1] = strtod(pEnd, &pEnd);
        phProfUp[3*ii + 2] = strtod(pEnd, NULL);
        //printf("phProfUp %d %f %f %f\n",ii,phProfUp[3*ii], phProfUp[3*ii+1] ,phProfUp[3*ii+2]); getchar();

    }

    // read pinhole line
    flag = 0;
    while (1)
    {
		CharRetCd = fgets (oneline, maxlength, fp1);
        ptr1 = strstr(oneline,"profile2");

        if(ptr1)
        {
            flag = 1;
            break;
        }
    }

    if(flag == 0)
    {
        printf("error could find in ph profile2 %s\n",fileName);
        getchar();
        exit(-1);
    }

    for(ii = 0;ii<nPoints;ii++)
    {
        CharRetCd = fgets(oneline, maxlength, fp1);
        pointIndex = (int)strtod(oneline,&pEnd);
        if (ii != pointIndex)
        {
            printf("error in read ph profile ii %d pointInde %d, those two value should be equal",ii,pointIndex);
            getchar();
            exit(-1);
        }

        phProfCenter[3*ii] = strtod(pEnd, &pEnd);
        phProfCenter[3*ii + 1] = strtod(pEnd, &pEnd);
        phProfCenter[3*ii + 2] = strtod(pEnd, NULL);
        //printf("phProfcenter %d %f %f %f\n",ii,phProfCenter[3*ii], phProfCenter[3*ii+1] ,phProfCenter[3*ii+2]); getchar();

    }

     // read profile close to detector
    flag = 0;
    while (1)
    {
		CharRetCd = fgets (oneline, maxlength, fp1);
        ptr1 = strstr(oneline,"profile3");

        if(ptr1)
        {
            flag = 1;
            break;
        }
    }

    if(flag == 0)
    {
        printf("error could find in ph profile3 %s\n",fileName);
        getchar();
        exit(-1);
    }

    for(ii = 0;ii<nPoints;ii++)
    {
        CharRetCd = fgets(oneline, maxlength, fp1);
        pointIndex = (int)strtod(oneline,&pEnd);
        if (ii != pointIndex)
        {
            printf("error in read ph profile ii %d pointInde %d, those two value should be equal",ii,pointIndex);
            getchar();
            exit(-1);
        }

        phProfDown[3*ii] = strtod(pEnd, &pEnd);
        phProfDown[3*ii + 1] = strtod(pEnd, &pEnd);
        phProfDown[3*ii + 2] = strtod(pEnd, NULL);
        //printf("phProfdown %d %f %f %f\n",ii,phProfDown[3*ii], phProfDown[3*ii+1] ,phProfDown[3*ii+2]); getchar();

    }



    fclose(fp1);
}
