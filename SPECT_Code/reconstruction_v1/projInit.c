#include"headFile.h"

struct projection *projInit()
{
    FILE *fp1;
    FILE *fp;
    char fileName[1000];
    int maxlength, CharRetCd, arglen;
	char  oneline[256], TagValue[256], TagName[256];
	char* ptr2;
	char* ptr1;
	int nDet;
	int nSubDet;

	int NX;
	int NY;

	int nSou;

    struct projection *proj;
    int flag;

    // read the parallel sequecne setting files
    maxlength=256;

    sprintf(fileName,"%s/setting.txt",DATA_DIR);
    fp1=fopen(fileName,"r");

    checkFile(fp1,fileName);

    flag = 0;
    while (1)
    {
		CharRetCd = fgets (oneline, maxlength, fp1);
        ptr1 = strstr(oneline,"detector setting");

        if(ptr1)
        {
            flag = 1;
            break;
        }
    }

    if(flag == 0)
    {
        printf("error could find the detector setting in %s\n",fileName);
        getchar();
        exit(-1);
    }


    while (1)
    {
		CharRetCd = fgets (oneline, maxlength, fp1);
		if(!CharRetCd) break;

		ptr1 = strstr(oneline,"detector setting end");
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

		if(strstr(TagName, "nDet"))	nDet=(int)atof(TagValue);
		if(strstr(TagName, "NX"))   NX=(int)atof(TagValue);
		if(strstr(TagName, "NY"))	NY=(int)atof(TagValue);
        if(strstr(TagName, "nSubDet"))   nSubDet=(int)atof(TagValue);

    }

    // read numer of source position;
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

		if(strstr(TagName, "nSou"))	nSou=(int)atof(TagValue);
    }

    fclose(fp1);


    printf("projection parameters\n");
    printf("nDet %d NX %d NY %d nSubDet %d nSou %d\n",nDet,NX,NY,nSubDet,nSou);
    proj = subProjInit(nDet,NX,NY,nSubDet,nSou);

    printf("press enter key to continue\n");getchar();


    return(proj);
}

