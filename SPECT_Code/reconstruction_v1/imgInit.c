#include"headFile.h"

void imgInit(struct source *img, struct source *imgBuf)
{

    FILE *fp1;
    FILE *fp;
    char fileName[1000];
    int maxlength, CharRetCd, arglen;
	char  oneline[256], TagValue[256], TagName[256];
	char* ptr2;
	char* ptr1;


	int NX;
	int NY;
	int NZ;

	int flag;
	int sn;


    // read the source file
    maxlength = 256;

    sprintf(fileName,"%s/setting.txt",DATA_DIR);
    fp1 = fopen(fileName,"r");


    checkFile(fp1,fileName);

    flag = 0;

    while (1) // find the section of source setting in setting file
    {
		CharRetCd = fgets (oneline, maxlength, fp1);
        ptr1 = strstr(oneline,"source setting");

        if(ptr1)
        {
            flag = 1;
            break;
        }
    }

    if(flag == 0)
    {
        printf("error could find the source setting in %s\n",fileName);
        getchar();
        exit(-1);
    }


    while (1)
    {
		CharRetCd = fgets (oneline, maxlength, fp1);
		if(!CharRetCd) break;

		ptr1 = strstr(oneline,"source setting end");

		if(ptr1)
		{
            break;
		}


		ptr2 = oneline;
		ptr1 = strstr(oneline,"=");
		//printf("%s\n",oneline);getchar();

		if(ptr1)
		{
			ptr2 = ptr1+1;
			arglen = strlen (ptr2);
			if(arglen<=0 ||  arglen>=256) continue;
			*ptr1 = 0;
		}
		strcpy(TagName, oneline);

		strcpy(TagValue, ptr2);

		if(strstr(TagName, "SNX"))   NX  = (int)atof(TagValue);
		if(strstr(TagName, "SNY"))	 NY  = (int)atof(TagValue);
		if(strstr(TagName, "SNZ"))	 NZ  = (int)atof(TagValue);


    }
    fclose(fp1);

    (*img).dx = -1;
    (*img).dy = -1;
    (*img).dz = -1;
    (*img).iPosi = -1;
    (*img).nPosi = -1;
    (*img).NX = NX;
    (*img).NY = NY;
    (*img).NZ = NZ;

    (*img).imgReady = 1;

    (*imgBuf).dx = -1;
    (*imgBuf).dy = -1;
    (*imgBuf).dz = -1;
    (*imgBuf).iPosi = -1;
    (*imgBuf).nPosi = -1;
    (*imgBuf).NX = NX;
    (*imgBuf).NY = NY;
    (*imgBuf).NZ = NZ;
    (*imgBuf).imgReady = 0;

    (*img).image =  (double *) malloc(sizeof(double)*NX*NY*NZ);

    (*imgBuf).image =  (double *) malloc(sizeof(double)*NX*NY*NZ);

    sprintf(fileName,"%s//img0",DATA_DIR);
    fp1 = fopen(fileName,"rb");

    checkFile(fp1,fileName);

    fread(&( (*img).image[0] ), sizeof(double), NX*NY*NZ, fp1);
    fclose(fp1);

    for(sn = 0; sn < NX*NY*NZ; sn++)
    {
        imgBuf[0].image[sn] = 0;
    }


    (*img).souMap = NULL;
    (*imgBuf).souMap = NULL;

    souMaskInit(img,1,0);


    for(sn = 0; sn < NX*NY*NZ; sn++)
    {
        if(img[0].souMask[sn] == 0)
        {
            img[0].image[sn] = 0;
        }
    }

}
