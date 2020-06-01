#include"headFile.h"

void resetImg(struct source *img, struct source *imgBuf)
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
	int nRead;




    (*img).imgReady = 1;

    (*imgBuf).imgReady = 0;

    NX = img[0].NX;
    NY = img[0].NY;
    NZ = img[0].NZ;

    sprintf(fileName,"%s//img0",DATA_DIR);
    fp1 = fopen(fileName,"rb");

    checkFile(fp1,fileName);

    nRead = fread(&( (*img).image[0] ), sizeof(double), NX*NY*NZ, fp1);
    fclose(fp1);

    if(nRead!=NX*NY*NZ)
    {
        printf("error in reading file from %s\n",fileName);
        getchar();
        getchar();
        exit(-1);
    }



    for(sn = 0; sn < NX*NY*NZ; sn++)
    {
        imgBuf[0].image[sn] = 0;
    }


}
