#include"headFile.h"

struct detector *detInit()
{
    FILE *fp1;
    FILE *fp;
    char fileName[1000];
    int maxlength, CharRetCd, arglen;
	char  oneline[256], TagValue[256], TagName[256];
	char* ptr2;
	char* ptr1;
	int nDet;
	int nDivide;
	int nSubDet;
	double dx;
	double dy;
	double dz;
	int NX;
	int NY;
	int NZ;
	double attenCoef;
    struct detector *det;
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
		if(strstr(TagName, "NZ"))	NZ=(int)atof(TagValue);
        if(strstr(TagName, "nSubDet"))   nSubDet=(int)atof(TagValue);
		if(strstr(TagName, "nDivide"))	 nDivide=(int)atof(TagValue);
		if(strstr(TagName, "dx"))	     dx=(double)atof(TagValue);
		if(strstr(TagName, "dy"))	     dy=(double)atof(TagValue);
		if(strstr(TagName, "dz"))	     dz=(double)atof(TagValue);
		if(strstr(TagName, "attenCoef"))	     attenCoef=(double)atof(TagValue);
    }

    printf("detector parameters\n");
    printf("nDet %d NX %d NY %d NZ %d nSubDet %d nDivide %d \n dx %f dy %f dz%f\n attenuate coeffecient %f \n",nDet,NX,NY,NZ,nSubDet,nDivide,dx,dy,dz,attenCoef);
    det = subDetInit(nDet,NX,NY,NZ,nSubDet,nDivide,dx,dy,dz,attenCoef);



    return(det);
}

