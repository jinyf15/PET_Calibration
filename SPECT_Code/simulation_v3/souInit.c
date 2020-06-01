#include"headFile.h"

struct source *souInit()
{
    int iSou;
    FILE *fp1;
    FILE *fp;
    char fileName[1000];
    int maxlength, CharRetCd, arglen;
	char  oneline[256], TagValue[256], TagName[256];
	char* ptr2;
	char* ptr1;
	int nSou;
	double dx;
	double dy;
	double dz;
	int NX;
	int NY;
	int NZ;
	double attenCoef;
	int flag;
	struct source *sou;

    // read the parallel sequecne setting files
    maxlength = 256;

    sprintf(fileName,"%s/setting.txt",DATA_DIR);
    fp1 = fopen(fileName,"r");

    printf("setting file %s\n;press anykey to continue",fileName);getchar();

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

    printf("source parameter\n");

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

		if(strstr(TagName, "nSou"))	nSou = (int)atof(TagValue);
		if(strstr(TagName, "SNX"))   NX  = (int)atof(TagValue);
		if(strstr(TagName, "SNY"))	 NY  = (int)atof(TagValue);
		if(strstr(TagName, "SNZ"))	 NZ  = (int)atof(TagValue);
		if(strstr(TagName, "sdx"))	 dx  = (double)atof(TagValue);
		if(strstr(TagName, "sdy"))	 dy  = (double)atof(TagValue);
		if(strstr(TagName, "sdz"))	 dz  = (double)atof(TagValue);


    }
    fclose(fp1);

    printf("nSou %d NX %d NY %d NZ %d \n dx %f dy %f dz%f\n",nSou,NX,NY,NZ,dx,dy,dz);


    //sou = malloc(sizeof(struct source)*nSou); // allocate the memory for the source vector

    sou = subSouInit(nSou,NX,NY,NZ,dx,dy,dz);
    //sou = NULL;

    return(sou);


}
