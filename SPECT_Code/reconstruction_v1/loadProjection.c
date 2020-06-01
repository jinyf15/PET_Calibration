#include"headFile.h"

void loadProjection(struct projection *proj,int projIndex)
{
    int NX;
    int NY;
    int iSou;
    int iDet;
    double buf;
    int nRead;
    int iProj;
    char fileName[1000];
    FILE *fp1;

    int maxlength, CharRetCd, arglen;
	char  oneline[256], TagValue[256], TagName[256];
	char* ptr2;
	char* ptr1;
	char projFolder[1000];

	// get projection folder;
    maxlength = 256;
    sprintf(fileName,"%s/setting.txt",DATA_DIR);
    fp1=fopen(fileName,"r");

    checkFile(fp1,fileName);

    iProj = 0;

    while (1)
    {
		CharRetCd = fgets (oneline, maxlength, fp1);
		if(!CharRetCd) break;

		ptr2 = oneline;
		ptr1 = strstr(oneline,"=");

		if(ptr1)
		{
			ptr2 = ptr1+2;
			arglen = strlen (ptr2);
			if(arglen<=0 ||  arglen>=256) continue;
			*ptr1 = 0;
		}
		strcpy(TagName, oneline);
		strcpy(TagValue, ptr2);

		if(strstr(TagName, "projFolder"))
		{
            strncpy(&(proj[projIndex].projFolder[0]),TagValue,arglen-1);
            //printf("loading projection from %s\n",TagValue);
		}
		if(strstr(TagName, "iProj"))   iProj=(int)atof(TagValue);
    }

    fclose(fp1);
    NX = proj[projIndex].NX;
    NY = proj[projIndex].NY;
    proj[projIndex].iProj = iProj;


    proj[projIndex].detImage = (double *)malloc(sizeof(double)*NX*NY);

    iDet = proj[projIndex].iDet;
    iSou = proj[projIndex].iPosi;

    if(iProj == 0) // projection is generated load from disk
    {
        sprintf(fileName,"%s//proj_iDet%d_iPosi%d", &(proj[projIndex].projFolder[0]), iDet, iSou);


        printf("loading projection %s\n",fileName);

        fp1 = fopen(fileName,"rb");
        checkFile(fp1,fileName);

        nRead = fread(&(proj[projIndex].detImage[0]),sizeof(double),NX*NY,fp1);
        if(nRead != NX*NY)
        {
            printf("error % reading file from %s,number of data request %d,actualy read %d\n",fileName,NX*NY,nRead);
            getchar();getchar();
            exit(-1);
        }



        fclose(fp1);

    }



}
