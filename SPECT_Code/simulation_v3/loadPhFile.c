#include"headFile.h"


void loadPhFile( int phIndex,int *iph,int *N1,int *N2,int *N3,int *N4,int *NPoint,double *radius,double *Rlimit,double *openAngle,double *angleLimit,double *incidentAngle,double *length,double *chanLength, double *muPlatinum )
{
    FILE *fp1;
    FILE *fp;
    char fileName[1000];
    int maxlength, CharRetCd, arglen;
	char  oneline[256], TagValue[256], TagName[256];
	char* ptr2;
	char* ptr1;



    maxlength=256;

    sprintf(fileName,"%s/phMap/%d_profile.txt",DATA_DIR,phIndex);
    fp1=fopen(fileName,"r");
    printf("fileName %s\n",fileName);

    checkFile(fp1,fileName);

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

        if(strstr(TagName, "iph"))	    *iph = (int)atof(TagValue);
        if(strstr(TagName, "radius"))	*radius = (double)atof(TagValue);
        if(strstr(TagName, "Rlimit"))	*Rlimit = (double)atof(TagValue);
        if(strstr(TagName, "openAngle"))	 *openAngle = (double)atof(TagValue);
        if(strstr(TagName, "angleLimit"))	 *angleLimit = (double)atof(TagValue);
        if(strstr(TagName, "incidentAngle")) *incidentAngle = (double)atof(TagValue);
        if(strstr(TagName, "length"))	     *length = (double)atof(TagValue);
        if(strstr(TagName, "chanLength"))	 *chanLength = (double)atof(TagValue);
        if(strstr(TagName, "muPlatinum"))	 *muPlatinum = (double)atof(TagValue);
        if(strstr(TagName, "N1"))	*N1=(int)atof(TagValue);
        if(strstr(TagName, "N2"))	*N2=(int)atof(TagValue);
        if(strstr(TagName, "N3"))	*N3=(int)atof(TagValue);
        if(strstr(TagName, "N4"))	*N4=(int)atof(TagValue);
        if(strstr(TagName, "NPoint"))	*NPoint=(int)atof(TagValue);

    }
    fclose(fp1);
}
