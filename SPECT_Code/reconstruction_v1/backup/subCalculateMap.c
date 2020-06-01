#include"headFile.h"
#include"globalVariable.h"

void subCalculateMap(void *threadId)
{
    FILE *fp1;
    FILE *fp;
    char fileName[1000];
    int maxlength, CharRetCd, arglen;
	char  oneline[256], TagValue[256], TagName[256];
	char* ptr2;
	char* ptr1;
    int flag;
    int nph;
    int phIndex;
    int iph;
    int N1,N2,N3,N4;
    int NPoint;
    int threadIndex;
    int iPrint;
    double tmp;



    double radius,Rlimit,openAngle,angleLimit,incidentAngle,length,chanLength;
    double angleRa;
    double dirt[4];
    double *phProfDown,*phProfUp,*phProfCenter;
    double phX[3],phY[3],phZ[3];
    double muPlatinum;
    double *phMap;

    threadIndex = (int) threadId;

    maxlength=256;
    // read file setting files to get iMap value;
    sprintf(fileName,"%s/phMap/0_profile.txt",DATA_DIR);
    fp1=fopen(fileName,"r");
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

		if(strstr(TagName, "nph"))	nph=(int)atof(TagValue);
    }
    fclose(fp1);

    if(threadIndex==0) printf("number of phinhole for calculation %d\n",nph);

    for(phIndex=threadIndex; phIndex<nph;phIndex += PSeq[0][0].nThread)
    {

        // loading other parameters


        loadPhFile( phIndex, &iph, &N1, &N2, &N3, &N4, &NPoint, &radius, &Rlimit, &openAngle, &angleLimit, &incidentAngle, &length, &chanLength, &muPlatinum);
        loadPhAxis( phIndex, phX, phY, phZ );

        phProfUp = (double*) malloc( sizeof (double) * 3 * NPoint );
        phProfDown = (double*) malloc( sizeof (double) * 3 * NPoint );
        phProfCenter = (double*) malloc( sizeof (double) * 3 * NPoint );

        loadPhProfile( phIndex, NPoint, phProfCenter, phProfDown, phProfUp );


        phMap = ( double* ) malloc( sizeof( double ) * N1 * N2 * N3 * N4 );

        dirt[ 0 ] = Rlimit / N1;
        dirt[ 1 ] = Rlimit / N2;
        dirt[ 2 ] = angleLimit * PI / 180.0 / N3;
        dirt[ 3 ] = 2 * PI / N4;

        if(threadIndex==0)
        {
            printf("iph %d : radius %f Rlimit %f openAngle %f \n angleLimit %f incident angle %f length %f\nchanLength %f muPlatinum %f\n",iph,radius,Rlimit,openAngle,angleLimit,incidentAngle,length,chanLength,muPlatinum);
            printf("N1 %d N2 %d N3 %d N4 %d NPoint %d \n",N1,N2,N3,N4,NPoint);
            printf("phX %f %f %f\n",phX[0],phX[1],phX[2]);
            printf("phY %f %f %f\n",phY[0],phY[1],phY[2]);
            printf("phZ %f %f %f\n",phZ[0],phZ[1],phZ[2]);

            printf("dirt0 %f dirt1 %f dirt2 %f dirt3 %f\n",dirt[0],dirt[1],dirt[2],dirt[3]);

        }

        if(threadIndex==0) iPrint = 1;

        else iPrint = 0;

        penetration_map_compound_eye (phX, phY, phZ, phProfUp, phProfCenter, phProfDown, (NPoint+3)/4, dirt, N1, N2, N3, N4, phIndex, nph, phMap, muPlatinum,iPrint);


        sprintf(fileName,"%s//phMap//map%d",DATA_DIR,phIndex);
        printf("fileName %s\n",fileName);
        fp1 = fopen(fileName,"wb");
        checkFile(fp1,fileName);

        tmp = (double) nph;
        fwrite(&tmp,sizeof(double),1,fp1);

        tmp = (double) phIndex;
        fwrite(&tmp,sizeof(double),1,fp1);

        fwrite(&radius,sizeof(double),1,fp1);
        fwrite(&Rlimit,sizeof(double),1,fp1);

        angleRa = openAngle*PI/180;
        fwrite(&angleRa,sizeof(double),1,fp1);

        angleRa = angleLimit*PI/180;
        fwrite(&angleRa,sizeof(double),1,fp1);

        angleRa = (incidentAngle*PI/180);
        fwrite(&angleRa,sizeof(double),1,fp1);

        fwrite(&length,sizeof(double),1,fp1);
        fwrite(&chanLength,sizeof(double),1,fp1);
        fwrite(&length,sizeof(double),1,fp1);

        tmp = (double) N1;
        fwrite(&tmp,sizeof(double),1,fp1);

        tmp = (double) N2;
        fwrite(&tmp,sizeof(double),1,fp1);

        tmp = (double) N3;
        fwrite(&tmp,sizeof(double),1,fp1);

        tmp = (double) N4;
        fwrite(&tmp,sizeof(double),1,fp1);

        fwrite(phMap,sizeof(double),N1*N2*N3*N4,fp1);

        fwrite(dirt,sizeof(double),4,fp1);

        fclose(fp1);


        free(phMap);
        free(phProfCenter);
        free(phProfDown);
        free(phProfUp);

    }

    pthread_exit(NULL);


}
