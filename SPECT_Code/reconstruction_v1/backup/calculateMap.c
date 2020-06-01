#include"headFile.h"


int calculateMap(int nThread)
{
    FILE *fp1;
    FILE *fp;
    char fileName[1000];
    int maxlength, CharRetCd, arglen;
	char  oneline[256], TagValue[256], TagName[256];
	char* ptr2;
	char* ptr1;
    int flag;
    int iMap;

    int iph;
    int N1,N2,N3,N4;
    int NPoint;

    double radius,Rlimit,openAngle,angleLimit,incidentAngle,length,chanLength;
    double muPlatinum;

    pthread_t *threadId;
    pthread_attr_t attr;

    size_t defStackSize;
    int threadIndex;
    void *state;



    maxlength=256;

    // read file setting files to get iMap value;
    sprintf(fileName,"%s/setting.txt",DATA_DIR);
    fp1=fopen(fileName,"r");

    checkFile(fp1,fileName);

    // import the value of iMap;
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

		if(strstr(TagName, "iMap"))	iMap=(int)atof(TagValue);

    }

    fclose(fp1);


    if(iMap==0) return(0);// iMap = 0; do not need to calcualated the map; 1 needs to calculat.


    pthread_attr_init(&attr);
    threadId = (pthread_t*) malloc( sizeof(pthread_t)*nThread);
    pthread_attr_setdetachstate(&attr,PTHREAD_CREATE_JOINABLE);


    for(threadIndex = 0; threadIndex < nThread ; threadIndex++)
    {
        flag = pthread_create( &threadId[threadIndex], &attr,subCalculateMap,(void *)threadIndex);
        //status = pthread_create(&threadId[iThread],&(attr),srfsForThread,(void *)iThread);
        if(flag)
        {
            printf("error: in thread creating in penetration map creation");
            exit(-1);
        }

    }

    for(threadIndex = 0; threadIndex < nThread ; threadIndex++)
    {
        flag = pthread_join( threadId[threadIndex], &state);
        //status = pthread_create(&threadId[iThread],&(attr),srfsForThread,(void *)iThread);
        if(flag)
        {
            printf("error: in thread joining error in penetration map creation");
            exit(-1);
        }

    }

    pthread_attr_destroy(&attr);
    //pthread_exit(NULL);
    return(1);






}
