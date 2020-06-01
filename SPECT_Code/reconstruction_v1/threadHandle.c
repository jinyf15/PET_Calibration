#include"headFile.h"
#include"globalVariable.h"

int threadHandle(int index)
{
    char fileName[1000];
    FILE *fp1;

    pthread_t threadID[MAX_CPU];
    pthread_attr_t attr;
    int nThread;
    int iThread;
    int status;
    void *tmp;
    int kk;


    nThread = PSeq[0][0].nThread;


    pthread_attr_init(&attr);
    pthread_attr_setdetachstate(&attr,PTHREAD_CREATE_JOINABLE);

    pthread_mutex_init(&comuMutex,NULL);
    pthread_cond_init(&sigalForCal,NULL);
    pthread_cond_init(&sigalForUpdate,NULL);

    iThread = nThread; // open thread for update;
    status = pthread_create(&(threadID[iThread]),&attr,reconUpdate,(void *)index);
    if(status)
    {
        printf("error in create thread %d\n",iThread);
        getchar();
        exit(-1);
    }

    // open thread for cal;
    for(iThread = 0; iThread < nThread; iThread++)
    {
        if(iThread<nThread)
        {
            status = pthread_create(&(threadID[iThread]),&attr,Recon_Cal,(void *) iThread);
            if(status)
            {
                printf("error in create thread %d\n",iThread);
                getchar();
                exit(-1);
            }
        }
    }

    // join all thread;
    pthread_attr_destroy(&attr);

    for(iThread = 0; iThread <= nThread; iThread++)
    {

        status = pthread_join((threadID[iThread]),&tmp);
        if(status)
        {
            printf("error in joining thread %d\n",iThread);
            getchar();
            exit(-1);
        }
    }


    pthread_mutex_destroy(&comuMutex);
    pthread_cond_destroy(&sigalForCal);
    pthread_cond_destroy(&sigalForUpdate);


    return(1);

}
