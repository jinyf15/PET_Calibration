#include "headFile.h"
#include "globalVariable.h"


int main()
{

    pthread_t *threadId;
    pthread_attr_t attr;
    int nThread;
    int iThread;
    int status;
    cpu_set_t cpus;
    void *state;



    // ininitilized source
    Sou = souInit();

    // ininitilized detector
    Det = detInit();

    // sequence initialized
    PSeq = seqInit(Det,Sou);

    // srf initialized
    Srf = sysSRFInit(PSeq);

    // system sensitivity;
    Sen = senInit(PSeq,Sou); // Sen stores sensitivity of each point;

    //nitilized pinhole
    Ph = phInit(PSeq[0][0].nThread);


    // start parallel thread
    nThread =  PSeq[0][0].nThread;
    threadId = (pthread_t *)malloc(sizeof(pthread_t)*(nThread+1));

    pthread_mutex_init(&mutexSRF,NULL); // mutesSRF is for talking between saveSRF thread and calculating SRF thread

    pthread_attr_init(&attr);
    pthread_attr_setdetachstate(&attr,PTHREAD_CREATE_JOINABLE); // set each thread joinable; so that main.c could wait for all the thread finish workload;

    for(iThread = 0; iThread < nThread + 1; iThread++)
    {
        if(0)
        {
            CPU_ZERO(&cpus); // specifiy the working in each cpu;
            CPU_SET(iThread,&cpus);
            pthread_attr_setaffinity_np(&attr, sizeof(cpu_set_t), &cpus);
        }

        if(iThread<nThread)
        status = pthread_create(&threadId[iThread],&(attr),srfsForThread,(void *)iThread);
        else
        status = pthread_create(&threadId[iThread],&(attr),saveSRF,(void *)iThread);


        if(status)
        {
            printf("ERRROR; return code from pthread_creat() is %d \n",status);
        }
    }

    // waiting all thread to be finished
    for(iThread=0; iThread < nThread +1; iThread++)
    {
        status = pthread_join(threadId[iThread],&state);

        if(status)
        {
            printf("ithread %d ERRROR; return code from pthread_join() is %d \n",iThread,status);

        }
    }

    printf("thread joined\n");

    pthread_attr_destroy(&attr);
    pthread_mutex_destroy(&mutexSRF);

    free(Srf);
    free(Sou);
    free(Det);
    free(Ph);
    free(PSeq);
    free(Sen);

    printf("memory free ok end\n");
    pthread_exit(NULL);


}
