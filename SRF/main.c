#include "headFile.h"
#include "GlobalVariables.h"

int main(){
	printf("start!\n");
	fflush(stdout);
	pthread_t *threadId;
	pthread_attr_t attr;
	nThread = 4;
	int iThread;
	Det = detInit();
	printf("Initialize detectors successfully!\npress any key to continue\n");
	fflush(stdout);
	getchar();
	
	Sou = souInit();
	printf("Initialize source space successfully!\npress any key to continue\n");
	fflush(stdout);
	getchar();
	
	PET = petInit();
	printf("Initialize PET panels successfully!\npress any key to continue\n");
	fflush(stdout);
	getchar();


	pthread_attr_init(&attr);
    	pthread_attr_setdetachstate(&attr, PTHREAD_CREATE_JOINABLE);
	threadId = (pthread_t *)malloc(sizeof(pthread_t)*nThread);
	if (pthread_mutex_init(&mutexSRF,NULL) != 0){
		free(threadId);
		return 0;
	}	
	for (iThread = 0;iThread < nThread;iThread++){
		/*if (0){
			CPU_ZERO(&cpus);
			CPU_SET(iThread, cpus);
			pthread_attr_setaffinity_np(&attr, sizeof(cpu_set_t), &cpus);
		}*/
		
		printf("start multi threads %d!\n",iThread);
		fflush(stdout);
		
		if (pthread_create(&threadId[iThread], &(attr), CalSen, (void *)iThread)){
			printf("Thread%d create failed", iThread);
			return 0;
		}
	}
	for (iThread = 0;iThread < nThread;iThread++){
		pthread_join(threadId[iThread],NULL);
	}
	pthread_mutex_destroy(&mutexSRF);
	free(threadId);
	return 1;
}
