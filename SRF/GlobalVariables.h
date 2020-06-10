#ifndef GLOBALVARIABLES_H_INCLUDED
#define GLOBALVARIABLES_H_INCLUDED
#include "pthread.h"
// global variable
struct detector *Det;
struct source Sou;
struct PETpanel *PET;
int nThread;
// for parallell sequence

pthread_mutex_t mutexSRF;




#endif // GLOBALVARIABLE_H_INCLUDED
