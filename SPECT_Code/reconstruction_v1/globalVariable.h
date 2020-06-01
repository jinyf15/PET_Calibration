#ifndef GLOBALVARIABLE_H_INCLUDED
#define GLOBALVARIABLE_H_INCLUDED
// global variable
long randSeed;
struct projection *Proj;

struct parellSequence **PSeq;
struct source Img;
struct source ImgBuf;
struct sysSRF *Srf;
struct source *Sen;

// for parallell sequence

pthread_mutex_t comuMutex;
pthread_cond_t sigalForUpdate;
pthread_cond_t sigalForCal;




#endif // GLOBALVARIABLE_H_INCLUDED
