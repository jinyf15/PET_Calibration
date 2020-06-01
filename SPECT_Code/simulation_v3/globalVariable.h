#ifndef GLOBALVARIABLE_H_INCLUDED
#define GLOBALVARIABLE_H_INCLUDED
// global variable
struct detector *Det;
struct pinMap *Ph;
struct parellSequence **PSeq;
struct source *Sou;
struct sysSRF *Srf;
struct source *Sen;

// for parallell sequence

pthread_mutex_t mutexSRF;




#endif // GLOBALVARIABLE_H_INCLUDED
