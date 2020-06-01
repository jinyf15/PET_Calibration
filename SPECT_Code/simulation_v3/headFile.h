#ifndef HEADFILE_H_INCLUDED
#define HEADFILE_H_INCLUDED
#define _GNU_SOURCE

#include <stdio.h>
#include<stdlib.h>
#include<string.h>
#include <pthread.h>
#include<sched.h>
#include <math.h>
#include <time.h>
//#include "./nr/nrutil.h"
#include "system.h"
#include "./lib/cblas.h"

#define DATA_DIR "../data/"

#define PI 3.14159265359

// function declaration
struct parellSequence **seqInit(struct detector *det,struct source *sou);
struct parellSequence **subSeqInit(int nDet,int nSouP,int nSubDet,int nIter,int nThread);

struct detector *detInit();
struct detector *subDetInit(int nDet,int NX,int NY,int NZ,int nSubDet,int nDivide,double dx,double dy,double dz,double attenCoef);
void loadDetFile(struct detector *det,int nDet,int iDet);
void loadDetPixPos(struct detector *det,int nDet,int iDet);
void loadDetPixMask(struct detector *dp,int nDet,int iDet);
void phGlobalPosi(struct detector *det,int iDet,double *sysmPara);
void calDetGlobalPos(struct detector *det,double *sysPar,int iDet);

struct pinMap *phInit(int nThread);
void subPhInit(struct pinMap *ph,int iph);

struct source *souInit();
struct source *subSouInit(int nPos,int NX, int NY, int NZ, double dx, double dy,double dz);
void calSouGlobalPos(struct source *sou,int nPos,int iPos,int NX, int NY, int NZ, double dx, double dy,double dz);
void souMaskInit(struct source *sou,int nPos, int iPos);

struct sysSRF *sysSRFInit(struct parellSequence **pSeq);

void checkFile(FILE *fp,char *fileName);

void srfsForThread(void * ID);
double pixelprob_det2ML_sum_ph_insert2_double10(double *buf, double *dpos, double *subdpos, double *spos, double *map, int *N_step, double *cos_D, int DM_DET, double mu_det, int fn );
double pixelprob_det2ML_sum_ph_insert2_double11(double *buf, double *dpos, double *subdpos, double *spos, double *map, int *N_step, double *cos_D, int DM_DET, double mu_det, int fn );
void LoadPhAxis(int phIndex,double *phX,double *phY,double *phZ);
int genSubSRF(struct detector * det,int iDet,int iSubDet,struct source *sou,int iSouP,struct pinMap *ph,struct sysSRF *srf,int iThread);
struct source *senInit(struct parellSequence **PSeq,struct source *Sou);

int tranSRF(double *bufSen,int iThread);
int saveSRF();

int calculateMap(int nThread);
void subCalculateMap();
void loadPhFile(int phIndex,int *iph,int *N1,int *N2,int *N3,int *N4,int *NPoint,double *radius,double *Rlimit,double *openAngle,double *angleLimit,double *incidentAngle,double *length,double *chanLength, double *muPlatinum);
void loadPhAxis(int phIndex,double *phX,double *phY,double *phZ);
void loadPhProfile(int phIndex,int nPoints,double *phProfCenter,double *phProfDown,double *phProfUp);
void penetration_map_compound_eye(double *phX,double *phY, double *phZ,double *line_up,double *line_p, double *line_dp,int N11,double *dirt,int N1,int N2,int N3,int N4,int i_nph,int nph,double *phMap,double muPlatinum,int iPrint);
double dot(double*,double*);
void  cross(double*,double*,double*);
void  vector_ge(double*,double*,double*);
int IsInTriangle(double *,double*);
int IsInPolygon(int,double*,double*);
void print_float(int,double*);
void print_int(int Num,int *v);
void bubbleSort(double *, int ,int,int *);


#endif // HEADFILE_H_INCLUDED
