#ifndef HEADFILE_H_INCLUDED
#define HEADFILE_H_INCLUDED

#include <stdio.h>
#include<stdlib.h>
#include<string.h>
#include <pthread.h>
#include <math.h>
#include <unistd.h>
#include <time.h>
#include "./system.h"
#include "./lib/cblas.h"
#define DATA_DIR "../data/"

#define PI 3.14159265359

// function declaration
struct parellSequence **seqInit(struct detector *det,struct source *sou);
struct parellSequence **subSeqInit(int nDet,int nSouP,int nSubDet,int nIter,int nThread);
struct parellSequence **seqInit_ForRecon(struct projection *proj);
int  **loadReconFile(int nDet,int nSubDet,int nSouP,int nSRF,int *nIter, int *nSubSet_OSEM,int *nThread);

struct projection *projInit();
struct projection *subProjInit(int nDet,int NX,int NY,int nSubDet,int nSou);
void loadProjection(struct projection *proj,int projIndex);
void loadDetPixMaskForRecon(struct projection *proj,int projIndex);

void imgInit(struct source *img, struct source *imgBuf);

struct sysSRF *loadSRF(struct projection *proj);
unsigned int srfread1(char *filename);
struct source *senInit_ForRecon(struct parellSequence **PSeq,struct source img);

void A_multiply_x(float sa[], unsigned int ija[], double x[], double b[],unsigned int n);
void TranA_multiply_x(float sa[], unsigned int ija[], double x[], double b[],unsigned int n);

void Recon_Cal(void * threadID);
void reconUpdate(void *index);

long *initRandSeed(struct projection *Proj);
float ran1(long *idum);
float poidev(float xm, long *idum);
float gammln(float xx);
void projGeneration(void * threadID); // generate the mean projection
void genNoiseProj(struct projection * Proj,int index); // generation noise project
void loadNoiseProj(struct projection * Proj,int index) ;; // load projection from file;
void addPulse2Proj(); // add pulse to projection;




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

int threadHandle(int index);
void resetSeq(struct parellSequence **seq);
void resetImg(struct source *img, struct source *imgBuf);




#endif // HEADFILE_H_INCLUDED
