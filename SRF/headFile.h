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
#include "system.h"
#include "./lib/cblas.h"

#define DATA_DIR "../data/"

#define PI 3.14159265359

// function declaration

struct detector *detInit();
struct source souInit();
struct PETpanel *petInit();
void checkFile(FILE *fp,char *fileName);

void *CalSen(void *iThread);
double fabs(double a);
double dotproduct(double *a, double *b);
double *scalevector(double a, double *b);
double *minusvector(double *a, double *b);
int checkinrectangle(double *x0, double *x1, double *x2, double *x3, double *x4);
int CalIntersect(double *a, double *b, double *c, double *normalvector, double *intersect);
#endif // HEADFILE_H_INCLUDED

