#include "headFile.h"

struct detector *detInit(){
	struct detector *det;
	FILE *fid;
	int i;
	int NDet = 16;
	char filename[100];
	sprintf(filename,"%s/DetDef",DATA_DIR);
	printf("loading data from%s\n",filename);
	fflush(stdout);
	fid = fopen(filename,"rb");
	checkFile(fid,filename);
	det = (struct detector *)malloc(sizeof(struct detector)*NDet);
	int NX = 32;
	int NY = 32;
	int NZ = 20;
	printf("loading detector geometry\n");
	fflush(stdout);
	for (i = 0;i < NDet;i++){
		det[i].Mtr = (double *)malloc(sizeof(double)*9);
		fread((void *)det[i].Mtr, sizeof(double),9,fid);
	}
	for (i = 0;i < NDet;i++){
		det[i].NormalVector = (double *)malloc(sizeof(double)*18);
		fread((void *)det[i].NormalVector, sizeof(double),18,fid);
	}
	for (i = 0;i < NDet;i++){
		det[i].VertexPos = (double *)malloc(sizeof(double)*72);
		fread((void *)det[i].VertexPos, sizeof(double),72,fid);
	}
	for (i = 0;i < NDet;i++){
		int total = NX*NZ*2+NY*NZ*2+NX*NY;
		det[i].PixelPos = (double *)malloc(sizeof(double)*total*3);
		fread((void *)det[i].PixelPos, sizeof(double),total*3,fid);
	}
	return det;
	fclose(fid);
}
