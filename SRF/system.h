#ifndef SYSTEM_H_INCLUDED
#define SYSTEM_H_INCLUDED

struct detector{
	double *Mtr;
	double *NormalVector;
	double *VertexPos;
	double *PixelPos;
};

struct source{
	int NX;
	int NY;
	int NZ;
	double dx;
	double dy;
	double dz;
	double *Pos;
};

struct PETpanel{
	double MU;
	int NX;
	int NY;
	int NZ;
	double DDX;
	double DDY;
	double DDZ;
	double *dp;
};
#endif
