#include "headFile.h"

struct source souInit(){
	struct source Sou;
	Sou.NX = 128;
	Sou.NY = 128;
	Sou.NZ = 128;
	Sou.dx = 0.25;
	Sou.dy = 0.25;
	Sou.dz = 0.25;
	int total = Sou.NX*Sou.NY*Sou.NZ;
	Sou.Pos = (double *)malloc(sizeof(double)*3*total);
	printf("NX:%d NY:%d NZ:%d dx:%lf dy:%lf dz:%lf\n",Sou.NX,Sou.NY,Sou.NZ,Sou.dx,Sou.dy,Sou.dz);
	fflush(stdout);
	int i,j,k,ind;
	for (k = 0;k < Sou.NZ;k++){
		for (j = 0;j < Sou.NY;j++){
			for (i = 0;i < Sou.NX;i++){
				ind = i+j*Sou.NX+k*Sou.NY*Sou.NX;
				Sou.Pos[3*ind] = Sou.dx * ((1-Sou.NX) / 2.0 + i);
				Sou.Pos[3*ind+1] = Sou.dy * ((1-Sou.NY) / 2.0 + j);
			        Sou.Pos[3*ind+2] = Sou.dz * ((1-Sou.NZ) / 2.0 + k);
			}
		}
	}
	//FILE *fp;
	//fp = fopen("../data/source","wb");
	//fwrite(Sou.Pos, sizeof(double),total*3,fp);
	//fclose(fp);
	return Sou;
}
