#include "headFile.h"

struct PETpanel *petInit(){
	struct PETpanel *pet;
	pet = (struct PETpanel *)malloc(sizeof(struct PETpanel)*4);
	
	int i;
	for (i = 0;i < 4;i++){
		pet[i].NX = 64;
		pet[i].NY = 64;
		pet[i].NZ = 20;
		pet[i].DDX = 2*20.9/pet[i].NX;
		pet[i].DDY = 2*20.9/pet[i].NY;
		pet[i].DDZ = 10.0/pet[i].NZ;
		pet[i].MU = 0.058;
	}
	return pet;
}
