#include "headFile.h"

double dotproduct(double *a, double *b){
	return a[0]*b[0]+a[1]*b[1]+a[2]*b[2];
}
double *scalevector(double a, double *b){
	double *y;
	y = (double *)malloc(sizeof(double)*3);
	y[0] = a*b[0];
	y[1] = a*b[1];
	y[2] = a*b[2];
	return y;
}
double *minusvector(double *a, double *b){
	double *y;
	y = (double *)malloc(sizeof(double)*3);
	y[0] = a[0]-b[0];
	y[1] = a[1]-b[1];
	y[2] = a[2]-b[2];
	return y;
}
int checkinrectangle(double *x0, double *x1, double *x2, double *x3, double *x4){
	double delta2[3], delta3[3],delta4[3];
	int i;
	for (i = 0;i < 3;i++){
		delta2[i] = x2[i]-x1[i];
		delta3[i] = x3[i]-x1[i];
		delta4[i] = x4[i]-x1[i];
	}

	if (fabs(dotproduct(delta2,delta3))<0.000001){
		if (dotproduct(x1,delta2)<=dotproduct(x0,delta2)&&dotproduct(x2,delta2)>=dotproduct(x0,delta2)){
			if (dotproduct(x1,delta3)<=dotproduct(x0,delta3)&&dotproduct(x3,delta3)>=dotproduct(x0,delta3)){
				return 1;	
			}	
		}
		return 0;
	}
	else{
		if (fabs(dotproduct(delta2,delta4))<0.000001){
			if (dotproduct(x1,delta2)<=dotproduct(x0,delta2)&&dotproduct(x2,delta2)>=dotproduct(x0,delta2)){
				if (dotproduct(x1,delta4)<=dotproduct(x0,delta4)&&dotproduct(x4,delta4)>=dotproduct(x0,delta4)){
					return 1;		
				}
			}
			return 0;
		}
		else{
			if (dotproduct(x1,delta4)<=dotproduct(x0,delta4)&&dotproduct(x4,delta4)>=dotproduct(x0,delta4)){
				if (dotproduct(x1,delta3)<=dotproduct(x0,delta3)&&dotproduct(x3,delta3)>=dotproduct(x0,delta3)){
					return 1;		
				}
			}
			return 0;
		}
	}

}

int CalIntersect(double *a,double *b, double *c,double *normalvector, double *intersect){
	double dirV[3],delta[3];
	int i;
	for (i = 0;i < 3;i++)
		dirV[i] = a[i]-b[i];
	double t;
	t = dotproduct(normalvector,dirV);
	if (t==0){
		return 0;
	}
	for (i = 0;i < 3;i++)
		delta[i] = c[i]-a[i];
	t = dotproduct(delta,normalvector)/t;
	intersect[0] = a[0] + dirV[0]*t;
	intersect[1] = a[1] + dirV[1]*t;
	intersect[2] = a[2] + dirV[2]*t;
	return 1;
}

