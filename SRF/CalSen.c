#include "GlobalVariables.h"
#include "headFile.h"

double fabs(double a){
	return (a>0)?a:(0-a);
}
int deterSurf(int indexS, int indexD, int ExitSide[], int IncomeSide[]){
	double s0[3],s[3], x1[3],x2[3],x3[3],x4[3],normalvector[3],diff[3],dist;
	int i,ind,j,k;
	double tmp[3];
	ind = 0;
	for (i = 0;i < 3;i++){
		s[i] = Sou.Pos[indexS*3+i];
		x1[i] = Det[indexD].VertexPos[48+i];
		x2[i] = Det[indexD].VertexPos[48+3+i];
		x3[i] = Det[indexD].VertexPos[48+6+i];
		x4[i] = Det[indexD].VertexPos[48+9+i];
		normalvector[i] = Det[indexD].NormalVector[12+i];
		//printf("s%lf x1%lf x2%lf x3%lf x4%lf\n",s[i],x1[i],x2[i],x3[i],x4[i]);
	}
	//fflush(stdout);
	diff[0] = s[0]-x1[0];
	diff[1] = s[1]-x1[1];
	diff[2] = s[2]-x1[2];
	dist = dotproduct(diff, normalvector);
	double d2[3],d4[3];
	for (i = 0;i < 3;i++){
		tmp[i] = dist*normalvector[i];
		s0[i] = s[i]-tmp[i];
		d2[i] = x2[i]-x1[i];
		d4[i] = x4[i]-x1[i];
	}
	if (!(dotproduct(x1,d4)<=dotproduct(s0,d4)&&dotproduct(x4,d4)>=dotproduct(s0,d4))){
		IncomeSide[ind] = ((indexD%2)==0)?1:0;
		ind++;
	}
	if (!(dotproduct(x1,d2)<=dotproduct(s0,d2)&&dotproduct(x2,d2)>=dotproduct(s0,d2))){
		IncomeSide[ind] = (((indexD%4)>>1)==0)?3:2;
		ind++;
	}
	IncomeSide[ind] = 5;
	ind++;
	j = 0;k = 0;
	
	for (i = 0;i < 6;i++){
		if (IncomeSide[j]==i){			
			j++;
		}
		else{
			ExitSide[k] = i;
			k++;
		}
	}
	
	return (6-ind);
}

double CalSubSen(int indexS, int indexD, int ExitSide, int indexP, int *IncomeSide, int n){
	double s0[3],x1[3],dist1,hs1,s1;
	double diff1[3];
	double y1[3], block1[3], block2[3],l1, mu1,costheta1;
	double normalvector[3];
	int i,j;
	double pr1,pr2;
	for (i = 0;i < 3;i++){
		s0[i] = Sou.Pos[indexS*3+i];
		x1[i] = Det[indexD].PixelPos[indexP*3+i];
		normalvector[i] = Det[indexD].NormalVector[ExitSide*3+i];
	}
	mu1 = PET[indexD/4].MU;
	switch (ExitSide/2){
		case 0:
			s1 = PET[indexD/4].DDY*PET[indexD/4].DDZ;
			break;
		case 1:
			s1 = PET[indexD/4].DDX*PET[indexD/4].DDZ;
			break;
		case 2:
			s1 = PET[indexD/4].DDX*PET[indexD/4].DDY;
			break;
	}
	diff1[0] = s0[0]-x1[0];
	diff1[1] = s0[1]-x1[1];
	diff1[2] = s0[2]-x1[2];
	dist1 = fabs(dotproduct(diff1,diff1));
	hs1 = fabs(dotproduct(normalvector,diff1));
	costheta1 = hs1/sqrt(dist1);
	pr1 = s1*costheta1/4/PI/dist1;
	//printf("source bin:%d and pixel %d on detector %d\t",indexS,indexP,indexD);
	//fflush(stdout);
	for  (i = 0;i < n;i++){ //search point on the income side
		double V1[3], V2[3], V3[3], V4[3], normalvectorIn[3];
		for (j = 0;j < 3;j++){
			V1[j] = Det[indexD].VertexPos[IncomeSide[i]*12+j];
			V2[j] = Det[indexD].VertexPos[IncomeSide[i]*12+3+j];
			V3[j] = Det[indexD].VertexPos[IncomeSide[i]*12+6+j];
			V4[j] = Det[indexD].VertexPos[IncomeSide[i]*12+9+j];
			normalvectorIn[j] = Det[indexD].NormalVector[IncomeSide[i]*3+j];
		}
		if (CalIntersect(s0,x1,V1,normalvectorIn,y1)){
			if (checkinrectangle(y1, V1, V2, V3, V4)){
				//printf("find intersection point on income side:%d\n",IncomeSide[i]);
				//printf("press any key to continue\n");
				//fflush(stdout);
				//getchar();
				break;
			}
		}
	}
	//printf("i=%d n=%d\n",i,n);
	//fflush(stdout);	
	diff1[0] = y1[0]-x1[0];
	diff1[1] = y1[1]-x1[1];
	diff1[2] = y1[2]-x1[2];
	l1 = sqrt(dotproduct(diff1,diff1));
	pr1 = pr1 * (1-exp(-mu1*l1));
	if (IncomeSide[i]!=5){
		int detInd,sideInd,k;
		switch (IncomeSide[i]){
			case 0:
				if ((indexD%4) == 1){
					detInd = indexD/4*4 + 0;
				}
				else{
					detInd = indexD/4*4 + 2;
				}
				sideInd = 1;
				break;
			case 1:
				if ((indexD%4) == 0){
					detInd = indexD/4*4 + 1;
				}
				else{
					detInd = indexD/4*4 + 3;
				}
				sideInd = 0;
				break;
			case 2:
				if ((indexD%4) == 0){
					detInd = indexD/4*4+ 2;
				}
				else{
					detInd = indexD/4*4+ 3;
				}
				sideInd = 3;
				break;
			case 3:
				if ((indexD%4) == 2){
					detInd = indexD/4*4+ 0;
				}
				else{
					detInd = indexD/4*4+ 1;
				}
				sideInd = 2;
				break;
		}
		double V1[3], V2[3], V3[3], V4[3], normalvectorIn[3];
		for (j = 0;j < 3;j++){
			V1[j] = Det[detInd].VertexPos[sideInd*12+j];
			V2[j] = Det[detInd].VertexPos[sideInd*12+3+j];
			V3[j] = Det[detInd].VertexPos[sideInd*12+6+j];
			V4[j] = Det[detInd].VertexPos[sideInd*12+9+j];
			normalvectorIn[j] = Det[detInd].NormalVector[sideInd*3+j];
		}
		if (CalIntersect(s0,x1,V1,normalvectorIn,block1)){
			if (checkinrectangle(block1, V1, V2, V3, V4)){
				double V1[3], V2[3], V3[3], V4[3], normalvectorIn[3];
				for (j = 0;j < 3;j++){
					V1[j] = Det[detInd].VertexPos[60+j];
					V2[j] = Det[detInd].VertexPos[60+3+j];
					V3[j] = Det[detInd].VertexPos[60+6+j];
					V4[j] = Det[detInd].VertexPos[60+9+j];
					normalvectorIn[j] = Det[detInd].NormalVector[15+j];
				}
				if (CalIntersect(s0,x1,V1,normalvectorIn, block2)){
					if (checkinrectangle(block2,V1,V2,V3,V4)){
						double diff_block[3],l_block,mu2;
						mu2 = PET[detInd/4].MU;
						diff_block[0] = block2[0]-block1[0];
						diff_block[1] = block2[1]-block1[1];
						diff_block[2] = block2[2]-block1[2];
						l_block = sqrt(dotproduct(diff_block,diff_block));
						pr1 = pr1 * exp(-mu2*l_block);
					}

				}
			}
		}
	}
		// start to calculate pr2
	int indexD2 = indexD/4*4 + 8;
	int k;
	int flag = 0;
	double z1[3],z2[3];
	double l2 = 0.0;
	for (i = indexD/4*4+4;i < indexD2;i++){	
		int ExitSide2[6];
		int IncomeSide2[6];
		int numSide = 6-deterSurf(indexS, i, ExitSide2, IncomeSide2);
		for (k = 0;k < numSide;k++){
			double V1[3], V2[3], V3[3], V4[3], normalvectorIn[3];
			for (j = 0;j < 3;j++){
				V1[j] = Det[i].VertexPos[IncomeSide2[k]*12+j];
				V2[j] = Det[i].VertexPos[IncomeSide2[k]*12+3+j];
				V3[j] = Det[i].VertexPos[IncomeSide2[k]*12+6+j];
				V4[j] = Det[i].VertexPos[IncomeSide2[k]*12+9+j];
				normalvectorIn[j] = Det[i].NormalVector[IncomeSide2[k]*3+j];
			}
			if(CalIntersect(s0,x1,V1,normalvectorIn,z1)){
				if (checkinrectangle(z1, V1, V2, V3, V4))
					break;
			}
		}
		if (k != numSide){
			numSide = 6 - numSide;
			for (k = 0;k < numSide;k++){
				double V1[3], V2[3], V3[3], V4[3], normalvectorEx[3];
				for (j = 0;j < 3;j++){
					V1[j] = Det[i].VertexPos[ExitSide2[k]*12+j];
					V2[j] = Det[i].VertexPos[ExitSide2[k]*12+3+j];
					V3[j] = Det[i].VertexPos[ExitSide2[k]*12+6+j];
					V4[j] = Det[i].VertexPos[ExitSide2[k]*12+9+j];
					normalvectorEx[j] = Det[i].NormalVector[ExitSide2[k]*3+j];
				}
				if (CalIntersect(s0,x1,V1,normalvectorEx,z2)){
					if (checkinrectangle(z2, V1, V2, V3, V4))
						break;
				}
			}
			if (k!=numSide){
				double diff2[3];
				diff2[0] = z1[0] - z2[0];
				diff2[1] = z1[1] - z2[1];
				diff2[2] = z1[2] - z2[2];
				l2 = l2 + sqrt(dotproduct(diff2,diff2));
				flag++;;
			}
		}
		if (flag==2)
			break;
	}
	double mu2 = PET[(indexD2-1)/4].MU;
	pr2 = 1 - exp(-mu2*l2);
	return pr1*pr2;
}
void *CalSen(void *ThreadId){
	int totalS = Sou.NX*Sou.NY*Sou.NZ;
	int start, end;
	int i,j,k,l;
        int type;	// type = 3,4,5 is the number of exit sides
	int pid = (int)ThreadId;
	FILE *fp;
	char filename[100];
	sprintf(filename,"%s/srf/sen_sub%d",DATA_DIR,pid);
	fp = fopen(filename,"wb");
	fclose(fp);
	fp = fopen(filename,"ab");
	for (i = pid;i < totalS;i+=nThread){
		double sen = 0;
		j = 0;
		while (j < 12){
			if (j == 4){
				j = 8;
			}
			int ExitSide[6];
			int IncomeSide[6];
			type = deterSurf(i, j, ExitSide, IncomeSide);
			for (k = 0;k < type;k++){
				switch((int)(ExitSide[k]/2)){
					case 0:
						start = PET[j/4].NY*PET[j/4].NZ*ExitSide[k]/2;
						end = PET[j/4].NY*PET[j/4].NZ*(ExitSide[k]+1)/2;
						break;
					case 1:
						start = PET[j/4].NY/2*PET[j/4].NZ*2+PET[j/4].NX/2*PET[j/4].NZ*(ExitSide[k]-2);
						end = PET[j/4].NY/2*PET[j/4].NZ*2+PET[j/4].NX/2*PET[j/4].NZ*(ExitSide[k]-1);
						break;
					case 2:
						start = PET[j/4].NY/2*PET[j/4].NZ*2+PET[j/4].NX/2*PET[j/4].NZ*2+PET[j/4].NX/2*PET[j/4].NY/2*(ExitSide[k]-4);
						end = PET[j/4].NY/2*PET[j/4].NZ*2+PET[j/4].NX/2*PET[j/4].NZ*2+PET[j/4].NX/2*PET[j/4].NY/2*(ExitSide[k]-3);
						break;
				}
				for (l = start;l < end;l++){
					sen = sen + CalSubSen(i,j,ExitSide[k],l,IncomeSide,6-type);
				}
			}
			j++;
		}
		if (pid == 0 && i%10000==0){
			printf("Sen%d is %lf\n",i,sen);
			fflush(stdout);
		}
		fwrite(&sen,sizeof(double),1,fp);
	}	
	fclose(fp);
	return;
}
