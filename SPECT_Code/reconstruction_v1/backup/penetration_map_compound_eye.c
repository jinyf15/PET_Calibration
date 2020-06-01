#include "headFile.h"
#include "./lib/cblas.h"


void penetration_map_compound_eye(double *phX, double *phY, double *phZ,double *line_up,double *line_p, double *line_dp,int N11,double *dirt,int N1,int N2,int N3,int N4,int i_nph,int nph,double *map,double muPlatinum,int iPrint)
{
	FILE *fp1,*fp2;
	int i,j,k,m,ii,jj,kk,ll;
	char file1[1000];

	double *p,*up,*dp,*op,*oup,*odp;

	double p_ph[3],tmp[4],dp_ph[3],up_ph[3];
	double *M;
	double alpha,phi,sinAlpha,cosAlpha,sinPhi,cosPhi,xx,yy,*local_orig,*local_dir,*line_orig,*line_dir;
	double *tmp2;
	double *TriVexU1,*TriVexD1,*TriVexU2,*TriVexD2,*TriVexU3,*TriVexD3;
	double *TriNormU,*TriNormD;
	double x1,x2,y1,y2,z1,z2;
	double t,tmp4[3],tmp5[3];
	double v0[3],v1[3],v2[3],v3[3],*InterSP,triangle[9],*InterSP1;
	int *Inside,*Inside1,*Index;
	double n[3],t1,t2;
	int size_1;
	double tmp0;
	//float *line_cross;
	time_t time1;

	p  = (double*) malloc((3*5)*sizeof(double));
	up = (double*) malloc((3*5)*sizeof(double));
	dp = (double*) malloc((3*5)*sizeof(double));

	M=(double*)malloc(9*sizeof(double));

	local_orig = (double*) malloc(3*N1*N2*N3*N4*sizeof(double));
	local_dir = (double*) malloc(3*N1*N2*N3*N4*sizeof(double));
	line_orig = (double*) malloc(3*N1*N2*N3*N4*sizeof(double));
	line_dir = (double*) malloc(3*N1*N2*N3*N4*sizeof(double));
	tmp2 = (double*) malloc(3*N1*N2*N3*N4*sizeof(double));

	TriVexU1 = (double*) malloc(3*N11*4*2*sizeof(double));
	TriVexD1 = (double*) malloc(3*N11*4*2*sizeof(double));
	TriVexU2 = (double*) malloc(3*N11*4*2*sizeof(double));
	TriVexD2 = (double*) malloc(3*N11*4*2*sizeof(double));
	TriVexU3 = (double*) malloc(3*N11*4*2*sizeof(double));
	TriVexD3 = (double*) malloc(3*N11*4*2*sizeof(double));
	TriNormU = (double*) malloc(3*N11*4*2*sizeof(double));
	TriNormD = (double*) malloc(3*N11*4*2*sizeof(double));
	InterSP = (double*) malloc(3*N11*4*2*sizeof(double));
	Inside  = (int*) malloc(N11*4*2*sizeof(int));
	Inside1 = (int*) malloc(N11*4*2*sizeof(int));
	Index  = (int*) malloc(N11*4*2*sizeof(int));
	InterSP1 = (double*) malloc(3*N11*4*2*sizeof(double));


	// get the votex of profile
	cblas_dcopy(3,&line_p[0],1,p,1);
	cblas_dcopy(3,&line_p[3*(N11-1)],1,&p[3],1);
	cblas_dcopy(3,&line_p[3*(2*N11-2)],1,&p[6],1);
	cblas_dcopy(3,&line_p[3*(3*N11-3)],1,&p[9],1);
	cblas_dcopy(3,&line_p[0],1,&p[12],1);

    cblas_dcopy(3,&line_up[0],1,up,1);
	cblas_dcopy(3,&line_up[3*(N11-1)],1,&up[3],1);
	cblas_dcopy(3,&line_up[3*(2*N11-2)],1,&up[6],1);
	cblas_dcopy(3,&line_up[3*(3*N11-3)],1,&up[9],1);
	cblas_dcopy(3,&line_up[0],1,&up[12],1);

    cblas_dcopy(3,&line_dp[0],1,dp,1);
	cblas_dcopy(3,&line_dp[3*(N11-1)],1,&dp[3],1);
	cblas_dcopy(3,&line_dp[3*(2*N11-2)],1,&dp[6],1);
	cblas_dcopy(3,&line_dp[3*(3*N11-3)],1,&dp[9],1);
	cblas_dcopy(3,&line_dp[0],1,&dp[12],1);

	// get pinhole position in global geometry;
	tmp[0] = 1;
	tmp[1] = 1;
	tmp[2] = 1;
	tmp[3] = 1;

	p_ph[0] = cblas_ddot(4,&p[0],3,tmp,1)/4.0;
	p_ph[1] = cblas_ddot(4,&p[1],3,tmp,1)/4.0;
	p_ph[2] = cblas_ddot(4,&p[2],3,tmp,1)/4.0;

	dp_ph[0] = cblas_ddot(4,&dp[0],3,tmp,1)/4.0;
	dp_ph[1] = cblas_ddot(4,&dp[1],3,tmp,1)/4.0;
	dp_ph[2] = cblas_ddot(4,&dp[2],3,tmp,1)/4.0;

	up_ph[0] = cblas_ddot(4,&up[0],3,tmp,1)/4.0;
	up_ph[1] = cblas_ddot(4,&up[1],3,tmp,1)/4.0;
	up_ph[2] = cblas_ddot(4,&up[2],3,tmp,1)/4.0;

	// ray is intersection with pinhole plane and direction
	m = 0;
	for(ii = 0; ii < N1; ++ii)
	{
		xx = (ii-N1/2.0)*dirt[0];

		for(jj = 0; jj < N2; ++jj)
		{
			yy = (jj-N1/2.0)*dirt[1];
			for(kk = 0;kk < N3;++kk)
			{
				alpha = kk*dirt[2];
				sinAlpha = sin(alpha);
				cosAlpha = cos(alpha);

				for(ll = 0;ll < N4;++ll)
				{
					phi = ll*dirt[3];
					sinPhi = sin(phi);
					cosPhi = cos(phi);
					local_orig[3*m] = xx; // the intersection of ray with pinhole plane,in aperture coordinate
					local_orig[3*m+1] = yy;
					local_orig[3*m+2] = 0;
					local_dir[3*m] = sinAlpha*cosPhi; // the ray direction in apertuer coordinate,attention! the baseline of phi is changed to x axis!!!
					local_dir[3*m+1] = sinAlpha*sinPhi;
					local_dir[3*m+2] = cosAlpha;
					m++;
				}
			}
		}
	}

	// transfer the ray intersection and diretion to global geometry
	M[0] = phX[0];M[1] = phX[1];M[2] = phX[2]; // transform matrix from aperture coordinate to global coordinate;
	M[3] = phY[0];M[4] = phY[1];M[5] = phY[2]; //
	M[6] = phZ[0];M[7] = phZ[1];M[8] = phZ[2]; //

	m = 0;
	for(i = 0;i < N1*N2*N3*N4;++i)
	{
		line_orig[3*m] = p_ph[0];
		line_orig[3*m+1] = p_ph[1];
		line_orig[3*m+2] = p_ph[2];
		m++;
	}

	cblas_dgemm(CblasColMajor,CblasNoTrans,CblasNoTrans,3,N1*N2*N3*N4,3,1.,M,3,local_orig,3,1.,line_orig,3); // the intersection of ray in global geometry;
	cblas_dgemm(CblasColMajor,CblasNoTrans,CblasNoTrans,3,N1*N2*N3*N4,3,1.0,M,3,local_dir,3,0.,line_dir,3); // the direction of ray in global geometry;

	free(local_orig);free(local_dir);

	// define the triangle surface of pinhole profile
	m = 0;
	for(ii = 0;ii < (4*N11-4);++ii)
	{ // define the vertice of triagnel surface consistiting the pinhole profile
		cblas_dcopy(3,&line_p[3*m],1,&TriVexU1[3*m],1);
		cblas_dcopy(3,&line_up[3*m],1,&TriVexU2[3*m],1);
		cblas_dcopy(3,&line_up[3*(m+1)],1,&TriVexU3[3*m],1);

		vector_ge(&TriVexU1[3*m],&TriVexU2[3*m],v1); //x1=B-A;
		vector_ge(&TriVexU1[3*m],&TriVexU3[3*m],v2); //y1=C-A;
		cross(v1,v2,&TriNormU[3*m]); //A X B

		cblas_dcopy(3,&line_up[3*(m+1)],1,&TriVexU1[3*m+3*(4*N11-4)],1);
		cblas_dcopy(3,&line_p[3*m],1,&TriVexU2[3*m+3*(4*N11-4)],1);
		cblas_dcopy(3,&line_p[3*(m+1)],1,&TriVexU3[3*m+3*(4*N11-4)],1);

		vector_ge(&TriVexU1[3*m+3*((4*N11-4))],&TriVexU2[3*m+3*((4*N11-4))],v1); //v1=B-A;
		vector_ge(&TriVexU1[3*m+3*((4*N11-4))],&TriVexU3[3*m+3*((4*N11-4))],v2); //v2=C-A;
		cross(v1,v2,&TriNormU[3*m+3*((4*N11-4))]); // norm=v1xv2; cross product

		cblas_dcopy(3,&line_p[3*m],1,&TriVexD1[3*m],1);
		cblas_dcopy(3,&line_dp[3*m],1,&TriVexD2[3*m],1);
		cblas_dcopy(3,&line_dp[3*(m+1)],1,&TriVexD3[3*m],1);

		vector_ge(&TriVexD1[3*m],&TriVexD2[3*m],v1); //x1=B-A;
		vector_ge(&TriVexD1[3*m],&TriVexD3[3*m],v2); //y1=C-A;
		cross(v1,v2,&TriNormD[3*m]);

		cblas_dcopy(3,&line_dp[3*(m+1)],1,&TriVexD1[3*m+3*((4*N11-4))],1);
		cblas_dcopy(3,&line_p[3*m],1,&TriVexD2[3*m+3*((4*N11-4))],1);
		cblas_dcopy(3,&line_p[3*(m+1)],1,&TriVexD3[3*m+3*((4*N11-4))],1);

		vector_ge(&TriVexD1[3*m+3*((4*N11-4))],&TriVexD2[3*m+3*((4*N11-4))],v1); //x1=B-A;
		vector_ge(&TriVexD1[3*m+3*((4*N11-4))],&TriVexD3[3*m+3*((4*N11-4))],v2); //y1=C-A;
		cross(v1,v2,&TriNormD[3*m+3*((4*N11-4))]);
		m++;
	}

	time(&time1);
	if(iPrint) printf("%s", asctime(localtime(&time1)));

	for(ii = 0;ii < N1*N2*N3*N4;++ii)
	{
		if(ii%(N2*N3*N4)==0 & iPrint)
		{
            time(&time1);
            printf("%d %s",ii/(N2*N3*N4),asctime(localtime(&time1)));
        }
		map[ii] = 0;
		k = 0;
		vector_ge(&line_orig[3*ii],up_ph,v1);

		n[0] = 0;
		n[1] = 0;
		n[2] = 1;

		t = dot(n,v1)/dot(n,&(line_dir[3*ii]));

		cblas_dcopy(3,&line_orig[3*ii],1,v2,1);
		cblas_daxpy(3,t,&(line_dir[3*ii]),1,v2,1); // y=a*x+y; intersection of ray with upper pinhole face

		if(IsInPolygon(4*N11-4,line_up,v2) == 0)
		{
			Inside[k] = -2;
			cblas_dcopy(3,v2,1,&InterSP[3*k],1);
			++k;
		}
		for(jj=0;jj<2*(4*N11-4);++jj)
		{
			vector_ge(&line_orig[3*ii],&TriVexU1[3*jj],v1);
			cblas_dcopy(3,&TriNormU[3*jj],1,n,1);

			t1 = dot(n,v1);
			t2 = dot(n,&(line_dir[3*ii])); // if tmp0==0 then the direction of line is parallel to the facet;

			if(fabs(t2)<1.0e-12)
			{
				if(t1==0)
				{
					printf("warning the ray is in one of triangel plane consisted of the pinhole profile plane, we treat it no intersection with this plane");
					printf("ii %d jj %d \n t1,t2,line_orig, line_dir,TriVexU1,TriNormU\n",ii,jj);
					continue;
				}
				else continue;  // parallel and not intersection;
			}
			else
			{
				t = t1/t2;
				cblas_dcopy(3,&line_orig[3*ii],1,v2,1);
				cblas_daxpy(3,t,&(line_dir[3*ii]),1,v2,1); // y=a*x+y; intersection with upper cone

				cblas_dcopy(3,&TriVexU1[3*jj],1,triangle,1);
				cblas_dcopy(3,&TriVexU2[3*jj],1,&triangle[3],1);
				cblas_dcopy(3,&TriVexU3[3*jj],1,&triangle[6],1);

				if(IsInTriangle(triangle,v2)==1)
				{ // whether the point is inside the triangle;
					Inside[k]=jj;
					cblas_dcopy(3,v2,1,&InterSP[3*k],1);
					++k;
				}
			}
		}
		for(jj=0;jj<2*(4*N11-4);++jj)
		{
			vector_ge(&line_orig[3*ii],&TriVexD1[3*jj],v1);

			cblas_dcopy(3,&TriNormD[3*jj],1,n,1);

			t1 = (dot(n,v1));
			t2 = dot(n,&(line_dir[3*ii]));

			if(fabs(t2)<1.0e-12)
			{
				if(t1==0)
				{
					printf("warning the ray is in one of triangel plane consisted of the pinhole profile plane, we treat it no intersection with this plane");
					printf("ii %d jj %d \n t1,t2,line_orig, line_dir,TriVexU1,TriNormU\n",ii,jj);
					continue;
				}
				else continue;
			}

			t = t1/t2;
			cblas_dcopy(3,&line_orig[3*ii],1,v2,1);
			cblas_daxpy(3,t,&(line_dir[3*ii]),1,v2,1); // y=a*x+y; intersection of ray with pinhole lower cone
			cblas_dcopy(3,&TriVexD1[3*jj],1,triangle,1);
			cblas_dcopy(3,&TriVexD2[3*jj],1,&triangle[3],1);
			cblas_dcopy(3,&TriVexD3[3*jj],1,&triangle[6],1);

			if(IsInTriangle(triangle,&v2)==1)
			{//
				Inside[k]=jj+2*(4*N11-4);
				cblas_dcopy(3,v2,1,&InterSP[3*k],1);
				++k;
			}
		}

		vector_ge(&line_orig[3*ii],dp_ph,v1);

		n[0] = 0;
		n[1] = 0;
		n[2] = 1;

		t = dot(n,v1)/dot(n,&(line_dir[3*ii]));

		cblas_dcopy(3,&line_orig[3*ii],1,v2,1);
		cblas_daxpy(3,t,&(line_dir[3*ii]),1,v2,1); // y=a*x+y; inserction of ray with lower pinhole face;

		if(IsInPolygon(4*N11-4,line_dp,v2)==0)
		{
			Inside[k] = -1;
			cblas_dcopy(3,v2,1,&InterSP[3*k],1);
			++k;
		}
		// the interection point with the pinhole upper surface;

		for(ll=0;ll<k;ll++) Index[ll]=ll;

		bubbleSort(&InterSP[2],k,3,Index);

		cblas_dcopy(3*k,InterSP,1,InterSP1,1);
		cblas_dcopy(k,Inside,1,Inside1,1);

		for(ll=0;ll<k;ll++) if(Index[ll]!=ll)
		{
			cblas_dcopy(2,&InterSP1[3*Index[ll]],1,&InterSP[3*ll],1);
			Inside[ll]=Inside1[Index[ll]];
		}

		for(i=0;i<k-1;++i)
		{
			if(Inside[i]>-3)
			{
				// intersect with up or lower surface;
				j=i+1;
				vector_ge(&InterSP[3*i],&InterSP[3*j],v3); // the vector between two inter section points

				if(Inside[i]==-2)
				{   // the inter section with aperture upper surface;
					map[ii]+=sqrt(dot(v3,v3))*muPlatinum;
				}
				else if(Inside[i]>=0&&Inside[i]<2*(4*N11-4))
				{
					cblas_dcopy(3,v3,1,v0,1);
					cblas_dscal(3,0.1,v0,1);
					cblas_daxpy(3,1,&InterSP[3*i],1,v0,1); //v0=(InterSP[j]-InterSP[i])*0.1+InterSP[i];
					kk=Inside[i];
					cblas_dcopy(3,&TriNormU[3*kk],1,n,1);
				    vector_ge(p_ph,&TriVexU1[3*kk],v1);
				    vector_ge(v0,&TriVexU1[3*kk],v2);

					if(dot(v2,n)*dot(v1,n)<0)
					{
						map[ii]+=sqrt(dot(v3,v3))*muPlatinum;
					}

				}
				else if(Inside[i]>=2*(4*N11-4))
				{
					cblas_dcopy(3,v3,1,v0,1);
					cblas_dscal(3,0.1,v0,1);
					cblas_daxpy(3,1,&InterSP[3*i],1,v0,1); //v0=(InterSP[j]-InterSP[i])*0.1+InterSP[i];

					kk = Inside[i]-2*(4*N11-4);

					cblas_dcopy(3,&TriNormD[3*kk],1,n,1);
					vector_ge(p_ph,&TriVexD1[3*kk],v1);
					vector_ge(v0,&TriVexD1[3*kk],v2);

					if(dot(v2,n)*dot(v1,n)<0)
					{
						map[ii]+=sqrt(dot(v3,v3))*muPlatinum;
					}

				}
				Inside[i]=-3;
			}

		}
	}
	time(&time1);
	if(iPrint) printf("%s", asctime(localtime(&time1)));



	free(p);free(up);free(dp);//free(op);free(oup);free(odp);
	free(M);
	free(line_orig);free(line_dir);
	free(tmp2);
	free(TriVexU1);free(TriVexD1);free(TriVexU2);free(TriVexD2);free(TriVexU3);free(TriVexD3);free(TriNormD);free(TriNormU);
	free(InterSP);free(Inside);free(Index);free(InterSP1);free(Inside1);

}

double dot(double v1[3],double v2[3]){
	double value;
	value=v1[0]*v2[0]+v1[1]*v2[1]+v1[2]*v2[2];
	return(value);
}
void cross(double v1[3],double v2[3],double v3[3]){
	v3[0]=v1[1]*v2[2]-v1[2]*v2[1];
	v3[1]=v1[2]*v2[0]-v1[0]*v2[2];
	v3[2]=v1[0]*v2[1]-v1[1]*v2[0];
}

void vector_ge(double v1[3],double v2[3],double v3[3]){
	v3[0]=v2[0]-v1[0];v3[1]=v2[1]-v1[1];v3[2]=v2[2]-v1[2];
}

int IsInTriangle(double*triangle,double*point){
	double v0[3],v1[3],v2[3];
	double dot00,dot01,dot02,dot11,dot12;
	double invDenom,v,u;

	// Compute vectors
	vector_ge(&triangle[0],&triangle[3],v0); //v0=C-A;
	vector_ge(&triangle[0],&triangle[6],v1); //V1=B-A;
	vector_ge(&triangle[0],point,v2);        //V2=P-A;

	// Compute dot products
	dot00 = dot(v0, v0);
	dot01 = dot(v0, v1);
	dot02 = dot(v0, v2);
	dot11 = dot(v1, v1);
	dot12 = dot(v1, v2);

	// Compute barycentric coordinates
	invDenom = 1.0 /(dot00 * dot11 - dot01 * dot01);
	u = (dot11 * dot02 - dot01 * dot12) * invDenom;
	v = (dot00 * dot12 - dot01 * dot02) * invDenom;

	// Check if point is in triangle
	return((u >=0) && (v >=0) && (u + v <=1.0));
}

int IsInPolygon(int NumVertex,double*polygon,double*point){ // this code is not finished....
	int ii,flag;
	double Triangle[9];
	Triangle[0]=polygon[0];Triangle[1]=polygon[1];Triangle[2]=polygon[2];
	flag=0;
	for(ii=0;ii<NumVertex-2;++ii){
		Triangle[3]=polygon[3*(ii+1)];Triangle[4]=polygon[3*(ii+1)+1];Triangle[5]=polygon[3*(ii+1)+2];
		Triangle[6]=polygon[3*(ii+2)];Triangle[7]=polygon[3*(ii+2)+1];Triangle[8]=polygon[3*(ii+2)+2];
		if(IsInTriangle(Triangle,point)==1){
			flag=1;
			break;
		}
	}
	return(flag);

}

void print_float(int Num,double *v){
	int i;
	for(i=0;i<Num;++i) printf("i%d,%15.14e\n",i,v[i]);
	getchar();
}

void print_int(int Num,int *v){
	int i;
	for(i=0;i<Num;++i) printf("i%d,%d\n",i,v[i]);
	getchar();
}
//
void bubbleSort(double arr[], int Count,int InN,int *Index){
	int i = Count,j;
	double temp;
	int tmp1;
	while(i > 0)
	{
		for(j = 0; j < i - 1; j++)
		{
			if(arr[j*InN] > arr[(j +1)*InN])
			{   temp = arr[j*InN];
			arr[j*InN] = arr[(j +1)*InN];
			arr[(j+1)*InN] = temp;
			tmp1=Index[j];
			Index[j]=Index[j+1];
			Index[j+1]=tmp1;
			}
		}
		i--;
	}

}

//int LineInTrianglePlane(double *Triangle,double *Norm,double *Line_orig,double *Line_dir,double *InterSP1,double *InterSP2){
//	double AxisX[3],AxisY[3];
//	double v0[1],v1[3],v2[3];
//	double x1,x2,x3,x4,y1,y2,y3,y4,tmp;
//	int ii,jj;
//
//	double M[9],Triangle_local[12],Triangle_global[12];
//	tmp=dot(Norm,Norm);
//	cblas_dscal(3,sqrt(tmp),Norm,1,Norm,1); // normolization of Norm:AxisZ;
//	vector_ge(Triangle,&Triangle[3],AxisX);
//	tmp=dot(AxisX,AxisX);
//	cblas_dscal(3,sqrt(tmp),AxisX,1,AxisX,1);// normolization of AxisX;
//	cross(Norm,AxisX,AxisY);
//	cblas_dcopy(3,AxisX,1,M,3);// tranformation matrix;
//	cblas_dcopy(3,AxisY,1,&M[3],3);
//	cblas_dcopy(3,Norm,1,&M[6],3);
//    vector_ge(Triangle,Line_orig,v0); // Triangle[0~2] the points for original of local geometry;
//	cblas_dgemm(CblasColMajor,CblasNoTrans,CblasNoTrans,3,1,3,1.,M,3,v0,3,1.,v1,3); // the line in local geometry;
//	cblas_dgemm(CblasColMajor,CblasNoTrans,CblasNoTrans,3,3,3,1.,M,3,Line_dir,3,1.,v2,3); // the line in local geometry;
//	x1=v1[0];y1=v1[1];
//	x2=v1[0]+v2[0];y2=v1[1]+v2[1];
//	cblas_dgemm(CblasColMajor,CblasNoTrans,CblasNoTrans,3,3,3,1.,M,3,Triangle,3,1.,Triangle_local,3);
//	vector_ge(Triangle_local,&Triangle_local[0],&Triangle_local[0]); // origin is set to be vex1;
//	vector_ge(Triangle_local,&Triangle_local[3],&Triangle_local[3]);
//	vector_ge(Triangle_local,&Triangle_local[6],&Triangle_local[6]);
//	cblas_dcopy(3,Triangle_local,1,&Triangle_local[9],1);
//	cblas_dcopy(9,Triangle,1,Triangle_global,1);
//	cblas_dcopy(3,Triangle,1,&Triangle_global[9],1);
//	for(ii=0;ii<3;ii++){
//       x3=Triangle_local[3*ii];y3=Triangle_local[3*ii+1];
//	   x4=Triangle_local[3*(ii+1)];y4=Triangle_local[3*(ii+1)+1];
//	   if(((x1-x2)*(y3-y4)-(y1-y2)*(x3-x4))==0){
//		   if(((x3-x1)*(y4-y1)-(x4-x1)*(y3-y1))==0){// the two lines are overlapping
//		     cblas_dcopy(3,&Triangle_global[3*ii],1,InterSP1,1);
//			 cblas_dcopy(3,&Triangle_global[3*(ii+1)],1,InterSP2,1);
//			 break;
//		   }
//		   else continue;
//	   }
//
//
//	}
//
//
//
//
//}
