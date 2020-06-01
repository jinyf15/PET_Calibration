void sprstx_float(float sa[], unsigned long ija[], float x[], float b[],
	unsigned long n)
{
	void nrerror(char error_text[]);
	long int i,j,k;

	if (ija[1] != n+2) nrerror("mismatched vector and matrix in sprstx");
	//for (i=1;i<=n;i++) if(sa[i]>0) b[i]=sa[i]*x[i];//b[i]=sa[i]*x[i];
	for (i=1;i<=n;i++){
		b[i] = sa[i]*x[i];
	}

	for (i=1;i<=n;i++) {
		for (k=ija[i]; k<=(long int)(ija[i+1]-1); ++k) {
			j = ija[k];
			//if (j!=14597) continue;
			/*if ((ija[k]>ija[ija[1]-1])||(ija[k]<1)) { 
				continue;
			}*/
			b[j] += sa[k]*x[i];
			//if(sa[k]>0) b[j] += sa[k]*x[i];
			//if (b[j]!=0) {printf("%d %d %d %e %f %f\n", i,j,k, sa[k], x[i], b[j]);getchar();}
			
//			if ((j==13477)&&(sa[k]*x[i]!=0)) {printf("%d %d %f %f %f %f\n", k, i, sa[k], x[i], sa[k]*x[i], b[j]);getchar();}
//			if (j==78) {printf("sa=%.12e\tx=%.12e\tb=%.12e\n", sa[k], x[i], b[j]);}
		}
	}
}
