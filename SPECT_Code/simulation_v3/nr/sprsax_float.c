void sprsax_float(float sa[], unsigned long ija[], float x[], float b[],
	unsigned long n)
{
	void nrerror(char error_text[]);
	unsigned long i,k;

	if (ija[1] != n+2) nrerror("sprsax: mismatched vector and matrix");
	for (i=1;i<=n;i++) {
		b[i]=sa[i]*x[i];
		//if(sa[i]>0) b[i]=sa[i]*x[i];
		for (k=ija[i];k<=ija[i+1]-1;k++){
			if ((ija[k]>ija[ija[1]-1])||(ija[k]<1)) { continue;}
			b[i] += sa[k]*x[ija[k]];
			//if(sa[k]>0) b[i]+=sa[k]*x[ija[k]];
//			if (i==100) printf("sa=%.12e\tx=%.12e\tb=%.12e\n", sa[k], x[ija[k]], b[i]);
		}
	}
}
