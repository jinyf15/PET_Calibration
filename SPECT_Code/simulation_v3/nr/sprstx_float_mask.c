void sprstx_float_mask(float sa[], unsigned long ija[], float x[], float b[],
	unsigned long n, float mask[])
{
	void nrerror(char error_text[]);
	long int i,j,k;

	if (ija[1] != n+2) nrerror("mismatched vector and matrix in sprstx");
	for (i=1;i<=n;i++) b[i]=sa[i]*x[i]*mask[i];

	for (i=1;i<=n;i++) {
		for (k=ija[i]; k<=(long int)(ija[i+1]-1); k++) {
			j=ija[k];
			if ((ija[k]>ija[ija[1]-1])||(ija[k]<1)) { continue;}
			b[j] += sa[k]*x[i]*mask[i];
		}
	}
}
