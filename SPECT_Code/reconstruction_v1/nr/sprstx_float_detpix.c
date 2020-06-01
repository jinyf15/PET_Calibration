void sprstx_float_detpix(float sa[], unsigned long ija[], float x[], float b[],
	unsigned long n, unsigned long pix[], unsigned long npix)
{
	void nrerror(char error_text[]);
	long int i,j,k, l;

	if (ija[1] != n+2) nrerror("mismatched vector and matrix in sprstx");

	for (i=1; i<=npix; i++) {
		l=pix[i];
		b[l]=sa[l]*x[l];
	}

	for (i=1;i<=npix;i++) {
		l=pix[i];
		for (k=ija[l]; k<=(long int)(ija[l+1]-1); k++) {
			j=ija[k];
			if ((j>n)||(j<1)) { continue;}
			b[j] += sa[k]*x[l];
		}
	}
}
