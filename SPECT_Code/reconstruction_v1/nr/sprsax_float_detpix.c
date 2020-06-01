void sprsax_float_detpix(float sa[], unsigned long ija[], float x[], float b[],
	unsigned long n, unsigned long pix[], unsigned long npix)
{
	void nrerror(char error_text[]);
	unsigned long i,k, j, ii;

	if (ija[1] != n+2) nrerror("sprsax: mismatched vector and matrix");
	
	for (i=1; i<=npix; i++) {
		j=pix[i];
		b[j]=sa[j]*x[j];
	}

	ii=0;
 	for (i=1;i<=npix;i++) {
		j=pix[i];
		for (k=ija[j];k<=ija[j+1]-1;k++){
			if ((ija[k]>n)||(ija[k]<1)) { continue;}
			b[j] += sa[k]*x[ija[k]];
			ii++;
		}
	}
	printf("ii=%d\n", ii);getchar();
}
