#include "headFile.h"

unsigned int srfread1(char *filename){
	FILE *fp1;
	unsigned int msize, tmp1, tmp2, i;

	tmp1=0;
	tmp2=0;

	fp1 = fopen(filename,"rb");
	checkFile(fp1,filename);

	if(fp1!=NULL) {
		fread(&tmp1, sizeof(unsigned int), 1, fp1);
		fclose(fp1);

		fp1=fopen(filename, "rb");
		for(i=1; i<=(tmp1-1); i++) {
			fread(&tmp2, sizeof(unsigned int), 1, fp1);
		}
		fclose(fp1);

		msize=tmp2-1;
		printf("loading SRF %s\n numbers of non-zero elements: %u, number of columns %u\n", filename, msize, tmp1-2);
		return msize;
	}
	else
		return 0;
}



