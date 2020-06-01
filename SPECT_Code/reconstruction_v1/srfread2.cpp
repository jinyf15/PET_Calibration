#include "structureSettings.h"
#include "functionsList.h"
#include "constants.h"


unsigned long srfread2(float sat[], unsigned long ijat[], unsigned long msize, char filename1[], char filename2[]){
	FILE *fp1, *fp2;
    time_t time1;

	//for sparse storage
	unsigned long int tmp1, tmp2, i, nnz;

	//finding msize from the file ijat_na0

	fp2=openFile(filename2, "rb");
	if(fp2!=NULL){
		fread(&(tmp1), sizeof(unsigned long), 1, fp2);
		fclose(fp2);
	}

	fp2=openFile(filename2, "rb");
	for(i=1; i<=(tmp1-1); i++) {
		fread(&tmp2, sizeof(unsigned long), 1, fp2);
	}
	fclose(fp2);
	//msize=tmp2-1;
	//testing
	if(msize!=(tmp2-1)) { printf("<srfread2.c>: msize error !! %ld, %ld\n", msize, tmp2-1); getchar(); }

	//reading SRF at from the file sat_na0 and ijat_na0
	fp1=openFile(filename1, "rb"); printf("<srfread2.c>: reading srf from %s\n", filename1);
	fp2=openFile(filename2, "rb"); printf("<srfread2.c>: reading srf from %s\n", filename2);
	time(&time1); printf("%s", asctime(localtime(&time1)));
	nnz=0;
	for(i=1; i<=msize; i++) {
		if(fread(&(ijat[i]), sizeof(unsigned long), 1, fp2)!=1) {
			printf("<srfread2.c>: error in reading srf 1 !!\n");
			return nnz;
			//getchar();
		}
		if(fread(&(sat[i]), sizeof(float), 1, fp1)!=1) {
			printf("<srfread2.c>: error in reading srf 2 !!\n");
			//getchar();
			return nnz;
		}
		if(sat[i]<0) { printf("i, sat[i], ijat[i]: %ld %.6e %ld\n", i, sat[i], ijat[i]); getchar(); }
		nnz++;
	}
	time(&time1); printf("%s", asctime(localtime(&time1)));
	fclose(fp1); fclose(fp2);
	msize=ijat[ijat[1]-1]-1;
	printf("<srfread2.c>: response fucntion--numbers of none zeros: %ld, number of col %ld\n", msize, ijat[1]-2);
	return msize;
}


