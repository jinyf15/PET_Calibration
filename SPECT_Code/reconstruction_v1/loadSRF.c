#include"headFile.h"

struct sysSRF *loadSRF(struct projection *proj)
{
	char fileName[1000];
	char fileNameIndex[1000];
	char fileNameValue[1000];
	FILE *fp1;
	unsigned int  mSize;
	unsigned int nRead;
	unsigned long int totalSize;
	int iSou;
	int iDet;
	int iSubDet;
	int counter;
	struct sysSRF *SRF;
	int nSou;
	int nDet;
	int nSubdet;
	int iSRF;
	time_t time1;


	nSou = proj[0].nPosi;
	nDet = proj[0].nDet;
	nSubdet = proj[0].nSubDet;

	SRF = (struct sysSRF*) malloc( sizeof(struct sysSRF) * nSou * nDet * nSubdet);

	time(&time1); printf("%s", asctime(localtime(&time1)));

	iSRF = 0;
	totalSize = 0;
	for(iSou = 0; iSou < nSou; iSou++)
	{
        for(iDet = 0; iDet < nDet; iDet++)
        {
            for(iSubDet = 0; iSubDet < nSubdet; iSubDet++)
            {
                SRF[iSRF].iDet = iDet;
                SRF[iSRF].iSou = iSou;
                SRF[iSRF].iSubDet = iSubDet;
                SRF[iSRF].iSRF = iSRF;
                SRF[iSRF].nSRF = nDet * nSou * nSubdet;

                sprintf(fileName,"%s//srf//ijat_iSou%d_iDet%d_iSubDet%d",DATA_DIR,iSou,iDet,iSubDet);
                mSize = srfread1(fileName);

                totalSize += 2*mSize; // total element of SRF;
                if(PC_MAX_MEMORY <= totalSize*sizeof(float) / (1.e9) )
                {
                    printf("error: the SRF is larger than system available memory,the program will exit");
                    getchar();
                    getchar();
                    exit(-1);
                }

                SRF[iSRF].mSize = mSize;
                SRF[iSRF].nMax = mSize;

                SRF[iSRF].index = (unsigned int *) malloc( sizeof(unsigned int) *(mSize+1));
                SRF[iSRF].value = (float *) malloc( sizeof(float) *(mSize+1));


                sprintf(fileName,"%s//srf//ijat_iSou%d_iDet%d_iSubDet%d",DATA_DIR,iSou,iDet,iSubDet);
                fp1 = fopen(fileName,"rb");
                checkFile(fp1,fileName);

                nRead = fread( &( SRF[iSRF].index[1] ), sizeof(unsigned int),mSize,fp1);

                if(nRead !=mSize)
                {
                    printf("error in reading file %s nRead %ud, request %ud",fileName,nRead,mSize);
                    getchar();
                    exit(-1);
                }

                fclose(fp1);


                sprintf(fileName,"%s//srf//sat_iSou%d_iDet%d_iSubDet%d",DATA_DIR,iSou,iDet,iSubDet);
                fp1 = fopen(fileName,"rb");
                checkFile(fp1,fileName);

                nRead = fread( &( SRF[iSRF].value[1] ), sizeof(float),mSize,fp1);

                if(nRead !=mSize)
                {
                    printf("error in reading file %s nRead %ud, request %ud",fileName,nRead,mSize);
                    getchar();
                    exit(-1);
                }

                fclose(fp1);


            iSRF++;
            }
        }
	}

	time(&time1); printf("%s", asctime(localtime(&time1)));

	printf("loading of SRFs done, press enter key to continue\n");
	getchar();

	return(SRF);

}
