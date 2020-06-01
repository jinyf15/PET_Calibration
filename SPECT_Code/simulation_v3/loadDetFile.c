#include"headFile.h"

void loadDetFile(struct detector *det,int nDet,int iDet)
{
    char fileName1[1000];
    FILE *fp1;
    int ii;
    int nphPerDet;
    int inph;
    double *sysPar;
    double buffer;




    sprintf(fileName1, "%s//a%d.txt", DATA_DIR,iDet);
    fp1=fopen(fileName1, "r");

    printf("using %s as the system geometry\n", fileName1);
    printf("press anykey to continue\n");
    //getchar();

    checkFile(fp1,fileName1);

    if(fp1!=NULL) // load number of pinhole and pinhole index;
    {
        fscanf(fp1,"%le\n",&(buffer));
        nphPerDet=(int)(buffer);
        printf("number of pinhole belong to this detector %d\n",buffer,nphPerDet);

        det[iDet].phIndex = (int *)malloc(sizeof(int)*nphPerDet);
        det[iDet].nphPerDet = nphPerDet;

        for(inph=0;inph<nphPerDet;inph++) // load the pinhole index for detector
        {
            fscanf(fp1,"%le\n",&(buffer));
            det[iDet].phIndex[inph] = (int)(buffer);
            printf("detector %d phIndex %d\n",iDet,det[iDet].phIndex[inph]);
            getchar();
        }

    }

    sysPar = (double *)malloc(sizeof(double)*(10*nphPerDet+16));

    ii = 0;
    while (!(feof(fp1)))
    {


        fscanf(fp1, "%le\n", &(sysPar[ii]));
        printf("%d %e %d\n", ii,sysPar[ii], feof(fp1));
        ii++;
    }
getchar();
    //printf("total data read in a file: ii =%d nph=%d ii=%d\n", ii, nphPerDet, nphPerDet*10+16); //getchar();
    fclose(fp1);

    calDetGlobalPos(det,sysPar,iDet); // calcualting detector center position and detector coordinate;

    printf("iDet %d detX %le %le %le \n",iDet,det[iDet].detX[0],det[iDet].detX[1],det[iDet].detX[2]);
    printf("iDet %d detY %le %le %le \n",iDet,det[iDet].detY[0],det[iDet].detY[1],det[iDet].detY[2]);
    printf("iDet %d detZ %le %le %le \n",iDet,det[iDet].detZ[0],det[iDet].detZ[1],det[iDet].detZ[2]);
    printf("iDet %d det center %le %le %le \n",iDet,det[iDet].detCent[0],det[iDet].detCent[1],det[iDet].detCent[2]);

    printf("press anykey to continue\n");
    //getchar();

    det[iDet].phCenter = (double *)malloc(sizeof(double*)*nphPerDet); // allocate the memory for pinhole center belong to the detector;
    det[iDet].phX = (double **)malloc(sizeof(double*)*nphPerDet);
    det[iDet].phY = (double **)malloc(sizeof(double*)*nphPerDet);
    det[iDet].phZ = (double **)malloc(sizeof(double*)*nphPerDet);


    phGlobalPosi(det,iDet,sysPar);



    free(sysPar);


}
