#include"headFile.h"

void loadDetPixMask(struct detector *det,int nDet,int iDet)
{
    int NX;
    int NY;
    int iPix;
    int iSubDet;
    double buf;
    char fileName[1000];
    FILE *fp1;

    NX = det[iDet].NX;
    NY = det[iDet].NY;


    det[iDet].detPixMask = (int *)malloc(sizeof(int)*NX*NY);

    sprintf(fileName,"%s//det_pix_mask%d.txt",DATA_DIR,iDet);
    fp1=fopen(fileName,"r");
    checkFile(fp1,fileName);


    iSubDet = 0;

    for(iPix=0;iPix<NX*NY;iPix++)
    {
        if(iSubDet>=det[iDet].nSubDet)
        {
            iSubDet = 0;
        }

        fscanf(fp1,"%le",&buf);

        if(buf>0.5)
        {
            det[iDet].detPixMask[iPix] = iSubDet;
            iSubDet++;
        }
        else
        {
            det[iDet].detPixMask[iPix] = -1000;
        }
    }

    fclose(fp1);

    // save copy of detpixelMask;

    sprintf(fileName,"%s//detPixMask_srfs/det_pix_mask%d.txt",DATA_DIR,iDet);
    fp1=fopen(fileName,"w");

    if(fp1==NULL) // can not open the file;
    {
        sprintf("error in writing file %s\n",fileName);
        getchar();
        getchar();
        exit(-1);
    }

    for(iPix=0;iPix<NX*NY;iPix++)
    {
        fprintf(fp1,"%d\n",(det[iDet].detPixMask[iPix]));
    }

    fclose(fp1);



}
