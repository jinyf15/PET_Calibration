#include"headFile.h"

void loadDetPixMaskForRecon(struct projection *proj,int projIndex)
{
    int NX;
    int NY;
    int iDet;
    int iPix;
    int iSou;
    int iSubDet;
    int buf;
    char fileName[1000];
    FILE *fp1;

    NX = proj[projIndex].NX;
    NY = proj[projIndex].NY;

    iDet = proj[projIndex].iDet;
    iSou = proj[projIndex].iPosi;


    proj[projIndex].detPixMask = (int *)malloc(sizeof(int)*NX*NY);

    sprintf(fileName,"%s//detPixMask_srfs//det_pix_mask%d.txt",DATA_DIR,iDet);
    fp1=fopen(fileName,"r");
    checkFile(fp1,fileName);



    for(iPix=0;iPix<NX*NY;iPix++)
    {

        fscanf(fp1,"%d\n",&(proj[projIndex].detPixMask[iPix]));
    }

    fclose(fp1);

    sprintf(fileName,"%s//det_pix_mask_iDet%d_iPos%d.txt",proj[projIndex].projFolder,iDet,iSou);
    fp1=fopen(fileName,"rt");
    checkFile(fp1,fileName);


    for(iPix=0;iPix<NX*NY;iPix++)
    {

        fscanf(fp1,"%d\n",&buf);

        if(buf<0)
        {
            proj[projIndex].detPixMask[iPix] = -1000;
        }

    }

    fclose(fp1);



    // save copy of detpixelMask;
//
//    sprintf(fileName,"%s//detPixMask_srfs/det_pix_mask%d.txt",DATA_DIR,iDet);
//    fp1=fopen(fileName,"w");
//
//    if(fp1==NULL) // can not open the file;
//    {
//        sprintf("error in writing file %s\n",fileName);
//        getchar();
//        getchar();
//        exit(-1);
//    }
//
//    for(iPix=0;iPix<NX*NY;iPix++)
//    {
//        fprintf(fp1,"%d\n",(det[iDet].detPixMask[iPix]));
//    }
//
//    fclose(fp1);



}
