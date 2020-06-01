#include "headFile.h"

struct source *senInit_ForRecon(struct parellSequence **PSeq,struct source img)
{

    int nSubSet;
    int iSubSet;
    int iDet;
    int iSubDet;
    int iSouP;
    int nLoad;
    int iLoad;
    struct source *sen;
    unsigned long int sn;
    unsigned long int NIMG;
    double *buf;
    int nRead;
    char fileName[1000];
    FILE *fp;


     printf("loading sensitivity for each subset of osem\n");
     nSubSet = PSeq[0][0].nSubSet;


    sen = (struct source*) malloc(sizeof(struct source)*nSubSet);

    NIMG = img.NX*img.NY*img.NZ;

    buf = (double*) malloc(sizeof(double)*NIMG);



    for(iSubSet = 0; iSubSet < nSubSet; iSubSet++)
    {

       sen[iSubSet].NX = img.NX;
       sen[iSubSet].NY = img.NY;
       sen[iSubSet].NZ = img.NZ;

       sen[iSubSet].image = (double*)malloc(sizeof(double)*NIMG);

        for(sn = 0; sn< NIMG; sn++)
        {
            sen[iSubSet].image[sn] = 0.0;
        }

       nLoad = PSeq[iSubSet][0].nLoad;

       for(iLoad = 0; iLoad < nLoad; iLoad++)
       {
            iDet = PSeq[iSubSet][iLoad].iDet;
            iSouP = PSeq[iSubSet][iLoad].iSouP;
            iSubDet = PSeq[iSubSet][iLoad].iSubDet;
            printf("loading sen for osem %d nLoad %d iLoad %d iDet %d iSou %d iSubDet%d\n",iSubSet,nLoad,iLoad,iDet,iSouP,iSubDet);

            sprintf(fileName,"%s//srf//sen_iSou%d_iDet%d_iSubDet%d",DATA_DIR,iSouP,iDet,iSubDet);
            fp = fopen(fileName,"rb");
            checkFile(fp,fileName);

            nRead = fread(buf,sizeof(double), NIMG,fp);

            if(nRead !=NIMG)
            {
                printf("error in read file%s, nRead %d,request %ld",fileName,nRead,NIMG);

                getchar(); getchar();
                exit(-1);
            }

            fclose(fp);

            for(sn = 0; sn< NIMG; sn++)
            {
                sen[iSubSet].image[sn] += buf[sn];
            }

       }



       //printf("iSubSet %d nSubSet %d NX NY NZ %d %d %d NIMG %d\n",iSubSet,nSubSet,sen[iSubSet].NX,sen[iSubSet].NY,sen[iSubSet].NZ,NIMG);


    }

    printf("loading sensitivity for each subset of osem done\n press enter key to  continue\n");
    getchar();


    return(sen);
}
