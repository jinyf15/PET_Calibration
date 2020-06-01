#include"headFile.h"
#include"globalVariable.h"

void reconUpdate(void *index)
{
    unsigned int sn;
    unsigned int nImg;
    int iIter;
    int nIter;
    int iSubSet;


    char fileName[1000];
    FILE *fp1;

    nIter = PSeq[0][0].nIter;
    nImg = Img.NX * Img.NY * Img.NZ;

    printf("update niter %d nImg %d\n",nIter,nImg);


    for(iIter = 0; iIter < nIter; iIter++)
    {
        iSubSet = PSeq[iIter][0].iSubSet;


        printf("update iIter %d \n",iIter);


        while(PSeq[iIter][0].nFinished < PSeq[iIter][0].nLoad || ImgBuf.imgReady !=1)
        {
            //printf("waiting update iter %d\n",iIter);
            usleep(1000);

        }

        pthread_mutex_lock(&comuMutex); // update img; needs to lock it;

        for(sn = 0; sn < nImg; sn++ )
        {
            if(Sen[iSubSet].image[sn] > 1.e-12 & Img.souMask[sn] > 0)
            {
                Img.image[sn] *= ImgBuf.image[sn] / Sen[iSubSet].image[sn];
                ImgBuf.image[sn] = 0;
            }

        }

        Img.imgReady = 1;
        ImgBuf.imgReady = 0;

        pthread_mutex_unlock(&comuMutex);


        sprintf(fileName,"%s//results//img_iter%d_iSet%d",DATA_DIR,iIter,(int)index);
        printf("%s\n",fileName);
        fp1 = fopen(fileName,"wb");

        checkFile(fp1,fileName);
        fwrite(&( Img.image[0] ),sizeof(double), nImg, fp1);
        fclose(fp1);


    }
    printf("here update end\n");


    pthread_exit(NULL);


}
