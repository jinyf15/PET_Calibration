#include"headFile.h"
#include"globalVariable.h"



int main()
{
    char fileName[1000];
    FILE *fp1;

    pthread_t *threadID;
    pthread_attr_t attr;
    int nThread;
    int iThread;
    int status;
    void *tmp;
    int kk;

    Proj = projInit();  // read the projection files

    imgInit(&Img,&ImgBuf);
    getchar();

    Srf = loadSRF(Proj);

    PSeq = seqInit_ForRecon(Proj);
    getchar();

    Sen = senInit_ForRecon(PSeq,Img);




    if(Proj[0].iProj) for(iThread = 0; iThread < 1; iThread++)
    {
        randSeed = -1;  // initiate the random see
        projGeneration((void *) iThread); // get meanprojection from source
    }



    for(kk = 0; kk<1; kk++)
    {
//        if(Proj[0].iProj&kk==0)
//
            genNoiseProj(Proj,kk); //add noise
//
//        else if(Proj[0].iProj&&kk==2)
//
//            addPulse2Proj();
//
//
//            //loadNoiseProj(Proj,kk);  // using to load noise projection from file



        threadHandle(kk); // handling process for reconstruction;

//        resetImag(&Img,&ImgBuf);
//        resetSeq(PSeq);
//
//        printf("kk %d\n",kk);
    }


    free(Proj);
    free(Img.image);
    free(Img.souMask);
    free(ImgBuf.image);
    free(ImgBuf.souMask);
    free(Srf);


}
