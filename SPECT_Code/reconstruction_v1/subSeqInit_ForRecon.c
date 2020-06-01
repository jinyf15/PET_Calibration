#include"headFile.h"

struct parellSequence **subSeqInit_ForRecon(int nDet,int nSouP,int nSubDet,int nIter,int nThread,int nSubSet_OSEM,int **osemPara)
{
    int ii;
    int jj;
    int *nLoad;
    int iIter;
    int iDet;
    int iSouP;
    int iSubDet;
    int iThread;
    int *nSequPerThread;
    int nSequ;
    struct parellSequence **pSeq;
    int nSRF;
    int iSRF;
    int nIterTrue;

    int nLoop,iLoop;
    int iSubSet;


    printf("intialized the Parallel Sequence \n nDet %d nSouce Position %d nSubDet %d nIter %d nThread %d  \n",nDet,nSouP,nSubDet,nIter,nThread);


    nSRF = nSouP * nDet * nSubDet;

    nIterTrue = fmax(nSubSet_OSEM,nIter/nSubSet_OSEM * nSubSet_OSEM);
    nLoop = nIterTrue/nSubSet_OSEM;

    if( ((nIter)%nSubSet_OSEM)!=0 )
    {
        printf("each subset of OSEM define as 1 iter, request iter %d,true iter %d",nIter,nIterTrue);
        getchar();
    }

    pSeq = (struct parellSequence**)malloc(sizeof(struct parellSequence * )*nIterTrue);


    for(iIter=0;iIter<nIterTrue;iIter++)
    {
        pSeq[iIter] = (struct parellSequence *)malloc(sizeof(struct parellSequence)*nSRF);
    }



    nLoad = (int *) malloc(sizeof(int)* nIterTrue);



    iThread = 0;
    iIter = 0;


    for(iLoop =  0; iLoop < nLoop; iLoop ++)
    {

        for(iSubSet = 0; iSubSet < nSubSet_OSEM; iSubSet++)
        {
            ii = 0;
            for(iSRF = 0; iSRF < nSRF; iSRF ++)
            {

                if(iSubSet == osemPara[iSRF][0]) // subset of osem match setting file;
                {

                    jj = osemPara[iSRF][1];
                    iSouP = osemPara[iSRF][2];
                    iDet = osemPara[iSRF][3];
                    iSubDet = osemPara[iSRF][4];

                    if(jj != ii)
                    {
                        printf("error in the OSEM setting file,in file jj %d in code ii %d\n",jj,ii);
                        getchar();
                        exit(-1);
                    }

                    pSeq[iIter][ii].nIter   = nIterTrue;
                    pSeq[iIter][ii].iIter = iIter;

                    pSeq[iIter][ii].nSouP   = nSouP;
                    pSeq[iIter][ii].iSouP = iSouP;

                    pSeq[iIter][ii].nDet    = nDet;
                    pSeq[iIter][ii].iDet  = iDet;

                    pSeq[iIter][ii].nSubDet = nSubDet;
                    pSeq[iIter][ii].iSubDet = iSubDet;


                    pSeq[iIter][ii].nThread = nThread;
                    pSeq[iIter][ii].iThread = iThread;

                    pSeq[iIter][ii].nSubSet = nSubSet_OSEM;
                    pSeq[iIter][ii].iSubSet = iSubSet;


                    pSeq[iIter][ii].totalSeq = -1;
                    pSeq[iIter][ii].iSeq  = -1;



                    ii++;
                    iThread++;
                    if(iThread==nThread) iThread = 0;

                }
            }

            nLoad[iIter] = ii; // numbe of SRF belong to coresponding iIter;
            iIter ++;
        }

    }



    for(iIter=0;iIter<nIterTrue;iIter++)
    {
        for(ii=0;ii<nLoad[iIter];ii++)
        {
            pSeq[iIter][ii].nLoad = nLoad[ii]; // # of srf belong to iIter;
            pSeq[iIter][ii].nFinished = 0;
            //printf("iIter %d %d \n",iIter,pSeq[iIter][ii].nLoad);
            //printf("iIter %d iDet %d iIter %d iSeq %d iSouP %d iThread %d  \n nDet %d nIter %d nSoup %d nSubDet %d iSubDet %d nTotalSeq %d\n",iIter,pSeq[iIter][ii].iDet,pSeq[iIter][ii].iIter,pSeq[iIter][ii].iSeq,pSeq[iIter][ii].iSouP,pSeq[iIter][ii].iThread,pSeq[iIter][ii].nDet,pSeq[iIter][ii].nIter,pSeq[iIter][ii].nSouP,pSeq[iIter][ii].nSubDet,pSeq[iIter][ii].iSubDet,pSeq[iIter][ii].totalSeq);
        }

    }

    //getchar();
    //iThread =0; ii=0;
    //printf("iThread %d iDet %d iIter %d iSeq %d iSouP %d iThread %d nDet %d nIter %d nSoup %d nSubDet %d iSubDet %d nTotalSeq %d\n",iThread,pSeq[iThread][ii].iDet,pSeq[iThread][ii].iIter,pSeq[iThread][ii].iSeq,pSeq[iThread][ii].iSouP,pSeq[iThread][ii].iThread,pSeq[iThread][ii].nDet,pSeq[iThread][ii].nIter,pSeq[iThread][ii].nSouP,pSeq[iThread][ii].nSubDet,pSeq[iThread][ii].iSubDet,pSeq[iThread][ii].totalSeq);

free(nLoad);

return(pSeq);
}

