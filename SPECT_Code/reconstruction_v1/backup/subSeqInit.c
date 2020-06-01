#include"headFile.h"

struct parellSequence **subSeqInit(int nDet,int nSouP,int nSubDet,int nIter,int nThread)
{
    int ii;
    int jj;
    int nLoad;
    int iIter;
    int iDet;
    int iSouP;
    int iSubDet;
    int iThread;
    int *nSequPerThread;
    int nSequ;
    struct parellSequence **pSeq;


    printf("intialized the Parallel Sequence \nnDet %d nSouce Postion %d nSbuDet %d nIter %d nThread %d  \n",nDet,nSouP,nSubDet,nIter,nThread);
    nSequPerThread = (int *)malloc(sizeof(int)*nThread);
    pSeq = (struct parellSequence**)malloc(sizeof(struct parellSequence * )*nThread);

    if( ((nDet*nSouP*nSubDet*nIter)%nThread)!=0 )
    {
        printf("warning  in parallel sequence arrang;  each thread work load is not the same; it will increase calculation time; please check your sequence file\n");
        getchar();
    }

    nLoad = nDet*nSouP*nSubDet*nIter;

    for(iThread=0;iThread<nThread;iThread++)
    {
        pSeq[iThread] = (struct parellSequence *)malloc(sizeof(struct parellSequence)*nLoad);

        for(ii=0;ii<nLoad;ii++)
        {
            //pSeq[iThread][ii]=malloc(sizeof(struct parellSequence));
        }
    }

    iThread=0;
    ii=0;
    jj=0;
    nSequ=0;

    for(iIter=0;iIter<nIter;iIter++)
    {
        for(iSouP=0;iSouP<nSouP;iSouP++)
        {
            for(iDet=0;iDet<nDet;iDet++)
            {
                for(iSubDet=0;iSubDet<nSubDet;iSubDet++)
                {
                    pSeq[iThread][ii].nIter   = nIter;
                    pSeq[iThread][ii].iIter = iIter;

                    pSeq[iThread][ii].nSouP   = nSouP;
                    pSeq[iThread][ii].iSouP = iSouP;

                    pSeq[iThread][ii].nDet    = nDet;
                    pSeq[iThread][ii].iDet  = iDet;

                    pSeq[iThread][ii].nSubDet = nSubDet;
                    pSeq[iThread][ii].iSubDet = iSubDet;


                    pSeq[iThread][ii].nThread = nThread;
                    pSeq[iThread][ii].iThread = iThread;


                    pSeq[iThread][ii].totalSeq = nIter*nSouP*nDet*nSubDet;
                    pSeq[iThread][ii].iSeq  = jj;

                    printf("iThread %d iDet %d iIter %d iSeq %d iSouP %d iThread %d nDet %d nIter %d nSoup %d nSubDet %d iSubDet %d nTotalSeq %d\n",iThread,pSeq[iThread][ii].iDet,pSeq[iThread][ii].iIter,pSeq[iThread][ii].iSeq,pSeq[iThread][ii].iSouP,pSeq[iThread][ii].iThread,pSeq[iThread][ii].nDet,pSeq[iThread][ii].nIter,pSeq[iThread][ii].nSouP,pSeq[iThread][ii].nSubDet,pSeq[iThread][ii].iSubDet,pSeq[iThread][ii].totalSeq);

                    nSequPerThread[iThread] = ii+1;
                    //printf("ithread %d,subSequence %d \n",iThread,nSequPerThread[iThread]);

                    iThread++;
                    jj++;


                    if(iThread==nThread)
                    {
                        iThread = 0;
                        ii++;
                    }

                }
            }
        }
    }


    for(iThread=0;iThread<nThread;iThread++)
    {
        for(ii=0;ii<nSequPerThread[iThread];ii++)
        {
            pSeq[iThread][ii].nLoad = nSequPerThread[iThread]; // # of sequence belong to iThread;
           // printf("iThread %d %d \n",iThread,pSeq[iThread][ii].iDet);
            //printf("iThread %d iDet %d iIter %d iSeq %d iSouP %d iThread %d  nDet %d nIter %d nSoup %d nSubDet %d iSubDet %d nTotalSeq %d\n",iThread,pSeq[iThread][ii].iDet,pSeq[iThread][ii].iIter,pSeq[iThread][ii].iSeq,pSeq[iThread][ii].iSouP,pSeq[iThread][ii].iThread,pSeq[iThread][ii].nDet,pSeq[iThread][ii].nIter,pSeq[iThread][ii].nSouP,pSeq[iThread][ii].nSubDet,pSeq[iThread][ii].iSubDet,pSeq[iThread][ii].totalSeq);
        }

    }
    //iThread =0; ii=0;
    //printf("iThread %d iDet %d iIter %d iSeq %d iSouP %d iThread %d nDet %d nIter %d nSoup %d nSubDet %d iSubDet %d nTotalSeq %d\n",iThread,pSeq[iThread][ii].iDet,pSeq[iThread][ii].iIter,pSeq[iThread][ii].iSeq,pSeq[iThread][ii].iSouP,pSeq[iThread][ii].iThread,pSeq[iThread][ii].nDet,pSeq[iThread][ii].nIter,pSeq[iThread][ii].nSouP,pSeq[iThread][ii].nSubDet,pSeq[iThread][ii].iSubDet,pSeq[iThread][ii].totalSeq);

return(pSeq);
}

