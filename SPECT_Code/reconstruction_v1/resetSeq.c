#include "headFile.h"


void resetSeq(struct parellSequence **seq)
{

    int iIter;
    int ii;


    for(iIter=0;iIter<seq[0][0].nIter;iIter++)
    {
        for(ii=0;ii<seq[iIter][0].nLoad;ii++)
        {
            seq[iIter][ii].nFinished = 0;

        }
    }


}
