#include "headFile.h"

long *initRandSeed(struct projection *Proj)
{
    int nPix;
    int nSou;
    int i;
    long *randSeed;

    nPix = Proj[0].NX * Proj[0].NY;
    nSou = Proj[0].nPosi;

    printf("nSou %d nPix %d\n",nPix,nSou);
    randSeed = (long *) malloc(sizeof(long) * nPix * nSou);

    for(i=0; i < nPix * nSou ; i++) randSeed[i] = -1; // negative value will initiate the random sequence;

    return(randSeed);



}
