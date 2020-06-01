#include"headFile.h"

struct pinMap *phInit(int nThread)
{
    char fileName[1000];
    FILE *fp1;
    double buffer;
    int iph;
    int nph;
    struct pinMap *ph;

    calculateMap(nThread); // calulating map;

    sprintf(fileName,"%s//phMap//map0",DATA_DIR);
    fp1=fopen(fileName,"rb");
    printf("%s\n",fileName);

    checkFile(fp1,fileName);

    fread(&buffer,sizeof(double),1,fp1);
    fclose(fp1);

    nph = (int) buffer;
    printf("numer of pinhole in penetration map %d",nph);getchar();
    ph = (struct pinMap *)malloc(sizeof(struct pinMap)*nph);

    for(iph=0;iph<nph;iph++)
    {
        subPhInit(ph,iph);
    }

    return(ph);
}
