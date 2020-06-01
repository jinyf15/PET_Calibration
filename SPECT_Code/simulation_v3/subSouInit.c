#include"headFile.h"

struct source *subSouInit(int nPos,int NX, int NY, int NZ, double dx, double dy,double dz)
{
    char fileName[1000];
    FILE *fp1;
    int iPos;
    int flag;
    int maxlength, CharRetCd, arglen;
    char oneline[256];
    char *ptr2;
	char *ptr1;
	struct source *sou;


    maxlength =  256;
    sou = (struct source *)malloc(sizeof(struct source)*nPos); // allocate the memory for the source vector

    sprintf(fileName,"%s/setting.txt",DATA_DIR);
    fp1 = fopen(fileName,"r");
    checkFile(fp1,fileName);

    flag = 0;
    while (1)
    {
		CharRetCd = fgets (oneline, maxlength, fp1);
		//printf("oneline %s \n",oneline);
        ptr1 = strstr(oneline,"source parameter");
        //printf("ptr1 %s \n",ptr1);getchar();

        if(ptr1)
        {
            flag = 1;
            CharRetCd = fgets (oneline, maxlength, fp1);
            break;
        }
    }

    if(flag == 0)
    {
        printf("error could NOT find the source setting in %s\n",fileName);
        getchar();
        exit(-1);
    }


    for(iPos = 0; iPos < nPos; iPos++)
    {

    fscanf(fp1,"%lf\t %lf\t %lf\t %lf\t %lf\t %lf\n",&(sou[iPos].eulerAg[2]),&(sou[iPos].eulerAg[1]),&(sou[iPos].eulerAg[0]),&(sou[iPos].souCenter[0]),&(sou[iPos].souCenter[1]),&(sou[iPos].souCenter[2]));
    // attention rotation sequence is z->y->x; this order is different with input file-setting.txt, which is x,y,z;
    printf("nSouce %d souce %dth \n phi     beta     alpha    centerX  centerY  centerZ\n",nPos,iPos);
    printf(" %lf %lf %lf %lf %lf %lf\n",(sou[iPos].eulerAg[2]),(sou[iPos].eulerAg[1]),(sou[iPos].eulerAg[0]),(sou[iPos].souCenter[0]),(sou[iPos].souCenter[1]),(sou[iPos].souCenter[2]));

    calSouGlobalPos(sou,nPos,iPos,NX,NY,NZ,dx,dy,dz); // calculated souce global position;
    souMaskInit(sou,nPos,iPos);

    }

    printf("press anykey to continue\n");
    getchar();

    fclose(fp1);

    return(sou);






}
