#include"headFile.h"

void checkFile(FILE *fp,char *fileName)
{

    if(fp==NULL) // can not open the file;
    {
        printf("error in reading file %s\n",fileName);
	fflush(stdout);
        getchar();
        getchar();
        exit(-1);
    }
}
