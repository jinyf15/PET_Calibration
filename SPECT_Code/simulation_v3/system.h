#ifndef SYSTEM_H_INCLUDED
#define SYSTEM_H_INCLUDED

// maxium pc memeory could be used by the pc
#define PC_MAX_MEMORY 480 // be cautious to change this value; if using too large value, pc will be dead;

struct detector
{
    int nDet; // number of detector in system
    int iDet; // which detector it is
    int nSubDet; // number of sub detector pixel
    int nDivide; // devide the detector pixel to smaller one, do not mass with nSubDet;
    int NX,NY,NZ; // number of detector pixel in each demension
    double dx,dy,dz; // pixel size
    double detX[3],detY[3],detZ[3];// detector coordinate in global
    double *dp; // detector position
    double detCent[3]; // detector center position
    double attenCoef; // attenuation coefficient
    int nphPerDet;
    int *phIndex;
    double **phCenter; // pinhhole center position
    double **phX; // pinhole X
    double **phY; //pinhole Y
    double **phZ; // pinhole Z;
    int *detPixMask; // detector pixelMask
};

struct source
{
    int nPosi; // num of source positon in system
    int iPosi; // which detector it is
    double eulerAg[3]; // rotation angle of source
    double souCenter[3];// source center Position
    double dx,dy,dz; // voxel size;
    int NX, NY, NZ; // demenision size of object space
    double *sp;// source position
    int *souMask; // source mask;
    int *souMap; // for rebin;
    double *image; // the value image;
};

struct parellSequence
{
    int nIter;// number of reconstruction iteration,in simulation set it to 1
    int iIter;

    int nDet; // number of detector;
    int iDet; // which detector is going to be calcuated in this sequence

    int nSouP; // number of source position;
    int iSouP; // which source positon is going to be calcuated in the sequence;

    int nSubDet;// number of subdetector
    int iSubDet; // which subdet is going to be calcualted in the sequence;

    int totalSeq;// total number of sequence = nIter*nSouP*nDet*nSubDet;
    int iSeq; // which sequence it is;

    int nThread; // number of thread;
    int iThread; // which thread this squence belong;

    int nLoad; // # of sequence belong to ithread

};

struct pinMap
{
    int  nph;//total number of pinhole
	int  iph;//pinhole index
	double radius;//pinhole radius
	double Rlimit; // radius limite for penetration map generation;
	double openAngle;//pinhole openangle ; this is half of total open area;
	double angleLimit;// open angle limit ;
	double incidenceAngle;//pinhole incidence angle
	double length;//pinhole length
	double channel;//pinhole center channel length
	double leadThickness;//lead thickness
	//double center;//pinhole center positioin
	//double *detAx;//pinhole x axis direction
	//double *detAy;//pinhole y axis direction
	//double *detAz;//pinhole z axis direction
	int NStep[4];//pinhole map 1~4st index

	double dirt[4];//pinhole map step in 1st index
	double *map;//pinhole map
};

struct sysSRF
{
    int nSRF; // total number of subSRFS;
    int iSRF; // which SRF it is belonged to
    int iDet; // which detector the response function belonged to
    int iSou; // which source postion the SRF belong;
    int iSubDet; // which sub detector the SRF belong to
    unsigned long int nMax; // size of memory could be size by vector value and index, nMax should be Larger than nSize;
    unsigned int mSize; // size of subSRF;
    float *value; // value of vector of subSFR;
    unsigned int *index; // index of vector for subSFR;
    int readyTran;
    int readySave;

};







#endif // SYSTEM_H_INCLUDED
