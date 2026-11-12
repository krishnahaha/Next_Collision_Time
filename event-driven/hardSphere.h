/*Aim: Defines the structure of the hardSphere class needed for the simulation. It is the backbone of the DSMC software.
*/

#ifndef HARDSPHERE_H_INCLUDED
#define HARDSPHERE_H_INCLUDED

#include "particle.h"
#include "neighList.h"
#include<iostream>
#include<cstdlib>
#include<math.h>
#include<fstream>
#include <vector>
#include <unistd.h>
#include <time.h>
#include <iomanip>
#include <algorithm>

struct parameter{
    int npart;   //precision = 1 for float, = 2 for double. 
    double packFrac;
    double temp;
    double maxTime;
    double sampleObsInterval; //In fraction of maxTime.
    double eps; //Coeff of restitution
    int posInitType;
    double shellSize;
    std::string fileName;
    int useThermostat;
    int trial;
};

void readPara( FILE *fp, struct parameter *p ); //Function to read the parameter.

template < typename real >
class hardSphere{
public:
    hardSphere(struct parameter *p); //Constructor

    int npart;
    real packFrac;
    real temp;
    real maxTime;
    real sampleObsInterval;
    real eps;
    int posInitType;
    real shellSize;
    std::string fileName;
    int useThermostat;
    int trial;

    particle<real> *op; //An object of 'particle' class. How does the scope fo the particle class is visible here.
    particle<real> *opRef;

    particle<real> **cellList; //An pointer to create the cellList.
    particle<real> **eventLists;
    neighList<real> *on;
    particle<real> *root;
    
    real boxVolm, boxLen;
    real rho;
    real simTime = 0.0, cpuTime = 0.0;
    unsigned int nColls;
    int nCell, nCellx;
    real cellLen;
    real cellLenI;
    real thermostatInterval;

    real avgOmega, avgOmega1, avgNCT, avgFFT, avgFPL, avgSpd, avgSpdSim; //Statistical quantity of interest.
    real ke = 0.0;

    real viscosityT1[DIM] = {0.0}, viscosityT2[DIM] = {0.0}, viscosityAxx[DIM]={0.0};
    real thermalCondT1[DIM] = {0.0}, thermalCondT2[DIM] = {0.0}, thermalCondAxx[DIM] = {0.0};
    real thermalCondAxxPrev[DIM] = {0.0}, viscosityAxxPrev[DIM] = {0.0};

    real brTot = 0.0, brTotPrev =  0.0, timePrev = 0.0, prevCollTime = 0.0, prevViscTime = 0.0, prevThermTime = 0.0;
    real mopTransfer = 0.0, prevmopTransfer = 0.0, mopprevTime = 0.0;
    real mopUcontr = 0.0, prevmopUcontr = 0.0;
    real wallTemp = 1, wallVel = 1.0, vmax = 1.0;
    real eventsDeltaT = 0.0;
    real maxScheduleTime = 0.0;
    real refTime = 0.0;
    int nEventBoxes = 0;
    int currentList = 0;
    int numOfConfigForGr;
    unsigned int totalCollEvents = 0;
    unsigned int totalEvents = 0, totalSamplingEvents = 0, totalNeighEvents = 0, totalThermoEvents = 0;

    //The tags to save to operate.
    int saveFinalConfig = 0;
    int saveSystemReport = 1;
    int makeEfficiencyReport = 1;
    int makeNextCollTime = 0;



    void printVariables();
    void initPosition( int ); //=0 for random, = 1 for FCC. 
    void initVelocity( real meanVel, real temp ); 
    void initCorrelations();
    void copyConfigurations();
    void initParticleClassOtherVariables();
    real gaussianRv( real mean, real std );
    void displayPosVel( int num );
    void printPosVel( particle<real>* pi);
    real kineticEnergy();
    real avgSpeed();
    void calcStatisticalQuantities( real packFrac, real temp );
    void initCellList();
    void updatePosition( particle<real> *pi );
    void addToCellList( particle<real> *pi, int cellX, int cellY, int cellZ );
    void removeFromCellList( particle<real> *pi );
    real findNeighborlistUpdate( particle<real> *pi );
    void createNeighborList( particle<real> *pi );
    void initEvents();
    void addEvent( particle<real>* piTypeVariable );
    void addEventToBinaryTree( particle<real>* piTypeVariable );
    void removeEvent(  particle<real>* piTypeVariable );
    void createEvent( particle<real> *pi, particle<real> * pj, real time, int type );
    void addNextEventList();
    void checkOverlap();

    int findPairCollisionTime( particle<real> *pi, particle<real> *pj, real *tmin );
    void findMinCollTimeEvent( particle<real> *pi );
    void executeCollision( particle<real> *pi );
    void edmdSimStep();


    inline void storeIntoCell( int ci, int i );
    inline int cellIndex( int cx, int cy, int cz );
    int collision(); //return the number of collisions executed by the subroutine.
    int numOfCollToAttempt( int npa );
    int pairCollision( int pi, int pj );
    real getWallCollTime( int pi );
    void boundaryConditionVelChange( int pi, real* dv );
    void sampleObservables( particle<real> *we );
    void sampleYvel();
    real samplePressure();
    real samplePressureMop();
    void sampleMeanSqDistance();
    void sampleConfigForGr();
    void sampleCorrelations( particle<real> *pi, particle<real> *pj, real delVij, real eventTime );
    void sampleShearViscosity();
    void sampleDiagonalThermalConductivity();

    void saveObservables();
    void saveYvel();
    void savePressure();
    void saveConfig();
    void saveMeanSqDistance();
    void saveViscosity();
    void saveDiagonalThermalConductivity();
    void systemReport();



    // Observables.
    std::vector<real> yVel;
    std::vector<real> simT;
    std::vector<real> pressArr;
    std::vector<real> pressArrMop;
    std::vector<real> nextCollTime; int nctCount = 0;
    std::vector<real> msdArr;
    std::vector<real> viscosityArr;
    std::vector<real> thermalCondArr;
    std::vector<real> pressSplits;

    ~hardSphere(); //Destructor.
};


//Read the parameter file from the commandline.
void readPara ( FILE *fp, struct parameter *p){
    char string[300];
    char *endptr;

    //npart
    fscanf(fp, "%s", string); fscanf(fp, "%s", string);
    p->npart = atoi(string); //std::cout<<"\nnpart = "<<p->npart;

    //packFrac
    fscanf(fp, "%s", string); fscanf(fp, "%s", string);
    p->packFrac = strtod(string, &endptr); //std::cout<<"\npackFrac = "<<p->packFrac;

    //temp
    fscanf(fp, "%s", string); fscanf(fp, "%s", string);
    p->temp = strtod(string, &endptr); //std::cout<<"\ntemp = "<<p->temp;

    //maxTime
    fscanf(fp, "%s", string); fscanf(fp, "%s", string);
    p->maxTime = strtod(string, &endptr); //std::cout<<"\nmaxTime = "<<p->maxTime;

    //sampleObsInterval
    fscanf(fp, "%s", string); fscanf(fp, "%s", string);
    p->sampleObsInterval = strtod( string, &endptr);

    //eps
    fscanf(fp, "%s", string); fscanf(fp, "%s", string);
    p->eps = strtod(string, &endptr);  

    //posInit
    fscanf( fp, "%s", string); fscanf( fp, "%s", string);
    p->posInitType = atoi( string ); 

    //shellSize
    fscanf( fp, "%s", string); fscanf( fp, "%s", string);
    p->shellSize = strtod(string, &endptr);

    //fileName
    fscanf( fp, "%s", string); fscanf( fp, "%s", string);
    p->fileName = string;


    //useThermostat
    fscanf( fp, "%s", string); fscanf( fp, "%s", string);
    p->useThermostat = atoi( string ); 


    //trial
    fscanf(fp, "%s", string); fscanf(fp, "%s", string);
    p->trial = atoi(string); //std::cout<<"\n trial: "<<p->trial;
}
#endif //HARDSPHERE_H_INCLUDED