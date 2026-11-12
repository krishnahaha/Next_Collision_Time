/*Aim: Defines the structure of the hardSphere class needed for the simulation. It is the backbone of the DSMC software.
*/

#ifndef HARDSPHERE_H_INCLUDED
#define HARDSPHERE_H_INCLUDED

#include "particle.h"
#include "cellList.h"
#include<iostream>
#include<cstdlib>
#include<math.h>
#include<fstream>
#include <vector>
#include <unistd.h>
#include <time.h>
#include <iomanip>




struct parameter{
    int npart;   //precision = 1 for float, = 2 for double. 
    double packFrac;
    double temp;
    double lambdaMult;
    double maxTime;
    double dtFrac;
    double wallTemp, wallVel;
    int trial;
};
void readPara( FILE *fp, struct parameter *p ); //Function to read the parameter.



template < typename real >
class hardSphere{
public:
    hardSphere(struct parameter *p); //Constructor

    unsigned int npart;
    real packFrac, rho;
    real temp;
    real lambdaMult;

    particle<real>  *op; //An object of 'particle' class. How does the scope fo the particle class is visible here.
    cellList<real>  *oc; //An object of 'cellList' class
    real boxVolm, boxLen;
    real maxTime = 0.0, simTime = 0.0, prevTime = 0.0;
    unsigned int nColls;
    int nCellx, nCell;
    real dt, cellLen;
    real avgOmega1, avgNct, avgFFT, avgFPL, avgSpd, avgSpdSim; //Statistical quantity of interest.
    real vmax, dtFrac;
    real eps;
    real ke; //Kinetic energy.
    real brtot = 0.0, timeLast = 0.0;
    real wallTemp, wallVel;
    real saveInterval;
    int trial;


    
    
    void initPosition( int ); //=0 for random, = 1 for FCC. 
    void initVelocity( real meanVel, real temp ); 
    real gaussianRv( real mean, real std );
    real getdt();
    void displayPosVel( int num );
    real kineticEnergy();
    real avgSpeed();
    void calcStatisticalQuantities( real packFrac, real temp );
    void updatePositions( real dt );
    void sortIntoCells();
    inline void storeIntoCell( int ci, int i );
    inline int cellIndex( int cx, int cy, int cz );
    int collision(); //return the number of collisions executed by the subroutine.
    int numOfCollToAttempt( int npa );
    int pairCollision( int pi, int pj );
    real getWallCollTime( int pi );
    void boundaryConditionVelChange( int pi, real* dv );
    void sampleObservables();
    void sampleYvel();
    void samplePressure();
    
    void saveObservables();
    void saveYvel();
    void savePressure();
    void saveConfig();



    /*Observables*/
    std::vector<real> yVel;
    std::vector<real> simT;
    std::vector<real> pressArr;

    ~hardSphere(); //Destructor.
};



/* Read the parameter file from the commandline */
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

    //lambdaMult.
    fscanf(fp, "%s", string); fscanf(fp, "%s", string);
    p->lambdaMult = strtod(string, &endptr); //std::cout<<"\nlambdaMult = "<<p->lambdaMult;

    //maxTime
    fscanf(fp, "%s", string); fscanf(fp, "%s", string);
    p->maxTime = strtod(string, &endptr); //std::cout<<"\nmaxTime = "<<p->maxTime;

    //dtFrac
    fscanf(fp, "%s", string); fscanf(fp, "%s", string);
    p->dtFrac = strtod(string, &endptr); //std::cout<<"\nmaxTime = "<<p->maxTime;

    //wallTemp
    fscanf(fp, "%s", string); fscanf(fp, "%s", string);
    p->wallTemp = strtod(string, &endptr);

    //wallVel
    fscanf(fp, "%s", string); fscanf(fp, "%s", string);
    p->wallVel = strtod(string, &endptr);

    //trial
    fscanf(fp, "%s", string); fscanf(fp, "%s", string);
    p->trial = atoi(string); //std::cout<<"\n trial: "<<p->trial;
}
#endif //HARDSPHERE_H_INCLUDED