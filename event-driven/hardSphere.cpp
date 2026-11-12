/*Aim: The definition file for hardSphere.h; Included definitions of the functions
*/
#include "hardSphere.h"

#include<iostream>
#include<fstream>
#include<mm_malloc.h>
#include<numeric>
#include <cstdlib> 

//#define M_PI 3.1415926535897932
#define EXTRAEVENTS 12


#define compileMode 1 //0: Debug mode; 1: Run; 2: Efficiency calculation

#define ObsSplit 1

//Class constructor
template<typename real>
hardSphere<real>::hardSphere( struct parameter *p ){
    //std::cout<<"\n Constructing hardSphere class variables";

    //Extrensic variables.
    npart             = p->npart;
    packFrac          = p->packFrac;
    temp              = p->temp;
    maxTime           = p->maxTime;
    sampleObsInterval = p->sampleObsInterval;
    eps               = p->eps;
    posInitType       = p->posInitType;
    shellSize         = p->shellSize;
    fileName          = p->fileName;
    useThermostat     = p->useThermostat;
    trial             = p->trial;
    
    //Intrinsic variables.
    unsigned long seed = time(NULL); //std::cout<<"\n seed "<<seed;
    srand48(seed); //Seed for the random number generator.
    //srand48(1); //For checking, the random number generator is initialized to 1. 
    
    boxVolm            = npart * M_PI / ( packFrac * 6);
    boxLen             = std::cbrt(boxVolm);
    simTime            = 0.0;
    eps                = 1.0; //Coefficient of restitution.
    rho                = 6 * packFrac / M_PI ; 
    numOfConfigForGr   = 50;
    

    //Initializing arrays to store data
    op = (particle<real>*)_mm_malloc( ( npart + EXTRAEVENTS ) * sizeof(particle<real>), 4096); //particle class array.
    opRef = (particle<real>*)_mm_malloc( ( npart + EXTRAEVENTS) * sizeof(particle<real>), 4096);
    on = (neighList<real>*)_mm_malloc( 20*npart * sizeof(real), 4096);
    //NeighList object. Check if sizeof(int) is coming out to be correct here.
    initPosition( posInitType );
    copyConfigurations(); //Needed for msd

    initVelocity( 0.0 , temp ); //Mean velocity = 0, and standard-deviation is of the temperature here.
    initParticleClassOtherVariables();
    calcStatisticalQuantities( packFrac, temp );
    //shellSize          = 1+avgFPL; //1 is here for the diameter of the Hard-Spheres.  
    //maxScheduleTime    = 5*avgFFT; //Modify it to see if putting it equal to 1 works.
    maxScheduleTime    = 1.0;
    thermostatInterval = 10*avgFFT;
    //eventsDeltaT       = 1 * avgNCT;
    eventsDeltaT   = 1.0/npart;
    ke = kineticEnergy();
    
    initEvents(); //Initialize the event list, root and other parameters.

    initCellList(); //Cell list.
    
    initCorrelations();

    particle<real> *pi;
    for( int i = 0; i < npart; i++ ){
        pi = op + i;
        //std::cout<<"\n Creating NeighborList pi: "<<pi->index;
        //std::cin.get();
        for( int k = 0; k < DIM; k++ ){
            pi->boxesTravelx[k] = 0;
            pi->rneig[k] = pi->pos(k);
        }
        pi->counter = 0;
        pi->tlast = 0.0;
        createNeighborList( pi );
    }
    //printVariables();
    //displayPosVel(10);
}



/*class destructor*/
template<typename real>
hardSphere<real>::~hardSphere(){
    //std::cout<<"\n A destructor variable called here";
    _mm_free(op);
    free(cellList);
    _mm_free(on);
}


/* Code the driver function */
template<typename real>
void hardSphere<real>::edmdSimStep(){
    //std::cout<<"\n\n ############################### edmdSimStep() begins ###############################";
    particle<real> *ei;
    ei = root->right;
    //std::cout<<"\n root->right->eventType: "<<ei->eventType<<"\t ei->index: "<<ei->index<<"\t ei->pos() "<<ei->pos(0)<<"\t"<<ei->pos(1)<<"\t "<<ei->pos(2) <<"\t simTime: "<<simTime;
    //std::cin.get();
    while( ei == NULL ){
        //std::cout<<"\n Adding new event list"; std::cin.get();
        addNextEventList(); //This connects the first event to the root event. Smart!!!
        ei = root->right; //Starting from the root. In the first step of adding events in the binary tree, a event is attached to the right of the root.
    }
    while( ei->left ){
        ei = ei->left;
    }

    simTime = ei->eventTime;
    //std::cout<<"\n The eventTime: "<<ei->eventTime<<"\t eventType: "<<ei->eventType<<"\t ei->index: "<<ei->index<<"\t ei->parent: "<<ei->parent->index;
    //std::cin.get();
    removeEvent( ei );
    switch( ei->eventType ){
        case 1:
            executeCollision( ei ); //Execute the collision, predict the new collisions for the particle, booking in the event calender.
            //std::cout<<"\n Executed collision";
            //std::cin.get();
            break;
        case 2:
            createNeighborList( ei );
            //std::cout<<"\n Creating neighbor list";
            break;
        case 3:
            sampleObservables( ei );
    }
    totalEvents++;
}


/*Calculate the mean Free Flight Time for the system of particles*/
template<typename real>
void hardSphere<real>::calcStatisticalQuantities( real packFrac, real temp ){
    real enskogCorr = (2 - packFrac) / ( 2 * pow( ( 1 - packFrac), 3) );
    //std::cout<<"\n enskogCorr: "<<enskogCorr;
    avgOmega1 = 24 * packFrac * enskogCorr * sqrt( temp / M_PI );
    avgOmega  = 12 * packFrac * npart * enskogCorr * sqrt( temp / M_PI );
    avgFFT    = 1.0 / avgOmega1; //std::cout<<"\n avgFFT: "<<avgFFT;
    avgNCT    = 1.0/avgOmega;
    avgFPL    = 1.0 / ( 6 * sqrt(2.0) * packFrac * enskogCorr );
    avgSpd    = sqrt ( 8 * temp / M_PI ); //std::cout<<"\n avgSpdTh: "<<avgSpd;
    avgSpdSim = avgSpeed(); //std::cout<<"\n avgSpdSim: "<<avgSpdSim;
    //double avgSpdSim   = kineticEnergy() / ( 1.5 * npart ); //std::cout<<"\n rmsSpdSim: "<<avgSpdSim;
}



/*initialize postitions randomly*/
template<typename real>
void hardSphere<real>::initPosition(int initType){
    if( initType == 0){
        //Random initialization.
        for( int i = 0; i < npart; i++){
            for( int k = 0; k < DIM; k++){
                op[i].pos(k) = boxLen * drand48();
            }
        }
    }
    else if( initType == 1 ){
        //Initialize on a FCC lattice.
        int subBox = round( cbrt( npart / 4 ) );
        if( 4* subBox * subBox * subBox != npart){
            std::cout<<"\n npart should be 4 * perfect cube: npart= "<<4*subBox*subBox*subBox;
            exit(3);
        }

        //std::cout<<"\n subBox: "<<subBox;
        real subBoxLen = boxLen / subBox;
        //std::cout<<"\n subBoxLen: "<<subBoxLen;
        int pi = 0;
        for( int a = 0; a < subBox; a++ ){
            for( int b = 0; b < subBox; b++ ){
                for( int c = 0; c < subBox; c++){
                    op[pi].pos(0) = (a + 0.25) * subBoxLen;
                    op[pi].pos(1) = (b + 0.25) * subBoxLen;
                    op[pi].pos(2) = (c + 0.25) * subBoxLen;
                    op[pi].index  = pi;  
                    pi++;

                    op[pi].pos(0) = (a + 0.75) * subBoxLen;
                    op[pi].pos(1) = (b + 0.75) * subBoxLen;
                    op[pi].pos(2) = (c + 0.25) * subBoxLen;
                    op[pi].index  = pi;
                    pi++;

                    op[pi].pos(0) = (a + 0.75) * subBoxLen;
                    op[pi].pos(1) = (b + 0.25) * subBoxLen;
                    op[pi].pos(2) = (c + 0.75) * subBoxLen;
                    op[pi].index  = pi;
                    pi++;

                    op[pi].pos(0) = (a + 0.25) * subBoxLen;
                    op[pi].pos(1) = (b + 0.75) * subBoxLen;
                    op[pi].pos(2) = (c + 0.75) * subBoxLen;
                    op[pi].index  = pi;
                    pi++;
                }//c ends.
            }//b ends.
        }//a ends.        
        //std::cout<<"\n No. of particles initialized positions: "<<pi;
        //std::cout<<"\n initPosition() done  op[pi].index: "<<op[pi-1].index;
    }//initType == 1/ FCC ends.
    else if( initType == 2 ){ //Read from file.
        int initIgnore = 2;
        std::ifstream fiN;
        std::string line;
        char pType; real radius;
        //std::cout<<fileName<<"\n";
        int count=0;

        fiN.open(fileName);
        if( fiN.fail() ){
            std::cerr<<"\n Error opeaning the file.\n\n"; exit(1);
        }
        std::setprecision(16);
        while( count < initIgnore ){
            getline( fiN, line);
            count++;
        }
        //fiN>>npart;
        //fiN>>boxLen; fiN>>boxLen; fiN>>boxLen; std::cout<<"\n boxLen: "<<boxLen;
        int nlines = npart;

        particle<real> *pi = &op[0];
        int j = 0;
        for( int i = 0; i < nlines; i++ ){
            fiN>>pType;
            for( int k = 0; k < DIM; k++ ){
                fiN>>pi->pos(k);

                if( pi->pos(k) >= boxLen || pi->pos(k) < 0.0){
                    pi->pos(k) -= boxLen * floor( pi->pos(k) / boxLen );
                }
            }
            //pi->setPbc( boxLen );
            pi->index = j; j++;
            pi++;
            fiN>>radius;
        }
        fiN.close();   
    }
    //saveConfig(); //Save the initial configuration
    //checkOverlap();
}



/*initialize velocities according to MB distribution*/
template<typename real>
void hardSphere<real>::initVelocity( real meanVel, real temp ){
    real vcm[3] = {0.0, 0.0, 0.0};
    for( int i = 0; i < npart; i++){
        for( int k = 0; k < DIM; k++){
            op[i].vel(k) = gaussianRv( meanVel, sqrt(temp) ); //te
            vcm[k] += op[i].vel(k);
        }
    }
    //Get the cm velocity by dividing it with npart;
    vcm[0] /= npart; 
    vcm[1] /= npart; 
    vcm[2] /= npart; 
    //Subtract from each particle velocity to get velocity to be equal to zero.
    real vcm2 = 0.0;
    for( int i = 0; i < npart; i++){
        for( int k = 0; k < DIM; k++){
            op[i].vel(k) -= vcm[k];
            vcm2 += op[i].vel(k) * op[i].vel(k);
        }
    }
    //Scale to match the temperature and the mean value.
    double localTemp = vcm2 / (npart*3.0);
    localTemp = (temp / localTemp);
    for( int i = 0; i < npart; i++){
        for( int k = 0; k < DIM; k++){
            op[i].vel(k) = meanVel + sqrt(localTemp) * op[i].vel(k);
        }
    }
    //std::cout<<"\n Kinetic energy/npart: "<<kineticEnergy()/npart;
    
}



/* Store the initial configurations for MSD calcualations */
template< typename real >
void hardSphere<real>::copyConfigurations(){
    for( int i = 0; i < npart; i++ ){
        //std::cout<<"\ni "<<i<<"\t";
        for( int k = 0; k < DIM; k++ ){
            opRef[i].pos(k) = op[i].pos(k); //std::cout<<"\t"<<opRef[i].pos(k);
            opRef[i].vel(k) = op[i].vel(k);
        }
    }
}



/*Generate gaussian random variable*/
template<typename real>
real hardSphere<real>::gaussianRv( real mean, real std ){
    static int out = 0;
    real v1, v2, r2;
    real y1; static real y2;
    double sqrtLnx1;
    if( out == 0){
        do{
            v1 = -1 + 2 * drand48();
            v2 = -1 + 2 * drand48();
            r2 = v1 * v1 + v2 * v2;
        }while( r2 >= 1.0 || r2 == 0);

        sqrtLnx1 = sqrt( -2 * log( r2 ) / r2 );

        y1 = mean + std * sqrtLnx1 * v1;

        y2 = mean + std * sqrtLnx1 * v2;
        out = 1;
        return y1;
    }
    else{
        out = 0; 
        return y2;
    }
}


/*Initialize other variables of the particle class*/
template<typename real>
void hardSphere<real>::initParticleClassOtherVariables(){
    particle<real> *pi = op; //Pointer to the first element of the op array of particles.
    for( int i = 0; i < npart; i++ ){
        for( int k = 0; k < DIM; k++ ){
            pi->rneig[k] = pi->pos(k);
            pi->boxesTravelx[k] = 0;
        }
        pi->tlast = 0.0;
        pi++;
    }
}



// Initialize thermalConductivity calculations.
template<typename real>
void hardSphere<real>::initCorrelations(){
    particle<real> *pi = op;
    real KE = 0.0;
    for( int i = 0; i < npart; i++ ){
        //Get the kinetic energy of particle i
        KE = 0.0;
        for( int k = 0 ; k < DIM; k++ ){
            KE += 0.5 * ( pi->vel(k) * pi->vel(k) );
        }
        
        //Multiply it with the velocity in different directions.
        for( int k = 0; k < DIM; k++ ){
            thermalCondT1[k] += ( KE * pi->vel(k) ) ; //Along one direction.
            viscosityT1[k]   += (pi->vel(k) * pi->vel( (k+1)%DIM )); //Shear direction.
        }

        //Go to the next particle.
        pi++;
    }
}


/*Display position and velocity of first 'num' number of particles on screen*/
template<typename real>
void hardSphere<real>::displayPosVel( int num){
    std::cout<<"\n #particle \t pos(x) \t pos (y) \t pos(z) \t vel (x) \t vel(y) \t vel(z)";
    std::cout<<"\n #boxLe: "<<boxLen;
    for( int i = 0; i < num; i++ ){
        std::cout<<"\n"<<i<<"\t";
        for( int k = 0; k < DIM; k++ ){
            std::cout<<op[i].pos(k)<<"\t";
        }
        for( int k = 0; k < DIM; k++ ){
            std::cout<<op[i].vel(k)<<"\t";
        }
    }
}



/* Display position-velocity of the specified particle*/
template<typename real>
void hardSphere<real>::printPosVel( particle<real>* pi){
    updatePosition(pi);

    std::cout<<"\nsimTime: "<<simTime<<"\tpi: "<<pi->index<<"\tPos: ";
    for( int k = 0; k < DIM; k++ ){
        std::cout<<pi->pos(k)<<"\t";
    }
    std::cout<<"Vel: ";
    for( int k = 0; k < DIM; k++ ){
        std::cout<<pi->vel(k)<<"\t";
    }
    std::cout<<"\t boxLen: "<<boxLen;
}




/*Calculate the total kinetic energy*/
template<typename real>
real hardSphere<real>::kineticEnergy(){
    ke = 0.0;
    for( int i = 0; i < npart; i++ ){
        for( int k = 0; k < DIM; k++ ){
            ke += op[i].vel(k) * op[i].vel(k); 
        }
    }
    ke *= 0.5; //ke is a hardSphere class variable.
    return ke;
}


/*Calculate the average speed*/
template<typename real>
real hardSphere<real>::avgSpeed(){
    real as = 0.0, s = 0.0;
    for( int i = 0; i < npart; i++ ){
        s = 0.0;
        for( int k = 0; k < DIM; k++ ){
            s += op[i].vel(k) * op[i].vel(k);
        }
        as += sqrt(s);
    }
    return ( as/(1.0 * npart) );
}



/*Initialize the cell list*/
template<typename real>
void hardSphere<real>::initCellList(){
    nCellx = (int)( boxLen - 0.0001 ) / shellSize;

    while( nCellx * nCellx * nCellx > 8 * npart ){
        nCellx--; 
    }

    //std::cout<<"\n shellSize: "<<shellSize;
    nCell = nCellx * nCellx * nCellx; //std::cout<<"\n nCell: "<<nCell<<"\t";
    cellList = (particle<real>**)calloc( nCell , sizeof(particle<real>*) );
    for( int i = 0; i < nCell; i++ ){
        cellList[i] = NULL;
    }
    cellLenI = nCellx / boxLen;

    //Add to cellList
    for( int i = 0; i < npart; i++ ){
        particle<real> *pi = op + i;
        //std::cout<<"\npi: "<<pi->index<<"\t"<<pi->pos(0)*cellLenI<<"\t "<<pi->pos(1)*cellLenI<<"\t "<<pi->pos(2)*cellLenI ;
        addToCellList( pi, pi->pos(0)*cellLenI, pi->pos(1)*cellLenI, pi->pos(2)*cellLenI );
    }
    //std::cout<<"\n Cell List initialized ";
    //std::cout<<"\n Members of cell "<<(nCell-1);
    //particle<real> *pi = cellList[nCell-1];
    //std::cout<<"\n pi->index: "<<pi->index;
    //while(pi->next){
    //    std::cout<<"\t"<<pi->index;
    //    pi = pi->next;
    //}
    
}


/*Get the cell index*/
template<typename real>
inline int hardSphere<real>::cellIndex( int cx, int cy, int cz ){
    //std::cout<<"\n cx: "<<cx<<"\t cy: "<<cy<<"\t "<<cz;
    return ( (int)(cx + cy*nCellx + cz * nCellx * nCellx) );
}



/**** addToCellList the particle op[i] ***
 * Add particle to the cell list at cell index cellX, cellY, cellZ
 * Each of the particle in the cell list is appended to a doubly linked list.
 * Why this is needed to be doubly linked list???? 
*/
template<typename real>
void hardSphere<real>::addToCellList( particle<real> *pi, int cellX, int cellY, int cellZ){
    pi->cell = cellIndex( cellX, cellY, cellZ ); //Get cellIndex.
    //std::cout<<"\n pi->cell: "<<pi->cell;
    //std::cin.get();
    pi->next  = cellList[pi->cell]; //The address of the first element at site p->cell is linked the the next address of pi.
    //std::cout<<"\n pi->cell: "<<pi->cell<<"\t";
    if( pi->next ){ 
        pi->next->prev = pi; //The reverse linking.
    }
    cellList[pi->cell] = pi; 
    pi->prev = NULL; //Reset the prev link of pi to NULL so that it is helpful in the if() condition.
    pi->nearBoxEdge = ( cellX == 0 || cellY == 0 || cellZ == 0 || cellX == nCellx - 1 || cellY == nCellx - 1 || cellZ == nCellx - 1 );
    //std::cout<<"\n pi: "<<pi->index<<" Added to cell list "<<pi->cell;
}


/* removeFromCellList() */
template<typename real>
void hardSphere<real>::removeFromCellList(particle<real> *pi){
    if(pi->prev)
        pi->prev->next = pi->next;  //Remove particle from cellList.
    else
        cellList[pi->cell] = pi->next;
    if(pi->next) 
        pi->next->prev = pi->prev; //Reverse linking.
        
        
    //std::cout<<"\n pi: "<<pi->index<<" removed from cell list";
}


// Update the position of a the particle.
template<typename real>
void hardSphere<real>::updatePosition( particle<real> *pi ){
    real dt = simTime - pi->tlast;
    //std::cout<<"\n dt: "<<dt;
    pi->tlast = simTime; //Update the local time to t.
    
    //For pressure using method of planes.
    real aftPosY, prevPosY = pi->pos(1);

    for( int k = 0; k < DIM; k++ ){
        pi->pos(k) += dt * pi->vel(k);
        
        if( k == 1 ){
            aftPosY = pi->pos(1);
        }

        if( pi->pos(k) >= boxLen ){
            pi->pos(k) -= boxLen;
            pi->boxesTravelx[k] += 1;
        }
        else if( pi->pos(k) < 0.0 ){
            pi->pos(k) += boxLen;
            pi->boxesTravelx[k] -= 1;
        }
    }

    if( prevPosY < boxLen*0.5 && aftPosY > boxLen*0.5 )
        mopTransfer += abs(pi->vel(1)); 
    else if( prevPosY > boxLen*0.5 && aftPosY < boxLen*0.5 )
        mopTransfer += abs(pi->vel(1)); 

    //pi->updatePos(dt); //Update the positions with velocity and by dt amount.
    //pi->setPbc( boxLen );
}


// checkoverlap between the pair of particles.
template<typename real>
void hardSphere<real>::checkOverlap(){
    particle<real> *pi, *pj;
    real r2 = 0.0, tmp = 0.0;
    for( int i = 1; i < npart; i++ ){
        pi = &op[i];
        updatePosition( pi );
        //pi->setPbc(boxLen);
        for( int j = 0; j < i; j++ ){
            pj = &op[j];
            r2 = 0.0;
            updatePosition( pj );
            //pj->setPbc( boxLen ); //This is not going to work here, as the particles have a fixed diameter, if the center does not crosses, then also there can be a overlap to the particle positions.
            for( int k = 0 ; k < DIM; k++ ){
                tmp = pi->pos(k) - pj->pos(k);
                r2 += tmp * tmp;
            }
            if( ( round(r2 * pow(10, 8) )/pow( 10, 8 ) ) < 0.99999990){
                //std::setprecision(16);
                std::cout<<std::setprecision(16)<<"\n r2: "<<r2;
                printPosVel(pi); printPosVel( pj );
                //std::cin.get();
            }
        }
    }
}




/***findNeighborlistUpdate**
 * Assumes pi is up to date in position. Note that the particle is always in the same box as its neighborlist position (p->xn)
 * Returns the time when the particle will cross the allowed value of the shellsize if moves with its velocity.
*/
template<typename real>
real hardSphere<real>::findNeighborlistUpdate( particle<real> *pi){
    real dr[DIM] = {0.0, 0.0, 0.0}, dv[DIM]={0.0, 0.0, 0.0}, b=0.0, dv2=0.0, dr2=0.0;
    //updatePosition( pi );
    for( int k = 0; k < DIM; k++ ){
        dr[k] = pi->pos(k) - pi->rneig[k];
        dr2 += dr[k] * dr[k];

        dv[k] = pi->vel(k);
        dv2 += dv[k] * dv[k];

        b += dr[k] * dv[k];
    }
    real md = (shellSize - 1) * 0.5; // Whenever the particle crosses the specified path length the neighborlist will be updated.
    real disc = b * b - dv2 * ( dr2 - md * md );
    if (disc < 0  || -b + sqrt(disc) < 0 ){
        return 0.0; //No real root. Means, the particle has already crossed the allowed shell thickness, immediately needs the neighbourlist update.
        std::cout<<"\n disc < 0 or -b+sqrt(disc)<=0 "<<disc;
    }
    real tupd1 = ( ( -b + sqrt( disc) )/dv2 );
    real tupd2 = ( ( -b - sqrt( disc) )/dv2 );
    if ( tupd1 > 0 && tupd2 > 0){
        return std::min( tupd1, tupd2 );
    }

    return tupd1;
}



/**** makeNeighborList for a given particle pi
 * Update the position of the particle, use PBC for its positions.
 * Remove the particle from its cell list as its position is updated and add it to a correct cell list site.
 * Set the current neighbor list position to the current position of the particle.
*/
template<typename real>
void hardSphere<real>::createNeighborList( particle<real> *pi ){
    particle<real> *pj;
    //std::cout<<"\n Making neighborlist of pi: "<<pi->index;
    int cellX[DIM];
    real dr[DIM];
    real r2, rm2 = shellSize * shellSize;
    //std::cout<<"\n Variables initialized";
    //std::cin.get();

    updatePosition( pi );
    //pi->setPbc( boxLen );
    removeFromCellList( pi );

    //Remove pi^th particle from the neighborlist of its old neighbors. This step can be slightly costly here to get into.
    //std::cout<<"\n Remove from ith particle from neighborlist of others";
    //std::cin.get();
    int bi, ai = pi->index;
    int nneigh = on[ai].neigh.size();
    //std::cout<<"\n ai: "<<ai;
    //std::cout<<"\t No. of neighbours: "<<nneigh;
    //Printing the index of all the neighbors.
    //for(  int ji = 0; ji < nneigh; ji++){
        //std::cout<<"\t bi: "<<on[ai].neigh[ji] ;
    //}

    for( int i = 0; i < nneigh; i++ ){
        bi = on[ai].neigh[i]; //std::cout<<"\n bi: "<<bi;
        //std::cout<<"\n No. of neighbours: "<<on[bi].neigh.size()<<"\t bi: "<<bi;
        auto& neighVector = on[bi].neigh;
        neighVector.erase( std::remove(neighVector.begin(), neighVector.end(), ai),  neighVector.end() );
    }

    for( int k = 0 ; k < DIM; k++ ){
        pi->rneig[k] = pi->pos(k); //Set the center of the neighbourlist to the particle's current position. This update is for the findNeighborlistUpdate() function to keep track of how much the particle has moved from the time when its neighbor list was updated there.
        cellX[k]     = pi->pos(k) * cellLenI;
        //std::cout<<"\t cellX[k]: "<<cellX[k];
    }

    addToCellList( pi, cellX[0], cellX[1], cellX[2] );
    cellX[0] += nCellx; cellX[1] += nCellx; cellX[2] += nCellx; //This step is useful if in case cellX[k] == 0, then in the following loop it will not create the problem.
    
    
    on[pi->index].neigh.clear(); //Clear the neighbourlist array.
    //std::cout<<"\n cellX[0]"<<cellX[0]<<"  "<<cellX[1]<<"  "<<cellX[2];
    //std::cin.get();
    int c1, c2, c3;
    for( int ci = cellX[0] - 1; ci < cellX[0] + 2; ci++ ){
        c1 = ci%nCellx;
        for( int cj = cellX[1] - 1; cj < cellX[1] + 2; cj++ ){
            c2 = cj%nCellx; c2 *= nCellx;
            for( int ck = cellX[2] - 1; ck < cellX[2] + 2; ck++ ){
                c3 = ck%nCellx; c3 *= (nCellx * nCellx); //Method to find the cell list index here.
                //std::cout<<"\n ci \% nCellx: "<<(ci % nCellx); std::cout<<"\t cj \% nCellx: "<<(cj % nCellx); std::cout<<"\t ck \% nCellx: "<<(ck % nCellx);

                //pj = cellList[cellIndex( ci % nCellx, cj % nCellx, ck % nCellx )];
                pj = cellList[ c1 + c2 + c3 ];
                //std::cout<<"\t pj->index: "<<pj->index;

                //Above will give the address of the first element sitting in neighbouring cell. Loop over all the elements in the cell to determine the neighbourlist. 
                int piInd = pi->index;
                while(pj){ //This demands that the last element of the cellList has to be NULL. 
                    if( pj->index != piInd ){
                        //std::cout<<"\n got pj: "<<pj->index <<"\t pi: "<<pi->index;
                        //std::cout<<"\n pkk: "<<pkk;
                        //std::cin.get();
                        
                        updatePosition( pj ); //It might cross the given cell because of the position update here. But this is not getting used here...
                        r2 = 0.0;
                        for( int k = 0; k < DIM; k++ ){  //dr[k] = pi->rneig[k] - pj->rneig[k];    //the rniegh[k] is not updated for the pj particle ?????
                            dr[k] = pi->rneig[k] - pj->rneig[k]; //THIS IS DIFFERENT FROM THE FRANK'S CODE.
                            //if(pi->nearBoxEdge){    //Only for the pi because its neighbourlist is getting calculated, if it not the last/first cell, its next neighbour will not be that far away. This is how the nearBoxEdge is determined.
                                if( dr[k] > 0.5*boxLen )    //If greater than half the box length here.
                                    dr[k] -= boxLen;
                                else if( dr[k] < -0.5*boxLen )
                                    dr[k] += boxLen;                        
                            //}//nearBoxEdge ends.
                            r2 += dr[k] * dr[k];
                        }//k ends.
                        if( r2 < rm2 ){
                            on[pi->index].neigh.push_back( pj->index );
                            bool isPresent = std::find(on[pj->index].neigh.begin(), on[pj->index].neigh.end(), pi->index) != on[pj->index].neigh.end();
                            if (!isPresent){
                                on[pj->index].neigh.push_back(pi->index);
                            }
                            //on[pj->index].neigh.push_back( pi->index );
                        }
                    }//pj!=pi ends.
                    pj = pj->next; //At one place this is getting a wierd value, not NULL, can not be addressed here as such.
                }//while pj != NULL ends.
            }//ck ends.
        }//cj ends. 
    }//ci ends.     
    //std::cout<<"\n no. of neighbours: "<<on[pi->index].neigh.size();
    findMinCollTimeEvent( pi ); //After neighborlist update, this is needed.
    //std::cout<<"\n neighborlist created ...";
    //std::cin.get();
    totalNeighEvents++;
}


/************* initEvents ******************
 * Determine the time for which a event is scheduled in the binary tree ( eventsDeltaT ).
 * Get the number of event buckets. Allocate the memory.
 * initialize the root event and the write event to a file. 
 * If thermostating is switched on, add it to the event tree.
*/
template<typename real>
void hardSphere<real>::initEvents(){
    //std::cout<<"\n eventsDeltaT: "<<eventsDeltaT;
    nEventBoxes = std::ceil( maxScheduleTime / eventsDeltaT );
    maxScheduleTime = nEventBoxes * eventsDeltaT;

    eventLists = (particle<real> **)calloc( nEventBoxes+1, sizeof(particle<real>*) );
    if( !eventLists ){
        std::cout<<"\n Memory allocation to eventLists failed"; exit(3);
    }

    root = op + npart; //Root event created. Pointing to the last+1 particle in system. This address is assigned to it.
    root->eventTime = -10.11; //The root will be empty, its childs will have events.
    root->eventType = 9; 
    root->index = npart+1;
    root->parent = NULL; //Its Brahma, the originator!!!
    root->right = NULL;
    root->left  = NULL;

    particle<real> *we = op + npart + 1; //Setting up the write event.
    we->eventTime = sampleObsInterval;
    we->eventType = 3;
    we->pj = NULL;
    addEvent( we );  //we is short of writeEvent.
    

    /*if( useThermostat ){
        particle<real> *thermostatEvent = op + npart + 2;
        thermostatEvent->eventTime = thermostatInterval;
        thermostatEvent->eventType = 3;
        thermostatEvent->pj = NULL;
        addEvent( thermostatEvent );
    }*/
}


/*Add a event to the event tree *
 * Check if the event time is in the current list, if so add it to the binary tree
 * If not, add it to a doubly linked list, which is for those events which occur far in the future.
*/
template<typename real>
void hardSphere<real>::addEvent( particle<real> *newEvent ){
    real dt = newEvent->eventTime - refTime;
    //std::cout<<"\n dt: "<<dt<<"\t refTime: "<<refTime<<"\t eventTime: "<<newEvent->eventTime<<"\t eventType: "<<newEvent->eventType;
    if( dt < eventsDeltaT ){
        newEvent->queue = currentList;
        addEventToBinaryTree( newEvent );
        //std::cout<<"\n Added event to binary tree";
    }
    else{
        //Put it in the appropriate event bucket.
        int listId = ( currentList + dt/eventsDeltaT );
        //std::cout<<"\t currentList: "<<currentList<<"\t dt: "<<dt<<"\t eventsDeltaT: "<<eventsDeltaT;
        if( listId >= nEventBoxes ){
            //Means, belongs to the overflow list.
            listId -= nEventBoxes;
            if( listId > currentList -1 ){ //Even after subtraction, it is away in time, make it in the overflow list. Current list gets reset in the addNextEventList() function.
                listId = nEventBoxes;
            }
        }
        //Doubly linked list is getting created.
        newEvent->queue = listId;
        //std::cout<<"\t listId: "<<listId;
        //std::cin.get();
        newEvent->right = eventLists[listId];
        newEvent->left  = NULL; //Making the newEvent the first elemment of the eventLists[][] array/pointer.
        if( newEvent->right ){
            newEvent->right->left = newEvent; //Doubly linking the list.
            //std::cout<<"\n doublely linked";
        }
        eventLists[listId] = newEvent;
        //std::cout<<"\t back linked";
    }
}


/***Add the event to the binary tree if its in the current list**
 * Find the location to add event into the B tree (loc)
    * Start from the root
    * Goto left, if event time is less than the current location event time otherwise right
 * 
*/
template<typename real>
void hardSphere<real>::addEventToBinaryTree( particle<real>* newEvent){
    real time = newEvent->eventTime; //std::cout<<"\t time: "<<time;
    particle<real> *loc = root; //Start from the root element.
    int occupied = 1;
    while( occupied ){ //Find the location to add event into the B tree (loc). Loops have a very strong shifting properties.
        //std::cout<<"\n loc->eventTime: "<<loc->eventTime<<"\tdt: "<<(loc->eventTime - refTime);
        //std::cin.get();
        if( time < loc->eventTime ){ //Then add it to the event tree by going to the left most leaf.
            if( loc->left){
                loc = loc->left; //In loop, go always to the left.
            }
            else{
                loc->left = newEvent; 
                occupied = 0;
            }
        }//time < loc->eventTime ends.
        else{ //Go to the right of the current location to find the possible fit.
            if( loc->right ){ //If not NULL.
                loc = loc->right;
                //std::cout<<"\n From addEventToBinaryTree, loc->eventTime: "<<loc->eventTime;
                //std::cin.get();
            }
            else{  //If the current right location is NULL, happens for the root event.
                loc->right = newEvent; //In the first runs, attaches the event to the root.
                occupied = 0; //Now the correct position found.
            }
        }//time > loc->eventTime ends.  
    }//while ends.

    newEvent->parent = loc; //What ever the location obtained, it will be the parent of the newEvent and will be at the leaf node. Hence, making its left and right nodes NULL.
    newEvent->left = NULL;
    newEvent->right = NULL;
}


/* Remove an event ( from the Binary tree if its going to happen in near future/from the doubly linked list, if its going to happen far in the future )*/
template<typename real>
void hardSphere<real>::removeEvent(particle<real> *oldEvent){
    //std::cout<<"\n removeEvent, time: "<<oldEvent->eventTime<<"\t type: "<<oldEvent->eventType;
    if (oldEvent->queue != currentList){ //The event is not in the Binary Tree. No sorting is required, just removal form the unsorted doubly linked list.
        if (oldEvent->right)
            oldEvent->right->left = oldEvent->left; //If it is in the right extreme, the next if will make the oldEvent->left->right = NULL.
        
        if (oldEvent->left){
            oldEvent->left->right = oldEvent->right; //simply pulling and attaching.
        }
        else{
            eventLists[oldEvent->queue] = oldEvent->right; //If left side is NULL, means that event is the first event of the queue, make it the first.
            //eventList[], then make its first successor the first event:)
        }
        return; //Deletion completed from the unsorted doubly linked list, go out of the function removeEvent. 
    }

    /* If the The event was in the current list / Binary Tree.
        3 Steps:
        a) The oldnode does not have a left element; just attach its right element to the parent of the oldnode.
        b) The old node have a right element, no left branch; just attach the left element to the parent of oldnode.
        c) The old node has branches on both sides. Get the left of oldEvent->right and swap the two + do required linkings.
    */
    particle<real> *parent = oldEvent->parent; //Start from the particle parent initialization
    particle<real> *node;  //A temporary object.

    //a) If the event to be removed does not have a element lower than that attached to it.
    if ( oldEvent->left == NULL ){ 
        node = oldEvent->right; //just remove it and attach the right event to the parent of the current node.
        if (node){ // If the right node is not null, then perform this operation here.
            node->parent = parent; //Shift the parent to the node. But what about the reverse linking... It is done at the last
        //It is not needed, because it is not a doubly linked list, only one way the information needs to be traversed.
        }
    } //b) The old event has just the right branch, no left branch.
    else if (oldEvent->right == NULL){ //Only one child again.
        node = oldEvent->left; //Child1 is attached to parent
        node->parent = parent; //if( node){  node->parent = parent} is not needed here...it is already been tested in the above if loop.
    } // c) The old element has both left ad right branches.
    else{ 
        /*if neither of left and right are NULL, the oldEvent has branches in both the sides.
        a) Find first node having no left in tree after oldEvent. (This node will have no left)
        b) Remove this node from the tree     (Attach node->right to node->parent)
        c) Put this node in place of oldEvent (Oldevent's children are adopted by node)
        */
        node = oldEvent->right; //just remove it and attach the right event to the parent of the current node.
        while (node->left) 
            node = node->left; //Keep on going to the leftmost element untill you find a NULL. 

        //Get the parent of that leftmost node.    
        particle<real> *pnode = node->parent; //It will be greater than the old event but less than the right branch of the old event.

        //Check if the rightnode of the oldEvent does not have an event to the left, then just replace it to oldEvent.
        if (pnode == oldEvent){                         
            oldEvent->left->parent = node;
            node->left = oldEvent->left;
        } //If there are deep branches to the right element of the oldEvent.
        else{ 
            //node is not a child of oldEvent.
            //Both of oldEvent's children should be connected to node.
            pnode->left = node->right; //Make the right event of the old node the left event of parent ( logical as ).
            if (node->right) 
                node->right->parent = pnode;

            //Making the connections.
            oldEvent->left->parent = node;   //Forward.
            node->left = oldEvent->left;     //Backward.
            oldEvent->right->parent = node;  //Forward.
            node->right = oldEvent->right;   //Backward.
        }
        //But the parent of the new node replacing the oldNode has to be allocated as well. Let us do it now.
        node->parent = parent; //No need to go the the last right of the event to be deleted, just go to the first Right only
        //Good intuition, and will be efficient.
    }
    //The parent of the old event is liked to the new member at the place of the old event similarly as it was placed for old event.
    if (parent->left == oldEvent) 
        parent->left = node;
    else
        parent->right = node;
}


/*Create and Event */
template<typename real>
void hardSphere<real>::createEvent( particle<real> *pi, particle<real> * pj, real etime, int type){
    if( time > 0){
        pi->eventTime = etime;
        //std::cout<<"\n pi->eventTime: "<<pi->eventTime;
        pi->eventType = type;
        pi->pj = pj;
        addEvent( pi );
    }
    else{
        std::cout<<"\n pi->eventTime < 0.0: "<<etime;
        //std::cin.get();
    }
}



/*Add a next event list to the event calender *
 * Sorts the events from the next event list.
 * If it is overflowed checks the overflow list, store it in the current list 
 * Push the events in current list into the Binary Tree
 */
template< typename real>
void hardSphere<real>::addNextEventList(){
    //std::cout<<"\n Adding new event list";
    currentList++; //Increase the current list variable by 1.
    //std::cout<<"\n currentList: "<<currentList<<"\t SimTime: "<<simTime;
    refTime += eventsDeltaT; //Progress in the reference time.

    if (currentList == nEventBoxes){ //Reached the end of the array of the event list, start from begining as the events are stored in a PBC manner :)
        currentList = 0; //Reset the current list to 0.

        //Now process the overflow event queue and add it to the current list.
        particle<real> *ei = eventLists[nEventBoxes]; //Processing the overflow list
        eventLists[nEventBoxes] = NULL; //Set the overflow event box to NULL.
        //listCounter2 = 0;
        particle<real> *nextev;
        while (ei){
            nextev = ei->right;
            addEvent(ei);
            ei = nextev;
            //listCounter2 += 1; // Count how many events there were in last event list.
        }
    }

    particle<real> *ei = eventLists[currentList];
    particle<real> *nextev;
    while (ei){
        nextev = ei->right;
        addEventToBinaryTree(ei);
        ei = nextev;
        //listCounter1++; // Count how many events there were in event list
    }

    eventLists[currentList] = NULL; //The elements have been added to the binary tree, so  reset it to zero.
    //mergeCounter++;//The number of times the merge operation/The event list got empited, has been performed.
}


/* Finding a collision time between two particles*/
template< typename real>
int hardSphere<real>::findPairCollisionTime( particle<real> *pi, particle<real> *pj, real *tmin ){
    real dtj = simTime - pj->tlast;
    real dr[DIM] = {0.0}, dv[DIM] = {0.0}, dr2=0.0, b=0.0, dv2=0.0;
    for( int k = 0; k < DIM; k++ ){
        dr[k] = pi->pos(k) - pj->pos(k) - dtj * pj->vel(k);
        //if( pi->nearBoxEdge ){
            if( dr[k] > 0.5 * boxLen )
                dr[k] -= boxLen;
            else if( dr[k] < -0.5 * boxLen )
                dr[k] += boxLen;
        //}
        dv[k] = pi->vel(k) - pj->vel(k);
        dr2 += dr[k] * dr[k];
        dv2 += dv[k] * dv[k];
        b   += dr[k] * dv[k];
    }//k ends.

    if( b >= 0.0){ //Particles flying apart.
        return 0;
    }

    real A = 1 - dr2; //1 is the diameter of the particles.
    if( 2 * b * *tmin > A ){
        return 0; //check to discard hits at times > tmin. Note than b<0 here.
    }

    real tcol;
    if( dr2 <= 0.99999999 ){ //Overlapping particles, hence the collision time will come out to be negative here.
        std::setprecision(16);
        std::cout<<"\n dr_ij2 = "<<dr2;
        printPosVel( pi ); printPosVel( pj );
        std::cout<<"\n Press enter to continue"; //std::cin.get();
        tcol = 0.0;
        *tmin = tcol;
        return 1;
        
    }

    real disc = b * b + dv2 * A;
    if( disc < 0 ){
        return 0; //Particles have closest distance approach greater than the sum of radii.
    }

    //real tcol = ( ( -b - sqrt( disc ) ) / dv2 );
    tcol = -A / ( -b + sqrt( disc ) );
    //std::cout<<"\n tcol1 : "<<tcol <<"\t *tneig: "<<*tmin;

    if( tcol <= *tmin ){
        *tmin = tcol; //Particle collision time is greater than the wall other update times.
        //std::cout<<"\n tmin : "<<*tmin;
        if( *tmin < 0.0  ){
            std::cout<<"\n *tmin < 0 encountered: "<<*tmin;
            printPosVel( pi ); printPosVel( pj );
            //std::cin.get();
        }
        return 1;
    }

    return 0; //If non of these is happening, make sure that the event is not scheduled.
}


/* Find the minimum collision time for particle pi */
template<typename real>
void hardSphere<real>::findMinCollTimeEvent( particle<real> *pi ){
    real tmin = findNeighborlistUpdate( pi ); 
    //std::cout<<"\n tmin: "<<tmin;
    particle<real> *pj, *partner = pi; //Default partner initialized.
    int etype = 2; //Neighborlist update type.
    int iI = pi->index;
    int iJ;
    for( unsigned int i = 0; i < on[iI].neigh.size(); i++ ){
        iJ = on[iI].neigh[i];
        //std::cout<<"\n neig count i: "<<i<<"\t iJ: "<<iJ;
        pj = &op[iJ]; //std::cout<<"\t pj->index: "<<pj->index;
        if( findPairCollisionTime( pi, pj, &tmin ) ){
            partner = pj;
            etype = 1; //Pair collision type.
        }
    }
    //std::cout<<"\n tmin: "<<tmin;
    //std::cout<<"\n Found minCollTimeEvent for pi: "<<pi->index;
    createEvent( pi, partner, tmin+simTime, etype );
    pi->counterPj = partner->counter;
}


// Execute a collision event of a particle.
template< typename real >
void hardSphere<real>::executeCollision( particle<real> *pi ){
    particle<real> *pj = pi->pj;
    updatePosition( pi );

    if( pi->counterPj != pj->counter ){
        //The partner of pi has undergone collision with other particles, the event is invalid, find all collisions for the 'pi' particle.
        findMinCollTimeEvent( pi );
        return;
    }
    //If partner has not undergone collision with other particle. The update position includes the PBC terms as well. 
    updatePosition( pj );

    //Increase collision counter.
    //std::cout<<"\n Executing collisions: ";
    pi->counter++; pj->counter++;

    real dr[DIM] = {0.0, 0.0, 0.0}, dv[DIM] = {0.0, 0.0, 0.0};
    real b=0.0, dr2=0.0;

    for( int k = 0;  k < DIM; k++ ){
        dr[k] = pi->pos(k) - pj->pos(k);
        //if( pi->nearBoxEdge ){
            if( dr[k] > 0.5 * boxLen )
                dr[k] -= boxLen;
            else if( dr[k] < -0.5 * boxLen )
                dr[k] += boxLen;
        //}
        dv[k] = pi->vel(k) - pj->vel(k);
        dr2 += dr[k] * dr[k];

        b += dr[k] * dv[k];
    }

    //If the particles are overlapping, pop up a warning message.
    if( dr2 < 0.999999999){
        std::cout<<"\n Overlapping particles, but colliding now. i: "<<pi->index<<"\t j: "<<pj->index<<"\t dr2: "<<dr2;
    }

    real delVij = 0.5 * ( 1 + eps ) * b; 

    #if compileMode==0
    std::cout<<"\n Pre-collision velocities from execute collision, pi: ";
    for( int k=0; k<DIM; k++ ){
        std::cout<<pi->vel(k)<<"\t";
    }
    std::cout<<"\n pj: ";
    for( int k=0; k<DIM; k++ ){
        std::cout<<pj->vel(k)<<"\t";
    }

    std::cout<<"\n Relative position from executeCollision: ";
    for( int k=0; k<DIM; k++ ){
        std::cout<<dr[k]<<"\t";
    }

    #endif

    //Update velocities of particles.
    for( int k = 0; k < DIM; k++ ){
        pi->vel(k) -= delVij * dr[k];
        pj->vel(k) += delVij * dr[k];
    }

    #if compileMode==0
    std::cout<<"\n Post-collision velocities from execute collision, pi\n";
    for( int k=0; k<DIM; k++ ){
        std::cout<<pi->vel(k)<<"\t";
    }
    std::cout<<"\n pj \n";
    for( int k=0; k<DIM; k++ ){
        std::cout<<pj->vel(k)<<"\t";
    }
    #endif

    // Save the virial for the pressure calculation.
    brTot += delVij;
    totalCollEvents++;


    // Sampling MOP potential energy contribution term.
    real halfBoxLen = 0.5 * boxLen;
    if( pi->pos(1) < halfBoxLen && pj->pos(1) > halfBoxLen )
        mopUcontr += abs( 2 * delVij*dr[1] ); //The absolute value is added latter.
    else if( pi->pos(1) > halfBoxLen && pj->pos(1) < halfBoxLen )
        mopUcontr += abs( (-2 * delVij*dr[1]) );

    //Sample for the Viscosity and Thermal conductivity.
    // sampleCorrelations( pi, pj, delVij, simTime );

    //Save the next collision time value.
    if(makeNextCollTime){
        nextCollTime.push_back( pi->eventTime ); //Store the next collision time here.
        nextCollTime.push_back( pi->index );
        //Bring the particles back in box before saving their positions. Rapped coordinates.
        
        nextCollTime.push_back( pi->pos(0) - boxLen * floor(pi->pos(0) / boxLen) ); 
        nextCollTime.push_back( pi->pos(1) - boxLen * floor(pi->pos(1) / boxLen) );
        nextCollTime.push_back( pi->pos(2) - boxLen * floor(pi->pos(1) / boxLen) );

        nextCollTime.push_back( pj->index );
        nextCollTime.push_back( pj->pos(0) - boxLen * floor(pj->pos(0) / boxLen) );
        nextCollTime.push_back( pj->pos(1) - boxLen * floor(pj->pos(1) / boxLen) );
        nextCollTime.push_back( pj->pos(2) - boxLen * floor(pj->pos(2) / boxLen) );

        nctCount += 9;
    }


    //Remove the event from the event calender. pi is already removed in the driver code.
    removeEvent( pj ); 
    findMinCollTimeEvent( pi); findMinCollTimeEvent( pj );
    //std::cout<<"\n Executed Collision between pi: "<<pi->index<<"\t pj: "<<pj->index<<"\t eventTime: "<<pi->eventTime;
}


/*Function to save the average of y component of velocity*/
template<typename real>
void hardSphere<real>::saveYvel(){
    std::fstream fiO;
    char fileName[300];
    sprintf( fileName, "meanOf_velInYdirection_N%d_pf%.4lf_v%.4lf_eNormal%.2lf_cellLen%.2lf_maxTime%.1lf.dat", npart, packFrac, eps, cellLen, boxLen*boxLen*boxLen, maxTime );
    fiO.open( fileName, std::ios::app );
    fiO<< std::fixed << std::setprecision(16);
    if( trial == 0){
        fiO<<"#The first row is simulation time, the second row onwards are the average of the y componenet velocity of cells near z = 0 at the specified simulation time in the first row\n";
        for( int i = 0; i < simT.size(); i++ ){
            fiO<<simT[i]<<"\t";
        }
        fiO<<std::endl;
        for( int i = 0; i < yVel.size(); i++ ){
            fiO<<yVel[i]<<"\t";
        }
        fiO<<std::endl;
    }
    else{
        for( int i = 0; i < yVel.size(); i++ ){
            fiO<<yVel[i]<<"\t";
        }
        fiO<<std::endl;
    }
    fiO.close();
    //std::cout<<"\n yVel.size(): "<<yVel.size();
    //yVel.clear();
    //simT.clear();
}



// Function to sample the observables.
template<typename real>
void hardSphere<real>::sampleObservables( particle<real> *we ){
    //sampleMeanSqDistance();
    //sampleConfigForGr();
    real ptot = samplePressure();
    real ptotMop = samplePressureMop(); //std::cout<<"\n"<<simTime<<"\t"<<ptotMop;
    pressArr.push_back(ptot);
    pressArrMop.push_back( ptotMop );
    //sampleShearViscosity();
    //sampleDiagonalThermalConductivity();
    simT.push_back(simTime);


    //Write the Next Collision Time data to a file.
    static bool counter = 1;
    std::fstream fiO;
    char fileName[300];

    if(makeNextCollTime){
        sprintf( fileName, "nextCollTime_N%d_pf%.4lf_eNormal%.2lf_cellLen%.2lf_maxTime%.1lf.dat", npart, packFrac, eps, cellLen, maxTime );
        
        if(counter){ 
            fiO.open( fileName, std::ios::out );
            fiO<<"#NextCollTime\t pi->index\t pi->x\t pi->y \tpi->z\t pj->index\t pj->x\t pj->y\t pj->z\n";  
            fiO<<std::endl;      
        }
        else fiO.open( fileName, std::ios::app );

        fiO<< std::fixed<<std::setprecision(16);
        for( int i = 0; i<nctCount; i = i+9 ){
            for( int k = 0; k<9; k++){
                fiO<<nextCollTime[i+k]<<"\t";
            }
            fiO<<std::endl;
        }
        fiO.close();
        nextCollTime.clear();
        nctCount=0; //resetting it to zero.
    }
    counter = 0;

    
    //Schedule the next write event.
    we->eventTime = simTime + sampleObsInterval;
    we->eventType = 3; //This is the event type for write intervals.
    we->pj = NULL;
    addEvent( we ); 
    totalSamplingEvents++;

    //Periodically save the sampled observables.
    if( totalSamplingEvents % 5 == 0){
        saveObservables();
        systemReport();

        //Reset the simulation time array.
        simT.clear(); pressArr.clear(); msdArr.clear(); 
        pressArrMop.clear(); //This is added latter in this code, check if its creating problem.
        viscosityArr.clear(); thermalCondArr.clear();
        pressSplits.clear(); //This was there, hence safe with the pressure splits.
    }
}


//Save the observables
template<typename real>
void hardSphere<real>::saveObservables(){
    savePressure();
    //saveConfig(); //It will save only the last configuration.
    //saveMeanSqDistance();
    //saveViscosity();
    //saveDiagonalThermalConductivity();
    //sampleConfigForGr();
}



/*Sample the pressure from system*/
template<typename real>
real hardSphere<real>::samplePressure(){
    
    double v2 = 0.0; 
    for( int i = 0; i < npart; i++ ){
        for( int k = 0; k < DIM; k++ ){
            v2 += op[i].vel(k) * op[i].vel(k);
        }
    }

    double pid = ( 1/(3.0*boxVolm) ) * v2; //Mass of the particle is taken to be equal to 1.
    // real pid = rho; //Next time calcualte with the above mentioned procedure.
    real pvir =  -(brTot - brTotPrev) / ( 3.0 * boxVolm * ( simTime - timePrev ) );
    real ptot = pid + pvir;
    timePrev = simTime;
    brTotPrev = brTot;
    //std::cout<<"\n pid: "<<pid<<"\tpvir: "<<pvir<<"\tptot: "<<ptot<<"\tbrTot: "<<brTot;
    pressSplits.push_back( pid ); 
    pressSplits.push_back( pvir );
    return ptot;
}


template<typename real>
real hardSphere<real>::samplePressureMop(){
    real Atau = ( boxLen*boxLen*( simTime - mopprevTime ) );
    real pKin = ( mopTransfer - prevmopTransfer ) /  Atau ;
    real pPot = ( mopUcontr - prevmopUcontr )/( 4*Atau );
    
    prevmopUcontr = mopUcontr;
    prevmopTransfer = mopTransfer;
    mopprevTime = simTime;
    pressSplits.push_back( pKin );
    pressSplits.push_back( pPot );

    return (pKin+pPot);
}


// Sample mean square distance for diffusivity calculations.
template<typename real>
void hardSphere<real>::sampleMeanSqDistance(){
    int i; int k;
    real disp = 0.0, disp2 = 0.0;
    particle<real> *pi = op;
    for( i = 0; i < npart; i++ ){
        updatePosition( pi ); //pi->setPbc( boxLen ); pi++;
        
        for( k = 0; k < DIM; k++ ){
            disp = op[i].pos(k) - opRef[i].pos(k) + ( op[i].boxesTravelx[k] * boxLen );
            disp2 += disp * disp;
        }
    }
    disp2 /= (1.0*npart);
    msdArr.push_back( disp2 );
}


/*sample configurations for g(r)*/
template< typename real>
void hardSphere<real>::sampleConfigForGr(){
    real grDelT = maxTime / numOfConfigForGr;
    static int numOfConfigSaved = 1;
    if( simTime > numOfConfigSaved * grDelT ){
        saveConfig();
        numOfConfigSaved += 1;
    }
}


//Sample the different correlation components for viscosity and thermal conductivity. The updated position and velocity is here post collision.
template<typename real>
void hardSphere<real>::sampleCorrelations( particle<real> *pi, particle<real> *pj, real delVi, real eventTime ){
    //Calculate the second term of the correlation coefficient.
    real dviBeta[DIM], drijAlpha[DIM]; //These are post collision.
    real viPrev[DIM], vjPrev[DIM], delVi2 = 0.0; //For the difference in the velocity.
    real keIi = 0.0, keIj = 0.0, keFi = 0.0, keFj = 0.0;

    #if compileMode==0
    real drij2 = 0.0;
    #endif 

    real deltaTau = eventTime - prevCollTime;
    prevCollTime = eventTime;

    int k;
    for( k = 0; k < DIM; k++ ){
        //Viscosity.
        drijAlpha[k]  = pi->pos(k) - pj->pos(k); //Position of the particles does not changes.
        if( drijAlpha[k] > 0.5 * boxLen )
            drijAlpha[k] -= boxLen;
        else if( drijAlpha[k] < -0.5 * boxLen )
            drijAlpha[k] += boxLen;
        //drijAlpha[k] -= boxLen * lround( drijAlpha[k]/boxLen ); //The PBC being employed here.

        dviBeta[k]   = -delVi * drijAlpha[k]; //Change in ith particle velocity due to collision.
        
        //Thermal conductivity.
        viPrev[k]    = pi->vel(k) + delVi * drijAlpha[k]; //- changes to + when taken the other side. pi->vel() gives the updated velocity of particle i.
        keIi        += 0.5 * viPrev[k] * viPrev[k]; 
        keFi        += 0.5 * pi->vel(k) * pi->vel(k);
        vjPrev[k]    = pj->vel(k) - delVi * drijAlpha[k]; 
        keIj        += 0.5 * vjPrev[k] * vjPrev[k];
        keFj        += 0.5 * pj->vel(k) * pj->vel(k);
        delVi2      += ( pi->vel(k) * pi->vel(k) - viPrev[k] * viPrev[k] );
        
        #if compileMode==0
        drij2 += drijAlpha[k]*drijAlpha[k];
        #endif
    }

    #if compileMode==0
    std::cout<<"\n Pre-collision velocities from sampleCorrelations, pi: ";
    for( int k=0; k<DIM; k++ ){
        std::cout<<viPrev[k]<<"\t";
    }
    std::cout<<"\n pj: ";
    for( int k=0; k<DIM; k++ ){
        std::cout<<vjPrev[k]<<"\t";
    }
    std::cout<<"\n Relative position from sampleCorrelations: ";
    for( int k=0; k<DIM; k++ ){
        std::cout<<drijAlpha[k]<<"\t";
    }
    #endif


    #if compileMode==0
    std::cout<<"\n Post-collision velocities from sampleCorrelations, pi: ";
    for( int k=0; k<DIM; k++ ){
        std::cout<<pi->vel(k)<<"\t";
    }
    std::cout<<"\n pj: ";
    for( int k=0; k<DIM; k++ ){
        std::cout<<pj->vel(k)<<"\t";
    }
    #endif

    #if compileMode==0
    if( drij2 < 0.999999999999 || drij2> 1.00000000001){
        std::cout<<"\n drij2 != unity "<<std::setprecision(16)<<drij2;
    }
    std::cout<<"\n Checking if the percentage operation is implemented correctly";
    std::cout<<"\n The initial value: ";
    for( int k=0; k<DIM; k++ ){
        std::cout<<"\t"<<dviBeta[k];
    }
    std::cout<<"\n One shifted value";
    for( int k=0; k<DIM; k++ ){
        std::cout<<"\t"<<dviBeta[(k+1)%DIM];
    }

    #endif

    //The three cross components.
    for( k = 0; k < DIM; k++ ){
        //Viscosity.
        viscosityT2[k]   = drijAlpha[k] * dviBeta[(k+1)%DIM];
        //Thermal conductivity.
        thermalCondT2[k] = 0.5 * delVi2 * drijAlpha[k];
    }

    //Getting the Hefland moments.
    for( k = 0; k < DIM; k++ ){
        //Viscosity
        viscosityAxx[k]   += viscosityT1[k] * deltaTau + viscosityT2[k];
        //Thermal Conductivity
        thermalCondAxx[k] += thermalCondT1[k]*deltaTau + thermalCondT2[k];
    }
    
    //Check for the viscosityAxx data, if the integration of the stress-stress correlations function comes out to be an exponentially decaying function.
    //real corrOut = viscosityAxx[0];
    //std::cout<<simTime<<"\t"<<( corrOut*corrOut/( 2 * simTime * boxVolm * boxVolm ) )  <<"\n";




    //Update T1 arrays. Here the velocities are the updated velocities after collision, be careful.
    for( k = 0; k < DIM; k++ ){
        //Viscosity.
        //viscosityT1[k] += -delVi * ( delVi * drijAlpha[k] * drijAlpha[(k+1)%DIM] + pi->vel(k) * drijAlpha[(k+1)%DIM] + pi->vel((k+1)%DIM) * drijAlpha[k] );  //Done for ith particle.
        //viscosityT1[k] += delVi * ( -delVi * drijAlpha[k] * drijAlpha[(k+1)%DIM] + pj->vel(k) * drijAlpha[(k+1)%DIM] + pj->vel((k+1)%DIM) * drijAlpha[k] ); //For the jth particle.
        
        viscosityT1[k] += ( pi->vel(k) * pi->vel( (k+1)%DIM ) + pj->vel(k) * pj->vel( (k+1)%DIM ) - viPrev[k] * viPrev[ (k+1)%DIM ] - vjPrev[k] * vjPrev[ (k+1)%DIM ] );

        //Thermal conductivity.
        thermalCondT1[k] += ( keFi * pi->vel(k) + keFj * pi->vel(k) - keIi * viPrev[k] - keIj * vjPrev[k] );   
    }

    //The sampling of the viscosity and thermalConductivity will be done in the sampleViscosity function and the sampleThermalConductivity.
}


//Sample the viscosity.
template<typename real>
void hardSphere<real>::sampleShearViscosity(){
    //Do the integration of the viscosity from the Green-Kubo relation.
    real deltaAalphaBeta;
    for( int k = 0; k < DIM; k++ ){
        deltaAalphaBeta = viscosityAxx[k] - viscosityAxxPrev[k];
        viscosityAxxPrev[k] = viscosityAxx[k];
        
        //From the Einstein's relation there.
        deltaAalphaBeta = ( deltaAalphaBeta*deltaAalphaBeta/( 2 * boxVolm * temp * ( simTime - prevViscTime ) ) );

        viscosityArr.push_back( deltaAalphaBeta );
    }
    //Push the last element to the viscosity calculation.
    prevViscTime = simTime;
}


//Sample thermal conductivity of the particles along the xx, yy, and the zz direction ( the diagonal directions )
template<typename real>
void hardSphere<real>::sampleDiagonalThermalConductivity(){
    real deltaAalphaBeta;
    for( int k = 0; k < DIM; k++ ){
        deltaAalphaBeta       = thermalCondAxx[k] - thermalCondAxxPrev[k];
        thermalCondAxxPrev[k] = thermalCondAxx[k];
        
        //From the Einstein's relation.
        deltaAalphaBeta = ( deltaAalphaBeta*deltaAalphaBeta/( 2 * boxVolm * temp * temp * ( simTime - prevThermTime ) ) );

        thermalCondArr.push_back( deltaAalphaBeta );
    }

    //Push the last element to the viscosity calculation.
    prevThermTime = simTime;
}


//Save the pressure into a data file.
template<typename real>
void hardSphere<real>::savePressure(){
    std::fstream fiO;
    char fileName[300];
    static bool pcounter = 1; static bool psplitcounter = 1;
    sprintf( fileName, "pressure_N%d_pf%.4lf_maxTime%.0lf_trial%d.dat", npart, packFrac, maxTime, trial );
    fiO.open( fileName, std::ios::app );
    fiO<< std::fixed<<std::setprecision(16);
    //if( trial == 0){
        //fiO<<"#The first row is simulation time, the second row onwards are the pressure at different time instances of the given trial"<<std::endl;
        if( pcounter == 1 ){
            fiO<<"#SimTime(t) \t Pressure(t) \t PressureMop(t)"<<std::endl;
            pcounter = 0;
        }
        
        for( unsigned int i = 0; i < simT.size(); i++ ){
            fiO<<simT[i]<<"\t"<<pressArr[i]<<"\t"<<pressArrMop[i]<<"\n";
        }
    fiO.close();

    #if ObsSplit==1
    sprintf( fileName, "pressureSplit_N%d_pf%.4lf_maxTime%.0lf_trial%d.dat", npart, packFrac, maxTime, trial );
    fiO.open( fileName, std::ios::app );
    fiO<< std::fixed<<std::setprecision(16);
    //if( trial == 0){
        //fiO<<"#The first row is simulation time, the second row onwards are the pressure at different time instances of the given trial"<<std::endl;
        if( psplitcounter == 1 ){
            fiO<<"#SimTime(t) \t Pressure(t) Kin \t Pressure(t) Pot \t Pressure(t) \t PressureMop(t) Kin \t PressureMop(t) Pot \t PressureMop(t)"<<std::endl;
            psplitcounter = 0;
        }
        
        for( unsigned int i = 0; i < simT.size(); i++ ){
            fiO<<simT[i];
            for( unsigned int k = 0; k < 2; k++ ){
            fiO<<"\t"<<pressSplits[i*4+k];
            }
            fiO<<"\t"<<( pressSplits[i*4+0] + pressSplits[i*4+1] );

            for( unsigned int k = 2; k < 4; k++ ){
            fiO<<"\t"<<pressSplits[i*4+k];
            }
            fiO<<"\t"<<( pressSplits[i*4+2] + pressSplits[i*4+3] );

            fiO<<"\n";
        }
    fiO.close();
    #endif
}



/*Save the configuration of particle configurations*/
template<typename real>
void hardSphere<real>::saveConfig(){
    //std::cout<<"\n Saving the configurations";
    std::fstream fiO;
    char fileName[300];
    sprintf( fileName, "config_N%d_pf%.4lf_simTime%.1lf_trial%d.xyz", npart, packFrac, simTime, trial );
    fiO.open( fileName, std::ios::out );
    fiO<< std::fixed<<std::setprecision(16);
    
    //fiO<<"#npart, boxLx, boxLy, boxLz \n#x\ty\tz\tvx\tvy\tvz"<<std::endl;

    fiO<<npart<<std::endl;
    fiO<<boxLen<<"\t"<<boxLen<<"\t"<<boxLen<<std::endl;

    for( int i = 0; i < npart; i++ ){
        for( int k = 0; k < DIM; k++ )
            fiO<<op[i].pos(k)<<"\t";
        
        for( int k = 0; k < DIM; k++ )
            fiO<<op[i].vel(k)<<"\t";
        
        fiO<<std::endl;
    }
    fiO.close();
}



/*Save the pressure into a data file*/
template<typename real>
void hardSphere<real>::saveMeanSqDistance(){
    std::fstream fiO;
    char fileName[300];
    static bool msdcounter = 1;
    sprintf( fileName, "msd_N%d_pf%.4lf_maxTime%.0lf_trial%d.dat", npart, packFrac, maxTime, trial );
    fiO.open( fileName, std::ios::app );
    fiO<< std::fixed<<std::setprecision(16);

    if( msdcounter == 1 ){
        fiO<<"#SimTime(t) \t Msd(t)"<<std::endl;
        msdcounter = 0;
    }

    for( unsigned int i = 0; i < simT.size(); i++ ){
        fiO<<simT[i]<<"\t"<<msdArr[i]<<"\n";
    }
    fiO.close();
}


//Save the time series of viscosity 
template<typename real>
void hardSphere<real>::saveViscosity(){
    std::fstream fiO;
    char fileName[300];
    static bool visccounter = 1;
    sprintf( fileName, "viscosity_N%d_pf%.4lf_maxTime%.0lf_trial%d.dat", npart, packFrac, maxTime, trial );
    fiO.open( fileName, std::ios::app );
    fiO<<std::fixed<<std::setprecision(16);
    real meanVisc;
    if( visccounter == 1 ){
        fiO<<"#SimTime(t) \t Viscosity(xy) \t Viscosity(yz) \t Viscosity(zx) \t meanViscosity( (xy+yz+zx)/3)"<<std::endl;
        visccounter = 0; 
    }
    
    for( unsigned int i = 0; i < simT.size(); i++ ){
        fiO<<simT[i];
        meanVisc = 0.0;
        for(int k = 0; k < DIM; k++ ){
            fiO<<"\t"<<viscosityArr[DIM*i+k];
            meanVisc += viscosityArr[DIM*i+k];
            //Printing on screen.
            //std::cout<<"\t"<<viscosityArr[DIM*i+k];
        }
        meanVisc /= DIM;
        fiO<<"\t"<<meanVisc;
        fiO<<std::endl;
    }
    fiO.close();
}


//Save the time series of viscosity. 
template<typename real>
void hardSphere<real>::saveDiagonalThermalConductivity(){
    std::fstream fiO;
    char fileName[300];
    static bool kapacounter = 1;
    sprintf( fileName, "thermalCond_N%d_pf%.4lf_maxTime%.0lf_trial%d.dat", npart, packFrac, maxTime, trial );
    fiO.open( fileName, std::ios::app );
    fiO<<std::fixed<<std::setprecision(16);
    real meanthermalCond = 0.0;

    if( kapacounter == 1 ){
        fiO<<"#SimTime(t) \tthermalCond(xx) \t thermalCond(yy) \t thermalCond(zz) \t meanthermalCond( (xx+yy+zz)/3)"<<std::endl;
        kapacounter = 0;
    }

    for( unsigned int i = 0; i < simT.size(); i++ ){
        fiO<<simT[i]; 
        meanthermalCond = 0.0;
        for(int k = 0; k < DIM; k++ ){
            fiO<<"\t"<<thermalCondArr[DIM*i+k];
            meanthermalCond +=thermalCondArr[DIM*i+k];
            //Printing on screen.
            //std::cout<<"\t"<<viscosityArr[DIM*i+k];
        }
        meanthermalCond /= DIM;
        fiO<<"\t"<<meanthermalCond;
        fiO<<std::endl;
    }
    fiO.close();
}


//Print the variables.
template<typename real>
void hardSphere<real>::printVariables(){
    std::cout<<"\n Displaying the varaibles of the simulations ";
    std::cout<<"\n npart: "<<npart;
    std::cout<<"\n packFrac: "<<packFrac;
    std::cout<<"\n temp: "<<temp;
    std::cout<<"\n eps: "<<eps;
    std::cout<<"\n maxTime: "<<maxTime;
    std::cout<<"\n sampleObsInterval: "<<sampleObsInterval;
    std::cout<<"\n shellSize: "<<shellSize;
    std::cout<<"\n posInit: "<<posInitType;
    std::cout<<"\n fileName: "<<fileName;
    std::cout<<"\n useThermostat: "<<useThermostat;
    std::cout<<"\n trial: "<<trial;

    std::cout<<"\n avgFFT: "<<avgFFT;
    std::cout<<"\n avgFPL: "<<avgFPL;
    std::cout<<"\n avgNCT: "<<avgNCT;
    std::cout<<"\n boxVolm: "<<boxVolm;
    std::cout<<"\n boxLen: "<<boxLen;
    std::cout<<"\n simTime: "<<simTime;
    //std::cout<<"\n nCellx: "<<nCellx;
    //std::cout<<"\n nCell: "<<nCell;
    //std::cout<<"\n cellLen: "<<cellLen;
    std::cout<<"\n ke/npart: "<<ke/npart;
    std::cout<<"\n brTot: "<<brTot;
    std::cout<<"\n timePrev: "<<timePrev;
}


/* Taking the simulation report */
template<typename real>
void hardSphere<real>::systemReport(){
    //particle<real> *p = &op[0]; //Just for the sake of initializations.
    double v2tot = 0;


    //The final write event before closing the simulation.
    //printf("\n Writing final ");
    //sampleObservables(p);
    //totalSamplingEvents++; totalEvents++;
    //Calculate the total kinetic energy.
    double v2 = 0.0; 
    for( int i = 0; i < npart; i++ ){
        for( int k = 0; k < DIM; k++ ){
            v2 += op[i].vel(k) * op[i].vel(k);
        }
    }

    double pressid = ( 1/(3.0*boxVolm) ) * v2; //Mass of the particle is taken to be equal to 1.
    //double pid = rho * temp; //Next time calcualte with the above mentioned procedure.
    double pvir =  -(brTot ) / ( 3.0 * boxVolm * ( simTime ) );

    double presstot = pvir + pressid;

    double nCollPerSec = totalCollEvents/cpuTime;

    //Shear Viscosity
    real deltaAalphaBeta, meanShearVisc = 0.0;
    for( int k = 0; k < DIM; k++ ){
        deltaAalphaBeta = viscosityAxx[k] - 0.0;
        //viscosityAxxPrev[k] = viscosityAxx[k];
        
        //From the Einstein's relation there.
        meanShearVisc += ( deltaAalphaBeta*deltaAalphaBeta/( 2 * boxVolm * temp * ( simTime - 0.0 ) ) );

        //viscosityArr.push_back( deltaAalphaBeta );
    }
    meanShearVisc /= 3;

    //Mean thermal conductivity.
    real meanThermalCond = 0.0;
    for( int k = 0; k < DIM; k++ ){
        deltaAalphaBeta       = thermalCondAxx[k] - 0.0;
        //thermalCondAxxPrev[k] = thermalCondAxx[k];
        
        //From the Einstein's relation.
        meanThermalCond += ( deltaAalphaBeta*deltaAalphaBeta/( 2 * boxVolm * temp * temp * ( simTime - 0.0 ) ) );

        //thermalCondArr.push_back( deltaAalphaBeta );
    }
    meanThermalCond /= 3;


    if(makeEfficiencyReport){
        FILE *file;
        char fileName[300];

        sprintf(fileName, "systemEfficiency_N%d_pf%.4lf_eNormal%.2lf.dat", npart, packFrac, eps );
        file = fopen(fileName, "a");
        fprintf(file, "%d \t %lf \t %lf \t %lf \t %d  \t %d \t %lf\n", npart, ceil(nCollPerSec),simTime , cpuTime, totalCollEvents, totalEvents, packFrac);
        fclose(file);
    }

    if(saveSystemReport){
        FILE *file;

        char fileName[300];
        //File to save the system information.
        sprintf( fileName, "systemInfo_N%d_pf%.4lf_eNormal%.2lf.dat", npart, packFrac, eps );
        file = fopen( fileName, "w");
        fprintf(file, "Number of particles      : %d\n", npart);
        fprintf(file, "Packing fraction         : %lf\n", packFrac);
        fprintf(file, "Density                  : %lf\n", (double) npart / boxVolm);
        fprintf(file, "Volume                   : %lf\n", boxVolm );
        fprintf(file, "Box Len                  : %lf\n", cbrt(boxVolm));
        fprintf(file, "shellSize                : %lf\n", shellSize);
        fprintf(file, "Average kinetic energy   : %lf\n", 0.5 * v2tot / npart);
        fprintf(file, "Total time simulated     : %lf\n", maxTime);
        fprintf(file, "Measured pressure        : %lf + %lf = %lf\n", pvir, pressid, presstot);
        fprintf(file, "Average shear viscosity  : %lf\n", meanShearVisc );
        fprintf(file, "Average thermal conduct  : %lf\n", meanThermalCond );
        fprintf(file, "Collision/s              : %lf\n", nCollPerSec);
        fprintf(file, "CpuRunTime(s)            : %lf\n", cpuTime);
        fprintf(file, "simRunTime               : %lf\n", maxTime);
        fprintf(file, "Total Events             : %d\n", totalEvents);
        fprintf(file, "Total collision events   : %d\n", totalCollEvents );
        fprintf(file, "Theoretical Collisions   : %d\n", int(avgOmega*simTime) );
        fprintf(file, "total neighborhood events: %d\n", totalNeighEvents);
        fprintf(file, "total thermostat events  : %d\n", totalThermoEvents);
        fprintf(file, "total sampling events    : %d\n", totalSamplingEvents);
    }

    if(saveFinalConfig){
        saveConfig();
    }
}


