/*Aim: The definition file for hardSphere.h; Included definitions of the functions
    **
    **
    **
*/
#include "hardSphere.h"

#include<iostream>
#include<fstream>
#include<mm_malloc.h>
#include<numeric>

//#define M_PI 3.1415926535897932

/*class constructor*/
template<typename real>
hardSphere<real>::hardSphere( struct parameter *p ){
    //std::cout<<"\n Constructing hardSphere class variables";

    npart     = p->npart;
    packFrac  = p->packFrac;
    temp      = p->temp;
    

    lambdaMult    = p->lambdaMult;
    op = (particle<real>*)_mm_malloc( npart * sizeof(particle<real>), 4096);

    //std::cout<<"\n sizeof(particle<real>): "<<sizeof(particle<real>);
    boxVolm = npart * M_PI / ( packFrac * 6);
    
    dtFrac  = p->dtFrac;
    simTime = 0.0;
    boxLen  = std::cbrt(boxVolm);
    eps     = 1.0; //Coefficient of restitution.
    wallTemp= p->wallTemp; //std::cout<<"\n wallTemp: "<<wallTemp;
    wallVel = p->wallVel; //std::cout<<"\n wallVel: "<<wallVel;
    trial = p->trial;

    rho = 6 * packFrac / M_PI ; 

    unsigned long seed = time(NULL); //std::cout<<"\n seed "<<seed;
    srand48(seed); //Seed for the random number generator.
    
    initPosition( 0 );
    initVelocity( 0.0 , temp ); //Mean velocity is 0, and standard-deviation is of the temperature here.
    //displayPosVel(10);
    
    calcStatisticalQuantities( packFrac, temp );
    vmax = 5.0 * avgSpdSim;

    dt    = dtFrac * avgFFT; //std::cout<<"\n dt: "<<dt;
    //dt    = 0.01 * boxLen / avgSpd; //std::cout<<"\ndt: "<<dt;
    //std::cout<<"\n avgFFT: "<<avgFFT;
    //dt = 0.2 * vmax * cellLen; //WHY THIS CHOICE IS MADE HERE.
    cellLen = lambdaMult * avgFPL; //std::cout<<"\n cellLen: "<<cellLen<<"\t boxLen: "<<boxLen;
    nCellx  = boxLen / cellLen; //std::cout<<"\n nCellx: "<<nCellx;
    nCell   = nCellx * nCellx * nCellx;
    oc = (cellList<real>*)_mm_malloc( nCell * sizeof(particle<real>), 4096 );

    maxTime = p->maxTime * avgFFT; //std::cout<<"\n maxTime: "<<maxTime;  //The number times the free flight time for each particle there.
    saveInterval = 0.01 * maxTime;
}


/*class destructor*/
template<typename real>
hardSphere<real>::~hardSphere(){
    //std::cout<<"\n A destructor variable called here";
    _mm_free(op);
    _mm_free(oc);
}


/*initialize postitions randomly*/
template<typename real>
void hardSphere<real>::initPosition(int initType){
    for( unsigned int i = 0; i < npart; i++){
        for( int k = 0; k < dim; k++){
            op[i].pos(k) = boxLen * drand48();
        }
    }
    //std::cout<<"\n Positions initialized";
}


/*initialize velocities according to MB distribution*/
template<typename real>
void hardSphere<real>::initVelocity( real meanVel, real temp ){
    real vcm[3] = {0.0, 0.0, 0.0};
    for( unsigned int i = 0; i < npart; i++){
        for( unsigned int k = 0; k < dim; k++){
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
    for( unsigned int i = 0; i < npart; i++){
        for( unsigned int k = 0; k < dim; k++){
            op[i].vel(k) -= vcm[k];
            vcm2 += op[i].vel(k) * op[i].vel(k);
        }
    }
    //Scale to match the temperature and the mean value.
    double localTemp = vcm2 / (npart*3.0);
    localTemp = (temp / localTemp);
    for( unsigned int i = 0; i < npart; i++){
        for( unsigned int k = 0; k < dim; k++){
            op[i].vel(k) = meanVel + sqrt(localTemp) * op[i].vel(k);
        }
    }
    //std::cout<<"\n Kinetic energy: "<<kineticEnergy();
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


/* Getting the value for dt, from an exponential random number generator */
template< typename real >
real hardSphere< real > ::getdt(){
    return -log( drand48()  ) * avgNct;
}



/*Display position and velocity of first 'num' number of particles on screen*/
template<typename real>
void hardSphere<real>::displayPosVel( int num){
    std::cout<<"\n #particle \t pos(x) \t pos (y) \t pos(z) \t vel (x) \t vel(y) \t vel(z)";
    std::cout<<"\n #boxLe: "<<boxLen;
    for( int i = 0; i < num; i++ ){
        std::cout<<"\n"<<i<<"\t";
        for( int k = 0; k < dim; k++ ){
            std::cout<<op[i].pos(k)<<"\t";
        }
        for( int k = 0; k < dim; k++ ){
            std::cout<<op[i].vel(k)<<"\t";
        }
    }
}



/*Calculate the total kinetic energy*/
template<typename real>
real hardSphere<real>::kineticEnergy(){
    ke = 0.0;
    for( unsigned int i = 0; i < npart; i++ ){
        for( int k = 0; k < dim; k++ ){
            ke += op[i].vel(k) * op[i].vel(k); 
        }
    }
    ke *= 0.5; //ke is a hardSphere class variable.
    return ke;
}



/*Calculate the average speed*/
template<typename real>
real hardSphere<real>::avgSpeed(){
    real as = 0.0;
    for( unsigned int i = 0; i < npart; i++ ){
        real s = 0.0;
        for( int k = 0; k < dim; k++ ){
            s += op[i].vel(k) * op[i].vel(k);
        }
        as += sqrt( s );
    }
    return as /( 1.0 * npart );
}



/*Calculate the mean Free Flight Time for the system of particles*/
template<typename real>
void hardSphere<real>::calcStatisticalQuantities( real packFrac, real temp ){
    real enskogCorr = (2 - packFrac) / ( 2 * pow( ( 1 - packFrac), 3) );
    //std::cout<<"\n enskogCorr: "<<enskogCorr;
    avgOmega1 = 24 * packFrac * enskogCorr * sqrt( temp / M_PI );
    avgFFT    = 1.0 / avgOmega1; //std::cout<<"\n avgFFT: "<<avgFFT;
    avgFPL    = 1.0 / ( 6 * sqrt(2.0) * packFrac * enskogCorr );
    avgSpd    = sqrt ( 8 * temp / M_PI ); //std::cout<<"\n avgSpdTh: "<<avgSpd;
    avgSpdSim = avgSpeed();               //std::cout<<"\n avgSpdSim: "<<avgSpdSim;
    avgNct    = ( 2 * avgFFT )/npart;
    //double avgSpdSim   = kineticEnergy() / ( 1.5 * npart ); //std::cout<<"\n rmsSpdSim: "<<avgSpdSim;
}



/*Sort the particles into different cells*/
template<typename real>
void hardSphere<real>::sortIntoCells(){
    int c[3], ci;
    //Clear the cell contents.
    for( int ci = 0; ci < nCell; ci++ ){
        //std::cout<<" before oc[ci].cell.size(): "<<oc[ci].cell.size();
        oc[ci].cell.clear();
        //std::cout<<" after oc[ci].cell.size(): "<<oc[ci].cell.size()<<"\n";
    }
    //usleep(3000000);


    //Get the cell index and sort the particles into cells.
    for( unsigned int i = 0; i < npart; i++ ){
        for( int k = 0; k < dim; k++ ){
            //std::cout<<"\n op[i].pos(k): "<<op[i].pos(k);
            //std::cout<<"\n cellLen: "<<cellLen;
            c[k] = op[i].pos(k) / cellLen;
        }
        //usleep(3000000); //In microseconds units.
        //std::cout<<"\n cellIndex: "<<c[0]<<"\t"<<c[1]<<"\t"<<c[2];
        ci = cellIndex( c[0], c[1], c[2] );
        //std::cout<<"\n ci: "<<ci;
        oc[ci].cell.push_back(i);
        //storeIntoCell(ci, vi);
    }
}



/*Get the cell index*/
template<typename real>
inline int hardSphere<real>::cellIndex( int cx, int cy, int cz ){
    return cx + cy*nCellx + cz * nCellx * nCellx;
}



template<typename real>
inline void hardSphere<real>::storeIntoCell( int ci, int i ){
    oc[ci].cell.push_back(i);
}


/*Update the position of all the particles*/
template<typename real>
void hardSphere<real>::updatePositions( real dt ){
    
    for( unsigned int i = 0; i < npart; i++ ){

        //double wt = getWallCollTime(i); //Get the yth direction wall collision time for the particle i.
        //std::cout<<"\n wt: "<<wt<<"\tdt: "<<dt;
        //if( wt > dt ){
            op[i].updatePos( dt ); //It is a function in the class particle.h, This is mostly the case.
            op[i].setPbc(boxLen);
        /*}
        else{
            op[i].updatePos( wt );
            real dv[3] = { 0.0, 0.0, 0.0 };
            boundaryConditionVelChange( i, dv );
            //std::cout<<"\n dv: "<<dv[0]<<"\t"<<dv[1]<<"\t"<<dv[2];
            op[i].resetVel( dv ); //Update the velocity.

            op[i].updatePos( (dt - wt) );
            op[i].setPbc( boxLen );
        }
        */
    }
}



/*Get wall collision time for particle i in y direction, top specular, bottom is thermal wall */
template<typename real>
real hardSphere<real>::getWallCollTime( int pi ){
    double dist; //The wall is in z direction.
    int k = 2;
    if( op[pi].vel(k) > 0 )
        dist = boxLen - op[pi].pos(k); //Moving towards the positive side.
    else 
        dist = op[pi].pos(k); //The co-ordinate of the bottom of the boxLen is 0. The position is always positive, check if it's not here.
    //std::cout<<"\n op[pi].vel(k): "<<op[pi].vel(k) ;
    //std::cout<<"\n abs( op[pi].vel(k) ): "<<abs( op[pi].vel(k) );
    return ( dist / abs( op[pi].vel(k) ) );
}



/*Implement the specular reflection from the top wall, and a gaussian change in veloicty from the bottom wall*/
template<typename real>
void hardSphere<real>::boundaryConditionVelChange( int pi, real* dv ){
    int k = 2; //Direction is z direction.
    if( op[pi].vel(k) > 0.0 ){
        //Specular collision.
        dv[k] = -op[pi].vel(k); 
        dv[0] = op[pi].vel(0);
        dv[1] = op[pi].vel(1);
    }
    else{
        
        //The Gaussian initialization to the velocity variable.
        dv[0] = sqrt( wallTemp ) * gaussianRv( 0, 1 ); //mean = 0, std = 1
        dv[1] = sqrt( wallTemp ) * gaussianRv( 0, 1 ) + wallVel;
        dv[2] = sqrt(-2 * wallTemp * log( drand48() ) ); //always posititve.
        //yVel.push_back( simTime );
        //yVel.push_back( dv[1] );
    }
}



/* Perform the collision */
template<typename real>
int hardSphere<real>::collision(){
    //Loop over all the cells, select a pair of particles, if the collision between them is accepted, call the pair collision operation here.
    int numCellParticle = 0, nCollP = 0, nCollE = 0;
    int i, j, pi, pj; 

    //std::cout<<"\n nCell: "<<nCell;
    for( int ci = 0; ci < nCell; ci++){
        //For each cell.
        numCellParticle = oc[ci].cell.size(); //std::cout<<"\n numElem: "<<numElem;
        nCollP = numOfCollToAttempt( numCellParticle );

        for( int np = 0; np < nCollP; np++ ){
            //For particles in each cell.
            i = numCellParticle * drand48(); 
            j = numCellParticle * drand48();
            //i = rand( numCellParticle );
            //j = rand( numCellParticle );
            pi = oc[ci].cell[i];
            pj = oc[ci].cell[j];
            //Check if particles are approaching towards each other.
            
            nCollE += pairCollision( pi, pj );
            //std::cout<<"\n i, j, pi, pj: "<<i<<"\t"<<j<<"\t"<<pi<<"\t"<<pj;
        }
    }
    return nCollE;
}



/*Number of collisions to attempt*/
template<typename real>
int hardSphere<real>::numOfCollToAttempt( int npa ){
    return ( 4 * M_PI * npa * npa * vmax * dt ) / ( cellLen*cellLen*cellLen );
}



/*Perform the pair collision between the two particles, Select a random angle, check for the collision to be executed and then do it.*/
template<typename real>
int hardSphere<real>::pairCollision( int pi, int pj ){
    real dvr[dim], ehat[dim];
    real vnorm = 0.0, vrij = 0.0;
    real drij2 = 0.0;
    real drx = 0.0;

    //Get the relative velocity vector.
    for( int k = 0; k < dim; k++ ){
        dvr[k] = op[pi].vel(k) - op[pj].vel(k);
        drx = op[pi].pos(k) - op[pj].pos(k);
        vrij += dvr[k] * drx;
        drij2 += drx * drx;
    }

    //Check if they are not approching then get out of the loop here.
    if( vrij >=0.0 ){
        //std::cout<<"\n vnorm: "<<vnorm;
        return 0;
    }
    //If not, it means the particles are approaching hence proceed from here onwards. Get the random line of alps now.
    
    real phi = 2 * M_PI * drand48(); //Azimuthal angle.
    real cosTheta = 2 * drand48() - 1; //The p.
    real sinTheta = sqrt( 1 - cosTheta * cosTheta );
    ehat[0] = sinTheta * cos(phi);
    ehat[1] = sinTheta * sin(phi);
    ehat[2] = cosTheta;

    //Get the dot product of line of alps and the relative velocity.
    vnorm = 0.0;
    for( int k = 0; k < dim; k++ ){
        vnorm += dvr[k] * ehat[k];
    }
    //Assigning probability for the collision between the pair of particles.
    if( fabs( vnorm) < vmax * drand48() ) return 0;

    //Collision successful.
    double h = ( 1 + eps ) * vnorm * 0.5; 
    real dvi[3] = { -1.0 * h * ehat[0], -1.0 * h * ehat[1], -1.0 * h * ehat[2] };
    real dvj[3] = { h * ehat[0], h * ehat[1], h * ehat[2] };
    op[pi].updateVel( dvi );
    op[pj].updateVel( dvj );
    ke -= ( 1- eps * eps ) * 0.25 * vnorm * vnorm; //Updating the kinetic energy here.

    //Collect the virial for the pressure.
    brtot += vrij/sqrt(drij2); 
    return 1;
}

/*Function to sample the y component of the velocity near the z = 0 plate. This needs a different boundary condition for the Rayleigh problem.*/
template<typename real>
void hardSphere<real>::sampleYvel(){
    int count = 1; real sum = 0.0;
    int y = 1;
    real stride = 0.5 * avgFFT; //std::cout<<"\nstride: "<<stride;
    unsigned int pi; 
    int ci, cz = 0; //cz si zero because the bottom cell is to be considered.
    //Sample the y component of velocity periodically at 0.5 avgFFT;
    double rem = remainder(simTime, stride );
    if( 0 < rem && rem < dt ){
        //std::cout<<"\n simTime mod stride: "<<remainder(simTime, stride );
        //Loop over the bottom cell to get the velocity.
        //std::cout<<"\n simTime from sample velocity: "<<simTime<<"\tdt: "<<dt;
        
        count = 0, sum = 0.0; //Resetting them to zero.
        for( int cx = 0 ; cx < nCellx; cx++ ){
            for( int cy = 0; cy < nCellx; cy++ ){
                ci = cellIndex( cx, cy, cz); //std::cout<<"\n ci: "<<ci;
                for( unsigned int p = 0; p < oc[ci].cell.size(); p++ ){
                    pi = oc[ci].cell[p];
                    //std::cout<<"\n pi: "<<pi<<"\toc[ci].cell.size(): "<<oc[ci].cell.size();
                    sum += op[pi].vel(y);
                    count++;
                }
            }
        }
        sum /= (1.0*count); //std::cout<<"\nY vel avg: "<<sum;
        yVel.push_back(sum);
        
    }
}



/*Function to save the average of y component of velocity*/
template<typename real>
void hardSphere<real>::saveYvel(){
    std::fstream fiO;
    char fileName[300];
    sprintf( fileName, "meanOf_velInYdirection_N%d_pf%.4lf_v%.4lf_eNormal%.2lf_cellLen%.2lf_maxTime%.1lf.dat", npart, packFrac, eps, cellLen, boxLen*boxLen*boxLen, maxTime );
    fiO.open( fileName, std::ios::app );
    fiO<< std::fixed << std::setprecision(16) << std::endl;
    if( trial == 0){
        fiO<<"#The first row is simulation time, the second row onwards are the average of the y componenet velocity of cells near z = 0 at the specified simulation time in the first row\n";
        for( unsigned int i = 0; i < simT.size(); i++ ){
            fiO<<simT[i]<<"\t";
        }
        fiO<<std::endl;
        for( unsigned int i = 0; i < yVel.size(); i++ ){
            fiO<<yVel[i]<<"\t";
        }
        fiO<<std::endl;
    }
    else{
        for( unsigned int i = 0; i < yVel.size(); i++ ){
            fiO<<yVel[i]<<"\t";
        }
        fiO<<std::endl;
    }
    fiO.close();
    //std::cout<<"\n yVel.size(): "<<yVel.size();
    //yVel.clear();
    //simT.clear();
}


/*Function to save the yth component of the velocity near z = 0 for any particle near the wall*/
/*
template <typename real>
void hardSphere<real>::saveYvel(){
    std::fstream fiO;
    char fileName[300];
    sprintf( fileName, "mean_velInYdirection_N%d_pf%.4lf_v%.4lf_eNormal%.2lf_cellLen%.2lf_maxTime%.1lf.dat", npart, packFrac, eps, cellLen, boxLen*boxLen*boxLen, maxTime );
*/
    /*
    double sum = std::accumulate( yVel.begin(), yVel.end(), 0.00);
    double mean = sum / yVel.size(); std::cout<<"\n mean: "<<mean;

    double sq_sum = std::inner_product( yVel.begin(), yVel.end(), yVel.begin(), 0.0 );
    double stdDev = std::sqrt( sq_sum / yVel.size() - mean * mean ); 
    std::cout<<"\n stdDev: "<<stdDev;
    */
    
 /*   

    //for( unsigned int i = 0; i < yVel.size(); i+=2 ){
    //    fiO<<yVel[i]<<"\t"<<yVel[i+1] << std::endl;
    //}
    real sum = 0.0, stdDev = 0.0;
    for( unsigned int i = 1; i < yVel.size(); i+=2 ){
        sum += yVel[i];
        stdDev += yVel[i] * yVel[i];
        //fiO<<yVel[i]<<"\t"<<yVel[i+1] << std::endl;
    }
    sum /= (0.5 * ( yVel.size() )); std::cout<<"\n yVel.size(): "<<yVel.size();
    stdDev /= (0.5 * ( yVel.size() ));std::cout<<"\nstdDev: "<<stdDev;
    stdDev -= (sum * sum); 
    stdDev = sqrt( stdDev );
    std::cout<<"\n"<<yVel[2]<<"\t"<<sum<<"\t"<<stdDev<<std::endl;



    fiO.open( fileName, std::ios::app );
    //std::vector<int>::iterator ptr;
    static int count = 0;
    if( count == 0){
        //fiO<<"#SimTime(t)   v_y( z = 0, t )"<<std::endl;
        fiO<<"#SimTime(t)   <v_y( z = 0, t )>    stdDev "<<std::endl;
        count++;
    }

    fiO<<yVel[2]<<"\t"<<sum<<"\t"<<stdDev<<std::endl;
    fiO.close();
    yVel.clear();

}
*/

/*Function to sample the observables*/
template<typename real>
void hardSphere<real>::sampleObservables(){
    real stride = 0.2 * avgFFT; //std::cout<<"\nstride: "<<stride;
    double rem = remainder(simTime, stride );
    if( 0 < rem && rem < dt ){
        //sampleYvel();
        samplePressure();
        //std::cout<<"\n Sampling pressure: ";
        //saveConfig(); //There is no sampling required, just save it periodically. Just a Jugaad.
        simT.push_back(simTime);
    }
}

/*Save the observables*/
template<typename real>
void hardSphere<real>::saveObservables(){
    //saveYvel();
    savePressure();
    saveConfig(); //It will save only the last configuration.
}



/*Sample the pressure from system*/
template<typename real>
void hardSphere<real>::samplePressure(){
    double v2 = 0.0; 
    for( unsigned int i = 0; i < npart; i++ ){
        for( int k = 0; k < dim; k++ ){
            v2 += op[i].vel(k) * op[i].vel(k);
        }
    }

    double pid = ( 1/(3.0*boxVolm) ) * v2; //Mass of the particle is taken to be equal to 1.
    //double pid = rho * temp; //Next time calcualte with the above mentioned procedure.
    double pvir =  -brtot / ( 3.0 * boxVolm * ( simTime - prevTime ) );
    double ptot = pid + pvir;
    prevTime = simTime;
    brtot = 0.0;
    pressArr.push_back(ptot);
    //std::cout<<"\n pid: "<<pid<<"\tpvir: "<<pvir<<"\tptot: "<<ptot;
}


/*Save the pressure into a data file*/
template<typename real>
void hardSphere<real>::savePressure(){
    std::fstream fiO;
    char fileName[300];
    sprintf( fileName, "pressure_N%d_pf%.4lf_eNormal%.2lf_cellLen%.2lf_maxTime%.1lf_simTime%.2lf.dat", npart, packFrac, eps, cellLen, maxTime, simTime - dt );
    fiO.open( fileName, std::ios::app );
    fiO<< std::fixed<<std::setprecision(16)<<std::endl;
    if( trial == 0){
        fiO<<"#The first row is simulation time, the second row onwards are the pressure at different time instances of the given trial"<<std::endl;
        for( unsigned int i = 0; i < simT.size(); i++ ){
            fiO<<simT[i]<<"\t";
        }
        fiO<<std::endl;
        for( unsigned int i = 0; i < pressArr.size(); i++ ){
            fiO<<pressArr[i]<<"\t";
        }
        fiO<<std::endl;
    }
    else{
        for( unsigned int i = 0; i < pressArr.size(); i++ ){
            fiO<<pressArr[i]<<"\t";
        }
        fiO<<std::endl;
    }
    fiO.close();
}

/*Save the configuration of particle configurations*/
template<typename real>
void hardSphere<real>::saveConfig(){
    //std::cout<<"\n Saving the configurations";
    std::fstream fiO;
    char fileName[300];
    sprintf( fileName, "config_N%d_pf%.4lf_eNormal%.2lf_cellLen%.2lf_maxTime%.1lf_simTime%.2lf.dat", npart, packFrac, eps, cellLen, maxTime, simTime-dt );
    fiO.open( fileName, std::ios::out );
    fiO<< std::fixed<<std::setprecision(16)<<std::endl;
    
    fiO<<"#npart, boxLx, boxLy, boxLz \n#x\ty\tz\tvx\tvy\tvz"<<std::endl;

    fiO<<npart<<std::endl;
    fiO<<boxLen<<"\t"<<boxLen<<"\t"<<boxLen<<std::endl;

    for( unsigned int i = 0; i < npart; i++ ){
        for( int k = 0; k < dim; k++ )
            fiO<<op[i].pos(k)<<"\t";
        
        for( int k = 0; k < dim; k++ )
            fiO<<op[i].vel(k)<<"\t";
        
        fiO<<std::endl;
    }
    fiO.close();
}
