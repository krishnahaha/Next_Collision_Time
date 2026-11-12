/*Aim: A software for Direct Simulation Monte Carlo (DSMC) simulation for Rayleigh problem. The thermal and specular walls are in z-direction. The thermal wall at bottom is moving with some speed vy; in y-direction.

The following macroscopic properties can be obtained from this code.

    ** Fluid velcoity of gas near the wall as a function of time.
    ** Pressure
    ** 
Compile: g++7 main.cpp 
Run: ./a.out
*/
typedef double real;
#include<iostream>
#include<cstdlib>
#include<math.h>
#include<fstream>

#include "hardSphere.cpp"


int main( int argc, char *argv[] ){
    FILE *fpPara; struct parameter p;
    fpPara = fopen( argv[1], "r");
    if( fpPara == NULL ){
        std::cout<<"\n\n ############### Parameter file not found, simulation exit ###############"; 
        exit(0);
    }
    readPara( fpPara, &p );

    hardSphere<real>hs( &p ); //Initializing the hardSphere system. Initializes particle positions, velocity and other system variables.
    hs.simTime = 0.0;
    //std::cout<<"\n Before while loop";
    while( hs.simTime <= hs.maxTime ){
        hs.dt = hs.getdt(); //generate from the exponential distribution.
        //std::cout<<"\n hs.dt: "<<hs.dt;
        hs.simTime += hs.dt; //std::cout<<"\n"<<hs.simTime;
        hs.updatePositions( hs.dt );
        hs.sortIntoCells();
        hs.nColls += hs.collision();
        //hs.sampleYvel();
        hs.sampleObservables();
        //std::cout<<"\n In while loop";
        
        
    }
    //hs.saveYvel();
    hs.saveObservables();
    //std::cout<<"\n Particles propagated for: "<<hs.simTime<<" s";
    //hs.displayPosVel(10);
    //std::cout<<"\n Number of collisions executed: "<<hs.nColls;

    //std::cout<<"\n\n"; 
    return 0;
}

