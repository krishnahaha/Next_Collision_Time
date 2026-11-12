/*Aim: A software for Event Driven Molecular Dynamics simulation of Hard-Sphere, including event calender of binary tree, neighbourlist and cell list.
Compile: g++ main.cpp 
Run: ./a.out 'parameter.par'
*/
typedef double real;

#include<iostream>
#include<cstdlib>
#include<math.h>
#include<fstream>
#include <chrono>

#include "hardSphere.cpp"


int main( int argc, char *argv[] ){
    FILE *fpPara; struct parameter p;
    fpPara = fopen( argv[1], "r");
    if( fpPara == NULL ){
        std::cout<<"\n\n ############### Parameter file not found, simulation exit ###############"; 
        exit(0);
    }
    std::setprecision(12);
    readPara( fpPara, &p );

    hardSphere<real>hs( &p ); //Initializing the hardSphere system. Initializes particle positions, velocity and other system variables.
    hs.simTime = 0.0;
    
    auto start0 = std::chrono::system_clock::now();
    while( hs.simTime <= hs.maxTime ){
        hs.edmdSimStep();
        //hs.checkOverlap(); //Its a costly operation, use wisely.
    }
    auto end0 = std::chrono::system_clock::now();
    hs.cpuTime = std::chrono::duration_cast<std::chrono::milliseconds>(end0 - start0).count(); hs.cpuTime *= 0.0010;

    hs.systemReport();
    hs.saveObservables();   

    return 0;
}

