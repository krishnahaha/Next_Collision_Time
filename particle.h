/*Aim: Creating a particle class for DSMC simulations. Include data which are very specific and needed for each particle.
*/

#ifndef PARTICLE_H_INCLUDED
#define PARTICLE_H_INCLUDED

#include<iostream>
#include <cstdint>
#include<cstdlib>
#include<math.h>
#include<fstream>
#include <vector>

#define DIM 3

template<typename real>
class particle{
public:
    real data[2*DIM]; //Data for each particle. //Which are being used repeatedly there.
    real rneig[DIM]; //To sore center of neighbourlist.
    class particle *prev, *next; //Each pointer occupies 8 bytes.
    class particle *left, *right, *parent; //For Event calender.
    class particle *pj;
    //std::vector<particle<real> > *neighbors;
    real eventTime, tlast;
    int index, counter, counterPj, queue, cell;
    int boxesTravelx[DIM], eventType;
    uint8_t nearBoxEdge;

    real  pos( unsigned int k ) const{ return data[2*k]; } //Getter 
    real& pos( unsigned int k )      { return data[2*k]; } //Setter

    real  vel( unsigned int k ) const{ return data[2*k+1]; } //Getter 
    real& vel( unsigned int k )      { return data[2*k+1]; } //Setter

    inline void updatePos( real dt){
        for( int k = 0; k < DIM; k++){
            //std::cout<<"\npos(k): "<<pos(k)<<"\tvel(k): "<<vel(k)<<"\t dt: "<<dt;
            pos(k) += dt * vel(k);
            //std::cout<<"\tupdated: "<<pos(k);
        }
    }

    inline void setPbc( real boxL ){
        for( int k = 0; k < DIM; k++ ){
            if( pos(k) >= boxL ){
                pos(k) -= boxL;
                boxesTravelx[k] += 1; 
                //std::cout<<"\n boxesTravelx[k]: "<<boxesTravelx[k];
            }
            else if( pos(k) < 0.0 ){
                pos(k) += boxL;
                boxesTravelx[k] -= 1;
                //std::cout<<"\n boxesTravelx[k]: "<<boxesTravelx[k];
            }
        }
    }

    inline void updateVel( real* dv){
        for( int k = 0; k < DIM; k++){
            vel(k) += dv[k];
        }
    }
    
    inline void resetVel( real* dv){
        for( int k = 0; k < DIM; k++){
            vel(k) = dv[k];
        }
    }
};

#endif //PARTICLE_H_INCLUDED


