/*Aim: Creating a particle class for DSMC simulations. Include data which are very specific and needed for each particle.
*/

#ifndef PARTICLE_H_INCLUDED
#define PARTICLE_H_INCLUDED

#include<iostream>
#include<cstdlib>
#include<math.h>
#include<fstream>

#define dim 3

template<typename real>
class particle{
public:
    real data[6]; //Data for each particle. //Which are being used repeatedly there.

    real  pos( unsigned int k ) const{ return data[k]; } //Getter 
    real& pos( unsigned int k )      { return data[k]; } //Setter

    real  vel( unsigned int k ) const{ return data[dim+k]; } //Getter 
    real& vel( unsigned int k )      { return data[dim+k]; } //Setter
    
    inline void updatePos( real dt){
        for( int k = 0; k < dim; k++){
            pos(k) += dt * vel(k);
        }
        /*  
            if( k==2){
                if( pos(k) >= boxL ){
                    std::cout<<"\nIN Z DIRECTION greater than boxL pos(k): "<<pos(k)<<"\tboxL: "<<boxL;
                }
                else if( pos(k) < 0.0 )
                    std::cout<<"\nIN Z DIRECTION less than 0.0: pos(k): "<<pos(k)<<"\tboxL: "<<boxL;
            }
            //PBC
            if(pos(k) > boxL){
                pos(k) = fmod( pos(k), boxL ); //Setting up the PBC.
            }
            else if( pos(k) <= 0.0){
                pos(k) = boxL - fmod( fabs(pos(k) ), boxL);
            }
        }
        */
    }

    inline void setPbc( real boxL){
        for( int k = 0; k < dim; k++ ){
            if( pos(k) >= boxL){
                pos(k) = fmod( pos(k), boxL);
            }
            else if( pos(k) < 0.0){
                pos(k) = boxL - fabs( pos(k));
            }
            /*if( k == 2 ){
                if( pos(k) >= boxL || pos(k) <= 0.0 ){
                    std::cout<<"pos(k): "<<pos(k)<<"\tboxLen: "<<boxL;
                }
            }*/
        }
    }

    
    inline void updateVel( real* dv){
        for( int k = 0; k < dim; k++){
            vel(k) += dv[k];
        }
    }
    inline void resetVel( real* dv){
        for( int k = 0; k < dim; k++){
            vel(k) = dv[k];
        }
    }
    
};

#endif //PARTICLE_H_INCLUDED


