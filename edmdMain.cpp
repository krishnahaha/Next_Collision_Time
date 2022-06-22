/*********************************************************************************************
 *   Copyright (c) <2022>, <Santosh Ansumali@JNCASR> <Srikanth Sastry@JNCASR>                *
 *   All rights reserved.                                                                    *
 *   Redistribution and use in source and binary forms, with or without modification, are    *
 *   permitted provided that the following conditions are met:                               *
 *                                                                                           *
 *    1. Redistributions of source code must retain the above copyright notice, this list of *
 *       conditions and the following disclaimer.                                            *
 *    2. Redistributions in binary form must reproduce the above copyright notice, this list *
 *       of conditions and the following disclaimer in the documentation and/or other        *
 *       materials provided with the distribution.                                           *
 *    3. Neither the name of the <JNCASR> nor the names of its contributors may be used to   *
 *       endorse or promote products derived from this software without specific prior       *
 *       written permission.                                                                 *
 *                                                                                           *
 *       THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND     *
 *       ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED       *
 *       WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED.  *
 *       IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT,    *
 *       INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING,      *
 *       BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,       *
 *       DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF     *
 *       LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE     *
 *       OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED   *
 *       OF THE POSSIBILITY OF SUCH DAMAGE.                                                  *
 *                                                                                           *
 *       Contributors:         Krishna K Tiwari, Santosh Ansumali, Srikanth Sastry           *
 *       Suggestions:          ansumali@jncasr.ac.in, sastry@jncasr.ac.in                    *
 *       Bugs:                 ansumali@jncasr.ac.in, sastry@jncasr.ac.in                    *
 *                                                                                           *
 *********************************************************************************************/


//Include standard libraries
#define DIM 3
#define STRIDE 6
#define IMGN 27

#include<stdio.h>
#include<string.h>
#include <cstdio>
#include<iostream>
#include<cmath>
#include<fstream>
#include<math.h>
#define _USE_MATH_DEFINES
#include<time.h>
#include<iterator>
#include<map>
#include<vector>
#include<tuple>
#include <chrono>
using namespace std;

//#include "functions/pcg64.h" //random number gen.
//pcg64_rng *rng; seed_rng_with_entropy(&rng);
#include "functions/structure.h"
#include "functions/hardSphere.h"
#include "functions/readPara.h"
#include "functions/initArrays.h"
#include "functions/initPosVel.h"
#include "functions/eneForce.h"
#include "functions/saveData.h"
#include "functions/hardSphere.cpp"


int main(int argc, char *argv[])
{
  cout<<"\n Edmd begins...\n";

  //Create the class object for the hardSphere, a struct object for the parameter file and a FILE object for the input.
  hardSphere hs; struct PARAINPUT p; FILE *fpPara;
  srand48(time(0));
  double saveTime=0.0;
  fpPara = fopen(argv[1], "r");
  if(fpPara==NULL)
  { cout<<"\n### Parameter file not found ###"; exit(0);  }

  //Else read the parameter file and perform initialization.
  hs.readPara(fpPara, &p);
  hs.initHsVar(&p);
  hs.initArrays(&p);
  hs.initPosVel(&p);
  //hs.initPosVelScc(&p); //This is for initiliaing in high packing fraction on a cubic lattice.
  //hs.printPosVel();
  //hs.savePosVel();

  auto start0 = std::chrono::system_clock::now();

  hs.collTableInit(0.0);
  //hs.printCollTable(10);
  //cout<<"\n Initialization is sucessful!!!";
  //Start the event driven molecular dynamics.
  cout<<"\nInitial pos vel";
  hs.printPosVel();
  while(hs.nColl<p.nColl)
  {
    //auto start = std::chrono::system_clock::now();
    hs.edmd();
    /*//Periodic save of observables.
    if(hs.nColl % p.nCollSave == 0)
    {
      auto starts = std::chrono::system_clock::now();

      hs.savePosVel(); hs.saveMinCollTime();
      hs.saveFreeFlightTime(); hs.saveKineticEne();
      hs.saveImpactSeparation();
      hs.saveCollTable();

      auto ends = std::chrono::system_clock::now();

      saveTime += std::chrono::duration_cast<std::chrono::microseconds>(ends - starts).count();
    }
   */
  }
  auto end0 = std::chrono::system_clock::now();

  //hs.printCollTable(10);
  hs.printEfficiencyReport();
  cout <<"\n\nTotal run time: "<<std::chrono::duration_cast<std::chrono::microseconds>(end0 - start0).count() << " mu_s" <<endl;

  cout<<"\nTime spend in saving the files:"<<saveTime<<"mu_s";

  //FINAL SAVE OF OBSERVABLES.
  hs.savePosVel(); hs.saveMinCollTime();
  hs.saveFreeFlightTime(); hs.saveKineticEne();
  hs.saveImpactSeparation();
  hs.saveCollTable();

  cout<<"\n\n!!! Congratulations EDMD successfull !!!\n\n";
  return 0;
}//main ends.
