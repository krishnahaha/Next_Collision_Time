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


/*Definitions of the function for the hardSphere.h header file which are used in doing the edmd simulation for a set of hardSphere particles
Index:
  ** printEfficiencyReport()
  ** calcImpactSeparations()
  ** void particlePropagate()
  ** void particleToWallPropagate()
  ** pwCollTime()
  ** imageFunction()
  ** ppCollTime()
  ** void collTableInit()
  ** void velUpdate()
  ** void collTableUpdate()
  ** int checkOverlap()
  ** void edmd()
*/

void hardSphere::printEfficiencyReport()  {
  cout<<"\n#### ALL TIME IN MICROSECONDS #####";
  cout<<"\nparticlePropagateT: \t"<<particlePropagateT;
  cout<<"\nParticleToWallPropagateT:\t"<<particleToWallPropagateT;
  cout<<"\npwCollTimeT:\t"<<pwCollTimeT;
  cout<<"\nppCollTimeT:\t"<<ppCollTimeT;
  cout<<"\nppCollTime includes the time to the first initialization of the collision table";
  cout<<"\ncollTableUpdateT:\t"<<collTableUpdateT;
  cout<<"\nvelUpdateT:\t"<<velUpdateT;
  cout<<"\ncalcKineticEneT:\t"<<calcKineticEneT;
}

////////////// printCollTable /////////////////
void hardSphere::printCollTable(int limit) {
  float tmp21;
  int tmp22, tmp23, tmp24, tmp25;
  map< double, tuple<int, int, int> >::iterator it = cseq.begin();
  cout<<"\nSr.N0.\tCollTime\t i\t j\t imgJ\n";
  tmp25 = 1;
  while( it != cseq.end())  {
    tmp21 = it->first;
    tie(tmp22, tmp23, tmp24) = it->second;
    cout<<"\n"<<tmp25<<"\t"<<tmp21<<"\t\t"<<tmp22<<"\t"<<tmp23<<"\t"<<tmp24;
    tmp25++;
    if(tmp25>limit) {
      break; //come out of printing the collision table.
    }
    it++;
  }
}

/////////////// impactSeparations //////////////
double hardSphere::calcImpactSeparation(int i1, int j1) {
  double tmp20, dis2;
  tmp20 = 0.0, dis2=0.0;
  for(int k=0; k<DIM; k++)  {
    //Considering only the nearest neighbour distances.
    tmp20 = posVel[i1*STRIDE + k] - posVel[j1*STRIDE + k];
    tmp20 = tmp20 - boxLen*lround(tmp20*boxLenInv);
    dis2 += tmp20*tmp20;
  }
  return dis2;
}

/////////////////// particlePropagate /////////////////////
void hardSphere::particlePropagate(double delT) {
  auto start1 = std::chrono::system_clock::now();
  for(int i = 0; i<npart; i++)  {
    for(int k = 0; k<DIM; k++)  {
      posVel[i*STRIDE + k ] += delT*posVel[(i*STRIDE) + (DIM+k)]; //Velocity is stored after the positions.
    }//k ends.
  }//i ends.

  auto end1 = std::chrono::system_clock::now();
  particlePropagateT += std::chrono::duration_cast<std::chrono::microseconds>(end1 - start1).count();
}

////////////// particleToWallPropagate /////////
void hardSphere::particleToWallPropagate(int i13, int i14, int i15, double delT)
{
  //i13: Particle index; i14: wall index; i15: nearest/farthest wall index there.
  auto start2 = std::chrono::system_clock::now();
  for(int i=0; i<npart; i++)  {
    //because the particle crossing the box needs come out of the box from the next side, so...
    if(i!=i13) {
      for(int k = 0; k<DIM; k++)  {
        posVel[i*STRIDE + k ] += delT*posVel[(i*STRIDE) + (DIM+k)]; //Velocity is stored after the positions.
      }//k ends.
    }
    else  {
    //bring the particle to from the other side of the box there and then see what happens there. How to do this...
      for(int k=0; k<DIM; k++)  {
        if(k!=i14)  {
          posVel[i*STRIDE + k ] += delT*posVel[(i*STRIDE) + (DIM+k)];
        }
        else  {
          //If the particle is moving forward it will come from the backward to the negative side.
          if (posVel[(i*STRIDE) + (DIM+k)]>0) {
            posVel[i*STRIDE + k ] = -0.5*boxLen;
          }
          else  {
            posVel[i*STRIDE + k ] = 0.5*boxLen;
          }
        }
      }//k ends.
    }//else ends.
  }//i ends.

  auto end2 = std::chrono::system_clock::now();
  particleToWallPropagateT += std::chrono::duration_cast<std::chrono::microseconds>(end2 - start2).count();
}//particleToWallPropagate ends.

//////////////////// particleWallCollision ////////////////

tuple<double, int, int, int> hardSphere::pwCollTime(int i, double globTime) {
  double tmp6=0.0; int i10, i11, i12;//Temporary Var

  auto start3 = std::chrono::system_clock::now();
  wallColl.clear();
  for(int k=0; k<DIM; k++)  {
    //Approaching to the nearest wall.
    if((posVel[i*STRIDE + (DIM+k)]) > 0)  {
      tmp6 = (globTime + (0.5*boxLen - posVel[i*STRIDE + k])/posVel[i*STRIDE + (DIM +k)]);
      wallColl[tmp6] = tuple<int, int, int>(i,k,-1);
      //Particle, in which dimension, with the nearest wall there.
    }
    else if((posVel[i*STRIDE + (DIM+k)]) < 0) {
      //Approaching farthest wall, the distances gets added up there.
      tmp6 = (globTime + (-0.5*boxLen - posVel[i*STRIDE + k])/posVel[i*STRIDE + (DIM +k)]);
      wallColl[tmp6] = tuple<int, int, int>(i,k,-2);
      //Particle, in which dimension, with the next wall there.
    }
    else
    {
      tmp6 = infinity;
      wallColl[tmp6] = tuple<int, int, int>(i,k,-3);
    }
  }//k ends
  //Get the time of the minimum collision.
  tmp6 = wallColl.begin()->first;
  tie(i10, i11, i12) = wallColl.begin()->second;
  //cout<<"\ntmp6 from pwColl: "<<tmp6;
  pwReturn = make_tuple(tmp6, i10, i11, i12);

  auto end3 = std::chrono::system_clock::now();
  pwCollTimeT += std::chrono::duration_cast<std::chrono::microseconds>(end3 - start3).count();
  return pwReturn;
}

////////////////// particleParticleColl ///////////////////

tuple<double, int, int, int> hardSphere::ppCollTime(int i, int j, double globTime)  {
  double rij2, vij2, bij, xrij, yrij, zrij, xvij, yvij, zvij;
  double tmp, tmp1, tmp2, tmptmp=infinity;
  int imgCnt = 0, i3 = 0;

  auto start4 = std::chrono::system_clock::now();

  //Calculate the collision time for original position and the minimum image in all the three directions.
  for (int ii=0; ii<2; ii++)  {
    xrij = posVel[i*STRIDE + 0] - posVel[j*STRIDE + 0];
    xvij = posVel[i*STRIDE + DIM + 0] - posVel[j*STRIDE + DIM + 0];
    if(ii==1)
      xrij = xrij - boxLen*lround(xrij*boxLenInv);

    for(int jj=0; jj<2; jj++) {
      yrij = posVel[i*STRIDE + 1] - posVel[j*STRIDE + 1];
      yvij = posVel[i*STRIDE + DIM+ 1] - posVel[j*STRIDE + DIM + 1];
      if(jj==1)
        yrij = yrij - boxLen*lround(yrij*boxLenInv);

      for(int kk=0; kk<2; kk++) {
        zrij = posVel[i*STRIDE + 2] - posVel[j*STRIDE + 2];
        zvij = posVel[i*STRIDE + DIM + 2] - posVel[j*STRIDE + DIM + 2];
        if(kk==1)
          zrij = zrij - boxLen*lround(zrij*boxLenInv);

        bij  = xvij*xrij + yvij*yrij + zvij*zrij;
        rij2 = xrij*xrij + yrij*yrij + zrij*zrij;
        vij2 = xvij*xvij + yvij*yvij + zvij*zvij;

        //Check if they are approaching.
        if(bij >= 0) {
          tmp = infinity;
        }
        else if( rij2 <= 1-null) {
          //check if overlapping.Make them to exchange their velocities by assigining the global Time as a time of collision.
          tmp = globTime;
          cout<<"\n\nOverlapping particles\ti= "<<i<<"\tj= "<<j<<"\tOvp dist: "<<(1-sqrt(rij2));
        }
        else if( (bij*bij - vij2*(rij2 - 1)) < 0 )  {
          //Check if the discriminant is negative.
          tmp = infinity;
        }
        else  {
          //Calculate the smallest time for the two particles to have separation sigma.
          tmp1 = rij2 - 1;
          tmp2 = -bij + sqrt(bij*bij - vij2*tmp1);
          tmp = globTime + (tmp1/tmp2);
        }

        //Comparing the collision times.
        if(tmp<tmptmp)  {
          tmptmp = tmp;
          i3 = imgCnt;
        }
        imgCnt++;
      }//kk ends.
    }//jj ends.
  }//ii ends.

  //Make a tuple out of these three and then return the value there.
  ppReturn = make_tuple(tmptmp, i, j, i3 );
  auto end4 = std::chrono::system_clock::now();
  ppCollTimeT += std::chrono::duration_cast<std::chrono::microseconds>(end4 - start4).count();
  return ppReturn;
}

/////////////// collTableInit() /////////////////

void hardSphere::collTableInit(double globTime) {
  double tmp3, tmp4, tmp5;
  int i4, i5, i6;

  //Clear the content of the vector and the map container list. This is initialied after every 10^5 collisions. Keep the global time list here.
  cseq.clear();
  for (int i=1; i<npart; i++) {
    clist[i].clear();
    for(int j=0; j<i; j++)  {
      tie(tmp4, i4, i5, i6) = ppCollTime(i, j, globTime);
      //Storing it in a form of the sorted array.
      cseq[tmp4] = tuple<int, int, int>(i4, i5, i6);
      //store the collision times if the collision time is less than the infinity time there.
      if(tmp4<infinity) {
        clist[i].push_back(tmp4);
        clist[j].push_back(tmp4);
      }
    }
  }
  //Determinig the collision with the wall time there.
  for(int i=0; i<npart; i++)
  {
    tie(tmp4, i4, i5, i6) = pwCollTime(i, globTime);
    cseq[tmp4] = tuple<int, int, int>(i4, i5, i6);
    clist[i].push_back(tmp4);
  }
  tmp5 = cseq.begin()->first;
  ijnext = cseq.begin()->second;

  cout<<"\n\nThe first collision is in between, i: "<<get<0>(ijnext)<<"\tj: "<<get<1>(ijnext)<<"\tWith image of J index :"<<get<2>(ijnext) <<"\tAt: "<<tmp5;
}//collTableInit ends.

//////////////// velUpdate //////////////////////
void hardSphere::velUpdate(int i, int j, int imJ) {
  double rij2, vij2, bij, xrij, yrij, zrij, xvij, yvij, zvij;
  int imx, imy, imz, imtmp;

  auto start5 = std::chrono::system_clock::now();
  imz = imJ%2; imtmp = imJ/2;
  imy = imtmp%2; imtmp = imtmp/2;
  imx = imtmp%2;
  //cout<<"\nImg index: "<<imJ;
  //cout<<"\timz, imy, imx: "<<imz<<"\t"<<imy<<"\t"<<imx;

  //Velocity update between two particles i and j.
  xrij = posVel[i*STRIDE + 0] - posVel[j*STRIDE + 0];
  xvij = posVel[i*STRIDE + DIM + 0] - posVel[j*STRIDE + DIM + 0];
  if(imx==1)
    xrij = xrij - boxLen*lround(xrij*boxLenInv);

  yrij = posVel[i*STRIDE + 1] - posVel[j*STRIDE + 1];
  yvij = posVel[i*STRIDE + DIM + 1] - posVel[j*STRIDE + DIM + 1];
  if(imy==1)
    yrij = yrij - boxLen*lround(yrij*boxLenInv);

  zrij = posVel[i*STRIDE + 2] - posVel[j*STRIDE + 2];
  zvij = posVel[i*STRIDE + DIM + 2] - posVel[j*STRIDE + DIM + 2];
  if(imz==1)
    zrij = zrij - boxLen*lround(zrij*boxLenInv);

  bij  = xvij*xrij + yvij*yrij + zvij*zrij;
  rij2 = xrij*xrij + yrij*yrij + zrij*zrij;
  double temp1= bij/rij2;

  posVel[i*STRIDE + DIM + 0] -= (xrij*temp1);
  posVel[j*STRIDE + DIM + 0] +=  (xrij*temp1);
  posVel[i*STRIDE + DIM + 1] -= (yrij*temp1);
  posVel[j*STRIDE + DIM + 1] +=  (yrij*temp1);
  posVel[i*STRIDE + DIM + 2] -= (zrij*temp1);
  posVel[j*STRIDE + DIM + 2] +=  (zrij*temp1);

  auto end5 = std::chrono::system_clock::now();
  velUpdateT += std::chrono::duration_cast<std::chrono::microseconds>(end5 - start5).count() ;
}//vel update ends.

////////////////// collTableUpdate ///////////////

double ct=0.0;
int i7, i8, i9;

void hardSphere::collTableUpdate(int i, double globTime)  {
  auto start6 = std::chrono::system_clock::now();
  for(int ii=0; ii!=clist[i].size(); ii++)  {
    cseq.erase(clist[i][ii]); //clist[i][ii] is the key for the cseq there.
  }
  cseq.erase(infinity);
  clist[i].clear();
  for(int j = 0; j<npart; j++)  {
    if( j!=i)
    {
      tie(ct, i7, i8, i9) = ppCollTime(i, j, globTime); //Calculate the particle-particle collision time.
      if(ct < infinity)  {
        //Store it in the ordered list and the each particle collision array. ct is the key here.
        cseq[ct] = tuple<int, int, int>(i7, i8, i9);
        clist[i].push_back(ct);
        clist[j].push_back(ct);
      }//ct<box->infty ends.
    }//i!=j ends.
  }//j loop ends.
  tie(ct, i7, i8, i9) = pwCollTime(i, globTime);
  cseq[ct] = tuple<int, int, int>(i7, i8, i9);
  clist[i].push_back(ct);
  auto end6 = std::chrono::system_clock::now();
  collTableUpdateT+= std::chrono::duration_cast<std::chrono::microseconds>(end6 - start6).count();
}//collTableUpdate ends.

/////////////////// checkOverlap ///////////////

int hardSphere::checkOverlap()  {
  double dist2, tmp;
  int tmp19=0, ovp=0;
  //Overlap between the particles.
  for(int ii = 1; ii<npart; ii++)  {
    for( int jj = 0; jj<ii; jj++)  {
      dist2 = 0;
      for(int kk=0; kk<DIM; kk++)  {
        tmp19 = ( posVel[ii*STRIDE+kk]-posVel[jj*STRIDE+kk] );
        tmp19 = tmp - boxLen*lround(tmp*boxLenInv);
        dist2 += (tmp19*tmp19);
      }
      if(dist2<1)  {
        cout<<"\n### Particle Overlap ### i = "<<ii<<" j = "<<jj <<" Overlap= "<<dist2;
        tmp=1;
        ovp = 1;
      }
    }//j ends.
  }//i ends.

  //Overlap with the wall.
  for(int ii = 1; ii<npart; ii++)  {
    for( int kk=0; kk<DIM; kk++ )  {
      if( posVel[ii*STRIDE+kk] < -0.5*boxLen || posVel[ii*STRIDE+kk]  > 0.5*boxLen )  {
        cout<<"\n ### Wall overlap ### i= "<<ii<<"\tposVel: "<<posVel[ii*STRIDE+kk]<<"\tboxLen: "<<boxLen;
        ovp=1;
      }
    }//k ends.
  }//i ends.
  return ovp;
}//checkOverlap ends.

double globT0 = 0.0, globT = 0; //These are needed to be a global variable.
////////////////// edmd() //////////////////////
void hardSphere::edmd()
{
  double tmp7=0.0; int i13, i14, i15;
  double tmp8, tmp9, tmp10, dt=0.0;

  tmp7 = cseq.begin()->first;
  dt = tmp7 - globT;
  tie(i13, i14, i15) = cseq.begin()->second;
  minCollTime[0] = 0; //Initialization of the array
  //If wall collision.
  if(i15<0)
  {
    //If one particle crosses the box, bring it back from the other side, previously calculated time of collisions becomes invalid, recalculate the collision times for this particle.
    tie(tmp8, i13, i14, i15) = pwCollTime(i13, globT);
    dt = tmp8 - globT;
    globT = tmp8;

    particleToWallPropagate(i13, i14, i15, dt);
    //cout<<"\nUpdating time of collision...";
    collTableUpdate(i13, globT);
  }
  else {
    //if particle collision recalculate the time once again.
    tie(tmp9, i13, i14, i15) = ppCollTime(i13, i14, globT);
    dt = tmp9 - globT;
    globT = tmp9;
    particlePropagate(dt);
    velUpdate(i13, i14, i15);

    //Update the collision table.
    collTableUpdate(i13, globT);
    collTableUpdate(i14, globT);
    nColl++;

    /////////////// CALCULATING OBSERVABES //////////////

    //kineticEne[nColl] = calcKineticEne();
    //impactSeparation[nColl] = calcImpactSeparation(i13, i14);
    minCollTime[nColl] = globT + globT0;
    /*if(nColl >= nCollIgnore)  {
      freeFlightTime[i13].push_back(globT + globT0);
      freeFlightTime[i14].push_back(globT + globT0);
    }
    */
  }
  if(globT>infinity)  {
    //Reinitializing collision table, as precautionary measure.
    globT0 += globT;
    globT  = 0;
    collTableInit(0.0);
    cout<<"\nGlobal time reinitialized, globT0: "<<globT0;
  }
}//edmd ends.
