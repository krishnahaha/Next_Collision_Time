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

// Function to initialize position and velocity of hard particles randomly.
void hardSphere::initPosVel(struct PARAINPUT *p)  {
  //Initialize the first particle position.
  double tmp1=0;
  for(int k = 0; k<DIM; k++)  {
    posVel[k] = -0.5*boxLen + boxLen*drand48();
  }
  double rij2 = 0.0;  int cnt=0;
  for(int i = 1; i<npart; i++)  {
    rij2 = 0.0; cnt = 0;
    while(rij2<= 1 )  {
      for(int k=0; k<DIM; k++)  {
        posVel[i*STRIDE+k] = boxLen*(drand48() - 0.5);
      }
      for(int j = 0; j<i; j++)  {
        rij2 = 0.0;
        for (int k = 0; k<DIM; k++) {
          tmp1 = (posVel[j*STRIDE + k] - posVel[i*STRIDE+k]);
          tmp1 = tmp1 - boxLen*lround(tmp1*boxLenInv);
          rij2 += tmp1*tmp1;
        }
        if(rij2<=1) {
          break;
        }
      }
      cnt++;
      if(cnt>1000*npart)  {
        cout<<"\nPosition of i = "<<i<<" can't be specified after "<<cnt<<" steps####";
        exit(0);
      }
    }//while ends.
  }
  cout<<"\n Position initialized sucessfully!!!";

  //Initializing the velocities is tricky for the two and the three dimensions. The first dimesion is given random velocity, third dimension is adjusted to maintain a velocity of 1 with random direction.
  for (int i=0; i<npart; i++) {
    int ind=1; double tmp = 0, tmp1;
    float randS;
    while(ind==1) {
      tmp=0;
      for(int k=DIM; k<2*DIM-1; k++)  {
        posVel[i*STRIDE+k] = p->vInit*(2*drand48() - 1);
        tmp1 = posVel[i*STRIDE+k];
        tmp += tmp1*tmp1;
      }

      tmp = (p->vInit*p->vInit - tmp);
      if (tmp>=0) {
        ind = 0; //Makes the while loop terminates.
        posVel[i*STRIDE+2*DIM-1] = pow(tmp, 0.5);
        randS = (2*drand48() - 1);
        if (randS<0)  {
          posVel[i*STRIDE+2*DIM-1] = -posVel[i*STRIDE+2*DIM-1]; //randomly associate a velocity sign.
        }
      }
      else  {
        //cout<<"\n Invalid velocity encountered";
      }
    }//while ends.
  }//npart ends.
  cout<<"\n Velocity initialized!!!";
}//initPosVel ends.

void hardSphere::initPosVelScc(struct PARAINPUT *p) {
  int partCount = 0;
  double bSPp, smL, gridLen;
  bSPp = boxVolm/npart;
  if(bSPp < 1)  {
    cout<<"\n ### Caution overlapping while initializing in lattice ###";
  }
  smL = pow(bSPp, 0.3333333333333333333);
  cout<<"\nsmL: "<<smL;
  gridLen = pow(npart, 0.3333333333333333);
  for (int i=0; i< gridLen; i++)   {
    for (int j=0; j< gridLen; j++) {
      for (int k=0; k< gridLen; k++)   {
        posVel[partCount*STRIDE]     = (i+0.5)*smL - boxLen*0.5;
        posVel[partCount*STRIDE + 1] = (j+0.5)*smL - boxLen*0.5;
        posVel[partCount*STRIDE + 2] = (k+0.5)*smL - boxLen*0.5;
        partCount++;
      }//k ends
    }//j ends
  }//i ends.

    //Initializing the velocities is tricky for the two and the three dimensions. The first dimesion is given random velocity, third dimension is adjusted to maintain a velocity of 1 with random direction.
    for (int i=0; i<npart; i++) {
      int ind=1; double tmp = 0, tmp1;
      float randS;
      while(ind==1) {
        tmp=0;
        for(int k=DIM; k<2*DIM-1; k++)  {
          posVel[i*STRIDE+k] = p->vInit*(2*drand48() - 1);
          tmp1 = posVel[i*STRIDE+k];
          tmp += tmp1*tmp1;
        }

        tmp = (p->vInit*p->vInit - tmp);
        if (tmp>=0) {
          ind = 0; //Makes the while loop terminates.
          posVel[i*STRIDE+2*DIM-1] = pow(tmp, 0.5);
          randS = (2*drand48() - 1);
          if (randS<0)  {
            posVel[i*STRIDE+2*DIM-1] = -posVel[i*STRIDE+2*DIM-1]; //randomly associate a velocity sign.
          }
        }
        else  {
          //cout<<"\n Invalid velocity encountered";
        }
      }//while ends.
    }//npart ends.
    cout<<"\n Velocity initialized!!!";
}//initPosVelScc ends.
