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


/* Subrouting to initialize the simulation box.
INDEX
  ** initHsVar()
  ** initArrays()
  ** ~hardSphere()
  ** printPosVel()
  **
*/

//////////////////// initHsVar //////////////////
void hardSphere::initHsVar(struct PARAINPUT *p) {
  npart = p->npart; dim = p->dim; stride = 2*DIM;
  phi = p->phi; nCollSave = p->nCollSave;
  nCollIgnore = p->nCollIg;
  double expo;
  expo = 0.33333333333333333333;
  boxVolm = (p->npart*4*expo*M_PI*0.125/p->phi); cout<<"\n boxVolm = "<<boxVolm;
  cout<<"\n PI = "<<M_PI;
  boxLen = pow(boxVolm, expo); cout<<"\n boxLen = "<<boxLen;
  boxLenInv = 1.0/boxLen;
	cout<<"\n Allowed min and max positions for hard sphere center = "<<(0.5*boxLen);

  collCount = 0; null = 1e-10; infinity= 1e10;
  Time = 0.000000000000000; imgIndex = 0;

	//initPosVel(box, p); Call this in the Main.
  cout<<"\n Array initialization Successful!!!";
  /*rng = (pcg64_rng*)malloc(sizeof(pcg64_rng));
  if(seed_rng_with_entropy(rng) < 0)
  {
		//If the seeding function returns -1, then there is not enough entropy in the pool.
  	fprintf(stderr, "\nSeeding from /dev/urandom failed Exiting...\n");
		exit(0); // if so, do not continue!
  }*/
}

//////////////////// initArrays ////////////////
void hardSphere::initArrays(struct PARAINPUT *p) {
  //ARRAY INITIALIZATION.
  posVel = (double *)calloc(int(p->npart*DIM*2), sizeof(double));
  imgPos = (double *)calloc(IMGN*DIM, sizeof(double));
  edgeVec = (double *)calloc(DIM*DIM, sizeof(double));
  impactSeparation = (double *)calloc(p->nColl+1, sizeof(double));
  //In chuncks of nCollSave, the memory is allocated and the data is flushed after appending them in a file.
  kineticEne  = (float*)calloc(p->nColl+1, sizeof(float));
  minCollTime = (double *)calloc(p->nColl+1, sizeof(double));

  edgeVec[0] = boxLen; edgeVec[DIM+1] = boxLen;
  edgeVec[2*DIM+2] = boxLen;
  //vector<double> clist[npart];
  //vector<float>freeFlightTime[npart];
}

///////////////// ~hardSphere///////////////////
hardSphere::~hardSphere()
{
  free(posVel); free(imgPos); free(kineticEne);
  free(impactSeparation);
  free(minCollTime); free(edgeVec);
  //free(rng);
}
