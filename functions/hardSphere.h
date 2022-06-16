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


/*Create a class for the hardSphere object having functions to perform the edmd simulation and the data variables to store the description of the system.*/

class hardSphere
{
public:
  //Reading input data and initializing the class variables and printing.
  void readPara(FILE *fp, struct PARAINPUT *p);
  void initArrays(struct PARAINPUT *p);
  void initHsVar(struct PARAINPUT *p);
  void initPosVel(struct PARAINPUT *p);
  void initPosVelScc(struct PARAINPUT *p);

  //System properties
  void printPosVel();
  float calcKineticEne();
  void printCollTable(int);
  void printEfficiencyReport();

  //edmd functions.
  void particlePropagate(double delT);
  void particleToWallPropagate(int, int, int, double delT);
  tuple<double, int, int, int> pwCollTime(int, double globTime);
  void imageFunction(int);
  tuple<double, int, int, int>ppCollTime(int, int, double);
  void collTableInit(double globTime);
  void velUpdate(int, int, int);
  void collTableUpdate(int, double);
  int checkOverlap();
  double calcImpactSeparation(int, int);
  void edmd();

  //Data save functions.
  void savePosVel();
  void saveMinCollTime();
  void saveFreeFlightTime();
  void saveCollTable();
  void saveKineticEne();
  void saveImpactSeparation();


  //Class destructor.
  ~hardSphere();

  //Public variables for the class.
  int nInit=0, nColl=0;
  //float *kineticEne, *minCollTime;
private:
  int npart, dim, stride, nCollSave, nCollIgnore;
  double phi;
  double *posVel, *imgPos, *edgeVec, *impactSeparation, *minCollTime;
  float *kineticEne;
  double boxVolm, boxLen, boxLenInv, null, infinity, Time;
  int collCount, imgIndex;

  //Maps and container classes.
  map<double, tuple<int, int, int>>wallColl;
  tuple<double, int, int, int>pwReturn;
  tuple<int, int, int> ijnext;
  tuple<double, int, int, int> ppReturn;
  map <double, tuple<int, int, int> > imgCseq;
  map< double, tuple<int,int, int> >cseq;

  //vectors
  vector<double> clist[100000];
  vector<float>freeFlightTime[100000];

  //Efficiency test variables
  double particlePropagateT=0, pwCollTimeT=0, particleToWallPropagateT=0, collTableUpdateT=0, ppCollTimeT=0, velUpdateT=0, calcKineticEneT=0;
};
