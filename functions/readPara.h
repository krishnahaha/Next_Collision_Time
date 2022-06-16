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


//Subroutine to read the input parameters from script

void hardSphere::readPara(FILE *fp, struct PARAINPUT *p)  {
  char string[300]; char *endptr;
  //npart
  fscanf(fp, "%s", string);fscanf(fp, "%s", string);
  p->npart = atoi(string); cout<<"\n npart = "<<p->npart;

  //dim
  fscanf(fp, "%s", string);fscanf(fp, "%s", string);
  p->dim = atoi(string); cout<<"\n dim = "<<p->dim;

  //nColl
  fscanf(fp, "%s", string);fscanf(fp, "%s", string);
  p->nColl = atoi(string); cout<<"\n nColl = "<<p->nColl;

  //nInit
  fscanf(fp, "%s", string);fscanf(fp, "%s", string);
  p->nInit = atoi(string); cout<<"\n nInit = "<<p->nInit;

  //phi
  fscanf(fp, "%s", string);fscanf(fp, "%s", string);
  p->phi = strtod(string, &endptr); cout<<"\n phi = "<<p->phi;

  //vInit
  fscanf(fp, "%s", string);fscanf(fp, "%s", string);
  p->vInit = strtod(string, &endptr); cout<<"\n vInit = "<<p->vInit;

  //nCollSave
  fscanf(fp, "%s", string);fscanf(fp, "%s", string);
  p->nCollSave = atoi(string); cout<<"\n nCollSave = "<<p->nCollSave;

  //nInit
  fscanf(fp, "%s", string);fscanf(fp, "%s", string);
  p->nCollIg = atoi(string); cout<<"\n nCollIg = "<<p->nCollIg;

cout<<"\n Parameters read successfully !!!";
}
