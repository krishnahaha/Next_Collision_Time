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


/*Functions to save the data generated during simulation

INDEX
** saveCollTable()
** savePosVel()
** saveMinCollTime()
** saveFreeFlightTime()
** saveKineticEne()
** saveImpactSeparations()

*/
/////////////// saveCollTable /////////////////

void hardSphere::saveCollTable() {
  FILE *dataFile;
  char str[300] = "collTable";
  char ch[400];
  strcat(str,"_phi_" ); sprintf(ch, "%0.4f", phi); strcat(str, ch);
  strcat(str, "_npart_"); sprintf(ch, "%d", npart); strcat(str, ch);
  strcat(str, "_nInit_"); sprintf(ch, "%d", nInit); strcat(str, ch);
  strcat(str, "_nColl_"); sprintf(ch, "%d", nColl); strcat(str, ch);
  strcat(str, ".dat");
  dataFile = fopen(str, "w");
  if(!dataFile)  {
    cout<<"\n File not created###"; exit(0);
  }
  else  {
    map< double, tuple<int, int, int> >::iterator it1 = cseq.begin();
    double tmp26;
    int tmp27, tmp28, tmp29, tmp30=1;
    fprintf(dataFile, "#ppCollStep \t tcol \t\t i \t\t j \t\t imgJ \n");
    while( it1 != cseq.end())
    {
      tmp26 = it1->first;
      tie(tmp27, tmp28, tmp29) = it1->second;
      if (tmp29>0){
        //Its particle-particle collision.
        fprintf(dataFile, "%d \t %.15f \t %d \t %d \t %d\n", tmp30, tmp26, tmp27, tmp28, tmp29);
        tmp30++;
      }
      //fprintf(dataFile, "%d \t %.15f \n", tmp26, tmp27 )
      //cout<<"\n"<<tmp25<<"\t"<<tmp21<<"\t\t"<<tmp22<<"\t"<<tmp23<<"\t"<<tmp24;
      //tmp25++;
      //if(tmp25>limit){
      //  break; //come out of printing the collision table.
      //}
      it1++;
    }
    fclose(dataFile);
  }
}

/////////////// savePosVel ////////////////////
void hardSphere::savePosVel()
{
  FILE *dataFile1;
  char str[300] = "posVel";
  char ch[400];
  strcat(str,"_phi_" ); sprintf(ch, "%0.4f", phi); strcat(str, ch);
  strcat(str, "_npart_"); sprintf(ch, "%d", npart); strcat(str, ch);
  strcat(str, "_nInit_"); sprintf(ch, "%d", nInit); strcat(str, ch);
  strcat(str, "_nColl_"); sprintf(ch, "%d", nColl); strcat(str, ch);
  strcat(str, ".xyz");
  dataFile1 = fopen(str, "w");
  if(!dataFile1)
  {
    cout<<"\n File not created###"; exit(0);
  }
  else
  {
    fprintf(dataFile1, "%d \n", npart);
    for (int ii=0; ii<npart; ii++)
    {
      fprintf(dataFile1, "\n %s \t", "Ar");
      for(int kk=0; kk<stride; kk++)
      {
        fprintf(dataFile1, "%.15f \t", posVel[ii*stride + kk]);
      }
    }
  }	fclose(dataFile1);
}//savePosVel ends.

//////////////// saveMinCollTime //////////////
void hardSphere::saveMinCollTime()
{
  FILE *dataFile2;
	char str[300] = "minCollTime";
	char ch[400];
  strcat(str,"_phi_" ); sprintf(ch, "%0.4f", phi); strcat(str, ch);
  strcat(str, "_npart_"); sprintf(ch, "%d", npart); strcat(str, ch);
  strcat(str, "_nInit_"); sprintf(ch, "%d", nInit); strcat(str, ch);
  strcat(str, ".dat");

	dataFile2 = fopen(str, "w");
	if(!dataFile2)
	{	cout<<"File not created###"; exit(0);	}
	else
	{
    //int pad, indexS;
    //pad = int(nColl/nCollSave);
    //indexS = (pad-1)*nCollSave;
		fprintf(dataFile2, "#Step \t tcol \n");
		for(int kk=0; kk<=nColl; kk++)
		{
      fprintf(dataFile2, "%d \t %.15f \n",kk, minCollTime[kk]);
    }

	}	fclose(dataFile2);
}

////////////// saveFreeFlightTime ////////////
void hardSphere::saveFreeFlightTime()
{
  FILE *dataFile2;
	char str[300] = "freeFlightTime";
	char ch[400];
  strcat(str,"_phi_" ); sprintf(ch, "%0.4f", phi); strcat(str, ch);
  strcat(str, "_npart_"); sprintf(ch, "%d", npart); strcat(str, ch);
  strcat(str, "_nInit_"); sprintf(ch, "%d", nInit); strcat(str, ch);
  strcat(str, ch); strcat(str, ".csv");
	dataFile2 = fopen(str, "w");
	if(!dataFile2)
	{	cout<<"Velocity file not created###"; exit(0);	}
	else
	{
		fprintf(dataFile2, "ParticleIndex \t \t AbsoluteCollTime\n");
		for (int ik=0; ik<npart; ik++)
		{
			fprintf(dataFile2, "%d", ik);
			for (int kk=0; kk<freeFlightTime[ik].size(); kk++)
			{
				fprintf(dataFile2, ",%.7f", freeFlightTime[ik][kk]);
			}
			fprintf(dataFile2, "\n");
		}
	}
	fclose(dataFile2);
}

/////////////////// saveKineticEne ////////////

void hardSphere::saveKineticEne()
{	//function to save the positional configuration of the particles in the box.
	FILE *dataFile2;
  char str[300] = "kineticEne";
  char ch[400];
  strcat(str,"_phi_" ); sprintf(ch, "%0.4f", phi); strcat(str, ch);
  strcat(str, "_npart_"); sprintf(ch, "%d", npart); strcat(str, ch);
  strcat(str, "_nInit_"); sprintf(ch, "%d", nInit); strcat(str, ch);
  strcat(str, ".dat");

	dataFile2 = fopen(str, "w");
	if(!dataFile2)
	{	cout<<"File not created###"; exit(0);	}
	else
	{	//fprintf(dataFile1, "%s \t %s \t %s \t %s", "Particle", "x-cord", "y-cord","z-cord");
		//fprintf(dataFile2, "%d \n", p->npart);
		fprintf(dataFile2, "#Step \t kineticEne \n");
		//for (int ii=0; ii<p->npart; ii++)
		//{	//fprintf(dataFile1, "\n %d \t", ii+1);
			//fprintf(dataFile2, "\n %s \t", "Ar");
			for(int kk=0; kk<nColl; kk++)
			{	fprintf(dataFile2, "%d \t %4.7f \n",kk, kineticEne[kk]);	}
			//fprintf(dataFile1, "\n");
		//}
	}	fclose(dataFile2);
}
/////////////// saveImpactSeparations //////////
void hardSphere::saveImpactSeparation()
{
  FILE *dataFile2;
  char str[300] = "impactSeparations";
  char ch[400];
  strcat(str,"_phi_" ); sprintf(ch, "%0.4f", phi); strcat(str, ch);
  strcat(str, "_npart_"); sprintf(ch, "%d", npart); strcat(str, ch);
  strcat(str, "_nInit_"); sprintf(ch, "%d", nInit); strcat(str, ch);
  strcat(str, ".dat");

  dataFile2 = fopen(str, "w");
  if(!dataFile2)
  {	cout<<"File not created###"; exit(0);	}
  else
  {	//fprintf(dataFile1, "%s \t %s \t %s \t %s", "Particle", "x-cord", "y-cord","z-cord");
  //fprintf(dataFile2, "%d \n", p->npart);
  fprintf(dataFile2, "#Step \t ImpactSeparations \n");
  //for (int ii=0; ii<p->npart; ii++)
  //{	//fprintf(dataFile1, "\n %d \t", ii+1);
    //fprintf(dataFile2, "\n %s \t", "Ar");
    for(int kk=0; kk<nColl; kk++) {
      fprintf(dataFile2, "%d \t %4.15f \n",kk, impactSeparation[kk]);
    }
    //fprintf(dataFile1, "\n");
  //}
}	fclose(dataFile2);
}
