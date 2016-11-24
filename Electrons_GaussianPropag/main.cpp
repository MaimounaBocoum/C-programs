#include <cmath>
#include <limits>
#include <iostream>
#include <fstream>
#include <cstdlib>
#include <string>
#include <vector>
#include <list>
#include <string.h>
#include <omp.h>
#include <Variables.h>



using namespace std;


//**************
// MAIN PROGRAM
//**************
int main ()
{
	std::vector<LaserPulse> lpv;
	std::vector<Wake> wkv;
	std::vector<ParticlePopulation> ppv;

	double tmin,tmax,t,dt;
	int ppvsz,save,tnot;
	int i=0,j=0,k=0;

	if(!readparam(lpv,wkv,ppv,tmin,tmax,dt,save))
	{
  	ppvsz=ppv.size();

		for(k=0;k<ppvsz;k++)
			ppv.at(k).MoveParticles(lpv,wkv,tmin,0);
		savetime(tmin,0);
		saveppv(ppv,0);

		for(t=tmin; t<=tmax; t+=dt) {
                // %: modulo operator
			if( (j++ % save == 0) && j > 0 ) {
				savetime(t,i);
				saveppv(ppv,i++);
				cout << t << endl;
			}
			for(k=0;k<ppvsz;k++)
				ppv.at(k).MoveParticles(lpv,wkv,t,dt);
		}
	}
	return(0);
}
