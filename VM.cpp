/*C++ CODE - MANGEAT MATTHIEU - 2023*/
/*VICSEK MODEL*/

//////////////////////
///// LIBRAIRIES /////
//////////////////////

//Public librairies.
#include <cstdlib>
#include <stdlib.h>
#include <cmath>
#include <vector>
#include <map>
#include <string>
#include <string.h>
#include <ctime>
#include <iostream>
#include <fstream>
#include <sstream>

using namespace std;

//Personal libraries.
#include "lib/random.cpp"
#include "lib/special_functions.cpp"

void modulo(double &x, const int &L);
int modulo(const int &x, const int &L);

//////////////////////////
///// PARTICLE CLASS /////
//////////////////////////

class particle
{
	public:
	
	double x,y; //Position of the particle.
	double dx, dy; //Displacement of the particle.
	double theta; //Orientation of the particle.
	double theta_avg; //Average orientation in the neighborhood.
	
	particle(const int &LX, const int &LY, const int &init); //Initial position and orientation.
	void updateTheta(const double &eta); //Update the orientation.
	void move(const double &v0, const int &LX, const int &LY); //Motion of the particle.
};

particle::particle(const int &LX, const int &LY, const int &init)
{
	//Displacement.
	dx=0.;
	dy=0.;
	
	//Average orientation.
	theta_avg=0.;
	
	if (init==0) //Gas phase.
	{
		theta=M_PI*(2*ran()-1);
		x=LX*ran();
		y=LY*ran();
	}
	else if (init==1) //Band moving right.
	{
		theta=0;
		x=(2+ran())*0.2*LX;
		y=LY*ran();
	}
	else if (init==2) //Liquid phase.
	{
		//static const double phi=M_PI*(2*ran()-1);
		static const double phi=0; //right liquid phase.
		theta=phi;
		x=LX*ran();
		y=LY*ran();
	}
	else
	{
		cerr << "BAD INIT VALUE: " << init << endl;
		abort();
	}
}

//Update the orientation.
void particle::updateTheta(const double &eta)
{
	theta = theta_avg + 2*M_PI*eta*(ran()-0.5);
	if (theta<-M_PI)
	{
		theta+=2*M_PI;
	}
	else if (theta>M_PI)
	{
		theta-=2*M_PI;
	}
}

//Update the position and the displacement with the new angle.
void particle::move(const double &v0, const int &LX, const int &LY)
{
	//Update of the displacement.
	dx += v0*cos(theta);
	dy += v0*sin(theta);	
	
	//Update of the position (+periodicity).
	x += v0*cos(theta);
	y += v0*sin(theta);
	modulo(x,LX);
	modulo(y,LY);
}

///////////////////////////
///// BASIC FUNCTIONS /////
///////////////////////////

//Modulo function.
void modulo(double &x, const int &L)
{
	if (x<0)
	{
		x+=L;
	}
	else if (x>=L)
	{
		x-=L;
	}
}

//Modulo function.
int modulo(const int &x, const int &L)
{
	if (x<0)
	{
		return x+L;
	}
	else if (x>=L)
	{
		return x-L;
	}
	else
	{
		return x;
	}
}

//Distance between two particles in periodic space.
double distance2(const particle &part1, const particle &part2, const int &LX, const int &LY)
{
	const double DX=fabs(part1.x-part2.x);
	const double DY=fabs(part1.y-part2.y);	
	return square(min(DX,LX-DX)) + square(min(DY,LY-DY));
}

//Order parameter.
vector<double> mag(const int &Npart, const vector<particle> &PART)
{
	double MX=0, MY=0;
	for (int i=0; i<Npart; i++)
	{
		MX+=cos(PART[i].theta);
		MY+=sin(PART[i].theta);
	}
	
	vector<double> mag(2,0.);
	mag[0]=MX/Npart; //mx
	mag[1]=MY/Npart; //my
	return mag;
}

//Mean square displacement.
vector<double> MSD(const int &Npart, const vector<particle> &PART)
{
	double DX=0, DY=0;
	double DX2=0., DY2=0.;
	for (int i=0; i<Npart; i++)
	{
		DX+=PART[i].dx/Npart;
		DX2+=square(PART[i].dx);
		DY+=PART[i].dy;
		DY2+=square(PART[i].dy);
	}	
	vector<double> DR(2,0.);
	DR[0]=(DX2+DY2)/Npart;
	DR[1]=square(DX/Npart)+square(DY/Npart);
	return DR;
}

//Average angle near each particles at time t.
void THETA_AVG(const int &LX, const int &LY, const int &Npart, vector<particle> &PART)
{
	//Create a matrix with particle locations, box[i][j] regroups the particle indices with i<x<i+1 and j<y<j+1.
	vector< vector< vector<int> > > box(LX,vector< vector<int> >(LY,vector<int>(0)));
	for (int i=0; i<Npart; i++)
	{
		box[int(PART[i].x)][int(PART[i].y)].push_back(i);
	}
	
	//Calculate the average angle of neighbours.
	for (int i=0; i<Npart; i++)
	{
		const int X0=int(PART[i].x), Y0=int(PART[i].y);
		
		//Take the average magnetisation for particles in neighbour boxes and with a distance smaller than 1.
		double MX=0., MY=0.;		
		for (int XN=X0-1;XN<=X0+1;XN++)
		{
			for (int YN=Y0-1;YN<=Y0+1;YN++)
			{
				const vector<int> neighbours=box[modulo(XN,LX)][modulo(YN,LY)];
				for (int l=0; l<neighbours.size(); l++)
				{
					const int k=neighbours[l];
					if (distance2(PART[i],PART[k],LX,LY)<1)
					{
						MX+=cos(PART[k].theta);
						MY+=sin(PART[k].theta);
					}
				}
			}
		}
		PART[i].theta_avg=atan2(MY,MX);
	}
}

//Export density in a file
void exportDensity(const double &v0, const double &eta, const double &rho0, const int &LX, const int &LY, const int &init, const int &RAN, const int &t, const vector<particle> &PART)
{
	//Local density.
	vector< vector<int> > RHO(LX,vector<int>(LY,0));
	for (int i=0; i<PART.size(); i++)
	{
		const int X0=int(PART[i].x), Y0=int(PART[i].y);
		RHO[X0][Y0]++;
	}
	
	//Creation of the filename.
	static const int dossier=system("mkdir -p ./data_VM_dynamics/");
	stringstream ssRHO;
	
	ssRHO << "./data_VM_dynamics/VM_rho_v0=" << v0 << "_eta=" << eta << "_rho0=" << rho0 << "_LX=" << LX << "_LY=" << LY << "_init=" << init << "_ran=" << RAN << "_t=" << t << ".txt";
	string nameRHO = ssRHO.str();			
	ofstream fileRHO(nameRHO.c_str(),ios::trunc);	
	fileRHO.precision(6);
	
	//Write in the file.
	for (int Y0=0;Y0<LY;Y0++)
	{
		for (int X0=0; X0<LX; X0++)
		{
			fileRHO << RHO[X0][Y0] << "\t";
		}
		fileRHO << endl;
	}
	fileRHO.close();
}

//Export magnetisation in a file
void exportMagnetization(const double &v0, const double &eta, const double &rho0, const int &LX, const int &LY, const int &init, const int &RAN, const int &t, const vector<particle> &PART)
{
	//Local magnetisation.
	vector< vector<double> > MX(LX,vector<double>(LY,0.)), MY(LX,vector<double>(LY,0.));
	for (int i=0; i<PART.size(); i++)
	{
		const int X0=int(PART[i].x), Y0=int(PART[i].y);
		MX[X0][Y0]+=cos(PART[i].theta);
		MY[X0][Y0]+=sin(PART[i].theta);
	}
	
	//Creation of filenames.
	static const int dossier=system("mkdir -p ./data_VM_dynamics/");
	stringstream ssMX,ssMY;
	
	ssMX << "./data_VM_dynamics/VM_mx_v0=" << v0 << "_eta=" << eta << "_rho0=" << rho0 << "_LX=" << LX << "_LY=" << LY << "_init=" << init << "_ran=" << RAN << "_t=" << t << ".txt";
	string nameMX = ssMX.str();	
	ofstream fileMX(nameMX.c_str(),ios::trunc);
	fileMX.precision(6);
	
	ssMY << "./data_VM_dynamics/VM_my_v0=" << v0 << "_eta=" << eta << "_rho0=" << rho0 << "_LX=" << LX << "_LY=" << LY << "_init=" << init << "_ran=" << RAN << "_t=" << t << ".txt";
	string nameMY = ssMY.str();	
	ofstream fileMY(nameMY.c_str(),ios::trunc);
	fileMY.precision(6);	
	
	//Write in files.	
	for (int Y0=0;Y0<LY;Y0++)
	{
		for (int X0=0; X0<LX; X0++)
		{
			fileMX << MX[X0][Y0] << "\t";
			fileMY << MY[X0][Y0] << "\t";
		}
		fileMX << endl;
		fileMY << endl;
	}
	fileMX.close();
	fileMY.close();
}

//Read parameters from command line.
void ReadCommandLine(int argc, char** argv, double &v0, double &eta, double &rho0, int &LX, int &LY, int &tmax, int &init, int &RAN)
{
 	for( int i = 1; i<argc; i++ )
	{
		
		if (strstr(argv[i], "-v0=" ))
		{
			v0=atof(argv[i]+4);
		}
		else if (strstr(argv[i], "-eta=" ))
		{
			eta=atof(argv[i]+5);
		}
		else if (strstr(argv[i], "-rho0=" ))
		{
			rho0=atof(argv[i]+6);
		}
		else if (strstr(argv[i], "-LX=" ))
		{
			LX=atof(argv[i]+4);
		}
		else if (strstr(argv[i], "-LY=" ))
		{
			LY=atof(argv[i]+4);
		}
		else if (strstr(argv[i], "-tmax=" ))
		{
			tmax=atoi(argv[i]+6);
		}
		else if (strstr(argv[i], "-init=" ))
		{
			init=atoi(argv[i]+6);
		}
		else if (strstr(argv[i], "-ran=" ))
		{
			RAN=atoi(argv[i]+5);
		}
		else
		{
			cerr << "BAD ARGUMENT : " << argv[i] << endl;
			abort();
		}
	}
}

/////////////////////
///// MAIN CODE /////
/////////////////////

int main(int argc, char *argv[])
{
	//Physical parameters: v0=self-propulsion velocity, eta=noise parameter, rho0=average density, LX*LY=size of the box.
	double v0=1.0, eta=0.3, rho0=0.64;
	int LX=512, LY=512;
	
	//Numerical parameters: init=initial condition, tmax=maximal time, RAN=index of RNG.
	int init=0, tmax=100000, RAN=0;

	//Read imported parameters in command line.
	ReadCommandLine(argc,argv,v0,eta,rho0,LX,LY,tmax,init,RAN);
	
	//Number of particles.
	const int Npart=int(rho0*LX*LY);

	//Start the random number generator.
	init_gsl_ran();
	cout << "GSL index (initial time) = " << RAN << "\n";
	gsl_rng_set(GSL_r,RAN);
	
	//Creation of the file to export global averages.
	const int dossier=system("mkdir -p ./data_VM_averages/");
	stringstream ssAverages;
	ssAverages << "./data_VM_averages/VM_averages_v0=" << v0 << "_eta=" << eta << "_rho0=" << rho0 << "_LX=" << LX << "_LY=" << LY << "_init=" << init << "_ran=" << RAN << ".txt";
	string nameAverages = ssAverages.str();			
	ofstream fileAverages(nameAverages.c_str(),ios::trunc);
	fileAverages.precision(6);	
	
	//Initial state.
	vector<particle> PART;
	for (int k=0; k<Npart; k++)
	{
		particle part(LX,LY,init);
		PART.push_back(part);
	}
	
	//Time evolution.
	for(int t=0; t<=tmax; t++)
	{
		//Export data.
		if (t%25==0)
		{
			//Order parameter.
			const vector<double> MAG=mag(Npart,PART);
			const double MX=MAG[0], MY=MAG[1];
			const double VALPHA=sqrt(MX*MX+MY*MY), PHI=atan2(MY,MX);
			
			//Mean-square displacement.
			vector<double> msd=MSD(Npart,PART);
			
			fileAverages << t << "\t" << rho0 << "\t" << VALPHA << "\t" << PHI << "\t" << MX << "\t" << MY << "\t" << msd[0] << "\t" << msd[1] << endl;
			cout << "time=" << t << " -N/V=" << rho0 << " -M=" << VALPHA << " -PHI=" << PHI << " -MX=" << MX << " -MY=" << MY << " -MSD=" << msd[0] << " -R2=" << msd[1] << running_time.TimeRun(" ") << endl;
			
			exportDensity(v0,eta,rho0,LX,LY,init,RAN,t,PART);
			exportMagnetization(v0,eta,rho0,LX,LY,init,RAN,t,PART);
		}
		
		//Update orientations and positions sequentially.
		THETA_AVG(LX,LY,Npart,PART);
		
		//Update orientations and positions sequentially.
		for (int i=0; i<Npart; i++)
		{
			PART[i].updateTheta(eta);
			PART[i].move(v0,LX,LY);		
		}
	}
	return 0;
}
