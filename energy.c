// Calculates the energy of the system (kinetic energy excluded) (09/01/2018).
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "dcd.h"

#define LL 500

int main(int argc, char *argv[])
{
	float *x,*y,*z,lx,ly,lz,lx2,ly2,lz2,timestep,dx,dy,dz,fbin;
	int N,i,j,wcell,flag,nset,tbsave,t,tmax,whichfile;
	int nfiles,fileconfigs;
	long int pos=-1,n,norm,gt;
	char file[LL],filename[500][LL];
	FILE *input;
	long int read_dcd_head(FILE *input, int *N,int flag);
	int gdcdp(float *x, float *y, float *z, char file[],long int n,int flag,long int *pos, int wcell);
	long int dcd_info(char file[], int *N,int *nset, int *tbsave, float *timestep, int *wcell);
	double getboxlength(char filename[],int wcell,int config);

	//new variables
	double cut, cutsq, dRsq, idRsq, potential, potentialave;
	//double hist[maxbin];
	int k;
	FILE *handle;
	
	//Check to make sure there are enough arguments on the command line. 	
	//If not print an error message giving information on how to run the program. 
	if(argc != 2){
		fprintf(stderr,"analysis nfiles\n");
		exit(1);
	}
	nfiles = atoi(argv[1]);

	//Populate the list if file names if the trajectory is spread between more that one
	//file. May never need to do this these days. 
	for(i=0;i<nfiles;i++){
		sprintf(filename[i],"traj%i.dcd",i+1);
	}
	//Print out a check to make sure you are starting with the correct filename. 
	fprintf(stderr,"%s\n",filename[0]);
	//Read information stored in the header of the dcd file. It stores much useful information. 
	//The file dcd.h has the code needed to read the header, and you can see the format in that file. 
	dcd_info(filename[0],&N,&nset,&tbsave,&timestep,&wcell);

	//nset is the total number of configurations, which is the number of configurations in one file
	//times the number of files. 	
	fileconfigs = nset;
	nset = nset*nfiles;
	tmax = nset;

	//Read in the boxlength. This function only reads in the x value of the boxlength and 
	//I assume that the box is a cube. We can work on writing a function that reads in all
	//the box length parameters if we need. 	
	lx = getboxlength(filename[0],wcell,0);
	ly = lx;
	lz = lx;
	lx2 = lx/2;
	ly2 = lx2;
	lz2 = lx2;

	//Allocate the memory.
	x = calloc(N,sizeof(float));
	y = calloc(N,sizeof(float));
	z = calloc(N,sizeof(float));

	norm = 0;
	flag = 0;

	//Do the needed calcualtion. Here we are calculating an average over every 10th
	//configuration of the trajectory.
	
	//fprintf(stderr,"lx = %f\n",lx);

	cut = 2.5;
	cutsq = cut*cut;
	k = 1;
	potential = 0.0;
	double pairpot[nset];

	for (i=0;i<nset;i++) pairpot[i] = 0.0;
	
	for(n=0;n<nset;n+=k){
		//Determine which file we need to read from.
		whichfile = n/fileconfigs;
		//What configuration in the file.
		gt = n-fileconfigs*whichfile;
		//Read the x,y,z positions.
		gdcdp(x,y,z,filename[whichfile],gt,flag,&pos,wcell);
				
//		for (i=0;i<N;i++){
//			kinetic += v(n)*v(n)/2;
//		}
		for(i=0;i<N-1;i++){
			for(j=i+1;j<N;j++){
				dx = x[i]-x[j];
				dy = y[i]-y[j];
				dz = z[i]-z[j];

				dx -= lx*rint(dx/lx);
				dy -= ly*rint(dy/ly);
				dz -= lz*rint(dz/lz);

				dRsq = dx*dx + dy*dy + dz*dz;
				if (dRsq < cutsq){
					idRsq = 1.0/dRsq;
					pairpot[norm] += 4.0*(idRsq*idRsq*idRsq*idRsq*idRsq*idRsq - idRsq*idRsq*idRsq);
				}
			}
		}

		norm++;
		fprintf(stderr,"%lu\n",n);
	}
	
	handle = fopen("datafile", "w");
	if (!handle) {
		printf("Unable to open datafile!");
		exit(1);
	}

	fprintf(handle, "# <V> = %f\n", potentialave);	
	fprintf(handle, "# timestep\tV(r)\n");
	for(i=0; i<norm; i++){

		fprintf(handle, "%i\t\t%.7f\n", i*k, pairpot[i]/N);
	}
	fclose(handle);
}
