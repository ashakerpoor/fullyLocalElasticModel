// Self intermediate scattering function, calculated with a more efficient method (09/01/2018).
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

	//New variables
	double dt,k;
	float *rx,*ry,*rz;
	int m,len, index;
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
	printf("nset=%d\n", nset);
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
	dt = .002;

	//k = 2*pi/sigma	with	kx=(k,0,0); ky=(0,k,0); kz=(0,0,k)
	m = 1;
	k = 2*(4*atan(1));

	while(nset%m != 0){
		m++;
	}
	printf("m=%i\n",m);

	len = nset/m;
	rx = calloc((2*N),sizeof(float));
	ry = calloc((2*N),sizeof(float));
	rz = calloc((2*N),sizeof(float));

	if(rx==NULL){
                fprintf(stderr, "Not enough Xspace!\n");
                exit(1);
        }
	if(ry==NULL){
                fprintf(stderr, "Not enough Yspace!\n");
                exit(1);
        }
        if(rz==NULL){
                fprintf(stderr, "Not enough Zspace!\n");
                exit(1);
        }


	float relx[nset], rely[nset], relz[nset], imx[nset], imy[nset], imz[nset];
	int step[nset];

	for(i=0;i<len;i++){
		relx[i] = 0.0;
		rely[i] = 0.0;
		relz[i] = 0.0;
		imx[i] = 0.0;
		imy[i] = 0.0;
		imz[i] = 0.0;
		step[i] = 0;
        }

	index = 1;

	for(j=0;j<nset;j++){

		for(n=j;n<nset;n+=m){
			whichfile = n/fileconfigs;
			gt = n-fileconfigs*whichfile;
			gdcdp(x,y,z,filename[whichfile],gt,flag,&pos,wcell);

			if(index != j){
				for(i=0;i<N;i++){
					rx[i] = x[i];
					ry[i] = y[i];
					rz[i] = z[i];
				}
				index = j;
			}

			for(i=0;i<N;i++){
				rx[N+i] = x[i];
				ry[N+i] = y[i];
				rz[N+i] = z[i];
			}

			for(i=0;i<N;i++){
				relx[n-j] += cos(k*(rx[N+i] - rx[i]));
				rely[n-j] += cos(k*(ry[N+i] - ry[i]));
				relz[n-j] += cos(k*(rz[N+i] - rz[i]));

				imx[n-j] += sin(k*(rx[N+i] - rx[i]));
				imy[n-j] += sin(k*(ry[N+i] - ry[i]));
				imz[n-j] += sin(k*(rz[N+i] - rz[i]));
			}
			step[n-j] += 1;

		}
	}
	
	handle = fopen("datafile", "w");
	if (!handle) {
		printf("Unable to open datafile!");
		exit(1);
	}
	
	fprintf(handle, "# timestep(dt)\treal x\t\treal y\t\treal z\t\timaginary x\t\timaginary y\t\timaginary z\n");
	for(i=0; i<len; i++){

		fprintf(handle, "%f\t%f\t%f\t%f\t%f\t%f\t%f\n", i*m*dt, relx[i]/(N*step[i]), rely[i]/(N*step[i]), relz[i]/(N*step[i]), imx[i]/(N*step[i]), imy[i]/(N*step[i]), imz[i]/(N*step[i]));
	}
	fclose(handle);
}
