// Calculates the pressure tensor components: sigma_{i,j}, pressure, and the average pressure (09/01/2018).
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
	double cut,cut2,dr2,idr2,idr8,idr14,dpotential;
	double temp,vol,p11,p22,p33;

	FILE *handle_data;
	
	if(argc != 2){
		fprintf(stderr,"analysis nfiles\n");
		exit(1);
	}
	nfiles = atoi(argv[1]);

	for(i=0;i<nfiles;i++){
		sprintf(filename[i],"traj%i.dcd",i+1);
	}
	fprintf(stderr,"%s\n",filename[0]);
	dcd_info(filename[0],&N,&nset,&tbsave,&timestep,&wcell);

	fileconfigs = nset;
	nset = nset*nfiles;
	tmax = nset;

	lx = getboxlength(filename[0],wcell,0);
	ly = lx;
	lz = lx;
	lx2 = lx/2;
	ly2 = lx2;
	lz2 = lx2;

	x = calloc(N,sizeof(float));
	y = calloc(N,sizeof(float));
	z = calloc(N,sizeof(float));

	norm = 0;
	flag = 0;
	cut = 2.5;
	cut2 = cut*cut;
	dpotential = 0.0;

	p11 = 0.0;
	p22 = 0.0;
	p33 = 0.0;
	vol = lx*ly*lz;

	double sigma11[nset];
	double sigma22[nset];
	double sigma33[nset];
	double pressure[nset];

	
	for(n=100;n<nset;n++){

		whichfile = n/fileconfigs;
		gt = n-fileconfigs*whichfile;
		gdcdp(x,y,z,filename[whichfile],gt,flag,&pos,wcell);
				
		for(i=0;i<N;i++){

			for(j=i+1;j<N;j++){
				dx = x[i]-x[j];
				dy = y[i]-y[j];
				dz = z[i]-z[j];

				dx -= lx*rint(dx/lx);
				dy -= ly*rint(dy/ly);
				dz -= lz*rint(dz/lz);

				dr2 = dx*dx + dy*dy + dz*dz;

				if (dr2 < cut2){

					idr2 = 1.0/dr2;
					idr8 = idr2*idr2*idr2*idr2;
					idr14 = idr8*idr2*idr2*idr2;

					dpotential = 24.0*(idr8 - 2*idr14);	// =dpotential*(1/dr);

					p11 += (dx*dx)*(dpotential);
					p22 += (dy*dy)*(dpotential);
					p33 += (dz*dz)*(dpotential);
				}
			}
		}

		sigma11[norm] = -p11/vol;
		sigma22[norm] = -p22/vol;
		sigma33[norm] = -p33/vol;

		pressure[norm] = -(p11 + p22 + p33)/(3.0*vol);

		p11 = 0.0;
		p22 = 0.0;
		p33 = 0.0;
		norm++;
		fprintf(stderr,"%lu\n",n);
	}


	handle_data = fopen("datafile_press", "w");
	if (!handle_data) {
		printf("Unable to open datafile!");
		exit(1);
	}

	temp = 0.0;
	for(i=0; i<norm; i++){
		temp += pressure[i];
	}
	fprintf(handle_data, "<p> = %.7f\n\n",temp/norm);

	fprintf(handle_data, "# timestep\tpxx\t\tpyy\t\tpzz\t\tpress\n");
	for(i=0; i<norm; i++){

		fprintf(handle_data, "%i\t\t%.7f\t%.7f\t%.7f\t%.7f\n", i, sigma11[i], sigma22[i], sigma33[i], pressure[i]);
	}
	fclose(handle_data);

return 0;
}
