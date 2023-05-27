#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "dcd.h"

#define LL 500
#define ncellx 6
#define ncelly 6
#define ncellz 6

/***** Creates a linked list.  *****/
/***** Calculates the potential energy using a linked list algorithm. *****/

int IX,IY,IZ,celno;

int icell(IX, IY, IZ){
	celno = ( IX-1+ncellx) % ncellx \
		  + ((IY-1+ncelly) % ncelly)*ncellx \
		  + ((IZ-1+ncellz) % ncellz)*ncellx*ncelly;
        return celno;
}


int main(int argc, char *argv[])
{
	float *x,*y,*z,lx,ly,lz,lx2,ly2,lz2,timestep,dx,dy,dz,fbin;
	int N,i,j,k,wcell,flag,nset,tbsave,t,tmax,whichfile;
	int nfiles,fileconfigs;
	long int pos=-1,n,norm,gt;
	char file[LL],filename[500][LL];
	FILE *input;
	long int read_dcd_head(FILE *input, int *N,int flag);
	int gdcdp(float *x, float *y, float *z, char file[],long int n,int flag,long int *pos, int wcell);
	long int dcd_info(char file[], int *N,int *nset, int *tbsave, float *timestep, int *wcell);
	double getboxlength(char filename[],int wcell,int config);
	/******/
	double cut, cut2, dr2, idr2, potential;
	float xlo,ylo,zlo,xtemp,ytemp,ztemp,ilencelx,ilencely,ilencelz;
	int imap,ncell,mapsize,jcell,jcell0,nabor;

	FILE *handle;
	
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

	xlo = -lx2 + 0.000001;
	ylo = -ly2 + 0.000001;
	zlo = -lz2 + 0.000001;

	ilencelx = ((float) ncellx)/lx;
	ilencely = ((float) ncelly)/ly;
	ilencelz = ((float) ncellz)/lz;
	ncell = ncellx*ncelly*ncellz;
	mapsize = 13*ncell;

	x = calloc(N,sizeof(float));
	y = calloc(N,sizeof(float));
	z = calloc(N,sizeof(float));

	norm = 0;
	flag = 0;
	cut = 2.5;
	cut2 = cut*cut;
	potential = 0.0;
	double pairpot[nset];

	for (i=0;i<nset;i++) pairpot[i] = 0.0;

/***********************************************************************************/

	int head[ncell],list[N],mapp[mapsize];

        for(IZ=1;IZ<ncellz+1;IZ++){
                for(IY=1;IY<ncelly+1;IY++){
                        for(IX=1;IX<ncellx+1;IX++){

                                imap = (icell(IX, IY, IZ)) * 13;

                                mapp[imap + 0] = icell( IX + 1, IY    , IZ     );
                                mapp[imap + 1] = icell( IX + 1, IY + 1, IZ     );
                                mapp[imap + 2] = icell( IX    , IY + 1, IZ     );
                                mapp[imap + 3] = icell( IX - 1, IY + 1, IZ     );
                                mapp[imap + 4] = icell( IX + 1, IY    , IZ - 1 );
                                mapp[imap + 5] = icell( IX + 1, IY + 1, IZ - 1 );
                                mapp[imap + 6] = icell( IX    , IY + 1, IZ - 1 );
                                mapp[imap + 7] = icell( IX - 1, IY + 1, IZ - 1 );
                                mapp[imap + 8] = icell( IX + 1, IY    , IZ + 1 );
                                mapp[imap + 9] = icell( IX + 1, IY + 1, IZ + 1 );
                                mapp[imap + 10] = icell( IX    , IY + 1, IZ + 1 );
                                mapp[imap + 11] = icell( IX - 1, IY + 1, IZ + 1 );
                                mapp[imap + 12] = icell( IX    , IY    , IZ + 1 );
                        }
                }
        }


        if((1.0/ilencelx) < cut || (1.0/ilencely) < cut || (1.0/ilencelz) < cut){
                printf("cell size too small for cutoff.");
                exit(0);
        }

/***********************************************************************************/


	for(n=0;n<nset;n++){

		for(i=0;i<ncell;i++) head[i] = -1;

		whichfile = n/fileconfigs;
		gt = n-fileconfigs*whichfile;
		gdcdp(x,y,z,filename[whichfile],gt,flag,&pos,wcell);

	        for(i=0;i<N;i++){

                        x[i] -= lx*rint(x[i]/lx);
                        y[i] -= ly*rint(y[i]/ly);
                        z[i] -= lz*rint(z[i]/lz);

        	        k  = (int)((x[i]-xlo)*ilencelx);
			k += (int)((y[i]-ylo)*ilencely) * ncellx;
			k += (int)((z[i]-zlo)*ilencelz) * ncellx * ncelly;

			if(k > ncell || k < 0){
				fprintf(stderr,"cell index is out of bound.");
				exit(0);
			}

	                list[i] = head[k];
        	        head[k] = i;
	        }

		for(celno=0;celno<ncell;celno++){

			i = head[celno];

			while(i != -1){

	                        xtemp = x[i];
        	                ytemp = y[i];
                	        ztemp = z[i];
                        	j = list[i];

				while(j != -1){

	                                dx = x[j] - xtemp;
        	                        dy = y[j] - ytemp;
                	                dz = z[j] - ztemp;
	
        	                        dx -= lx*rint(dx/lx);
                	                dy -= ly*rint(dy/ly);
                        	        dz -= lz*rint(dz/lz);
                                	dr2 = dx*dx + dy*dy + dz*dz;

	                                if(dr2 < cut2){

						idr2 = 1.0/dr2;
						pairpot[norm] += 4.0*(idr2*idr2*idr2*idr2*idr2*idr2 - idr2*idr2*idr2);
					}

					j = list[j];
				}

				jcell0 = 13*celno;

				for(nabor=0;nabor<13;nabor++){

	                                jcell = mapp[jcell0+nabor];
        	                        j = head[jcell];

                	                while(j != -1){

                        	                dx = x[j] - xtemp;
                                	        dy = y[j] - ytemp;
                                        	dz = z[j] - ztemp;

	                                        dx -= lx*rint(dx/lx);
        	                                dy -= ly*rint(dy/ly);
                	                        dz -= lz*rint(dz/lz);
                        	                dr2 = dx*dx + dy*dy + dz*dz;

                                	        if(dr2 < cut2){

							idr2 = 1.0/dr2;
							pairpot[norm] += 4.0*(idr2*idr2*idr2*idr2*idr2*idr2 - idr2*idr2*idr2);
						}

						j = list[j];
					}
				}

				i = list[i];
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

	fprintf(handle, "# timestep\tV(r)\n");
	for(i=0; i<norm; i++){

		fprintf(handle, "%i\t\t%.7f\n", i, pairpot[i]/N);
	}
	fclose(handle);
}
