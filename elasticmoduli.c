/*
 * Calculates the local shear moduli in Mizuno and Mossa 2013 paper for smaller boxes.
 * Calculates one specified shear/bulk modulus matrix component at a time: C_ijkl
 * The code uses the lines and planes approach and the linked-list algorithm.
*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <stdbool.h>
#include "dcd.h"
#define LL 500

/*
 * Cell characteristics of the overlapping boxes.
 * ngrid-1 =  number of points inside big cube *** ngrid = number of boxes.
*/

#define nbox 3
#define ngrid 7
#define ngrid3 ngrid*ngrid*ngrid

/*
 * Cell characteristics of the linked-list algorithm.
*/

#define ncellx 6
#define ncelly 6
#define ncellz 6

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
	int N,i,j,k,l,wcell,flag,nset,tbsave,t,tmax,whichfile;
	int nfiles,fileconfigs;
	long int pos=-1,n,norm,gt;
	char file[LL],filename[500][LL];
	FILE *input;
	long int read_dcd_head(FILE *input, int *N,int flag);
	int gdcdp(float *x, float *y, float *z, char file[],long int n,int flag,long int *pos, int wcell);
	long int dcd_info(char file[], int *N,int *nset, int *tbsave, float *timestep, int *wcell);
	double getboxlength(char filename[],int wcell,int config);

	double cut,cut2,dr2,idr2,idr8,idr14,dpotential,d2potential;
	double vol,w,iw,w3,T,rho;
	double xlo,ylo,zlo,temp,factor,midpoint,move;
	float xImage,yImage,zImage;
	float xi,yi,zi,xj,yj,zj,xtemp,ytemp,ztemp;
	double sigIJboxsinglts[ngrid3][nbox][nbox][nbox],sigIJbox[ngrid3][nbox][nbox][nbox],multbox[ngrid3][nbox][nbox][nbox];
	double cBm[ngrid3][nbox][nbox][nbox],cNm[ngrid3][nbox][nbox][nbox],cKm[ngrid3][nbox][nbox][nbox];
        double Born_m[ngrid3][nbox][nbox][nbox];
	double rab[4],a[20],grid[20];
	double sigmaIJ,sigmaKL,sigIJcube,sigKLcube,multcube,Born,cB,cN,cK,cTotal;
	float sigsum[ngrid3],sigbig;
	int p,q,r,s,deltaIJ,deltaKL,deltaIK,deltaJL,deltaIL,deltaJK;
	int boxpop[ngrid3][nbox][nbox][nbox],rho_m[ngrid3][nbox][nbox][nbox];
	int alpha,beta,gamma,phi,chi,psi,an,bn,cn;
	int lnum,hnum,swap,count,m;

	/** LINKED LIST VARIABLES **/

	double ilencelx,ilencely,ilencelz;
	int imap,ncell,mapsize,jcell,jcell0,nabor;

	FILE *handle1,*handle2;

	if(argc != 6){
		fprintf(stderr,"input nfiles and 4 indices (i,j,k,l).\n");
		exit(1);
	}

	nfiles = atoi(argv[1]);
	p = atoi(argv[2]);
        q = atoi(argv[3]);
        r = atoi(argv[4]);
        s = atoi(argv[5]);

        printf("i j k l = %i %i %i %i\n",p,q,r,s);

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

	x = calloc(N,sizeof(float));
	y = calloc(N,sizeof(float));
	z = calloc(N,sizeof(float));

	norm = 0;
	m = 0;
	flag = 0;

	cut = 2.5;
	cut2 = cut*cut;
	T = 1e-3;
	rho = 1.015;

	w = lx/nbox;
	iw = 1.0/w;
	w3 = w*w*w;
	vol = lx*ly*lz;

	deltaIJ = 0;
	deltaKL = 0;
	deltaIK = 0;
	deltaJL = 0;
	deltaIL = 0;
	deltaJK = 0;

	if(p==q) deltaIJ = 1;
	if(r==s) deltaKL = 1;
	if(p==r) deltaIK = 1;
	if(q==s) deltaJL = 1;
	if(p==s) deltaIL = 1;
	if(q==r) deltaJK = 1;

	rab[0] = 0.0;
	a[0] = 0.0;
	grid[0] = xlo;
        move = lx/((float) ngrid);


/***********************************************************************************/
        /** LINKED LIST CELLS **/

        ilencelx = ((float) ncellx)/lx;
        ilencely = ((float) ncelly)/ly;
        ilencelz = ((float) ncellz)/lz;

        ncell = ncellx*ncelly*ncellz;
        mapsize = 13*ncell;
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
                printf("cell size too small for cutoff.\n");
                exit(0);
        }

/***********************************************************************************/


	for(l=0;l<ngrid3;l++){
		for(i=0;i<nbox;i++){
			for(j=0;j<nbox;j++){
				for(k=0;k<nbox;k++){

					boxpop[l][i][j][k] = 0;
					rho_m[l][i][j][k] = 0;
					Born_m[l][i][j][k] = 0.0;
					sigIJbox[l][i][j][k] = 0.0;
					sigIJboxsinglts[l][i][j][k] = 0.0;
					multbox[l][i][j][k] = 0.0;

				}
			}
		}
		sigsum[l] = 0.0;
	}

	Born = 0.0;
	sigmaIJ = 0.0;
	sigmaKL = 0.0;
	multcube = 0.0;
	sigIJcube = 0.0;
	sigKLcube = 0.0;

	for(i=1;i<ngrid;i++) grid[i] = grid[i-1] + move;


	for(n=100;n<nset;n++){

		whichfile = n/fileconfigs;
		gt = n-fileconfigs*whichfile;
		gdcdp(x,y,z,filename[whichfile],gt,flag,&pos,wcell);

		for(i=0;i<ncell;i++) head[i] = -1;

                xlo = grid[0];
                ylo = grid[0];
                zlo = grid[0];

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

						rab[1] = dx;
						rab[2] = dy;
						rab[3] = dz;

	                                        idr2 = 1.0/dr2;
        	                                idr8 = idr2*idr2 * idr2*idr2;
                	                        idr14 = idr8 * idr2*idr2*idr2;
                        	                dpotential = 24*(idr8 - 2*idr14); // this is dpotential*(1/dr)
                                	        d2potential = 24*(26*idr14 - 7*idr8);
                                        	factor = (d2potential - dpotential)*(rab[p]*rab[q]*rab[r]*rab[s])*idr2;

						Born += factor;
						sigmaIJ += dpotential*rab[p]*rab[q];
						sigmaKL += dpotential*rab[r]*rab[s];
					}


	                                for(an=0;an<ngrid;an++){

        	                                xlo = grid[an];
                	                        xi = xtemp;
                        	                xj = x[j];

                                	        if(xi < xlo) xi += lx;
        	                                if(xj < xlo) xj += lx;

	                                        for(bn=0;bn<ngrid;bn++){

        	                                        ylo = grid[bn];
                	                                yi = ytemp;
                        	                        yj = y[j];

	                                                if(yi < ylo) yi += ly;
                        	                        if(yj < ylo) yj += ly;

	                                                for(cn=0;cn<ngrid;cn++){

        	                                                zlo = grid[cn];
                	                                        zi = ztemp;
                        	                                zj = z[j];
	
        	                                                if(zi < zlo) zi += lz;
                                	                        if(zj < zlo) zj += lz;

	                                                        alpha = (int)((xi-xlo)*iw);
        	                                                beta = (int)((yi-ylo)*iw);
                	                                        gamma = (int)((zi-zlo)*iw);
	
        	                                                if(j == list[i]) boxpop[m][alpha][beta][gamma] += 1;

        	                                                if(dr2 < cut2){

	                                                                xImage = dx + xi;
                	                                                yImage = dy + yi;
                        	                                        zImage = dz + zi;

                                	                                phi = (int)((xImage-xlo)*iw);
                                        	                        chi = (int)((yImage-ylo)*iw);
                                                	                psi = (int)((zImage-zlo)*iw);

	                                                                if(xImage < xlo) phi--;
                        	                                        if(yImage < ylo) chi--;
                                                	                if(zImage < zlo) psi--;

	                                                                count = 1;

	                                                                if(phi>alpha){
        	                                                                hnum = phi+1;
                	                                                        lnum = alpha+1;
                        	                                                for(k=lnum;k<hnum;k++){
                                	                                                a[count] = (k*w-(xi-xlo))/dx;
                                        	                                        count++;
                                                	                        }
                                                        	        }
	                                                                else{
        	                                                                hnum = alpha+1;
                	                                                        lnum = phi+1;
                        	                                                for(k=lnum;k<hnum;k++){
                                	                                                a[count] = (k*w-(xi-xlo))/dx;
                                        	                                        count++;
                                                	                        }
                                                        	        }
	                                                                if(chi>beta){
        	                                                                hnum = chi+1;
                	                                                        lnum = beta+1;
                        	                                                for(k=lnum;k<hnum;k++){
                                	                                                a[count] = (k*w-(yi-ylo))/dy;
                                        	                                        count++;
                                                	                        }
                                                        	        }
	                                                                else{
        	                                                                hnum = beta+1;
                	                                                        lnum = chi+1;
                        	                                                for(k=lnum;k<hnum;k++){
                                	                                                a[count] = (k*w-(yi-ylo))/dy;
                                        	                                        count++;
                                                	                        }
									}
	                                                                if(psi>gamma){
        	                                                                hnum = psi+1;
                	                                                        lnum = gamma+1;
                        	                                                for(k=lnum;k<hnum;k++){
                                	                                                a[count] = (k*w-(zi-zlo))/dz;
                                        	                                        count++;
                                                	                        }
                                                        	        }
	                                                                else{
        	                                                                hnum = gamma+1;
                	                                                        lnum = psi+1;
                        	                                                for(k=lnum;k<hnum;k++){
                                	                                                a[count] = (k*w-(zi-zlo))/dz;
                                        	                                        count++;
                                                	                        }
                                                        	        }


	                                                                a[count] = 1.0;
        	                                                        swap = 1;

                	                                                while(swap==1){
                        	                                                swap = 0;
                                	                                        for(k=1;k<count;k++){
                                        	                                        if(a[k-1]>a[k]){
                                                	                                        temp = a[k];
                                                        	                                a[k] = a[k-1];
                                                                	                        a[k-1] = temp;
                                                                        	                swap = 1;
                                                                                	}
                                                                        	}
                                                                	}


        	                                                        count++;

	                                                                for(k=1;k<count;k++){

	                                                                        midpoint = a[k-1] + (a[k]-a[k-1])/2;
        	                                                                phi = (int)((midpoint*dx + (xi-xlo))*iw);
                	                                                        chi = (int)((midpoint*dy + (yi-ylo))*iw);
                        	                                                psi = (int)((midpoint*dz + (zi-zlo))*iw);
	
	                                                                        if((midpoint*dx + (xi-xlo))<0){
        	                                                                        phi--;
                	                                                                phi += nbox;
                        	                                                }
                                	                                        else if(phi>nbox-1){
                                        	                                        phi -= nbox;
                                                	                        }

                                                        	                if((midpoint*dy + (yi-ylo))<0){
                                                                	                chi--;
                                                                        	        chi += nbox;
	                                                                        }
        	                                                                else if(chi>nbox-1){
                	                                                                chi -= nbox;
                        	                                                }
	
	                                                                        if((midpoint*dz + (zi-zlo))<0){
        	                                                                        psi--;
                	                                                                psi += nbox;
                        	                                                }
	                                                                        else if(psi>nbox-1){
        	                                                                        psi -= nbox;
                	                                                        }
	
	                                                                        Born_m[m][phi][chi][psi] += factor*(a[k]-a[k-1]);
										sigIJboxsinglts[m][phi][chi][psi] += dpotential*rab[p]*rab[q]*(a[k]-a[k-1]);
                	                                                }
                        	                                }

                                	                        m++;
                                        	        }
	                                        }
        	                        }

                	                m = 0;
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

							rab[1] = dx;
							rab[2] = dy;
							rab[3] = dz;

		                                        idr2 = 1.0/dr2;
                		                        idr8 = idr2*idr2 * idr2*idr2;
                                		        idr14 = idr8 * idr2*idr2*idr2;
		                                        dpotential = 24*(idr8 - 2*idr14); // this is dpotential*(1/dr)
                		                        d2potential = 24*(26*idr14 - 7*idr8);
                                		        factor = (d2potential - dpotential)*(rab[p]*rab[q]*rab[r]*rab[s])*idr2;

							Born += factor;
							sigmaIJ += dpotential*rab[p]*rab[q];
							sigmaKL += dpotential*rab[r]*rab[s];
		                                }


		                                for(an=0;an<ngrid;an++){

                		                        xlo = grid[an];
                                		        xi = xtemp;
		                                        xj = x[j];
	
		                                        if(xi < xlo) xi += lx;
		                                        if(xj < xlo) xj += lx;
	
		                                        for(bn=0;bn<ngrid;bn++){

                		                                ylo = grid[bn];
                                		                yi = ytemp;
                                                		yj = y[j];
	
		                                                if(yi < ylo) yi += ly;
		                                                if(yj < ylo) yj += ly;

		                                                for(cn=0;cn<ngrid;cn++){

                		                                        zlo = grid[cn];
                                		                        zi = ztemp;
                                                		        zj = z[j];

		                                                        if(zi < zlo) zi += lz;
		                                                        if(zj < zlo) zj += lz;


		                                                        if(dr2 < cut2){

	                                                                        alpha = (int)((xi-xlo)*iw);
        	                                                                beta = (int)((yi-ylo)*iw);
                	                                                        gamma = (int)((zi-zlo)*iw);

                		                                                xImage = dx + xi;
                                		                                yImage = dy + yi;
                                                		                zImage = dz + zi;

		                                                                phi = (int)((xImage-xlo)*iw);
                		                                                chi = (int)((yImage-ylo)*iw);
                                		                                psi = (int)((zImage-zlo)*iw);
	
		                                                                if(xImage < xlo) phi--;
		                                                                if(yImage < ylo) chi--;
		                                                                if(zImage < zlo) psi--;

                                                		                count = 1;
	
		                                                                if(phi>alpha){
                		                                                        hnum = phi+1;
                                		                                        lnum = alpha+1;
                                                		                        for(k=lnum;k<hnum;k++){
                                                                		                a[count] = (k*w-(xi-xlo))/dx;
                                                                                		count++;
		                                                                        }
                		                                                }
                                		                                else{
                                                		                        hnum = alpha+1;
                                                                		        lnum = phi+1;
		                                                                        for(k=lnum;k<hnum;k++){
                		                                                                a[count] = (k*w-(xi-xlo))/dx;
                                		                                                count++;
		                                                                        }
                		                                                }
		                                                                if(chi>beta){
                		                                                        hnum = chi+1;
                                		                                        lnum = beta+1;
                                                		                        for(k=lnum;k<hnum;k++){
                                                                		                a[count] = (k*w-(yi-ylo))/dy;
                                                                                		count++;
		                                                                        }
                		                                                }
		                                                                else{
                		                                                        hnum = beta+1;
                                		                                        lnum = chi+1;
                                                		                        for(k=lnum;k<hnum;k++){
                                                                		                a[count] = (k*w-(yi-ylo))/dy;
                                                                                		count++;
		                                                                        }
										}
	        	                                                        if(psi>gamma){
        	        	                                                        hnum = psi+1;
                                		                                        lnum = gamma+1;
                                                		                        for(k=lnum;k<hnum;k++){
                                                                		                a[count] = (k*w-(zi-zlo))/dz;
                                                                                		count++;
		                                                                        }
                		                                                }
		                                                                else{
                		                                                        hnum = gamma+1;
                                		                                        lnum = psi+1;
                                                		                        for(k=lnum;k<hnum;k++){
                                                                		                a[count] = (k*w-(zi-zlo))/dz;
                                                                                		count++;
		                                                                        }
                		                                                }
	

		                                                                a[count] = 1.0;
                		                                                swap = 1;

                                		                                while(swap==1){
                                                		                        swap = 0;
                                                                		        for(k=1;k<count;k++){
                                                                                		if(a[k-1]>a[k]){
		                                                                                        temp = a[k];
                		                                                                        a[k] = a[k-1];
                                		                                                        a[k-1] = temp;
                                                		                                        swap = 1;
                                                                		                }
		                                                                        }
                		                                                }


		                                                                count++;

                		                                                for(k=1;k<count;k++){

                                		                                        midpoint = a[k-1] + (a[k]-a[k-1])/2;
                                                		                        phi = (int)((midpoint*dx + (xi-xlo))*iw);
                                                                		        chi = (int)((midpoint*dy + (yi-ylo))*iw);
		                                                                        psi = (int)((midpoint*dz + (zi-zlo))*iw);
	
		                                                                        if((midpoint*dx + (xi-xlo))<0){
                		                                                                phi--;
                                		                                                phi += nbox;
                                                		                        }
		                                                                        else if(phi>nbox-1){
                		                                                                phi -= nbox;
                                		                                        }
	
		                                                                        if((midpoint*dy + (yi-ylo))<0){
                		                                                                chi--;
                                		                                                chi += nbox;
                                                		                        }
		                                                                        else if(chi>nbox-1){
                		                                                                chi -= nbox;
                                		                                        }
	
		                                                                        if((midpoint*dz + (zi-zlo))<0){
                		                                                                psi--;
                                		                                                psi += nbox;
		                                                                        }
                		                                                        else if(psi>nbox-1){
                                		                                                psi -= nbox;
                                                		                        }
	
		                                                                        Born_m[m][phi][chi][psi] += factor*(a[k]-a[k-1]);
											sigIJboxsinglts[m][phi][chi][psi] += dpotential*rab[p]*rab[q]*(a[k]-a[k-1]);
                                		                                }
		                                                        }

                		                                        m++;
                                		                }
		                                        }
                		                }

                                		m = 0;
						j = list[j];
					}
				}

				i = list[i];
			}
		}


		for(l=0;l<ngrid3;l++){
			for(i=0;i<nbox;i++){
				for(j=0;j<nbox;j++){
					for(k=0;k<nbox;k++){

						rho_m[l][i][j][k] += boxpop[l][i][j][k];
					}
				}
			}
		}

		for(l=0;l<ngrid3;l++){
			for(i=0;i<nbox;i++){
				for(j=0;j<nbox;j++){
					for(k=0;k<nbox;k++){

						sigIJboxsinglts[l][i][j][k] = -(boxpop[l][i][j][k]/w3)*T*deltaIJ + sigIJboxsinglts[l][i][j][k]/w3;

					}
				}
			}
		}


                sigmaIJ = -rho*T*deltaIJ + sigmaIJ/vol;
                sigmaKL = -rho*T*deltaKL + sigmaKL/vol;

		multcube += sigmaIJ*sigmaKL;
                sigIJcube += sigmaIJ;
                sigKLcube += sigmaKL;


                for(l=0;l<ngrid3;l++){
                        for(i=0;i<nbox;i++){
                                for(j=0;j<nbox;j++){
                                        for(k=0;k<nbox;k++){

						multbox[l][i][j][k] += sigIJboxsinglts[l][i][j][k]*sigmaKL;
						sigIJbox[l][i][j][k] += sigIJboxsinglts[l][i][j][k];

					}
				}
			}
		}


		for(l=0;l<ngrid3;l++){
			for(i=0;i<nbox;i++){
				for(j=0;j<nbox;j++){
					for(k=0;k<nbox;k++){

						boxpop[l][i][j][k] = 0;
						sigIJboxsinglts[l][i][j][k] = 0.0;

					}
				}
			}
		}

		sigmaIJ = 0.0;
		sigmaKL = 0.0;

		norm++;
		fprintf(stderr,"%lu\n",n);
	}


	sigbig = sigIJcube/norm;

	cB = Born/(norm*vol);
	cN = (vol/T)*(multcube/norm - (sigIJcube/norm)*(sigKLcube/norm));
	cK = 2*rho*T*(deltaIK*deltaJL + deltaIL*deltaJK);
	cTotal = cB - cN + cK;


	for(l=0;l<ngrid3;l++){
		for(i=0;i<nbox;i++){
			for(j=0;j<nbox;j++){
				for(k=0;k<nbox;k++){

					sigsum[l] += sigIJbox[l][i][j][k]/norm;

				}
			}
		}

		sigsum[l] *= (w3/vol);
	}


	for(l=0;l<ngrid3;l++){
		for(i=0;i<nbox;i++){
			for(j=0;j<nbox;j++){
				for(k=0;k<nbox;k++){

					cBm[l][i][j][k] = Born_m[l][i][j][k]/(norm*w3);
					cNm[l][i][j][k] = (vol/T)*(multbox[l][i][j][k]/norm - (sigIJbox[l][i][j][k]/norm)*(sigKLcube/norm));
					cKm[l][i][j][k] = 2*(rho_m[l][i][j][k]/(norm*w3))*T*(deltaIK*deltaJL + deltaIL*deltaJK);

				}
			}
		}
	}


	handle1 = fopen("backup.dat", "w");
	if (!handle1) {
		printf("Unable to open datafile!");
		exit(0);
	}


	handle2 = fopen("extract.dat", "w");
	if (!handle2) {
		printf("Unable to open datafile!");
		exit(0);
	}


	fprintf(handle1, "# i j k l\tcB\t\tcN\t\tcK\t\tcTotal\n");
	fprintf(handle1, "  %i %i %i %i\t%f\t%f\t%f\t%f\n\n", p,q,r,s,cB,cN,cK,cTotal);


	for(l=0;l<ngrid3;l++){

		fprintf(handle1, "# grid no.\tsigmaIJ\t\tsigmaIJ_m sum\n");
		fprintf(handle1, "\t%i\t%e\t%e\n",l+1,sigbig,sigsum[l]);
		fprintf(handle1, "# box#\tcBm\t\tcNm\t\tcKm\t\tcTotalm\t\trho_m\n");
		m = 1;

		for(i=0;i<nbox;i++){
			for(j=0;j<nbox;j++){
				for(k=0;k<nbox;k++){

					fprintf(handle1, " %i\t%f\t%f\t%f\t%f\t%f\n",m,cBm[l][i][j][k],cNm[l][i][j][k],cKm[l][i][j][k],(cBm[l][i][j][k]-cNm[l][i][j][k]+cKm[l][i][j][k]),rho_m[l][i][j][k]/(norm*w3));
					fprintf(handle2, "%e\n",(cBm[l][i][j][k]-cNm[l][i][j][k]+cKm[l][i][j][k]));
					m++;
				}
			}
		}
	}

        fclose(handle1);
        fclose(handle2);
	return 0;
}
