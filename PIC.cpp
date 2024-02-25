#include "stdio.h"
#include "math.h"
#include <iostream>
#include <fstream>
#include <string>
#include <stdlib.h>
#include <complex>

using namespace std;
#define NP 512 //number of particles
#define NC 4 //number of cells in each direction
complex <double> NoC = 4.0;
complex <double> NoP = 512.0;

const double pi = std::acos(-1);
const std::complex<double> ii(0, 1);

typedef struct tagPART
{
 double x, y; //coordinates
 double vx, vy; //velocity
 double vx0, vy0; //initial velocity
 double vpar; //velocity component along initial
 double vperp; //velocity component perpendicular to initial
 double v,v0;
} PART;
double m = 9.11e-28; //g   
double k = 1.38e-16; //boltzmann 
complex<double> Ex[NC][NC];  //field components
complex<double> Ey[NC][NC]; 
complex<double> phif [NC][NC];
complex <double> phi [NC][NC];
complex<double> densCIC [NC][NC]; //density for cloud-in-cell  method
complex<double> densCICf [NC][NC];
complex<double> densTSC [NC][NC]; //density for TSC  method
complex<double> densTSCf [NC][NC];
PART electrons[NP];  
double T = 2.0;// 
double H = 1.0; //cell size in lambda_d
double t = 0.5;  //timestep in wpe^-1
double Vte = sqrt(T)*t/H; //initial electrons thermal speed
double vmax = 5.0*Vte;
double L = double(NC);  //area size in cells
double A = 0.05; //amplitude of electric field fluctuations

complex<double> W=L*L*t*t/NP*0.5;

//random value from min to max
double GetRand (double min, double max)
{
 double x = min+(rand()/ (RAND_MAX + 1.0))*(max-min);
 return (x);
}

//Maxwell distributin
double fmaxwell (double v)
{
 double fv = exp(-v*v/(2.0*T));
 return (fv);
}

//random value distributed by Maxwell, any method you prefer
double RandMaxwell (double Vt)
{
 int a=0.0;
 double f0;
 double v;
 do
	{
	 v = GetRand(-vmax,vmax);
	 f0 = GetRand(0,1.0);
	 if (f0 <= fmaxwell(v)) a = 1.0;
	}
 while (a != 1.0);
 return (v);
}

//initialize particles
void Init (void)
{
 int i, j;
 double v,v0;
 for (i = 0; i < NP; i++)
	{
	 electrons[i].x = GetRand (0,L);
	 electrons[i].y = GetRand (0,L);
	 electrons[i].vx = RandMaxwell (Vte);
	 electrons[i].vy = RandMaxwell (Vte);
	 electrons[i].vx0 = electrons[i].vx;
     electrons[i].vy0 = electrons[i].vy;
	 electrons[i].v = sqrt (electrons[i].vx*electrons[i].vx +
	 + electrons[i].vy*electrons[i].vy);
	 electrons[i].v0 = sqrt (electrons[i].vx0*electrons[i].vx0 +
	 + electrons[i].vy0*electrons[i].vy0);
	 electrons[i].vpar = v;
	 electrons[i].vperp = 0.0;
	}
}

void TSCmethod (void)
{
 for (int i = 0; i < NC; i++)
	{
	 for (int k = 0; k < NC; k++)
		{
		 densTSC [i][k]=W*NoP/NoC;
		}	
	}
 double a = 0.0; 
 double b = 0.0; 
 double c = 0.0; 
 double d = 0.0;
 double sx = 0.0; 
 double lx ;
 double rx ;
 double xx ;
 double ly;
 double ry ;
 double yy;
 double* cvalx = new double [NC]; //cells center coordinates
 double* cvaly = new double [NC];
 for(int i = 0; i<NC; i++)
	{
	 cvalx[i]=i+0.5;
	 cvaly[i]=i+0.5;
	}
 for(int i = 0; i<NP; i++)
	{	
	 b=electrons[i].x-(cvalx[int(floor(electrons[i].x))]-0.5);
	 a=(cvalx[int(floor(electrons[i].x))]+0.5)-electrons[i].x;
	 d=electrons[i].y-(cvaly[int(floor(electrons[i].y))]-0.5);
	 c=(cvaly[int(floor(electrons[i].y))]+0.5)-electrons[i].y;	
		
	 lx = (a*a)*0.5;
	 rx = (b*b)*0.5;
	 xx = 1.0-((a*a*0.5)+(b*b*0.5));
	 ly = (c*c)*0.5;
	 ry = (d*d)*0.5;
	 yy = 1.0-((c*c*0.5)+(d*d*0.5));
	 	
	 int hx1=(int(floor(electrons[i].x))-1+NC)%NC;
	 int hx2=(int(floor(electrons[i].x))+1)%NC;
	 int hx3=int(floor(electrons[i].x));
	 
	 int hy1=(int(floor(electrons[i].y))-1+NC)%NC;
	 int hy2=(int(floor(electrons[i].y))+1)%NC;
	 int hy3=int(floor(electrons[i].y));
	 
	 	
	 densTSC[hx1][hy1]-=lx*ly*W;
	 densTSC[hx1][hy2]-=lx*ry*W;
	 densTSC[hx1][hy3]-=lx*yy*W;
	
	 densTSC[hx2][hy1]-=rx*ly*W;
	 densTSC[hx2][hy2]-=rx*ry*W;
	 densTSC[hx2][hy3]-=rx*yy*W;
		
	 densTSC[hx3][hy1]-=xx*ly*W;
	 densTSC[hx3][hy2]-=xx*ry*W;
	 densTSC[hx3][hy3]-=xx*yy*W;
	}
}

void MFT (complex <double> dens[NC][NC])
{
 int l;
 complex<double> nsuum [NC];
 complex<double> densf [NC][NC];
 for ( int k = 0; k < NC; k++)
	{
	 for ( l =0; l < NC; l++)
		{
		 nsuum [l]=0;
		 for (int n = 0; n < NC; n++)
			{
			 nsuum [l] = nsuum [l] +
			 + dens[k][n]*(std::exp(-2.0*pi*l*n*ii/NoC));
			}    
		}
	 for (l = 0; l<NC;l++)
		{
		 densf [k][l]=nsuum [l];
		}
	}
	
 for ( int k = 0; k < NC; k++)
	{
	 for ( l =0; l < NC; l++)
		{
		 nsuum [l]=0;
		 for (int n = 0; n < NC; n++)
			{
			 nsuum [l] = nsuum [l] + 
			 +densf [n][k]*(exp(-2.0*pi*l*n*ii/NoC));
			}
		}
	 for (l = 0; l<NC;l++)
		{
		 densf [l][k]=nsuum [l];
		}
	}
	
 complex <double> G = 0.0;
	
 for ( int k = 0; k < NC; k++)
	{
	 for ( int l = 0; l < NC; l++)
		{
		 G=1.0 + 2.0*(cos(2.0*pi*k/NoC)-1.0) +
		 + 2.0*(cos(2.0*pi*l/NoC)-1.0);
		 phif[k][l] = densf[k][l]/G;
		}
	}

//inverted MFT
 for ( int k = 0; k < NC; k++)
	{
	 for ( l =0; l < NC; l++)
		{
		 nsuum [l]=0;
		 for (int n = 0; n < NC; n++)
			{
			 nsuum [l] = nsuum [l] + 
			 +phif[k][n]*(std::exp(2.0*pi*l*n*ii/NoC));
			}    
		}
	 for (l = 0; l<NC;l++)
		{
		 phi [k][l]=nsuum [l]/NoC;
		}
	}
	
 for ( int k = 0; k < NC; k++)
	{
	 for ( l =0; l < NC; l++)
		{
		 nsuum [l]=0;
		 for (int n = 0; n < NC; n++)
			{
			 nsuum [l] = nsuum [l] +
			 + phi [n][k]*(exp(2.0*pi*l*n*ii/NoC));
			}
		}
	 for (l = 0; l<NC;l++)
		{
		 phi [l][k]=nsuum [l]/NoC;
		}
	}
//computing_field_components_in_cells
 for ( int k = 0; k < NC; k++)
	{
	 for ( int l = 0; l < NC; l++)
		{			
		 if (k==NC) Ex[k][l]=-(phi[0][l]-phi[k][l])/H;
		 if (l==NC) Ey[k][l]=-(phi[k][0]-phi[k][l])/H;
		 else 
		 Ex[k][l]=-(phi[k+1][l]-phi[k][l])/H;
		 Ey[k][l]=-(phi[k][l+1]-phi[k][l])/H;
		}
	}
}
//particle_mover
void MoveParticles (void)
{		
 int j,k;
 for (int i = 0; i < NP; i++)
	{
	 electrons[i].x = electrons[i].x + electrons[i].vx*t;
	 while (electrons[i].x > L) electrons[i].x-=L;
	 while (electrons[i].x < 0) electrons[i].x+=L;
	 j = floor(electrons[i].x);
		
	 electrons[i].y = electrons[i].y + electrons[i].vy*t;
	 while (electrons[i].y > L) electrons[i].y-=L;
	 while (electrons[i].y < 0) electrons[i].y+=L;
	 k = floor(electrons[i].y);
	
	 electrons[i].vy = electrons[i].vy+Ey[j][k].real()*t;
	 electrons[i].vx = electrons[i].vx+Ex[j][k].real()*t;
		
	}
}

int main (int argc, char **argv)
{
 double step;
 FILE *fout = fopen ("output.dat", "wt");
 Init();
	
 for (step = t; step < 10.0; step+=t)
	{
	 double vpar = 0.0;
	 double avpar = 0.0;
	 double h = 0.0;
	 double ah = 0.0;
	 for ( int i = 0; i<NP ; i++)
		{
		 TSCmethod();
		 MFT(densTSC);
		 MoveParticles();
		 electrons[i].v = electrons[i].vx*electrons[i].vx + 
		 +electrons[i].vy*electrons[i].vy;
		 electrons[i].vpar = (electrons[i].vx * electrons[i].vx0 +
		 + electrons[i].vy*electrons[i].vy0) / electrons[i].v0;
		 vpar=vpar+electrons[i].vpar ;
		 h = h + (electrons[i].v - electrons[i].v0*electrons[i].v0);
		}
	 avpar=vpar/(NP);
	 ah= h/(2*NP);
	 cout <<endl << "t="<< step;
     fprintf (fout, "%f %f %f\n", step, avpar, ah);
	}
 fclose (fout);
}	
