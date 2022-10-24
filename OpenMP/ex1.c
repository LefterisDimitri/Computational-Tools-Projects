#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <omp.h>

int main()
{
  FILE *fp;
  fp = fopen("pde.dat","w");
 
  int i, j;
  int N=50;               /* default value */
  int M=50;               /* default value */
  
  double 
       **uold,
       **unew,
        *x,
        *y,
         dx,
         dy,
         dt,
         t = 0.0,
         start_time,
         end_time;

 

// Allocate memory

  uold = (double**) malloc(N*sizeof(double*));
  for (int i = 0; i < N; i++)
      uold[i] = (double*) malloc(M*sizeof(double));

  unew = (double**) malloc(N*sizeof(double*));
  for (int i = 0; i < N; i++)
      unew[i] = (double*) malloc(M*sizeof(double));

  x = (double*) malloc(N*sizeof(double));
  y = (double*) malloc(M*sizeof(double));
  
  // grid spacings
  dx = (2.0*M_PI-0.0)/(N-1);
  dy = (2.0*M_PI-0.0)/(M-1);
  dt = 0.005;                        //     dt < (5/16)*((2*pi)/(N-1))^2
  
  start_time = omp_get_wtime(); 
  #pragma omp parallel private(i,j) shared(dt,x,dx,y,dy,uold,unew,N,M) reduction(+:t) default(none)
  {
    // x and y grids
    #pragma omp for 
    for(i=0;i<N;i++)
    {
     x[i] = 0.0 + i*dx;
    }

    #pragma omp for 
    for(j=0;j<M;j++)
    {
     y[j] = 0.0 + j*dy;
    }

    // arxikopoihsh pinaka
    #pragma omp for collapse(2)
    for(i=0;i<N;i++)
    { 
       for(j=0;j<M;j++)
       {
          uold[i][j] = 0.0;
          unew[i][j] = 0.0;
       }
    }

    // sunoriakes 
    #pragma omp for 
    for(i=0;i<N;i++)
    {
       uold[i][0] = 0.0;
       uold[i][N-1] = 0.0;
       unew[i][0] = 0.0;
       unew[i][N-1] = 0.0;
    }

    #pragma omp for 
    for(j=0;j<M;j++)
    {
       uold[0][j] = 0.0;
       uold[M-1][j] = 0.0;
       unew[0][j] = 0.0;
       unew[M-1][j] = 0.0;
    }
    
    #pragma omp for collapse(2) 
    for(i=0;i<N;i++)
    {
     for(j=0;j<M;j++)
     {
      uold[i][j] = sin(x[i])*sin(y[j]/2.0);
     }
    }

    // iterative solution
    do
    {
       t += dt; // moving forward in time
       #pragma omp for collapse(2) 
       for(i=1;i<N-1;i++) 
       {
        for(j=1;j<M-1;j++) 
        {
         // ---------- oi sunoriakes sunthikes se kathe xrono ---------- //
         uold[0][j]=0.0;      // u_old(x=0,y,t) = 0
         uold[N-1][j]=0.0;    // u_old(x=2*pi,y,t) = 0
	  
         uold[i][0]=0.0;      // u_old(x,y=0,t) = 0
         uold[i][N-1]=0.0;    // u_old(x,y=2*pi,t) = 0
	  
         unew[0][j]=0.0;     // u_new(x=0,y,t) = 0
         unew[N-1][j]=0.0;   // u_new(x=2*pi,y,t) = 0
         
         unew[i][0]=0.0;     // u_new(x,y=0,t) = 0
         unew[i][N-1]=0.0;   // u_new(x,y=2*pi,t) = 0
	 // ------------------------------------------------------------ //
	 
	 // h anadromikh sxesh ths methodou //
	 unew[i][j] = uold[i][j]+(4.0/5.0)*dt*((uold[i+1][j] - 2.0*uold[i][j] + uold[i-1][j])/(dx*dx) + (uold[i][j+1] - 2.0*uold[i][j] + uold[i][j-1])/(dy*dy));          
        }
       }

       // reset uold for next iteration
       #pragma omp for collapse(2) 
       for(i=1;i<N-1;i++)
       {
          for(j=1;j<M-1;j++)
          { 
             uold[i][j] = unew[i][j]; 
          }
       }
       
    }while(t<2.0*M_PI);
 }
  
  
 end_time = omp_get_wtime();
 double d = end_time - start_time;
 
 printf("\n\n\nt = %lf sec", d);
 
  return 0; 
}

