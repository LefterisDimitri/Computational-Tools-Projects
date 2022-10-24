#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <omp.h>


double f(double x)
{
  return x*exp(2.0*x);
}
 

int main()
{
  int i, N = 1e9;
  double x_0, x_n, h, x;
  double sum, I;
  
  x_0 = 0.0;
  x_n = 3.0;
  h=(x_n-x_0)/N;
  sum = 0.0;
   
  //omp_set_num_threads(4);
    
  double start_time = omp_get_wtime(); 
  #pragma omp parallel for shared(x_0,x_n,h,N) private(x,i) reduction(+:sum) default(none)
  for(i=1;i<N;i++)
  {
    x = x_0 + i*h;
    
    if(i%2==0)
    {
      sum = sum + 2.0*f(x);
    }
    else
    {
      sum = sum + 4.0*f(x);
    }
  }
  double end_time = omp_get_wtime();

  printf("timer gave us: t = %lf sec\n", end_time - start_time);
 
  I = (h/3.0)*(f(x_0) + sum + f(x_n));

  printf("Integral: I = %lf\n", I);
return 0;
}
