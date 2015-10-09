//   Project 3 b) 
//   A second approach to calculating the 6D integral numerically
//   with Gauss-Laguerre quadrature for the r integrals combinded
//   with Gauss-Legendre quadrature for the other integrals in polar
//   coordinates.

#include <cmath>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <cmath>
#include <stdlib.h>
#include <stdio.h>
#define EPS 3.0e-14
#define MAXIT 10
#define ZERO 1.0E-10

using namespace std;
ofstream ofile;

//     Here we define various functions called by the main program:
void gauss_laguerre(double *, double *, int, double);
void gauleg(double, double, double *, double *, int);
double gammln(double);
double r12(double, double, double, double, double, double);
double Integrand(double, double, double, double, double, double);

//   Main function begins here
int main()
{
  // Initial read in of some numbers:
  int n;
  int A[] = {10,15,20,25,30,35,40}; // Values of n to be tested 
  double const  pi = 3.14159265359;
  double alf = 2.0; // the power of r in the integrals; comes from the Jacobi determinant.
  double int_gauss_laguerre = 0.0;
  double exact_integral = 5*pi*pi/256;
  double calculation_time, relatative_error;
  double t1,t2,t3,t4,t5,t6;
  string outfilename = "Gauss_Laguerre_table.txt";

  ofile.open(outfilename);
  ofile << setiosflags(ios::showpoint | ios::uppercase);
  ofile << " Integration in polar coordinates and same resolution in all dimensions. " << endl;
  ofile << " n:        Result with Gauss-Legendre:      Exact result:       Relative error:     Calculation time [s]:" << endl;
  // loop over all values of n:
  for (int l=0; l < 7; l++) {

  int_gauss_laguerre = 0.0;
  calculation_time = 0.0;

  n = A[l];

  double *rgl = new double [n+1];
  double *wgl = new double [n+1];
  gauss_laguerre(rgl, wgl, n, 2.0); // r integrals

  double *theta = new double [n];
  double *w1 = new double [n];
  gauleg(0.0, pi, theta, w1, n);  // theta integrals

  double *phi = new double [n];
  double *w2 = new double [n];
  gauleg(0.0, 2.0*pi, phi, w2, n);  // phi integrals

  clock_t start, finish;
  start = clock();
  for (int i=1; i <= n; i++) {
    t1 = wgl[i];
    for (int j=1; j <= n; j++) {
      t2 = t1*wgl[j];
      for (int k=0; k < n; k++) {
        t3 = t2*w1[k];
        for (int l=0; l < n; l++) {
          t4 = t3*w1[l];
          for (int m=0; m < n; m++) {
            t5 = t4*w2[m];
            for (int m2=0; m2 < n; m2++) {
              t6 = t5*w2[m2];
              int_gauss_laguerre += t6*Integrand(rgl[i],rgl[j],theta[k],theta[l],phi[m],phi[m2]);
            }
          }
        }
      }
    }
  }
  finish = clock();
  calculation_time = (finish - start)/(double)CLOCKS_PER_SEC;
  relatative_error = fabs(int_gauss_laguerre - exact_integral)/exact_integral;

  ofile << setw(5) << setprecision(5) << n;
  ofile << setw(25) << setprecision(10) << int_gauss_laguerre;
  ofile << setw(25) << setprecision(10) << exact_integral;
  ofile << setw(20) << setprecision(4) << relatative_error;
  ofile << setw(20) << setprecision(4) << calculation_time << endl;

  delete [] rgl;
  delete [] wgl;
  delete [] theta;
  delete [] phi;
  delete [] w1;
  delete [] w2;
  }
  ofile.close();

  return 0;
}  // end of main program

double Integrand(double u1, double u2, double theta1, double theta2, double phi1, double phi2) {
  double r_value = r12(u1, u2, theta1, theta2, phi1, phi2);
  double epsilon = 1.0E-5;
  if (r_value <= epsilon) { return 0; }
  else {
      double value = sin(theta1)*sin(theta2)/(pow(4.0,5)*r_value);
      return value;
  }
}

double r12(double u1, double u2, double theta1, double theta2, double phi1, double phi2) {
  double cos_val = cos(theta1)*cos(theta2) + sin(theta1)*sin(theta2)*cos(phi1-phi2);
  return sqrt(u1*u1 + u2*u2 - 2.0*u1*u2*cos_val);
}

/*
** The function
**              gauleg()
** takes the lower and upper limits of integration x1, x2, calculates
** and return the abcissas in x[0,...,n - 1] and the weights in w[0,...,n - 1]
** of length n of the Gauss--Legendre n--point quadrature formulae.
*/
void gauleg(double x1, double x2, double x[], double w[], int n)
{
  int         m,j,i;
  double      z1,z,xm,xl,pp,p3,p2,p1;
  double      const  pi = 3.14159265359;
  double      *x_low, *x_high, *w_low, *w_high;
  m  = (n + 1)/2;                             // roots are symmetric in the interval
  xm = 0.5 * (x2 + x1);
  xl = 0.5 * (x2 - x1);

  x_low  = x;                                       // pointer initialization
  x_high = x + n - 1;
  w_low  = w;
  w_high = w + n - 1;

  for(i = 1; i <= m; i++) {                             // loops over desired roots
    z = cos(pi * (i - 0.25)/(n + 0.5));
      /*
      ** Starting with the above approximation to the ith root
      ** we enter the main loop of refinement by Newtons method.
      */
      do {
        p1 =1.0;
        p2 =0.0;
        /*
        ** loop up recurrence relation to get the
        ** Legendre polynomial evaluated at x
        */
        for(j = 1; j <= n; j++) {
          p3 = p2;
          p2 = p1;
          p1 = ((2.0 * j - 1.0) * z * p2 - (j - 1.0) * p3)/j;
        }
        /*
        ** p1 is now the desired Legenrdre polynomial. Next compute
        ** ppp its derivative by standard relation involving also p2,
        ** polynomial of one lower order.
        */
        pp = n * (z * p1 - p2)/(z * z - 1.0);
        z1 = z;
        z  = z1 - p1/pp;                   // Newton's method
      } while(fabs(z - z1) > ZERO);
      /*
      ** Scale the root to the desired interval and put in its symmetric
      ** counterpart. Compute the weight and its symmetric counterpart
      */
      *(x_low++)  = xm - xl * z;
      *(x_high--) = xm + xl * z;
      *w_low      = 2.0 * xl/((1.0 - z * z) * pp * pp);
      *(w_high--) = *(w_low++);
  }
} // End_ function gauleg()

void gauss_laguerre(double *x, double *w, int n, double alf)
{
        int i,its,j;
        double ai;
        double p1,p2,p3,pp,z,z1;

        for (i=1;i<=n;i++) {
                if (i == 1) {
                        z=(1.0+alf)*(3.0+0.92*alf)/(1.0+2.4*n+1.8*alf);
                } else if (i == 2) {
                        z += (15.0+6.25*alf)/(1.0+0.9*alf+2.5*n);
                } else {
                        ai=i-2;
                        z += ((1.0+2.55*ai)/(1.9*ai)+1.26*ai*alf/
                                (1.0+3.5*ai))*(z-x[i-2])/(1.0+0.3*alf);
                }
                for (its=1;its<=MAXIT;its++) {
                        p1=1.0;
                        p2=0.0;
                        for (j=1;j<=n;j++) {
                                p3=p2;
                                p2=p1;
                                p1=((2*j-1+alf-z)*p2-(j-1+alf)*p3)/j;
                        }
                        pp=(n*p1-(n+alf)*p2)/z;
                        z1=z;
                        z=z1-p1/pp;
                        if (fabs(z-z1) <= EPS) break;
                }
                if (its > MAXIT) cout << "too many iterations in gaulag" << endl;
                x[i]=z;
                w[i] = -exp(gammln(alf+n)-gammln((double)n))/(pp*n*p2);
        }
}
// end function gaulag

double gammln( double xx)
{
        double x,y,tmp,ser;
        static double cof[6]={76.18009172947146,-86.50532032941677,
                24.01409824083091,-1.231739572450155,
                0.1208650973866179e-2,-0.5395239384953e-5};
        int j;

        y=x=xx;
        tmp=x+5.5;
        tmp -= (x+0.5)*log(tmp);
        ser=1.000000000190015;
        for (j=0;j<=5;j++) ser += cof[j]/++y;
        return -tmp+log(2.5066282746310005*ser/x);
}

// end function gammln
//#undef EPS
//#undef MAXIT
