//   Project 3 a)
//   A first approach to calculating the 6D integral numerically
//   with Gauss-Legendre quadrature.

#include <cmath>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <cmath>
#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include <random>
#define EPS 3.0e-14
#define MAXIT 10
#define ZERO 1.0E-10

using namespace std;
ofstream ofile;

//     Here we define various functions called by the main program:
double func_6D(double, double, double, double, double, double);
double abs_dist(double, double, double, double, double, double);

//   Main function begins here
int main()
{
  // Initial read in of some numbers:
  int A[] = {1E5,1E6,1E7,1E8,1E9}; // Values of n to be tested 
  int N;
  double a,b;
  int d = 6; // number of dimension.
  double const  pi = 3.14159265359;
  double x1,x2,y1,y2,z1,z2;
  double exact_integral = 5*pi*pi/256;
  double calculation_time, relative_error, variance;
  double integral_MC, sum_sigma, sigma, f_val;

  string outfilename = "Bruteforce_MC_table.txt";

  cout << "Read in integration limits" << endl;
  cin >> a >> b;
  
  double jacobi_det = pow((b-a),d);
  random_device rd;
  mt19937 gen(rd());
  uniform_real_distribution<> ran0(a, b);

  ofile.open(outfilename);
  ofile << setiosflags(ios::showpoint | ios::uppercase);
  ofile << " Integration with brute force Monte Carlo  " << endl;
  ofile << " Integration limits in each dimension: a = " << a << " and b = " << b << endl;
  ofile << "        N:       MC result:     Exact result:  Relative error: Calculation time [s]:  sigma:" << endl;

  for (int l=0; l < 5; l++) {
  N = A[l];
  integral_MC = 0.0;
  sum_sigma = 0.0;
  sigma = 0.0;

  clock_t start, finish;
  start = clock();
  for (int i=0; i < N; i++) {
    x1 = ran0(gen); y1 = ran0(gen); z1 = ran0(gen);
    x2 = ran0(gen); y2 = ran0(gen); z2 = ran0(gen);
    f_val = func_6D(x1,y1,z1,x2,y2,z2);
    integral_MC += f_val;
    sum_sigma += f_val*f_val;
  } 
  finish = clock();

  integral_MC /= (double) N;
  sum_sigma /= (double) N;
  variance = sum_sigma - integral_MC*integral_MC;
  integral_MC *= jacobi_det;
  sigma = jacobi_det*sqrt(variance/((double) N));
  relative_error = fabs(integral_MC-exact_integral)/exact_integral;
  calculation_time = (finish - start)/(double)CLOCKS_PER_SEC;


  ofile << setw(15) << setprecision(5) << N;
  ofile << setw(15) << setprecision(10) << integral_MC;
  ofile << setw(15) << setprecision(10) << exact_integral;
  ofile << setw(15) << setprecision(4) << relative_error;
  ofile << setw(15) << setprecision(4) << calculation_time;
  ofile << setw(20) << setprecision(5) << sigma << endl;

  }
  ofile.close();

  return 0;
}  // end of main program

// Help function for func_6D.
// It calculates the reciprocal value of the distance between r_1 and r_2.
double abs_dist(double x1, double y1, double z1, double x2, double y2, double z2) {
  double x1_x2 = (x1-x2)*(x1-x2);
  double y1_y2 = (y1-y2)*(y1-y2);
  double z1_z2 = (z1-z2)*(z1-z2);
  double r1_r2 = sqrt( x1_x2 + y1_y2 + z1_z2 );
  return r1_r2;
}

// Six dimensional integrand to be integrated.
// We simply exclude the integration points where the integrand diverges.
double func_6D(double x1, double y1, double z1, double x2, double y2, double z2) {
  double alpha = 2.0;
  double f_val = 0.0;
  double r1 = sqrt(x1*x1 + y1*y1 + z1*z1);
  double r2 = sqrt(x2*x2 + y2*y2 + z2*z2);

  double r1_r2 = abs_dist(x1,y1,z1,x2,y2,z2);
  // Simply skip the parts where the integrand is singular:
  return exp(-2.0*alpha*(r1 + r2)) / r1_r2; 
}
