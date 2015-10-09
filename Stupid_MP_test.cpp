
//   Stupid test of openMP

#include <iostream>
#include <fstream>
#include <iomanip>
#include <cmath>
#include <stdlib.h>

#include <stdio.h>
// #include <libomp.h>
//#include <omp.h>

using namespace std;

//     Here we define various functions called by the main program:

double func_6D(double,double,double,double,double,double);
double abs_dist(double,double,double,double,double,double);

//   Main function begins here
int main()
{
  //#pragma omp parallel
  //printf("Hello from thread %d, nthreads %d\n", omp_get_thread_num(), omp_get_num_threads());
  cout << "Hello, world!" << endl;
  return 0;
}  // end of main program
