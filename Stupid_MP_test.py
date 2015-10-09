# Program for running Test of integration methods:

import os

# -I/usr/local/include -I/usr/local/lib/libiomp5.dylib -I/usr/local/bin/clang-omp

os.system('g++ Stupid_MP_test.cpp -o Stupid_MP_test.o -I/usr/local/include -I/usr/local/lib/libiomp5.dylib -I/usr/local/bin/clang-omp -fopenmp')
os.system('./Stupid_MP_test.o')