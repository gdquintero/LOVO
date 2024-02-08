rm -f covid

gfortran -c -O3 -Wall sort.f90
gfortran -c -O3 -Wall covid.f90
gfortran -L$PWD sort.o covid.o -o covid

./covid