rm -f main

gfortran -c -O3 -Wall sort.f90
gfortran -c -O3 -Wall main.f90
gfortran -L$PWD sort.o main.o -llapack -o main

./main