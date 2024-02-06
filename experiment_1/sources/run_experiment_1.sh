rm -f farrington

gfortran -c -O3 -Wall sort.f90
gfortran -c -O3 -Wall farrington.f90
gfortran -L$PWD sort.o farrington.o -o farrington

./farrington