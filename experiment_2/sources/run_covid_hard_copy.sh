rm -f covid_hard_copy

gfortran -c -O3 -Wall sort.f90
gfortran -c -O3 -Wall covid_hard_copy.f90
gfortran -L$PWD sort.o covid_hard_copy.o -llapack -o covid_hard_copy

./covid_hard_copy