rm -f covid_attempt1

gfortran -c -O3 -Wall sort.f90
gfortran -c -O3 -Wall covid_attempt1.f90
gfortran -L$PWD sort.o covid_attempt1.o -llapack -o covid_attempt1

./covid_attempt1