rm -f covid_attempt2

gfortran -c -O3 -Wall sort.f90
gfortran -c -O3 -Wall covid_attempt2.f90
gfortran -L$PWD sort.o covid_attempt2.o -llapack -o covid_attempt2

./covid_attempt2