rm -f covid_hard

gfortran -c -O3 -Wall -w sort.f90
gfortran -c -O3 -Wall -w covid_hard.f90
gfortran -L$PWD sort.o covid_hard.o -o covid_hard

./covid_hard