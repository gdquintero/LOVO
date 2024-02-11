rm -f covid_mixed

gfortran -c -O3 -Wall sort.f90
gfortran -c -O3 -Wall covid_mixed.f90
gfortran -L$PWD sort.o covid_mixed.o -o covid_mixed

./covid_mixed