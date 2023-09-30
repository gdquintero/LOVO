export ALGENCAN=/home/gustavo/algencan-4.0.0

rm -f lovo_covid_cubic_copy

# gfortran -c -O3 -w -fcheck=all -Wall -I$ALGENCAN/sources/algencan/inc lovo.f90 

gfortran -g -w -fcheck=all \
-I$ALGENCAN/sources/algencan/inc lovo_covid_cubic_copy.f90 \
-L$ALGENCAN/sources/algencan/lib -lalgencan -L$ALGENCAN/sources/hsl/lib -lhsl \
-L$ALGENCAN/sources/blas/lib -lblas sort.o subset.o -o lovo_covid_cubic_copy

./lovo_covid_cubic_copy

