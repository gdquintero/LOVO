export ALGENCAN=/home/gustavo/algencan-4.0.0

rm -f covid_logistic

gfortran -c -O3 $ALGENCAN/sources/blas/dgemm.f
gfortran -c -O3 $ALGENCAN/sources/blas/dgemv.f
gfortran -c -O3 $ALGENCAN/sources/blas/dtpmv.f
gfortran -c -O3 $ALGENCAN/sources/blas/dtpsv.f
gfortran -c -O3 $ALGENCAN/sources/blas/idamax.f
gfortran -c -O3 $ALGENCAN/sources/blas/lsame.f
gfortran -c -O3 $ALGENCAN/sources/blas/xerbla.f

ar rcs libblas.a dgemm.o dgemv.o dtpmv.o dtpsv.o idamax.o lsame.o xerbla.o

gfortran -c -O3 $ALGENCAN/sources/hsl/hsl_zd11d.f90
gfortran -c -O3 $ALGENCAN/sources/hsl/hsl_ma57d.f90
gfortran -c -O3 $ALGENCAN/sources/hsl/ma57ad.f
gfortran -c -O3 $ALGENCAN/sources/hsl/mc34ad.f
gfortran -c -O3 $ALGENCAN/sources/hsl/mc47ad.f
gfortran -c -O3 $ALGENCAN/sources/hsl/mc59ad.f
gfortran -c -O3 $ALGENCAN/sources/hsl/mc64ad.f
gfortran -c -O3 $ALGENCAN/sources/hsl/mc21ad.f
gfortran -c -O3 $ALGENCAN/sources/hsl/mc71ad.f
gfortran -c -O3 $ALGENCAN/sources/hsl/fakemetis.f

ar rcs libhsl.a hsl_zd11d.o hsl_ma57d.o ma57ad.o mc34ad.o mc47ad.o mc59ad.o mc64ad.o mc21ad.o mc71ad.o fakemetis.o

gfortran -c -O3 -Wall $ALGENCAN/sources/algencan/lss.f90
gfortran -c -O3 -Wall $ALGENCAN/sources/algencan/gencan.f90

ar rcs libgencan.a lss.o gencan.o

gfortran -c -O3 -Wall sort.f90
gfortran -c -O3 -Wall covid_logistic.f90
gfortran -L$PWD sort.o covid_logistic.o -lgencan -lhsl -lblas -o covid_logistic

./covid_logistic

