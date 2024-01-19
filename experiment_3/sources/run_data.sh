rm -f data_cubic

gfortran data_cubic.f90 -o data_cubic

./data_cubic

python3 plot_data.py