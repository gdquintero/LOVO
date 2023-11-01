import pandas as pd

# ind_excel = 122 # 2 de mayo 2020
# ind_excel = 305 # 1 de nov 2020 Colombia

country = "data/spain.xlsx"
ind_excel = 731 # 1 de ene 2022 España

df = pd.read_excel(country)

initial_date = ind_excel - 2
n_train = 60
n_test = 10
total_days = n_train + n_test

with open("data/covid_train.txt","w") as f:
    f.write("%i\n" % n_train)
    for i in range(initial_date,initial_date + n_train):
        x = df["new_deaths_smoothed_per_million"][i]

        if pd.isna(x) == True:
            f.write("%f\n" % 0.0)
        else:
            f.write("%f\n" % x)

with open("data/covid_test.txt","w") as f:
    f.write("%i\n" % n_test)
    for i in range(initial_date + n_train,initial_date + total_days):
        x = df["new_deaths_smoothed_per_million"][i]

        if pd.isna(x) == True:
            f.write("%f\n" % 0.0)
        else:
            f.write("%f\n" % x)
