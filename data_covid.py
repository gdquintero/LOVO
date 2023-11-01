import pandas as pd

# Dataframe of COVID-19 in Brazil
df = pd.read_excel("data/dados_italia.xlsx")

# ind_excel = 122 # 2 de mayo 2020
# ind_excel = 548 # 2 de julio 2021
ind_excel = 2 # 1 de octubre 2020 UK

initial_date = ind_excel - 2
total_days = 70
n_train = 60
n_test = 10

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
