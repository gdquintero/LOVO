import pandas as pd
countries = ["ar.xlsx","br.xlsx","co.xlsx","es.xlsx","it.xlsx","pa.xlsx","uk.xlsx","us.xlsx"]
country = countries[1]
# ind_excel = 100 # 10 abr 2020 desorden panama
# ind_excel = 122 # 2 de mayo 2020 subida brasil - bajada uk
# ind_excel = 171 # 20 de junio usa
# ind_excel = 166 # 15 jun 2020 estados unidos con outliers
# ind_excel = 216 # 4 de agosto italia 

ind_excel  = 91 # 01-04-2020


df = pd.read_excel("data/"+country)

initial_date = ind_excel - 2
n_train = 10
n_test = 5
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


