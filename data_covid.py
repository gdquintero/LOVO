import pandas as pd

df = pd.read_excel("dados_full_brasil.xlsx")

ind_excel = 122 # 2 de mayo 2020

initial_date = ind_excel - 2
total_days = 40

with open("output/data_covid.txt","w") as f:
    f.write("%i\n" % total_days)
    for i in range(initial_date,initial_date + total_days):
        x = df["new_deaths_smoothed_per_million"][i]

        if pd.isna(x) == True:
            f.write("%f\n" % 0.0)
        else:
            f.write("%f\n" % x)
            
print
