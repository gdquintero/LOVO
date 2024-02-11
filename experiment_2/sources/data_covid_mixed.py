import pandas as pd
import os

countries = ["ar","br","co","es","it","pa","uk","us"]
country = countries[1]

cwd = os.getcwd()
parent =  os.path.abspath(os.path.join(cwd,os.pardir))

df = pd.read_excel(parent+"/data/"+country+".xlsx")


ind_excel_init  = 91 # 01-04-2020
ind_excel_end   = 365 # 31-12-2020
init_date       = ind_excel_init - 2
end_date        = ind_excel_end - 2
n = ind_excel_end-ind_excel_init + 1


with open(parent+"/data/covid_mixed.txt","w") as f:
    f.write("%i\n" % n)
    for i in range(n):
        x = df["new_deaths_smoothed_per_million"][i+init_date]

        if pd.isna(x) == True:
            f.write("%f\n" % 0.0)
        else:
            f.write("%f\n" % x)