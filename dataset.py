import pandas as pd

countries = ["ar.xlsx","br.xlsx","co.xlsx","es.xlsx","it.xlsx","pa.xlsx","uk.xlsx","us.xlsx"]
country = countries[1]

df = pd.read_excel("data/"+country)

ind_excel_init  = 91
ind_excel_end   = 365
init_date       = ind_excel_init - 2
end_date        = ind_excel_end - 2

with open("data/br.txt","w") as f:
    for i in range(ind_excel_end-ind_excel_init+1):
        x = df["new_deaths_smoothed_per_million"][i+init_date]

        if pd.isna(x) == True:
            f.write("%f\n" % 0.0)
        else:
            f.write("%f\n" % x)
