import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import os

font = {'family': 'serif',
        'color':  'black',
        'weight': 'normal',
        'size': 16,
        }

countries = {
    "ar" : "ar.xlsx",
    "br" : "br.xlsx",
    "co" : "co.xlsx",
    "es" : "es.xlsx",
    "it" : "it.xlsx",
    "pa" : "pa.xlsx",
    "uk" : "uk.xlsx",
    "us" : "us.xlsx"
}

country = countries["br"]

cwd = os.getcwd()
parent =  os.path.abspath(os.path.join(cwd,os.pardir))

df = pd.read_excel(parent+"/data/"+country)

ind_excel_init  = 86
# ind_excel_init  = 547
ind_excel_end   = ind_excel_init + 99 
init_date       = ind_excel_init - 2
end_date        = ind_excel_end - 2
n = ind_excel_end-ind_excel_init + 1
data = np.empty(n)

with open(parent+"/data/covid_mixed.txt","w") as f:
    # f.write("%i\n" % n)
    for i in range(n):
        x = df["new_deaths_smoothed_per_million"][i+init_date]
        data[i] = x

        if pd.isna(x) == True:
            f.write("%f\n" % 0.0)
        else:
            f.write("%f\n" % x)

t = np.linspace(1,100,100)

# ind = 80
# ind = 100

plt.plot(t,data[:100],"--o",color="darkgreen")
plt.xlabel("Days",fontdict=font)
plt.ylabel("New deaths per million",fontdict=font)
plt.savefig(parent+"/images/data.pdf",bbox_inches = "tight")
plt.show()