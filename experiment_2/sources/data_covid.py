import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import os

cwd = os.getcwd()
parent =  os.path.abspath(os.path.join(cwd,os.pardir))

plt.rcParams.update({'font.size': 13})
plt.rc('text', usetex=True)
plt.rc('font', family='serif')

# font = {'family': 'serif',
#         'color':  'darkred',
#         'weight': 'normal',
#         'size': 16,
#         }

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
# ind_excel = 100 # 10 abr 2020 desorden panama
# ind_excel = 122 # 2 de mayo 2020 subida brasil - bajada uk
# ind_excel = 171 # 20 de junio usa
# ind_excel = 166 # 15 jun 2020 estados unidos con outliers
# ind_excel = 216 # 4 de agosto italia 

ind_excel = 476

df = pd.read_excel(parent+"/data/"+country)

initial_date = ind_excel - 2
n_train = 5
n_test = 5
total_days = n_train + n_test  
data = np.zeros(n_train)

with open(parent+"/data/covid.txt","w") as f:
    f.write("%i\n" % n_train)
    f.write("%i\n" % n_test)
    j = 0

    for i in range(initial_date,initial_date + total_days):
        x = df["new_deaths_smoothed_per_million"][i]

        if pd.isna(x) == True:
            f.write("%f\n" % 0.0)
        else:
            f.write("%f\n" % x)

        if j <= n_train - 1:
            data[j] = x
    
        j += 1

print(data)

plt.plot(np.linspace(1,n_train,n_train),data,"--o",color="darkgreen")
plt.xlabel("Days")
plt.ylabel("New deaths per million")
plt.savefig(parent+"/images/data_ita.pdf",bbox_inches = "tight")
plt.show()
