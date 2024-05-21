import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import os

cwd = os.getcwd()
parent =  os.path.abspath(os.path.join(cwd,os.pardir))
size_img = 0.6
plt.rcParams.update({'font.size': 11})
plt.rcParams['figure.figsize'] = [size_img * 6.4,size_img * 4.8]
plt.rc('text', usetex=True)
plt.rc('font', family='serif')

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

country = countries["it"]

# ind_excel = 475
ind_excel = 218

df = pd.read_excel(parent+"/data/"+country)

initial_date = ind_excel - 2
n_train = 27
n_test = 3
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


l = plt.plot(np.linspace(1,n_train,n_train),data,":o",color="darkgreen")
plt.setp(l,'markersize',4)

plt.xlabel("Days")
plt.ylabel("Deaths per million")
plt.xticks(np.arange(5,26,5))
plt.yticks(np.arange(0.1,0.55,0.1))
plt.tick_params(axis='both',direction='in')
plt.savefig(parent+"/images/image.pdf",bbox_inches="tight")
# plt.show()
plt.close()
