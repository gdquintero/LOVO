import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import os

cwd = os.getcwd()
parent =  os.path.abspath(os.path.join(cwd,os.pardir))

plt.rcParams.update({'font.size': 13})
plt.rc('text', usetex=True)
plt.rc('font', family='serif')

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

ind_excel = 475
# ind_excel = 218

df = pd.read_excel(parent+"/data/"+country)

initial_date = ind_excel - 2
n_train = 25
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


plt.plot(np.linspace(1,n_train,n_train),data,":o",color="darkgreen")
plt.xlabel("Days",fontdict=font)
plt.ylabel("Deaths per million",fontdict=font)
# plt.savefig(parent+"/images/image.pdf",bbox_inches="tight")
# plt.show()
plt.close()

def plot_data():
    for i in range(4):
        plt.plot(np.linspace(i*25+1,(i+1)*25,25),data[i*25:i*25+25],":o",color="darkgreen")
        plt.xlabel("Days",fontdict=font)
        plt.ylabel("Deaths per million",fontdict=font)
        plt.savefig(parent+"/images/image"+str(i+1)+".pdf",bbox_inches="tight")
        plt.close()

# plot_data()