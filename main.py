import csv
import numpy as np
import pandas
import matplotlib.pyplot as plt




data = pandas.read_csv('utslappRCP45.csv')
print(data.head())


beta = 0.35
# beta = 0.6 #Bästa fit

b0 = [600, 600, 1500]
b = [600, 600, 1500]

f = [
    [0, 60, 0],
    [15, 0, 45],
    [45, 0, 0]
]

alpha = [
    [value / b[0] for value in f[0]],
    [value / b[1] for value in f[1]],
    [value / b[2] for value in f[2]]
]

NPP = 60

# Define the equations
def dB1_dt(b, NPP, U):
    return alpha[2][0] * b[2] + alpha[1][0] * b[1] - NPP + U

def dB2_dt(b, NPP):
    return NPP - alpha[2][1] * b[1] - alpha[0][1] * b[1]

def dB3_dt(b):
    return alpha[1][2] * b[1] - alpha[2][0] * b[2]



year = data["Time (year)"].values
u = data["CO2 Emissions  (CO2 GtC/yr)"].values

co2 = []

for i in range(len(year)):
    npp = NPP * (1+beta*(np.log(b[0]/b0[0])))
    dB1 = dB1_dt(b, npp, u[i])
    dB2 = dB2_dt(b, npp)
    dB3 = dB3_dt(b)    

    b = [b[0] + dB1, b[1] + dB2, b[2] + dB3]

    co2.append(b[0]*0.469)



# Read the 'koncentrationerRCP45.csv' file
koncentrationer_data = pandas.read_csv('koncentrationerRCP45.csv')
print(koncentrationer_data.head())


print(koncentrationer_data.columns)
co2conc = koncentrationer_data["CO2ConcRCP45 (ppm CO2)"].values


plt.plot(year, co2conc, label='Observed CO2 Concentration', color='blue')
plt.plot(year, co2, label='Model CO2 Concentration', color='red')    
plt.show()