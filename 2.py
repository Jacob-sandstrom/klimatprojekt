import pandas as pd
import numpy as np
import matplotlib.pyplot as plt


pco20 = 278.05158


data = pd.read_csv('radiativeForcingRCP45.csv')
s = 1
year = data["Time (year)"].values   
rfdata = data["RF CO2 (W/m2)"].values
rfdata += data["RF aerosols (W/m2)"].values
rfdata += data["RF other than CO2 and aerosols (W/m2)"].values * s


concentration_data = pd.read_csv('koncentrationerRCP45.csv')
year_concentration = concentration_data["Time (year)"].values
co2_concentration = concentration_data["CO2ConcRCP45 (ppm CO2)"].values


rfmodel = []
for i in range(len(year_concentration)):
    rfmodel.append(5.35*np.log(co2_concentration[i]/pco20))


plt.plot(year_concentration, rfmodel, label='model', color='blue')
plt.plot(year, rfdata, label='Observed', color='red')
plt.show()