import pandas as pd
import numpy as np
import matplotlib.pyplot as plt


pco20 = 278.05158


data = pd.read_csv('radiativeForcingRCP45.csv')
s = 1
year = data["Time (year)"].values   
rfdata = data["RF CO2 (W/m2)"].values
rfdata += data["RF aerosols (W/m2)"].values * s
rfdata += data["RF other than CO2 and aerosols (W/m2)"].values 


concentration_data = pd.read_csv('koncentrationerRCP45.csv')
year_concentration = concentration_data["Time (year)"].values
co2_concentration = concentration_data["CO2ConcRCP45 (ppm CO2)"].values


rfmodel = []
for i in range(len(year_concentration)):
    rfmodel.append(5.35*np.log(co2_concentration[i]/pco20))


plt.plot(year_concentration, rfmodel, label='model', color='blue')
plt.plot(year, rfdata, label='Observed', color='red')
plt.show()



dT1 = 0
dT2 = 0

lambd = 0.8 #0.5-1.3
k = 0.5 # 0.2-1
RF = 1
c = 4186
p = 1020
h = 50
c1 = c*p*h/1314000
d = 2000
c2 = c*p*d/1314000

def delt1(c1, rf, lambd, t1, t2, k):
    return (rf-t1/lambd-k*(t1-t2))/c1

def delt2(c2, t1, t2, k):
    return k* (t1-t2)/c2



for i in range(100000):
    ddT1 = delt1(c1, RF, lambd, dT1, dT2, k)
    ddT2 = delt2(c2, dT1, dT2, k)

    dT1 += ddT1
    dT2 += ddT2


    if dT1 > (1-np.exp(-1))*RF*lambd and dT2 > (1-np.exp(-1))*RF*lambd:
        break

print(i)
# print(dT1- dT2)
print(dT1, dT2)
print(RF*lambd)

#b) higher lambda => higher equilibrium temp
#b) higher k => faster heat transfer


dT1 = 0
dT2 = 0
oHeatUptake = []
spaceRadiation = []
rf = []
for i in range(200):
    rf.append(RF)
    oHeatUptake.append(k* (dT1-dT2))
    spaceRadiation.append(dT1/lambd)

    ddT1 = delt1(c1, RF, lambd, dT1, dT2, k)
    ddT2 = delt2(c2, dT1, dT2, k)

    dT1 += ddT1
    dT2 += ddT2

plt.plot(oHeatUptake, label='ocean heat uptake', color='blue')
plt.plot(spaceRadiation, label='space radiation', color='red')
plt.plot(rf, label='radiative forcing', color='green')
plt.xlabel('Time (years)')
plt.legend()
plt.title('Ocean heat uptake, space radiation and radiative forcing')
plt.show()

