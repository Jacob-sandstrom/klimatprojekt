import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

pco20 = 278.05158

data = pd.read_csv('radiativeForcingRCP45.csv')
s = 1
year = data["Time (year)"].values   
# rfdata = data["RF CO2 (W/m2)"].values
rfdata = data["RF aerosols (W/m2)"].values * s
rfdata += data["RF other than CO2 and aerosols (W/m2)"].values 


concentration_data = pd.read_csv('koncentrationerRCP45.csv')
year_concentration = concentration_data["Time (year)"].values
co2_concentration = concentration_data["CO2ConcRCP45 (ppm CO2)"].values


rfmodel = []
for i in range(len(year_concentration)):
    rfmodel.append(5.35*np.log(co2_concentration[i]/pco20))

rfdata += rfmodel





dT1 = 0
dT2 = 0

lambd = 0.8 #0.5-1.3
k = 0.5 # 0.2-1
RF = 1
c = 4186
p = 1020
h = 50
c1 = c*p*h/31556952
d = 2000
c2 = c*p*d/31556952

def delt1(c1, rf, lambd, t1, t2, k):
    return (rf-t1/lambd-k*(t1-t2))/c1

def delt2(c2, t1, t2, k):
    return k* (t1-t2)/c2


t1data = []
for i in range(len(rfdata)):
    ddT1 = delt1(c1, rfdata[i], lambd, dT1, dT2, k)
    ddT2 = delt2(c2, dT1, dT2, k)

    dT1 += ddT1
    dT2 += ddT2

    t1data.append(dT1)

    # if dT1 > (1-np.exp(-1))*RF*lambd and dT2 > (1-np.exp(-1))*RF*lambd:
    #     break


# Read the CSV file
file_path = 'NASA_GISS.csv'
data = pd.read_csv(file_path)

# Display the first few rows to understand the structure

# Remove the first row from the data
data = data.iloc[1:]
print(data.head())

years = data['Land-Ocean Temperature Index (C)'].values
temperature_anomalies = data['Unnamed: 1'].values
# Convert the years to a numeric format
years = pd.to_numeric(years, errors='coerce')

temperature_anomalies = pd.to_numeric(temperature_anomalies, errors='coerce')

# print(years)
# print(year)

tempaverage = np.mean(t1data[186:216])
print(year[186:216])
# t1data = np.subtract(t1data, tempaverage)

temperature_anomalies = np.add(temperature_anomalies, tempaverage)

year = year[:260]
rfdata = rfdata[:260]
t1data = t1data[:260]

# t1data[196:216]

# Plot the data
plt.figure(figsize=(10, 6))
plt.plot(years, temperature_anomalies, label='Observed', color='red')
plt.plot(year, t1data, label='Model', color='blue')
# plt.plot(year, rfdata, label='Observed', color='green')
plt.xlabel('Year')
plt.ylabel('Temperature Anomaly (°C)')
plt.title('Global Temperature Anomalies Over Time (NASA GISS)')
# plt.grid(True)
plt.legend()
plt.hlines(0, 1765, 2024, colors='black', linestyles='dashed')


plt.show()