import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import scipy.optimize as opt

data = pd.read_csv('radiativeForcingRCP45.csv')
year = data["Time (year)"].values   


# concentration_data = pd.read_csv('koncentrationerRCP45.csv')
# year_concentration = concentration_data["Time (year)"].values
# co2_concentration = concentration_data["CO2ConcRCP45 (ppm CO2)"].values


# rfmodel = []
# for i in range(len(year_concentration)):
#     rfmodel.append(5.35*np.log(co2_concentration[i]/pco20))


utslapp = pd.read_csv('utslappRCP45.csv')

# Scenario 1

# Scenario 2

# Scenario 3


beta = 0.35

B_0 = np.array([600, 600, 1500])
F_0 = np.array([[0, 60, 0],
                [15, 0, 45],
                [45, 45, 0]])

Alpha = np.divide(F_0.T, B_0).T

B = B_0

A = np.array([0.113, 0.213, 0.258, 0.273, 0.143])
tau_0 = np.array([2, 12.2, 50.4, 243.3, np.inf])

k = 3.06*10**(-3)

def Mt6(t, U):
    tau = np.multiply(tau_0, 1 + k*sum(U[:t]))
    I = lambda t_hat: np.sum(np.multiply(A, np.exp(np.divide(-t_hat, tau))))
    return B_0[0] + sum(I(t-t_tilde)*U[t_tilde] for t_tilde in range(t+1))

def utslappscen(U):

    B_array = np.zeros((len(U), 3))
    B_array[0] = B_0
    B = B_0

    U_array = np.zeros(len(U))


    for i in range(len(U)):
        NPP = F_0[0,1]*(1+beta*np.log(B[0]/B_0[0]))

        U_array[i] = U[i] - NPP + Alpha[2,0]*B[2] + Alpha[1,0]*B[1]

        B = np.array([Mt6(i, U_array), 
                    NPP - Alpha[1,2]*B[1] - Alpha[1,0]*B[1] + B[1],
                    Alpha[1,2]*B[1] - Alpha[2,0]*B[2] + B[2]])
        
        B_array[i] = B

    pco20 = 278.05158

    
    s = 1.20110178
    
    # rfdata = data["RF CO2 (W/m2)"].values
    data = pd.read_csv('radiativeForcingRCP45.csv')
    rfdata = data["RF aerosols (W/m2)"].values * s
    rfdata += data["RF other than CO2 and aerosols (W/m2)"].values 

    rfmodel = 5.35*(np.log(B_array[:,0] * 0.469/pco20))
    rfdata = np.add(rfdata, rfmodel)

    dT1 = 0
    dT2 = 0

    lambd = 0.8 #0.5-1.3
    kap = 0.70137457 # 0.2-1
    RF = 1
    c = 4186
    p = 1020
    h = 50
    c1 = c*p*h/31556952
    d = 2000
    c2 = c*p*d/31556952

    t1data = []
    for i in range(len(rfdata)):
        ddT1 = delt1(c1, rfdata[i], lambd, dT1, dT2, kap)
        ddT2 = delt2(c2, dT1, dT2, kap)

        dT1 += ddT1
        dT2 += ddT2

        t1data.append(dT1)
    
    return t1data



def delt1(c1, rf, lambd, t1, t2, kap):
    return (rf-t1/lambd-kap*(t1-t2))/c1

def delt2(c2, t1, t2, kap):
    return kap* (t1-t2)/c2


# Define scenarios
UtslappScenarios = np.zeros((3, len(utslapp['CO2 Emissions  (CO2 GtC/yr)'])))
UtslappScenarios[0,:260] = utslapp['CO2 Emissions  (CO2 GtC/yr)'][:260]
UtslappScenarios[1,:260] = utslapp['CO2 Emissions  (CO2 GtC/yr)'][:260]
UtslappScenarios[2,:260] = utslapp['CO2 Emissions  (CO2 GtC/yr)'][:260]

UtslappScenarios[0,259:306] =  np.linspace(UtslappScenarios[0,259], 0, 47)
UtslappScenarios[1,259:] = [UtslappScenarios[1,259] for _ in range(len(UtslappScenarios[1,259:]))]
UtslappScenarios[2,259:336] = np.linspace(UtslappScenarios[0,259], 2*UtslappScenarios[0,259], 77)
UtslappScenarios[2,335:] = [UtslappScenarios[2,335] for _ in range(len(UtslappScenarios[2,335:]))]


T1Scenarios = np.zeros((3, len(utslapp['CO2 Emissions  (CO2 GtC/yr)'])))

year = year[:435]
print(UtslappScenarios[0,259:306])



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




# t1data[196:216]

# Plot the data
plt.figure(figsize=(15, 6))
for i in range(3):
    T1Scenarios[i] = utslappscen(UtslappScenarios[i])

tempaverage = np.mean(T1Scenarios[0,186:216])
# print(year[186:216])
# t1data = np.subtract(t1data, tempaverage)

temperature_anomalies = np.add(temperature_anomalies, tempaverage)


plt.plot(year, T1Scenarios[0,:435], label='Scenario 1', color='blue')
plt.plot(year, T1Scenarios[1,:435], label='Scenario 2', color='green')
plt.plot(year, T1Scenarios[2,:435], label='Scenario 3', color='orange')

plt.plot(years, temperature_anomalies, label='Observed', color='red')
# plt.plot(year, t1data, label='Model', color='blue')
# plt.plot(year, rfdata, label='Observed', color='green')
plt.xlabel('Year')
plt.ylabel('Temperature Anomaly (°C)')
plt.title('Global Temperature Anomalies Over Time (NASA GISS)')
# plt.grid(True)
plt.legend()
plt.hlines(0, 1765, 2200, colors='black', linestyles='dashed')

plt.xlim(1880, 2200)

plt.show()