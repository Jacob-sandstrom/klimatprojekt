import numpy as np
import pandas as pd
#import matplotlib.pyplot as plt

utslapp = pd.read_csv('utslappRCP45.csv')


# TASK 1
beta = 0.35

B_0 = np.array([600, 600, 1500])
F_0 = np.array([[0, 60, 0],
                [15, 0, 45],
                [45, 45, 0]])

Alpha = np.divide(F_0.T, B_0).T

B = B_0
# for i in range(len(utslapp['CO2 Emissions  (CO2 GtC/yr)'])):
#     NPP = F_0[0,1]*(1+beta*np.log(B[0]/B_0[0]))

#     dB = np.array([Alpha[2,0]*B[2] + Alpha[1,0]*B[1] - NPP + utslapp['CO2 Emissions  (CO2 GtC/yr)'][i], 
#                    NPP - Alpha[1,2]*B[1] - Alpha[1,0]*B[1],
#                    Alpha[1,2]*B[1] - Alpha[2,0]*B[2]])
    
#     B = B + dB

#     print(B*0.469)


# TASK 3
A = np.array([0.113, 0.213, 0.258, 0.273, 0.143])
tau_0 = np.array([2, 12.2, 50.4, 243.3, np.inf])

k = 3.06*10**(-3)

def M(t):
    tau = np.multiply(tau_0, 1 + k*sum(utslapp['CO2 Emissions  (CO2 GtC/yr)'][:t]))
    I = lambda t_hat: np.sum(np.multiply(A, np.exp(np.divide(-t_hat, tau))))
    return B[0] + np.sum(I(t-t_tilde)*utslapp['CO2 Emissions  (CO2 GtC/yr)'][t_tilde] for t_tilde in range(t+1))


# TASK 4
# for t in range(len(utslapp['CO2 Emissions  (CO2 GtC/yr)'])):
#     print(M(t)*0.469)


# TASK 6
def Mt6(t, U):
    tau = np.multiply(tau_0, 1 + k*sum(U[:t]))
    I = lambda t_hat: np.sum(np.multiply(A, np.exp(np.divide(-t_hat, tau))))
    return B[0] + np.sum(I(t-t_tilde)*U[t_tilde] for t_tilde in range(t+1))

# B_array = np.zeros((len(utslapp['CO2 Emissions  (CO2 GtC/yr)']) + 1, 3))
# B_array[0] = B_0
B = B_0

U_array = np.zeros(len(utslapp['CO2 Emissions  (CO2 GtC/yr)']) + 1)


for i in range(len(utslapp['CO2 Emissions  (CO2 GtC/yr)'])):
    NPP = F_0[0,1]*(1+beta*np.log(B[0]/B_0[0]))

    U_array[i] = utslapp['CO2 Emissions  (CO2 GtC/yr)'][i] - NPP + Alpha[2,0]*B[2] + Alpha[1,0]*B[1]

    B = np.array([Mt6(i, U_array), 
                   NPP - Alpha[1,2]*B[1] - Alpha[1,0]*B[1] + B[1],
                   Alpha[1,2]*B[1] - Alpha[2,0]*B[2] + B[2]])
    
    # B_array[i] = B

    print(utslapp['Time (year)'][i], B*0.469)