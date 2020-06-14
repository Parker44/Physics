"""
Created on Thu Feb 21 17:43:53 2019
@author: Parker Lloyd

The following script processes temperature-dependent voltage data to determine the temperature-dependence 
of the London penetration depth for a superconducting slab of YBCO.
"""
import numpy as np
import matplotlib.pyplot as plt

def is_float(string):
    try:
        return float(string)
    except ValueError:
        return False

# Open .dat file
data = []
with open('YBCO_heating.dat', 'r') as f:
    d = f.readlines()
    for i in d:
        k = i.rstrip().split("\t")
        data.append([float(i) if is_float(i) else i for i in k]) 

data = np.array(data, dtype='O')

# Remove header
data = np.delete(data, (0), 0)

# Remove artifacts
for i in range(176, 170, -1):
    data = np.delete(data, (i), 0)
    
for i in range(8, -1, -1):
    data = np.delete(data, (i), 0)

# Extract temperature data
T = [row[1] for row in data]

# Extract voltage data
V = [row[2] for row in data]

# Linear fit data above critical temperature (~92K)
lin_fit = np.polyfit(T[107::], V[107::], 1)
T = np.asarray(T)
V_linear = lin_fit[1] + lin_fit[0]*T[::]

# Subtract linear fit from entire dataset for voltage to be zero above critical temperature
V = np.subtract(V, V_linear)

# Normalize voltage/magnetization
maxV = np.amax(np.absolute(V))
maxV_index = np.argmax(np.absolute(V))
maxV = maxV * np.sign(V[maxV_index])

V = V[::]/maxV

# Determine London penetration depth as a function of temperature

t = 0.0013  # slab thickness
d0 = 3e-7  # penetration depth at ~80K

d = (1 - V[::]*(1 - (2*d0)/t))*(t/2)

# Take points within Tc - 10K
Tc = T[107]

d = d[17:107]
T = T[17:107]

x = (Tc - T)/Tc
y = d

# Plot log-log
plt.loglog(x, y, 'o')
plt.show()
