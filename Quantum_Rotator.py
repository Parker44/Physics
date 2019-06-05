"""
Created on Thu Feb 21 17:43:53 2019
@author: Parker Lloyd

This script models a heteronuclear quantum rotator and numerically calculates the partition 
function and heat capacity. The plot axes are transformed into dimensionless quantities.
"""
import math
import mpmath
import numpy as np
import matplotlib.pyplot as plt

# Initialize variables
x = np.linspace(0, 3, 50)
Z = []

# Determine partition function for range of temperatures
# Heat capacity will become constant at high temperatures
for i in range(0, 50):
    try:
        Z.append(mpmath.nsum(lambda j: (2*j + 1)*math.exp(-j*(j+1)/x[i]), [0, mpmath.inf]))
    except ZeroDivisionError:
        continue

# Calculate heat capacity
dZdx = np.diff(Z) / np.diff(x[1::])
dZdB = -np.power(x[2::], 2)*(dZdx[::])

U = -np.power(Z[1::], -1)*(dZdB[::])
C = np.diff(U[::]) / np.diff(x[2::])

plt.plot(x[3::], C, 'r-')
plt.show()
