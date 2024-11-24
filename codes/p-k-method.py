import numpy as np
import matplotlib.pyplot as plt

# Model validation for p-k-method
# Method is validated by comparings results with Hodges & Pierces (Fig. 5.20 & 5.21 [p.229-230])

# Dimensionless coefficients 
a = -1/5
e = -1/10
mu = 20
r = np.sqrt(6/25)
sigma = 2/5
x_theta = e-a

