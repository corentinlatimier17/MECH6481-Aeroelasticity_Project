import matplotlib.pyplot as plt
import numpy as np
from scipy.linalg import eig

# Model validation for p-method using eigenvalue analysis (see Handout 1 for further informations)
# Method is validated by comparings results with Hodges & Pierces (Fig. 5.3 & 5.4 [p.186-187])

# -> This code has been validated and gives the same flutter speed and flutter frequency than the Hodges & Pierces book 


def identify_flutter_speed(V_vec, roots):
    # IDENTIFY REDUCED FLUTTER SPEED 
    flutter_speed = None
    flutter_frequency = None
    for i, V in enumerate(V_vec):
        for rt in roots[i]:
            if rt.real > 0.01 and rt.imag > 0: # the 0.01 can be modified depending on the case
                flutter_speed = V 
                flutter_frequency = rt.imag*V
                break
        if flutter_speed is not None:
            break 

    if flutter_speed is not None:
        print(f"The recuced flutter speed is approximately: {flutter_speed} and reduced flutter frequency is {flutter_frequency}")
    else:
        print("No flutter speed found within the specified range.")
    return flutter_speed


# Dimensionless coefficients 
a = -1/5
e = -1/10
mu = 20
r = np.sqrt(6/25)
sigma = 2/5
x_theta = e-a

# Range of dimensionless velocity to study
V_vec = np.linspace(0.01, 3, 200)

# Matrix to store the complex-roots -> for a given reduced velocity, the 4 roots are stored in the corresponding row
roots = np.zeros((len(V_vec), 4), dtype=complex)


for i,V  in enumerate(V_vec):
    # Construct mass matrix
    M = np.zeros((2,2)) # Mass matrix
    M[0,0] = 1
    M[0,1] = x_theta
    M[1,0] = x_theta
    M[1,1] = r**2
    # Construct stifness matrix
    K = np.zeros((2,2)) # Stifness matrix
    K[0,0] = sigma**2/V**2
    K[0,1] = 2/mu
    K[1,1] = r**2/V**2 -2/mu*(1/2+a)

    # Construct state-spaces matrix (see Handout 1 for definition)
    B = np.zeros((4,4))
    B[0:2, 0:2] = M
    B[2:4, 2:4] = np.identity(2)
    A = np.zeros((4,4))
    A[0:2, 2:4] = -K
    A[2:4, 0:2] = np.identity(2)

    # Solve the generalized eigenvalue problem
    roots[i, :], _ = eig(A, B)


# Create subplots
fig, axs = plt.subplots(1, 2, figsize=(16, 6))  # 1 row, 2 columns

# Plot modal damping versus reduced velocity
axs[0].plot(V_vec, roots[:, 0].real * V_vec, label=r"$\Gamma_1/(\omega_{\theta})$", color="tab:blue", linewidth=2)
axs[0].plot(V_vec, roots[:, 1].real * V_vec, label=r"$\Gamma_2/(\omega_{\theta})$", color="tab:orange", linewidth=2)
axs[0].plot(V_vec, roots[:, 2].real * V_vec, label=r"$\Gamma_3/(\omega_{\theta})$", color="tab:purple", linewidth=2)
axs[0].plot(V_vec, roots[:, 3].real * V_vec, label=r"$\Gamma_4/(\omega_{\theta})$", color="tab:green", linewidth=2)
axs[0].set_xlabel(r"$V/(b\omega_{\theta})$ - Reduced velocity", fontsize=12)
axs[0].set_ylabel(r"$\Gamma/\omega_{\theta}$", fontsize=12)
axs[0].legend()
axs[0].set_title("Modal Damping vs Reduced Velocity", fontsize=14)

# Plot modal frequency versus reduced velocity
axs[1].plot(V_vec, roots[:, 0].imag * V_vec, label=r"$\Omega_1/(\omega_{\theta})$", color="tab:blue", linewidth=2)
axs[1].plot(V_vec, roots[:, 1].imag * V_vec, label=r"$\Omega_2/(\omega_{\theta})$", color="tab:orange", linewidth=2)
axs[1].plot(V_vec, roots[:, 2].imag * V_vec, label=r"$\Omega_3/(\omega_{\theta})$", color="tab:purple", linewidth=2)
axs[1].plot(V_vec, roots[:, 3].imag * V_vec, label=r"$\Omega_4/(\omega_{\theta})$", color="tab:green", linewidth=2)
axs[1].set_xlabel(r"$V/(b\omega_{\theta})$ - Reduced velocity", fontsize=12)
axs[1].set_ylabel(r"$\Omega/\omega_{\theta}$", fontsize=12)
axs[1].legend()
axs[1].set_title("Modal Frequency vs Reduced Velocity", fontsize=14)

# Adjust layout
plt.tight_layout()
plt.show()

identify_flutter_speed(V_vec, roots)












