import numpy as np
import matplotlib.pyplot as plt

def identify_flutter_speed(V_vec, k_sol, gamma_sol):
    # IDENTIFY REDUCED FLUTTER SPEED 
    flutter_speed = None
    flutter_frequency = None
    for i, V in enumerate(V_vec):
        for rt in range(0, 4):
            if gamma_sol[i,rt]*k_sol[i,rt]* V>0.01 and k_sol[i,rt]* V>0: # the 0.01 can be modified depending on the case
                flutter_speed = V 
                flutter_frequency = k_sol[i,rt]* V
                break
        if flutter_speed is not None:
            break 

    if flutter_speed is not None:
        print(f"The recuced flutter speed is approximately: {flutter_speed} and reduced flutter frequency is {flutter_frequency}")
    else:
        print("No flutter speed found within the specified range.")
    return flutter_speed



# Model validation for p-k-method
# Method is validated by comparings results with Hodges & Pierces (Fig. 5.20 & 5.21 [p.229-230])

def C(k): # Approximation of Theodorsen's function
    return (0.01365 + 0.2808j*k - k**2/2)/(0.01365 + 0.3455j*k -k**2)

def f11(k, V): # see handout 2 for further informations
    return sigma**2/V**2 - k**2/mu + 2j*k*C(k)/mu

def f12(k, V): # see handout 2 for further informations
    return (k*(1j+a*k) + (2+ k*1j*(1-2*a))*C(k))/mu

def f21(k,V): # see handout 2 for further informations
    return (a*k**2 -1j*k*(1+2*a)*C(k))/mu

def f22(k, V): # see handout 2 for further informations
    return ((8*mu*r**2)/V**2 + 4j*(1+2*a)*(2j-k*(1-2*a))*C(k)-k*(k-4j+8*a*(1j+a*k)))/(8*mu)

# Dimensionless coefficients 
a = -1/5
e = -1/10
mu = 20
r = np.sqrt(6/25)
sigma = 2/5
x_theta = e-a

# Range of dimensionless velocity to study
V_vec = np.linspace(0.01, 3, 200)

# Matrix to store the solutions -> for a given reduced velocity, the 4 solutions are stored in the corresponding row
k_sol = np.zeros((len(V_vec), 4))
gamma_sol = np.zeros((len(V_vec), 4))


# Numerical parameters
eps = 10**(-5) # Convergence criterion on k and gamma
itermax = 50000  # Number of maximum iterations
k0 = 0.1 # Initial guess for k

for i, V in enumerate(V_vec):  # Outer loop over V_vec
    print(f"\nV = {V}")
    for root in range(4):  # Loop over roots
        k = [k0]  # Initial guess for k
        gamma = [0]  # Initial guess for gamma
        converged = False  # Flag to check convergence
        for iter in range(itermax):  # Iterative solver loop
            coeff = [r**2 - x_theta**2, 
                     0,
                     f22(k[-1], V) + f11(k[-1], V) * r**2 - f21(k[-1], V) * x_theta - f12(k[-1], V) * x_theta,
                     0,
                     f11(k[-1], V) * f22(k[-1], V) - f12(k[-1], V) * f21(k[-1], V)]
            
            sol = np.roots(coeff)
            # Sort the roots based on their imaginary parts
            sorted_sol = sorted(sol, key=lambda x: x.imag)
            k.append(np.abs(sorted_sol[root].imag))
            gamma.append(sorted_sol[root].real/ k[-1])
            
            # Check for convergence
            if np.abs(gamma[-1] - gamma[-2]) < eps and np.abs(k[-1] - k[-2]) < eps:
                print(f"Solver has converged for root n° {root+1} considering V={V} in {iter+1} iterations!")
                k_sol[i, root] = sorted_sol[root].imag
                gamma_sol[i, root] = sorted_sol[root].real/sorted_sol[root].imag
                converged = True  # Mark as converged
                break  # Exit the `iter` loop
        
        if not converged:
            print(f"Maximum number of iterations reached, solver has not converged for root n° {root+1} considering V={V}.")

# Create subplots
fig, axs = plt.subplots(1, 2, figsize=(16, 6))  # 1 row, 2 columns

# Plot modal damping versus reduced velocity
axs[0].plot(V_vec, gamma_sol[:,0]*k_sol[:,0]* V_vec, label=r"$\Gamma_1/(\omega_{\theta})$", color="tab:blue", linewidth=2)
axs[0].plot(V_vec, gamma_sol[:,1]*k_sol[:,1]* V_vec, label=r"$\Gamma_2/(\omega_{\theta})$", color="tab:orange", linewidth=2)
axs[0].plot(V_vec, gamma_sol[:,2]*k_sol[:,2]* V_vec, label=r"$\Gamma_3/(\omega_{\theta})$", color="tab:purple", linewidth=2)
axs[0].plot(V_vec, gamma_sol[:,3]*k_sol[:,3]* V_vec, label=r"$\Gamma_4/(\omega_{\theta})$", color="tab:green", linewidth=2)
axs[0].set_xlabel(r"$V/(b\omega_{\theta})$ - Reduced velocity", fontsize=12)
axs[0].set_ylabel(r"$\Gamma/\omega_{\theta}$", fontsize=12)
axs[0].legend(fontsize=15)
axs[0].set_title("Modal Damping vs Reduced Velocity", fontsize=14)

# Plot modal frequency versus reduced velocity
axs[1].plot(V_vec, k_sol[:,0]* V_vec, label=r"$\Omega_1/(\omega_{\theta})$", color="tab:blue", linewidth=2)
axs[1].plot(V_vec, k_sol[:,1]* V_vec, label=r"$\Omega_2/(\omega_{\theta})$", color="tab:orange", linewidth=2)
axs[1].plot(V_vec, k_sol[:,2]* V_vec, label=r"$\Omega_3/(\omega_{\theta})$", color="tab:purple", linewidth=2)
axs[1].plot(V_vec, k_sol[:,3]* V_vec, label=r"$\Omega_4/(\omega_{\theta})$", color="tab:green", linewidth=2)
axs[1].set_xlabel(r"$V/(b\omega_{\theta})$ - Reduced velocity", fontsize=12)
axs[1].set_ylabel(r"$\Omega/\omega_{\theta}$", fontsize=12)
axs[1].legend(fontsize=15)
axs[1].set_title("Modal Frequency vs Reduced Velocity", fontsize=14)

# Adjust layout
plt.tight_layout()
plt.show()


identify_flutter_speed(V_vec, k_sol, gamma_sol)


