import numpy as np
import matplotlib.pyplot as plt

class pkmethod():
    def __init__(self, a, e, mu, r, sigma, V_vec, threshold_detection=0.01, itermax=1000, eps=10**(-5), k0 = 0.1):
        self.itermax = itermax
        self.eps = eps
        self.k0 = k0
        self.thresold = threshold_detection

        self.a = a
        self.e = e
        self.mu = mu
        self.r = r
        self.sigma = sigma
        self.V_vec = V_vec
        self.x_theta = e-a

        self.k_sol = None
        self.gamma_sol = None

        self.flutter_speed = None
        self.flutter_frequency = None
    
    def run(self):
        self.k_sol = np.zeros((len(self.V_vec), 4))
        self.gamma_sol = np.zeros((len(self.V_vec), 4))

        for i, V in enumerate(self.V_vec):  # Outer loop over V_vec
            print(f"\nV = {V}")
            for root in range(4):  # Loop over roots
                k = [self.k0]  # Initial guess for k
                gamma = [0]  # Initial guess for gamma
                converged = False  # Flag to check convergence
                for iter in range(self.itermax):  # Iterative solver loop
                    coeff = [self.r**2 - self.x_theta**2, 
                            0,
                            self.f22(k[-1], V) + self.f11(k[-1], V) * self.r**2 - self.f21(k[-1], V) * self.x_theta - self.f12(k[-1], V) * self.x_theta,
                            0,
                            self.f11(k[-1], V) * self.f22(k[-1], V) - self.f12(k[-1], V) * self.f21(k[-1], V)]                    
                    
                    sol = np.roots(coeff)
                    # Sort the roots based on their imaginary parts
                    sorted_sol = sorted(sol, key=lambda x: x.imag)
                    k.append(np.abs(sorted_sol[root].imag))
                    gamma.append(sorted_sol[root].real/ k[-1])
                    
                    # Check for convergence
                    if np.abs(gamma[-1] - gamma[-2]) < self.eps and np.abs(k[-1] - k[-2]) < self.eps:
                        print(f"Solver has converged for root n° {root+1} considering V={V} in {iter+1} iterations!")
                        self.k_sol[i, root] = sorted_sol[root].imag
                        self.gamma_sol[i, root] = sorted_sol[root].real/sorted_sol[root].imag
                        converged = True  # Mark as converged
                        break  # Exit the `iter` loop
                
                if not converged:
                    print(f"Maximum number of iterations reached, solver has not converged for root n° {root+1} considering V={V}.")
        
    def plot_results(self):
        if self.k_sol is not None:
            # Create subplots
            fig, axs = plt.subplots(1, 2, figsize=(16, 6))  # 1 row, 2 columns

            # Plot modal damping versus reduced velocity
            axs[0].plot(self.V_vec, self.gamma_sol[:,0]*self.k_sol[:,0]*self.V_vec, label=r"$\Gamma_1/(\omega_{\theta})$", color="tab:blue", linewidth=2)
            axs[0].plot(self.V_vec, self.gamma_sol[:,1]*self.k_sol[:,1]*self.V_vec, label=r"$\Gamma_2/(\omega_{\theta})$", color="tab:orange", linewidth=2)
            axs[0].plot(self.V_vec, self.gamma_sol[:,2]*self.k_sol[:,2]*self.V_vec, label=r"$\Gamma_3/(\omega_{\theta})$", color="tab:purple", linewidth=2)
            axs[0].plot(self.V_vec, self.gamma_sol[:,3]*self.k_sol[:,3]*self.V_vec, label=r"$\Gamma_4/(\omega_{\theta})$", color="tab:green", linewidth=2)
            axs[0].set_xlabel(r"$V/(b\omega_{\theta})$ - Reduced velocity", fontsize=12)
            axs[0].set_ylabel(r"$\Gamma/\omega_{\theta}$", fontsize=12)
            axs[0].legend(fontsize=15)
            axs[0].set_title("Modal Damping vs Reduced Velocity", fontsize=14)

            # Plot modal frequency versus reduced velocity
            axs[1].plot(self.V_vec, self.k_sol[:,0]*self.V_vec, label=r"$\Omega_1/(\omega_{\theta})$", color="tab:blue", linewidth=2)
            axs[1].plot(self.V_vec, self.k_sol[:,1]*self.V_vec, label=r"$\Omega_2/(\omega_{\theta})$", color="tab:orange", linewidth=2)
            axs[1].plot(self.V_vec, self.k_sol[:,2]*self.V_vec, label=r"$\Omega_3/(\omega_{\theta})$", color="tab:purple", linewidth=2)
            axs[1].plot(self.V_vec, self.k_sol[:,3]*self.V_vec, label=r"$\Omega_4/(\omega_{\theta})$", color="tab:green", linewidth=2)
            axs[1].set_xlabel(r"$V/(b\omega_{\theta})$ - Reduced velocity", fontsize=12)
            axs[1].set_ylabel(r"$\Omega/\omega_{\theta}$", fontsize=12)
            axs[1].legend(fontsize=15)
            axs[1].set_title("Modal Frequency vs Reduced Velocity", fontsize=14)

            # Adjust layout
            plt.tight_layout()
            plt.show()
        else:
            print("Error : you have to use first .run() method ! ")
        
    def find_flutter(self):
        if self.k_sol is None: 
            print("Error : you have to use first .run() method ! ")
            return -1
        else:
            self.flutter_speed = None
            self.flutter_frequency = None
            for i, V in enumerate(self.V_vec):
                for rt in range(0, 4):
                    if self.gamma_sol[i,rt]*self.k_sol[i,rt]*V>self.thresold and self.k_sol[i,rt]*V>0: # the 0.01 can be modified depending on the case
                        self.flutter_speed = V 
                        self.flutter_frequency = self.k_sol[i,rt]* V
                        break
                if self.flutter_speed is not None:
                    break 

            if self.flutter_speed is not None:
                print(f"The reduced flutter speed is approximately: {self.flutter_speed} and reduced flutter frequency is {self.flutter_frequency}")
            else:
                print("No flutter speed found within the specified range.")
            return self.flutter_speed, self.flutter_frequency
        
    def C(self, k): # Approximation of Theodorsen's function
        return (0.01365 + 0.2808j*k - k**2/2)/(0.01365 + 0.3455j*k -k**2)

    def f11(self, k, V): # see handout 2 for further informations
        return self.sigma**2/V**2 - k**2/self.mu + 2j*k*self.C(k)/self.mu

    def f12(self, k, V): # see handout 2 for further informations
        return (k*(1j+self.a*k) + (2+ k*1j*(1-2*self.a))*self.C(k))/self.mu

    def f21(self, k,V): # see handout 2 for further informations
        return (self.a*k**2 -1j*k*(1+2*self.a)*self.C(k))/self.mu

    def f22(self, k, V): # see handout 2 for further informations
        return ((8*self.mu*self.r**2)/V**2 + 4j*(1+2*self.a)*(2j-k*(1-2*self.a))*self.C(k)-k*(k-4j+8*self.a*(1j+self.a*k)))/(8*self.mu)


