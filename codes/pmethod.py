import matplotlib.pyplot as plt
import numpy as np
from scipy.linalg import eig

class pmethod():
    def __init__(self, a, e, mu, r, sigma, V_vec, threshold_detection=0.01):
        self.a = a
        self.e = e
        self.mu = mu
        self.r = r
        self.sigma = sigma
        self.x_theta = e-a
        self.V_vec = V_vec
        self.roots = None
        self.flutter_speed = None
        self.flutter_frequency = None
        self.threshold = threshold_detection

    def run(self):
        # Matrix to store the complex-roots -> for a given reduced velocity, the 4 roots are stored in the corresponding row
        self.roots = np.zeros((len(self.V_vec), 4), dtype=complex)
        for i,V  in enumerate(self.V_vec):
            # Construct mass matrix
            M = np.zeros((2,2)) # Mass matrix
            M[0,0] = 1
            M[0,1] = self.x_theta
            M[1,0] = self.x_theta
            M[1,1] = self.r**2
            # Construct stifness matrix
            K = np.zeros((2,2)) # Stifness matrix
            K[0,0] = self.sigma**2/V**2
            K[0,1] = 2/self.mu
            K[1,1] = self.r**2/V**2 -2/self.mu*(1/2+self.a)

            # Construct state-spaces matrix (see Handout 1 for definition)
            B = np.zeros((4,4))
            B[0:2, 0:2] = M
            B[2:4, 2:4] = np.identity(2)
            A = np.zeros((4,4))
            A[0:2, 2:4] = -K
            A[2:4, 0:2] = np.identity(2)

            # Solve the generalized eigenvalue problem
            self.roots[i, :], _ = eig(A, B)
    
    def plot_results(self):
        if self.roots is None: 
            print("Error : you have to use first .run() method ! ")
        else:
            # Create subplots
            fig, axs = plt.subplots(1, 2, figsize=(16, 6))  # 1 row, 2 columns

            # Plot modal damping versus reduced velocity
            axs[0].plot(self.V_vec, self.roots[:, 0].real * self.V_vec, label=r"$\Gamma_1/(\omega_{\theta})$", color="tab:blue", linewidth=2)
            axs[0].plot(self.V_vec, self.roots[:, 1].real * self.V_vec, label=r"$\Gamma_2/(\omega_{\theta})$", color="tab:orange", linewidth=2)
            axs[0].plot(self.V_vec, self.roots[:, 2].real * self.V_vec, label=r"$\Gamma_3/(\omega_{\theta})$", color="tab:purple", linewidth=2)
            axs[0].plot(self.V_vec, self.roots[:, 3].real * self.V_vec, label=r"$\Gamma_4/(\omega_{\theta})$", color="tab:green", linewidth=2)
            axs[0].set_xlabel(r"$V/(b\omega_{\theta})$ - Reduced velocity", fontsize=12)
            axs[0].set_ylabel(r"$\Gamma/\omega_{\theta}$", fontsize=12)
            axs[0].legend()
            axs[0].set_title("Modal Damping vs Reduced Velocity", fontsize=14)

            # Plot modal frequency versus reduced velocity
            axs[1].plot(self.V_vec, self.roots[:, 0].imag * self.V_vec, label=r"$\Omega_1/(\omega_{\theta})$", color="tab:blue", linewidth=2)
            axs[1].plot(self.V_vec, self.roots[:, 1].imag * self.V_vec, label=r"$\Omega_2/(\omega_{\theta})$", color="tab:orange", linewidth=2)
            axs[1].plot(self.V_vec, self.roots[:, 2].imag * self.V_vec, label=r"$\Omega_3/(\omega_{\theta})$", color="tab:purple", linewidth=2)
            axs[1].plot(self.V_vec, self.roots[:, 3].imag * self.V_vec, label=r"$\Omega_4/(\omega_{\theta})$", color="tab:green", linewidth=2)
            axs[1].set_xlabel(r"$V/(b\omega_{\theta})$ - Reduced velocity", fontsize=12)
            axs[1].set_ylabel(r"$\Omega/\omega_{\theta}$", fontsize=12)
            axs[1].legend()
            axs[1].set_title("Modal Frequency vs Reduced Velocity", fontsize=14)

            # Adjust layout
            plt.tight_layout()
            plt.show()
    
    def find_flutter(self):
        if self.roots is None: 
            print("Error : you have to use first .run() method ! ")
            return -1
        else:
            self.flutter_speed = None
            self.flutter_frequency = None
            # IDENTIFY REDUCED FLUTTER SPEED 
            for i, V in enumerate(self.V_vec):
                for rt in self.roots[i]:
                    if rt.real > 0.01 and rt.imag > 0: # the 0.01 can be modified depending on the case
                        self.flutter_speed = V 
                        self.flutter_frequency = rt.imag*V
                        break
                if self.flutter_speed is not None:
                    break 

            if self.flutter_speed is not None:
                print(f"The recuced flutter speed is approximately: {self.flutter_speed} and reduced flutter frequency is {self.flutter_frequency}")
            else:
                print("No flutter speed found within the specified range.")
            return self.flutter_speed, self.flutter_frequency
    











