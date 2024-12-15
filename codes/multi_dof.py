import matplotlib.pyplot as plt
import numpy as np
import sympy as sp
from scipy.special import hankel1

# Constants definition
RHO_0 = 1.225 # Density at sea level [kg/m³]
R = 287.058 # Temperature of at sea level [K]
T0 = 288.15 # Specific gas constant [J/(kg.K)]
G = 9.80665  # Gravitational constant [m/s²]

# Theodorsen function C(k)
def theodorsen_function(k): # Approximation of Theodorsen's function
    return (0.01365 + 0.2808j*k - k**2/2)/(0.01365 + 0.3455j*k -k**2)

def flight_envelope(h):
    # This function returns the minimum and maximum true speed for a given altitude h [km]
    if h < 0.8:
        vmin = 0
        vmax = 252 + 62.5 * h
    elif h > 1.2:
        vmin = 62.5 * (h - 1.2)
        vmax = 302
    else:
        vmin = 0
        vmax = 302
    return vmin, vmax

def rho_fun(h):
   # This function returns the density of air [kg/m³] depending on a given altitude h [km]
   # International Standard Atmosphere is used
   T = T0 - 0.0065*h*10**(3)
   rho = RHO_0*(T/T0)**(G/(0.0065*R)-1)
   return rho


def plot_flight_envelope():
    # Define h vector
    h_vec = np.linspace(0, 3, 100)
    V = np.zeros((len(h_vec), 2)) # for each row (h constant) -> store vmin, vmax

    # Calculate vmin and vmax
    for i, h in enumerate(h_vec):
        V[i, 0], V[i, 1] = flight_envelope(h)

    # First subplot: True airspeed flight envelope
    h_closed = np.concatenate([h_vec, h_vec[::-1]])  # Combine h and its reverse
    V_closed = np.concatenate([V[:, 1], V[::-1, 0]])  # Combine vmax and reversed vmin
    plt.plot(V[:, 0], h_vec, color="blue", label=r"Min true airspeed $V_{min}^{true}$")
    plt.plot(V[:, 1], h_vec, color="red", label=r"Max true airspeed $V_{max}^{true}$")
    plt.fill(V_closed, h_closed, color="lightblue", alpha=0.6, label="Flight Envelope")
    
    # Adding the yellow zone from Vmax to 1.15*Vmax
    V_1_15_max = 1.15 * V[:, 1]
    plt.fill_betweenx(h_vec, V[:, 1], V_1_15_max, color="yellow", alpha=0.5, label=r"Margin zone")
    plt.plot(V_1_15_max, h_vec, color="black")

    plt.hlines(0, flight_envelope(0)[0], 1.15*flight_envelope(0)[1], colors="black", linestyles="--", linewidth=1.5)
    plt.hlines(3, flight_envelope(3)[0], 1.15*flight_envelope(3)[1], colors="black", linestyles="--", linewidth=1.5)


    
    plt.xlabel("True airspeed [km/h]", fontsize=12)
    plt.ylabel("Altitude (h) [km]", fontsize=12)
    plt.title("True Airspeed Flight Envelope", fontsize=14)
    plt.xlim(-10, 400)
    plt.legend(fontsize=12)
    plt.grid(True)

    # Show plot
    plt.show()
    

# --------------------------------------- Parameters ------------------------------------- #

# Regulation parameters
alpha = 1.15

# -------------------- Wing parameters ---------------------- #

# Team parameters (team 5 - see Excel sheet)
S = 17.0 # wing area (whole wing) [m²]
AR = 7.55 # Aspect ratio of the wing
l = np.sqrt(AR*S)/2 # Span (length) of the right/left wing

me = 26.91*S # Empty mass of the whole wing
m_fuel = 80 # mass of fuel carried in each wing (i.e right and left)

Ip_max = 7 # Max inertia (max weight) of the right/left at point P [kg.m²]
e_max = -0.3 # Position of the center of mass at max weight situation
m_max = me/2 + m_fuel # max weight of the right/left wing

Ip_min = 4 # Min inertia (min weight) of the right/left at point P [kg.m²]
e_min = -0.1 # Position of the center of mass at min weight situation
m_min = me/2 # min weight of the right/left wing

a = -0.2 # Position of the elastic axis in all situations

c = S/(2*l) # Chord length
b = c/2 # See course

r_min = np.sqrt(Ip_min/(m_min*b**2))
r_max = np.sqrt(Ip_max/(m_max*b**2))

EI = 2*10**5 # Bending rigidity of the wing [N.m²]
GJ = 10**5 # Torsional rigidity of the wing [N.m²]

# See appendix for these calculations
omega_h_min = (1.8751)**2*np.sqrt(EI/(m_min*l**3))
omega_h_max = (1.8751)**2*np.sqrt(EI/(m_max*l**3))

omega_theta_min = np.pi/2*np.sqrt(GJ/(Ip_min*l))
omega_theta_max = np.pi/2*np.sqrt(GJ/(Ip_max*l))

# Print results
print(f"{'Parameter':<30} {'Symbol':<15} {'Value':<15} {'Unit'}")
print("-" * 65)
print(f"{'Wing Area (whole wing)':<30} {'S':<15} {S:<15.2f} {'m²'}")
print(f"{'Aspect Ratio':<30} {'AR':<15} {AR:<15.2f} {'-'}")
print(f"{'Span (right/left wing)':<30} {'l':<15} {l:<15.2f} {'m'}")
print(f"{'Chord Length':<30} {'c':<15} {c:<15.2f} {'m'}")
print(f"{'Half-Chord':<30} {'b':<15} {b:<15.2f} {'m'}")
print(f"{'Empty Mass (whole wing)':<30} {'m_e':<15} {me:<15.2f} {'kg'}")
print(f"{'Mass of Fuel (each wing)':<30} {'m_fuel':<15} {m_fuel:<15.2f} {'kg'}")
print(f"{'Max Mass (each wing)':<30} {'m_max':<15} {m_max:<15.2f} {'kg'}")
print(f"{'Min Mass (each wing)':<30} {'m_min':<15} {m_min:<15.2f} {'kg'}")
print(f"{'Elastic Axis Location':<30} {'a':<15} {a:<15.2f} {'-'}")
print(f"{'Center of Mass (max)':<30} {'e_max':<15} {e_max:<15.2f} {'-'}")
print(f"{'Center of Mass (min)':<30} {'e_min':<15} {e_min:<15.2f} {'-'}")
print(f"{'Bending Rigidity':<30} {'EI':<15} {EI:<15.2e} {'N·m²'}")
print(f"{'Torsional Rigidity':<30} {'GJ':<15} {GJ:<15.2e} {'N·m²'}")
print(f"{'Inertia (max)':<30} {'Ip_max':<15} {Ip_max:<15.2f} {'kg·m²'}")
print(f"{'Inertia (min)':<30} {'Ip_min':<15} {Ip_min:<15.2f} {'kg·m²'}")
print(f"{'Plunging natural frequency (min)':<30} {'omega_h_min':<15} {omega_h_min:<15.2f} {'rad/s'}")
print(f"{'Plunging natural frequency (max)':<30} {'omega_h_max':<15} {omega_h_max:<15.2f} {'rad/s'}")
print(f"{'Pitching natural frequency (min)':<30} {'omega_theta_min':<15} {omega_theta_min:<15.2f} {'rad/s'}")
print(f"{'Pitching natural frequency (max)':<30} {'omega_theta_max':<15} {omega_theta_max:<15.2f} {'rad/s'}")

def get_alphaL4(n_w):
    alphaL4 = sp.zeros(n_w)
    for i in range(n_w):
        if i == 0:
            alphaL4[i, i] = 1.875**4
        elif i == 1:
            alphaL4[i, i] = 4.694**4
        elif i == 2 : 
            alphaL4[i, i] = 7.855**4
        else : 
            alphaL4[i, i] = ((2*(i+1)-1)*np.pi/2)**4
    return alphaL4

def get_gammaL2(n_theta):
    gammaL2 = sp.zeros(n_theta)
    for i in range(n_theta):
        gammaL2[i, i] =  ((2*(i+1)-1)*np.pi/2)**2
    return gammaL2

def get_matrix_A(n_w, n_theta, l):
    A = sp.zeros(n_w, n_theta)
    y = sp.symbols("y")
    for i in range(n_w):
        for j in range(n_theta):
            if i==0:
                alpha_i = 1.875/l
            elif i == 1:
                alpha_i = 4.694/l
            elif i == 2 :
                alpha_i = 7.855/l
            else:
                alpha_i = (2*(i+1)-1)*np.pi/(2*l)
            beta_i = (-sp.sin(alpha_i*l) + sp.sinh(alpha_i*l))/(sp.cos(alpha_i*l)+sp.cosh(alpha_i*l))
            f = ((sp.cosh(alpha_i*y)-sp.cos(alpha_i*y))-beta_i*(sp.sinh(alpha_i*y)-sp.sin(alpha_i*y)))*sp.sqrt(2)*sp.sin(y*(2*(j+1)-1)*np.pi/(2*l))
            A[i,j] = 1/l*sp.integrate(f, (y, 0, l))
    return A


def get_det(n_w, n_theta, m, l, b, r, EI, GJ, x_theta, a, rho):
    p, U_v, k_v = sp.symbols('p U_v k_v')

    # Initialize a matrix of size (n_w + n_theta) x (n_w + n_theta) with zeros
    matrix_1 = sp.zeros(n_w+n_theta)
    matrix_2 = sp.zeros(n_w+n_theta)
    matrix_3 = sp.zeros(n_w+n_theta)
    matrix_4 = sp.zeros(n_w+n_theta)

    # Create an identity matrix of size n_w and one of size n_theta
    Inw = sp.eye(n_w)
    Intheta = sp.eye(n_theta)

    # Construct matrices alphaL⁴ and gammaL² and A 
    alphaL4 = get_alphaL4(n_w)
    gammaL2 = get_gammaL2(n_theta)
    A = get_matrix_A(n_w, n_theta, l)

    # Construct matrix 1
    matrix_1[0:n_w, 0:n_w] = p**2*(m*l*U_v**2/b**2)*Inw + EI/l**3 * alphaL4
    matrix_1[n_w:, n_w:] = r**2*p**2*(m*l*U_v**2/b)*Intheta + GJ/l*gammaL2
    matrix_1[0:n_w, n_w:] = -2*p**2*x_theta*(m*l*U_v**2/b)*A
    matrix_1[n_w:, 0:n_w] = -2*p**2*x_theta*(m*l*U_v**2/b)*sp.transpose(A)

    # Construct matrix 2
    matrix_2 [0:n_w, 0:n_w] = Inw
    matrix_2[n_w:, n_w:] = (b**2*a**2 + 1/8)*Intheta
    matrix_2[0:n_w, n_w:] = b*a*A
    matrix_2[n_w:, 0:n_w] = b*a*sp.transpose(A)
    matrix_2 = matrix_2*k_v**2*U_v**2/b**2*np.pi*rho*l

    # Construct matrix 3
    matrix_3[0:n_w, 0:n_w] = 2*theodorsen_function(k_v)*Inw
    matrix_3[n_w:, n_w:] = b**2*(1/2-a)*(1-2*(1/2+a)*theodorsen_function(k_v))*Intheta
    matrix_3[0:n_w, n_w:] = -b*(1+2*(1/2-a)*theodorsen_function(k_v))*A
    matrix_3[n_w:, 0:n_w] = 2*b*(1/2+a)*theodorsen_function(k_v)*sp.transpose(A)
    matrix_3 = matrix_3*1j*np.pi*rho*U_v**2*l*k_v

    # Construct matrix 4
    matrix_4[n_w:, n_w:] = -b*(1+2*a)*theodorsen_function(k_v)*Intheta
    matrix_4[0:n_w, n_w:] = -2*theodorsen_function(k_v)*A
    matrix_4 = matrix_4*np.pi*rho*b*U_v**2*l

    mat = matrix_1-matrix_2+matrix_3+matrix_4
    mat = mat.expand()
    det = mat.det(method="laplace") 
    return det


# p-k method for multi degree of freedom model

# Numerical parameters
N_w = 3 # number of modes conserved for heaving motion 
N_theta = 3 # number of modes conserved for pitching motion
itermax = 200  # Maximum number of iterations
eps = 10**(-3) # Convergence criterion
k0 = 1 # Initial guess for k

# Other parameters
n_roots = (N_w+N_theta)*2
V_vec = np.linspace(0.1, 2, 100) # dimensional speed in m/s
U_vec = V_vec*b*omega_theta_min

k_sol = np.zeros((len(U_vec), n_roots))
gamma_sol = np.zeros((len(U_vec), n_roots))

# Analysis for h=1.5
rho = rho_fun(1.5)
# Determinant is computed symbolically
det = get_det(N_w, N_theta, m_min, l, b, r_min, EI, GJ, e_min-a, a, rho)

for i, U in enumerate(U_vec):  # Outer loop over V_vec
    print(f"\nU = {U}")
    for root in range(n_roots):  # Loop over roots
        k = [k0]  # Initial guess for k
        gamma = [0]  # Initial guess for gamma
        converged = False  # Flag to check convergence
        for iter in range(itermax):  # Iterative solver loop

            p, U_v, k_v = sp.symbols('p U_v k_v')
            values = {U_v: U, k_v: k[-1]} 
            det_new = det.subs(values)
            coeff = sp.Poly(det_new, p).all_coeffs()

            coeff_np = np.fromiter(coeff, dtype=complex)
            sol = np.roots(coeff_np)
            
        
            # Sort the roots based on their imaginary parts
            sorted_sol = sorted(sol, key=lambda x: (np.abs(x.imag)!=0, np.abs(x.real)!=0, x.imag, np.sign(x.imag)))
            k.append(np.abs(sorted_sol[root].imag))
            gamma.append(sorted_sol[root].real/ k[-1])
                    
            # Check for convergence
            if np.abs(gamma[-1] - gamma[-2]) < eps and np.abs(k[-1] - k[-2]) < eps:
                print(f"Solver has converged for root n° {root+1} considering U={U} in {iter+1} iterations!")
                k_sol[i, root] = sorted_sol[root].imag
                gamma_sol[i, root] = sorted_sol[root].real/sorted_sol[root].imag
                converged = True  # Mark as converged
                break  # Exit the `iter` loop
                
        if not converged:
            print(f"Maximum number of iterations reached, solver has not converged for root n° {root+1} considering U={U}.")


# Create subplots
fig, axs = plt.subplots(1, 2, figsize=(16, 6))  # 1 row, 2 columns

# Colors for the plots
colors = plt.cm.tab10.colors  # Tab10 color map (supports up to 10 colors)

# Plot modal damping versus reduced velocity
for i in range(n_roots):
    axs[0].plot(
        U_vec / (b * omega_theta_min),
        gamma_sol[:, i] * k_sol[:, i] * (U_vec / (b * omega_theta_min)),
        label=rf"$\Gamma_{i+1}/(\omega_\theta)$",
        color=colors[i % len(colors)],
        linewidth=2
    )

axs[0].set_xlabel(r"$V/(b\omega_{\theta})$ - Reduced velocity", fontsize=12)
axs[0].set_ylabel(r"$\Gamma/\omega_{\theta}$", fontsize=12)
axs[0].legend(fontsize=15)
axs[0].set_title("Modal Damping vs Reduced Velocity", fontsize=14)

# Plot modal frequency versus reduced velocity
for i in range(n_roots):
    axs[1].plot(
        U_vec / (b * omega_theta_min),
        k_sol[:, i] * (U_vec / (b * omega_theta_min)),
        label=rf"$\Omega_{i+1}/(\omega_\theta)$",
        color=colors[i % len(colors)],
        linewidth=2
    )

axs[1].set_xlabel(r"$V/(b\omega_{\theta})$ - Reduced velocity", fontsize=12)
axs[1].set_ylabel(r"$\Omega/\omega_{\theta}$", fontsize=12)
axs[1].legend(fontsize=15)
axs[1].set_title("Modal Frequency vs Reduced Velocity", fontsize=14)

# Adjust layout
plt.tight_layout()
plt.show()




