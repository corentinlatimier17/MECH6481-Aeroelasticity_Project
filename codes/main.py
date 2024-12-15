import matplotlib.pyplot as plt
import numpy as np
from pmethod import *
from pkmethod import *

# ------------------------------------------------------------ Constants definition -----------------------------------------------------------------# 
RHO_0 = 1.225 # Density at sea level [kg/m³]
R = 287.058 # Temperature of at sea level [K]
T0 = 288.15 # Specific gas constant [J/(kg.K)]
G = 9.80665  # Gravitational constant [m/s²]

# ------------------------------------------------------------ Functions definition -----------------------------------------------------------------# 

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


# ------------------------------------------------------------ Main code  -----------------------------------------------------------------# 

# 1. Visualisation of the flight envelope
plot_flight_envelope() # plot the flight envelope of the aircraft

# 2. Definition of the variables

# a) Regulation parameters
alpha = 1.15

# b) Wing parameters
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

EI = 2*10**5 # Bending rigidity of the wing [N.m²]
GJ = 10**5 # Torsional rigidity of the wing [N.m²]

# See appendix for these calculations
omega_h_min = (1.8751)**2*np.sqrt(EI/(m_min*l**3))
omega_h_max = (1.8751)**2*np.sqrt(EI/(m_max*l**3))

omega_theta_min = np.pi/2*np.sqrt(GJ/(Ip_min*l))
omega_theta_max = np.pi/2*np.sqrt(GJ/(Ip_max*l))

# c) Print parameters 
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

# 3. Divergence analysis of the initial design

q_d = (np.pi/(2*l))**2*(GJ)/(c*4.5*b*(a+0.5)) # Divergence dynamic pressure (flexible beam -clamped-free modelisation)

def plot_flight_envelope_with_divergence_speed(q_d):
    # Define h vector
    h_vec = np.linspace(0, 3, 100)
    V = np.zeros((len(h_vec), 2)) # for each row (h constant) -> store vmin, vmax

    V_d = []

    # Calculate vmin and vmax
    for i, h in enumerate(h_vec):
        V[i, 0], V[i, 1] = flight_envelope(h)
        V_d.append(np.sqrt(2*q_d/rho_fun(h))*3.6) # divergence dynamic pressure in km/h


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

    plt.plot(V_d, h_vec, color="tab:purple", linewidth=3, linestyle="--", label=r"Divergence speed - $V_d$ ")

    plt.xlabel("True airspeed [km/h]", fontsize=12)
    plt.ylabel("Altitude (h) [km]", fontsize=12)
    plt.title("True Airspeed Flight Envelope with divergence speed", fontsize=14)
    plt.xlim(-10, 400)
    plt.legend(fontsize=12)
    plt.grid(True)

    # Show plot
    plt.show()

# Plot results 
plot_flight_envelope_with_divergence_speed(q_d)

# 4. Flutter analysis of the initial design

# For each altitude h, we compute the flutter velocity using both p and p-k method
h_vec = np.linspace(0,3,50)

V_min_max = np.zeros((len(h_vec), 2))

V_flutter_pk_min = []
V_flutter_pk_max = []


for i, h in enumerate(h_vec) : 
    print(h)
    V_min_max[i, 0], V_min_max[i, 1] = flight_envelope(h)

    rho = rho_fun(h) # get the density at this level

    # Case 1 - min weight
    mu_min = m_min/(np.pi*rho*b**2*l)
    r_min = np.sqrt(Ip_min/(m_min*b**2))
    sigma_min = omega_h_min/omega_theta_min

    # Case 2 - max weight
    mu_max = m_max/(np.pi*rho*b**2*l)
    r_max = np.sqrt(Ip_max/(m_max*b**2))
    sigma_max = omega_h_max/omega_theta_max

    # Definition of the reduced velocity to explore
    V_vec_min = np.linspace(0.01,3*V_min_max[i,1]/(3.6*b*omega_theta_min),400)
    V_vec_max = np.linspace(0.01,3*V_min_max[i,1]/(3.6*b*omega_theta_max),400)

    # Results with p-k method

    pkmethod_min = pkmethod(a=a, e=e_min, mu=mu_min, r=r_min, sigma=sigma_min, V_vec=V_vec_min)
    pkmethod_max = pkmethod(a=a, e=e_max, mu=mu_max, r=r_max, sigma=sigma_max, V_vec=V_vec_max)

    pkmethod_min.run()
    pkmethod_max.run()

    V_flutter_min_pk, _ = pkmethod_min.find_flutter()
    V_flutter_max_pk, _ = pkmethod_max.find_flutter()

    V_flutter_pk_min.append(V_flutter_min_pk)
    V_flutter_pk_max.append(V_flutter_max_pk)


V_flutter_pk_max = np.array(V_flutter_pk_max)*3.6*b*omega_theta_max
V_flutter_pk_min = np.array(V_flutter_pk_min)*3.6*b*omega_theta_min



def plot_flight_envelope_with_flutter(h_vec, V, V_flutter_pk_min=None, V_flutter_pk_max=None):
    # Combine h and its reverse
    h_closed = np.concatenate([h_vec, h_vec[::-1]])
    V_closed = np.concatenate([V[:, 1], V[::-1, 0]])

    # Plot the flight envelope
    plt.plot(V[:, 0], h_vec, color="blue", label=r"Min true airspeed $V_{min}^{true}$")
    plt.plot(V[:, 1], h_vec, color="red", label=r"Max true airspeed $V_{max}^{true}$")
    plt.fill(V_closed, h_closed, color="lightblue", alpha=0.6, label="Flight Envelope")

    # Adding the yellow zone from Vmax to 1.15*Vmax
    V_1_15_max = 1.15 * V[:, 1]
    plt.fill_betweenx(h_vec, V[:, 1], V_1_15_max, color="yellow", alpha=0.5, label=r"Margin zone")
    plt.plot(V_1_15_max, h_vec, color="black")


    # Add flutter velocities if provided
    if V_flutter_pk_min is not None:
        plt.plot(V_flutter_pk_min, h_vec, color="green", linestyle="--", label=r"$V_{flutter}^{pk}$ - Min weight scenario", linewidth=2)
    if V_flutter_pk_max is not None:
        plt.plot(V_flutter_pk_max, h_vec, color="purple", linestyle="--", label=r"$V_{flutter}^{pk}$ - Max weight scenario", linewidth=2)

    # Add horizontal boundaries for the flight envelope
    plt.hlines(0, flight_envelope(0)[0], 1.15*flight_envelope(0)[1], colors="black", linestyles="--", linewidth=1.5)
    plt.hlines(3, flight_envelope(3)[0], 1.15*flight_envelope(3)[1], colors="black", linestyles="--", linewidth=1.5)

    # Plot settings
    plt.xlabel("True airspeed [km/h]", fontsize=15)
    plt.ylabel("Altitude (h) [km]", fontsize=15)
    plt.xlim(-10, 450)
    plt.legend(fontsize=12)
    plt.grid(True)

    # Show plot
    plt.show()

# Plot results 
plot_flight_envelope_with_flutter(h_vec, V_min_max, V_flutter_pk_min, V_flutter_pk_max)

# 5. Parametric analysis at h=1.5 km of mu, x_theta and r²

h = 1.5
rho = rho_fun(h)
V_min, V_max = flight_envelope(h)

# Case 1 - min weight
mu_min = m_min/(np.pi*rho*b**2*l)
r_min = np.sqrt(Ip_min/(m_min*b**2))
sigma_min = omega_h_min/omega_theta_min

# Mu analysis
mu_vec = np.linspace(0.1, 30, 200)
Vf = [] # list to store flutter speeds
V_vec = np.linspace(0.01, 2, 500)
for mu in mu_vec:
    pkmethod_ = pkmethod(a=a, e=e_min, mu=mu, r=r_min, sigma=sigma_min, V_vec=V_vec, itermax=1000, eps=10**(-4))
    pkmethod_.run()
    V_flutter_pk, _ = pkmethod_.find_flutter()
    Vf.append(V_flutter_pk)

plt.scatter(mu_vec, Vf, s=10, marker='o')
mu_min = m_min/(np.pi*rho*b**2*l)
pkmethod_ = pkmethod(a=a, e=e_min, mu=mu_min, r=r_min, sigma=sigma_min, V_vec=V_vec, itermax=1000, eps=10**(-4))
pkmethod_.run()
V_flutter_pk, _ = pkmethod_.find_flutter()
plt.scatter(mu_min, V_flutter_pk, s=30, marker='^', color="tab:orange", label="Empty wing scenario")
plt.xlabel(r"$\mu $ - Mass ratio", fontsize=15)
plt.ylabel(r"$\frac{U}{b\omega_{\theta}}$ - Reduced flutter velocity", fontsize=15)
plt.grid(True)
plt.legend(fontsize=15)
plt.show()



# r² analysis
r_vec = np.linspace(0.15, 0.25, 100)
Vf = [] # list to store flutter speeds
V_vec = np.linspace(0.01, 2, 500)
for r in r_vec:
    pkmethod_ = pkmethod(a=a, e=e_min, mu=mu_min, r=r, sigma=sigma_min, V_vec=V_vec, itermax=1000, eps=10**(-4))
    pkmethod_.run()
    V_flutter_pk, _ = pkmethod_.find_flutter()
    Vf.append(V_flutter_pk)

plt.scatter(r_vec**2, Vf, s=10, marker='o')
pkmethod_ = pkmethod(a=a, e=e_min, mu=mu_min, r=r_min, sigma=sigma_min, V_vec=V_vec, itermax=1000, eps=10**(-4))
pkmethod_.run()
V_flutter_pk, _ = pkmethod_.find_flutter()
plt.scatter(r_min**2, V_flutter_pk, s=30, marker='^', color="tab:orange", label="Empty wing scenario")
plt.xlabel(r"$r^2 $ - Dimensionless radius of gyration", fontsize=15)
plt.ylabel(r"$\frac{U}{b\omega_{\theta}}$ - Reduced flutter velocity", fontsize=15)
plt.grid(True)
plt.legend(fontsize=15)
plt.show()


# X-theta analysis

e_vec = np.linspace(-0.30, -0.05, 150)
Vf = [] # list to store flutter speeds
V_vec = np.linspace(0.01, 2, 500)
for e in e_vec:
    pkmethod_ = pkmethod(a=a, e=e, mu=mu_min, r=r_min, sigma=sigma_min, V_vec=V_vec, itermax=1000, eps=10**(-4))
    pkmethod_.run()
    V_flutter_pk, _ = pkmethod_.find_flutter()
    Vf.append(V_flutter_pk)

plt.scatter(e_vec-a, Vf, s=10, marker='o')
pkmethod_ = pkmethod(a=a, e=e_min, mu=mu_min, r=r_min, sigma=sigma_min, V_vec=V_vec, itermax=1000, eps=10**(-4))
pkmethod_.run()
V_flutter_pk, _ = pkmethod_.find_flutter()
plt.scatter(e_min-a, V_flutter_pk, s=30, marker='^', color="tab:orange", label="Empty wing scenario")
plt.xlabel(r"$x_{\theta}", fontsize=15)
plt.ylabel(r"$\frac{U}{b\omega_{\theta}}$ - Reduced flutter velocity", fontsize=15)
plt.grid(True)
plt.legend(fontsize=15)
plt.show()

# 6. Analysis of design modifications

# a) Wing parameters (modified design) 

# Team parameters (team 5 - see Excel sheet)
S = 17.0 # wing area (whole wing) [m²]
AR = 7.55 # Aspect ratio of the wing
l = np.sqrt(AR*S)/2 # Span (length) of the right/left wing

me = 26.91*S*0.8 # Empty mass of the whole wing
m_fuel = 80 # mass of fuel carried in each wing (i.e right and left)

Ip_max = 7*0.8 # Max inertia (max weight) of the right/left at point P [kg.m²]
e_max = -0.44 # Position of the center of mass at max weight situation
m_max = me/2 + m_fuel # max weight of the right/left wing

Ip_min = 4*0.8 # Min inertia (min weight) of the right/left at point P [kg.m²]
e_min = -0.24 # Position of the center of mass at min weight situation
m_min = me/2 # min weight of the right/left wing

a = -0.2 # Position of the elastic axis in all situations

c = S/(2*l) # Chord length
b = c/2 # See course

EI = 2*10**5*1.15 # Bending rigidity of the wing [N.m²]
GJ = 10**5*1.15 # Torsional rigidity of the wing [N.m²]

# See appendix for these calculations
omega_h_min = (1.8751)**2*np.sqrt(EI/(m_min*l**3))
omega_h_max = (1.8751)**2*np.sqrt(EI/(m_max*l**3))

omega_theta_min = np.pi/2*np.sqrt(GJ/(Ip_min*l))
omega_theta_max = np.pi/2*np.sqrt(GJ/(Ip_max*l))

# Print parameters
print(f"{'Parameter (modified design)':<30} {'Symbol':<15} {'Value':<15} {'Unit'}")
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

# b) Divergence analysis (modified design) 

q_d = (np.pi/(2*l))**2*(GJ)/(c*4.5*b*(a+0.5))
plot_flight_envelope_with_divergence_speed(q_d)

# c) Flutter analysis (modified design) [only for min weight scenario here]

# For each altitude h, we compute the flutter velocity using both p and p-k method
h_vec = np.linspace(0,3,50)

V_min_max = np.zeros((len(h_vec), 2))

V_flutter_pk_min = []

for i, h in enumerate(h_vec) : 
    V_min_max[i, 0], V_min_max[i, 1] = flight_envelope(h)

    rho = rho_fun(h) # get the density at this level

    # Case 1 - min weight
    mu_min = m_min/(np.pi*rho*b**2*l)
    r_min = np.sqrt(Ip_min/(m_min*b**2))
    sigma_min = omega_h_min/omega_theta_min

    # Definition of the reduced velocity to explore
    V_vec_min = np.linspace(0.01,3*V_min_max[i,1]/(3.6*b*omega_theta_min),400)

    pkmethod_min = pkmethod(a=a, e=e_min, mu=mu_min, r=r_min, sigma=sigma_min, V_vec=V_vec_min)
    pkmethod_min.run()
    V_flutter_min_pk, _ = pkmethod_min.find_flutter()

    V_flutter_pk_min.append(V_flutter_min_pk)

V_flutter_pk_min = np.array(V_flutter_pk_min)*3.6*b*omega_theta_min


def plot_flight_envelope_with_flutter(h_vec, V, V_flutter_pk_min=None, V_flutter_pk_max=None):
    # Combine h and its reverse
    h_closed = np.concatenate([h_vec, h_vec[::-1]])
    V_closed = np.concatenate([V[:, 1], V[::-1, 0]])

    # Plot the flight envelope
    plt.plot(V[:, 0], h_vec, color="blue", label=r"Min true airspeed $V_{min}^{true}$")
    plt.plot(V[:, 1], h_vec, color="red", label=r"Max true airspeed $V_{max}^{true}$")
    plt.fill(V_closed, h_closed, color="lightblue", alpha=0.6, label="Flight Envelope")

    # Adding the yellow zone from Vmax to 1.15*Vmax
    V_1_15_max = 1.15 * V[:, 1]
    plt.fill_betweenx(h_vec, V[:, 1], V_1_15_max, color="yellow", alpha=0.5, label=r"Margin zone")
    plt.plot(V_1_15_max, h_vec, color="black")


    # Add flutter velocities if provided
    if V_flutter_pk_min is not None:
        plt.plot(V_flutter_pk_min, h_vec, color="green", linestyle="--", label=r"$V_{flutter}^{pk}$ - Min weight scenario", linewidth=2)
    if V_flutter_pk_max is not None:
        plt.plot(V_flutter_pk_max, h_vec, color="purple", linestyle="--", label=r"$V_{flutter}^{pk}$ - Max weight scenario", linewidth=2)

    # Add horizontal boundaries for the flight envelope
    plt.hlines(0, flight_envelope(0)[0], 1.15*flight_envelope(0)[1], colors="black", linestyles="--", linewidth=1.5)
    plt.hlines(3, flight_envelope(3)[0], 1.15*flight_envelope(3)[1], colors="black", linestyles="--", linewidth=1.5)

    # Plot settings
    plt.xlabel("True airspeed [km/h]", fontsize=15)
    plt.ylabel("Altitude (h) [km]", fontsize=15)
    plt.xlim(-10, 450)
    plt.legend(fontsize=12)
    plt.grid(True)

    # Show plot
    plt.show()

plot_flight_envelope_with_flutter(h_vec, V_min_max, V_flutter_pk_min)
















