import matplotlib.pyplot as plt
import numpy as np
from pmethod import *
from pkmethod import *

# Constants definition
RHO_0 = 1.225 # Density at sea level [kg/m³]
R = 287.058 # Temperature of at sea level [K]
T0 = 288.15 # Specific gas constant [J/(kg.K)]
G = 9.80665  # Gravitational constant [m/s²]

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

def rho(h):
   # This function returns the density of air [kg/m³] depending on a given altitude h [km]
   # International Standard Atmosphere is used
   T = T0 - 0.0065*h*10**(3)
   rho = RHO_0*(T/T0)**(G/(0.0065*R)-1)
   return rho

def get_V_eq(h):
    # This function returns the maximum equivalent airspeed for a given altitude h [km]
    vmin, vmax = flight_envelope(h)
    return vmin*np.sqrt(rho(h)/RHO_0), vmax*np.sqrt(rho(h)/RHO_0)

def plot_flight_envelopes():
    # Define h vector
    h_vec = np.linspace(0, 3, 100)
    V = np.zeros((len(h_vec), 4)) # for each row (h constant) -> store vmin, vmax, vmin_eq, vmax_eq

    # Calculate vmin and vmax
    for i, h in enumerate(h_vec):
        V[i, 0], V[i, 1] = flight_envelope(h)
        V[i, 2], V[i, 3] = get_V_eq(h)

    # Create figure with two subplots
    fig, axes = plt.subplots(2, 1, figsize=(10, 16), constrained_layout=True)

    # First subplot: True airspeed flight envelope
    ax = axes[0]
    h_closed = np.concatenate([h_vec, h_vec[::-1]])  # Combine h and its reverse
    V_closed = np.concatenate([V[:, 1], V[::-1, 0]])  # Combine vmax and reversed vmin
    ax.plot(V[:, 0], h_vec, color="blue", label=r"Min true airspeed $V_{min}^{true}$")
    ax.plot(V[:, 1], h_vec, color="red", label=r"Max true airspeed $V_{max}^{true}$")
    ax.fill(V_closed, h_closed, color="lightblue", alpha=0.6, label="Flight Envelope")
    ax.hlines(0, flight_envelope(0)[0], flight_envelope(0)[1], colors="black", linestyles="--", linewidth=1.5)
    ax.hlines(3, flight_envelope(3)[0], flight_envelope(3)[1], colors="black", linestyles="--", linewidth=1.5)
    ax.set_xlabel("True airspeed [km/h]", fontsize=12)
    ax.set_ylabel("Altitude (h) [km]", fontsize=12)
    ax.set_title("True Airspeed Flight Envelope", fontsize=14)
    ax.set_xlim(-10, 350)
    ax.legend(fontsize=12)
    ax.grid(True)

    # Second subplot: Equivalent airspeed flight envelope
    ax = axes[1]
    h_closed = np.concatenate([h_vec, h_vec[::-1]])  # Combine h and its reverse
    V_closed = np.concatenate([V[:, 3], V[::-1, 2]])  # Combine vmax and reversed vmin
    ax.plot(V[:, 2], h_vec, color="blue", label=r"Min equivalent airspeed $V_{min}^{eq}$")
    ax.plot(V[:, 3], h_vec, color="red", label=r"Max equivalent airspeed $V_{max}^{eq}$")
    ax.plot(V[:, 3] * alpha, h_vec, linestyle="-", color="tab:orange", linewidth=1.5,
            label=r"Minimum flutter velocity required $V_{SL}^{flutter}$")
    ax.fill(V_closed, h_closed, color="lightblue", alpha=0.6, label="Flight Envelope")
    ax.fill_betweenx(h_vec, V[:, 3], V[:, 3] * alpha, color="yellow", alpha=0.5, label="Margin Zone for flutter")
    ax.fill_betweenx(h_vec, V[:, 3] * alpha, 350 * np.ones(len(h_vec)), color="green", alpha=0.5, label="Safety zone for flutter")
    ax.hlines(0, get_V_eq(0)[0], get_V_eq(0)[1], colors="black", linestyles="--", linewidth=1.5)
    ax.hlines(3, get_V_eq(3)[0], get_V_eq(3)[1], colors="black", linestyles="--", linewidth=1.5)
    ax.set_xlabel("Equivalent airspeed [km/h]", fontsize=12)
    ax.set_ylabel("Altitude (h) [km]", fontsize=12)
    ax.set_xlim(-10, 350)
    ax.set_title("Equivalent Airspeed Flight Envelope", fontsize=14)
    ax.legend(fontsize=12, loc="upper center", bbox_to_anchor=(0.5, -0.15), ncol=2)
    ax.grid(True)

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


# --- Computation of dimensionless parameters for SEA LEVEL flutter velocity ----- #
V_eq_max = get_V_eq(0.8)[1]/3.6 # [m/s]

# Case 1 - min weight
mu_min = m_min/(np.pi*RHO_0*b**2*l)
r_min = np.sqrt(Ip_min/(m_min*b**2))
sigma_min = omega_h_min/omega_theta_min

# Case 2 - max weight
mu_max = m_max/(np.pi*RHO_0*b**2*l)
r_max = np.sqrt(Ip_max/(m_max*b**2))
sigma_max = omega_h_max/omega_theta_max

# Definition of the reduced velocity to explore
V_vec = np.linspace(0.01,2,200)

# Results with p-method and p-k method
pmethod_SL_min = pmethod(a=a, e=e_min, mu=mu_min, r=r_min, sigma=sigma_min, V_vec=V_vec)
pmethod_SL_max = pmethod(a=a, e=e_max, mu=mu_max, r=r_max, sigma=sigma_max, V_vec=V_vec)
pkmethod_SL_min = pkmethod(a=a, e=e_min, mu=mu_min, r=r_min, sigma=sigma_min, V_vec=V_vec)
pkmethod_SL_max = pkmethod(a=a, e=e_max, mu=mu_max, r=r_max, sigma=sigma_max, V_vec=V_vec, eps=10**(-8), k0=0.5)

pmethod_SL_min.run()
pmethod_SL_max.run()

pkmethod_SL_min.run()
pkmethod_SL_max.run()

Vrf_min_pmethod, _ = pmethod_SL_min.find_flutter()
Vrf_max_pmethod, _ = pmethod_SL_max.find_flutter()
Vrf_max_pmethod = 1.290 # approximately, algorithm find_flutter does not work for this case

Vrf_min_pkmethod, _ = pkmethod_SL_min.find_flutter()
Vrf_max_pkmethod, _ = pkmethod_SL_max.find_flutter()

# Convert reduced flutter speed into real speed 
Vf_min_pmethod = Vrf_min_pmethod*b*omega_theta_min
Vf_max_pmethod = Vrf_max_pmethod*b*omega_theta_max

Vf_min_pkmethod = Vrf_min_pkmethod*b*omega_theta_min
Vf_max_pkmethod = Vrf_max_pkmethod*b*omega_theta_max


# Check and print for p-method
print("\n")
if Vf_min_pmethod > alpha * V_eq_max:
    print(f"PMETHOD: Flutter velocity for MIN weight ({Vf_min_pmethod:.2f} m/s) respects regulation (>{alpha * V_eq_max:.2f} m/s).\n")
else:
    print(f"PMETHOD: Flutter velocity for MIN weight ({Vf_min_pmethod:.2f} m/s) does NOT respect regulation (>{alpha * V_eq_max:.2f} m/s).\n")

if Vf_max_pmethod > alpha * V_eq_max:
    print(f"PMETHOD: Flutter velocity for MAX weight ({Vf_max_pmethod:.2f} m/s) respects regulation (>{alpha * V_eq_max:.2f} m/s).\n")
else:
    print(f"PMETHOD: Flutter velocity for MAX weight ({Vf_max_pmethod:.2f} m/s) does NOT respect regulation (>{alpha * V_eq_max:.2f} m/s).\n")

# Check and print for pk-method
if Vf_min_pkmethod > alpha * V_eq_max:
    print(f"PKMETHOD: Flutter velocity for MIN weight ({Vf_min_pkmethod:.2f} m/s) respects regulation (>{alpha * V_eq_max:.2f} m/s).\n")
else:
    print(f"PKMETHOD: Flutter velocity for MIN weight ({Vf_min_pkmethod:.2f} m/s) does NOT respect regulation (>{alpha * V_eq_max:.2f} m/s).\n")

if Vf_max_pkmethod > alpha * V_eq_max:
    print(f"PKMETHOD: Flutter velocity for MAX weight ({Vf_max_pkmethod:.2f} m/s) respects regulation (>{alpha * V_eq_max:.2f} m/s).\n")
else:
    print(f"PKMETHOD: Flutter velocity for MAX weight ({Vf_max_pkmethod:.2f} m/s) does NOT respect regulation (>{alpha * V_eq_max:.2f} m/s).\n")
























