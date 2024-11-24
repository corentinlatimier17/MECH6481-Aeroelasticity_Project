import matplotlib.pyplot as plt
import numpy as np

# Constants definition
RHO_0 = 1.225 # Density at sea level [kg/m³]
R = 287.058 # Temperature of at sea level [K]
T0 = 288.15 # Specific gas constant [J/(kg.K)]
G = 9.80665  # Gravitational constant [m/s²]

# Regulation 
alpha = 1.15

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











