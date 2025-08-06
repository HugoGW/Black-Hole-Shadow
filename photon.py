import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
from numba import njit, prange

# ----------------------------
# Constants and Initial Setup
# ----------------------------
Rs = 1.0               # Schwarzschild radius (in arbitrary units)
x0 = 20.0              # Initial x-position of all photons
vx0, vy0 = -1.0, 0.0   # Initial velocity vector: moving left
t_max = 100            # Total simulation time
npts = 5000            # Number of time steps
t_vals = np.linspace(0, t_max, npts)
dt = t_vals[1] - t_vals[0]

# Photon wall and animation parameters
y_lim = 10             # Vertical extent of the photon wall
n_photons = 1000       # Number of photons in the wall
fade_duration = 30     # Number of frames during which photons fade after being absorbed

# ---------------------------------------------
# Geodesic equations and RK4 integrator (Numba)
# ---------------------------------------------

@njit
def geodesic_eq_numba(y, Rs):
    """
    Returns the derivatives for the geodesic equations in Schwarzschild coordinates.
    y = [r, ψ, dr/dt, dψ/dt]
    """
    r, psi, drdt, dpsidt = y
    d2rdt2 = r * dpsidt**2 - 1.5 * Rs * dpsidt**2
    d2psidt2 = -2 * drdt * dpsidt / r
    return np.array([drdt, dpsidt, d2rdt2, d2psidt2])

@njit
def rk4_step(y, dt, Rs):
    """
    Runge-Kutta 4th order integrator step for solving ODEs.
    """
    k1 = geodesic_eq_numba(y, Rs)
    k2 = geodesic_eq_numba(y + 0.5 * dt * k1, Rs)
    k3 = geodesic_eq_numba(y + 0.5 * dt * k2, Rs)
    k4 = geodesic_eq_numba(y + dt * k3, Rs)
    return y + (dt / 6.0) * (k1 + 2*k2 + 2*k3 + k4)

@njit(parallel=True)
def compute_trajectories(x0, vx0, vy0, y_offsets, Rs, t_max, npts):
    """
    Compute the trajectories of photons using the geodesic equation,
    parallelized with Numba.
    """
    dt = t_max / (npts - 1)
    trajectories = np.zeros((len(y_offsets), npts, 2))
    cut_indices = np.zeros(len(y_offsets), dtype=np.int32)

    for i in prange(len(y_offsets)):
        y0 = y_offsets[i]
        r0 = np.sqrt(x0**2 + y0**2)
        psi0 = np.arctan2(y0, x0)
        vr0 = (x0 * vx0 + y0 * vy0) / r0
        vpsi0 = (x0 * vy0 - y0 * vx0) / (r0**2)
        y_vec = np.array([r0, psi0, vr0, vpsi0])

        cut_idx = npts + 1000  # Default value: photon not absorbed

        for j in range(npts):
            r, psi = y_vec[0], y_vec[1]
            x = r * np.cos(psi)
            y = r * np.sin(psi)
            trajectories[i, j, 0] = x
            trajectories[i, j, 1] = y

            if r < Rs and cut_idx == npts + 1000:
                cut_idx = j  # Mark where photon crosses the event horizon

            y_vec = rk4_step(y_vec, dt, Rs)

        cut_indices[i] = cut_idx

    return trajectories, cut_indices

# -----------------------------
# Compute all photon trajectories
# -----------------------------
print("Computing photon trajectories... this may take a few seconds.")
y_offsets = np.linspace(-y_lim, y_lim, n_photons)
trajectories, cut_indices = compute_trajectories(x0, vx0, vy0, y_offsets, Rs, t_max, npts)

# -----------------------------
# Plot and animation setup
# -----------------------------
fig, ax = plt.subplots(figsize=(8, 8))
fig.patch.set_facecolor('gray')
ax.set_facecolor('gray')
ax.set_xlim(-10, 20)
ax.set_ylim(-y_lim - 1, y_lim + 1)
ax.set_aspect('equal')
ax.set_title("Photon Deflection by a Black Hole (Schwarzschild)", color='white')

# Hide axes
ax.axis('off')

# Draw event horizon (black disk)
event_horizon = plt.Circle((0, 0), Rs, color='black')
ax.add_artist(event_horizon)

# Draw black hole shadow (in dim gray)
R_shadow = (np.sqrt(27) / 2) * Rs
shadow = plt.Circle((0, 0), R_shadow, facecolor='none', edgecolor='dimgrey', linewidth=1.5)
ax.add_patch(shadow)

# Draw photon sphere (in white)
R_photon = 1.5 * Rs
photon_sphere = plt.Circle((0, 0), R_photon, facecolor='none', edgecolor='white', linewidth=1.2)
ax.add_patch(photon_sphere)

# Parallel limit lines showing the shadow boundary
x_line_start = x0
x_line_end = 0
ax.plot([x_line_start, x_line_end], [R_shadow, R_shadow], color='dimgrey', linewidth=1)
ax.plot([x_line_start, x_line_end], [-R_shadow, -R_shadow], color='dimgrey', linewidth=1)

# Photon trajectory alpha value
alpha = 0.2

# Initialize photon lines (in yellow)
lines = [ax.plot([], [], color='yellow', alpha=alpha, linewidth=1)[0] for _ in range(n_photons)]

# -----------------------------
# Animation functions
# -----------------------------

def init():
    for line in lines:
        line.set_data([], [])
        line.set_alpha(0)
    return lines

def animate(frame):
    for i, line in enumerate(lines):
        traj = trajectories[i]
        cut_idx = cut_indices[i]

        if frame < cut_idx:
            line.set_data(traj[:frame, 0], traj[:frame, 1])
            line.set_alpha(alpha)
        elif cut_idx <= frame < cut_idx + fade_duration:
            fade_frame = frame - cut_idx
            new_alpha = max(0, alpha * (1 - fade_frame / fade_duration))
            line.set_data(traj[:cut_idx, 0], traj[:cut_idx, 1])
            line.set_alpha(new_alpha)
        else:
            line.set_data([], [])
            line.set_alpha(0)
    return lines

# -----------------------------
# Create and display animation
# -----------------------------
anim = FuncAnimation(fig, animate, init_func=init,
                     frames=npts, interval=1, blit=True)

# Optional: Save the animation as MP4 using ffmpeg
output_path = "/home/hugo-alexandre/pCloudDrive/Python/code_python_1m2p/Gravitational simulations/General Relativity and Cosmology/black hole/BH_shadow/bh_shadow.mp4"
anim.save(output_path, writer='ffmpeg', fps=30, dpi=200)

plt.show()
