import numpy as np
import matplotlib.pyplot as plt

# Parameter values
R0 = 5  # m (radius)
omega0 = np.pi  # rad/s (angular velocity)
dt = 0.01  # s (time step)
N = 1000  # number of time steps

# Exact solution
te = np.arange(0, 2 * np.pi, dt)  # s (time)
psie = omega0 * te  # rad (angular position)
xA1e = R0 * np.cos(psie)  # m (x-coordinate)
xA2e = R0 * np.sin(psie)  # m (y-coordinate)

# Plot exact trajectory
plt.figure(figsize=(8, 6))
plt.plot(xA1e, xA2e, '-r', linewidth=2)
plt.xlabel('$x_{A_1}$ (m)')
plt.ylabel('$x_{A_2}$ (m)')
plt.title('Exact Trajectory')
plt.axis('equal')
plt.grid(True)
plt.show()

# IMU simulation with noise
SNR_dB = 30  # Signal-to-noise ratio (dB)
noise_std = np.sqrt(10 ** (-SNR_dB / 10))  # Standard deviation of noise

psidote = omega0 * np.ones(len(te)) + np.random.normal(0, noise_std, len(te))  # rad/s (angular velocity)
a1e = (-R0 * omega0 ** 2) * np.ones(len(te)) + np.random.normal(0, noise_std, len(te))  # m/s^2 (acceleration)
a2e = np.zeros(len(te)) + np.random.normal(0, noise_std, len(te))  # m/s^2 (acceleration)

# Plot IMU data
plt.figure(figsize=(12, 8))
plt.subplot(3, 1, 1)
plt.plot(te, psidote, '-b', linewidth=2)
plt.xlabel('Time (s)')
plt.ylabel('Angular velocity $\omega_3$ (rad/s)')
plt.title('IMU Data - Noisy')
plt.grid(True)

plt.subplot(3, 1, 2)
plt.plot(te, a1e, '-b', linewidth=2)
plt.xlabel('Time (s)')
plt.ylabel('Specific force $f_1$ (m/s$^2$)')
plt.grid(True)

plt.subplot(3, 1, 3)
plt.plot(te, a2e, '-b', linewidth=2)
plt.xlabel('Time (s)')
plt.ylabel('Specific force $f_2$ (m/s$^2$)')
plt.grid(True)

plt.tight_layout()
plt.show()
