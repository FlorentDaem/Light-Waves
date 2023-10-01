import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation

# Paramètres de la simulation
n_steps = 400  # Nombre d'itérations temporelles
n_points = 10**2  # Nombre de points d'espace
c = 10**1  # Vitesse de la lumière
dx = 10**-1  # Espacement spatial
dt = 10**-2  # Pas de temps

# Position de la source
source_position = n_points // 3

# Fréquence de la source
source_frequency = 5.0  # Ajustez cette valeur selon vos besoins

# Initialisation du champ électromagnétique
E = np.zeros((n_steps, n_points))
E[0, source_position] = 0  # Source, oscillation constante

# Indices optiques
n_values = np.ones(n_points)
n_values[n_points // 2:] = 1.5  # Deuxième moitié avec un indice différent

# Coefficients d'absorption
alpha_values = np.ones(n_points) * 10**-3
alpha_values[n_points // 2:] = 10**-2  # Deuxième moitié avec une absorption forte

# Création de la figure
fig, ax = plt.subplots()
line, = ax.plot(np.linspace(0, n_points * dx, n_points), E[0, :])

# Ajout d'une ligne verticale pour la délimitation
ax.axvline(x=n_points // 2 * dx, color='gray', linestyle='--', linewidth=0.5)

# Fonction d'animation
def animate(i):
    global E
    E[i + 1, 1:-1] = (c**2 * dt**2 / (dx**2 * n_values[1:-1]**2)) * (E[i, 2:] - 2 * E[i, 1:-1] + E[i, :-2]) + 2 * E[i, 1:-1] - E[i - 1, 1:-1] - alpha_values[1:-1] * E[i, 1:-1]
    E[i + 1, source_position] += 0.1*np.sin(2 * np.pi * source_frequency * i * dt)  # Source oscillante avec la fréquence

    line.set_ydata(E[i, :])
    return line,

# Animation
ani = FuncAnimation(fig, animate, frames=n_steps - 1, interval=50, blit=True, repeat=False)

plt.ylim(-2, 2)
plt.xlabel('Position')
plt.ylabel('Champ électromagnétique')
plt.show()
