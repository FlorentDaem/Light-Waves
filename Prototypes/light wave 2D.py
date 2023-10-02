import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation

mode = 'vecteur'  # 'vecteur' ou 'intensite'

# Paramètres de la simulation
n_steps = 200  # Nombre d'itérations temporelles
n_points_x = 100  # Nombre de points en x
n_points_y = n_points_x  # Nombre de points en y
c = 14.0  # Vitesse de la lumière
dx = 0.1  # Espacement spatial en x
dy = dx  # Espacement spatial en y
dt = 0.005  # Pas de temps


# Calcul du paramètre de stabilité CFL
stability_param = c * dt / (dx)
print(f"Paramètre de stabilité CFL : {stability_param}")

# Vérification de la condition de stabilité
if stability_param <= 1.0:

    # Position de la source
    source_position_x = n_points_x // 3
    source_position_y = n_points_y // 3

    # Fréquence de la source
    source_frequency = 6.0  # Ajustez cette valeur selon vos besoins

    # Initialisation du champ électromagnétique (2D)
    Ex = np.zeros((n_steps, n_points_x, n_points_y))
    Ey = np.zeros((n_steps, n_points_x, n_points_y))

    Ex[0, source_position_x, source_position_y] = 1.0  # Source, oscillation constante

    # Indices optiques (2D)
    n_values = np.ones((n_points_x, n_points_y))
    n_values[:n_points_x // 2, :] = 1.5  # Première moitié en x avec un indice différent
    n_values[n_points_x // 2:, :] = 1.0  # Deuxième moitié en x avec un indice normal

    # Coefficients d'absorption (2D)
    alpha_values = np.zeros((n_points_x, n_points_y))
    alpha_values[:n_points_x // 2, :] = 0.2  # Première moitié en x avec une absorption forte

    # Création de la figure
    fig, ax = plt.subplots(figsize=(8, 8))
    fleche = ax.quiver(np.arange(n_points_x), np.arange(n_points_y), Ex[0], Ey[0], scale=2, color='b')

    # Fonction d'animation
    def animate(i):
        global Ex, Ey
        Ex[i + 1, 1:-1, 1:-1] = (
            (c**2 * dt**2 / (dx**2 * n_values[1:-1, 1:-1]**2)) * (Ex[i, 2:, 1:-1] - 2 * Ex[i, 1:-1, 1:-1] + Ex[i, :-2, 1:-1]) +
            (c**2 * dt**2 / (dy**2 * n_values[1:-1, 1:-1]**2)) * (Ex[i, 1:-1, 2:] - 2 * Ex[i, 1:-1, 1:-1] + Ex[i, 1:-1, :-2]) +
            2 * Ex[i, 1:-1, 1:-1] - Ex[i - 1, 1:-1, 1:-1] - alpha_values[1:-1, 1:-1] * Ex[i, 1:-1, 1:-1]
        )
        Ey[i + 1, 1:-1, 1:-1] = (
            (c**2 * dt**2 / (dx**2 * n_values[1:-1, 1:-1]**2)) * (Ey[i, 2:, 1:-1] - 2 * Ey[i, 1:-1, 1:-1] + Ey[i, :-2, 1:-1]) +
            (c**2 * dt**2 / (dy**2 * n_values[1:-1, 1:-1]**2)) * (Ey[i, 1:-1, 2:] - 2 * Ey[i, 1:-1, 1:-1] + Ey[i, 1:-1, :-2]) +
            2 * Ey[i, 1:-1, 1:-1] - Ey[i - 1, 1:-1, 1:-1] - alpha_values[1:-1, 1:-1] * Ey[i, 1:-1, 1:-1]
        )
        Ex[i + 1, source_position_x, source_position_y] = np.sin(2 * np.pi * source_frequency * i * dt)  # Source oscillante avec la fréquence

        # Calcul de la norme carée du vecteur E
        I = Ex[i + 1]**2 + Ey[i + 1]**2

        # Mise à jour de la figure
        ax.clear()
        # Affichage du champ de vecteurs
        if mode == 'vecteur':
            fleche = ax.quiver(
                np.arange(n_points_x),
                np.arange(n_points_y),
                Ex[i + 1],
                Ey[i + 1],
                color='r',
                scale=2,
                minlength=0
            )
            return fleche,
        # Affichage de l'intensité
        elif mode == 'intensite':
            intensite = ax.imshow(I, extent=[0, n_points_x, 0, n_points_y], origin='lower', cmap='viridis', vmin=0, vmax=10**-3)
            return intensite,

    # Animation
    ani = FuncAnimation(fig, animate, frames=n_steps - 1, interval=50, blit=True)

    ax.set_xlabel('Position en x')
    ax.set_ylabel('Position en y')
    ax.set_title('Champ de vecteurs Ex, Ey')

    plt.show()

else:
    print("La condition de stabilité CFL n'est pas satisfaite. Réduisez le pas de temps (dt) ou augmentez l'espacement spatial (dx).")