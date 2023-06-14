import numpy as np
import matplotlib.pyplot as plt

# leggi i dati dai file
data1 = np.loadtxt("best_final_path_mpi.dat")

# Estrai le coordinate X e Y dai dati
x_coords = data1[:,0]
y_coords = data1[:,1]

# Crea il grafico
for i in range(1, len(x_coords)):
    dx = x_coords[i] - x_coords[i-1]
    dy = y_coords[i] - y_coords[i-1]
    plt.arrow(x_coords[i-1], y_coords[i-1], dx, dy, head_width=0.06, head_length=0.06, fc='blue', ec='blue')


plt.xlabel("X")
plt.ylabel("Y")
plt.title("Percorso dei punti")
plt.grid(True)

# Mostra il grafico
plt.show()
