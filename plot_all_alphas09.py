import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1.inset_locator import inset_axes, mark_inset
import matplotlib.colors as colors
import matplotlib.cm as mplcm

size_img = 1
plt.rcParams.update({'font.size': 12})
plt.rcParams['figure.figsize'] = [size_img * 10, size_img * 6]
plt.rc('text', usetex=True)
plt.rc('font', family='serif')
plt.tick_params(axis='both', direction='in')
NUM_COLORS = 10
cm = plt.get_cmap('rainbow')
cNorm  = colors.Normalize(vmin=0, vmax=NUM_COLORS+1)
scalarMap = mplcm.ScalarMappable(norm=cNorm, cmap=cm)


# Cargar datos
data_files = ["alpha09/output_090.txt", "alpha091/output_091.txt",
              "alpha092/output_092.txt", "alpha093/output_093.txt",
              "alpha094/output_094.txt", "alpha095/output_095.txt",
              "alpha096/output_096.txt", "alpha097/output_097.txt", 
              "alpha098/output_098.txt", "alpha099/output_099.txt"]

labels = [r"$\alpha = 0.90$", r"$\alpha = 0.91$", 
          r"$\alpha = 0.92$", r"$\alpha = 0.93$", 
          r"$\alpha = 0.94$", r"$\alpha = 0.95$",
          r"$\alpha = 0.96$", r"$\alpha = 0.97$",
          r"$\alpha = 0.98$", r"$\alpha = 0.99$"]

# Crear figura y eje principal
fig, ax = plt.subplots()

ax.set_prop_cycle(color=[scalarMap.to_rgba(NUM_COLORS-i) for i in range(NUM_COLORS)])

# Dibujar las curvas en el eje principal
for i, file in enumerate(data_files):
    data = np.loadtxt(file, delimiter=" ")
    ax.plot(data[:, 0], data[:, 2], lw=0.7, label=labels[i])

# Insertar el zoom en la esquina superior derecha
axins = inset_axes(ax, width="40%", height="30%", loc="upper center")

# Dibujar las mismas curvas en el zoom (sin etiquetas para evitar duplicar la leyenda)
for file in data_files:
    data = np.loadtxt(file, delimiter=" ")
    axins.plot(data[:, 0], data[:, 2], lw=0.8)  # Línea ligeramente más gruesa en el zoom

# Definir la región de zoom
axins.set_xlim(73000, 90000)
axins.set_ylim(0, 0.0002)

# Quitar etiquetas y ticks del zoom
axins.set_xticklabels([])
axins.set_yticklabels([])
axins.set_xticks([])
axins.set_yticks([])

# Dibujar el cuadro rojo con transparencia
mark_inset(ax, axins, loc1=3, loc2=1, fc="red", alpha=0.2, ec="none")

# Dibujar las líneas de conexión en negro sólido
mark_inset(ax, axins, loc1=3, loc2=1, fc="none", ec="k", lw=0.5,alpha=0.5)

# Etiquetas y límites de los ejes
ax.set_xlabel(r"$n$", fontsize=16)
ax.set_ylabel(r"$s_n$", fontsize=16)
ax.set_xlim(-500, 100500)
ax.set_ylim(-0.00005, 0.0023)

# Agregar grid
ax.grid(True, linestyle="--", alpha=0.6)

# Agregar leyenda mejorada
ax.legend(loc="upper right", fontsize=10, frameon=True, edgecolor="black",fancybox=False)

# Guardar la figura
plt.savefig("alpha09.pdf", bbox_inches="tight")
plt.show()
