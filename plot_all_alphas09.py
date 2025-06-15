import numpy as np
import matplotlib.pyplot as plt

size_img = 1
plt.rcParams.update({'font.size': 12})
plt.rcParams['figure.figsize'] = [size_img * 10, size_img * 6]
plt.rc('text', usetex=True)
plt.rc('font', family='serif')
plt.tick_params(axis='both', direction='in')

# Cargar datos (aquí se asume que cada archivo tiene datos tipo columnas con x, y, etc.)
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

# Graficar
for i, file in enumerate(data_files):
    data = np.loadtxt(file, delimiter=" ")
    plt.plot(data[:, 0], data[:, 2], lw=0.5, label=labels[i])

plt.grid(True, linestyle="--", lw=0.5, alpha=0.6)

plt.xlabel(r"$n$", fontsize=16)
plt.ylabel(r"$s_n$", fontsize=16)
plt.xlim(-500, 100500)
plt.ylim(-0.00005, 0.0023)

# Agregar leyenda mejorada
plt.legend(loc="upper right", fontsize=10, frameon=True, edgecolor="black",fancybox=False)
plt.savefig("alpha09.pdf", bbox_inches="tight")
plt.show()