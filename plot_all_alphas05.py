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
data_files = ["alpha05/output_050.txt", "alpha051/output_051.txt",
              "alpha052/output_052.txt", "alpha053/output_053.txt",
              "alpha054/output_054.txt", "alpha055/output_055.txt",
              "alpha056/output_056.txt", "alpha057/output_057.txt", 
              "alpha058/output_058.txt", "alpha059/output_059.txt"]

labels = [r"$\alpha = 0.50$", r"$\alpha = 0.51$", 
          r"$\alpha = 0.52$", r"$\alpha = 0.53$", 
          r"$\alpha = 0.54$", r"$\alpha = 0.55$",
          r"$\alpha = 0.56$", r"$\alpha = 0.57$",
          r"$\alpha = 0.58$", r"$\alpha = 0.59$"]

# Graficar
for i, file in enumerate(data_files):
    data = np.loadtxt(file, delimiter=" ")
    plt.plot(data[:, 0], data[:, 2], lw=0.5, label=labels[i])

plt.grid(True, linestyle="--", lw=0.5, alpha=0.6)

plt.xlabel(r"$n$", fontsize=16)
plt.ylabel(r"$s_n$", fontsize=16)
plt.xlim(-500, 100500)
# plt.ylim(-0.00005, 0.0023)

# Agregar leyenda mejorada
plt.legend(loc="upper left", fontsize=10, frameon=True, edgecolor="black",fancybox=False)
plt.savefig("alpha05.pdf", bbox_inches="tight")
plt.show()