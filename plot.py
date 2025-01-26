import numpy as np
import matplotlib.pyplot as plt

data = np.loadtxt("output.txt",delimiter=" ")

x = data[:,0]
y = data[:,2]

plt.plot(x,y)
plt.savefig("fig.pdf")