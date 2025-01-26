import numpy as np
import matplotlib.pyplot as plt

data1 = np.loadtxt("alpha2.txt",delimiter=" ")
data2 = np.loadtxt("alpha09.txt",delimiter=" ")

x_alpha1 = data1[:,0]
y_alpha1 = data1[:,2]

x_alpha09 = data2[:,0]
y_alpha09 = data2[:,2]

plt.plot(x_alpha1,y_alpha1,label="alpha = 2")
plt.plot(x_alpha09,y_alpha09,label="alpha = 0.9")
plt.legend()
plt.show()