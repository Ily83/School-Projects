import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits import mplot3d

fig = plt.figure(figsize=(8,6))
ax3d = plt.axes(projection="3d")

xdata = np.linspace(-5.12,3,100)
ydata = np.linspace(-3,3,100)
x,y = np.meshgrid(xdata,ydata)
a = 1; b = 1; c= 0.8; d= 1;
Z =  a*np.exp(-b*((x - 1)**2 + (y - 2)**2)) + c*np.exp(-d*((x + 1)**2 + (y + 3)**2))

ax3d = plt.axes(projection='3d')
ax3d.plot_surface(x, y, Z,cmap='plasma')
ax3d.set_title('Surface Plot in Matplotlib')
ax3d.set_xlabel('X')
ax3d.set_ylabel('Y')
ax3d.set_zlabel('Z')

plt.show()
