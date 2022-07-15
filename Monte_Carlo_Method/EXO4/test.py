import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from mpl_toolkits import mplot3d

 

f2 = lambda x,y: 20+x**2+y**2-10*(np.cos(2*np.pi*x)+np.cos(2*np.pi*y))
xdata = np.linspace(-5.12,5.12,100)
ydata = np.linspace(-5.12,5.12,100)
x,y = np.meshgrid(xdata,ydata)
z =  f2(x,y)

fig4,ax4=plt.subplots(1,1)
# cp2 = ax4.contourf(x, y, z,cmap='plasma')
cp2 = ax4.pcolormesh(x, y, z,cmap='plasma',shading='gouraud')
fig4.colorbar(cp2) 
plt.show()