# geometry_desc.py
"""

provide a convenient set of tools for describing channel obstruction objects
for use with NFC.  Define a 3D lattice of prescribed density, select
an object and find lattice points within the object.  Additionally, 
determine which lattice points are on the inlet, outlet and channel solid 
boundary surfaces



"""

import math
import argparse
import numpy as np
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
from vtkHelper import saveStructuredPointsVTK_ascii as writeVTK


def brick3D(Lx_p,Ly_p,Lz_p,Lo,Ny_divs):
   """
   return the x, y and z coordinates for all lattice points
   in a uniform, regular, 3D parallelpiped lattice.
   """
   
   Ny = math.ceil((Ny_divs-1)*(Ly_p/Lo))+1
   Nx = math.ceil((Ny_divs-1)*(Lx_p/Lo))+1
   Nz = math.ceil((Ny_divs-1)*(Lz_p/Lo))+1
   nnodes = Nx*Ny*Nz
   
   # compute geometric data only once
   x = np.linspace(0.,Lx_p,Nx).astype(np.float32);
   y = np.linspace(0.,Ly_p,Ny).astype(np.float32);
   z = np.linspace(0.,Lz_p,Nz).astype(np.float32);
   
   Y,Z,X = np.meshgrid(y,z,x);
    
   XX = np.reshape(X,nnodes)
   YY = np.reshape(Y,nnodes)
   ZZ = np.reshape(Z,nnodes)
   
   return (np.array(XX),np.array(YY), np.array(ZZ))

def fluid_properties(fluid_str):  # add more fluids as desired
   """
   Return the physical density and kinematic viscosity for the prescribed
   fluid.
   
   """
   fluid_lib = {'water':(1000., 1.0e-6), 
                'glycol':(965.3,6.216e-4),
                'glycerin':(1260,1.18e-3)}
   if fluid_str in fluid_lib.keys():
     return fluid_lib[fluid_str]
   else:
     print 'valid fluids are:'
     for keys in fluid_lib:
       print " '%s' " % keys
     raise KeyError('invalid fluid specified')


# -- define the channel outer dimensions --- 
Lx_p = 1. # "thickness"
Ly_p = 2. # "height"
Lz_p = 3. # "length"

Ny_divs = 82. # number of divisions for the characteristic length
Lo = 1.
Ny = int(math.ceil((Ny_divs-1)*(Ly_p/Lo))+1)
Nx = int(math.ceil((Ny_divs-1)*(Lx_p/Lo))+1)
Nz = int(math.ceil((Ny_divs-1)*(Lz_p/Lo))+1)
nnodes = Nx*Ny*Nz

rho_p, nu_p = fluid_properties('water')

print 'fluid selected: %s, rho_p = %g, nu_d = %g. \n' % ('water',rho_p, nu_p)



(x, y, z) = brick3D(Lx_p, Ly_p, Lz_p,Lo,Ny_divs)

#for i in range(len(x)):
#   print '%g, %g, %g' % (x[i],y[i],z[i])

x_c = Lx_p/2.; y_c = Ly_p/2.; z_c = Lz_p/2.; r = Lx_p/5.;

# do it the bone-headed way first:
obst_list = []
for i in range(len(x)):
   if ((x[i] - x_c)**2 + (y[i] - y_c)**2 + (z[i] - z_c)**2 < r**2):
     obst_list.append(i)

print 'total nodes = %d ' % nnodes

dims = [Nx,Ny,Nz]
origin = [0,0,0]
dx = x[1]-x[0]; dy = dx; dz = dx;
spacing = [dx,dy,dz]


print 'preparing output for vtk'
obst_array = np.ones(len(x)); 
#for i in obst_list:
#  obst_array[i] = 100.0
#obst_array[obst_list]= 100.0

print 'writing to vtk file'
writeVTK(obst_array,'obstacle','obst.vtk',dims,origin,spacing);

inlet_list = np.where(z==0)
inlet_array = np.ones(len(x))
#inlet_array[inlet_list] = 200.0
#writeVTK(inlet_array,'inlet','inlet.vtk',dims,origin,spacing)


outlet_list = np.where(z==Lz_p)
outlet_array = np.ones(len(x))
#outlet_array[outlet_list] = 300.0
#writeVTK(outlet_array,'outlet','outlet.vtk',dims,origin,spacing)


# make the solid nodes on the x == 0, x == Lx_p, y == 0, and y == Ly_p boundaries

#solid_list = np.where((x == 0) or (x == Lx_p) or (y == 0) or (y == Ly_p))
solid_a = np.where((x==0.))
solid_b = np.where((x == Lx_p))
solid_c = np.where((y == 0.))
solid_d = np.where((y == Ly_p))
solid_list = np.union1d(solid_a,solid_b)
solid_list = np.union1d(solid_list,solid_c)
solid_list = np.union1d(solid_list,solid_d)

solid_array = np.ones(len(x))
solid_array[solid_list] = 500.0
writeVTK(solid_array,'solid','solid.vtk',dims,origin,spacing)

# remove solid nodes from the inlet node list:
inlet_list = np.setxor1d(inlet_list,np.intersect1d(inlet_list,solid_list))
outlet_list = np.setxor1d(outlet_list,np.intersect1d(outlet_list,solid_list))
obst_list = np.setxor1d(obst_list,np.intersect1d(obst_list,solid_list))

inlet_array[inlet_list] = 200.0
outlet_array[outlet_list] = 300.0
obst_array[obst_list] = 100.0

writeVTK(inlet_array,'inlet','inlet.vtk',dims,origin,spacing)
writeVTK(outlet_array,'outlet','outlet.vtk',dims,origin,spacing)
writeVTK(obst_array,'obst','obst.vtk',dims,origin,spacing)


fig = plt.figure()
ax = fig.add_subplot(111,projection='3d',aspect='equal')
ax.scatter(x[obst_list],y[obst_list],z[obst_list],c='r',marker='o')

ax.set_xlabel('X Label')
ax.set_ylabel('Y Label')
ax.set_zlabel('Z Label')
ax.set_xlim3d(0.,Lx_p)
ax.set_ylim3d(0.,Ly_p)
ax.set_zlim3d(0.,Lz_p)
#ax.pbaspect = [Lx_p,Ly_p,Lz_p]
#ax.auto_scale_xyz([0, Lx_p], [0, Ly_p], [0, Lz_p])
plt.show()


