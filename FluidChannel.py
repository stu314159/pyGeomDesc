#FluidChannel.py
"""
Class implementation file for the Python class FluidChannel
Depends on vtkHelper module for geometry visualization functionality

"""
import math
import argparse
import numpy as np
from vtkHelper import saveStructuredPointsVTK_ascii as writeVTK

class EmptyChannel:  
    """
     a channel with nothing in it
    """
    def __init__(self,Lo):
        """
         constructor
        """
        self.Lo = Lo

    def get_Lo(self):
        """
         set Lo if need be ?
        """
        return self.Lo

    def get_obstList(self,X,Y,Z):
        """
          for an empty channel - no obstacles 
        """
        return []

class SphereObstruction(EmptyChannel):
    """
     a channel with a sphere obstruction
    """

    def __init__(self,r,x_c,y_c,z_c):
        """
          just need to define the radius and position of the center of the obstacle.
          it is up to the caller to verify that the object will fit within the intended
          channel.  If it does not fit, the obstacle will effectively be
          truncated at the channel boundaries
         
        """
        self.r = r
        self.x_c = x_c
        self.y_c = y_c
        self.z_c = z_c
 
    def get_Lo(self):
        return self.r*2.

    def get_obstList(self,X,Y,Z):
        """
            return a list of all indices all indices within boundary of sphere 
        """

        x = np.array(X); y = np.array(Y); z = np.array(Z);
        dist = (x - self.x_c)**2 + (y - self.y_c)**2 + (z - self.z_c)**2
       
        return list(np.where(dist < self.r**2))
       
class EllipticalScourPit(EmptyChannel):
    """
     a channel with an elliptical scour pit with prescribed properties
     corresponds to case 3 of Bryan's geometry_desc.m
    """

    def __init__(self,x_c,z_c,cyl_rad):
        """
          constructor giving the x and z coordinates of the scour pit along with
          the radius of the cylindrical piling
        """
        self.x_c = x_c
        self.z_c = z_c
        self.cyl_rad = cyl_rad

    def get_Lo(self):
        return self.cyl_rad*2.

    def get_obstList(self,X,Y,Z):
        """
         return a list of all indices of lattice points within the boundaries of the
         scour pit obstacle

        """
       
        ellip_a = 2.*2.*self.cyl_rad
        ellip_b = 2.*self.cyl_rad
        ellip_c = 8.*self.cyl_rad
        ellip_x = self.x_c
        ellip_z = self.z_c + self.cyl_rad
        ellip_y = ellip_b 

        floor_part = np.array(np.where(Y < ellip_b)).flatten()

        dist = (X - self.x_c)**2 + (Z - self.z_c)**2;
        cyl_part = list(np.array(np.where( dist < self.cyl_rad**2)).flatten())

       scour_pit = np.array(np.where( (X - ellip_x)**2/(ellip_a**2) + 
                        (Y - ellip_y)**2/(ellip_b**2) +
                        (Z - ellip_z)**2/(ellip_c**2) <= 1.)).flatten()

        # remove the scour pit from the floor
        obst_list = np.setxor1d(floor_part[:], 
                        np.intersect1d(floor_part[:],scour_pit[:]))


        # then add the cylinder
        obst_list = np.union1d(obst_list[:],cyl_part[:])
        
        return list(obst_list[:])



def fluid_properties(fluid_str):  
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

class FluidChannel:
    def __init__(self,Lx_p=1.,
        Ly_p=1.,
        Lz_p=6.,
        fluid='water', 
        obst=EmptyChannel(1.),
        N_divs = 5,
        wallList=['left','right','top','bottom']):
        """
         class constructor

        """
        self.Lx_p = Lx_p
        self.Ly_p = Ly_p
        self.Lz_p = Lz_p
        self.N_divs = N_divs
        self.fluid = fluid
        self.obst = obst

        # generate the geometry

        Lo = obst.get_Lo()

        self.Ny = math.ceil((N_divs-1)*(Ly_p/Lo))+1
        self.Nx = math.ceil((N_divs-1)*(Lx_p/Lo))+1
        self.Nz = math.ceil((N_divs-1)*(Lz_p/Lo))+1
        self.nnodes = self.Nx*self.Ny*self.Nz
        print "Creating channel with %g lattice points." % self.nnodes
        x = np.linspace(0.,Lx_p,self.Nx).astype(np.float32);
        y = np.linspace(0.,Ly_p,self.Ny).astype(np.float32);
        z = np.linspace(0.,Lz_p,self.Nz).astype(np.float32);
   
        Y,Z,X = np.meshgrid(y,z,x);
    
        self.x = np.reshape(X,self.nnodes)
        self.y = np.reshape(Y,self.nnodes)
        self.z = np.reshape(Z,self.nnodes)

        # get fluid properties from the included fluid library
        self.rho_p, self.nu_p = fluid_properties(fluid)

        # identify inlet and outlet nodes - 
        # require the user to set solid boundaries separately
        self.inlet_list = np.where(self.z==0)
        self.outlet_list = np.where(self.z==Lz_p)
        
        print "Getting obstacle list"
        # get obstacle list
        self.obst_list = self.obst.get_obstList(self.x[:],self.y[:],self.z[:])
        

        print "Generating channel solid boundaries"
        # set channel walls
        self.set_channel_walls(wallList)

        # now eliminate overlap between node lists

        self.inlet_list = np.setxor1d(self.inlet_list[:],
            np.intersect1d(self.inlet_list[:],self.solid_list[:]))
        self.inlet_list = np.setxor1d(self.inlet_list[:],
            np.intersect1d(self.inlet_list[:],self.obst_list[:]))
        
        self.outlet_list = np.setxor1d(self.outlet_list[:],
            np.intersect1d(self.outlet_list[:],self.solid_list[:]))
        self.outlet_list = np.setxor1d(self.outlet_list[:],
            np.intersect1d(self.outlet_list[:],self.obst_list[:]))

        self.obst_list = np.setxor1d(self.obst_list[:],
            np.intersect1d(self.obst_list[:],self.solid_list[:]))
       
        
    def write_bc_vtk(self):
        """
         write node lists to properly formatted VTK files
        """
        print "Creating boundary condition arrays"
        obst_array = np.zeros(self.nnodes)
        obst_array[list(self.obst_list)] = 100.

        #print type(self.inlet_list)
        inlet_array = np.zeros(self.nnodes)
        inlet_array[list(self.inlet_list)] = 200.

        outlet_array = np.zeros(self.nnodes)
        outlet_array[list(self.outlet_list)] = 300.

        solid_array = np.zeros(self.nnodes)
        solid_array[list(self.solid_list)] = 500.
        
        dims = [int(self.Nx), int(self.Ny), int(self.Nz)]
        origin = [0., 0., 0.]
        dx = self.x[1] - self.x[0]
        spacing = [dx, dx, dx] #uniform lattice
        
        print "Writing boundary conditions to VTK files"
        writeVTK(inlet_array,'inlet','inlet.vtk',dims,origin,spacing)
        writeVTK(outlet_array,'outlet','outlet.vtk',dims,origin,spacing)
        writeVTK(obst_array,'obst','obst.vtk',dims,origin,spacing)
        writeVTK(solid_array,'solid','solid.vtk',dims,origin,spacing)


     # must have geometry set first
    def set_channel_walls(self,walls=['left','right','top','bottom']): 
        """
         set up to 4 walls as solid walls for the simulation
        """
        solid_list_a = np.empty(0).flatten()
        solid_list_b = np.empty(0).flatten()
        solid_list_c = np.empty(0).flatten()
        solid_list_d = np.empty(0).flatten()

        for w in walls:
            if w=='right':
                solid_list_a = np.array(np.where((self.x==0.))).flatten()
            elif w=='left':
                solid_list_b = np.array(np.where((self.x == self.Lx_p))).flatten()
            elif w=='top':
                solid_list_d = np.array(np.where((self.y == self.Ly_p))).flatten()
            elif w=='bottom':
                solid_list_c = np.array(np.where((self.y == 0.))).flatten()

        solid_list = np.array(np.union1d(solid_list_a,solid_list_b)); 
        solid_list = np.array(np.union1d(solid_list,solid_list_c))
        self.solid_list = np.array(np.union1d(solid_list,solid_list_d))


   







