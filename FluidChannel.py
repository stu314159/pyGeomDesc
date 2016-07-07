#FluidChannel.py
"""
Class implementation file for the Python class FluidChannel
Depends on vtkHelper module for geometry visualization functionality

"""
import math
import argparse
import numpy as np
from vtkHelper import saveStructuredPointsVTK_ascii as writeVTK

Class EmptyChannel():  # if nothing else, it's a worthwhile default
    """

    """
    def __init__(self,Lo):
        """

        """
        self.Lo = Lo

    def get_Lo(self):
        "

        "
        return self.Lo

    def get_obst(self,X,Y,Z):
        "
        
        "
        return []



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

Class FluidChannel:
    def __init__(self,Lx_p=1.,Ly_p=1.,Lz_p=6.,fluid='water', obst=EmptyChannel(1.),N_divs = 5):
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

        self.Ny = math.ceil((Ny_divs-1)*(Ly_p/Lo))+1
        self.Nx = math.ceil((Ny_divs-1)*(Lx_p/Lo))+1
        self.Nz = math.ceil((Ny_divs-1)*(Lz_p/Lo))+1
        nnodes = Nx*Ny*Nz
   
        
        x = np.linspace(0.,Lx_p,self.Nx).astype(np.float32);
        y = np.linspace(0.,Ly_p,self.Ny).astype(np.float32);
        z = np.linspace(0.,Lz_p,self.Nz).astype(np.float32);
   
        Y,Z,X = np.meshgrid(y,z,x);
    
        self.x = np.reshape(X,nnodes)
        self.y = np.reshape(Y,nnodes)
        self.z = np.reshape(Z,nnodes)

        # get fluid properties from the included fluid library
        self.rho_p, self.nu_p = fluid_properties(fluid)

        # identify inlet and outlet nodes - require the user to set solid boundaries separately
        self.inlet_list = np.where(z==0)
        self.outlet_list = np.where(z==Lz_p)
        
        


    def set_channel_walls(self,walls=['left','right','top','bottom']): # must have geometry set first
        """
         set up to 4 walls as solid walls for the simulation
        """
        solid_list_a = []
        solid_list_b = []
        solid_list_c = []
        solid_list_d = []

        for w in walls:
            if w=='right':
                solid_list_a = np.array(np.where((self.x==0.))).flatten()
            elif w=='left':
                solid_list_b = np.array(np.where((self.x == self.Lx_p))).flatten()
            elif w=='top':
                solid_list_d = np.array(np.where((self.y == self.Ly_p))).flatten()
            elif w=='bottom':
                solid_list_c = np.array(np.where((y == 0.))).flatten()

        solid_list = np.array(np.union1d(solid_list_a,solid_list_b)); 
        solid_list = np.array(np.union1d(solid_list,solid_list_c))
        self.solid_list = np.array(np.union1d(solid_list,solid_list_d))


   







