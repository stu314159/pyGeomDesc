# fc_test.py 
"""
 Test script for implementation of the FluidChannel.py object library

"""

import FluidChannel as fc


#openChannel = fc.FluidChannel() # basic empty default channel
#openChannel.write_bc_vtk()

#openChannel2 = fc.FluidChannel(N_divs = 50)
#openChannel2.write_bc_vtk()

#openChannel3 = fc.FluidChannel(Lx_p = 2.,Ly_p = 3., Lz_p = 14.,
#    obst = fc.EmptyChannel(3.))
#openChannel3.write_bc_vtk()

openChannel4 = fc.FluidChannel(wallList = ['left','right','top'], N_divs = 50)
openChannel4.write_bc_vtk()
