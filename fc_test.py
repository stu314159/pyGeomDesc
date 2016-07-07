# fc_test.py 
"""
 Test script for implementation of the FluidChannel.py object library

"""

import FluidChannel as fc


#openChannel = fc.FluidChannel() # basic empty default channel
#openChannel.write_bc_vtk()

openChannel2 = fc.FluidChannel(N_divs = 50)
openChannel2.write_bc_vtk()
