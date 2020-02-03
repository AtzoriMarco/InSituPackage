
#--------------------------------------------------------------

# Global timestep output options
timeStepToStartOutputAt=0
forceOutputAtFirstCall=False

# Global screenshot output options
imageFileNamePadding=5
rescale_lookuptable=False

# Whether or not to request specific arrays from the adaptor.
requestSpecificArrays=False

# a root directory under which all Catalyst output goes
rootDirectory=''

# makes a cinema D index table
make_cinema_table=False

#--------------------------------------------------------------
# Code generated from cpstate.py to create the CoProcessor.
# paraview version 5.6.2
#--------------------------------------------------------------

from paraview.simple import *
from paraview import coprocessing

import vmslice, vmisosurface, vmvolrender
from camutils import *

import os
print(os.getcwd()+"\n")

camPath = camutils.CameraPath()

# ---------------------- Data Selection method ----------------------

def RequestDataDescription(datadescription):
    "Callback to populate the request for current timestep"
    global vmslice
    global vmisosurface
    global vmvolrender

    # setup requests for all inputs based on the requirements of the
    # pipeline.
    vmslice.RequestDataDescription(datadescription)
    vmisosurface.RequestDataDescription(datadescription)
    vmvolrender.RequestDataDescription(datadescription)


# ------------------------ Processing method ------------------------

def DoCoProcessing(datadescription):
    "Callback to do co-processing for current timestep"
    global vmslice
    global vmisosurface
    global vmvolrender

    timestep = datadescription.GetTimeStep()
    print("A####### %i", timestep)

    if (timestep <= 20):
      #print("VM Slice Camera:", vmslice.coprocessor.Pipeline.renderView1.CameraPosition)
      vmslice.DoCoProcessing(datadescription)
    elif (timestep < 40):
      #print("VM Isosurface Camera:", vmisosurface.coprocessor.Pipeline.renderView1.CameraPosition)
      vmisosurface.DoCoProcessing(datadescription)
    elif (timestep < 60):
      #print("VM Slice Camera:", vmisosurface.coprocessor.Pipeline.renderView1.CameraPosition)
      vmvolrender.DoCoProcessing(datadescription)
