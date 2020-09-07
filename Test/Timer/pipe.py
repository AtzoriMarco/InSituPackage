
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
import time
from mpi4py import MPI


# ----------------------- CoProcessor definition -----------------------

def CreateCoProcessor():
  def _CreatePipeline(coprocessor, datadescription):
    class Pipeline:
      # state file generated using paraview version 5.6.2

      # ----------------------------------------------------------------
      # setup views used in the visualization
      # ----------------------------------------------------------------

      # trace generated using paraview version 5.6.2
      #
      # To ensure correct image size when batch processing, please search 
      # for and uncomment the line `# renderView*.ViewSize = [*,*]`

      #### disable automatic camera reset on 'Show'
      paraview.simple._DisableFirstRenderCameraReset()

      # get the material library
      materialLibrary1 = GetMaterialLibrary()

      # Create a new 'Render View'
      renderView1 = CreateView('RenderView')
      renderView1.ViewSize = [1327, 829]
      renderView1.AxesGrid = 'GridAxes3DActor'
      renderView1.OrientationAxesVisibility = 0
      renderView1.CenterOfRotation = [0.5000734329223633, 0.0, 0.05000000447034836]
      renderView1.StereoType = 0
      renderView1.CameraPosition = [2.167763566500773, 0.34826958060563895, 0.6160839114657211]
      renderView1.CameraFocalPoint = [-7.354350171537613, -3.3470051309966404, -6.363903005657959]
      renderView1.CameraViewUp = [-0.2898743543208173, 0.9509519829749983, -0.10799622578131131]
      renderView1.CameraParallelScale = 3.2018951952285404
      renderView1.Background = [0.31999694819562063, 0.3400015259021897, 0.4299992370489052]
      renderView1.OSPRayMaterialLibrary = materialLibrary1

      # init the 'GridAxes3DActor' selected for 'AxesGrid'
      renderView1.AxesGrid.XTitleFontFile = ''
      renderView1.AxesGrid.YTitleFontFile = ''
      renderView1.AxesGrid.ZTitleFontFile = ''
      renderView1.AxesGrid.XLabelFontFile = ''
      renderView1.AxesGrid.YLabelFontFile = ''
      renderView1.AxesGrid.ZLabelFontFile = ''

      # register the view with coprocessor
      # and provide it with information such as the filename to use,
      # how frequently to write the images, etc.
      coprocessor.RegisterView(renderView1,
          filename='WingLambda2Iso_%t.png', freq=50, fittoscreen=0, magnification=1, width=1920, height=1280, cinema={})
      renderView1.ViewTime = datadescription.GetTime()

      # Create a new 'Render View'
      renderView2 = CreateView('RenderView')
      renderView2.ViewSize = [1327, 829]
      renderView2.AxesGrid = 'GridAxes3DActor'
      renderView2.CenterOfRotation = [0.7403411054983735, 0.03242843225598335, 0.05000000447034836]
      renderView2.StereoType = 0
      renderView2.CameraPosition = [1.8284104995356767, 0.6774269348570306, 0.5479496071052073]
      renderView2.CameraFocalPoint = [-0.26942860250081746, -0.9801616353957846, -0.7121910220624629]
      renderView2.CameraViewUp = [-0.31851098383607845, 0.7952677677483379, -0.5158487479449276]
      renderView2.CameraParallelScale = 0.7650054918926504
      renderView2.Background = [0.32, 0.34, 0.43]
      renderView2.OSPRayMaterialLibrary = materialLibrary1

      # init the 'GridAxes3DActor' selected for 'AxesGrid'
      renderView2.AxesGrid.XTitleFontFile = ''
      renderView2.AxesGrid.YTitleFontFile = ''
      renderView2.AxesGrid.ZTitleFontFile = ''
      renderView2.AxesGrid.XLabelFontFile = ''
      renderView2.AxesGrid.YLabelFontFile = ''
      renderView2.AxesGrid.ZLabelFontFile = ''

      # register the view with coprocessor
      # and provide it with information such as the filename to use,
      # how frequently to write the images, etc.

      # ----------------------------------------------------------------
      # restore active view
      SetActiveView(renderView1)
      # ----------------------------------------------------------------

      # ----------------------------------------------------------------
      # setup the data processing pipelines
      # ----------------------------------------------------------------

      # create a new 'XML Partitioned Unstructured Grid Reader'
      # create a producer from a simulation input
      input = coprocessor.CreateProducer(datadescription, 'input')

      # create a new 'Calculator'
      magnitudeCalc = Calculator(Input=input)
      magnitudeCalc.ResultArrayName = 'Magnitude'
      magnitudeCalc.Function = 'mag(velocity)'


      # create a new 'Contour'
      lambda2Iso = Contour(Input=input)
      lambda2Iso.ContourBy = ['POINTS', 'temperature']
      lambda2Iso.Isosurfaces = [-200.0]
      lambda2Iso.PointMergeMethod = 'Uniform Binning'

      # create a new 'Contour'
      wingGeo = Contour(Input=magnitudeCalc)
      wingGeo.ContourBy = ['POINTS', 'Magnitude']
      wingGeo.Isosurfaces = [1e-05]
      wingGeo.PointMergeMethod = 'Uniform Binning'

      # ----------------------------------------------------------------
      # setup the visualization in view 'renderView1'
      # ----------------------------------------------------------------

      # show data from lambda2Iso
      lambda2IsoDisplay = Show(lambda2Iso, renderView1)

      # trace defaults for the display properties.
      lambda2IsoDisplay.Representation = 'Surface'
      lambda2IsoDisplay.ColorArrayName = ['POINTS', '']
      lambda2IsoDisplay.DiffuseColor = [0.8196078431372549, 1.0, 0.4823529411764706]
      lambda2IsoDisplay.OSPRayScaleArray = 'Normals'
      lambda2IsoDisplay.OSPRayScaleFunction = 'PiecewiseFunction'
      lambda2IsoDisplay.SelectOrientationVectors = 'None'
      lambda2IsoDisplay.ScaleFactor = 0.3001077304594219
      lambda2IsoDisplay.SelectScaleArray = 'None'
      lambda2IsoDisplay.GlyphType = 'Arrow'
      lambda2IsoDisplay.GlyphTableIndexArray = 'None'
      lambda2IsoDisplay.GaussianRadius = 0.015005386522971094
      lambda2IsoDisplay.SetScaleArray = ['POINTS', 'Normals']
      lambda2IsoDisplay.ScaleTransferFunction = 'PiecewiseFunction'
      lambda2IsoDisplay.OpacityArray = ['POINTS', 'Normals']
      lambda2IsoDisplay.OpacityTransferFunction = 'PiecewiseFunction'
      lambda2IsoDisplay.DataAxesGrid = 'GridAxesRepresentation'
      lambda2IsoDisplay.SelectionCellLabelFontFile = ''
      lambda2IsoDisplay.SelectionPointLabelFontFile = ''
      lambda2IsoDisplay.PolarAxes = 'PolarAxesRepresentation'

      # init the 'GridAxesRepresentation' selected for 'DataAxesGrid'
      lambda2IsoDisplay.DataAxesGrid.XTitleFontFile = ''
      lambda2IsoDisplay.DataAxesGrid.YTitleFontFile = ''
      lambda2IsoDisplay.DataAxesGrid.ZTitleFontFile = ''
      lambda2IsoDisplay.DataAxesGrid.XLabelFontFile = ''
      lambda2IsoDisplay.DataAxesGrid.YLabelFontFile = ''
      lambda2IsoDisplay.DataAxesGrid.ZLabelFontFile = ''

      # init the 'PolarAxesRepresentation' selected for 'PolarAxes'
      lambda2IsoDisplay.PolarAxes.PolarAxisTitleFontFile = ''
      lambda2IsoDisplay.PolarAxes.PolarAxisLabelFontFile = ''
      lambda2IsoDisplay.PolarAxes.LastRadialAxisTextFontFile = ''
      lambda2IsoDisplay.PolarAxes.SecondaryRadialAxesTextFontFile = ''

      # show data from wingGeo
      wingGeoDisplay = Show(wingGeo, renderView1)

      # trace defaults for the display properties.
      wingGeoDisplay.Representation = 'Surface'
      wingGeoDisplay.ColorArrayName = [None, '']
      wingGeoDisplay.OSPRayScaleArray = 'Normals'
      wingGeoDisplay.OSPRayScaleFunction = 'PiecewiseFunction'
      wingGeoDisplay.SelectOrientationVectors = 'None'
      wingGeoDisplay.ScaleFactor = 0.10003012699598912
      wingGeoDisplay.SelectScaleArray = 'None'
      wingGeoDisplay.GlyphType = 'Arrow'
      wingGeoDisplay.GlyphTableIndexArray = 'None'
      wingGeoDisplay.GaussianRadius = 0.005001506349799456
      wingGeoDisplay.SetScaleArray = ['POINTS', 'Normals']
      wingGeoDisplay.ScaleTransferFunction = 'PiecewiseFunction'
      wingGeoDisplay.OpacityArray = ['POINTS', 'Normals']
      wingGeoDisplay.OpacityTransferFunction = 'PiecewiseFunction'
      wingGeoDisplay.DataAxesGrid = 'GridAxesRepresentation'
      wingGeoDisplay.SelectionCellLabelFontFile = ''
      wingGeoDisplay.SelectionPointLabelFontFile = ''
      wingGeoDisplay.PolarAxes = 'PolarAxesRepresentation'

      # init the 'GridAxesRepresentation' selected for 'DataAxesGrid'
      wingGeoDisplay.DataAxesGrid.XTitleFontFile = ''
      wingGeoDisplay.DataAxesGrid.YTitleFontFile = ''
      wingGeoDisplay.DataAxesGrid.ZTitleFontFile = ''
      wingGeoDisplay.DataAxesGrid.XLabelFontFile = ''
      wingGeoDisplay.DataAxesGrid.YLabelFontFile = ''
      wingGeoDisplay.DataAxesGrid.ZLabelFontFile = ''

      # init the 'PolarAxesRepresentation' selected for 'PolarAxes'
      wingGeoDisplay.PolarAxes.PolarAxisTitleFontFile = ''
      wingGeoDisplay.PolarAxes.PolarAxisLabelFontFile = ''
      wingGeoDisplay.PolarAxes.LastRadialAxisTextFontFile = ''
      wingGeoDisplay.PolarAxes.SecondaryRadialAxesTextFontFile = ''

      # ----------------------------------------------------------------
      # setup the visualization in view 'renderView2'
      # ----------------------------------------------------------------

      # show data from wingGeo
      wingGeoDisplay_1 = Show(wingGeo, renderView2)

      # trace defaults for the display properties.
      wingGeoDisplay_1.Representation = 'Surface'
      wingGeoDisplay_1.ColorArrayName = [None, '']
      wingGeoDisplay_1.OSPRayScaleArray = 'Normals'
      wingGeoDisplay_1.OSPRayScaleFunction = 'PiecewiseFunction'
      wingGeoDisplay_1.SelectOrientationVectors = 'None'
      wingGeoDisplay_1.ScaleFactor = 0.10003012699598912
      wingGeoDisplay_1.SelectScaleArray = 'None'
      wingGeoDisplay_1.GlyphType = 'Arrow'
      wingGeoDisplay_1.GlyphTableIndexArray = 'None'
      wingGeoDisplay_1.GaussianRadius = 0.005001506349799456
      wingGeoDisplay_1.SetScaleArray = ['POINTS', 'Normals']
      wingGeoDisplay_1.ScaleTransferFunction = 'PiecewiseFunction'
      wingGeoDisplay_1.OpacityArray = ['POINTS', 'Normals']
      wingGeoDisplay_1.OpacityTransferFunction = 'PiecewiseFunction'
      wingGeoDisplay_1.DataAxesGrid = 'GridAxesRepresentation'
      wingGeoDisplay_1.SelectionCellLabelFontFile = ''
      wingGeoDisplay_1.SelectionPointLabelFontFile = ''
      wingGeoDisplay_1.PolarAxes = 'PolarAxesRepresentation'

      # init the 'GridAxesRepresentation' selected for 'DataAxesGrid'
      wingGeoDisplay_1.DataAxesGrid.XTitleFontFile = ''
      wingGeoDisplay_1.DataAxesGrid.YTitleFontFile = ''
      wingGeoDisplay_1.DataAxesGrid.ZTitleFontFile = ''
      wingGeoDisplay_1.DataAxesGrid.XLabelFontFile = ''
      wingGeoDisplay_1.DataAxesGrid.YLabelFontFile = ''
      wingGeoDisplay_1.DataAxesGrid.ZLabelFontFile = ''

      # init the 'PolarAxesRepresentation' selected for 'PolarAxes'
      wingGeoDisplay_1.PolarAxes.PolarAxisTitleFontFile = ''
      wingGeoDisplay_1.PolarAxes.PolarAxisLabelFontFile = ''
      wingGeoDisplay_1.PolarAxes.LastRadialAxisTextFontFile = ''
      wingGeoDisplay_1.PolarAxes.SecondaryRadialAxesTextFontFile = ''

      # ----------------------------------------------------------------
      # finally, restore active source
      SetActiveSource(wingGeo)
      # ----------------------------------------------------------------
    return Pipeline()

  class CoProcessor(coprocessing.CoProcessor):
    def CreatePipeline(self, datadescription):
      self.Pipeline = _CreatePipeline(self, datadescription)

  coprocessor = CoProcessor()
  # these are the frequencies at which the coprocessor updates.
  freqs = {'input': [1]}
  coprocessor.SetUpdateFrequencies(freqs)
  if requestSpecificArrays:
    arrays = [['temperature', 0], ['velocity', 0]]
    coprocessor.SetRequestedArrays('input', arrays)
  coprocessor.SetInitialOutputOptions(timeStepToStartOutputAt,forceOutputAtFirstCall)

  if rootDirectory:
      coprocessor.SetRootDirectory(rootDirectory)

  if make_cinema_table:
      coprocessor.EnableCinemaDTable()

  return coprocessor


#--------------------------------------------------------------
# Global variable that will hold the pipeline for each timestep
# Creating the CoProcessor object, doesn't actually create the ParaView pipeline.
# It will be automatically setup when coprocessor.UpdateProducers() is called the
# first time.
timee1 = time.time()
coprocessor = CreateCoProcessor()
timee2 = time.time()
comm = MPI.COMM_WORLD
rank=comm.Get_rank()
deltalT1 = timee2-timee1
if rank==0:
    with open("rank0timer.txt", "a+") as fileR0:
        fileR0.write("t_crtCop="+str(deltalT1)+"\n")
elif rank==1:
    with open("rank1timer.txt", "a+") as fileR1:
        fileR1.write("t_crtCop="+str(deltalT1)+"\n")        
        #print("t_createCop=",deltalT1)


#--------------------------------------------------------------
# Enable Live-Visualizaton with ParaView and the update frequency
coprocessor.EnableLiveVisualization(False, 1)

# ---------------------- Data Selection method ----------------------

def RequestDataDescription(datadescription):
    "Callback to populate the request for current timestep"
    global coprocessor

    # setup requests for all inputs based on the requirements of the
    # pipeline
    timel1=time.time()
    coprocessor.LoadRequestedData(datadescription)
    timel2=time.time()
    coprocessor.DoLiveVisualization(datadescription, "localhost", 22222)
    step=datadescription.GetTimeStep()
    comm = MPI.COMM_WORLD
    rank=comm.Get_rank()
    deltalT1 = timel2-timel1
    if rank==0:
        with open("rank0timer.txt", "a+") as fileR0:
            fileR0.write("Stp="+str(step)+", t_ldData="+str(deltalT1)+"\n")
    elif rank==1:
        with open("rank1timer.txt", "a+") as fileR1:
            fileR1.write("Stp="+str(step)+", t_ldData="+str(deltalT1)+"\n")
        #print("Step ",step,", t_ldData=",deltalT1)




# ------------------------ Processing method ------------------------

def DoCoProcessing(datadescription):
    "Callback to do co-processing for current timestep"
    global coprocessor
    time0 = time.time()
    # Update the coprocessor by providing it the newly generated simulation data.
    # If the pipeline hasn't been setup yet, this will setup the pipeline.
    coprocessor.UpdateProducers(datadescription)
    time1=time.time()
    # Write output data, if appropriate.
    coprocessor.WriteData(datadescription);
    time2=time.time()
    # Write image capture (Last arg: rescale lookup table), if appropriate.
    coprocessor.WriteImages(datadescription, rescale_lookuptable=rescale_lookuptable,
        image_quality=0, padding_amount=imageFileNamePadding)
    time3=time.time()

    deltaT = time3-time0
    deltaT1 = time1-time0
    deltaT2 = time2-time1
    deltaT3 = time3-time2
    # Live Visualization, if enabled.
    coprocessor.DoLiveVisualization(datadescription, "localhost", 22222)
    step=datadescription.GetTimeStep()
    comm = MPI.COMM_WORLD
    rank=comm.Get_rank()
    if rank==0:
        with open("rank0timer.txt", "a+") as fileR0:
            fileR0.write("Stp="+str(step)+", t_cop="+str(deltaT)+', t_sv='+str(deltaT3)+
                ', t_updt='+str(deltaT1)+', t_wrt='+str(deltaT2)+"\n")
    elif rank==1:
        with open("rank1timer.txt", "a+") as fileR1:
            fileR1.write("Stp="+str(step)+", t_cop="+str(deltaT)+', t_sv='+str(deltaT3)+
                ', t_updt='+str(deltaT1)+', t_wrt='+str(deltaT2)+"\n")                
        #print("Step ",step,", t_cop=",deltaT,', t_sv=',deltaT3,', t_updt=',
        #        deltaT1,', t_wrt=',deltaT2)



