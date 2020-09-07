
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
      renderView1.ViewSize = [973, 537]
      renderView1.AnnotationColor = [0.0, 0.0, 0.0]
      renderView1.AxesGrid = 'GridAxes3DActor'
      renderView1.OrientationAxesLabelColor = [0.0, 0.0, 0.0]
      renderView1.OrientationAxesOutlineColor = [0.0, 0.0, 0.0]
      renderView1.CenterOfRotation = [0.0, 0.0, 0.5]
      renderView1.StereoType = 0
      renderView1.CameraPosition = [3.0, 0.5, 1.0]
      renderView1.CameraFocalPoint = [-6.0, -3.0, -6.0]
      renderView1.CameraViewUp = [-0.32240347756710497, 0.938373684704941, -0.12454246466932267]
      renderView1.CameraParallelScale = 3.2018951952285404
      renderView1.Background = [1.0, 1.0, 1.0]
      renderView1.OSPRayMaterialLibrary = materialLibrary1

      # init the 'GridAxes3DActor' selected for 'AxesGrid'
      renderView1.AxesGrid.XTitleColor = [0.0, 0.0, 0.0]
      renderView1.AxesGrid.XTitleFontFile = ''
      renderView1.AxesGrid.YTitleColor = [0.0, 0.0, 0.0]
      renderView1.AxesGrid.YTitleFontFile = ''
      renderView1.AxesGrid.ZTitleColor = [0.0, 0.0, 0.0]
      renderView1.AxesGrid.ZTitleFontFile = ''
      renderView1.AxesGrid.GridColor = [0.0, 0.0, 0.0]
      renderView1.AxesGrid.XLabelColor = [0.0, 0.0, 0.0]
      renderView1.AxesGrid.XLabelFontFile = ''
      renderView1.AxesGrid.YLabelColor = [0.0, 0.0, 0.0]
      renderView1.AxesGrid.YLabelFontFile = ''
      renderView1.AxesGrid.ZLabelColor = [0.0, 0.0, 0.0]
      renderView1.AxesGrid.ZLabelFontFile = ''

      # register the view with coprocessor
      # and provide it with information such as the filename to use,
      # how frequently to write the images, etc.
      coprocessor.RegisterView(renderView1,
          filename='fig/Wing_VelocityMagnitude_VolRendering_%t.png', freq=10, fittoscreen=0, magnification=1, width=1920, height=1280, cinema={})
      renderView1.ViewTime = datadescription.GetTime()

      # ----------------------------------------------------------------
      # restore active view
      SetActiveView(renderView1)
      # ----------------------------------------------------------------

      # ----------------------------------------------------------------
      # setup the data processing pipelines
      # ----------------------------------------------------------------

      # create a new 'XML Partitioned Unstructured Grid Reader'
      # create a producer from a simulation input
      input = coprocessor.CreateProducer(datadescription, 'Input')

      # create a new 'Wavefront OBJ Reader'
      wingSurface = WavefrontOBJReader("wingsurface.obj")

      # create a new 'Clip'
      clipX = Clip(Input=input)
      clipX.ClipType = 'Plane'
      clipX.Scalars = ['POINTS', 'Magnitude']
      clipX.Value = 0.776977069017983

      # init the 'Plane' selected for 'ClipType'
      clipX.ClipType.Origin = [-0.2, 0.0, 0.05000000447034836]
      clipX.ClipType.Normal = [-1.0, 0.0, 0.0]

      # create a new 'Clip'
      clipYTop = Clip(Input=clipX)
      clipYTop.ClipType = 'Plane'
      clipYTop.Scalars = ['POINTS', 'Magnitude']
      clipYTop.Value = 0.776977069017983

      # init the 'Plane' selected for 'ClipType'
      clipYTop.ClipType.Origin = [0.0, 0.3, 0.05000000447034836]
      clipYTop.ClipType.Normal = [0.0, 1.0, 0.0]

      # create a new 'Clip'
      clipYBottom = Clip(Input=clipYTop)
      clipYBottom.ClipType = 'Plane'
      clipYBottom.Scalars = ['POINTS', 'Magnitude']
      clipYBottom.Value = 0.776977069017983

      # init the 'Plane' selected for 'ClipType'
      clipYBottom.ClipType.Origin = [0.0, -0.2, 0.05000000447034836]
      clipYBottom.ClipType.Normal = [0.0, -1.0, 0.0]

      # create a new 'Resample To Image'
      resample = ResampleToImage(Input=clipYBottom)
      resample.SamplingDimensions = [1024, 512, 64]
      resample.SamplingBounds = [-0.20000000298023224, 3.0, -0.20000000298023224, 0.30000001192092896, 0.0, 0.10000000894069672]

      # ----------------------------------------------------------------
      # setup the visualization in view 'renderView1'
      # ----------------------------------------------------------------

      # show data from wingSurface
      wingSurfaceDisplay = Show(wingSurface, renderView1)

      # trace defaults for the display properties.
      wingSurfaceDisplay.Representation = 'Surface'
      wingSurfaceDisplay.AmbientColor = [0.0, 0.0, 0.0]
      wingSurfaceDisplay.ColorArrayName = [None, '']
      wingSurfaceDisplay.OSPRayScaleFunction = 'PiecewiseFunction'
      wingSurfaceDisplay.SelectOrientationVectors = 'None'
      wingSurfaceDisplay.ScaleFactor = 0.10002639870217536
      wingSurfaceDisplay.SelectScaleArray = 'None'
      wingSurfaceDisplay.GlyphType = 'Arrow'
      wingSurfaceDisplay.GlyphTableIndexArray = 'None'
      wingSurfaceDisplay.GaussianRadius = 0.005001319935108768
      wingSurfaceDisplay.SetScaleArray = [None, '']
      wingSurfaceDisplay.ScaleTransferFunction = 'PiecewiseFunction'
      wingSurfaceDisplay.OpacityArray = [None, '']
      wingSurfaceDisplay.OpacityTransferFunction = 'PiecewiseFunction'
      wingSurfaceDisplay.DataAxesGrid = 'GridAxesRepresentation'
      wingSurfaceDisplay.SelectionCellLabelFontFile = ''
      wingSurfaceDisplay.SelectionPointLabelFontFile = ''
      wingSurfaceDisplay.PolarAxes = 'PolarAxesRepresentation'

      # init the 'GridAxesRepresentation' selected for 'DataAxesGrid'
      wingSurfaceDisplay.DataAxesGrid.XTitleColor = [0.0, 0.0, 0.0]
      wingSurfaceDisplay.DataAxesGrid.XTitleFontFile = ''
      wingSurfaceDisplay.DataAxesGrid.YTitleColor = [0.0, 0.0, 0.0]
      wingSurfaceDisplay.DataAxesGrid.YTitleFontFile = ''
      wingSurfaceDisplay.DataAxesGrid.ZTitleColor = [0.0, 0.0, 0.0]
      wingSurfaceDisplay.DataAxesGrid.ZTitleFontFile = ''
      wingSurfaceDisplay.DataAxesGrid.GridColor = [0.0, 0.0, 0.0]
      wingSurfaceDisplay.DataAxesGrid.XLabelColor = [0.0, 0.0, 0.0]
      wingSurfaceDisplay.DataAxesGrid.XLabelFontFile = ''
      wingSurfaceDisplay.DataAxesGrid.YLabelColor = [0.0, 0.0, 0.0]
      wingSurfaceDisplay.DataAxesGrid.YLabelFontFile = ''
      wingSurfaceDisplay.DataAxesGrid.ZLabelColor = [0.0, 0.0, 0.0]
      wingSurfaceDisplay.DataAxesGrid.ZLabelFontFile = ''

      # init the 'PolarAxesRepresentation' selected for 'PolarAxes'
      wingSurfaceDisplay.PolarAxes.PolarAxisTitleColor = [0.0, 0.0, 0.0]
      wingSurfaceDisplay.PolarAxes.PolarAxisTitleFontFile = ''
      wingSurfaceDisplay.PolarAxes.PolarAxisLabelColor = [0.0, 0.0, 0.0]
      wingSurfaceDisplay.PolarAxes.PolarAxisLabelFontFile = ''
      wingSurfaceDisplay.PolarAxes.LastRadialAxisTextColor = [0.0, 0.0, 0.0]
      wingSurfaceDisplay.PolarAxes.LastRadialAxisTextFontFile = ''
      wingSurfaceDisplay.PolarAxes.SecondaryRadialAxesTextColor = [0.0, 0.0, 0.0]
      wingSurfaceDisplay.PolarAxes.SecondaryRadialAxesTextFontFile = ''

      # show data from resample
      resampleDisplay = Show(resample, renderView1)

      # get color transfer function/color map for 'velocity'
      velocityLUT = GetColorTransferFunction('velocity')
      velocityLUT.AutomaticRescaleRangeMode = 'Never'
      velocityLUT.RGBPoints = [0.0, 0.001462, 0.000466, 0.013866, 0.005883, 0.002267, 0.00127, 0.01857, 0.011764499999999999, 0.003299, 0.002249, 0.024239, 0.0176475, 0.004547, 0.003392, 0.030909, 0.023528999999999998, 0.006006, 0.004692, 0.038558, 0.029412, 0.007676, 0.006136, 0.046836, 0.035293500000000005, 0.009561, 0.007713, 0.055143, 0.0411765, 0.011663, 0.009417, 0.06346, 0.0470595, 0.013995, 0.011225, 0.071862, 0.052941, 0.016561, 0.013136, 0.080282, 0.058824, 0.019373, 0.015133, 0.088767, 0.0647055, 0.022447, 0.017199, 0.097327, 0.07058849999999998, 0.025793, 0.019331, 0.10593, 0.07647, 0.029432, 0.021503, 0.114621, 0.082353, 0.033385, 0.023702, 0.123397, 0.088236, 0.037668, 0.025921, 0.132232, 0.09411749999999999, 0.042253, 0.028139, 0.141141, 0.1000005, 0.046915, 0.030324, 0.150164, 0.105882, 0.051644, 0.032474, 0.159254, 0.111765, 0.056449, 0.034569, 0.168414, 0.1176465, 0.06134, 0.03659, 0.177642, 0.12352949999999999, 0.066331, 0.038504, 0.186962, 0.12941250000000004, 0.071429, 0.040294, 0.196354, 0.135294, 0.076637, 0.041905, 0.205799, 0.14117699999999997, 0.081962, 0.043328, 0.215289, 0.14705850000000004, 0.087411, 0.044556, 0.224813, 0.1529415, 0.09299, 0.045583, 0.234358, 0.158823, 0.098702, 0.046402, 0.243904, 0.164706, 0.104551, 0.047008, 0.25343, 0.1705875, 0.110536, 0.047399, 0.262912, 0.1764705, 0.116656, 0.047574, 0.272321, 0.1823535, 0.122908, 0.047536, 0.281624, 0.18823499999999999, 0.129285, 0.047293, 0.290788, 0.194118, 0.135778, 0.046856, 0.299776, 0.1999995, 0.142378, 0.046242, 0.308553, 0.20588249999999997, 0.149073, 0.045468, 0.317085, 0.211764, 0.15585, 0.044559, 0.325338, 0.217647, 0.162689, 0.043554, 0.333277, 0.22353, 0.169575, 0.042489, 0.340874, 0.2294115, 0.176493, 0.041402, 0.348111, 0.23529450000000002, 0.183429, 0.040329, 0.354971, 0.241176, 0.190367, 0.039309, 0.361447, 0.24705899999999997, 0.197297, 0.0384, 0.367535, 0.25294049999999996, 0.204209, 0.037632, 0.373238, 0.25882350000000004, 0.211095, 0.03703, 0.378563, 0.26470649999999996, 0.217949, 0.036615, 0.383522, 0.270588, 0.224763, 0.036405, 0.388129, 0.276471, 0.231538, 0.036405, 0.3924, 0.2823525, 0.238273, 0.036621, 0.396353, 0.28823550000000003, 0.244967, 0.037055, 0.400007, 0.2941170000000001, 0.25162, 0.037705, 0.403378, 0.3, 0.258234, 0.038571, 0.406485, 0.305883, 0.26481, 0.039647, 0.409345, 0.3117645, 0.271347, 0.040922, 0.411976, 0.3176475, 0.27785, 0.042353, 0.414392, 0.323529, 0.284321, 0.043933, 0.416608, 0.329412, 0.290763, 0.045644, 0.418637, 0.3352935000000001, 0.297178, 0.04747, 0.420491, 0.3411765, 0.303568, 0.049396, 0.422182, 0.3470595, 0.309935, 0.051407, 0.423721, 0.352941, 0.316282, 0.05349, 0.425116, 0.358824, 0.32261, 0.055634, 0.426377, 0.3647055, 0.328921, 0.057827, 0.427511, 0.3705885, 0.335217, 0.06006, 0.428524, 0.37646999999999997, 0.3415, 0.062325, 0.429425, 0.38235300000000005, 0.347771, 0.064616, 0.430217, 0.388236, 0.354032, 0.066925, 0.430906, 0.3941175, 0.360284, 0.069247, 0.431497, 0.4000005, 0.366529, 0.071579, 0.431994, 0.40588199999999997, 0.372768, 0.073915, 0.4324, 0.41176499999999994, 0.379001, 0.076253, 0.432719, 0.4176465, 0.385228, 0.078591, 0.432955, 0.4235295, 0.391453, 0.080927, 0.433109, 0.4294125, 0.397674, 0.083257, 0.433183, 0.435294, 0.403894, 0.08558, 0.433179, 0.441177, 0.410113, 0.087896, 0.433098, 0.4470585, 0.416331, 0.090203, 0.432943, 0.4529415, 0.422549, 0.092501, 0.432714, 0.458823, 0.428768, 0.09479, 0.432412, 0.46470600000000006, 0.434987, 0.097069, 0.432039, 0.47058749999999994, 0.441207, 0.099338, 0.431594, 0.4764705, 0.447428, 0.101597, 0.43108, 0.4823535, 0.453651, 0.103848, 0.430498, 0.488235, 0.459875, 0.106089, 0.429846, 0.49411799999999995, 0.4661, 0.108322, 0.429125, 0.4999995, 0.472328, 0.110547, 0.428334, 0.5058825, 0.478558, 0.112764, 0.427475, 0.5117640000000001, 0.484789, 0.114974, 0.426548, 0.5176470000000001, 0.491022, 0.117179, 0.425552, 0.52353, 0.497257, 0.119379, 0.424488, 0.5294115, 0.503493, 0.121575, 0.423356, 0.5352944999999999, 0.50973, 0.123769, 0.422156, 0.541176, 0.515967, 0.12596, 0.420887, 0.547059, 0.522206, 0.12815, 0.419549, 0.5529405, 0.528444, 0.130341, 0.418142, 0.5588235, 0.534683, 0.132534, 0.416667, 0.5647065, 0.54092, 0.134729, 0.415123, 0.570588, 0.547157, 0.136929, 0.413511, 0.5764710000000001, 0.553392, 0.139134, 0.411829, 0.5823525, 0.559624, 0.141346, 0.410078, 0.5882354999999999, 0.565854, 0.143567, 0.408258, 0.594117, 0.572081, 0.145797, 0.406369, 0.6, 0.578304, 0.148039, 0.404411, 0.605883, 0.584521, 0.150294, 0.402385, 0.6117644999999999, 0.590734, 0.152563, 0.40029, 0.6176475000000001, 0.59694, 0.154848, 0.398125, 0.623529, 0.603139, 0.157151, 0.395891, 0.629412, 0.60933, 0.159474, 0.393589, 0.6352935, 0.615513, 0.161817, 0.391219, 0.6411764999999999, 0.621685, 0.164184, 0.388781, 0.6470595, 0.627847, 0.166575, 0.386276, 0.652941, 0.633998, 0.168992, 0.383704, 0.658824, 0.640135, 0.171438, 0.381065, 0.6647055, 0.64626, 0.173914, 0.378359, 0.6705884999999999, 0.652369, 0.176421, 0.375586, 0.67647, 0.658463, 0.178962, 0.372748, 0.682353, 0.66454, 0.181539, 0.369846, 0.6882360000000001, 0.670599, 0.184153, 0.366879, 0.6941175000000002, 0.676638, 0.186807, 0.363849, 0.7000005, 0.682656, 0.189501, 0.360757, 0.705882, 0.688653, 0.192239, 0.357603, 0.7117649999999999, 0.694627, 0.195021, 0.354388, 0.7176465, 0.700576, 0.197851, 0.351113, 0.7235294999999999, 0.7065, 0.200728, 0.347777, 0.7294125, 0.712396, 0.203656, 0.344383, 0.735294, 0.718264, 0.206636, 0.340931, 0.741177, 0.724103, 0.20967, 0.337424, 0.7470585, 0.729909, 0.212759, 0.333861, 0.7529414999999999, 0.735683, 0.215906, 0.330245, 0.758823, 0.741423, 0.219112, 0.326576, 0.7647060000000001, 0.747127, 0.222378, 0.322856, 0.7705875, 0.752794, 0.225706, 0.319085, 0.7764705, 0.758422, 0.229097, 0.315266, 0.7823534999999999, 0.76401, 0.232554, 0.311399, 0.788235, 0.769556, 0.236077, 0.307485, 0.794118, 0.775059, 0.239667, 0.303526, 0.7999995, 0.780517, 0.243327, 0.299523, 0.8058825000000001, 0.785929, 0.247056, 0.295477, 0.8117639999999999, 0.791293, 0.250856, 0.29139, 0.817647, 0.796607, 0.254728, 0.287264, 0.8235299999999999, 0.801871, 0.258674, 0.283099, 0.8294115, 0.807082, 0.262692, 0.278898, 0.8352945, 0.812239, 0.266786, 0.274661, 0.8411759999999999, 0.817341, 0.270954, 0.27039, 0.847059, 0.822386, 0.275197, 0.266085, 0.8529405, 0.827372, 0.279517, 0.26175, 0.8588235, 0.832299, 0.283913, 0.257383, 0.8647064999999999, 0.837165, 0.288385, 0.252988, 0.870588, 0.841969, 0.292933, 0.248564, 0.876471, 0.846709, 0.297559, 0.244113, 0.8823524999999999, 0.851384, 0.30226, 0.239636, 0.8882355000000001, 0.855992, 0.307038, 0.235133, 0.894117, 0.860533, 0.311892, 0.230606, 0.8999999999999999, 0.865006, 0.316822, 0.226055, 0.905883, 0.869409, 0.321827, 0.221482, 0.9117645000000001, 0.873741, 0.326906, 0.216886, 0.9176475, 0.878001, 0.33206, 0.212268, 0.9235289999999999, 0.882188, 0.337287, 0.207628, 0.9294120000000001, 0.886302, 0.342586, 0.202968, 0.9352935, 0.890341, 0.347957, 0.198286, 0.9411765, 0.894305, 0.353399, 0.193584, 0.9470594999999999, 0.898192, 0.358911, 0.18886, 0.952941, 0.902003, 0.364492, 0.184116, 0.958824, 0.905735, 0.37014, 0.17935, 0.9647055, 0.90939, 0.375856, 0.174563, 0.9705885000000001, 0.912966, 0.381636, 0.169755, 0.97647, 0.916462, 0.387481, 0.164924, 0.982353, 0.919879, 0.393389, 0.16007, 0.9882359999999999, 0.923215, 0.399359, 0.155193, 0.9941174999999999, 0.92647, 0.405389, 0.150292, 1.0000005, 0.929644, 0.411479, 0.145367, 1.0058819999999997, 0.932737, 0.417627, 0.140417, 1.011765, 0.935747, 0.423831, 0.13544, 1.0176465000000001, 0.938675, 0.430091, 0.130438, 1.0235295, 0.941521, 0.436405, 0.125409, 1.0294124999999998, 0.944285, 0.442772, 0.120354, 1.0352940000000002, 0.946965, 0.449191, 0.115272, 1.041177, 0.949562, 0.45566, 0.110164, 1.0470585, 0.952075, 0.462178, 0.105031, 1.0529415, 0.954506, 0.468744, 0.099874, 1.058823, 0.956852, 0.475356, 0.094695, 1.064706, 0.959114, 0.482014, 0.089499, 1.0705875, 0.961293, 0.488716, 0.084289, 1.0764705, 0.963387, 0.495462, 0.079073, 1.0823535, 0.965397, 0.502249, 0.073859, 1.088235, 0.967322, 0.509078, 0.068659, 1.094118, 0.969163, 0.515946, 0.063488, 1.0999995, 0.970919, 0.522853, 0.058367, 1.1058824999999999, 0.97259, 0.529798, 0.053324, 1.111764, 0.974176, 0.53678, 0.048392, 1.117647, 0.975677, 0.543798, 0.043618, 1.1235300000000001, 0.977092, 0.55085, 0.03905, 1.1294115000000002, 0.978422, 0.557937, 0.034931, 1.1352944999999997, 0.979666, 0.565057, 0.031409, 1.141176, 0.980824, 0.572209, 0.028508, 1.147059, 0.981895, 0.579392, 0.02625, 1.1529405, 0.982881, 0.586606, 0.024661, 1.1588235, 0.983779, 0.593849, 0.02377, 1.1647065, 0.984591, 0.601122, 0.023606, 1.170588, 0.985315, 0.608422, 0.024202, 1.1764709999999998, 0.985952, 0.61575, 0.025592, 1.1823525, 0.986502, 0.623105, 0.027814, 1.1882355, 0.986964, 0.630485, 0.030908, 1.1941169999999999, 0.987337, 0.63789, 0.034916, 1.2, 0.987622, 0.64532, 0.039886, 1.205883, 0.987819, 0.652773, 0.045581, 1.2117645, 0.987926, 0.66025, 0.05175, 1.2176474999999995, 0.987945, 0.667748, 0.058329, 1.2235289999999999, 0.987874, 0.675267, 0.065257, 1.229412, 0.987714, 0.682807, 0.072489, 1.2352934999999998, 0.987464, 0.690366, 0.07999, 1.2411765000000001, 0.987124, 0.697944, 0.087731, 1.2470595, 0.986694, 0.70554, 0.095694, 1.2529409999999999, 0.986175, 0.713153, 0.103863, 1.258824, 0.985566, 0.720782, 0.112229, 1.2647054999999998, 0.984865, 0.728427, 0.120785, 1.2705885000000001, 0.984075, 0.736087, 0.129527, 1.2764699999999998, 0.983196, 0.743758, 0.138453, 1.2823529999999999, 0.982228, 0.751442, 0.147565, 1.288236, 0.981173, 0.759135, 0.156863, 1.2941175000000003, 0.980032, 0.766837, 0.166353, 1.3000005, 0.978806, 0.774545, 0.176037, 1.305882, 0.977497, 0.782258, 0.185923, 1.311765, 0.976108, 0.789974, 0.196018, 1.3176464999999997, 0.974638, 0.797692, 0.206332, 1.3235295, 0.973088, 0.805409, 0.216877, 1.3294124999999999, 0.971468, 0.813122, 0.227658, 1.335294, 0.969783, 0.820825, 0.238686, 1.3411769999999998, 0.968041, 0.828515, 0.249972, 1.3470585000000002, 0.966243, 0.836191, 0.261534, 1.3529415000000002, 0.964394, 0.843848, 0.273391, 1.3588230000000001, 0.962517, 0.851476, 0.285546, 1.364706, 0.960626, 0.859069, 0.29801, 1.3705875, 0.95872, 0.866624, 0.31082, 1.3764705000000002, 0.956834, 0.874129, 0.323974, 1.3823534999999998, 0.954997, 0.881569, 0.337475, 1.3882350000000003, 0.953215, 0.888942, 0.351369, 1.394118, 0.951546, 0.896226, 0.365627, 1.3999994999999998, 0.950018, 0.903409, 0.380271, 1.4058825, 0.948683, 0.910473, 0.395289, 1.411764, 0.947594, 0.917399, 0.410665, 1.417647, 0.946809, 0.924168, 0.426373, 1.4235299999999997, 0.946392, 0.930761, 0.442367, 1.4294115, 0.946403, 0.937159, 0.458592, 1.4352945000000001, 0.946903, 0.943348, 0.47497, 1.441176, 0.947937, 0.949318, 0.491426, 1.4470589999999999, 0.949545, 0.955063, 0.50786, 1.4529405, 0.95174, 0.960587, 0.524203, 1.4588235000000003, 0.954529, 0.965896, 0.540361, 1.4647065000000001, 0.957896, 0.971003, 0.556275, 1.470588, 0.961812, 0.975924, 0.571925, 1.476471, 0.966249, 0.980678, 0.587206, 1.4823524999999997, 0.971162, 0.985282, 0.602154, 1.4882355, 0.976511, 0.989753, 0.61676, 1.494117, 0.982257, 0.994109, 0.631017, 1.5, 0.988362, 0.998364, 0.644924]
      velocityLUT.NanColor = [0.0, 1.0, 0.0]
      velocityLUT.ScalarRangeInitialized = 1.0

      # get opacity transfer function/opacity map for 'velocity'
      velocityPWF = GetOpacityTransferFunction('velocity')
      velocityPWF.Points = [0.0, 0.0, 0.5, 0.0, 0.0, 0.0, 0.5, 0.0, 0.425, 1.0, 0.5, 0.0, 0.85, 0.0, 0.5, 0.0, 1.5, 0.0, 0.5, 0.0]
      velocityPWF.ScalarRangeInitialized = 1

      # trace defaults for the display properties.
      resampleDisplay.Representation = 'Volume'
      resampleDisplay.AmbientColor = [0.0, 0.0, 0.0]
      resampleDisplay.ColorArrayName = ['POINTS', 'velocity']
      resampleDisplay.LookupTable = velocityLUT
      resampleDisplay.OSPRayScaleArray = 'pressure'
      resampleDisplay.OSPRayScaleFunction = 'PiecewiseFunction'
      resampleDisplay.SelectOrientationVectors = 'None'
      resampleDisplay.ScaleFactor = 0.3200000002980232
      resampleDisplay.SelectScaleArray = 'None'
      resampleDisplay.GlyphType = 'Arrow'
      resampleDisplay.GlyphTableIndexArray = 'None'
      resampleDisplay.GaussianRadius = 0.01600000001490116
      resampleDisplay.SetScaleArray = ['POINTS', 'pressure']
      resampleDisplay.ScaleTransferFunction = 'PiecewiseFunction'
      resampleDisplay.OpacityArray = ['POINTS', 'pressure']
      resampleDisplay.OpacityTransferFunction = 'PiecewiseFunction'
      resampleDisplay.DataAxesGrid = 'GridAxesRepresentation'
      resampleDisplay.SelectionCellLabelFontFile = ''
      resampleDisplay.SelectionPointLabelFontFile = ''
      resampleDisplay.PolarAxes = 'PolarAxesRepresentation'
      resampleDisplay.ScalarOpacityUnitDistance = 0.010109172854385859
      resampleDisplay.ScalarOpacityFunction = velocityPWF
      resampleDisplay.Slice = 31

      # init the 'GridAxesRepresentation' selected for 'DataAxesGrid'
      resampleDisplay.DataAxesGrid.XTitleColor = [0.0, 0.0, 0.0]
      resampleDisplay.DataAxesGrid.XTitleFontFile = ''
      resampleDisplay.DataAxesGrid.YTitleColor = [0.0, 0.0, 0.0]
      resampleDisplay.DataAxesGrid.YTitleFontFile = ''
      resampleDisplay.DataAxesGrid.ZTitleColor = [0.0, 0.0, 0.0]
      resampleDisplay.DataAxesGrid.ZTitleFontFile = ''
      resampleDisplay.DataAxesGrid.GridColor = [0.0, 0.0, 0.0]
      resampleDisplay.DataAxesGrid.XLabelColor = [0.0, 0.0, 0.0]
      resampleDisplay.DataAxesGrid.XLabelFontFile = ''
      resampleDisplay.DataAxesGrid.YLabelColor = [0.0, 0.0, 0.0]
      resampleDisplay.DataAxesGrid.YLabelFontFile = ''
      resampleDisplay.DataAxesGrid.ZLabelColor = [0.0, 0.0, 0.0]
      resampleDisplay.DataAxesGrid.ZLabelFontFile = ''

      # init the 'PolarAxesRepresentation' selected for 'PolarAxes'
      resampleDisplay.PolarAxes.PolarAxisTitleColor = [0.0, 0.0, 0.0]
      resampleDisplay.PolarAxes.PolarAxisTitleFontFile = ''
      resampleDisplay.PolarAxes.PolarAxisLabelColor = [0.0, 0.0, 0.0]
      resampleDisplay.PolarAxes.PolarAxisLabelFontFile = ''
      resampleDisplay.PolarAxes.LastRadialAxisTextColor = [0.0, 0.0, 0.0]
      resampleDisplay.PolarAxes.LastRadialAxisTextFontFile = ''
      resampleDisplay.PolarAxes.SecondaryRadialAxesTextColor = [0.0, 0.0, 0.0]
      resampleDisplay.PolarAxes.SecondaryRadialAxesTextFontFile = ''

      # ----------------------------------------------------------------
      # setup color maps and opacity mapes used in the visualization
      # note: the Get..() functions create a new object, if needed
      # ----------------------------------------------------------------

      # ----------------------------------------------------------------
      # finally, restore active source
      SetActiveSource(resample)
      # ----------------------------------------------------------------
    return Pipeline()

  class CoProcessor(coprocessing.CoProcessor):
    def CreatePipeline(self, datadescription):
      self.Pipeline = _CreatePipeline(self, datadescription)

  coprocessor = CoProcessor()
  # these are the frequencies at which the coprocessor updates.
  freqs = {'Input': [10, 10, 10, 10, 10]}
  coprocessor.SetUpdateFrequencies(freqs)
  if requestSpecificArrays:
    arrays = [['pressure', 0], ['temperature', 0], ['velocity', 0]]
    coprocessor.SetRequestedArrays('Input', arrays)
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
coprocessor = CreateCoProcessor()

#--------------------------------------------------------------
# Enable Live-Visualizaton with ParaView and the update frequency
coprocessor.EnableLiveVisualization(False, 1)

# ---------------------- Data Selection method ----------------------

def RequestDataDescription(datadescription):
    "Callback to populate the request for current timestep"
    global coprocessor

    # setup requests for all inputs based on the requirements of the
    # pipeline.
    coprocessor.LoadRequestedData(datadescription)

# ------------------------ Processing method ------------------------

def DoCoProcessing(datadescription):
    "Callback to do co-processing for current timestep"
    global coprocessor

    # Update the coprocessor by providing it the newly generated simulation data.
    # If the pipeline hasn't been setup yet, this will setup the pipeline.
    coprocessor.UpdateProducers(datadescription)

    # Write output data, if appropriate.
    coprocessor.WriteData(datadescription);

    # Write image capture (Last arg: rescale lookup table), if appropriate.
    coprocessor.WriteImages(datadescription, rescale_lookuptable=rescale_lookuptable,
        image_quality=0, padding_amount=imageFileNamePadding)

    # Live Visualization, if enabled.
    coprocessor.DoLiveVisualization(datadescription, "localhost", 22222)
