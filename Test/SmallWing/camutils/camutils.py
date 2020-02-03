"""Camera utility function for creating camera animations.

This module contains camera utility functions for interactively creating
camera animations in ParaView that can be saved and later used in 
catalyst scripts. paraview.simple and vtk modules need to be accessible.

"""

from paraview.simple import *
import vtk
import os

def CameraToString(cam, name="cam"):
	"""Generates a string with Python code to set up a given VTK Camera. 

The set attributes are position, focal point, view up vector, parallel 
scale and view angle. The output looks like this::

	cam = vtk.vtkCamera()
	cam.SetPosition([...])
	...
	cam.SetParallelScale(...)

Parameters:
    cam (vtkCamera): Camera to be set up.
    name (string): Variable name for the camera. Defaults to "cam"

Returns:
	string: Python code to setup *cam* as variable with *name*.
"""

	pos = cam.GetPosition()
	focalPos = cam.GetFocalPoint()
	viewUp = cam.GetViewUp()
	parallel = cam.GetParallelScale()
	output = ""
	output += "%s = vtk.vtkCamera()\n"%name
	output += "%s.SetPosition([%f, %f, %f])\n"% \
		(name, pos[0], pos[1], pos[2])
	output += "%s.SetFocalPoint([%f, %f, %f])\n"% \
		(name, focalPos[0], focalPos[1], focalPos[2])
	output += "%s.SetViewUp([%f, %f, %f])\n"% \
		(name, viewUp[0], viewUp[1], viewUp[2])
	output += "%s.SetViewAngle(%f)\n"% (name, cam.GetViewAngle())
	output += "%s.SetParallelScale(%f)\n"% \
		(name, cam.GetParallelScale())
	return output

def PrintCamera(cam):
	"""Prints Python code to set up a given VTK Camera.

See function *CameraToString(cam)*.

Parameters:
	cam(vtkCamera): Camera to be set up.
"""

	print(CameraToString(cam))

def PrintCurrentCamera():
	"""Prints Python code to set up the VTK Camera of the current view.

Only makes sense from within a Paraview pipeline, when a view is set up.
See function *CameraToString(cam)*.
"""

	view = GetActiveView()
	cam = view.GetActiveCamera()
	PrintCamera(cam)

def GetCurrentCamera():
	"""Returns a copy of the camera of the current view.

Only makes sense from within a Paraview pipeline, when a view is set up.

Returns:
	vtkCamera: Copy of the camera of the current view.
"""

	view = GetActiveView()
	cam = view.GetActiveCamera()
	newCam = vtk.vtkCamera()
	newCam.DeepCopy(cam)
	return newCam


def SetViewCamera(view, cam):
	"""Sets the camera of the given view to a given VTK Camera.

Parameters:
	view (RenderView): View for which camera is modified.
	cam (vtkCamera): Camera to switch to.
"""
	view.CameraPosition = cam.GetPosition()
	view.CameraFocalPoint = cam.GetFocalPoint()
	view.CameraParallelScale = cam.GetParallelScale()
	view.CameraViewUp = cam.GetViewUp()
	view.CameraViewAngle = cam.GetViewAngle()

def SetCurrentCamera(cam):
	"""Sets the camera of the current view to a given VTK Camera.

See function *SetViewCamera(cam)*

Parameters:
	cam (vtkCamera): Camera to switch to. 
"""
	view = GetActiveView()
	SetViewCamera(view, cam)


class CameraPath:
	"""Class for a camera path to be animated or loaded from file.

A camera path is made up out of cameras at different time points. They
are interpolated through splines and garantueed to go though the given 
camera settings at the given points. Different interpolation schemes 
(e.g. Linear) are theoreticaly possible, but not implemented right now. 
CameraPath objects are meant to be created interactively in Paraview, 
to then be exported and later used in a Catalyst Script. In a Catalyst 
Script one can load one or more CameraPaths with::
	
	from camutils import SetViewCamera, CameraPath
	from camconfig import *

This assumes the camutils folder is the same folder as the Catalyst 
pipeline. All variables in the camconfig file are included. In Catalyst 
the camera path can be used for interpolation in the CoProcessor class 
with::

	animationStart = 50
	animationEnd = 200
	timestep = datadescription.GetTimeStep()
	if (timestep >= animationStart and timestep <= animationEnd):
		cam = camPath.Interpolate(timeStep, animationStart, animationEnd)
		view = coprocessor.Pipeline.renderView1
		SetViewCamera(view, cam)

Attributes:
	interpolator (vtkCameraInterpolator): Interpolater for the cameras.
	camerasByTime (dict): Time and camera information for the path.

"""

	def __init__(self, cameras=[], times=[]):
		"""Initializes a camera path. 

Parameters:
	cameras (vtkCamera[]): List of cameras to be interpolated.
	times (float[]): List of times at which the path assumes the 
		respective camera setting.

"""
		self.interpolator = vtk.vtkCameraInterpolator()
		self.camerasByTime = dict()
		if (len(cameras) != len(times)):
			print("Number of times an cameras do not match")
		for i in range(len(cameras)):
			camera = vtk.vtkCamera()
			camera.DeepCopy(cameras[i])
			self.camerasByTime[times[i]] = cameras[i]
			self.interpolator.AddCamera(times[i], camera)

	def AddCamera(self, time, camera):
		"""Adds the given camera to the path at the given time. 

The camera is added to the iterpolator and to the camera dictionary. If 
there is already a camera setting present at the given time, that one is 
overwritten.

Parameters:
	time (float): Time at which to add the camera. 
	camera (vtkCamera): Camera to be added to the path.

"""

		cam = vtk.vtkCamera()
		# Create a copy of the camera the camera variable could be 
		# reused outside without the values in this path being modified.
		cam.DeepCopy(camera)
		self.interpolator.AddCamera(time, cam)
		self.camerasByTime[time] = cam
		print(self.camerasByTime)

	def AddCurrentCamera(self, time):
		"""Adds the current camera to the path at the given time. 

See function *AddCamera(cam, time)*

Parameters:
	time (float): Time at which to add the camera. 

"""
		view = GetActiveView()
		camera = view.GetActiveCamera()
		self.AddCamera(time, camera)

	def RemoveCamera(self, time):
		"""Removes the camera corresponding to given time from the path. 

The camera is removed from the iterpolator and the camera dictionary. If 
no camera present at the given time, nothing happens.

Parameters:
	time (float): Time at which to remove the camera. 

"""
		self.interpolator.RemoveCamera(time)
		self.camerasByTime.pop(time, None)

	def InterpolateCamera(self, t, tMin=0, tMax=1):
		"""Returns an interpolated camera along the path. 

The given *t* is the parameter for interpolating the camera. By default,
*t* should range between 0 and 1, but other ranges are supported by 
setting *tMin* and *tMax* (e.g. *tMax* could be set to numFrames-1). 
The given parameter is scaled to the range of the minimum and maximum 
times given on the camera path. 

Parameters:
	t (float): Time for which to interpolate the camera. 
	tMin (float): Minimum time for interpolation. Defaults to 0.
	tMax (float): Maximum time for interpolation. Defaults to 1.

Returns:
	vtkCamera: Camera that is interpolated along this path.

"""
		tMaxInterpolator = self.interpolator.GetMaximumT()
		tMinInterpolator = self.interpolator.GetMinimumT()
		time = (t - tMin) / (tMax - tMin) * \
			(tMaxInterpolator - tMinInterpolator) + tMinInterpolator
		camera = vtk.vtkCamera()
		self.interpolator.InterpolateCamera(time, camera)
		return camera

	def ToString(self, name="camPath"):
		"""Generates a string with Python code to set up this path. 

The times and cameras are put into arrays and later assembled by passing
these to a new CameraPath. This setup assumes that the CameraPath is 
directly imported from camutils(by *from camutils import CameraPath*). 
The output looks like this::

	camPath_times = [...]
	camPath_cam0 = vtk.vtkCamera()
	camPath_cam0.SetPosition([-24.066216, 31.326384, -20.081373])
	...
	camPath_cams = [camPath_cam0, camPath_cam1]
	camPath = CameraPath(camPath_cams, camPath_times)

Parameters:
	name (string): Variable name for the camera. Defaults to "camPath".

Returns:
	string: Python code to set up this camera as a variable with *name*.
		Members of this path are set up with suffixes.
"""

		outputTimes = "%s_times = [" % name
		outputCams = "%s_cams = [" % name
		outputCamSettings = ""
		counter = 0
		for time, camera in self.camerasByTime.items():
			outputTimes += "%f, " % time
			camName = "%s_cam%i" % (name, counter)
			outputCams += "%s, " % camName
			outputCamSettings += CameraToString(camera, camName)
			counter += 1
		# Chop of the last ", "
		if counter >= 1:
			outputTimes = outputTimes[:-2]
			outputCams = outputCams[:-2]
		outputTimes += "]\n"
		outputCams += "]\n"
		output = outputTimes + outputCamSettings + outputCams
		output += "%s = CameraPath(%s_cams, %s_times)\n" % \
			(name, name, name)
		return output

	def Print(self):
		"""Prints Python code to set up this camera path.

See function *ToString(self, name)*.
"""
		print(self.ToString())

	def Animate(self, numFrames=100):
		"""Animates the camera path with a given number of frames.

Only makes sense from within a Paraview pipeline, when a view is set up. 
Adds an additional PythonAnmationCue to the animation scene. Any other 
existing PythonAnimationCue is removed so that consecutively calling 
this function does not include older confliction camera path animations.
Then the animation scene is played in *Sequence* mode.

Parameters:
	numFrames (int): Number of Frames the animation runs with. Defaults
		to 100.
"""

		filelocation = os.path.dirname(os.path.abspath(__file__))
		camPathString = self.ToString()
		# Make sure to use the correct indent
		camPathString = camPathString.replace("\n", "\n    ")
		scene = GetAnimationScene()
		pythonCue = PythonAnimationCue()
		pythonCue.Script = """
# Camera interpolation inception
import sys
sys.path.append(r"%s")
from camutils import CameraPath, SetCurrentCamera
import vtk

# Executed once at the beginning of the animation
def start_cue(self):
    self.counter = 0
    self.numFrames = %i
    %s
    self.camPath = camPath

# Executed every time step
def tick(self):
    time = float(self.counter) / (self.numFrames-1)
    cam = self.camPath.InterpolateCamera(time, 0, 1)
    SetCurrentCamera(cam)
    self.counter += 1

# Executed once at the end of the animation
def end_cue(self): pass
""" % (filelocation, numFrames, camPathString)
		# If there are any other Pyhthon cues remove them
		cues = []
		for cue in scene.Cues:
			if not isinstance(cue, type(pythonCue)):
				cues.append(cue)
		scene.Cues = cues
		scene.Cues.append(pythonCue)
		timeTrack = GetTimeTrack()
		timeTrack.UseAnimationTime = 0
		scene.NumberOfFrames = numFrames
		scene.PlayMode = 'Sequence'
		scene.Play()

	def ExportToFile(self, filename="camconfig", varname="camPath", 
			exportMode=0):
		"""Exports this path to a file to be loaded by another script.

The file is saved in the directory of this file. Multiple export modes 
are supported for overwriting existing files, appending to existing 
files (enabling export of multiple paths to the same file) or doing 
neither. No check for duplicate variables is made.

Parameters:
	filename (string): Filename to export to. Defaults to "camconfig". 
	varname (string): Variable name this path is exported with. Defaults
		to "camPath".
	exportMode (int): If 0, create a new file if none exists or append
		if the file exsists. If 1, overwrite the file if it exists. 
		If 2, neither append nor overwrite, only create a new file if 
		none exists yet. 
"""
		
		directory = os.path.dirname(os.path.abspath(__file__))
		fullpath = directory + "/" + filename + ".py"
		if exportMode==2 and os.path.exists(fullpath):
			print("File already exists and overwriting is disabled.",
				"Set exportMode=0 for appending and =1 for overwriting.")
			return
		
		openMode = "a+" if exportMode==0 else "w+"
		print(openMode)
		with open(fullpath, openMode) as file:
			file.seek(0, os.SEEK_END)
			size = file.tell()
			# Only write the import of we are starting a new file 
			if (exportMode==1 or file.tell() == 0):
				imports = """import vtk
from camutils import CameraPath
"""
				for line in imports.splitlines():
					file.write("%s\n" % line)

			file.write("\n")

			camerapath = self.ToString(varname)
			for line in camerapath.splitlines():
				file.write("%s\n" % line)

		print("Written camera path to %s" %(fullpath))
		