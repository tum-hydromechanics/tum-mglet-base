# script-version: 2.0
# Catalyst state generated using paraview version 5.13.20240909
import paraview
paraview.compatibility.major = 5
paraview.compatibility.minor = 13

#### import the simple module from the paraview
from paraview.simple import *
#### disable automatic camera reset on 'Show'
paraview.simple._DisableFirstRenderCameraReset()

# ----------------------------------------------------------------
# setup views used in the visualization
# ----------------------------------------------------------------

# Create a new 'Render View'
renderView1 = CreateView('RenderView')
renderView1.ViewSize = [1858, 1120]
renderView1.AxesGrid = 'Grid Axes 3D Actor'
renderView1.CenterOfRotation = [1.5, 0.5, 1.5]
renderView1.StereoType = 'Crystal Eyes'
renderView1.CameraPosition = [6.361720043489708, 5.361720043489707, 6.361720043489709]
renderView1.CameraFocalPoint = [1.5, 0.5, 1.5]
renderView1.CameraViewUp = [-0.4082482904638631, 0.816496580927726, -0.40824829046386296]
renderView1.CameraFocalDisk = 1.0
renderView1.CameraParallelScale = 2.179449471770337
renderView1.LegendGrid = 'Legend Grid Actor'
renderView1.PolarGrid = 'Polar Grid Actor'

SetActiveView(None)

# ----------------------------------------------------------------
# setup view layouts
# ----------------------------------------------------------------

# create new layout object 'Layout #1'
layout1 = CreateLayout(name='Layout #1')
layout1.AssignView(0, renderView1)
layout1.SetSize(1858, 1120)

# ----------------------------------------------------------------
# restore active view
SetActiveView(renderView1)
# ----------------------------------------------------------------

# ----------------------------------------------------------------
# setup the data processing pipelines
# ----------------------------------------------------------------

# create a new 'XML UniformGrid AMR Reader'
grid = XMLUniformGridAMRReader(registrationName='grid', FileName=['/students/equintana/tum-mglet-base/paraview-examples/CxxOverlappingAMRExample/build/out.vthb'])
grid.CellArrayStatus = ['procid', 'vtkGhostType']
grid.PointArrayStatus = ['otherfield']

# ----------------------------------------------------------------
# setup the visualization in view 'renderView1'
# ----------------------------------------------------------------

# show data from grid
gridDisplay = Show(grid, renderView1, 'AMRRepresentation')

# trace defaults for the display properties.
gridDisplay.Representation = 'Outline'
gridDisplay.ColorArrayName = [None, '']
gridDisplay.SelectNormalArray = 'None'
gridDisplay.SelectTangentArray = 'None'
gridDisplay.SelectTCoordArray = 'None'
gridDisplay.TextureTransform = 'Transform2'
gridDisplay.OSPRayScaleArray = 'otherfield'
gridDisplay.OSPRayScaleFunction = 'Piecewise Function'
gridDisplay.Assembly = 'Hierarchy'
gridDisplay.SelectedBlockSelectors = ['']
gridDisplay.SelectOrientationVectors = 'None'
gridDisplay.ScaleFactor = 0.30000000000000004
gridDisplay.SelectScaleArray = 'None'
gridDisplay.GlyphType = 'Arrow'
gridDisplay.GlyphTableIndexArray = 'None'
gridDisplay.GaussianRadius = 0.015
gridDisplay.SetScaleArray = ['POINTS', 'otherfield']
gridDisplay.ScaleTransferFunction = 'Piecewise Function'
gridDisplay.OpacityArray = ['POINTS', 'otherfield']
gridDisplay.OpacityTransferFunction = 'Piecewise Function'
gridDisplay.DataAxesGrid = 'Grid Axes Representation'
gridDisplay.PolarAxes = 'Polar Axes Representation'
gridDisplay.ScalarOpacityUnitDistance = 2.095540042778318

# init the 'Piecewise Function' selected for 'ScaleTransferFunction'
gridDisplay.ScaleTransferFunction.Points = [0.0, 0.0, 0.5, 0.0, 2.9950041652780257, 1.0, 0.5, 0.0]

# init the 'Piecewise Function' selected for 'OpacityTransferFunction'
gridDisplay.OpacityTransferFunction.Points = [0.0, 0.0, 0.5, 0.0, 2.9950041652780257, 1.0, 0.5, 0.0]

# ----------------------------------------------------------------
# setup animation scene, tracks and keyframes
# note: the Get..() functions create a new object, if needed
# ----------------------------------------------------------------

# get time animation track
timeAnimationCue1 = GetTimeTrack()

# initialize the animation scene

# get the time-keeper
timeKeeper1 = GetTimeKeeper()

# initialize the timekeeper

# initialize the animation track

# get animation scene
animationScene1 = GetAnimationScene()

# initialize the animation scene
animationScene1.ViewModules = renderView1
animationScene1.Cues = timeAnimationCue1
animationScene1.AnimationTime = 0.0

# ----------------------------------------------------------------
# restore active source
SetActiveSource(grid)
# ----------------------------------------------------------------

# ------------------------------------------------------------------------------
# Catalyst options
from paraview import catalyst
options = catalyst.Options()
options.GlobalTrigger = 'Time Step'
options.CatalystLiveTrigger = 'Time Step'

# ------------------------------------------------------------------------------
if __name__ == '__main__':
    from paraview.simple import SaveExtractsUsingCatalystOptions
    # Code for non in-situ environments; if executing in post-processing
    # i.e. non-Catalyst mode, let's generate extracts using Catalyst options
    SaveExtractsUsingCatalystOptions(options)
