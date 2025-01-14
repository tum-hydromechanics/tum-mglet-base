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
renderView1.CameraPosition = [-3.2315992482676434, 4.380893334494689, 7.284427441318329]
renderView1.CameraFocalPoint = [1.4999999999999996, 0.4999999999999998, 1.4999999999999993]
renderView1.CameraViewUp = [0.15189194965635433, 0.8737805861126201, -0.46199169144290636]
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
grid.DefaultNumberOfLevels = 6

# ----------------------------------------------------------------
# setup the visualization in view 'renderView1'
# ----------------------------------------------------------------

# show data from grid
gridDisplay = Show(grid, renderView1, 'AMRRepresentation')

# get 2D transfer function for 'timestep'
timestepTF2D = GetTransferFunction2D('timestep')

# get color transfer function/color map for 'timestep'
timestepLUT = GetColorTransferFunction('timestep')
timestepLUT.TransferFunction2D = timestepTF2D
timestepLUT.RGBPoints = [1.0, 0.0564, 0.0564, 0.47, 1.000041892578125, 0.243, 0.46035, 0.81, 1.0000728737792968, 0.356814, 0.745025, 0.954368, 1.0001055002441406, 0.6882, 0.93, 0.91791, 1.0001220703125, 0.899496, 0.944646, 0.768657, 1.0001436098632812, 0.957108, 0.833819, 0.508916, 1.000172397705078, 0.927521, 0.621439, 0.315357, 1.000206943359375, 0.8, 0.352, 0.16, 1.000244140625, 0.59, 0.0767, 0.119475]
timestepLUT.ScalarRangeInitialized = 1.0

# get opacity transfer function/opacity map for 'timestep'
timestepPWF = GetOpacityTransferFunction('timestep')
timestepPWF.Points = [1.0, 0.0, 0.5, 0.0, 1.000244140625, 1.0, 0.5, 0.0]
timestepPWF.ScalarRangeInitialized = 1

# trace defaults for the display properties.
gridDisplay.Representation = 'Surface With Edges'
gridDisplay.ColorArrayName = ['FIELD', 'timestep']
gridDisplay.LookupTable = timestepLUT
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
gridDisplay.ScalarOpacityFunction = timestepPWF

# init the 'Piecewise Function' selected for 'ScaleTransferFunction'
gridDisplay.ScaleTransferFunction.Points = [0.0, 0.0, 0.5, 0.0, 2.9950041652780257, 1.0, 0.5, 0.0]

# init the 'Piecewise Function' selected for 'OpacityTransferFunction'
gridDisplay.OpacityTransferFunction.Points = [0.0, 0.0, 0.5, 0.0, 2.9950041652780257, 1.0, 0.5, 0.0]

# setup the color legend parameters for each legend in this view

# get color legend/bar for timestepLUT in view renderView1
timestepLUTColorBar = GetScalarBar(timestepLUT, renderView1)
timestepLUTColorBar.Title = 'timestep'
timestepLUTColorBar.ComponentTitle = ''

# set color bar visibility
timestepLUTColorBar.Visibility = 1

# show color legend
gridDisplay.SetScalarBarVisibility(renderView1, True)

# ----------------------------------------------------------------
# setup color maps and opacity maps used in the visualization
# note: the Get..() functions create a new object, if needed
# ----------------------------------------------------------------

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
