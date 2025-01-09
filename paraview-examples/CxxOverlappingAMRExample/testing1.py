# script-version: 2.0
# Catalyst state generated using paraview version 5.13.1-627-g56f006c7fe
import paraview
paraview.compatibility.major = 5
paraview.compatibility.minor = 13   

#### import the simple module from the paraview
from paraview.simple import *
#### disable automatic camera reset on 'Show'
paraview.simple._DisableFirstRenderCameraReset()
print('changed')
# ----------------------------------------------------------------
# setup views used in the visualization
# ----------------------------------------------------------------

# Create a new 'Render View'
renderView1 = CreateView('RenderView')
renderView1.ViewSize = [1396, 726]
renderView1.AxesGrid = 'Grid Axes 3D Actor'
renderView1.CenterOfRotation = [2.5, 0.5, 1.5]
renderView1.StereoType = 'Crystal Eyes'
renderView1.CameraPosition = [-4.104634517044713, 7.964624244632029, 7.092849546306956]
renderView1.CameraFocalPoint = [2.4999999999999996, 0.5000000000000046, 1.500000000000001]
renderView1.CameraViewUp = [0.6633805229653105, 0.7251823913894793, -0.18449059859219807]
renderView1.CameraFocalDisk = 1.0
renderView1.CameraParallelScale = 2.958039891549808
renderView1.LegendGrid = 'Legend Grid Actor'
renderView1.PolarGrid = 'Polar Grid Actor'

SetActiveView(None)

# ----------------------------------------------------------------
# setup view layouts
# ----------------------------------------------------------------

# create new layout object 'Layout #1'
layout1 = CreateLayout(name='Layout #1')
layout1.AssignView(0, renderView1)
layout1.SetSize(1396, 726)

# ----------------------------------------------------------------
# restore active view
SetActiveView(renderView1)
# ----------------------------------------------------------------

# ----------------------------------------------------------------
# setup the data processing pipelines
# ----------------------------------------------------------------

# create a new 'XML UniformGrid AMR Reader'
grid = XMLUniformGridAMRReader(registrationName='grid', FileName=['/home/eduardoq/git/ParaviewCatalyst/catalyst2/Examples/Catalyst2/CxxOverlappingAMRExample/build/data.vthb'])
grid.CellArrayStatus = ['procid', 'vtkGhostType']
grid.PointArrayStatus = ['otherfield']
grid.DefaultNumberOfLevels = 1

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
timestepLUT.RGBPoints = [4.0, 0.0564, 0.0564, 0.47, 4.0001675703125, 0.243, 0.46035, 0.81, 4.000291495117187, 0.356814, 0.745025, 0.954368, 4.000422000976562, 0.6882, 0.93, 0.91791, 4.00048828125, 0.899496, 0.944646, 0.768657, 4.000574439453125, 0.957108, 0.833819, 0.508916, 4.000689590820312, 0.927521, 0.621439, 0.315357, 4.0008277734375, 0.8, 0.352, 0.16, 4.0009765625, 0.59, 0.0767, 0.119475]
timestepLUT.ScalarRangeInitialized = 1.0

# get opacity transfer function/opacity map for 'timestep'
timestepPWF = GetOpacityTransferFunction('timestep')
timestepPWF.Points = [4.0, 0.0, 0.5, 0.0, 4.0009765625, 1.0, 0.5, 0.0]
timestepPWF.ScalarRangeInitialized = 1

# trace defaults for the display properties.
gridDisplay.Representation = 'Wireframe'
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
gridDisplay.ScaleFactor = 0.5
gridDisplay.SelectScaleArray = 'None'
gridDisplay.GlyphType = 'Arrow'
gridDisplay.GlyphTableIndexArray = 'None'
gridDisplay.GaussianRadius = 0.025
gridDisplay.SetScaleArray = ['POINTS', 'otherfield']
gridDisplay.ScaleTransferFunction = 'Piecewise Function'
gridDisplay.OpacityArray = ['POINTS', 'otherfield']
gridDisplay.OpacityTransferFunction = 'Piecewise Function'
gridDisplay.DataAxesGrid = 'Grid Axes Representation'
gridDisplay.PolarAxes = 'Polar Axes Representation'
gridDisplay.ScalarOpacityUnitDistance = 2.398852817515996
gridDisplay.ScalarOpacityFunction = timestepPWF

# init the 'Piecewise Function' selected for 'ScaleTransferFunction'
gridDisplay.ScaleTransferFunction.Points = [0.0, 0.0, 0.5, 0.0, 4.921060994002885, 1.0, 0.5, 0.0]

# init the 'Piecewise Function' selected for 'OpacityTransferFunction'
gridDisplay.OpacityTransferFunction.Points = [0.0, 0.0, 0.5, 0.0, 4.921060994002885, 1.0, 0.5, 0.0]

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
# setup extractors
# ----------------------------------------------------------------

# create extractor
pNG1 = CreateExtractor('PNG', renderView1, registrationName='PNG1')
# trace defaults for the extractor.
pNG1.Trigger = 'Time Step'

# init the 'PNG' selected for 'Writer'
pNG1.Writer.FileName = 'RenderView1_{timestep:06d}{camera}.png'
pNG1.Writer.ImageResolution = [1396, 726]
pNG1.Writer.Format = 'PNG'

# ----------------------------------------------------------------
# restore active source
SetActiveSource(pNG1)
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
