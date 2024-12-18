# script-version: 2.0
# Catalyst state generated using paraview version 5.13.1-749-g7d48519d05
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
renderView1.ViewSize = [1487, 1144]
renderView1.InteractionMode = '2D'
renderView1.AxesGrid = 'Grid Axes 3D Actor'
renderView1.CenterOfRotation = [0.0, 3.725290298461914e-09, 3.725290298461914e-09]
renderView1.StereoType = 'Crystal Eyes'
renderView1.CameraPosition = [0.0, 3.725290298461914e-09, 1.8351886133500763]
renderView1.CameraFocalPoint = [0.0, 3.725290298461914e-09, 3.725290298461914e-09]
renderView1.CameraFocalDisk = 1.0
renderView1.CameraParallelScale = 0.324418935541363
renderView1.LegendGrid = 'Legend Grid Actor'
renderView1.PolarGrid = 'Polar Grid Actor'

SetActiveView(None)

# ----------------------------------------------------------------
# setup view layouts
# ----------------------------------------------------------------

# create new layout object 'Layout #1'
layout1 = CreateLayout(name='Layout #1')
layout1.AssignView(0, renderView1)
layout1.SetSize(1487, 1144)

# ----------------------------------------------------------------
# restore active view
SetActiveView(renderView1)
# ----------------------------------------------------------------

# ----------------------------------------------------------------
# setup the data processing pipelines
# ----------------------------------------------------------------

# create a new 'XML MultiBlock Data Reader'
grid = XMLMultiBlockDataReader(registrationName='grid', FileName=['/home/jakob/uni/hydro-project/tum-mglet-base/tests/Parker_micro/dataout.vtm'])
grid.CellArrayStatus = ['u']

# create a new 'Slice'
slice1 = Slice(registrationName='Slice1', Input=grid)
slice1.SliceType = 'Plane'
slice1.HyperTreeGridSlicer = 'Plane'
slice1.SliceOffsetValues = [0.0]
slice1.PointMergeMethod = 'Uniform Binning'

# init the 'Plane' selected for 'SliceType'
slice1.SliceType.Normal = [0.0, 0.0, 1.0]

# init the 'Plane' selected for 'HyperTreeGridSlicer'
slice1.HyperTreeGridSlicer.Origin = [0.0, 3.725290298461914e-09, 3.725290298461914e-09]

# ----------------------------------------------------------------
# setup the visualization in view 'renderView1'
# ----------------------------------------------------------------

# show data from grid
gridDisplay = Show(grid, renderView1, 'UniformGridRepresentation')

# get 2D transfer function for 'u'
uTF2D = GetTransferFunction2D('u')
uTF2D.ScalarRangeInitialized = 1
uTF2D.Range = [-8.723669052124023, 35.11651611328125, 0.0, 1.0]

# get color transfer function/color map for 'u'
uLUT = GetColorTransferFunction('u')
uLUT.TransferFunction2D = uTF2D
uLUT.RGBPoints = [-8.723669052124023, 0.0564, 0.0564, 0.47, -1.2010439992218016, 0.243, 0.46035, 0.81, 4.362231658082962, 0.356814, 0.745025, 0.954368, 10.22094632321739, 0.6882, 0.93, 0.91791, 13.196423530578613, 0.899496, 0.944646, 0.768657, 17.06426770698166, 0.957108, 0.833819, 0.508916, 22.233683140760423, 0.927521, 0.621439, 0.315357, 28.4370255014801, 0.8, 0.352, 0.16, 35.11651611328125, 0.59, 0.0767, 0.119475]
uLUT.ScalarRangeInitialized = 1.0

# get opacity transfer function/opacity map for 'u'
uPWF = GetOpacityTransferFunction('u')
uPWF.Points = [-8.723669052124023, 0.0, 0.5, 0.0, 35.11651611328125, 1.0, 0.5, 0.0]
uPWF.ScalarRangeInitialized = 1

# trace defaults for the display properties.
gridDisplay.Representation = 'Outline'
gridDisplay.ColorArrayName = ['CELLS', 'u']
gridDisplay.LookupTable = uLUT
gridDisplay.SelectNormalArray = 'None'
gridDisplay.SelectTangentArray = 'None'
gridDisplay.SelectTCoordArray = 'None'
gridDisplay.TextureTransform = 'Transform2'
gridDisplay.OSPRayScaleFunction = 'Piecewise Function'
gridDisplay.Assembly = 'Hierarchy'
gridDisplay.SelectedBlockSelectors = ['']
gridDisplay.SelectOrientationVectors = 'None'
gridDisplay.ScaleFactor = 0.05759999752044678
gridDisplay.SelectScaleArray = 'None'
gridDisplay.GlyphType = 'Arrow'
gridDisplay.GlyphTableIndexArray = 'None'
gridDisplay.GaussianRadius = 0.0028799998760223387
gridDisplay.SetScaleArray = [None, '']
gridDisplay.ScaleTransferFunction = 'Piecewise Function'
gridDisplay.OpacityArray = [None, '']
gridDisplay.OpacityTransferFunction = 'Piecewise Function'
gridDisplay.DataAxesGrid = 'Grid Axes Representation'
gridDisplay.PolarAxes = 'Polar Axes Representation'
gridDisplay.ScalarOpacityUnitDistance = 0.009989573422720137
gridDisplay.ScalarOpacityFunction = uPWF
gridDisplay.TransferFunction2D = uTF2D
gridDisplay.OpacityArrayName = ['CELLS', 'u']
gridDisplay.ColorArray2Name = ['CELLS', 'u']
gridDisplay.SliceFunction = 'Plane'
gridDisplay.Slice = 12
gridDisplay.SelectInputVectors = [None, '']
gridDisplay.WriteLog = ''

# init the 'Plane' selected for 'SliceFunction'
gridDisplay.SliceFunction.Origin = [0.0, 3.725290298461914e-09, 3.725290298461914e-09]

# show data from slice1
slice1Display = Show(slice1, renderView1, 'GeometryRepresentation')

# trace defaults for the display properties.
slice1Display.Representation = 'Surface'
slice1Display.ColorArrayName = ['CELLS', 'u']
slice1Display.LookupTable = uLUT
slice1Display.SelectNormalArray = 'None'
slice1Display.SelectTangentArray = 'None'
slice1Display.SelectTCoordArray = 'None'
slice1Display.TextureTransform = 'Transform2'
slice1Display.OSPRayScaleFunction = 'Piecewise Function'
slice1Display.Assembly = 'Hierarchy'
slice1Display.SelectedBlockSelectors = ['']
slice1Display.SelectOrientationVectors = 'None'
slice1Display.ScaleFactor = 0.05759999752044678
slice1Display.SelectScaleArray = 'None'
slice1Display.GlyphType = 'Arrow'
slice1Display.GlyphTableIndexArray = 'None'
slice1Display.GaussianRadius = 0.0028799998760223387
slice1Display.SetScaleArray = [None, '']
slice1Display.ScaleTransferFunction = 'Piecewise Function'
slice1Display.OpacityArray = [None, '']
slice1Display.OpacityTransferFunction = 'Piecewise Function'
slice1Display.DataAxesGrid = 'Grid Axes Representation'
slice1Display.PolarAxes = 'Polar Axes Representation'
slice1Display.SelectInputVectors = [None, '']
slice1Display.WriteLog = ''

# setup the color legend parameters for each legend in this view

# get color legend/bar for uLUT in view renderView1
uLUTColorBar = GetScalarBar(uLUT, renderView1)
uLUTColorBar.Title = 'u'
uLUTColorBar.ComponentTitle = ''

# set color bar visibility
uLUTColorBar.Visibility = 1

# show color legend
gridDisplay.SetScalarBarVisibility(renderView1, True)

# show color legend
slice1Display.SetScalarBarVisibility(renderView1, True)

# ----------------------------------------------------------------
# setup color maps and opacity maps used in the visualization
# note: the Get..() functions create a new object, if needed
# ----------------------------------------------------------------

# ----------------------------------------------------------------
# setup animation scene, tracks and keyframes
# note: the Get..() functions create a new object, if needed
# ----------------------------------------------------------------

# get the time-keeper
timeKeeper1 = GetTimeKeeper()

# initialize the timekeeper

# get time animation track
timeAnimationCue1 = GetTimeTrack()

# initialize the animation track

# get animation scene
animationScene1 = GetAnimationScene()

# initialize the animation scene
animationScene1.ViewModules = renderView1
animationScene1.Cues = timeAnimationCue1
animationScene1.AnimationTime = 10.0
animationScene1.StartTime = 10.0
animationScene1.EndTime = 11.0
animationScene1.PlayMode = 'Snap To TimeSteps'

# initialize the animation scene

# ----------------------------------------------------------------
# setup extractors
# ----------------------------------------------------------------

# create extractor
pNG1 = CreateExtractor('PNG', renderView1, registrationName='PNG1')
# trace defaults for the extractor.
pNG1.Trigger = 'Time Step'

# init the 'PNG' selected for 'Writer'
pNG1.Writer.FileName = 'RenderView1_{timestep:06d}{camera}.png'
pNG1.Writer.ImageResolution = [1487, 1144]
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

# init the 'Time Step' selected for 'GlobalTrigger'
options.GlobalTrigger.Frequency = 5

# ------------------------------------------------------------------------------
if __name__ == '__main__':
    from paraview.simple import SaveExtractsUsingCatalystOptions
    # Code for non in-situ environments; if executing in post-processing
    # i.e. non-Catalyst mode, let's generate extracts using Catalyst options
    SaveExtractsUsingCatalystOptions(options)
