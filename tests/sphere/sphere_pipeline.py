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
renderView1.ViewSize = [1900, 1120]
renderView1.AxesGrid = 'Grid Axes 3D Actor'
renderView1.CenterOfRotation = [2.7505929470062256, 2.9802322387695312e-08, -0.7502499967813492]
renderView1.StereoType = 'Crystal Eyes'
renderView1.CameraPosition = [2.7505929470062256, 2.9802322387695312e-08, 14.762661567712483]
renderView1.CameraFocalPoint = [2.7505929470062256, 2.9802322387695312e-08, -0.7502499967813492]
renderView1.CameraFocalDisk = 1.0
renderView1.CameraParallelScale = 4.0150369578821445
renderView1.LegendGrid = 'Legend Grid Actor'
renderView1.PolarGrid = 'Polar Grid Actor'

SetActiveView(None)

# ----------------------------------------------------------------
# setup view layouts
# ----------------------------------------------------------------

# create new layout object 'Layout #1'
layout1 = CreateLayout(name='Layout #1')
layout1.AssignView(0, renderView1)
layout1.SetSize(1900, 1120)

# ----------------------------------------------------------------
# restore active view
SetActiveView(renderView1)
# ----------------------------------------------------------------

# ----------------------------------------------------------------
# setup the data processing pipelines
# ----------------------------------------------------------------

# create a new 'XML MultiBlock Data Reader'
grid = XMLMultiBlockDataReader(registrationName='grid', FileName=['/home/jakob/uni/hydro-project/tum-mglet-base/tests/sphere/dataout.vtm'])
grid.CellArrayStatus = ['u', 'v', 'w']

# create a new 'Clip'
clip1 = Clip(registrationName='Clip1', Input=grid)
clip1.ClipType = 'Plane'
clip1.HyperTreeGridClipper = 'Plane'
clip1.Scalars = ['CELLS', 'u']
clip1.Value = 0.7382767051458359

# init the 'Plane' selected for 'ClipType'
clip1.ClipType.Origin = [4.000000029802322, 2.9802322387695312e-08, 2.9802322387695312e-08]
clip1.ClipType.Normal = [0.0, 0.0, 1.0]

# init the 'Plane' selected for 'HyperTreeGridClipper'
clip1.HyperTreeGridClipper.Origin = [4.000000029802322, 2.9802322387695312e-08, 2.9802322387695312e-08]

# create a new 'STL Reader'
spherestl = STLReader(registrationName='sphere.stl', FileNames=['/home/jakob/uni/hydro-project/tum-mglet-base/tests/sphere/sphere.stl'])

# ----------------------------------------------------------------
# setup the visualization in view 'renderView1'
# ----------------------------------------------------------------

# show data from clip1
clip1Display = Show(clip1, renderView1, 'UnstructuredGridRepresentation')

# get 2D transfer function for 'u'
uTF2D = GetTransferFunction2D('u')
uTF2D.ScalarRangeInitialized = 1
uTF2D.Range = [0.3922745883464813, 1.0842788219451904, 0.0, 1.0]

# get color transfer function/color map for 'u'
uLUT = GetColorTransferFunction('u')
uLUT.TransferFunction2D = uTF2D
uLUT.RGBPoints = [0.3922745883464813, 0.0564, 0.0564, 0.47, 0.511016978798151, 0.243, 0.46035, 0.81, 0.5988316240375936, 0.356814, 0.745025, 0.954368, 0.6913096858072578, 0.6882, 0.93, 0.91791, 0.7382767051458359, 0.899496, 0.944646, 0.768657, 0.7993294706593156, 0.957108, 0.833819, 0.508916, 0.8809271498641074, 0.927521, 0.621439, 0.315357, 0.9788450569140911, 0.8, 0.352, 0.16, 1.0842788219451904, 0.59, 0.0767, 0.119475]
uLUT.ScalarRangeInitialized = 1.0

# get opacity transfer function/opacity map for 'u'
uPWF = GetOpacityTransferFunction('u')
uPWF.Points = [0.3922745883464813, 0.0, 0.5, 0.0, 1.0842788219451904, 1.0, 0.5, 0.0]
uPWF.ScalarRangeInitialized = 1

# trace defaults for the display properties.
clip1Display.Representation = 'Surface'
clip1Display.ColorArrayName = ['CELLS', 'u']
clip1Display.LookupTable = uLUT
clip1Display.SelectNormalArray = 'None'
clip1Display.SelectTangentArray = 'None'
clip1Display.SelectTCoordArray = 'None'
clip1Display.TextureTransform = 'Transform2'
clip1Display.OSPRayScaleFunction = 'Piecewise Function'
clip1Display.Assembly = 'Hierarchy'
clip1Display.SelectedBlockSelectors = ['']
clip1Display.SelectOrientationVectors = 'None'
clip1Display.ScaleFactor = 0.4000000059604645
clip1Display.SelectScaleArray = 'None'
clip1Display.GlyphType = 'Arrow'
clip1Display.GlyphTableIndexArray = 'None'
clip1Display.GaussianRadius = 0.020000000298023225
clip1Display.SetScaleArray = [None, '']
clip1Display.ScaleTransferFunction = 'Piecewise Function'
clip1Display.OpacityArray = [None, '']
clip1Display.OpacityTransferFunction = 'Piecewise Function'
clip1Display.DataAxesGrid = 'Grid Axes Representation'
clip1Display.PolarAxes = 'Polar Axes Representation'
clip1Display.ScalarOpacityFunction = uPWF
clip1Display.ScalarOpacityUnitDistance = 0.18898816030037396
clip1Display.OpacityArrayName = ['CELLS', 'u']
clip1Display.SelectInputVectors = [None, '']
clip1Display.WriteLog = ''

# show data from spherestl
spherestlDisplay = Show(spherestl, renderView1, 'GeometryRepresentation')

# trace defaults for the display properties.
spherestlDisplay.Representation = 'Surface'
spherestlDisplay.ColorArrayName = [None, '']
spherestlDisplay.SelectNormalArray = 'None'
spherestlDisplay.SelectTangentArray = 'None'
spherestlDisplay.SelectTCoordArray = 'None'
spherestlDisplay.TextureTransform = 'Transform2'
spherestlDisplay.OSPRayScaleFunction = 'Piecewise Function'
spherestlDisplay.Assembly = ''
spherestlDisplay.SelectedBlockSelectors = ['']
spherestlDisplay.SelectOrientationVectors = 'None'
spherestlDisplay.ScaleFactor = 0.09995689988136292
spherestlDisplay.SelectScaleArray = 'None'
spherestlDisplay.GlyphType = 'Arrow'
spherestlDisplay.GlyphTableIndexArray = 'None'
spherestlDisplay.GaussianRadius = 0.004997844994068146
spherestlDisplay.SetScaleArray = [None, '']
spherestlDisplay.ScaleTransferFunction = 'Piecewise Function'
spherestlDisplay.OpacityArray = [None, '']
spherestlDisplay.OpacityTransferFunction = 'Piecewise Function'
spherestlDisplay.DataAxesGrid = 'Grid Axes Representation'
spherestlDisplay.PolarAxes = 'Polar Axes Representation'
spherestlDisplay.SelectInputVectors = [None, '']
spherestlDisplay.WriteLog = ''

# setup the color legend parameters for each legend in this view

# get color legend/bar for uLUT in view renderView1
uLUTColorBar = GetScalarBar(uLUT, renderView1)
uLUTColorBar.Title = 'u'
uLUTColorBar.ComponentTitle = ''

# set color bar visibility
uLUTColorBar.Visibility = 1

# show color legend
clip1Display.SetScalarBarVisibility(renderView1, True)

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
animationScene1.AnimationTime = 90.0
animationScene1.StartTime = 90.0
animationScene1.EndTime = 91.0

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
pNG1.Writer.ImageResolution = [1900, 1120]
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
options.GlobalTrigger.Frequency = 10

# ------------------------------------------------------------------------------
if __name__ == '__main__':
    from paraview.simple import SaveExtractsUsingCatalystOptions
    # Code for non in-situ environments; if executing in post-processing
    # i.e. non-Catalyst mode, let's generate extracts using Catalyst options
    SaveExtractsUsingCatalystOptions(options)
