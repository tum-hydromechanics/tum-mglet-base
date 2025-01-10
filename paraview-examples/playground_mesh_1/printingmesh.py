from paraview.simple import *

# Greeting to ensure that ctest knows this script is being imported
print("executing catalyst_pipeline")

# ------------------------------------------------------------------------------
# Catalyst options
from paraview import catalyst
options = catalyst.Options()
options.GlobalTrigger = 'TimeStep'
options.CatalystLiveTrigger = 'TimeStep'

# registrationName must match the channel name used in the
# 'CatalystAdaptor'.
producer = TrivialProducer(registrationName="grid")

def catalyst_execute(info):
    global producer
    producer.UpdatePipeline()

    # assert producer.GetDataInformation().GetNumberOfDataSets() == 2

    gridInfo = producer.GetSubsetDataInformation(0, "//grid", "Hierarchy");
    # particlesInfo = producer.GetSubsetDataInformation(0, "//particles", "Hierarchy");

    if info.cycle ==1:
        print("Grid info:")
        print(gridInfo)



        print("-----------------------------------")
        print("executing (cycle={}, time={})".format(info.cycle, info.time))
        print("bounds:", producer.GetDataInformation().GetBounds())
        print("procid:", producer.CellData["procid"].GetRange(-1))
        print("vtkGhostType:", producer.CellData["vtkGhostType"].GetRange(-1))
