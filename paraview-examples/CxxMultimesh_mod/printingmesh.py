from paraview.simple import *

# Greeting to ensure that ctest knows this script is being imported
print("executing catalyst_pipeline")

# registrationName must match the channel name used in the
# 'CatalystAdaptor'.
producer = TrivialProducer(registrationName="input")


def catalyst_execute(info):
    global grid, particles
    producer.UpdatePipeline()
#     print(producer.GetDataInformation().GetDataAssembly())
    assert producer.GetDataInformation().GetNumberOfDataSets() == 2

    gridInfo = producer.GetSubsetDataInformation(0, "//grid", "Hierarchy");
    particlesInfo = producer.GetSubsetDataInformation(0, "//particles", "Hierarchy");

    print("-----------------------------------")
    print("executing (cycle={}, time={})".format(info.cycle, info.time))
    print("grid:")
    print("  bounds:", gridInfo.GetBounds())
    if info.cycle ==1:
        print("Grid info:")
        print(gridInfo)
