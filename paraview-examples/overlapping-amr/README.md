# Overlapping AMR Example

This is a slightly modified version of `CxxOverlappingAMRExample` from Paraview's Catalyst2 examples. It also allows the user to output a representative dataset by specifying an output path.

## Build Instructions

Export environment variables for catalyst and the paraview implementation
```
export CATALYST_IMPLEMENTATION_PATHS=<Path to Catalyst Library>
export CATALYST_IMPLEMENTATION_NAME=paraview
```
Build with
```
mkdir build
cd build
cmake ..
cmake --build .
```

## Generating a Representative Dataset

From the overlapping-amr directory root, run the example binary while outputting a representative dataset with
```
./build/bin/CxxOverlappingAMRExampleV2 --output amr-data.vthb
```
A new `amr-data` directory and `amr-data.vthb` file will be created. Import the `amr-data.vthb` into paraview and set `Default Number of Levels` to the desired depth to see the different levels of the tree.

## Running a Simulation with Catalyst

From the overlapping-amr directory root, run the example binary along with catalyst
```
./build/bin/CxxOverlappingAMRExampleV2 catalyst_pipeline.py
```
This will output the simulation details with `print` statements in the `catalyst_execute` function.