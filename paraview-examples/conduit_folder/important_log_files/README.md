# saves log files

this files show the output from different conduit basic examples found on this documentation page
https://llnl-conduit.readthedocs.io/en/latest/tutorial_cpp.html

The important of this files is to be able to identify how the information is being sent from conduit to catalyst. To do so, we must be able to understand the structure of the mesh blueprint.


## Obtaining mesh information
There are multiple ways. The best way is to obtain a '.yaml' file containing the information. But this is not always possible due to installation or knowledge-of-functions limitations. 

On the examples we were able to write the mesh by using:
```
conduit::relay::io::blueprint::save_mesh(mesh,
                                        "complete_uniform_mesh_example2",
                                        "json");
```
Note that the library used is not conduit_cpp but conduit. This is because conduit was build and installed separately in order to access the blueprint-examples and explore all capabilities without the dependency of catalyst::conduit_cpp.



For the AMROverlappingExample file
```
std::cout<<mesh.to_yaml()<<std::endl;
```
Note that here we were using catalyst::conduit_cpp

## List of files

CxxOverlappingAMRExample_N2.txt
CxxOverlappingAMRExample_N4.txt
data1.yaml
FEDriver2_yaml_2.txt
FEDriver2_yaml_4_3.txt
FEDriver2_yaml_4.txt
julia_nest_mesh_try1.root
    this is the extention we get for writing the mesh. In another subfolder (named similar as this file) we are going to be able to see 


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
