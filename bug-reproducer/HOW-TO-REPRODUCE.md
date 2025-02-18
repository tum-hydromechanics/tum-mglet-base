## How to reproduce
- Install Catalyst: https://catalyst-in-situ.readthedocs.io/en/latest/build_and_install.html
- Compile Paraview with PARAVIEW_ENABLE_CATALYST:BOOL=ON
    - full configuration we used in paraview-config.txt
- Go to /minimal-example and follow the commands:
    ```
    mkdir build
    cd build
    cmake ..
    make
    mpirun -n 2 ./bin/reproducer
    ```
    - This will trigger an assertion and segfault later due to some invalid memory reference 
- Expected output:
    - dataout.vthb file correctly interpreting the grid structure
    - the conversion from a conduit node to a vtkOverlappingAMR should handle the case that some levels may not be represented on some MPI rank