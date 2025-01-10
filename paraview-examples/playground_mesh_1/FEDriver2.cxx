// SPDX-FileCopyrightText: Copyright (c) Kitware Inc.
// SPDX-License-Identifier: BSD-3-Clause
#ifdef USE_CATALYST
#include "CatalystAdaptor.h"
#endif
#include "FEDataStructures.h"
#include <iostream>
#include <mpi.h>



// using Node =   conduit_cpp::Node;
// using conduit = conduit_cpp;

int main(int argc, char* argv[]){

    /*
// Call on each of 4 MPI ranks.
CatalystAdaptor::conduit_cpp::Node mesh, bp_index;
CatalystAdaptor::conduit_cpp::blueprint::mesh::examples::braid("uniform", 10, 10, 10, mesh);
char domainFile[1024];
sprintf(domainFile, "./bp/bp_%04d.hdf5", rank);
conduit_cpp::relay::io::save(mesh, domainFile, "hdf5");

CatalystAdaptor::conduit_cpp::blueprint::mpi::mesh::generate_index(mesh,
                                              "",
                                              bp_index["blueprint_index/mesh"],
                                              MPI_COMM_WORLD);
bp_index["file_pattern"] = "./bp/bp_%04d.hdf5";
bp_index["number_of_files"] = 4;
bp_index["number_of_trees"] = 4;
bp_index["protocol/name"] = "hdf5";clear
bp_index["protocol/version"] = "0.4.0";
bp_index["tree_pattern"] = "/";
if(rank == 0)
    CatalystAdaptor::conduit_cpp::relay::io::save(bp_index, "bp.root", "hdf5");
}
*/

// Example of a C++ adaptor for a simulation code
// where the simulation code has an overlapping AMR
// grid. The grid in this case is a vtkOverlappingAMR
// data set with the MPI process id specified
// as cell data. Note that in order to see the AMR
// cells (i.e. Surface With Edges representation of
// the data) that the .vtk reader needs to increase
// the `Default Number Of Levels` parameter to
// greater than the default of 1.
  MPI_Init(&argc, &argv);
  int numRanks(1), myRank(0);
  MPI_Comm_size(MPI_COMM_WORLD, &numRanks);
  MPI_Comm_rank(MPI_COMM_WORLD, &myRank);

  unsigned int numberOfAMRLevels = 3;
  AMR amr(numberOfAMRLevels, myRank, numRanks);

  // The first argument is the program name
#ifdef USE_CATALYST
  CatalystAdaptor::Initialize(argc, argv);
#endif
  // keep the number of time steps small since nothing about the grid or fields is changing
  unsigned int numberOfTimeSteps = 3;
  for (unsigned int timeStep = 0; timeStep < numberOfTimeSteps; timeStep++)
  {
    // use a time step length of 0.1
    double time = timeStep * 0.1;
    // the number of AMR levels that each process will generate
    // the grid is one cell "wide" in the X-direction for each MPI process,
    // one cell deep in the Y-direction and numberOfAMRLevels high in the
    // Z-direction, relative to the 0th level cell. Each 0th level cell in the
    // Z-direction is subsequently refined.
#ifdef USE_CATALYST
    CatalystAdaptor::Execute(timeStep, time, amr);
#endif
  }
  CatalystAdaptor::Finalize();

#ifdef USE_CATALYST
  MPI_Finalize();
#endif
  return 0;
}