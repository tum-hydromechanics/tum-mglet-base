// SPDX-FileCopyrightText: Copyright (c) Kitware Inc.
// SPDX-License-Identifier: BSD-3-Clause
#ifdef USE_CATALYST
#include "CatalystAdaptor.h"
#endif
#include "FEDataStructures.h"
#include <iostream>
// #include <mpi.h>

int main(int argc, char* argv[])
{
  // MPI_Init(&argc, &argv);
  // int numRanks(1), myRank(0);
  int myRank = 5;
  int numRanks = 3;
  // MPI_Comm_size(MPI_COMM_WORLD, &numRanks);
  // MPI_Comm_rank(MPI_COMM_WORLD, &myRank);

  unsigned int numberOfAMRLevels = 4;
  std::vector<AMR*> amr;
  for (int rank=0;rank<myRank;rank++){
    amr.push_back(new AMR(numberOfAMRLevels, rank, numRanks));
  }
  std::cout<<"size: "<<amr.size()<<std::endl;

  // The first argument is the program name
#ifdef USE_CATALYST
    // for (int rank=0;rank<myRank;rank++){
  CatalystAdaptor::Initialize(argc, argv);
    // }
#endif
  // keep the number of time steps small since nothing about the grid or fields is changing
  unsigned int numberOfTimeSteps = 2;
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
    //for (int rank=0;rank<myRank;rank++){
    //  std::cout<<"Rank:"<<rank<<",address:"<<&amr[rank]<<std::endl;
    //  CatalystAdaptor::Execute(timeStep, time, *amr[rank], rank);
    //}
    CatalystAdaptor::Execute(timeStep, time, amr, myRank);

#endif
  }
  CatalystAdaptor::Finalize();

#ifdef USE_CATALYST
  // MPI_Finalize();
#endif
  return 0;
}
