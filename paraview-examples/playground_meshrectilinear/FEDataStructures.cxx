// SPDX-FileCopyrightText: Copyright (c) Kitware Inc.
// SPDX-License-Identifier: BSD-3-Clause
#include "FEDataStructures.h"
#include <math.h>
#include <iostream>
#include <vector>
#include <algorithm>
///////////////////////////////////////////////////////////////////////////////
//AMR starts here

AMR::AMR(int numberOfAMRLevels, int myRank, int numRanks)
  : NumberOfAMRLevels(numberOfAMRLevels)
{
  int numberOfCells = 0;
  std::vector<int> numberOfCellsperLevel(numberOfAMRLevels,0);
  for (int level = 0; level < numberOfAMRLevels; level++)
  {
    std::array<int, 6> levelIndices;
    levelIndices[0] = 0;                                            // smallest i
    levelIndices[1] = std::pow(2, level);                           // largest i
    levelIndices[2] = 0;                                            // smallest j
    levelIndices[3] = std::pow(2, level);                           // largest j
    levelIndices[4] = level * std::pow(2, level);                   // smallest k
    levelIndices[5] = this->NumberOfAMRLevels * std::pow(2, level); // largest k
    this->LevelIndices.push_back(levelIndices);
    int cellsPerLevel = (levelIndices[1] - levelIndices[0]) * (levelIndices[3] - levelIndices[2]) *
      (levelIndices[5] - levelIndices[4]);
    this->CellsPerLevel.push_back(cellsPerLevel);
    std::array<double, 3> levelOrigin;
    levelOrigin[0] = myRank;
    levelOrigin[1] = 0;
    levelOrigin[2] = level;
    this->LevelOrigin.push_back(levelOrigin);
    numberOfCells += cellsPerLevel;
    numberOfCellsperLevel[level]+=cellsPerLevel;
  }
  std::for_each(numberOfCellsperLevel.begin(), numberOfCellsperLevel.end(), [](int value) {
      std::cout << value << " ";
  });

  // std::cout << std::endl;
  this->BlockId.resize(this->NumberOfAMRLevels, -1);
  for (int level = 0; level < this->NumberOfAMRLevels; level++)
  {
    this->BlockId[level] = myRank * this->NumberOfAMRLevels + level;
  }
}
// AMR::AMR(AMR&  amr2copy){
//   //writedown the code for the copy-constructor
// }

std::array<int, 6> AMR::GetLevelIndices(int level)
{
  return this->LevelIndices[level];
}

std::array<double, 3> AMR::GetLevelOrigin(int level)
{
  return this->LevelOrigin[level];
}
//AMR ends here
///////////////////////////////////////////////////////////////////////////////

///////////////////////////////////////////////////////////////////////////////
//Jul starts here
Jul::Jul(int numberOfAMRLevels, int myRank, int numRanks,
      std::array<double, 3> origin,std::array<double, 3> lengths,//);
      int nx, int ny, std::string topo, int domain, int childrens) :
       origin_(origin), lengths_(lengths),
        domain_(domain),childrens_(childrens)
         {
  int numberOfCells = 0;
  type_ = "rectilinear";//uniform
  topo_ = "topo";
  coordset_ = "coords";
  dx_ = lengths_[0] / nx;
  dy_ = lengths_[1] / ny;
  dz_ = 1;
  int ntotal = nx*ny;
  coordsX = linspace(nx,origin[0],dx_);
  coordsY = linspace(ny,origin[1],dy_);
  coordsZ = linspace(2 ,origin[2], 1 );//[0,1];
  
}

std::vector<double> linspace(int ns, double o, double ds){
  std::vector<double> Vout(ns,0.0);
  Vout[0]=o; 
    for (int i=1;i<ns;++i){
      Vout[i]=Vout[i-1]+ds;
  }
  return Vout;
}


// std::array<int, 6> AMR_new::GetLevelIndices(int level)
// {
//   return this->LevelIndices[level];
// }

// std::array<double, 3> Jul::GetLevelOrigin(int level)
// {
//   return this->LevelOrigin[level];
// }

//Jul ends here
///////////////////////////////////////////////////////////////////////////////

///////////////////////////////////////////////////////////////////////////////
//Lia starts here
//create the childrens!!!!
  Lia::Lia(int numberOfAMRLevels, int myRank, int numRanks,
      std::array<double, 3> origin,std::array<double, 3> lengths,//);
      int nx, int ny, std::string topo, int domain, int childrens) : 
      origin_(origin), lengths_(lengths),
        domain_(domain),parentid_(parentid), 
        childrens_(childrens)
         {
  int numberOfCells = 0;
  type_ = "rectilinear";//uniform
  topo_ = "topo";
  coordset_ = "coords";
  dx_ = lengths_[0] / nx;
  dy_ = lengths_[1] / ny;
  dz_ = 1;
  int ntotal = nx*ny;
  coordsX = linspace(nx,origin[0],dx_);
  coordsY = linspace(ny,origin[1],dy_);
  coordsZ = linspace(2 ,origin[2], 1 );//[0,1];
  parentid_ = parentid;
  
}

std::vector<double> linspace(int ns, double o, double ds){
  std::vector<double> Vout(ns,0.0);
  Vout[0]=o; 
    for (int i=1;i<ns;++i){
      Vout[i]=Vout[i-1]+ds;
  }
  return Vout;
}


// std::array<int, 6> AMR_new::GetLevelIndices(int level)
// {
//   return this->LevelIndices[level];
// }

// std::array<double, 3> Jul::GetLevelOrigin(int level)
// {
//   return this->LevelOrigin[level];
// }

//Lia ends here
///////////////////////////////////////////////////////////////////////////////


///////////////////////////////////////////////////////////////////////////////
//Grid starts here


// Grid::Grid(const unsigned int numPoints[3], const double spacing[3])
// {
//   if (numPoints[0] == 0 || numPoints[1] == 0 || numPoints[2] == 0)
//   {
//     std::cerr << "Must have a non-zero amount of points in each direction.\n";
//   }
//   this->GlobalBounds[0] = this->GlobalBounds[2] = this->GlobalBounds[4] = 0;
//   for (int i = 0; i < 3; i++)
//   {
//     this->GlobalBounds[1 + 2 * i] = (numPoints[i] - 1) * spacing[i];
//   }
//   // in parallel, we do a simple partitioning in the x-direction.
//   int mpiSize = 1;
//   int mpiRank = 0;
//   MPI_Comm_rank(MPI_COMM_WORLD, &mpiRank);
//   MPI_Comm_size(MPI_COMM_WORLD, &mpiSize);

//   unsigned int startXPoint = mpiRank * numPoints[0] / mpiSize;
//   unsigned int endXPoint = (mpiRank + 1) * numPoints[0] / mpiSize;
//   if (mpiSize != mpiRank + 1)
//   {
//     endXPoint++;
//   }

//   // create the points -- slowest in the x and fastest in the z directions
//   double coord[3] = { 0, 0, 0 };
//   for (unsigned int i = startXPoint; i < endXPoint; i++)
//   {
//     coord[0] = i * spacing[0];
//     for (unsigned int j = 0; j < numPoints[1]; j++)
//     {
//       coord[1] = j * spacing[1];
//       for (unsigned int k = 0; k < numPoints[2]; k++)
//       {
//         coord[2] = k * spacing[2];
//         // add the coordinate to the end of the vector
//         std::copy(coord, coord + 3, std::back_inserter(this->Points));
//       }
//     }
//   }
//   // create the hex cells
//   unsigned int cellPoints[8];
//   unsigned int numXPoints = endXPoint - startXPoint;
//   for (unsigned int i = 0; i < numXPoints - 1; i++)
//   {
//     for (unsigned int j = 0; j < numPoints[1] - 1; j++)
//     {
//       for (unsigned int k = 0; k < numPoints[2] - 1; k++)
//       {
//         cellPoints[0] = i * numPoints[1] * numPoints[2] + j * numPoints[2] + k;
//         cellPoints[1] = (i + 1) * numPoints[1] * numPoints[2] + j * numPoints[2] + k;
//         cellPoints[2] = (i + 1) * numPoints[1] * numPoints[2] + (j + 1) * numPoints[2] + k;
//         cellPoints[3] = i * numPoints[1] * numPoints[2] + (j + 1) * numPoints[2] + k;
//         cellPoints[4] = i * numPoints[1] * numPoints[2] + j * numPoints[2] + k + 1;
//         cellPoints[5] = (i + 1) * numPoints[1] * numPoints[2] + j * numPoints[2] + k + 1;
//         cellPoints[6] = (i + 1) * numPoints[1] * numPoints[2] + (j + 1) * numPoints[2] + k + 1;
//         cellPoints[7] = i * numPoints[1] * numPoints[2] + (j + 1) * numPoints[2] + k + 1;
//         std::copy(cellPoints, cellPoints + 8, std::back_inserter(this->Cells));
//       }
//     }
//   }
// }

// size_t Grid::GetNumberOfPoints()
// {
//   return this->Points.size() / 3;
// }

// size_t Grid::GetNumberOfCells()
// {
//   return this->Cells.size() / 8;
// }

// double* Grid::GetPointsArray()
// {
//   if (this->Points.empty())
//   {
//     return nullptr;
//   }
//   return &(this->Points[0]);
// }

// double* Grid::GetPoint(size_t pointId)
// {
//   if (pointId >= this->Points.size())
//   {
//     return nullptr;
//   }
//   return &(this->Points[pointId * 3]);
// }

// unsigned int* Grid::GetCellPoints(size_t cellId)
// {
//   if (cellId >= this->Cells.size())
//   {
//     return nullptr;
//   }
//   return &(this->Cells[cellId * 8]);
// }

// void Grid::GetGlobalBounds(double globalBounds[6])
// {
//   for (int i = 0; i < 6; i++)
//   {
//     globalBounds[i] = this->GlobalBounds[i];
//   }
// }

//Grid ends here
///////////////////////////////////////////////////////////////////////////////