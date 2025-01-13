// SPDX-FileCopyrightText: Copyright (c) Kitware Inc.
// SPDX-License-Identifier: BSD-3-Clause
#ifndef FEDataStructures_h
#define FEDataStructures_h

#include <array>
#include <vector>

class AMR
{
public:
  // on each MPI proc the coarse grid is 1 first AMR level cell in the i-dir, 1 first AMR level
  // cell in the j-dir, and numberOfAMRLevels first AMR level cells in the k-dir.
  // The grid gets finer as we go in the k-dir.
  AMR(int numberOfAMRLevels, int myRank, int numRanks);
  // AMR(AMR& amr2copy);

  std::array<int, 6> GetLevelIndices(int level);
  std::array<double, 3> GetLevelOrigin(int level);

  int NumberOfAMRLevels;

  std::vector<std::array<int, 6>> LevelIndices; // inclusive min and max for point indices
  std::vector<int> CellsPerLevel;
  std::vector<int> BlockId; // We only have one child block under each parent block
  std::vector<std::array<double, 3>> LevelOrigin;
};

/*
class AMR_new
{
public:
  // on each MPI proc the coarse grid is 1 first AMR level cell in the i-dir, 1 first AMR level
  // cell in the j-dir, and numberOfAMRLevels first AMR level cells in the k-dir.
  // The grid gets finer as we go in the k-dir.
  AMR_new(int numberOfAMRLevels, int myRank, int numRanks,
      std::array<double, 3> AMR_origin,int size1, int size2);

  std::array<int, 6> GetLevelIndices(int level);
  std::array<double, 3> GetLevelOrigin(int level);

  int NumberOfAMRLevels;
  std::vector<std::array<int, 6>> LevelIndices; // inclusive min and max for point indices
  std::vector<int> CellsPerLevel;
  std::vector<int> BlockId; // We only have one child block under each parent block
  std::vector<std::array<double, 3>> LevelOrigin;
  private:
  std::array<double, 3> AMR_origin_;
  int size1_,size2_;
};
*/
// class Grid
// {
// public:
//   Grid(const unsigned int numPoints[3], const double spacing[3]);
//   size_t GetNumberOfPoints();
//   size_t GetNumberOfCells();
//   double* GetPointsArray();
//   double* GetPoint(size_t pointId);
//   unsigned int* GetCellPoints(size_t cellId);
//   void GetGlobalBounds(double globalBounds[6]);

// private:
//   std::vector<double> Points;
//   std::vector<unsigned int> Cells;
//   double GlobalBounds[6];
// };

#endif
