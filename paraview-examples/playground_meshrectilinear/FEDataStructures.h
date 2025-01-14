// SPDX-FileCopyrightText: Copyright (c) Kitware Inc.
// SPDX-License-Identifier: BSD-3-Clause
#ifndef FEDataStructures_h
#define FEDataStructures_h

#include <array>
#include <vector>
#include <catalyst.hpp>

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


class Jul
{
public:
  // on each MPI proc the coarse grid is 1 first AMR level cell in the i-dir, 1 first AMR level
  // cell in the j-dir, and numberOfAMRLevels first AMR level cells in the k-dir.
  // The grid gets finer as we go in the k-dir.
  Jul(int numberOfAMRLevels, int myRank, int numRanks,
      std::array<double, 3> origin,std::array<double, 3> lengths,//);
      int nx, int ny, std::string topo, std::string domain,
      int childrens);

  std::array<int, 6> GetLevelIndices(int level);
  std::array<double, 3> GetLevelOrigin(int level){return origin_;}
  std::array<double, 3> GetLengths(int level){return lengths_;}
  std::vector<double> coords(int i);
  int dx(){return dx_;}
  int dy(){return dy_;}
  int nx(){return nx_;}
  int ny(){return ny_;}
  int domain(){return domain_;}
  int level(){return level_;}

  std::string topo() {return topo_;}
  std::string type() {return type_;}

  std::array<double, 3> lengths(){return lengths_;}
  int NumberOfAMRLevels;
  std::vector<std::array<int, 6>> LevelIndices; // inclusive min and max for point indices
  std::vector<int> CellsPerLevel;
  std::vector<int> BlockId; // We only have one child block under each parent block
  // std::vector<std::array<double, 3>> LevelOrigin;
  // private:
  std::array<double, 3> origin_, lengths_,ds;
  int dx_,dy_,dz_,domain_;
  int nx_,ny_,nz_;
  std::string topo_,assotiation_,type_,coordset_;
  int level_;
  std::vector<double> coordsX, coordsY, coordsZ;
  int childrens_;
  // void make_fields(conduit_cpp::Node& proc_id_field);
  // void make_otherfields(conduit_cpp::Node& proc_id_field);

};

class Lia
{
public:
  // on each MPI proc the coarse grid is 1 first AMR level cell in the i-dir, 1 first AMR level
  // cell in the j-dir, and numberOfAMRLevels first AMR level cells in the k-dir.
  // The grid gets finer as we go in the k-dir.
  Lia(int numberOfAMRLevels, int myRank, int numRanks,
      std::array<double, 3> origin,std::array<double, 3> lengths,//);
      int nx, int ny, std::string topo, std::string domain,
      int childrens);

  std::array<int, 6> GetLevelIndices(int level);
  std::array<double, 3> GetLevelOrigin(int level){return origin_;}
  std::array<double, 3> GetLengths(int level){return lengths_;}
  std::vector<double> coords(int i);
  int dx(){return dx_;}
  int dy(){return dy_;}
  int nx(){return nx_;}
  int ny(){return ny_;}
  int domain(){return domain_;}
  int level(){return level_;}

  std::string topo() {return topo_;}
  std::string type() {return type_;}

  std::array<double, 3> lengths(){return lengths_;}
  int NumberOfAMRLevels;
  std::vector<std::array<int, 6>> LevelIndices; // inclusive min and max for point indices
  std::vector<int> CellsPerLevel;
  std::vector<int> BlockId; // We only have one child block under each parent block
  // std::vector<std::array<double, 3>> LevelOrigin;
  // private:
  std::array<double, 3> origin_, lengths_,ds;
  int dx_,dy_,dz_,domain_;
  int nx_,ny_,nz_;
  std::string topo_,assotiation_,type_,coordset_;
  int level_;
  std::vector<double> coordsX, coordsY, coordsZ;
  int childrens_parentid_;


};

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
