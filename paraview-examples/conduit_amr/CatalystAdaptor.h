// SPDX-FileCopyrightText: Copyright (c) Kitware Inc.
// SPDX-License-Identifier: BSD-3-Clause
#ifndef CatalystAdaptor_h
#define CatalystAdaptor_h

#include <catalyst.hpp>

#include "FEDataStructures.h"

#include <iostream>
#include <math.h>
#include <mpi.h>
#include <string>
#include <vector>
#include <cstring>
#include <fstream>
#include <string>

namespace CatalystAdaptor
{
  static std::vector<std::string> filesToValidate;
  bool trigger =true;

/**
 * In this example, we show how to pass overlapping AMR data
 * into Conduit for ParaView Catalyst processing.
 */
void Initialize(int argc, char* argv[])
{
  conduit_cpp::Node node;
  for (int cc = 1; cc < argc; ++cc)
  {
    if (strcmp(argv[cc], "--output") == 0 && (cc + 1) < argc)
    {
      node["catalyst/pipelines/0/type"].set("io");
      node["catalyst/pipelines/0/filename"].set(argv[cc + 1]);
      node["catalyst/pipelines/0/channel"].set("grid");
      ++cc;
    }
    else if (strcmp(argv[cc], "--exists") == 0 && (cc + 1) < argc)
    {
      filesToValidate.push_back(argv[cc + 1]);
      ++cc;
    }
    else
    {
      const auto path = std::string(argv[cc]);
      // note: one can simply add the script file as follows:
      // node["catalyst/scripts/script" + std::to_string(cc - 1)].set_string(path);

      // alternatively, use this form to pass optional parameters to the script.
      const auto name = "catalyst/scripts/script" + std::to_string(cc - 1);
      node[name + "/filename"].set_string(path);
      node[name + "/args"].append().set_string("argument0");
      node[name + "/args"].append().set_string("argument1=12");
      node[name + "/args"].append().set_string("--argument3");
      node[name + "/args"].append().set_string("--channel-name=grid");
    }
  }

  // indicate that we want to load ParaView-Catalyst
  node["catalyst_load/implementation"].set_string("paraview");
  node["catalyst_load/search_paths/paraview"] = PARAVIEW_IMPL_DIR;

  catalyst_status err = catalyst_initialize(conduit_cpp::c_node(&node));
  if (err != catalyst_status_ok)
  {
    std::cerr << "ERROR: Failed to initialize Catalyst: " << err << std::endl;
  }
}

void Execute(unsigned int cycle, double time, AMR& amr)
{
  int numRanks(1), myRank(0);
  MPI_Comm_size(MPI_COMM_WORLD, &numRanks);
  MPI_Comm_rank(MPI_COMM_WORLD, &myRank);

  conduit_cpp::Node exec_params;

  // add time/cycle information
  auto state = exec_params["catalyst/state"];
  state["timestep"].set(cycle);
  state["time"].set(time);

  // Add channels.
  // We only have 1 channel here. Let's name it 'grid'.
  auto channel = exec_params["catalyst/channels/grid"];

  // Since this example is using Conduit Mesh Blueprint to define the mesh,
  // we set the channel's type to "amrmesh".
  channel["type"].set("amrmesh");

  // now create the mesh.
  conduit_cpp::Node mesh = channel["data"];

  for (unsigned int level = 0; level < amr.NumberOfAMRLevels; level++)
  {
    std::string patch_name = "domain_" + std::to_string(level + amr.NumberOfAMRLevels * myRank);
    conduit_cpp::Node patch = mesh[patch_name];
    // add basic state info
    patch["state/domain_id"] = level + amr.NumberOfAMRLevels * myRank;
    patch["state/cycle"] = cycle;
    patch["state/time"] = time;
    patch["state/level"] = level;

    patch["coordsets/coords/type"] = "uniform";
    std::array<int, 6> levelIndices = amr.GetLevelIndices(level);
    patch["coordsets/coords/dims/i"] = levelIndices[1] - levelIndices[0] + 1;
    patch["coordsets/coords/dims/j"] = levelIndices[3] - levelIndices[2] + 1;
    patch["coordsets/coords/dims/k"] = levelIndices[5] - levelIndices[4] + 1;

    patch["coordsets/coords/spacing/dx"] = 1. / std::pow(2, level);
    patch["coordsets/coords/spacing/dy"] = 1. / std::pow(2, level);
    patch["coordsets/coords/spacing/dz"] = 1. / std::pow(2, level);

    std::array<double, 3> levelOrigin = amr.GetLevelOrigin(level);
    patch["coordsets/coords/origin/x"] = levelOrigin[0];
    patch["coordsets/coords/origin/y"] = levelOrigin[1];
    patch["coordsets/coords/origin/z"] = levelOrigin[2];

    // create a rectilinear topology that refs our coordset
    patch["topologies/topo/type"] = "uniform";
    patch["topologies/topo/coordset"] = "coords";

    // add logical elements origin
    patch["topologies/topo/elements/origin/i0"] = levelIndices[0];
    patch["topologies/topo/elements/origin/j0"] = levelIndices[2];
    patch["topologies/topo/elements/origin/k0"] = levelIndices[4];

    conduit_cpp::Node nest_set;
    nest_set["association"] = "element";
    nest_set["topology"] = "topo";
    // If level is not on the root level, parent_id is the level above
    if (level > 0)
    {
      int parent_id = amr.BlockId[level - 1];
      std::string parent_name = "windows/window_" + std::to_string(parent_id);
      conduit_cpp::Node parent = nest_set[parent_name];
      parent["domain_id"] = parent_id;
      parent["domain_type"] = "parent";
      std::array<int, 6> parentLevelIndices = amr.GetLevelIndices(level - 1);
      parent["origin/i"] = levelIndices[0] / 2;
      parent["origin/j"] = parentLevelIndices[2];
      parent["origin/k"] = parentLevelIndices[4];
      
      parent["dims/i"] = parentLevelIndices[1] - levelIndices[0] / 2 + 1;
      parent["dims/j"] = parentLevelIndices[3] - parentLevelIndices[2] + 1;
      ;
      parent["dims/k"] = parentLevelIndices[5] - parentLevelIndices[4] + 1;
      ;
      parent["ratio/i"] = 2;
      parent["ratio/j"] = 2;
      parent["ratio/k"] = 2;
    }
    // If level is not on the leaf level, child_id gets set?
    if (level < amr.NumberOfAMRLevels - 1)
    {
      int child_id = amr.BlockId[level];
      std::string child_name = "windows/window_" + std::to_string(child_id);
      conduit_cpp::Node child = nest_set[child_name];
      child["domain_id"] = child_id;
      child["domain_type"] = "child";

      child["origin/i"] = levelIndices[0];
      child["origin/j"] = levelIndices[2];
      child["origin/k"] = levelIndices[4];

      child["dims/i"] = levelIndices[1] - levelIndices[0] + 1;
      child["dims/j"] = levelIndices[3] - levelIndices[2] + 1;
      child["dims/k"] = levelIndices[5] - levelIndices[4] + 1;

      child["ratio/i"] = 2;
      child["ratio/j"] = 2;
      child["ratio/k"] = 2;
    }
    patch["nestsets/nest"].set(nest_set);
    // add fields
    conduit_cpp::Node fields = patch["fields"];

    // cell data corresponding to MPI process id
    conduit_cpp::Node cell_vals_field = fields["cellvals"];
    cell_vals_field["association"] = "element";
    cell_vals_field["topology"] = "topo";
    int num_cells = (levelIndices[1] - levelIndices[0]) * (levelIndices[3] - levelIndices[2]) *
      (levelIndices[5] - levelIndices[4]);
    std::vector<double> cellValues(num_cells, 0.0);

    for (size_t k = 0; k < levelIndices[5] - levelIndices[4]; k++)
    {
      for (size_t j = 0; j < levelIndices[3] - levelIndices[2]; j++)
      {
        for (size_t i = 0; i < levelIndices[1] - levelIndices[0]; i++)
        {
          size_t l = i + j * (levelIndices[1] - levelIndices[0]) +
            k * (levelIndices[1] - levelIndices[0]) * (levelIndices[3] - levelIndices[2]);
          cellValues[l] = i;
        }
      }
    }

    // we copy the data since cellValues will get deallocated
    cell_vals_field["values"] = cellValues;

    // point data that varies in time and X location.
    conduit_cpp::Node other_field = fields["otherfield"];
    other_field["association"] = "vertex";
    other_field["topology"] = "topo";
    int num_points = (levelIndices[1] - levelIndices[0] + 1) *
      (levelIndices[3] - levelIndices[2] + 1) * (levelIndices[5] - levelIndices[4] + 1);
    std::vector<double> point_values(num_points, 0);
    for (size_t k = 0; k < levelIndices[5] - levelIndices[4] + 1; k++)
    {
      for (size_t j = 0; j < levelIndices[3] - levelIndices[2] + 1; j++)
      {
        for (size_t i = 0; i < levelIndices[1] - levelIndices[0] + 1; i++)
        {
          size_t l = i + j * (levelIndices[1] - levelIndices[0] + 1) +
            k * (levelIndices[1] - levelIndices[0] + 1) * (levelIndices[3] - levelIndices[2] + 1);
          double spacing = 1. / std::pow(2, level);
          double xOrigin = levelOrigin[0];
          point_values[l] = xOrigin + i * spacing * std::cos(time);
        }
      }
    }
    // we copy the data since point_values will get deallocated
    other_field["values"] = point_values;
  
  }
  // if(trigger){

  // std::cout<<"\n\n------------------------------writing mesh. rank"<<myRank<<std::endl;
  // // mesh.print();

  // trigger = false;
  // }
  // MPI_Barrier(MPI_COMM_WORLD);

  // this is using conduit.hpp but we are using this embedded conduit_cpp
  // std::string path1="amr_mesh_1";

  // conduit_cpp::relay::io::blueprint::write_mesh(mesh, path1,"json");
  // so let's try this
  // std::cout<<mesh.to_yaml()<<std::endl;
  // and
  // std::cout<<mesh.write_mesh(path1,"json")<<std::endl;

  // exec_params.print(); // for viewing the Conduit node information

  catalyst_status err = catalyst_execute(conduit_cpp::c_node(&exec_params));
  if (err != catalyst_status_ok)
  {
    std::cerr << "Failed to execute Catalyst: " << err << std::endl;
  }
}

void Finalize()
{
  conduit_cpp::Node node;
  catalyst_status err = catalyst_finalize(conduit_cpp::c_node(&node));
  if (err != catalyst_status_ok)
  {
    std::cerr << "ERROR: Failed to finalize Catalyst: " << err << std::endl;
  }

  for (const auto& fname : filesToValidate)
  {
    std::ifstream istrm(fname.c_str(), std::ios::binary);
    if (!istrm.is_open())
    {
      std::cerr << "ERROR: Failed to open file '" << fname.c_str() << "'." << std::endl;
    }
  }
}
}

#endif
