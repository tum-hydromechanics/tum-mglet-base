#pragma once

#include <catalyst.hpp>
#include "grid.h"

namespace catalyst_adaptor {

void initialize(int argc, char* argv[]) {
  conduit_cpp::Node node;
  for (int cc = 1; cc < argc; ++cc) {
    if (strcmp(argv[cc], "--output") == 0 && (cc + 1) < argc) {
      node["catalyst/pipelines/0/type"].set("io");
      node["catalyst/pipelines/0/filename"].set(argv[cc + 1]);
      node["catalyst/pipelines/0/channel"].set("grid");
      ++cc;
    } else {
      const auto path = std::string(argv[cc]);
      node["catalyst/scripts/script/filename"] = path;
    }
  }

  node["catalyst_load/implementation"].set_string("paraview");

  catalyst_status err = catalyst_initialize(conduit_cpp::c_node(&node));
  if (err != catalyst_status_ok) {
    std::cerr << "ERROR: Failed to initialize Catalyst: " << err << std::endl;
  }
}



void execute(unsigned timestep, float time, const std::vector<Grid>& grids) {
  conduit_cpp::Node node;
  auto state = node["catalyst/state"];
  state["timestep"].set(timestep);
  state["time"].set(time);

  auto channel = node["catalyst/channels/grid"];
  channel["type"].set("amrmesh");

  auto data = channel["data"];

  for (unsigned i = 0; i < grids.size(); ++i) {
    const Grid& grid = grids[i];

    std::string name = "grid_" + std::to_string(i);
    auto grid_node = data[name];

    grid_node["state/domain_id"] = i;
    grid_node["state/cycle"] = timestep;
    grid_node["state/time"] = time;
    grid_node["state/level"] = grid.level;

    auto coords = grid_node["coordsets/coords"];
    coords["type"] = "uniform";
    coords["dims/i"] = grid.dims[0];
    coords["dims/j"] = grid.dims[1];
    coords["dims/k"] = grid.dims[2];
    coords["spacing/dx"] = grid.spacing[0];
    coords["spacing/dy"] = grid.spacing[1];
    coords["spacing/dz"] = grid.spacing[2];
    coords["origin/x"] = grid.origin[0];
    coords["origin/y"] = grid.origin[1];
    coords["origin/z"] = grid.origin[2];

    auto topo = grid_node["topologies/topo"];
    topo["type"] = "uniform";
    topo["coordset"] = "coords";

    auto fields = grid_node["fields"];
    auto cell_values = fields["cell_values"];
    cell_values["association"] = "element";
    cell_values["topology"] = "topo";
    cell_values["values"].set_float32_vector(grid.values);

    if (grid.parent_id >= 0) {
      // SPECIFY PARENT
      const Grid& parent_grid = grids[grid.parent_id];
      conduit_cpp::Node nestset;
      nestset["association"] = "element";
      nestset["topology"] = "topo";

      std::string parent_name = "windows/window_" + std::to_string(grid.parent_id);
      auto parent = nestset[parent_name];
      parent["domain_id"] = 0; // THIS IS THE PARENT ID
      parent["domain_type"] = "parent";
      parent["origin/i"] = 0; // ON WHICH PARENT ID DOES THE CHILD BEGIN
      parent["origin/j"] = 0;
      parent["origin/k"] = 0;
      parent["dims/i"] = grid.dims[0]; // DIMENSIONS OF CHILD
      parent["dims/j"] = grid.dims[1];
      parent["dims/k"] = grid.dims[2];
      parent["ratio/i"] = 4;
      parent["ratio/j"] = 4;
      parent["ratio/k"] = 4;
      grid_node["nestsets/nest"].set(nestset);
    } 
    
    if (!grid.child_ids.empty()) {
      // WE ARE THE PARENT (on the parent level)
      conduit_cpp::Node nestset;
      nestset["association"] = "element";
      nestset["topology"] = "topo";
      for (auto& child_id : grid.child_ids) {
        const Grid& child_grid = grids[child_id];
        std::string child_name = "windows/window_" + std::to_string(child_id);
        auto child = nestset[child_name];
        child["domain_id"] = child_id; // THIS IS THE CHILD ID
        child["domain_type"] = "child";
        child["origin/i"] = 0; // Origin of overlap
        child["origin/j"] = 0;
        child["origin/k"] = 0;
        child["dims/i"] = child_grid.dims[0]; // DIMENSION OF THE OVERLAP
        child["dims/j"] = child_grid.dims[1];
        child["dims/k"] = child_grid.dims[2];
        child["ratio/i"] = 4;
        child["ratio/j"] = 4;
        child["ratio/k"] = 4;
      }
      grid_node["nestsets/nest"].set(nestset);
    }
  }

  std::cout << node.to_yaml() << std::endl;

  catalyst_status err = catalyst_execute(conduit_cpp::c_node(&node));
  if (err != catalyst_status_ok) {
    std::cerr << "Failed to execute Catalyst: " << err << std::endl;
  }
}

void finalize() {
  conduit_cpp::Node node;
  catalyst_status err = catalyst_finalize(conduit_cpp::c_node(&node));
  if (err != catalyst_status_ok) {
    std::cerr << "ERROR: Failed to finalize Catalyst: " << err << std::endl;
  }
}

} // namespace catalyst_adaptor