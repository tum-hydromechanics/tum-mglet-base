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
    grid_node["state/level"] = grid.get_level();

    auto coords = grid_node["coordsets/coords"];
    coords["type"] = "uniform";
    auto grid_dims = grid.get_dims();
    auto grid_spacing = grid.get_spacing();
    auto origin_global = grid.get_global_origin();
    coords["dims/i"] = grid_dims[0];
    coords["dims/j"] = grid_dims[1];
    coords["dims/k"] = grid_dims[2];
    coords["spacing/dx"] = grid_spacing[0];
    coords["spacing/dy"] = grid_spacing[1];
    coords["spacing/dz"] = grid_spacing[2];
    coords["origin/x"] = origin_global[0];
    coords["origin/y"] = origin_global[1];
    coords["origin/z"] = origin_global[2];

    auto topo = grid_node["topologies/topo"];
    topo["type"] = "uniform";
    topo["coordset"] = "coords";

    auto fields = grid_node["fields"];
    auto cell_values = fields["cell_values"];
    cell_values["association"] = "element";
    cell_values["topology"] = "topo";
    cell_values["values"].set_float32_vector(grid.get_values());

    int parent_id = grid.get_parent_id();
    if (parent_id >= 0) {
      // This grid has a parent
      const Grid& parent_grid = grids[parent_id];
      conduit_cpp::Node nestset;
      nestset["association"] = "element";
      nestset["topology"] = "topo";

      std::string parent_name = "windows/window_" + std::to_string(parent_id);
      auto parent = nestset[parent_name];
      parent["domain_id"] = parent_id;
      parent["domain_type"] = "parent";
      // Origin of parent-child overlap box, local to the parent
      parent["origin/i"] = 9999;  // Placeholder
      parent["origin/j"] = 9999;  // Placeholder
      parent["origin/k"] = 9999;  // Placeholder
      // Dimensions of parent-child overlap box = child dimensions (since it is always fully contained)
      parent["dims/i"] = grid_dims[0];
      parent["dims/j"] = grid_dims[1];
      parent["dims/k"] = grid_dims[2];
      // Refinement ratio
      parent["ratio/i"] = 9999;  // Placeholder
      parent["ratio/j"] = 9999;  // Placeholder
      parent["ratio/k"] = 9999;  // Placeholder
      grid_node["nestsets/nest"].set(nestset);
    }
    auto child_ids = grid.get_child_ids();
    if (!child_ids.empty()) {
      // This grid has at least one child
      conduit_cpp::Node nestset;
      nestset["association"] = "element";
      nestset["topology"] = "topo";
      for (auto& child_id : child_ids) {
        const Grid& child_grid = grids[child_id];
        std::string child_name = "windows/window_" + std::to_string(child_id);
        auto child = nestset[child_name];
        child["domain_id"] = child_id; // THIS IS THE CHILD ID
        child["domain_type"] = "child";
        // Origin of parent-child overlap box, local to the parent
        child["origin/i"] = 9999;  // Placeholder
        child["origin/j"] = 9999;  // Placeholder
        child["origin/k"] = 9999;  // Placeholder
        // Dimensions of parent-child overlap box = child dimensions (since it is always fully contained)
        auto child_dims = child_grid.get_dims();
        child["dims/i"] = child_dims[0];
        child["dims/j"] = child_dims[1];
        child["dims/k"] = child_dims[2];
        // Refinement ratio
        child["ratio/i"] = 9999;  // Placeholder
        child["ratio/j"] = 9999;  // Placeholder
        child["ratio/k"] = 9999;  // Placeholder
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