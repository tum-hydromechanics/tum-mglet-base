#pragma once

#include <catalyst.hpp>
#include <mpi.h>

namespace catalyst_adaptor {

void initialize(int argc, char* argv[]) {
  conduit_cpp::Node node;
  
  node["catalyst/pipelines/0/type"].set("io");
  node["catalyst/pipelines/0/filename"].set("dataout.vthb");
  node["catalyst/pipelines/0/channel"].set("grid");

  node["catalyst_load/implementation"].set_string("paraview");

  catalyst_status err = catalyst_initialize(conduit_cpp::c_node(&node));
  if (err != catalyst_status_ok) {
    std::cerr << "ERROR: Failed to initialize Catalyst: " << err << std::endl;
  }
}

void execute(unsigned timestep, float time) {
  int rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  conduit_cpp::Node node;
  auto state = node["catalyst/state"];
  state["timestep"].set(timestep);
  state["time"].set(time);

  auto channel = node["catalyst/channels/grid"];
  channel["type"].set("amrmesh");
  auto data = channel["data"];

  // hardcode the minimal example
  if (rank == 0) {
    // grid node
    auto grid_node = data["grid_1"];
    grid_node["state/domain_id"] = 1;
    grid_node["state/cycle"] = timestep;
    grid_node["state/time"] = time;
    grid_node["state/level"] = 0;
    // coords node
    auto coords = grid_node["coordsets/coords"];
    coords["type"] = "uniform";
    int ii = 33;
    int jj = 33;
    int kk = 2;
    coords["dims/i"] = 33;
    coords["dims/j"] = 33;
    coords["dims/k"] = 2;
    coords["spacing/dx"] = 1.0;
    coords["spacing/dy"] = 1.0;
    coords["spacing/dz"] = 1.0;
    coords["origin/x"] = 0.0;
    coords["origin/y"] = 0.0;
    coords["origin/z"] = 0.0;
    // topo node
    auto topo = grid_node["topologies/topo"];
    topo["type"] = "uniform";
    topo["coordset"] = "coords";
    // fields node
    auto fields = grid_node["fields"];
    auto cell_values = fields["values"];
    cell_values["association"] = "element";
    cell_values["topology"] = "topo";
    std::vector<float> values(ii * jj * kk, 1.0);
    cell_values["values"].set_float32_vector(values);
  } else if (rank == 1) {
    // grid node
    auto grid_node = data["grid_2"];
    grid_node["state/domain_id"] = 2;
    grid_node["state/cycle"] = timestep;
    grid_node["state/time"] = time;
    grid_node["state/level"] = 1;
    // coords node
    auto coords = grid_node["coordsets/coords"];
    coords["type"] = "uniform";
    int ii = 33;
    int jj = 33;
    int kk = 3;
    coords["dims/i"] = 33;
    coords["dims/j"] = 33;
    coords["dims/k"] = 3;
    coords["spacing/dx"] = 0.5;
    coords["spacing/dy"] = 0.5;
    coords["spacing/dz"] = 0.5;
    coords["origin/x"] = 0.0;
    coords["origin/y"] = 0.0;
    coords["origin/z"] = 0.0;
    // topo node
    auto topo = grid_node["topologies/topo"];
    topo["type"] = "uniform";
    topo["coordset"] = "coords";
    // fields node
    auto fields = grid_node["fields"];
    auto cell_values = fields["values"];
    cell_values["association"] = "element";
    cell_values["topology"] = "topo";
    std::vector<float> values(ii * jj * kk, 2.0);
    cell_values["values"].set_float32_vector(values);
  }

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