#include <iostream>
#include "conduit.hpp"

#include <conduit_blueprint_mesh_examples_julia.hpp>
#include <conduit_relay_io_blueprint.hpp>



#include "conduit_blueprint.hpp"
#include <conduit_relay_io_blueprint.hpp>
#include <conduit_blueprint_mesh.hpp>
#include <conduit_blueprint_mesh_examples_generate.hpp>
#include <conduit_relay_io.hpp>

using namespace conduit;

int main(){

Node mesh;
// generate simple structured 2d 'basic' mesh
conduit::blueprint::mesh::examples::basic("structured", 3, 3, 10, mesh);
// print out results
std::cout << mesh.to_yaml() << std::endl;
}