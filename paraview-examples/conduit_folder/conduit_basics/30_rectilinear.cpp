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
    // conduit::blueprint::mesh::examples::julia_nestsets_simple
    Node mesh;
    // conduit::blueprint::mesh::examples::basic("polyhedra",3,3,3,mesh);
    conduit::blueprint::mesh::examples::basic("rectilinear", 4, 4, 2, mesh);

    // blueprint::mesh::examples::basic("polyhedra",3,3,3,mesh)
    // blueprint::mesh::examples::julia(8,8,0,10,0,10,1,1,mesh);
    std::cout<<mesh.to_yaml()<<std::endl;


    std::string path1="rectilinearmes";

    relay::io::blueprint::write_mesh(mesh, path1,"json");
    // conduit::relay::io::save(mesh, "iosave", "csv");

}