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
    // conduit::blueprint::mesh::examples::basic("polyhedra",3,3,3,mesh)
    // blueprint::mesh::examples::basic("polyhedra",3,3,3,mesh)
    index_t l = 4;
    blueprint::mesh::examples::julia_nestsets_complex(12,12,0,10,0,10,1,1,l,mesh);
    // std::cout<<mesh.to_yaml()<<std::endl;
    // std::cout<<"\n\n------------info"<<std::endl;

    // std::cout<<mesh.info()<<std::endl;

    // std::cout<<"\n\n------------print"<<std::endl;

    // std::cout<<mesh.print()<<std::endl;
    mesh.print();



    std::string path1="julia_complex_l4";

    relay::io::blueprint::write_mesh(mesh, path1,"json");
    // conduit::relay::io::save(mesh, "iosave", "hdf5");

}