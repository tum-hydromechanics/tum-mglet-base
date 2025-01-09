#include <iostream>
#include "conduit.hpp"
#include <conduit_relay_io_blueprint.hpp>
#include <conduit_blueprint_mesh.hpp>
#include <conduit_blueprint_mesh_examples_generate.hpp>
using namespace conduit;


int main(){
    // create a Conduit node to hold our mesh data
    Node mesh;

    // create the coordinate set
    mesh["coordsets/coords/type"] = "uniform";
    mesh["coordsets/coords/dims/i"] = 3;
    mesh["coordsets/coords/dims/j"] = 3;
    // add origin and spacing to the coordset (optional)
    mesh["coordsets/coords/origin/x"] = -10.0;
    mesh["coordsets/coords/origin/y"] = -10.0;
    mesh["coordsets/coords/spacing/dx"] = 10.0;
    mesh["coordsets/coords/spacing/dy"] = 10.0;

    // add the topology
    // this case is simple b/c it's implicitly derived from the coordinate set
    mesh["topologies/topo/type"] = "uniform";
    // reference the coordinate set by name
    mesh["topologies/topo/coordset"] = "coords";

    // add a simple element-associated field 
    mesh["fields/ele_example/association"] =  "element";
    // reference the topology this field is defined on by name
    mesh["fields/ele_example/topology"] =  "topo";
    // set the field values, for this case we have 4 elements
    mesh["fields/ele_example/values"].set(DataType::float64(4));

    float64 *ele_vals_ptr = mesh["fields/ele_example/values"].value();

    for(int i=0;i<4;i++)
    {
        ele_vals_ptr[i] = float64(i);
    }

    // add a simple vertex-associated field 
    mesh["fields/vert_example/association"] =  "vertex";
    // reference the topology this field is defined on by name
    mesh["fields/vert_example/topology"] =  "topo";
    // set the field values, for this case we have 9 vertices
    mesh["fields/vert_example/values"].set(DataType::float64(9));

    float64 *vert_vals_ptr = mesh["fields/vert_example/values"].value();

    for(int i=0;i<9;i++)
    {
        vert_vals_ptr[i] = float64(i);
    }

    // make sure we conform:
    Node verify_info;
    if(!blueprint::mesh::verify(mesh, verify_info))
    {
        std::cout << "Verify failed!" << std::endl;
        verify_info.print();
    }

    // print out results
    std::cout << mesh.to_yaml() << std::endl;

    // save our mesh to a file that can be read by VisIt
    //
    // this will create the file: complete_uniform_mesh_example.root
    // which includes the mesh blueprint index and the mesh data
    conduit::relay::io::blueprint::save_mesh(mesh,
                                            "complete_uniform_mesh_example2",
                                            "json");
}