#include <iostream>

#include "conduit.hpp"
#include "conduit_blueprint.hpp"
#include <conduit_relay_io_blueprint.hpp>
#include <conduit_blueprint_mesh.hpp>
#include <conduit_blueprint_mesh_examples_generate.hpp>
#include <conduit_relay_io.hpp>

using namespace conduit;

int main(){
    // Node io_about;
    // relay::io::about(io_about);
    // std::cout << "\nRelay I/O Info and Default Options:" << std::endl;
    // std::cout << io_about.to_yaml() << std::endl;

    // Node &hdf5_opts = io_about["options/hdf5"];
    // // change the default chunking threshold to 
    // // a smaller number to enable compression for
    // // a small array
    // hdf5_opts["chunking/threshold"]  = 2000;
    // hdf5_opts["chunking/chunk_size"] = 2000;

    // std::cout << "\nNew HDF5 I/O Options:" << std::endl;
    // hdf5_opts.print();
    // // set options
    // relay::io::hdf5_set_options(hdf5_opts);

    // int num_vals = 5000;
    // Node n;
    // n["my_values"].set(DataType::float64(num_vals));

    // float64 *v_ptr = n["my_values"].value();
    // for(int i=0; i< num_vals; i++)
    // {
    //     v_ptr[i] = float64(i);
    // }

    // // save using options
    // std::cout << "\nsaving data to 'myoutput_chunked.hdf5' " << std::endl;
    // relay::io::hdf5_save(n,"myoutput_chunked.hdf5");   



    // // Call on each of 4 MPI ranks.
    // conduit::Node mesh, bp_index;
    // conduit::blueprint::mesh::examples::braid("uniform", 10, 10, 10, mesh);
    // char domainFile[1024];
    // sprintf(domainFile, "./bp/bp_%04d.hdf5", rank);
    // conduit::relay::io::save(mesh, domainFile, "hdf5");

    // conduit::blueprint::mpi::mesh::generate_index(mesh,
    //                                             "",
    //                                             bp_index["blueprint_index/mesh"],
    //                                             MPI_COMM_WORLD);
    // bp_index["file_pattern"] = "./bp/bp_%04d.hdf5";
    // bp_index["number_of_files"] = 4;
    // bp_index["number_of_trees"] = 4;
    // bp_index["protocol/name"] = "hdf5";
    // bp_index["protocol/version"] = "0.4.0";
    // bp_index["tree_pattern"] = "/";
    // if(rank == 0)
    //     conduit::relay::io::save(bp_index, "bp.root", "hdf5");


    }