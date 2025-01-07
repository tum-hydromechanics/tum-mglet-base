#include <conduit.hpp>
#include <iostream>
#include <conduit_relay_io_blueprint.hpp>
#include <conduit_blueprint_mesh.hpp>
#include <conduit_blueprint_mesh_examples_generate.hpp>

#include <conduit_relay_io.hpp>

using namespace conduit;
//trying to save somehting with hdf5
int main(){
    //https://llnl-conduit.readthedocs.io/en/latest/relay_io.html#relay-i-o-hdf5-interface
    
    // ------------------------------------------------------------------
    // Create a 2D array and show it off.
    // int constexpr rank = 2;
    // int constexpr nrows = 3;
    // int constexpr ncols = 4;
    // int constexpr eltcount = nrows * ncols;
    // double data[eltcount];
    // for (int i = 0; i < eltcount; ++i)
    // {
    //     data[i] = i;
    // }

    // std::cout << "Array, in memory:\n";
    // for (int j = 0; j < nrows; ++j)
    // {
    //     for (int i = 0; i < ncols; ++i)
    //     {
    //         std::cout << std::right << std::setw(4) << data[j * ncols + i];
    //     }
    //     std::cout << std::endl;
    // }

    // // Create an HDF5 file with a 2D array.
    // herr_t status = 0;
    // // HDF5 dimensions are ordered from slowest- to fastest-varying.
    // // This is the same as C and C++ nested arrays and opposite from
    // // many people's geometric intuition.
    // hsize_t hdims[rank]{ nrows, ncols };

    // const char* fname = "t_relay_io_hdf5_read_ndarray.hdf5";
    // hid_t file = H5Fcreate(fname, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);

    // // Create, initialize a dataspace for the dataset
    // hid_t    dataset, dataspace;
    // dataspace = H5Screate_simple(rank, hdims, NULL);

    // // Create, initialize the dataset.  Element type is double.
    // const char* dsname = "twoDarray";
    // dataset = H5Dcreate(file, dsname, H5T_NATIVE_DOUBLE, dataspace, 
    //     H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    // status = H5Dwrite(dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, data);
    // status = H5Dclose(dataset);

    // // close the dataspace and file
    // status = H5Sclose(dataspace);
    // status = H5Fclose(file);

    // std::cout << "\nsaved array to '" << fname << ":" << dsname << "'" << std::endl;

    // // ------------------------------------------------------------------
    // // Now read a subset of that 2D array from the HDF5 file.
    // // Two rows, two columns; total of four elements.
    // int constexpr rnrows = 2;
    // int constexpr rncols = 2;
    // int constexpr reltcount = rnrows * rncols;
    // // As noted earlier, HDF5 orders all dimensions from slowest- to
    // // fastest-varying.  In this two-dimensional example, row (or y-index)
    // // always comes before column (or x-index).  If working with a 3D
    // // dataset, level (or z-index) would come before row.
    // int p_sizes[rank]{ rnrows, rncols };
    // // offset to row 0, column 1
    // int p_offsets[rank]{ 0, 1 };
    // // read every row, every other column
    // int p_strides[rank]{ 1, 2 };
    // // Store pointers to these parameters in the read_opts Node
    // Node read_opts;
    // read_opts["sizes"].set_external(p_sizes, rank);
    // read_opts["offsets"].set_external(p_offsets, rank);
    // read_opts["strides"].set_external(p_strides, rank);

    // std::cout << "\nHDF5 Options for reading the array:" << std::endl;
    // read_opts.print();

    // // Read some of the 2D array in the HDF5 file into an array of doubles
    // Node read_data;
    // double p_data_out[reltcount];
    // read_data.set_external(p_data_out, reltcount);
    // std::string in_path;
    // in_path.append(fname).append(":").append(dsname);
    // conduit::relay::io::hdf5_read(in_path.c_str(), read_opts, read_data);

    // // Show what we read
    // std::cout << "Subset of array, read from '" << in_path << "'" << std::endl;
    // for (int j = 0; j < rnrows; ++j)
    // {
    //     for (int i = 0; i < rncols; ++i)
    //     {
    //         std::cout << std::right << std::setw(8) << p_data_out[j * rncols + i];
    //     }
    //     std::cout << std::endl;
    // }

}