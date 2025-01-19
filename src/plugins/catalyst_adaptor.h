#ifndef CatalystAdaptor_h
#define CatalystAdaptor_h

// This Macro is required due to the deprecation of MPI C++ bindings
// Without it, some function argument casting of MPI::Op::Init(...) between the C and C++ version may fail unexpectedly
// Note: This is somehow not required on the hyd44 system, on my home system this is the only way to make it work
#define OMPI_SKIP_MPICXX 1

#include <catalyst.hpp>
#include <mpi.h>
#include <iostream>
#include <string>
#include <chrono>
#include <cmath>

#ifdef _MGLET_DOUBLE_PRECISION_
typedef double real;
#else
typedef float real;
#endif

struct CatalystConfig {
    const char* file;
    const char* impl;
    const char* path;
    bool is_repr;
    int myid;
};

// defining a struct of all information
struct MgletDataLink {
    // function pointers
    void (*cp_mgdims)(int*,int*,int*,const int*);
    void (*cp_mgbasb)(int*,int*,int*,int*,int*,int*,const int*);
    void (*cp_iterate_grids_lvl)(int*,const int*,const int*);
    void (*cp_get_bbox)(float*,float*,float*,float*,float*,float*,const int*);
    void (*cp_get_parent)(int*, int*);
    void (*cp_get_arrptr)(void*,void*,const int*);
    void (*cp_get_xyzptr)(void*,void*,void*,const int*);
    void (*cp_get_dxyzptr)(void*,void*,void*,const int*);
    void (*cp_get_ddxyzptr)(void*,void*,void*,const int*);
    // data value (not pointers)
    int myid;
    int numprocs;
    int istep;
    int nscal;
    int lvlmin;
    int lvlmax;
};


int get_ngrids_lvl(const MgletDataLink& args, int ilevel)
{
    int lvlcounter = 1;
    int igrid = -1;
    while ( true ) {
        args.cp_iterate_grids_lvl( &igrid, &lvlcounter, &ilevel );
        if (igrid > 0) { 
            lvlcounter++; 
        } else { 
            break;
        }
    }
    return lvlcounter-1;
};



namespace catalyst_adaptor
{

constexpr const char* REPR_DATASET_NAME = "dataout.vthb";
constexpr const char* PIPELINE_TYPE = "io";
constexpr const char* PIPELINE_CHANNEL = "grid";
constexpr unsigned GRID_NAME_SIZE = 8;

/// @brief Print to stdout if the specified rank is the MPI root
/// @tparam ...Args Variadic template arguments
/// @param rank_id MPI rank id
/// @param ...args Arguments to send to stdout
template <typename... Args>
void print_if_root(int rank_id, Args&&... args) {
    if (rank_id != 0) {
        return;
    }
    std::ostringstream oss;
    (oss << ... << args);
    std::cout << oss.str() << std::endl;
}

/// @brief Creates a string with a maximum number of characters which contains a number and is left padded with 0s
/// @param input Number that is contained within the result
/// @param num_chars Maximum number of characters in string
/// @return Padded string with 0s on the left and input on the right
std::string format_with_zeros(unsigned input, unsigned num_chars) {
    std::string input_str = std::to_string(input);
    if (input_str.length() > num_chars) {
        auto num_chars_str = std::to_string(num_chars);
        auto error_message = "Input '" + input_str + "' is larger than allowed number of characters (" + num_chars_str + ")";
        throw std::invalid_argument(error_message);
    }
    int leading_zeros = num_chars - input_str.length();
    std::string result(leading_zeros, '0');
    result += input_str;
    return result;
}

void get_field_ptr(const char* field_name, const int igrid, const MgletDataLink& data_link, real* field_ptr) {
    data_link.cp_get_arrptr( &field_ptr, &field_name, &igrid );
}


void initialize(const CatalystConfig& config)
{
    print_if_root(config.myid, "CATALYST:");
    conduit_cpp::Node init_node;

    // Representative dataset setup
    if (!config.is_repr) {
        print_if_root(config.myid, "  No representative dataset to be created");
    } else {
        print_if_root(config.myid, "  Creating representative dataset '/", REPR_DATASET_NAME, "'");
        init_node["catalyst/pipelines/0/type"].set(PIPELINE_TYPE);
        init_node["catalyst/pipelines/0/filename"].set(REPR_DATASET_NAME);
        init_node["catalyst/pipelines/0/channel"].set(PIPELINE_CHANNEL);
    }

    // Only set catalyst script if the path is not empty
    // The path may be empty if we create a representative dataset since a script is optional then
    if (strcmp(config.file, "") != 0) {
        init_node["catalyst/scripts/script/filename"] = config.file;
    }
    init_node["catalyst_load/implementation"] = config.impl;
    init_node["catalyst_load/search_paths/paraview"] = config.path;    

    catalyst_status err = catalyst_initialize(conduit_cpp::c_node(&init_node));
    if (err == catalyst_status_ok) {
        print_if_root(config.myid, "  catalyst/scripts/script/filename = '", config.file, "'");
        print_if_root(config.myid, "  catalyst_load/implementation = '", config.impl, "'");
        print_if_root(config.myid, "  catalyst_load/search_paths/paraview = '", config.path, "'");
        if (config.is_repr) {
            print_if_root(config.myid, "  catalyst/pipelines/0/type = '", PIPELINE_TYPE, "'");
            print_if_root(config.myid, "  catalyst/pipelines/0/filename = '", REPR_DATASET_NAME, "'");
            print_if_root(config.myid, "  catalyst/pipelines/0/channel = '", PIPELINE_CHANNEL, "'");
        }
        print_if_root(config.myid, "");
    }  else {
        std::cerr << "Failed to initialize Catalyst: " << err << std::endl;
    }
}


// Consistency requirement accross MPI ranks
//  - each MPI rank must manage the same number of grids
//  - otherwise conduit will append an empty .vtp entry which will trigger an import error in paraview
void execute(const MgletDataLink& args)
{
    print_if_root(args.myid, "Executing Catalyst (timestep=", args.istep, ")");

    // Conduit node for the execution
    conduit_cpp::Node exec_node;

    // Add Catalyst state
    auto state = exec_node["catalyst/state"];
    state["timestep"].set(args.istep);
    state["time"].set( (float) args.istep);

    // Opening grid channel as multimesh
    auto channel = exec_node["catalyst/channels/grid"];
    channel["type"].set("amrmesh");
    auto data = channel["data"];

    print_if_root(args.myid, "    Looping grids");
    // Loop levels independent of specific grid
    for (int ilvl = args.lvlmin; ilvl <= args.lvlmax; ilvl++) {
        // Get number of grids on specific ilvl
        int ngridlvl = get_ngrids_lvl(args, ilvl);
        // Loop through all grids on level ilvl
        for (int igrdlvl = 1; igrdlvl <= ngridlvl; igrdlvl++) {
            // Get grid id from ngridlvl'th grid on each level
            int igrid;
            args.cp_iterate_grids_lvl( &igrid, &igrdlvl, &ilvl );
            if (igrid <= 0) {
                continue;
            }

            // Get grid properties
            int kk; int jj; int ii;
            args.cp_mgdims( &kk, &jj, &ii, &igrid );

            // Grid bounding box
            real minx, maxx, miny, maxy, minz, maxz;
            args.cp_get_bbox( &minx, &maxx, &miny, &maxy, &minz, &maxz, &igrid );

            // Pointers to arrays (3D and 1D)
            void *ptr_arr = nullptr;
            // Velocity U
            char const *u_name = "U_C";
            args.cp_get_arrptr( &ptr_arr, &u_name, &igrid );
            real *vel_u = (real*) ptr_arr;
            // Velocity V
            char const *v_name = "V_C";
            args.cp_get_arrptr( &ptr_arr, &v_name, &igrid );
            real *vel_v = (real*) ptr_arr;
            // Velocity W
            char const *w_name = "W_C";
            args.cp_get_arrptr( &ptr_arr, &w_name, &igrid );
            real *vel_w = (real*) ptr_arr;          

            std::string grid_str = format_with_zeros(igrid, GRID_NAME_SIZE);
            std::string block_name = "grid_" + grid_str;
            auto grid_node = data[block_name];
            grid_node["state/domain_id"] = igrid;
            grid_node["state/cycle"] = args.istep;
            grid_node["state/time"] = args.istep;
            grid_node["state/level"] = ilvl;

            // Calculate grid spacing
            real mdx = ( maxx - minx ) / ( ii - 4 );
            real mdy = ( maxy - miny ) / ( jj - 4 );
            real mdz = ( maxz - minz ) / ( kk - 4 );

            // Set grid coordsets
            auto coords = grid_node["coordsets/coords"];
            coords["type"].set("uniform");
            coords["dims/i"].set(ii + 1 - 4);
            coords["dims/j"].set(jj + 1 - 4);
            coords["dims/k"].set(kk + 1 - 4);
            coords["origin/x"].set(minx);
            coords["origin/y"].set(miny);
            coords["origin/z"].set(minz); 
            coords["spacing/dx"].set(mdx);
            coords["spacing/dy"].set(mdy);
            coords["spacing/dz"].set(mdz);

            // Next, add topology
            auto topology = grid_node["topologies/mesh"];
            topology["type"].set("uniform");
            topology["coordset"].set("coords");

            // Finally, add fields.
            auto fields = grid_node["fields"];
            // Number of value entries without boundary layer
            const int nval = (ii - 4) * (jj - 4) * (kk - 4);
            // Velocity in x is cell-data
            fields["u/association"].set("element");
            fields["u/topology"].set("mesh");
            fields["u/volume_dependent"].set("false");
            fields["u/values"].set_external(vel_u, nval, 0, sizeof(real));
            // Velocity in y is cell-data
            fields["v/association"].set("element");
            fields["v/topology"].set("mesh");
            fields["v/volume_dependent"].set("false");
            fields["v/values"].set_external(vel_v, nval, 0, sizeof(float));
            // Velocity in z is cell-data
            fields["w/association"].set("element");
            fields["w/topology"].set("mesh");
            fields["w/volume_dependent"].set("false");
            fields["w/values"].set_external(vel_w, nval, 0, sizeof(float));
        }
    }
    print_if_root(args.myid, "    Passing node to Catalyst");
    catalyst_status err = catalyst_execute(conduit_cpp::c_node(&exec_node));
    if (err != catalyst_status_ok) {
        std::cerr << "Failed to execute Catalyst: " << err << std::endl;
    }
    print_if_root(args.myid, "    Catalyst Status Ok\n");
}

} // namespace catalyst_adaptor

#endif