
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


// defining a struct of all information
struct TransferFromMGLET{
    // function pointers
    void (*cp_mgdims)(int*,int*,int*,const int*);
    void (*cp_mgbasb)(int*,int*,int*,int*,int*,int*,const int*);
    void (*cp_iterate_grids_lvl)(int*,const int*,const int*);
    void (*cp_get_bbox)(float*,float*,float*,float*,float*,float*,const int*);
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


int get_ngrids_lvl( TransferFromMGLET* args, int ilevel )
{
    int lvlcounter = 1;
    int igrid = -1;
    while ( true ) {
        args->cp_iterate_grids_lvl( &igrid, &lvlcounter, &ilevel );
        if ( igrid > 0 ){ lvlcounter++; } else { break; }
    }
    return (lvlcounter-1);
};



namespace CatalystAdaptor
{

void Initialize( const char* file, const char* impl, const char* path )
{
    // Conduit node for the initialization
    conduit_cpp::Node init_node;

    
    /*
    init_node["catalyst/pipelines/0/type"].set("io");
    init_node["catalyst/pipelines/0/filename"].set("dataout.vtm");
    init_node["catalyst/pipelines/0/channel"].set("grid");
    
    */
    // Passing the parameters from the JSON
    init_node["catalyst/scripts/script/filename"] = file;
    init_node["catalyst_load/implementation"] = impl;
    init_node["catalyst_load/search_paths/paraview"] = path;

    catalyst_status err = catalyst_initialize(conduit_cpp::c_node(&init_node));
    if (err == catalyst_status_ok)
    {
        std::cerr << "Successfully initialized Catalyst" << std::endl;
        std::cerr << " - catalyst/scripts/script/filename = '" << file << "'" << std::endl;
        std::cerr << " - catalyst_load/implementation = '" << impl << "'" << std::endl;
        std::cerr << " - catalyst_load/search_paths/paraview = '" << path << "'" << std::endl;
        std::cerr << std::endl;
    } 
    else
    {
        std::cerr << "Failed to initialize Catalyst: " << err << std::endl;
    }
}



void Execute(TransferFromMGLET* args)
{
    // Conduit node for the execution
    conduit_cpp::Node exec_params;

    // testing MPI
    // int numRanks(1), myRank(0);
    // MPI_Comm_size(MPI_COMM_WORLD, &numRanks);
    // MPI_Comm_rank(MPI_COMM_WORLD, &myRank);
    // std::cout << numRanks << myRank << std::endl;

    // add cycle information
    auto state = exec_params["catalyst/state"];
    state["timestep"].set(args->istep);
    state["time"].set( (float) args->istep);
    state["multiblock"].set(1);

    // opening one channel named "grid"
    auto channel = exec_params["catalyst/channels/grid"];
    channel["type"].set("multimesh");

    std::cout << "---------- Execute ----------" << std::endl;
    std::cout << "lvlmin=" << args->lvlmin << std::endl;
    std::cout << "lvlmax=" << args->lvlmax << std::endl;
    std::cout << "itstep=" << args->istep << std::endl;
    std::cout << "myid=" << args->myid << std::endl;


    for ( int ilvl = args->lvlmin; ilvl <= args->lvlmax; ilvl++ )
    {
        std::cout << "  ilvl=" << ilvl << std::endl;
        int ngridlvl = get_ngrids_lvl( args, ilvl );

        std::cout << "  ngridlvl=" << ngridlvl << std::endl;

        for ( int igrdlvl = 0; igrdlvl <= ngridlvl; igrdlvl++ )
        {
            // std::cout << igrdlvl << std::endl;

            // grid properties
            int igrid; int kk; int jj; int ii;

            // pointers to arrays (3D and 1D)
            void *ptr_arr = nullptr;
            void *ptr_x = nullptr;
            void *ptr_y = nullptr;
            void *ptr_z = nullptr;
            void *ptr_dx = nullptr;
            void *ptr_dy = nullptr;
            void *ptr_dz = nullptr;
            void *ptr_ddx = nullptr;
            void *ptr_ddy = nullptr;
            void *ptr_ddz = nullptr;

            args->cp_iterate_grids_lvl( &igrid, &igrdlvl, &ilvl );

            std::cout << "    igrid=" << igrid << std::endl;
            std::cout << "    ilvl=" << ilvl << std::endl;


            if ( igrid > 0 )
            {
                // calls to MGLET routines

                args->cp_mgdims( &kk, &jj, &ii, &igrid );
                std::cout << "    kk=" << kk << ", jj=" << jj << ", ii=" << ii << std::endl;

                // casting of arrays

#ifdef _MGLET_DOUBLE_PRECISION_
                // casting arrc[ii][jj][kk] from arrf(kk,jj,ii)
                // double (*arr)[jj][kk] = (double (*)[jj][kk]) ptr_arr;
                double (***val_arr) = (double***) ptr_arr;
                double (*x_arr) = (double*) x_arr;  // [0 - (ii-1)]
                double (*y_arr) = (double*) y_arr;  // [0 - (jj-1)]
                double (*z_arr) = (double*) z_arr;  // [0 - (kk-1)]
                double (*dx_arr) = (double*) ptr_dx;  // [0 - (ii-1)]
                double (*dy_arr) = (double*) ptr_dy;  // [0 - (jj-1)]
                double (*dz_arr) = (double*) ptr_dz;  // [0 - (kk-1)]
                double (*ddx_arr) = (double*) ptr_ddx;  // [0 - (ii-1)]
                double (*ddy_arr) = (double*) ptr_ddy;  // [0 - (jj-1)]
                double (*ddz_arr) = (double*) ptr_ddz;  // [0 - (kk-1)]
                // grid bounding box
                double minx; double maxx;
                double miny; double maxy;
                double minz; double maxz;
                double mdx; double mdy; double mdz;
#else
                // casting arrc[ii][jj][kk] from arrf(kk,jj,ii)
                // float (*val_arr)[jj][kk] = (float (*)[jj][kk]) ptr_arr;
                // float (***val_arr) = (float***) ptr_arr;

                // grid bounding box
                float minx; float maxx;
                float miny; float maxy;
                float minz; float maxz;
                float mdx; float mdy; float mdz;

                args->cp_get_bbox( &minx, &maxx, &miny, &maxy, &minz, &maxz, &igrid );
                std::cout << "    minx=" << minx << ", maxx=" << maxx << std::endl;
                std::cout << "    miny=" << miny << ", maxy=" << maxy << std::endl;
                std::cout << "    minz=" << minz << ", maxz=" << maxz << std::endl;

                char const *u_name = "U_C";
                args->cp_get_arrptr( &ptr_arr, &u_name, &igrid );
                float *vel_u = (float*) ptr_arr;

                char const *v_name = "V_C";
                args->cp_get_arrptr( &ptr_arr, &v_name, &igrid );
                float *vel_v = (float*) ptr_arr;

                char const *w_name = "W_C";
                args->cp_get_arrptr( &ptr_arr, &w_name, &igrid );
                float *vel_w = (float*) ptr_arr;          

                args->cp_get_xyzptr( &ptr_x, &ptr_y, &ptr_z, &igrid );
                args->cp_get_dxyzptr( &ptr_dx, &ptr_dy, &ptr_dz, &igrid );
                args->cp_get_ddxyzptr( &ptr_ddx, &ptr_ddy, &ptr_ddz, &igrid );

                std::cout << "    vel_u=" << vel_u << ", vel_v=" << vel_v << ", vel_w=" << vel_w << std::endl;
                std::cout << "    ptr_x=" << ptr_x << ", ptr_y=" << ptr_y << ", ptr_z=" << ptr_z << std::endl;
                std::cout << "    ptr_dx=" << ptr_dx << ", ptr_dy=" << ptr_dy << ", ptr_dz=" << ptr_dz << std::endl;
                std::cout << "    ptr_ddx=" << ptr_ddx << ", ptr_ddy=" << ptr_ddy << ", ptr_ddz=" << ptr_ddz << std::endl;

                /*
                float (*x_arr) = (float*) ptr_x;  // [0 - (ii-1)]
                float (*y_arr) = (float*) ptr_y;  // [0 - (jj-1)]
                float (*z_arr) = (float*) ptr_z;  // [0 - (kk-1)]
                float (*dx_arr) = (float*) ptr_dx;  // [0 - (ii-1)]
                float (*dy_arr) = (float*) ptr_dy;  // [0 - (jj-1)]
                float (*dz_arr) = (float*) ptr_dz;  // [0 - (kk-1)]
                float (*ddx_arr) = (float*) ptr_ddx;  // [0 - (ii-1)]
                float (*ddy_arr) = (float*) ptr_ddy;  // [0 - (jj-1)]
                float (*ddz_arr) = (float*) ptr_ddz;  // [0 - (kk-1)]
                */
#endif

                // now create the mesh.
                int n_zero = 6;
                std::string dataStr = "data/";
                std::string numStr = std::to_string(igrid);
                numStr = std::string(n_zero - std::min(n_zero, (int) numStr.length()), '0') + numStr;
                std::string blockName = dataStr + numStr;
                std::cout << "    blockName=" << blockName << std::endl;

                auto mesh = channel[blockName];

                mdx = ( maxx - minx ) / ( ii - 4 );
                mdy = ( maxy - miny ) / ( jj - 4 );
                mdz = ( maxz - minz ) / ( kk - 4 );
                std::cout << "    mdx=" << mdx << ", mdy=" << mdy << ", mdz=" << mdz << std::endl;

                // start with coordsets
                mesh["coordsets/coords/type"].set("uniform");

                mesh["coordsets/coords/dims/i"].set(ii+1-4);
                mesh["coordsets/coords/dims/j"].set(jj+1-4);
                mesh["coordsets/coords/dims/k"].set(kk+1-4);

                mesh["coordsets/coords/origin/x"].set(minx);
                mesh["coordsets/coords/origin/y"].set(miny);
                mesh["coordsets/coords/origin/z"].set(minz); 

                mesh["coordsets/coords/spacing/dx"].set(mdx);
                mesh["coordsets/coords/spacing/dy"].set(mdy);
                mesh["coordsets/coords/spacing/dz"].set(mdz);

                // Next, add topology
                mesh["topologies/mesh/type"].set("uniform");
                mesh["topologies/mesh/coordset"].set("coords");

                // Finally, add fields.
                auto fields = mesh["fields"];

                // number of value entries without boundary layer
                const int nval = (ii-4)*(jj-4)*(kk-4);

                std::cout << "    nval=" << nval << std::endl;

                // velocity_x is cell-data
                fields["u/association"].set("element");
                fields["u/topology"].set("mesh");
                fields["u/volume_dependent"].set("false");
                fields["u/values"].set_external(vel_u, nval, 0, sizeof(float) );

                // velocity_y is cell-data
                fields["v/association"].set("element");
                fields["v/topology"].set("mesh");
                fields["v/volume_dependent"].set("false");
                fields["v/values"].set_external(vel_v, nval, 0, sizeof(float) );

                // velocity_z is cell-data
                fields["w/association"].set("element");
                fields["w/topology"].set("mesh");
                fields["w/volume_dependent"].set("false");
                fields["w/values"].set_external(vel_w, nval, 0, sizeof(float) );

                std::cout << "  -- end of grid --" << std::endl;
            }
        }
    }


    std::cout << "TEST" << std::endl;
    // submitting the execution node
    catalyst_status err = catalyst_execute(conduit_cpp::c_node(&exec_params));
    if (err != catalyst_status_ok)
    {
        std::cerr << "Failed to execute Catalyst: " << err << std::endl;
    }
}



void Finalize()
{
    conduit_cpp::Node node;
    


    // finalization by passing an empty node
    catalyst_status err = catalyst_finalize(conduit_cpp::c_node(&node));
    if (err != catalyst_status_ok)
    {
        std::cerr << "Failed to finalize Catalyst: " << err << std::endl;
    }
}

}

#endif