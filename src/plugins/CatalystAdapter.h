
#ifndef CatalystAdaptor_h
#define CatalystAdaptor_h

#include <catalyst.hpp>

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

    // Passing the parameters from the JSON
    init_node["catalyst/scripts/script/filename"] = file;
    init_node["catalyst_load/implementation"] = impl;
    init_node["catalyst_load/search_paths/paraview"] = path;

    catalyst_status err = catalyst_initialize(conduit_cpp::c_node(&init_node));
    if (err == catalyst_status_ok)
    {
        std::cerr << "Successfully initialized Catalyst" << std::endl;
        std::cerr << " - catalyst/scripts/script/filename = 'paraview'" << std::endl;
        std::cerr << " - catalyst_load/implementation = 'paraview'" << std::endl;
        std::cerr << " - catalyst_load/search_paths/paraview = " << std::endl;
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

    // declaring data vector (must be here - otherwise out-of-scope)
    std::vector<float> Pressure;

    // add cycle information
    auto state = exec_params["catalyst/state"];
    state["timestep"].set(args->istep);
    state["time"].set( (float) args->istep);
    // state["multiblock"].set(1);

    // opening one channel named "grid"
    auto channel = exec_params["catalyst/channels/grid"];
    channel["type"].set("mesh");

    for ( int ilvl = 1; ilvl <= 1; ilvl++ )
    {
        int ngridlvl = get_ngrids_lvl( args, ilvl );
        for ( int igrdlvl = 1; igrdlvl <= ngridlvl; igrdlvl++ )
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
            char const *name = "U";

            args->cp_iterate_grids_lvl( &igrid, &igrdlvl, &ilvl );
            if ( igrid > 0 )
            {
                // calls to MGLET routines

                args->cp_mgdims( &kk, &jj, &ii, &igrid );

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

                args->cp_get_arrptr( &ptr_arr, &name, &igrid );
                args->cp_get_xyzptr( &ptr_x, &ptr_y, &ptr_z, &igrid );
                args->cp_get_dxyzptr( &ptr_dx, &ptr_dy, &ptr_dz, &igrid );
                args->cp_get_ddxyzptr( &ptr_ddx, &ptr_ddy, &ptr_ddz, &igrid );

                float (*val_arr) = (float*) ptr_arr;
                float (*x_arr) = (float*) ptr_x;  // [0 - (ii-1)]
                float (*y_arr) = (float*) ptr_y;  // [0 - (jj-1)]
                float (*z_arr) = (float*) ptr_z;  // [0 - (kk-1)]
                float (*dx_arr) = (float*) ptr_dx;  // [0 - (ii-1)]
                float (*dy_arr) = (float*) ptr_dy;  // [0 - (jj-1)]
                float (*dz_arr) = (float*) ptr_dz;  // [0 - (kk-1)]
                float (*ddx_arr) = (float*) ptr_ddx;  // [0 - (ii-1)]
                float (*ddy_arr) = (float*) ptr_ddy;  // [0 - (jj-1)]
                float (*ddz_arr) = (float*) ptr_ddz;  // [0 - (kk-1)]
#endif

                // now create the mesh.
                int n_zero = 6;
                std::string dataStr = "data/";
                std::string numStr = std::to_string(igrid);
                numStr = std::string(n_zero - std::min(n_zero, (int) numStr.length()), '0') + numStr;
                std::string blockName = dataStr + numStr;

                auto mesh = channel["data"];

                // populate the data node following the Mesh Blueprint [4]
                // [4] https://llnl-conduit.readthedocs.io/en/latest/blueprint_mesh.html

                // start with coordsets
                mesh["coordsets/coords/type"].set("uniform");

                mesh["coordsets/coords/dims/i"].set(ii+1);
                mesh["coordsets/coords/dims/j"].set(jj+1);
                mesh["coordsets/coords/dims/k"].set(kk+1);

                mesh["coordsets/coords/origin/x"].set(minx);
                mesh["coordsets/coords/origin/y"].set(miny);
                mesh["coordsets/coords/origin/z"].set(minz);

                mdx = ( maxx - minx ) / ( ii - 4 );
                mesh["coordsets/coords/spacing/dx"].set(mdx);
                mdy = ( maxy - miny ) / ( jj - 4 );
                mesh["coordsets/coords/spacing/dy"].set(mdy);
                mdz = ( maxz - minz ) / ( kk - 4 );
                mesh["coordsets/coords/spacing/dz"].set(mdz);

                // Next, add topology
                mesh["topologies/mesh/type"].set("uniform");
                mesh["topologies/mesh/coordset"].set("coords");

                // Finally, add fields.
                auto fields = mesh["fields"];

                Pressure.resize( (ii*jj*kk) );
                for ( int a = 0; a < (ii*jj*kk); a++ ){
                    Pressure[a] = (float) a;
                }

                // pressure is cell-data
                fields["pressure/association"].set("element");
                fields["pressure/topology"].set("mesh");
                fields["pressure/volume_dependent"].set("false");
                fields["pressure/values"].set_external(&Pressure[0], (ii*jj*kk), 0, sizeof(float) );

            }
        }
    }

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