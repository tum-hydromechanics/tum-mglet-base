
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

    // add cycle information
    auto state = exec_params["catalyst/state"];
    state["timestep"].set(args->istep);
    state["time"].set( (float) args->istep);
    state["multiblock"].set(1);

    // opening one channel named "grid"
    auto channel = exec_params["catalyst/channels/grid"];
    channel["type"].set("multimesh");

    for ( int ilvl = 2; ilvl <= 2; ilvl++ )
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

                char const *u_name = "U";
                args->cp_get_arrptr( &ptr_arr, &u_name, &igrid );
                float (*vel_u) = (float*) ptr_arr;

                char const *v_name = "V";
                args->cp_get_arrptr( &ptr_arr, &v_name, &igrid );
                float (*vel_v) = (float*) ptr_arr;

                char const *w_name = "W";
                args->cp_get_arrptr( &ptr_arr, &w_name, &igrid );
                float (*vel_w) = (float*) ptr_arr;

                char const *ua_name = "U_AVG";
                args->cp_get_arrptr( &ptr_arr, &ua_name, &igrid );
                float (*vel_u_avg) = (float*) ptr_arr;

                char const *va_name = "V_AVG";
                args->cp_get_arrptr( &ptr_arr, &va_name, &igrid );
                float (*vel_v_avg) = (float*) ptr_arr;

                char const *wa_name = "W_AVG";
                args->cp_get_arrptr( &ptr_arr, &wa_name, &igrid );
                float (*vel_w_avg) = (float*) ptr_arr;                

                args->cp_get_xyzptr( &ptr_x, &ptr_y, &ptr_z, &igrid );
                args->cp_get_dxyzptr( &ptr_dx, &ptr_dy, &ptr_dz, &igrid );
                args->cp_get_ddxyzptr( &ptr_ddx, &ptr_ddy, &ptr_ddz, &igrid );

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

                auto mesh = channel[blockName];

                mdx = ( maxx - minx ) / ( ii - 4 );
                mdy = ( maxy - miny ) / ( jj - 4 );
                mdz = ( maxz - minz ) / ( kk - 4 );

                // start with coordsets
                mesh["coordsets/coords/type"].set("uniform");

                mesh["coordsets/coords/dims/i"].set(kk+1);
                mesh["coordsets/coords/dims/j"].set(jj+1);
                mesh["coordsets/coords/dims/k"].set(ii+1);

                mesh["coordsets/coords/origin/x"].set(minz - 2.0 * mdz);
                mesh["coordsets/coords/origin/y"].set(miny - 2.0 * mdy);
                mesh["coordsets/coords/origin/z"].set(minx - 2.0 * mdx);

                mesh["coordsets/coords/spacing/dx"].set(mdz);
                mesh["coordsets/coords/spacing/dy"].set(mdy);
                mesh["coordsets/coords/spacing/dz"].set(mdx);

                // Next, add topology
                mesh["topologies/mesh/type"].set("uniform");
                mesh["topologies/mesh/coordset"].set("coords");

                // Finally, add fields.
                auto fields = mesh["fields"];

                // velocity_x is cell-data
                fields["u/association"].set("element");
                fields["u/topology"].set("mesh");
                fields["u/volume_dependent"].set("false");
                fields["u/values"].set_external(vel_u, (ii*jj*kk), 0, sizeof(float) );

                // velocity_x is cell-data
                fields["u_avg/association"].set("element");
                fields["u_avg/topology"].set("mesh");
                fields["u_avg/volume_dependent"].set("false");
                fields["u_avg/values"].set_external(vel_u_avg, (ii*jj*kk), 0, sizeof(float) );

                // velocity_y is cell-data
                fields["v/association"].set("element");
                fields["v/topology"].set("mesh");
                fields["v/volume_dependent"].set("false");
                fields["v/values"].set_external(vel_v, (ii*jj*kk), 0, sizeof(float) );

                // velocity_y is cell-data
                fields["v_avg/association"].set("element");
                fields["v_avg/topology"].set("mesh");
                fields["v_avg/volume_dependent"].set("false");
                fields["v_avg/values"].set_external(vel_v_avg, (ii*jj*kk), 0, sizeof(float) );

                // velocity_z is cell-data
                fields["w/association"].set("element");
                fields["w/topology"].set("mesh");
                fields["w/volume_dependent"].set("false");
                fields["w/values"].set_external(vel_w, (ii*jj*kk), 0, sizeof(float) );

                // velocity_z is cell-data
                fields["w_avg/association"].set("element");
                fields["w_avg/topology"].set("mesh");
                fields["w_avg/volume_dependent"].set("false");
                fields["w_avg/values"].set_external(vel_w_avg, (ii*jj*kk), 0, sizeof(float) );

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