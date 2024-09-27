
// #include "mglet_cdef.h"


#include "CatalystAdapter.h"


// int get_ngrids_lvl( TransferFromMGLET* args, int ilevel )
// {
//     int lvlcounter = 1;
//     int igrid = -1;
//     while ( true ) {
//         args->cp_iterate_grids_lvl( &igrid, &lvlcounter, &ilevel );
//         if ( igrid > 0 ){ lvlcounter++; } else { break; }
//     }
//     return (lvlcounter-1);
// };


template<class prec>
void show_arrays( int ii, int jj, int kk,
                  prec*** val_arr,
                  prec* x, prec* y, prec* z,
                  prec* dx, prec* dy, prec* dz,
                  prec* ddx, prec* ddy, prec* ddz )
{
    // casting of the array dimensions
    prec (*arr)[jj][kk] = (prec (*)[jj][kk]) val_arr;

    for ( int k = 2; k < kk-2; k++ )
    {
        std::cout << "ARR = " << arr[3][3][k] << std::endl;
    }

};


// void process_arguments( TransferFromMGLET* args )
// {
//     for ( int ilvl = args->lvlmin; ilvl <= args->lvlmax; ilvl++ )
//     {
//         int ngridlvl = get_ngrids_lvl( args, ilvl );
//         for ( int igrdlvl = 1; igrdlvl <= ngridlvl; igrdlvl++ )
//         {
//             // std::cout << igrdlvl << std::endl;

//             // grid properties
//             int igrid; int kk; int jj; int ii;

//             // grid bounding box
//             float minx; float maxx;
//             float miny; float maxy;
//             float minz; float maxz;

//             // pointers to arrays (3D and 1D)
//             void *ptr_arr = nullptr;
//             void *ptr_x = nullptr;
//             void *ptr_y = nullptr;
//             void *ptr_z = nullptr;
//             void *ptr_dx = nullptr;
//             void *ptr_dy = nullptr;
//             void *ptr_dz = nullptr;
//             void *ptr_ddx = nullptr;
//             void *ptr_ddy = nullptr;
//             void *ptr_ddz = nullptr;
//             char const *name = "U";

//             args->cp_iterate_grids_lvl( &igrid, &igrdlvl, &ilvl );
//             if ( igrid > 0 )
//             {
//                 // calls to MGLET routines

//                 args->cp_mgdims( &kk, &jj, &ii, &igrid );
//                 args->cp_get_bbox( &minx, &maxx, &miny, &maxy, &minz, &maxz, &igrid );

//                 args->cp_get_arrptr( &ptr_arr, &name, &igrid );
//                 args->cp_get_xyzptr( &ptr_x, &ptr_y, &ptr_z, &igrid );
//                 args->cp_get_dxyzptr( &ptr_dx, &ptr_dy, &ptr_dz, &igrid );
//                 args->cp_get_ddxyzptr( &ptr_ddx, &ptr_ddy, &ptr_ddz, &igrid );

//                 // casting of arrays

// #ifdef _MGLET_DOUBLE_PRECISION_
//                 // casting arrc[ii][jj][kk] from arrf(kk,jj,ii)
//                 // double (*arr)[jj][kk] = (double (*)[jj][kk]) ptr_arr;
//                 double (***val_arr) = (double***) ptr_arr;
//                 double (*x_arr) = (double*) x_arr;  // [0 - (ii-1)]
//                 double (*y_arr) = (double*) y_arr;  // [0 - (jj-1)]
//                 double (*z_arr) = (double*) z_arr;  // [0 - (kk-1)]
//                 double (*dx_arr) = (double*) ptr_dx;  // [0 - (ii-1)]
//                 double (*dy_arr) = (double*) ptr_dy;  // [0 - (jj-1)]
//                 double (*dz_arr) = (double*) ptr_dz;  // [0 - (kk-1)]
//                 double (*ddx_arr) = (double*) ptr_ddx;  // [0 - (ii-1)]
//                 double (*ddy_arr) = (double*) ptr_ddy;  // [0 - (jj-1)]
//                 double (*ddz_arr) = (double*) ptr_ddz;  // [0 - (kk-1)]
// #else
//                 // casting arrc[ii][jj][kk] from arrf(kk,jj,ii)
//                 // float (*val_arr)[jj][kk] = (float (*)[jj][kk]) ptr_arr;
//                 float (***val_arr) = (float***) ptr_arr;
//                 float (*x_arr) = (float*) ptr_x;  // [0 - (ii-1)]
//                 float (*y_arr) = (float*) ptr_y;  // [0 - (jj-1)]
//                 float (*z_arr) = (float*) ptr_z;  // [0 - (kk-1)]
//                 float (*dx_arr) = (float*) ptr_dx;  // [0 - (ii-1)]
//                 float (*dy_arr) = (float*) ptr_dy;  // [0 - (jj-1)]
//                 float (*dz_arr) = (float*) ptr_dz;  // [0 - (kk-1)]
//                 float (*ddx_arr) = (float*) ptr_ddx;  // [0 - (ii-1)]
//                 float (*ddy_arr) = (float*) ptr_ddy;  // [0 - (jj-1)]
//                 float (*ddz_arr) = (float*) ptr_ddz;  // [0 - (kk-1)]
// #endif

//                 show_arrays( ii, jj, kk, val_arr,
//                              x_arr, y_arr, z_arr,
//                              dx_arr, dy_arr, dz_arr,
//                              ddx_arr, ddy_arr, ddz_arr );

//             }
//         }
//     }
// };


// Main function that is called from MGLET
// (written in C++ but appears as C to the outside)

extern "C" void catalyst_trigger(
    void (*cp_mgdims)(int*,int*,int*,const int*),
    void (*cp_iterate_grids_lvl)(int*,const int*,const int*),
    void (*cp_mgbasb)(int*,int*,int*,int*,int*,int*,const int*),
    void (*cp_get_bbox)(float*,float*,float*,float*,float*,float*,const int*),
    void (*cp_get_arrptr)(void*,void*,const int*),
    void (*cp_get_xyzptr)(void*,void*,void*,const int*),
    void (*cp_get_dxyzptr)(void*,void*,void*,const int*),
    void (*cp_get_ddxyzptr)(void*,void*,void*,const int*),
    int* myid, int* numprocs, int* istep,
    int* nscal, int* lvlmin, int* lvlmax )
{

    // function body -------------------------------------------------

    TransferFromMGLET args;

    // function pointers
    args.cp_mgdims = cp_mgdims;
    args.cp_mgbasb = cp_mgbasb;
    args.cp_iterate_grids_lvl = cp_iterate_grids_lvl;
    args.cp_get_bbox = cp_get_bbox;
    args.cp_get_arrptr = cp_get_arrptr;
    args.cp_get_xyzptr = cp_get_xyzptr;
    args.cp_get_dxyzptr = cp_get_dxyzptr;
    args.cp_get_ddxyzptr = cp_get_ddxyzptr;

    // data pointers (converted to values)
    args.myid = *myid;
    args.numprocs = *numprocs;
    args.istep = *istep;
    args.nscal = *nscal;
    args.lvlmin = *lvlmin;
    args.lvlmax = *lvlmax;

    // std::cout << "Greetings from C++" << std::endl;
    CatalystAdaptor::Execute( &args );

    // function body -------------------------------------------------

    return;
}



extern "C" void catalyst_init( 
    const char* file, 
    const char* impl, 
    const char* path )
{

    // function body -------------------------------------------------

    std::cout << file << std::endl;
    std::cout << impl << std::endl;
    std::cout << path << std::endl;

    CatalystAdaptor::Initialize( file, impl, path );

    // function body -------------------------------------------------

    return;
}
