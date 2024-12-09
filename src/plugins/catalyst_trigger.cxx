
#include "CatalystAdapter.h"


// template<class prec>
// void show_arrays( int ii, int jj, int kk,
//                   prec*** val_arr,
//                   prec* x, prec* y, prec* z,
//                   prec* dx, prec* dy, prec* dz,
//                   prec* ddx, prec* ddy, prec* ddz )
// {
//     // casting of the array dimensions
//     prec (*arr)[jj][kk] = (prec (*)[jj][kk]) val_arr;

//     for ( int k = 2; k < kk-2; k++ )
//     {
//         std::cout << "ARR = " << arr[3][3][k] << std::endl;
//     }

// };


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
