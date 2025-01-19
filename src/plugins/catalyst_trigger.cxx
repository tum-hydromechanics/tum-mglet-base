
#include "catalyst_adaptor.h"


extern "C" void catalyst_init(const char* file, const char* impl, const char* path, bool* is_repr, int* myid ) {
    CatalystConfig config;
    config.file = file;
    config.impl = impl;
    config.path = path;
    config.is_repr = *is_repr;
    config.myid = *myid;

    catalyst_adaptor::initialize(config);
}

extern "C" void catalyst_finish() {
    conduit_cpp::Node node;

    catalyst_status err = catalyst_finalize(conduit_cpp::c_node(&node));
    if (err != catalyst_status_ok) {
        std::cerr << "Failed to finalize Catalyst: " << err << std::endl;
    }
}

extern "C" void catalyst_trigger(
    void (*cp_mgdims)(int*,int*,int*,const int*),
    void (*cp_iterate_grids_lvl)(int*,const int*,const int*),
    void (*cp_mgbasb)(int*,int*,int*,int*,int*,int*,const int*),
    void (*cp_get_bbox)(float*,float*,float*,float*,float*,float*,const int*),
    void (*cp_get_parent)(int*, int*),
    void (*cp_get_arrptr)(void*,void*,const int*),
    void (*cp_get_xyzptr)(void*,void*,void*,const int*),
    void (*cp_get_dxyzptr)(void*,void*,void*,const int*),
    void (*cp_get_ddxyzptr)(void*,void*,void*,const int*),
    int* myid, int* numprocs, int* istep,
    int* nscal, int* lvlmin, int* lvlmax) {

    MgletDataLink data;
    // Function pointers
    data.cp_mgdims = cp_mgdims;
    data.cp_mgbasb = cp_mgbasb;
    data.cp_iterate_grids_lvl = cp_iterate_grids_lvl;
    data.cp_get_bbox = cp_get_bbox;
    data.cp_get_parent = cp_get_parent;
    data.cp_get_arrptr = cp_get_arrptr;
    data.cp_get_xyzptr = cp_get_xyzptr;
    data.cp_get_dxyzptr = cp_get_dxyzptr;
    data.cp_get_ddxyzptr = cp_get_ddxyzptr;

    // Data pointers (converted to values)
    data.myid = *myid;
    data.numprocs = *numprocs;
    data.istep = *istep;
    data.nscal = *nscal;
    data.lvlmin = *lvlmin;
    data.lvlmax = *lvlmax;

    catalyst_adaptor::execute(data);
}
