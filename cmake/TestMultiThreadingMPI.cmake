message(STATUS "Checking if MPI version supports multi-threading")

find_package(MPI QUIET REQUIRED)
set(CMAKE_REQUIRED_FLAGS ${MPI_C_COMPILE_FLAGS})
set(CMAKE_REQUIRED_INCLUDES ${MPI_C_INCLUDE_PATH})
set(CMAKE_REQUIRED_LIBRARIES ${MPI_C_LIBRARIES})

check_c_source_runs("
#include \"mpi.h\"
int main( int argc, char** argv )
{
    int required = MPI_THREAD_MULTIPLE;
    int provided;

    MPI_Init_thread(&argc, &argv, required, &provided);

    MPI_Finalize();

    if( provided == MPI_THREAD_MULTIPLE ) {
        return 0;
    } else {
        return 1;
    }
}" MPI_MULTITHREAD_OK)

if(NOT MPI_MULTITHREAD_OK)
    message(WARNING "Your MPI version does not support multi-threading (optional)")
endif()
