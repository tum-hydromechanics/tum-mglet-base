# ##
## Add Conduit's install path to CMAKE_PREFIX_PATH 
## and use following find_package call to import conduit's targets
##
#
cmake_minimum_required(VERSION 3.13)

find_package(Conduit REQUIRED)
# Check if Conduit was found and print its location
if(Conduit_FOUND)
    message(STATUS "Conduit package found!")
    message(STATUS "Conduit include directories: ${Conduit_INCLUDE_DIRS}")
    message(STATUS "Conduit libraries: ${Conduit_LIBRARIES}")
else()
    message(FATAL_ERROR "Conduit package not found!")
endif()


# find_package(Ascent REQUIRED)
# # Check if Conduit was found and print its location
# if(Ascent_FOUND)
#     message(STATUS "Conduit package found!")
#     message(STATUS "Conduit include directories: ${Conduit_INCLUDE_DIRS}")
#     message(STATUS "Conduit libraries: ${Conduit_LIBRARIES}")
# else()
#     message(FATAL_ERROR "Conduit package not found!")
# endif()


######
# Conduit requires c++14 support
######
set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED ON)


#######################################################
add_executable(conduit_example conduit_example.cpp)
target_link_libraries(conduit_example conduit::conduit)# link to conduit targets
# create our example 
#######################################################

add_executable(hierarchical_const hierarchical_const.cpp)
target_link_libraries(hierarchical_const conduit::conduit)# link to conduit targets
# target_link_libraries(hierarchical_const conduit::conduit ascent::ascent)# link to conduit targets

#######################################################

add_executable(julia julia.cpp)
target_link_libraries(julia conduit::conduit)# link to conduit targets

#######################################################

add_executable(mesh101 mesh101.cpp)
target_link_libraries(mesh101 conduit::conduit)# link to conduit targets

#######################################################

add_executable(mesh_io_1 mesh_io_1.cpp)
target_link_libraries(mesh_io_1 conduit::conduit)# link to conduit targets


#######################################################

add_executable(mesh_io_2 mesh_io_2.cpp)
target_link_libraries(mesh_io_2 conduit::conduit)# link to conduit targets


#######################################################

add_executable(polyhedra polyhedra.cpp)
target_link_libraries(polyhedra conduit::conduit)# link to conduit targets