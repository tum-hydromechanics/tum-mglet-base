# ##
## Add Conduit's install path to CMAKE_PREFIX_PATH 
## and use following find_package call to import conduit's targets
##
#
cmake_minimum_required(VERSION 3.13)
project( Conduit_basics)
find_package(Conduit REQUIRED)
# Check if Conduit was found and print its location
if(Conduit_FOUND)
    message(STATUS "Conduit package found!")
    message(STATUS "Conduit include directories: ${Conduit_INCLUDE_DIRS}")
    message(STATUS "Conduit libraries: ${Conduit_LIBRARIES}")
else()
    message(FATAL_ERROR "Conduit package not found!")
endif()

######
# Conduit requires c++14 support
######
set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
# #######################################################
# Get all .cpp files in the current directory
file(GLOB CPP_FILES "${CMAKE_SOURCE_DIR}/*.cpp")

# Loop through each .cpp file and create an executable
foreach(CPP_FILE ${CPP_FILES})
    # Extract the file name without extension
    get_filename_component(EXEC_NAME ${CPP_FILE} NAME_WE)
    
    # Add the executable
    add_executable(${EXEC_NAME} ${CPP_FILE})
    target_link_libraries(${EXEC_NAME} conduit::conduit)# link to conduit targets
endforeach()
