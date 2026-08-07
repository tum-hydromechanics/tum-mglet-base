# Distributed under the OSI-approved BSD 3-Clause License.  See accompanying
# file Copyright.txt or https://cmake.org/licensing for details.

cmake_minimum_required(VERSION 3.5)

file(MAKE_DIRECTORY
  "/home/pott/tum-mglet-base/src/_deps/exprtk-src"
  "/home/pott/tum-mglet-base/src/_deps/exprtk-build"
  "/home/pott/tum-mglet-base/src/_deps/exprtk-subbuild/exprtk-populate-prefix"
  "/home/pott/tum-mglet-base/src/_deps/exprtk-subbuild/exprtk-populate-prefix/tmp"
  "/home/pott/tum-mglet-base/src/_deps/exprtk-subbuild/exprtk-populate-prefix/src/exprtk-populate-stamp"
  "/home/pott/tum-mglet-base/src/_deps/exprtk-subbuild/exprtk-populate-prefix/src"
  "/home/pott/tum-mglet-base/src/_deps/exprtk-subbuild/exprtk-populate-prefix/src/exprtk-populate-stamp"
)

set(configSubDirs )
foreach(subDir IN LISTS configSubDirs)
    file(MAKE_DIRECTORY "/home/pott/tum-mglet-base/src/_deps/exprtk-subbuild/exprtk-populate-prefix/src/exprtk-populate-stamp/${subDir}")
endforeach()
