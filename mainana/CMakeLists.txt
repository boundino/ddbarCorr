cmake_minimum_required(VERSION 3.12.2)
project(ddbar CXX)
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED True)

# add_compile_options(-Wall -Wextra -march=native -O2)
add_compile_options(-Wall -Wextra -g -O0)

list(APPEND CMAKE_PREFIX_PATH $ENV{ROOTSYS})
find_package(ROOT 6.16 REQUIRED)

project(ddbar_savetree)
include_directories(../includes)
add_executable(ddbar_savetree ddbar_savetree.cc)
target_link_libraries(ddbar_savetree PUBLIC ROOT::TreePlayer)

project(lib2dfit)
include_directories(../includes)
add_library(lib2dfit fit2d.cc)
target_link_libraries(lib2dfit PUBLIC ROOT::TreePlayer ROOT::RooFit)

project(ddbar_2dfit)
include_directories(../includes)
add_executable(ddbar_2dfit ddbar_2dfit.cc)
target_link_libraries(ddbar_2dfit PUBLIC lib2dfit)
