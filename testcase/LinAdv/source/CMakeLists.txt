cmake_minimum_required(VERSION 3.3)
project(JAUMIN_example CXX C Fortran)
set(CMAKE_C_FLAGS "-std=gnu++0x")
set(CMAKE_CXX_FLAGS "-std=gnu++0x")

if(NOT JAUMIN_ROOT)
  set(JAUMIN_ROOT "/home/zouyg/.local/jasmin-mid")
endif()
if (NOT JAUMIN_NDIM)
  set(JAUMIN_NDIM 3)
endif()
if (NOT JAUMIN_CREL)
  set(JAUMIN_CREL opt)
endif()

list(APPEND CMAKE_MODULE_PATH "${JAUMIN_ROOT}/share/cmake/Modules/scns")
find_package(JAUMIN REQUIRED)
include_directories(${JAUMIN_INCLUDE_DIRS})
add_definitions(${JAUMIN_DEFINITIONS})
add_definitions(-Wall)

file(GLOB_RECURSE src_files ./*.cpp ./*.C ./fortran/${JAUMIN_NDIM}d/*.f ./fortran/*.f)
add_executable(main${JAUMIN_NDIM}d ${src_files})
target_link_libraries(main${JAUMIN_NDIM}d ${JAUMIN_LIBRARIES})
