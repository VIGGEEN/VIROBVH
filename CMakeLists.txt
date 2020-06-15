cmake_minimum_required(VERSION 3.16)

project(virobvh)

set(CMAKE_WINDOWS_EXPORT_ALL_SYMBOLS ON)
set(CMAKE_CXX_STANDARD 17)

add_compile_options(/EHsc /W2 /c /O2 /MTd)#/GL

add_subdirectory(ext/libmorton)
include_directories(ext/libmorton/libmorton/include)

add_library(virobvh SHARED virobvh.cpp virobvh.h ZippIterator.hpp capi.h capi.cpp)

target_link_libraries(virobvh PRIVATE libmorton)