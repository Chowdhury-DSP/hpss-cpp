cmake_minimum_required(VERSION 3.20)
project(hpss)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/cmake/")

## HPSS Library ##
add_library(hpss STATIC)
target_sources(hpss
    PRIVATE
        hpss-cpp/hpss/hpss.hpp
        hpss-cpp/hpss/hpss.cpp
)
target_include_directories(hpss PUBLIC hpss-cpp)
target_compile_definitions(hpss PRIVATE _USE_MATH_DEFINES=1)

## Accelerated Math Library ##
include(RuntimeSIMDLib)
make_lib_simd_runtime(math_accelerated
    SOURCES
        hpss-cpp/hpss/util/simd_math.hpp
        hpss-cpp/hpss/util/simd_math.cpp
)
target_link_libraries(hpss PRIVATE math_accelerated)

## Tests ##
option(HPSS_BUILD_TESTS OFF)
if(HPSS_BUILD_TESTS)
    include(CPM)
    CPMAddPackage(
        NAME libsndfile
        VERSION 1.0.30
        GITHUB_REPOSITORY libsndfile/libsndfile
        OPTIONS
        "BUILD_PROGRAMS OFF" "BUILD_EXAMPLES OFF" "BUILD_TESTING OFF"
    )
    add_executable(hpss_test test/test.cpp)
    target_link_libraries(hpss_test PRIVATE hpss sndfile)
endif()
