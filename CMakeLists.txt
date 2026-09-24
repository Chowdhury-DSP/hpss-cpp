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

if(NOT TARGET chowdsp_fft)
    include(CPM)
    CPMAddPackage(
        NAME chowdsp_fft
        GIT_TAG main
        GITHUB_REPOSITORY Chowdhury-DSP/chowdsp_fft
    )
endif()
target_link_libraries(hpss PUBLIC chowdsp_fft)

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

    add_executable(median_test test/median_test.cpp)
    target_link_libraries(median_test PRIVATE hpss)

    add_executable(mask_test test/mask_test.cpp)
    target_link_libraries(mask_test PRIVATE hpss)
endif()
