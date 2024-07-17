cmake_minimum_required(VERSION 3.20)
project(hpss)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/cmake/")

add_library(hpss STATIC)
target_sources(hpss
    PRIVATE
        hpss/hpss.hpp
        hpss/hpss.cpp
)
target_include_directories(hpss PUBLIC .)

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
